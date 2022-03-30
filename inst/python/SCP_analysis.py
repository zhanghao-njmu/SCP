def SCVELO(adata=None, h5ad=None, group_by=None, n_jobs=8,
          liner_reduction=None, nonliner_reduction=None,basis=None,
          mode=["deterministic","stochastic","dynamical"],fitting_by="stochastic",
          magic_impute=False,knn=5, t=2,
          min_shared_counts=30, n_pcs=30, n_neighbors=30, approx=True, 
          stream_smooth=None, stream_density=2,
          arrow_size=5, arrow_length=5,arrow_density=0.5,
          denoise=False,denoise_topn=3,kinetics=False,kinetics_topn=100,
          calculate_velocity_genes=False,
          s_genes=None, g2m_genes=None, axis="equal",
          show_plot=True, dpi=300, save=False, dirpath="./", fileprefix=""):
  import matplotlib.pyplot as plt
  import random
  random.seed(12)
  import scvelo as scv
  import scanpy as sc
  import pandas as pd
  import numpy as np
  
  import warnings
  warnings.simplefilter("ignore", category=UserWarning)
  warnings.simplefilter("ignore", category=FutureWarning)
  warnings.simplefilter("ignore", category=DeprecationWarning)

  import os
  prevdir = os.getcwd()
  os.chdir(os.path.expanduser(dirpath))

  try:
    if adata is None and h5ad is None:
      print("adata or h5ad must be provided.")
      exit()
      
    if adata is None:
      adata = scv.read(h5ad)

    if group_by is None :
      print("group_by must be provided.")
      exit()
      
    if liner_reduction is None and nonliner_reduction is None:
      print("liner_reduction or nonliner_reduction must be provided at least one.")
      exit()
      
    if liner_reduction is None:
      sc.pp.pca(adata, n_comps = n_pcs)
      liner_reduction="X_pca"
      
    if basis is None:
      if nonliner_reduction is not None:
        basis=nonliner_reduction
      else:
        basis="basis"
        adata.obsm["basis"]=adata.obsm[liner_reduction][:,0:2]

    mode.append(fitting_by)
    if kinetics is True or denoise is True:
      mode.append("dynamical")
      
    mode=list(set(mode))
    if "dynamical" in mode:
      mode.sort(key = "dynamical".__eq__)
      
    if not fitting_by in ["deterministic","stochastic"]:
      print("'fitting_by' must be one of 'deterministic' and 'stochastic'.")
      exit()
      
    if not all([m in ["deterministic","stochastic","dynamical"] for m in mode]):
      print("Invalid mode name! Must be the 'deterministic', 'stochastic' or 'dynamical'.")
      exit()

    adata.obs[group_by] = adata.obs[group_by].astype(dtype="category")
    ax=scv.pl.proportions(adata,groupby=group_by,save=False,show=False)
    if show_plot is True:
      plt.show() 
    if save:
        plt.savefig('.'.join(filter(None, [fileprefix, "sp_usp_proportions.png"])), dpi=dpi)

    scv.pp.filter_and_normalize(adata, min_shared_counts = min_shared_counts)
         
    if magic_impute is True:
      import magic
      magic_operator = magic.MAGIC(knn=knn, t=t)
      adata.layers["spliced_raw"]=adata.layers["spliced"]
      adata.layers["unspliced_raw"]=adata.layers["unspliced"]
      adata.layers["spliced"] = magic_operator.fit_transform(adata.layers["spliced_raw"])
      adata.layers["unspliced"] = magic_operator.transform(adata.layers["unspliced_raw"])

    scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors = n_neighbors, use_rep=liner_reduction)

    for m in mode:
      vkey_list=[m]
      dk_list=[False]
      gene_subset_list=[None]
      autoscale_list=[True]
      
      if m == "dynamical":
        adata2 = adata[:, adata.var[fitting_by+"_genes"]].copy()
        Ms = adata2.layers["Ms"]
        Mu = adata2.layers["Mu"]
        adata2.layers.clear()
        adata2.layers["Ms"] = Ms
        adata2.layers["Mu"] = Mu
        connectivities = adata2.obsp["connectivities"]
        adata2.obsp.clear()
        adata2.obsp["connectivities"] = connectivities
  
        scv.tl.recover_dynamics(adata2, var_names=fitting_by+"_genes", use_raw=False, n_jobs=n_jobs)
        
        var_add = [i for i in list(adata2.var.columns) if not i in list(adata.var.columns)]
        adata.var=adata.var.merge(adata2.var[var_add], how='left', left_index=True,right_index=True)
        adata.uns["recover_dynamics"]=adata2.uns["recover_dynamics"]
        
        adata.varm["loss"]=np.empty((adata.shape[1],adata2.varm["loss"].shape[1]))
        adata.varm["loss"][:] = np.nan
        adata.varm["loss"][adata.var[fitting_by+"_genes"],:]=adata2.varm["loss"]
        
        empty_layer=np.empty((adata.layers["spliced"].shape))
        empty_layer[:] = np.nan
        adata.layers["fit_t"]=adata.layers["fit_tau"]=adata.layers["fit_tau_"]=empty_layer
        adata.layers["fit_t"][:,adata.var[fitting_by+"_genes"]]=adata2.layers["fit_t"]
        adata.layers["fit_tau"][:,adata.var[fitting_by+"_genes"]]=adata2.layers["fit_tau"]
        adata.layers["fit_tau_"][:,adata.var[fitting_by+"_genes"]]=adata2.layers["fit_tau_"]
        
        if kinetics is True:
          vkey_list.append("dynamical_kinetics")
          dk_list.append(True)
          gene_subset_list.append(None)
          autoscale_list.append(True)
          top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:kinetics_topn]
          scv.tl.differential_kinetic_test(adata, var_names=top_genes, groupby=group_by)
          
        if denoise is True:
          vkey_list.append("dynamical_denoise")
          dk_list.append(False)
          gene_subset_list.append(adata.var['fit_likelihood'].sort_values(ascending=False).index[:denoise_topn])
          autoscale_list.append(False)
          adata.layers[vkey] = adata.layers[m] + np.random.normal(adata.layers[m], scale=adata.layers['Ms'].std(0))
          
      for i in range(len(vkey_list)):
        vkey=vkey_list[i]
        dk=dk_list[i]
        gene_subset=gene_subset_list[i]
        autoscale=autoscale_list[i]
      
        # Velocity graph
        scv.tl.velocity(adata, mode=m, vkey=vkey,diff_kinetics=dk)
        scv.tl.velocity_graph(adata, vkey=vkey, gene_subset=gene_subset, n_neighbors=n_neighbors, approx=approx, n_jobs=n_jobs)
        if m=="dynamical":
          adata.var["velocity_genes"]=adata.var[m+"_genes"]
          adata.layers["velocity"]=adata.layers[m]
          adata.layers["variance_u"]=adata.layers[m+"_u"]
        else:
          adata.var["velocity_gamma"]=adata.var[m+"_gamma"]
          adata.var["velocity_r2"]=adata.var[m+"_r2"]
          adata.var["velocity_genes"]=adata.var[m+"_genes"]
          adata.layers["velocity"]=adata.layers[m]
          adata.layers["variance_velocity"]=adata.layers["variance_"+m]
        
        scv.pl.velocity_graph(adata, vkey=vkey,basis=basis,title=vkey,color=group_by, save=False, show=False)
        plt.axis(axis)
        if show_plot is True:
          plt.show()
        if save:
          plt.savefig('.'.join(filter(None, [fileprefix, vkey+"_graph.png"])), dpi=dpi)
        
        # Velocity embedding
        scv.tl.velocity_embedding(adata, basis=basis, vkey=vkey, autoscale=autoscale)
        scv.pl.velocity_embedding_stream(adata,vkey=vkey,basis=basis, title=vkey, color=group_by, smooth=stream_smooth, density=stream_density, save=False, show=False)
        plt.axis(axis) 
        if show_plot is True:
          plt.show()
        if save:
          plt.savefig('.'.join(filter(None, [fileprefix, vkey+"_stream.png"])), dpi=dpi)
          
        scv.pl.velocity_embedding(adata,vkey=vkey,basis=basis, title=vkey,color=group_by, size=100, arrow_length=arrow_length, arrow_size=arrow_size, density=arrow_density, linewidth=0.2, save=False, show=False)
        plt.axis(axis) 
        if show_plot is True:
          plt.show()
        if save:
          plt.savefig('.'.join(filter(None, [fileprefix, vkey+"_arrow.png"])), dpi=dpi)
          
        scv.pl.velocity_embedding_grid(adata,vkey=vkey,basis=basis, title=vkey, color=group_by, size=100, arrow_length=arrow_length/2, arrow_size=arrow_size/2, density = arrow_density*2, save=False, show=False)
        plt.axis(axis) 
        if show_plot is True:
          plt.show()
        if save:
          plt.savefig('.'.join(filter(None, [fileprefix, vkey+"_embedding_grid.png"])), dpi=dpi)
        
        # Velocity confidence
        scv.tl.velocity_confidence(adata, vkey=vkey)
        scv.pl.scatter(adata, basis=basis, title=vkey+" length", color=vkey+"_length",cmap="coolwarm", perc=[5, 95], save=False, show=False)
        plt.axis(axis) 
        if show_plot is True:
          plt.show()
        if save:
          plt.savefig('.'.join(filter(None, [fileprefix, vkey+"_length.png"])), dpi=dpi)
          
        scv.pl.scatter(adata, basis=basis, title=vkey+" confidence",color=vkey+"_confidence",cmap="coolwarm", perc=[5, 95], save=False, show=False)
        plt.axis(axis) 
        if show_plot is True:
          plt.show()
        if save:
          plt.savefig('.'.join(filter(None, [fileprefix, vkey+"_confidence.png"])), dpi=dpi)
          
        # Terminal states
        for term in ["root_cells","end_points",vkey+'_root_cells',vkey+'_end_points']:
          if term in adata.obs.columns:
            adata.obs.drop(term, axis=1, inplace=True)
          
        scv.tl.terminal_states(adata,vkey=vkey)
        for term in ["root_cells","end_points"]:
          adata.obs[vkey+"_"+term]= adata.obs[term]
          adata.obs.drop(term, axis=1, inplace=True)
        
        # scv.pl.scatter(adata,basis=basis,title=vkey+" terminal_states",color_gradients=[vkey+'_root_cells', vkey+'_end_points'], perc=[5, 95], legend_loc="best", save=False, show=False)
        # plt.axis(axis);plt.show() 
        # if show_plot is True:
        #   plt.show()
        # if save:
        #   plt.savefig('.'.join(filter(None, [fileprefix, vkey+"_terminal_states.png"])), dpi=dpi)
        
        # Pseudotime
        scv.tl.velocity_pseudotime(adata, vkey=vkey,root_key=vkey+'_root_cells',end_key=vkey+'_end_points')
        scv.pl.scatter(adata, basis=basis,title=vkey+" pseudotime", color=vkey+"_pseudotime",cmap="gnuplot", save=False, show=False)
        plt.axis(axis) 
        if show_plot is True:
          plt.show()
        if save:
          plt.savefig('.'.join(filter(None, [fileprefix, vkey+"_pseudotime.png"])), dpi=dpi)
          
        # Latent time
        if m=="dynamical":
          scv.tl.latent_time(adata, vkey=vkey,root_key=vkey+'_root_cells',end_key=vkey+'_end_points')
          scv.pl.scatter(adata, basis=basis,title=vkey+" latent time", color='latent_time',color_map='gnuplot', save=False, show=False)
          plt.axis(axis)
          if show_plot is True:
            plt.show()
          if save:
            plt.savefig('.'.join(filter(None, [fileprefix, vkey+"_latent_time.png"])), dpi=dpi)
        
        # PAGA
        adata.uns["neighbors"]["distances"] = adata.obsp["distances"]
        adata.uns["neighbors"]["connectivities"] = adata.obsp["connectivities"]
        scv.tl.paga(adata, groups=group_by, vkey=vkey,root_key=vkey+'_root_cells',end_key=vkey+'_end_points')
        scv.pl.paga(adata,title=vkey+" PAGA ("+group_by+")", basis=basis, size=50, alpha=0.05,min_edge_width=2, node_size_scale=1.5,legend_loc="best",legend_fontsize=10, save=False, show=False)
        plt.axis(axis) 
        if show_plot is True:
          plt.show()
        if save:
          plt.savefig('.'.join(filter(None, [fileprefix, vkey+"_paga.png"])), dpi=dpi)
    
        # Velocity genes
        if calculate_velocity_genes is True:
          if m != "dynamical":
            scv.tl.rank_velocity_genes(adata, vkey=vkey, groupby=group_by)
            adata.var[vkey+"_score"]=adata.var["spearmans_score"]
            df = scv.get_df(adata.uns["rank_velocity_genes"]["names"])
            adata.uns["rank_"+vkey+"_genes"]=df
          else:
            scv.tl.rank_dynamical_genes(adata, groupby=group_by)
            df = scv.get_df(adata.uns['rank_dynamical_genes']['names'])
            adata.uns["rank_"+vkey+"_genes"]=df
            

          for cluster in df.columns:
            #df[0:1].values.ravel()[:12] ### by row
            ax=scv.pl.scatter(adata, color=group_by, basis=df[cluster].values[:6],vkey=vkey, size=20, linewidth=2, alpha=1, ylabel="cluster: "+cluster+"\nunspliced",
                              add_linfit=True, add_rug=True, add_outline=True, ncols=3, frameon=True, save=False, show=False)
            #[a.axis(axis) for a in ax]
            #if show_plot is True:
            # plt.show() 
            if save:
              plt.savefig('.'.join(filter(None, [fileprefix, cluster, vkey+"_genes1.png"])), dpi=dpi)
              
            ax=scv.pl.velocity(adata, color=group_by, var_names=df[cluster].values[:6],vkey=vkey, size=10, linewidth=2, alpha=1, ylabel="cluster: "+cluster+"\nunspliced",
                               add_outline=True,basis=basis, color_map=["Blues", "YlOrRd"], ncols=2, save=False, show=False)
            #[a.axis(axis) for a in ax]
            #if show_plot is True:
            # plt.show() 
            if save:
              plt.savefig('.'.join(filter(None, [fileprefix, cluster, vkey+"_genes2.png"])), dpi=dpi)
    
    # Cell cycle  
    # if s_genes is not None and g2m_genes is not None:
    scv.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    # scv.pl.scatter(adata, basis=basis, color_gradients=('S_score', 'G2M_score'), perc=[5, 95], smooth=True, legend_loc="best", save=False, show=False)
    # plt.axis(axis);
    # if show_plot is True:
    #  plt.show() 
    # if save:
    #  plt.savefig('.'.join(filter(None, [fileprefix, "cycle_score.png"])), dpi=dpi)
          
    scv.pl.scatter(adata, basis=basis, color='phase',legend_loc="best",save=False, show=False)
    plt.axis(axis) 
    if show_plot is True:
      plt.show()
    if save:    
      plt.savefig('.'.join(filter(None, [fileprefix, "cycle_phase.png"])), dpi=dpi)

  finally:
      os.chdir(prevdir)

  try:
    adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__[
        '_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
  except:
    pass

  return adata

def CellRank(adata=None, h5ad=None, group_by=None, n_jobs=8,
             liner_reduction=None, nonliner_reduction=None,basis=None,
             mode=["deterministic","stochastic","dynamical"],fitting_by="stochastic",
             magic_impute=False,knn=5, t=2,
             min_shared_counts=30, n_pcs=30, n_neighbors=30, approx=True, 
             stream_smooth=None, stream_density=2,
             arrow_size=5, arrow_length=5,arrow_density=0.5,
             s_genes=None, g2m_genes=None,calculate_velocity_genes=False,
             denoise=False,kinetics=False,axis="equal",
             show_plot=True, dpi=300, save=False, dirpath="./", fileprefix=""):
  import matplotlib.pyplot as plt
  import random
  random.seed(12)
  import scvelo as scv
  import cellrank as cr
  import pandas as pd
  import numpy as np
  
  import warnings
  warnings.simplefilter("ignore", category=UserWarning)
  warnings.simplefilter("ignore", category=FutureWarning)
  warnings.simplefilter("ignore", category=DeprecationWarning)

  import os
  prevdir = os.getcwd()
  os.chdir(os.path.expanduser(dirpath))

  try:
    if adata is None and h5ad is None:
      print("adata or h5ad must be provided.")
      exit()
      
    if adata is None:
      adata = scv.read(h5ad)
    # del adata.uns

    if group_by is None :
      print("group_by must be provided.")
      exit()
      
    if liner_reduction is None and nonliner_reduction is None:
      print("liner_reduction or nonliner_reduction must be provided at least one.")
      exit()
      
    if liner_reduction is None:
      sc.pp.pca(adata, n_comps = n_pcs)
      liner_reduction="X_pca"
      
    if basis is None:
      if nonliner_reduction is not None:
        basis=nonliner_reduction
      else:
        basis="basis"
        adata.obsm["basis"]=adata.obsm[liner_reduction][:,0:2]

    mode.append(fitting_by)
    if kinetics is True or denoise is True:
      mode.append("dynamical")
      
    mode=list(set(mode))
    if "dynamical" in mode:
      mode.sort(key = "dynamical".__eq__)
      
    if not fitting_by in ["deterministic","stochastic"]:
      print("'fitting_by' must be one of 'deterministic' and 'stochastic'.")
      exit()
      
    if not all([m in ["deterministic","stochastic","dynamical"] for m in mode]):
      print("Invalid mode name! Must be the 'deterministic', 'stochastic' or 'dynamical'.")
      exit()

    adata.obs[group_by] = adata.obs[group_by].astype(dtype="category")

    if mode[-1]+"_graph" not in adata.obs.keys():
      adata=SCVELO(adata=adata,group_by=group_by, n_jobs=n_jobs,
                   liner_reduction=liner_reduction, nonliner_reduction=nonliner_reduction,basis=basis,
                   mode=mode,fitting_by=fitting_by,magic_impute=magic_impute,knn=knn, t=t,
                   min_shared_counts=min_shared_counts, n_pcs=n_pcs, n_neighbors=n_neighbors, approx=approx, 
                   stream_smooth=stream_smooth, stream_density=stream_density,
                   arrow_size=arrow_size, arrow_length=arrow_length,arrow_density=arrow_density,
                   denoise=denoise,kinetics=kinetics,
                   calculate_velocity_genes=calculate_velocity_genes,
                   s_genes=s_genes, g2m_genes=g2m_genes,
                   show_plot=show_plot, dpi=dpi, save=save, dirpath=dirpath, fileprefix=fileprefix)
    adata.layers["velocity"]= adata.layers[mode[-1]]
    cr.tl.terminal_states(adata, cluster_key=group_by)
    cr.pl.terminal_states(adata)
    cr.tl.initial_states(adata, cluster_key=group_by)
    cr.pl.initial_states(adata)
    cr.tl.lineages(adata)
    cr.pl.lineages(adata, same_plot=False)
    cr.pl.lineages(adata, same_plot=True)
    
    scv.tl.recover_latent_time(adata, root_key="initial_states_probs", end_key="terminal_states_probs")
    scv.tl.paga(adata,groups=group_by,root_key="initial_states_probs",end_key="terminal_states_probs",use_time_prior="velocity_pseudotime")
    cr.pl.cluster_fates(adata,mode="paga_pie",cluster_key=group_by,basis=basis,legend_kwargs={"loc": "top right out"},legend_loc="top left out",node_size_scale=5,edge_width_scale=1,max_edge_width=4,title="directed PAGA")
    if show_plot is True:
      plt.show()
    
    cr.tl.lineage_drivers(adata,cluster_key=group_by)
    cr.pl.lineage_drivers(adata, lineage=adata.obs[group_by].unique()[1], n_genes=4)
    if show_plot is True:
      plt.show()
    
  finally:
      os.chdir(prevdir)

  try:
    adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__[
        '_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
  except:
    pass

  return adata

def PAGA(adata=None, h5ad=None, group_by=None, liner_reduction=None, nonliner_reduction=None,basis=None,
            n_pcs=30,n_neighbors=30, use_rna_velocity=False,vkey="stochastic",
            embedded_with_PAGA=False,paga_layout="fr", threshold=0.1, point_size=20, axis="equal",
            show_plot=True, dpi=300, save=False, dirpath="./", fileprefix=""):
  import matplotlib.pyplot as plt
  import random
  random.seed(11)
  import scanpy as sc

  import warnings
  warnings.simplefilter("ignore", category=UserWarning)
  warnings.simplefilter("ignore", category=FutureWarning)
  warnings.simplefilter("ignore", category=DeprecationWarning)
  
  import os
  prevdir = os.getcwd()
  os.chdir(os.path.expanduser(dirpath))
  
  try:
    if adata is None and h5ad is None:
      print("adata or h5ad must be provided.")
      exit()
      
    if adata is None:
      adata = sc.read(h5ad)
      
    if group_by is None:
      print("group_by must be provided.")
      exit()
      
    if liner_reduction is None and nonliner_reduction is None:
      print("liner_reduction or nonliner_reduction must be provided at least one.")
      exit()
    
    if liner_reduction is None:
      sc.pp.pca(adata, n_comps = n_pcs)
      liner_reduction="X_pca"
      
    if basis is None:
      if nonliner_reduction is not None:
        basis=nonliner_reduction
      else:
        basis="basis"
        adata.obsm["basis"]=adata.obsm[liner_reduction][:,0:2]
        
    if point_size is None:
      point_size = min(100000 / adata.shape[0],20)  
      
    if use_rna_velocity is True:
      adata.uns["velocity_graph"]=adata.uns[vkey+"_graph"]
    # del adata.uns
    
    adata.obs[group_by] = adata.obs[group_by].astype(dtype = "category")
  
    sc.pp.neighbors(adata, n_pcs = n_pcs, use_rep = liner_reduction, n_neighbors = n_neighbors)
    sc.tl.paga(adata, groups = group_by, use_rna_velocity = use_rna_velocity)
    
    if use_rna_velocity is True:
      ax=sc.pl.paga_compare(adata, basis=basis, threshold = threshold, size = point_size, min_edge_width = 1, node_size_scale = 1,
      dashed_edges='connectivities', solid_edges='transitions_confidence', transitions='transitions_confidence',
      title = basis,frameon = False, edges = True, save = False, show = False)
    else:
      ax=sc.pl.paga_compare(adata, basis=basis, threshold = threshold, size = point_size, title = basis,frameon = False, edges = True, save = False, show = False)
    [a.axis(axis) for a in ax] 
    if show_plot is True:
      plt.show()
    if save:
      plt.savefig('.'.join(filter(None, [fileprefix, "paga_compare.png"])), dpi=dpi)
    
    if embedded_with_PAGA is True:
      sc.pl.paga(adata, threshold = threshold, layout = paga_layout, title = "PAGA layout: " + paga_layout, frameon = False, save = False, show = False)
      plt.axis(axis) 
      if show_plot is True:
        plt.show()
      if save:
        plt.savefig('.'.join(filter(None, [fileprefix, "paga_layout.png"])), dpi=dpi)
        
      sc.tl.draw_graph(adata, init_pos = "paga", layout = paga_layout)
      sc.pl.draw_graph(adata, color=group_by, title = "PAGA layout: " + paga_layout, layout = paga_layout,  frameon = False, legend_loc="on data",show = False)
      plt.axis(axis)
      if show_plot is True:
        plt.show()
      if save:
        plt.savefig('.'.join(filter(None, [fileprefix, "paga_graph.png"])), dpi=dpi)
      
      umap2d = sc.tl.umap(adata, init_pos = "paga", n_components=2, copy = True)
      adata.obsm["PAGAUMAP2D"] = umap2d.obsm["X_umap"]
      ax=sc.pl.paga_compare(adata, basis = "PAGAUMAP2D", threshold = threshold, size = point_size, title = "PAGA-initialized UMAP", edges = True, save = False, show = False)
      [a.axis(axis) for a in ax] 
      if show_plot is True:
        plt.show()
      if save:
        plt.savefig('.'.join(filter(None, [fileprefix, "paga_umap.png"])), dpi=dpi)

  finally:
      os.chdir(prevdir)

  try:
    adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__[
        '_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
  except:
    pass

  return adata


def Palantir(adata=None, h5ad=None,group_by=None,
             liner_reduction=None,nonliner_reduction=None,basis=None,
             n_pcs=30,n_neighbors=30,dm_n_components=10,dm_alpha=0,dm_n_eigs=None,
             early_group = None, terminal_groups = None,early_cell=None,terminal_cells=None,
             num_waypoints=1200,scale_components=True,use_early_cell_as_start=False,
             max_iterations=25,n_jobs=8,
             point_size=20,axis="equal",
             show_plot=True, dpi=300, save=False, dirpath="./", fileprefix=""):
  import matplotlib.pyplot as plt
  import matplotlib
  import statistics
  from math import hypot
  import random
  random.seed(11)
  import scanpy as sc
  import numpy as np
  import pandas as pd
  import palantir
  import harmony
  
  import warnings
  warnings.simplefilter("ignore", category=UserWarning)
  warnings.simplefilter("ignore", category=FutureWarning)
  warnings.simplefilter("ignore", category=DeprecationWarning)
  
  import os
  prevdir = os.getcwd()
  os.chdir(os.path.expanduser(dirpath))

  try:
    if adata is None and h5ad is None:
      print("adata or h5ad must be provided.")
      exit()
    
    if adata is None:
      adata = sc.read(h5ad)
    # del adata.uns
    
    if group_by is None and (early_group is not None or terminal_groups is not None):
      print("group_by must be provided.")
      exit()
    
    if liner_reduction is None and nonliner_reduction is None:
      print("liner_reduction or nonliner_reduction must be provided at least one.")
      exit()
      
    if liner_reduction is None:
      sc.pp.pca(adata, n_comps = n_pcs)
      liner_reduction="X_pca"
      
    if basis is None:
      if nonliner_reduction is not None:
        basis=nonliner_reduction
      else:
        basis=liner_reduction
      
    if point_size is None:
      point_size = min(100000 / adata.shape[0],20)  
    
    if early_group is not None and early_cell is None:
      cell = adata.obs[group_by].index.values[adata.obs[group_by]==early_group]
      early_group_cell=adata.obsm[basis][adata.obs[group_by]==early_group,][:, [0, 1]]
      x=statistics.median(early_group_cell[:,0])
      y=statistics.median(early_group_cell[:,1])
      diff=np.array((x - early_group_cell[:,0], y - early_group_cell[:,1]))
      dist=[]
      for i in range(diff.shape[1]):
        dist.append(hypot(diff[0,i],diff[1,i]))
      early_cell=cell[dist.index(min(dist))]
      
    if early_cell is None:
      print("early_cell must be provided.")
      exit()
    else:  
      print("early_cell: ",early_cell)
          
    if terminal_groups is not None and terminal_cells is None:
      terminal_cells=[]
      for n in range(len(terminal_groups)):
        terminal_group=terminal_groups[n]
        cell = adata.obs[group_by].index.values[adata.obs[group_by]==terminal_group]
        terminal_group_cell=adata.obsm[basis][adata.obs[group_by]==terminal_group,][:, [0, 1]]
        x=statistics.median(terminal_group_cell[:,0])
        y=statistics.median(terminal_group_cell[:,1])
        diff=np.array((x - terminal_group_cell[:,0], y - terminal_group_cell[:,1]))
        dist=[]
        for i in range(diff.shape[1]):
          dist.append(hypot(diff[0,i],diff[1,i]))
        terminal_cells.append(cell[dist.index(min(dist))])

    if terminal_cells is None:
      print("terminal_cells: None")
    else:  
      print("terminal_cells: ",terminal_cells)
      
    # if HVF is None:
    #   sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    # else:
    #   df = pd.DataFrame([False] * adata.X.shape[1],columns=["highly_variable"])
    #   df = df.set_index(adata.var_names)
    #   df.highly_variable.iloc[:n_top_genes] = True
    #   df.loc[df['channel'].isin(['sale','fullprice'])]
    #   df.highly_variable.iloc[df.index.isin(HVF)] = True
    #   if "highly_variable" in adata.var.columns:
    #     adata.var.drop('highly_variable', axis=1, inplace=True)
    #   adata.var=adata.var.join(df)
      
    # adata.uns['pca']['variance_ratio']
    # pca_projections=n_comps = np.where(np.cumsum(ad.uns['pca']['variance_ratio']) > 0.85)[0][0]
    # pca_projections, _ = palantir.utils.run_pca(adata, use_hvg=True)

    pca_projections = pd.DataFrame(adata.obsm[liner_reduction][:,:n_pcs], index=adata.obs_names)
    dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=dm_n_components, knn=n_neighbors, alpha=dm_alpha)
    ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=dm_n_eigs)
    pr_res = palantir.core.run_palantir(ms_data=ms_data,early_cell=early_cell,terminal_states=terminal_cells, knn=n_neighbors,num_waypoints=num_waypoints,
    scale_components=scale_components,use_early_cell_as_start=use_early_cell_as_start,max_iterations=max_iterations,n_jobs=n_jobs)

    adata.obsm["palantir_dm"]=dm_res["T"].toarray()
    adata.uns["dm_kernel"]=dm_res["kernel"]
    for term in np.append(pr_res.branch_probs.columns.values, np.array(["palantir_pseudotime","palantir_diff_potential"])):
      if term in adata.obs.columns:
        adata.obs.drop(term, axis=1, inplace=True)
    adata.obs=adata.obs.join(pr_res.pseudotime.to_frame("palantir_pseudotime"))
    adata.obs=adata.obs.join(pr_res.entropy.to_frame("palantir_diff_potential"))
    adata.obs=adata.obs.join(pr_res.branch_probs)
    
    sc.pl.embedding(adata,basis=basis,color="palantir_pseudotime",size = point_size)
    if save:
      plt.savefig('.'.join(filter(None, [fileprefix, "palantir_pseudotime.png"])), dpi=dpi)
      
    sc.pl.embedding(adata,basis=basis,color="palantir_diff_potential",size = point_size)
    if save:
      plt.savefig('.'.join(filter(None, [fileprefix, "palantir_diff_potential.png"])), dpi=dpi)
      
    sc.pl.embedding(adata,basis=basis,color=pr_res.branch_probs.columns.values,size = point_size)
    if save:
      plt.savefig('.'.join(filter(None, [fileprefix, "palantir_probs.png"])), dpi=dpi)

  finally:
      os.chdir(prevdir)

  try:
    adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__[
        '_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
  except:
    pass

  return adata

def Dynamo(adata=None, h5ad=None,group_by=None,
           liner_reduction=None,nonliner_reduction=None,basis=None,
           n_pcs=30,n_neighbors=30,dm_n_components=10,dm_alpha=0,dm_n_eigs=None,
           early_group = None, terminal_groups = None,early_cell=None,terminal_cells=None,
           num_waypoints=1200,scale_components=True,use_early_cell_as_start=False,
           max_iterations=25,n_jobs=8,
           point_size=20,axis="equal",
           show_plot=True, dpi=300,save=False, dirpath="./", fileprefix=""):
  import matplotlib.pyplot as plt
  import matplotlib
  import statistics
  from math import hypot
  import random
  random.seed(11)
  import scanpy as sc
  import numpy as np
  import pandas as pd
  import palantir
  import harmony
  
  import warnings
  warnings.simplefilter("ignore", category=UserWarning)
  warnings.simplefilter("ignore", category=FutureWarning)
  warnings.simplefilter("ignore", category=DeprecationWarning)
  
  import os
  prevdir = os.getcwd()
  os.chdir(os.path.expanduser(dirpath))

  try:
    if adata is None and h5ad is None:
      print("adata or h5ad must be provided.")
      exit()
    
    if adata is None:
      adata = sc.read(h5ad)
    # del adata.uns
    
    if group_by is None and (early_group is not None or terminal_groups is not None):
      print("group_by must be provided.")
      exit()
    
    if liner_reduction is None and nonliner_reduction is None:
      print("liner_reduction or nonliner_reduction must be provided at least one.")
      exit()
      
    if liner_reduction is None:
      sc.pp.pca(adata, n_comps = n_pcs)
      liner_reduction="X_pca"
      
    if basis is None:
      if nonliner_reduction is not None:
        basis=nonliner_reduction
      else:
        basis=liner_reduction
      
    if point_size is None:
      point_size = min(100000 / adata.shape[0],20)  
    
    if early_group is not None and early_cell is None:
      cell = adata.obs[group_by].index.values[adata.obs[group_by]==early_group]
      early_group_cell=adata.obsm[basis][adata.obs[group_by]==early_group,][:, [0, 1]]
      x=statistics.median(early_group_cell[:,0])
      y=statistics.median(early_group_cell[:,1])
      diff=np.array((x - early_group_cell[:,0], y - early_group_cell[:,1]))
      dist=[]
      for i in range(diff.shape[1]):
        dist.append(hypot(diff[0,i],diff[1,i]))
      early_cell=cell[dist.index(min(dist))]
      
    if early_cell is None:
      print("early_cell must be provided.")
      exit()
    else:  
      print("early_cell: ",early_cell)
          
    if terminal_groups is not None and terminal_cells is None:
      terminal_cells=[]
      for n in range(len(terminal_groups)):
        terminal_group=terminal_groups[n]
        cell = adata.obs[group_by].index.values[adata.obs[group_by]==terminal_group]
        terminal_group_cell=adata.obsm[basis][adata.obs[group_by]==terminal_group,][:, [0, 1]]
        x=statistics.median(terminal_group_cell[:,0])
        y=statistics.median(terminal_group_cell[:,1])
        diff=np.array((x - terminal_group_cell[:,0], y - terminal_group_cell[:,1]))
        dist=[]
        for i in range(diff.shape[1]):
          dist.append(hypot(diff[0,i],diff[1,i]))
        terminal_cells.append(cell[dist.index(min(dist))])

    if terminal_cells is None:
      print("terminal_cells: None")
    else:  
      print("terminal_cells: ",terminal_cells)

    pca_projections = pd.DataFrame(adata.obsm[liner_reduction][:,:n_pcs], index=adata.obs_names)
    dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=dm_n_components, knn=n_neighbors, alpha=dm_alpha)
    ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=dm_n_eigs)
    pr_res = palantir.core.run_palantir(ms_data=ms_data,early_cell=early_cell,terminal_states=terminal_cells, knn=n_neighbors,num_waypoints=num_waypoints,
    scale_components=scale_components,use_early_cell_as_start=use_early_cell_as_start,max_iterations=max_iterations,n_jobs=n_jobs)

    adata.obsm["palantir_dm"]=dm_res["T"].toarray()
    adata.uns["dm_kernel"]=dm_res["kernel"]
    for term in np.append(pr_res.branch_probs.columns.values, np.array(["palantir_pseudotime","palantir_diff_potential"])):
      if term in adata.obs.columns:
        adata.obs.drop(term, axis=1, inplace=True)
    adata.obs=adata.obs.join(pr_res.pseudotime.to_frame("palantir_pseudotime"))
    adata.obs=adata.obs.join(pr_res.entropy.to_frame("palantir_diff_potential"))
    adata.obs=adata.obs.join(pr_res.branch_probs)
    
    sc.pl.embedding(adata,basis=basis,color="palantir_pseudotime",size = point_size)
    if save:
      plt.savefig('.'.join(filter(None, [fileprefix, "palantir_pseudotime.png"])), dpi=dpi)
      
    sc.pl.embedding(adata,basis=basis,color="palantir_diff_potential",size = point_size)
    if save:
      plt.savefig('.'.join(filter(None, [fileprefix, "palantir_diff_potential.png"])), dpi=dpi)
      
    sc.pl.embedding(adata,basis=basis,color=pr_res.branch_probs.columns.values,size = point_size)
    if save:
      plt.savefig('.'.join(filter(None, [fileprefix, "palantir_probs.png"])), dpi=dpi)

  finally:
      os.chdir(prevdir)

  try:
    adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__[
        '_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
  except:
    pass

  return adata
