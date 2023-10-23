def SCVELO(adata=None, h5ad=None, group_by=None, palette=None,
          linear_reduction=None, nonlinear_reduction=None,basis=None,
          mode=["deterministic","stochastic","dynamical"],fitting_by="stochastic",
          magic_impute=False,knn=5, t=2,
          min_shared_counts=30, n_pcs=30, n_neighbors=30, 
          stream_smooth=None, stream_density=2,
          arrow_size=5, arrow_length=5,arrow_density=0.5,
          denoise=False,denoise_topn=3,kinetics=False,kinetics_topn=100,
          calculate_velocity_genes=False,top_n=6,n_jobs=1,
          show_plot=True, dpi=300, save=False, dirpath="./", fileprefix=""):
  import matplotlib.pyplot as plt
  import scvelo as scv
  import scanpy as sc
  import numpy as np
  
  import warnings
  warnings.simplefilter("ignore", category=UserWarning)
  warnings.simplefilter("ignore", category=FutureWarning)
  warnings.simplefilter("ignore", category=DeprecationWarning)

  import os
  prevdir = os.getcwd()
  os.chdir(os.path.expanduser(dirpath))
  
  import platform
  if platform.system()=="Windows":
    import sys, multiprocessing, re
    if re.match(pattern=".*pythonw.exe$",string=sys.executable):
      pythonw=sys.executable
    else:
      pythonw=sys.executable.replace("python.exe", "pythonw.exe")
    sys.executable = pythonw
    sys._base_executable = pythonw
    multiprocessing.set_executable(pythonw) 

  try:
    if adata is None and h5ad is None:
      print("adata or h5ad must be provided.")
      exit()
      
    if adata is None:
      adata = scv.read(h5ad)

    if group_by is None :
      print("'group_by' must be provided.")
      exit()
      
    if linear_reduction is None and nonlinear_reduction is None:
      print("linear_reduction or nonlinear_reduction must be provided at least one.")
      exit()
      
    if linear_reduction is None:
      sc.pp.pca(adata, n_comps = n_pcs)
      linear_reduction="X_pca"
      
    if basis is None:
      if nonlinear_reduction is not None:
        basis=nonlinear_reduction
      else:
        basis="basis"
        adata.obsm[linear_reduction+"_basis"]=adata.obsm[linear_reduction][:,0:2]
    scv.pl.utils.check_basis(adata, basis)
     
    if "spliced" not in adata.layers.keys():
      print("'spliced' data must be provided.")
      exit()
      
    if "unspliced" not in adata.layers.keys():
      print("'unspliced' data must be provided.")
      exit()
      
    if type(mode) is str:
      mode=[mode]

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
    scv.pl.proportions(adata,groupby=group_by,save=False,show=False)
    if show_plot is True:
      plt.show() 
    if save:
      plt.savefig('.'.join(filter(None, [fileprefix, "proportions.png"])), dpi=dpi)

    scv.pp.filter_and_normalize(adata, min_shared_counts = min_shared_counts)
         
    if magic_impute is True:
      import magic
      magic_operator = magic.MAGIC(knn=knn, t=t)
      adata.layers["spliced_raw"]=adata.layers["spliced"]
      adata.layers["unspliced_raw"]=adata.layers["unspliced"]
      adata.layers["spliced"] = magic_operator.fit_transform(adata.layers["spliced_raw"])
      adata.layers["unspliced"] = magic_operator.transform(adata.layers["unspliced_raw"])

    scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors = n_neighbors, use_rep=linear_reduction)

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
        scv.tl.velocity_graph(adata, vkey=vkey, gene_subset=gene_subset, n_neighbors=n_neighbors, n_jobs=n_jobs)
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
        
        # scv.pl.velocity_graph(adata, vkey=vkey,basis=basis,title=vkey,color=group_by, save=False, show=False)
        # plt.axis(axis)
        # if show_plot is True:
        #   plt.show()
        # if save:
        #   plt.savefig('.'.join(filter(None, [fileprefix, vkey+"_graph.png"])), dpi=dpi)
        
        # Velocity embedding
        scv.tl.velocity_embedding(adata, basis=basis, vkey=vkey, autoscale=autoscale)
        scv.pl.velocity_embedding_stream(adata,vkey=vkey,basis=basis, title=vkey, color=group_by,palette=palette, smooth=stream_smooth, density=stream_density,legend_loc="none",save=False, show=False)
        if show_plot is True:
          plt.show()
        if save:
          plt.savefig('.'.join(filter(None, [fileprefix, vkey+"_stream.png"])), dpi=dpi)
          
        scv.pl.velocity_embedding(adata,vkey=vkey,basis=basis, title=vkey,color=group_by,palette=palette, arrow_length=arrow_length, arrow_size=arrow_size, density=arrow_density, linewidth=0.3, save=False, show=False)
        if show_plot is True:
          plt.show()
        if save:
          plt.savefig('.'.join(filter(None, [fileprefix, vkey+"_arrow.png"])), dpi=dpi)
          
        scv.pl.velocity_embedding_grid(adata,vkey=vkey,basis=basis, title=vkey, color=group_by,palette=palette, arrow_length=arrow_length/2, arrow_size=arrow_size/2, density = arrow_density*2, save=False, show=False)
        if show_plot is True:
          plt.show()
        if save:
          plt.savefig('.'.join(filter(None, [fileprefix, vkey+"_embedding_grid.png"])), dpi=dpi)
        
        # Velocity confidence
        scv.tl.velocity_confidence(adata, vkey=vkey)
        scv.pl.scatter(adata, basis=basis, title=vkey+" length", color=vkey+"_length",cmap="coolwarm", save=False, show=False)
        if show_plot is True:
          plt.show()
        if save:
          plt.savefig('.'.join(filter(None, [fileprefix, vkey+"_length.png"])), dpi=dpi)
          
        scv.pl.scatter(adata, basis=basis, title=vkey+" confidence",color=vkey+"_confidence",cmap="magma", save=False, show=False)
        if show_plot is True:
          plt.show()
        if save:
          plt.savefig('.'.join(filter(None, [fileprefix, vkey+"_confidence.png"])), dpi=dpi)
          
        # Terminal states
        for term in ["root_cells","end_points",vkey+'_root_cells',vkey+'_end_points']:
          if term in adata.obs.columns:
            adata.obs.drop(term, axis=1, inplace=True)

        scv.tl.terminal_states(adata,vkey=vkey,)
        for term in ["root_cells","end_points"]:
          adata.obs[vkey+"_"+term]= adata.obs[term]
          adata.obs.drop(term, axis=1, inplace=True)

        # scv.pl.scatter(adata,basis=basis,title=vkey+" terminal_states",color_gradients=[vkey+'_root_cells', vkey+'_end_points'], legend_loc="best", save=False, show=False)
        # if show_plot is True:
        #   plt.show()
        # if save:
        #   plt.savefig('.'.join(filter(None, [fileprefix, vkey+"_terminal_states.png"])), dpi=dpi)
        
        # Pseudotime
        scv.tl.velocity_pseudotime(adata, vkey=vkey,root_key=vkey+'_root_cells',end_key=vkey+'_end_points')
        scv.pl.scatter(adata, basis=basis,title=vkey+" pseudotime", color=vkey+"_pseudotime",cmap="cividis", save=False, show=False)
        if show_plot is True:
          plt.show()
        if save:
          plt.savefig('.'.join(filter(None, [fileprefix, vkey+"_pseudotime.png"])), dpi=dpi)
          
        # Latent time
        if m=="dynamical":
          scv.tl.latent_time(adata, vkey=vkey,root_key=vkey+'_root_cells',end_key=vkey+'_end_points')
          scv.pl.scatter(adata, basis=basis,title=vkey+" latent time", color='latent_time',color_map='cividis', save=False, show=False)
          if show_plot is True:
            plt.show()
          if save:
            plt.savefig('.'.join(filter(None, [fileprefix, vkey+"_latent_time.png"])), dpi=dpi)
        
        # PAGA
        adata.uns["neighbors"]["distances"] = adata.obsp["distances"]
        adata.uns["neighbors"]["connectivities"] = adata.obsp["connectivities"]
        scv.tl.paga(adata, groups=group_by, vkey=vkey,root_key=vkey+'_root_cells',end_key=vkey+'_end_points')
        scv.pl.paga(adata,title=vkey+" PAGA ("+group_by+")",node_colors=palette, basis=basis, alpha=0.5,min_edge_width=2, node_size_scale=1.5,legend_loc="none",save=False, show=False)
        if show_plot is True:
          plt.show()
        if save:
          plt.savefig('.'.join(filter(None, [fileprefix, vkey+"_paga.png"])), dpi=dpi)
    
        # Velocity genes
        if calculate_velocity_genes is True:
          if m != "dynamical":
            scv.tl.rank_velocity_genes(adata, vkey=vkey, groupby=group_by)
            adata.var[vkey+"_score"]=adata.var["spearmans_score"]
            df1 = scv.get_df(adata.uns["rank_velocity_genes"]["names"])
            adata.uns["rank_"+vkey+"_genenames"]=df1
            df2 = scv.get_df(adata.uns["rank_velocity_genes"]["scores"])
            adata.uns["rank_"+vkey+"_genescores"]=df2
            del adata.uns["rank_velocity_genes"]
          else:
            scv.tl.rank_dynamical_genes(adata, groupby=group_by)
            df1 = scv.get_df(adata.uns['rank_dynamical_genes']['names'])
            adata.uns["rank_"+vkey+"_genenames"]=df1
            df2 = scv.get_df(adata.uns['rank_dynamical_genes']['scores'])
            adata.uns["rank_"+vkey+"_genescores"]=df2
            del adata.uns["rank_dynamical_genes"]

          for cluster in df1.columns:
            #df1[0:1].values.ravel()[:12] ### by row

            scv.pl.scatter(adata, color=group_by,palette=palette, basis=df1[cluster].values[:top_n],vkey=vkey, size=10, linewidth=2, alpha=1, ylabel="cluster: "+cluster+"\nunspliced",
                            add_linfit=True, add_rug=True, add_outline=True, ncols=3, frameon=True, save=False, show=False)
            if show_plot is True:
              plt.show()
            if save:
              plt.savefig('.'.join(filter(None, [fileprefix, cluster, vkey+"_genes1.png"])), dpi=dpi)
              
            scv.pl.velocity(adata, color=group_by, var_names=df1[cluster].values[:top_n],vkey=vkey, size=10, linewidth=2, alpha=1, ylabel="cluster: "+cluster+"\nunspliced", 
                             add_outline=True,basis=basis, color_map=["Blues", "YlOrRd"], ncols=2, save=False, show=False)
            if show_plot is True:
              plt.show()
            if save:
              plt.savefig('.'.join(filter(None, [fileprefix, cluster, vkey+"_genes2.png"])), dpi=dpi)
    
    # Cell cycle  
    # if s_genes is not None and g2m_genes is not None:
    # scv.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    # scv.pl.scatter(adata, basis=basis, color_gradients=('S_score', 'G2M_score'), smooth=True, legend_loc="best", save=False, show=False)
    # if show_plot is True:
    #  plt.show() 
    # if save:
    #  plt.savefig('.'.join(filter(None, [fileprefix, "cycle_score.png"])), dpi=dpi)
    # scv.pl.scatter(adata, basis=basis, color='phase',legend_loc="best",save=False, show=False)
    # if show_plot is True:
    #   plt.show()
    # if save:    
    #   plt.savefig('.'.join(filter(None, [fileprefix, "cycle_phase.png"])), dpi=dpi)

  finally:
      os.chdir(prevdir)

  try:
    adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__[
        '_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
  except:
    pass

  return adata

def CellRank(adata=None, h5ad=None, group_by=None,palette=None,
             linear_reduction=None, nonlinear_reduction=None,basis=None,
             mode=["deterministic","stochastic","dynamical"],fitting_by="stochastic",
             magic_impute=False,knn=5, t=2,
             min_shared_counts=30, n_pcs=30, n_neighbors=30, 
             stream_smooth=None, stream_density=2,
             arrow_size=5, arrow_length=5,arrow_density=0.5,
             s_genes=None, g2m_genes=None,calculate_velocity_genes=False,
             denoise=False,kinetics=False, n_jobs=1,
             show_plot=True, dpi=300, save=False, dirpath="./", fileprefix=""):
  import matplotlib.pyplot as plt
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
  
  import platform
  if platform.system()=="Windows":
    import sys, multiprocessing, re
    if re.match(pattern=".*pythonw.exe$",string=sys.executable):
      pythonw=sys.executable
    else:
      pythonw=sys.executable.replace("python.exe", "pythonw.exe")
    sys.executable = pythonw
    sys._base_executable = pythonw
    multiprocessing.set_executable(pythonw) 

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
      
    if linear_reduction is None and nonlinear_reduction is None:
      print("linear_reduction or nonlinear_reduction must be provided at least one.")
      exit()
      
    if linear_reduction is None:
      sc.pp.pca(adata, n_comps = n_pcs)
      linear_reduction="X_pca"
      
    if basis is None:
      if nonlinear_reduction is not None:
        basis=nonlinear_reduction
      else:
        basis="basis"
        adata.obsm["basis"]=adata.obsm[linear_reduction][:,0:2]

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
                   linear_reduction=linear_reduction, nonlinear_reduction=nonlinear_reduction,basis=basis,
                   mode=mode,fitting_by=fitting_by,magic_impute=magic_impute,knn=knn, t=t,
                   min_shared_counts=min_shared_counts, n_pcs=n_pcs, n_neighbors=n_neighbors,
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

def PAGA(adata=None, h5ad=None, group_by=None,palette=None, linear_reduction=None, nonlinear_reduction=None,basis=None,
            n_pcs=30,n_neighbors=30, use_rna_velocity=False,vkey="stochastic",
            embedded_with_PAGA=False,paga_layout="fr", threshold=0.1, point_size=20, 
            infer_pseudotime = False,root_cell = None,root_group=None,n_dcs = 10, n_branchings = 0, min_group_size = 0.01,
            show_plot=True, dpi=300, save=False, dirpath="./", fileprefix=""):
  import matplotlib.pyplot as plt
  import scanpy as sc
  import numpy as np
  import statistics
  from math import hypot

  import warnings
  warnings.simplefilter("ignore", category=UserWarning)
  warnings.simplefilter("ignore", category=FutureWarning)
  warnings.simplefilter("ignore", category=DeprecationWarning)
  
  import os
  prevdir = os.getcwd()
  os.chdir(os.path.expanduser(dirpath))
  
  import platform
  if platform.system()=="Windows":
    import sys, multiprocessing, re
    if re.match(pattern=".*pythonw.exe$",string=sys.executable):
      pythonw=sys.executable
    else:
      pythonw=sys.executable.replace("python.exe", "pythonw.exe")
    sys.executable = pythonw
    sys._base_executable = pythonw
    multiprocessing.set_executable(pythonw) 
  
  try:
    if adata is None and h5ad is None:
      print("adata or h5ad must be provided.")
      exit()
      
    if adata is None:
      adata = sc.read(h5ad)
      
    if group_by is None:
      print("group_by must be provided.")
      exit()
      
    if linear_reduction is None and nonlinear_reduction is None:
      print("linear_reduction or nonlinear_reduction must be provided at least one.")
      exit()
    
    if linear_reduction is None:
      sc.pp.pca(adata, n_comps = n_pcs)
      linear_reduction="X_pca"
      
    if basis is None:
      if nonlinear_reduction is not None:
        basis=nonlinear_reduction
      else:
        basis="basis"
        adata.obsm["basis"]=adata.obsm[linear_reduction][:,0:2]
        
    if point_size is None:
      point_size = min(100000 / adata.shape[0],20)  
      
    if infer_pseudotime is True and root_cell is None and root_group is None:
      print("root_cell or root_group should be provided.")
      exit()

    if use_rna_velocity is True:
      adata.uns["velocity_graph"]=adata.uns[vkey+"_graph"]
    # del adata.uns
    
    adata.obs[group_by] = adata.obs[group_by].astype(dtype = "category")
    
    if "X_diffmap" in adata.obsm_keys():
      X_diffmap = adata.obsm['X_diffmap']
      del adata.obsm['X_diffmap']
      sc.pp.neighbors(adata, n_pcs = n_pcs, use_rep = linear_reduction, n_neighbors = n_neighbors)
      adata.obsm['X_diffmap'] = X_diffmap
    else:
      sc.pp.neighbors(adata, n_pcs = n_pcs, use_rep = linear_reduction, n_neighbors = n_neighbors)
    
    sc.tl.paga(adata, groups = group_by, use_rna_velocity = use_rna_velocity)
    
    if use_rna_velocity is True:
      sc.pl.paga_compare(adata, basis=basis,palette=palette, threshold = threshold, size = point_size, min_edge_width = 1, node_size_scale = 1,
      dashed_edges='connectivities', solid_edges='transitions_confidence', transitions='transitions_confidence',
      title = basis,frameon = False, edges = True, save = False, show = False)
    else:
      sc.pl.paga_compare(adata, basis=basis,palette=palette, threshold = threshold, size = point_size, title = basis,frameon = False, edges = True, save = False, show = False)
    if show_plot is True:
      plt.show()
    if save:
      plt.savefig('.'.join(filter(None, [fileprefix, "paga_compare.png"])), dpi=dpi)
    
    sc.pl.paga(adata, threshold = threshold, layout = paga_layout, title = "PAGA layout: " + paga_layout, frameon = False, save = False, show = False)
    if show_plot is True:
      plt.show()
    if save:
      plt.savefig('.'.join(filter(None, [fileprefix, "paga_layout.png"])), dpi=dpi)
      
    sc.tl.draw_graph(adata, init_pos = "paga", layout = paga_layout)
    sc.pl.draw_graph(adata, color=group_by,palette=palette, title = "PAGA layout: " + paga_layout, layout = paga_layout,  frameon = False, legend_loc="on data",show = False)
    if show_plot is True:
      plt.show()
    if save:
      plt.savefig('.'.join(filter(None, [fileprefix, "paga_graph.png"])), dpi=dpi)
      
    if embedded_with_PAGA is True:
      umap2d = sc.tl.umap(adata, init_pos = "paga", n_components=2, copy = True)
      adata.obsm["PAGAUMAP2D"] = umap2d.obsm["X_umap"]
      sc.pl.paga_compare(adata, basis = "PAGAUMAP2D",palette=palette, threshold = threshold, size = point_size, title = "PAGA-initialized UMAP", edges = True, save = False, show = False)
      if show_plot is True:
        plt.show()
      if save:
        plt.savefig('.'.join(filter(None, [fileprefix, "paga_umap.png"])), dpi=dpi)
        
    if infer_pseudotime is True:
      if root_group is not None and root_cell is None:
        cell = adata.obs[group_by].index.values[adata.obs[group_by]==root_group]
        root_group_cell=adata.obsm[basis][adata.obs[group_by]==root_group,][:, [0, 1]]
        x=statistics.median(root_group_cell[:,0])
        y=statistics.median(root_group_cell[:,1])
        diff=np.array((x - root_group_cell[:,0], y - root_group_cell[:,1]))
        dist=[]
        for i in range(diff.shape[1]):
          dist.append(hypot(diff[0,i],diff[1,i]))

        root_cell=cell[dist.index(min(dist))]
        
      sc.tl.diffmap(adata, n_comps = n_dcs)
      adata.uns['iroot'] = np.flatnonzero(adata.obs_names == root_cell)[0]
      sc.tl.dpt(adata, n_dcs = n_dcs, n_branchings = n_branchings, min_group_size = min_group_size)
      
      sc.pl.embedding(adata, basis = basis, color = "dpt_pseudotime", save = False, show = False)
      if show_plot is True:
          plt.show()
      if save:
          plt.savefig('.'.join(filter(None, [fileprefix, "dpt_pseudotime.png"])), dpi=dpi)
      

  finally:
      os.chdir(prevdir)

  try:
    adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__[
        '_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
  except:
    pass

  return adata


def Palantir(adata=None, h5ad=None,group_by=None,palette=None,
             linear_reduction=None,nonlinear_reduction=None,basis=None,
             n_pcs=30,n_neighbors=30,dm_n_components=10,dm_alpha=0,dm_n_eigs=None,
             early_group = None, terminal_groups = None,early_cell=None,terminal_cells=None,
             num_waypoints=1200,scale_components=True,use_early_cell_as_start=False,
             adjust_early_cell=False,adjust_terminal_cells=False,
             max_iterations=25,n_jobs=1,
             point_size=20,
             show_plot=True, dpi=300, save=False, dirpath="./", fileprefix=""):
  import matplotlib.pyplot as plt
  import matplotlib
  import statistics
  from math import hypot
  import scanpy as sc
  import numpy as np
  import pandas as pd
  import palantir
  
  import warnings
  warnings.simplefilter("ignore", category=UserWarning)
  warnings.simplefilter("ignore", category=FutureWarning)
  warnings.simplefilter("ignore", category=DeprecationWarning)
  
  import os
  prevdir = os.getcwd()
  os.chdir(os.path.expanduser(dirpath))
  
  import platform
  if platform.system()=="Windows":
    import sys, multiprocessing, re
    if re.match(pattern=".*pythonw.exe$",string=sys.executable):
      pythonw=sys.executable
    else:
      pythonw=sys.executable.replace("python.exe", "pythonw.exe")
    sys.executable = pythonw
    sys._base_executable = pythonw
    multiprocessing.set_executable(pythonw) 

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
    
    if linear_reduction is None and nonlinear_reduction is None:
      print("linear_reduction or nonlinear_reduction must be provided at least one.")
      exit()
      
    if linear_reduction is None:
      sc.pp.pca(adata, n_comps = n_pcs)
      linear_reduction="X_pca"
      
    if basis is None:
      if nonlinear_reduction is not None:
        basis=nonlinear_reduction
      else:
        basis=linear_reduction
      
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
          
    terminal_cells_dict=dict()
    if terminal_groups is not None and terminal_cells is None:
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
        terminal_cells_dict[cell[dist.index(min(dist))]]=terminal_group.replace(" ", ".") + "_diff_potential"
        
      terminal_cells=list(terminal_cells_dict.keys())

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

    pca_projections = pd.DataFrame(adata.obsm[linear_reduction][:,:n_pcs], index=adata.obs_names)
    dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=dm_n_components, knn=n_neighbors, alpha=dm_alpha)
    ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=dm_n_eigs)
    pr_res = palantir.core.run_palantir(ms_data = ms_data,early_cell = early_cell, terminal_states = terminal_cells, knn = n_neighbors, num_waypoints = num_waypoints,
    scale_components = scale_components, use_early_cell_as_start = use_early_cell_as_start, max_iterations = max_iterations, n_jobs = n_jobs)
    
    if adjust_early_cell is True or adjust_terminal_cells is True:
      if adjust_early_cell is True:
        early_cell_group=adata.obs[group_by][early_cell]
        cells = adata.obs[group_by].index.values[adata.obs[group_by]==early_cell_group]
        early_cell = pr_res.pseudotime[cells].index.values[pr_res.pseudotime[cells]==min(pr_res.pseudotime[cells])][0]
      if adjust_terminal_cells is True:
        terminal_cells_dict=dict()
        for n in range(len(terminal_cells)):
          terminal_cell=terminal_cells[n]
          terminal_cell_group=adata.obs[group_by][terminal_cell]
          cells = adata.obs[group_by].index.values[adata.obs[group_by]==terminal_cell_group]
          terminal_cells_dict[pr_res.pseudotime[cells].index.values[pr_res.pseudotime[cells]==max(pr_res.pseudotime[cells])][0]]=terminal_cell_group.replace(" ", ".") + "_diff_potential"
        terminal_cells=list(terminal_cells_dict.keys())
        
      pr_res = palantir.core.run_palantir(ms_data = ms_data,early_cell = early_cell, terminal_states = terminal_cells, knn = n_neighbors, num_waypoints = num_waypoints,
    scale_components = scale_components, use_early_cell_as_start = use_early_cell_as_start, max_iterations = max_iterations, n_jobs = n_jobs)

    adata.obsm["palantir_dm"]=dm_res["T"].toarray()
    adata.uns["dm_kernel"]=dm_res["kernel"]
    if len(terminal_cells_dict)>0:
      pr_res.branch_probs=pr_res.branch_probs.rename(columns = terminal_cells_dict)
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

def WOT(adata=None, h5ad=None, group_by=None,palette=None, 
        time_field = "Time", growth_iters = 3, tmap_out = "tmaps/tmap_out",
        time_from = None, time_to = None, get_coupling = False,recalculate=False,
        show_plot=True, dpi=300, save=False, dirpath="./", fileprefix=""):
  import matplotlib.pyplot as plt
  import scanpy as sc
  import numpy as np
  import statistics
  import pandas as pd
  from math import hypot
  import wot

  import warnings
  warnings.simplefilter("ignore", category=UserWarning)
  warnings.simplefilter("ignore", category=FutureWarning)
  warnings.simplefilter("ignore", category=DeprecationWarning)
  
  import os
  prevdir = os.getcwd()
  os.chdir(os.path.expanduser(dirpath))
  
  import platform
  if platform.system()=="Windows":
    import sys, multiprocessing, re
    if re.match(pattern=".*pythonw.exe$",string=sys.executable):
      pythonw=sys.executable
    else:
      pythonw=sys.executable.replace("python.exe", "pythonw.exe")
    sys.executable = pythonw
    sys._base_executable = pythonw
    multiprocessing.set_executable(pythonw) 
  
  try:
    if adata is None and h5ad is None:
      print("adata or h5ad must be provided.")
      exit()
      
    if adata is None:
      adata = sc.read(h5ad)
      
    if group_by is None:
      print("group_by must be provided.")
      exit()
      
    if time_field is None:
      print("time_field must be provided.")
      exit()

    adata.obs[group_by] = adata.obs[group_by].astype(dtype = "category")
    if pd.api.types.is_categorical_dtype(adata.obs[time_field]):
      adata.obs["time_field"] = adata.obs[time_field].cat.codes
    elif not pd.api.types.is_numeric_dtype(adata.obs[time_field]):
      try:
        adata.obs['time_field'] = adata.obs[time_field].astype("float")
      except ValueError:
        print("Unable to convert column '" + time_field + "' to float type.")
    else:
      adata.obs["time_field"] = adata.obs[time_field]
      
    time_dict = dict(zip(adata.obs[time_field],adata.obs["time_field"]))
    if time_from not in time_dict.keys():
      print("'time_from' is incorrect")
      exit()
    
    ot_model = wot.ot.OTModel(adata, growth_iters = growth_iters, day_field = "time_field")
    
    if recalculate is True:
      ot_model.compute_all_transport_maps(tmap_out = tmap_out)
      tmap_model = wot.tmap.TransportMapModel.from_directory(tmap_out)
    else:
      try:
        tmap_model = wot.tmap.TransportMapModel.from_directory(tmap_out)
      except (FileNotFoundError, ValueError):
        ot_model.compute_all_transport_maps(tmap_out = tmap_out)
        tmap_model = wot.tmap.TransportMapModel.from_directory(tmap_out)

    cell_sets = {}
    for k, v in zip(adata.obs[group_by], adata.obs_names):
        if k not in cell_sets:
          cell_sets[k] = []
        cell_sets[k].append(v)
    
    from_populations = tmap_model.population_from_cell_sets(cell_sets, at_time = time_dict[time_from])
    
    trajectory_ds = tmap_model.trajectories(from_populations)
    trajectory_df = pd.DataFrame(trajectory_ds.X, index=trajectory_ds.obs_names, columns=trajectory_ds.var_names)
    adata.uns["trajectory_"+str(time_from)]= trajectory_df.reindex(adata.obs_names)

    fates_ds = tmap_model.fates(from_populations)
    fates_df = pd.DataFrame(fates_ds.X, index=fates_ds.obs_names, columns=fates_ds.var_names)
    existing_rows = fates_df.index.tolist()
    new_rows = list(set(adata.obs_names) - set(existing_rows))
    new_df = pd.DataFrame(0, index=new_rows, columns=fates_df.columns)
    fates_df = pd.concat([fates_df, new_df])
    adata.uns["fates_"+str(time_from)]= fates_df.reindex(adata.obs_names)
    
    # obs_list = wot.tmap.trajectory_trends_from_trajectory(trajectory_ds = trajectory_ds, expression_ds = adata)
  
    if time_to is not None:
      if time_to not in time_dict.keys():
        print("'time_to' is incorrect")
        exit()
      
      to_populations = tmap_model.population_from_cell_sets(cell_sets, at_time = time_dict[time_to])
      transition_table = tmap_model.transition_table(from_populations, to_populations)
      transition_df = pd.DataFrame(transition_table.X, index=transition_table.obs_names, columns=transition_table.var_names)
      adata.uns["transition_"+str(time_from)+"_to_"+str(time_to)]= fates_df
      if get_coupling is True:
        coupling = tmap_model.get_coupling(time_dict[time_from], time_dict[time_to])
        coupling_df = pd.DataFrame(coupling.X, index=coupling.obs_names, columns=coupling.var_names)
        adata.uns["coupling_"+str(time_from)+"_to_"+str(time_to)]= coupling_df

  finally:
      os.chdir(prevdir)

  try:
    adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__[
        '_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
  except:
    pass

  return adata
