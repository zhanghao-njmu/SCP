# 测试覆盖率增强计划

**当前覆盖率**: 40-50%
**目标覆盖率**: 70-80%
**策略**: 优先测试新拆分模块和可快速验证的函数

---

## 当前状态分析

### 已有测试文件 (20个)
- ✅ test-basic-validation.R
- ✅ test-seurat-v5-compatibility.R
- ✅ test-ggplot2-compatibility.R
- ✅ test-utils.R
- ✅ test-utility-helpers.R
- ✅ test-color-functions.R
- ✅ test-data.R
- ✅ test-seurat-integration.R
- ✅ test-qc-functions.R
- ✅ test-workflow-functions.R
- ✅ test-dimensionality-reduction.R
- ✅ test-analysis-functions.R
- ✅ test-integration-functions.R
- ✅ test-annotation-functions.R
- ✅ test-plotting-helpers.R
- ✅ test-plotting-core.R
- ✅ test-sankey-alluvial.R
- ✅ test-gene-conversion.R
- ✅ test-trajectory-functions.R
- ✅ test-python-integration.R

---

## 需要增加测试的模块

### 🎨 绘图模块 (新拆分)

#### SCP-plot-themes.R
- ✅ theme_scp - 已测试
- ✅ theme_blank - 已测试
- ✅ palette_scp - 已测试
- ✅ show_palettes - 已测试
- 🟡 panel_fix - 需要更多测试
- 🟡 panel_fix_overall - 需要更多测试

#### SCP-plot-helpers.R
- ✅ drop_data - 已测试
- ✅ slim_data - 已测试
- ✅ get_vars - 已测试
- ✅ adjcolors - 已测试
- ✅ blendcolors - 已测试

#### SCP-plot-dimred.R
- 🟡 CellDimPlot - 需要更多edge cases
- 🟡 FeatureDimPlot - 需要更多edge cases
- ⚠️ CellDimPlot3D - 最小测试（3D依赖）
- ⚠️ FeatureDimPlot3D - 最小测试（3D依赖）

#### SCP-plot-stats.R
- 🟡 FeatureStatPlot - 需要更多测试
- 🟡 CellStatPlot - 需要更多测试
- 🟡 StatPlot - 需要更多测试
- 🟡 FeatureCorPlot - 需要更多测试
- 🟡 CellDensityPlot - 已基本测试

#### SCP-plot-volcano.R
- 🟡 VolcanoPlot - 需要更多测试

#### SCP-plot-heatmaps.R
- 🟠 GroupHeatmap - 基础测试
- 🟠 FeatureHeatmap - 基础测试
- 🟠 CellCorHeatmap - 基础测试

#### SCP-plot-dynamic.R
- ⚠️ DynamicHeatmap - 最小测试（交互式）
- ⚠️ DynamicPlot - 最小测试（交互式）
- 🟠 ProjectionPlot - 需要测试

#### SCP-plot-enrichment.R
- 🟡 EnrichmentPlot - 基础测试
- 🟡 GSEAPlot - 基础测试

### 🔬 分析模块 (新拆分)

#### SCP-analysis-gene.R
- ⚠️ GeneConvert - 最小测试（需要网络）
- ✅ CC_GenePrefetch - 已测试

#### SCP-analysis-scoring.R
- 🟡 CellScoring - 需要更多测试

#### SCP-analysis-de.R
- 🟡 FindExpressedMarkers - 基础测试
- 🟡 RunDEtest - 基础测试

#### SCP-analysis-database.R
- ⚠️ ListDB - 最小测试（网络依赖）
- ⚠️ PrepareDB - 最小测试（网络依赖，1523行）

#### SCP-analysis-enrichment.R
- 🟡 RunEnrichment - 基础测试
- 🟡 RunGSEA - 基础测试

#### SCP-analysis-trajectory.R
- 🟡 RunSlingshot - 基础测试
- 🟡 RunMonocle2 - 基础测试
- 🟡 RunMonocle3 - 基础测试
- 🟠 RunDynamicFeatures - 需要测试
- 🟠 RunDynamicEnrichment - 需要测试

#### SCP-analysis-python.R
- 🟡 srt_to_adata - 基础测试
- 🟡 adata_to_srt - 基础测试
- ⚠️ RunPAGA - 最小测试（Python依赖）
- ⚠️ RunSCVELO - 最小测试（Python依赖）
- ⚠️ RunPalantir - 最小测试（Python依赖）
- ⚠️ RunWOT - 最小测试（Python依赖）

### 🔄 工作流模块 (新拆分)

#### SCP-workflow-utils.R
- ✅ check_DataType - 已测试
- ✅ check_srtList - 已测试
- ✅ check_srtMerge - 已测试
- 🟡 RecoverCounts - 需要测试
- 🟡 RenameFeatures - 需要测试
- 🟡 RenameClusters - 需要测试
- ✅ SrtReorder - 已测试
- ✅ SrtAppend - 已测试

#### SCP-workflow-reduction.R
- 🟡 RunDimReduction - 基础测试
- ✅ DefaultReduction - 已测试

#### SCP-workflow-integration.R
- ✅ Uncorrected_integrate - 已测试
- 🟠 Seurat_integrate - 基础框架测试
- 🟠 scVI_integrate - 基础框架测试
- 🟠 MNN_integrate - 基础框架测试
- 🟠 fastMNN_integrate - 基础框架测试
- 🟠 Harmony_integrate - 基础框架测试
- 🟠 Scanorama_integrate - 基础框架测试
- 🟠 BBKNN_integrate - 基础框架测试
- 🟠 CSS_integrate - 基础框架测试
- 🟠 LIGER_integrate - 基础框架测试
- 🟠 Conos_integrate - 基础框架测试
- 🟠 ComBat_integrate - 基础框架测试

#### SCP-workflow-pipelines.R
- ⚠️ Standard_SCP - 最小测试（计算密集）
- ⚠️ Integration_SCP - 最小测试（计算密集）

---

## 优先增加的测试

### Phase 1: 高优先级（快速测试，无外部依赖）

1. **工作流工具函数**
   - RecoverCounts
   - RenameFeatures
   - RenameClusters

2. **统计图表增强**
   - FeatureStatPlot edge cases
   - CellStatPlot edge cases
   - StatPlot edge cases
   - FeatureCorPlot edge cases

3. **火山图**
   - VolcanoPlot 更多参数组合

4. **投影图**
   - ProjectionPlot

5. **面板函数**
   - panel_fix 更多场景
   - panel_fix_overall 更多场景

### Phase 2: 中优先级（需要数据但可测试）

1. **热图增强**
   - GroupHeatmap 更多参数
   - FeatureHeatmap 更多参数
   - CellCorHeatmap 更多参数

2. **降维图增强**
   - CellDimPlot 更多edge cases
   - FeatureDimPlot 更多edge cases

3. **分析函数增强**
   - CellScoring 更多基因集
   - FindExpressedMarkers 更多参数
   - RunDEtest 更多方法

4. **轨迹分析增强**
   - RunDynamicFeatures
   - RunDynamicEnrichment

### Phase 3: 低优先级（跳过或最小测试）

1. **网络依赖**
   - GeneConvert (skip)
   - PrepareDB (skip)
   - ListDB (skip)

2. **Python依赖**
   - RunPAGA (skip)
   - RunSCVELO (skip)
   - RunPalantir (skip)
   - RunWOT (skip)

3. **计算密集**
   - Standard_SCP (skip)
   - Integration_SCP (skip)
   - 大部分整合方法 (minimal)

4. **交互式/3D**
   - DynamicHeatmap (skip)
   - DynamicPlot (skip)
   - CellDimPlot3D (minimal)
   - FeatureDimPlot3D (minimal)

---

## 预期覆盖率提升

| 阶段 | 新增测试 | 预期覆盖率 |
|------|----------|------------|
| 当前 | 191 | 40-50% |
| Phase 1 | +50 | 55-65% |
| Phase 2 | +40 | 70-75% |
| Phase 3 | skip/minimal | 70-80% |

---

## 实施步骤

1. ✅ 分析当前覆盖
2. ⏳ Phase 1: 创建工作流工具函数测试
3. ⏳ Phase 1: 增强统计图表测试
4. ⏳ Phase 1: 增加火山图和投影图测试
5. ⏳ Phase 2: 增强热图测试
6. ⏳ Phase 2: 增强分析函数测试
7. ⏳ 验证所有新测试
8. ⏳ 提交代码

---

**目标**: 从40-50%提升到70-80%的实际可测试代码覆盖率
**方法**: 聚焦可快速验证的函数，跳过外部依赖
**时间**: 估计增加90个新测试用例
