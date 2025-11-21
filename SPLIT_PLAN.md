# SCP-plot.R 拆分方案

## 原文件统计
- **文件名**: SCP-plot.R
- **大小**: 714K
- **行数**: 14,794行
- **导出函数**: 41个

## 拆分策略

将SCP-plot.R拆分为9个功能模块文件，每个文件500-3100行，按功能逻辑组织。

---

## 拆分文件列表

### 1. SCP-plot-themes.R (832行)
**功能**: 主题、调色板、面板控制

**包含函数**:
- `theme_scp` (行16-83) - SCP默认主题
- `theme_blank` (行84-167) - 空白主题
- `palette_scp` (行168-297) - 调色板生成
- `show_palettes` (行298-452) - 展示调色板
- `panel_fix` (行453-587) - 面板修复
- `panel_fix_overall` (行588-832) - 整体面板修复

**行范围**: 1-832

---

### 2. SCP-plot-helpers.R (495行)
**功能**: 数据处理和颜色辅助函数

**包含函数**:
- `drop_data` (行833-840) - 数据丢弃通用方法
- `drop_data.ggplot` (行841-887) - ggplot方法
- `drop_data.patchwork` (行888-898) - patchwork方法
- `drop_data.default` (行899-917) - 默认方法
- `slim_data` (行918-925) - 数据精简通用方法
- `slim_data.ggplot` (行926-942) - ggplot方法
- `slim_data.patchwork` (行943-952) - patchwork方法
- `slim_data.default` (行953-962) - 默认方法
- `get_vars` (行963-999) - 获取变量
- `adjcolors` (行1000-1025) - 调整颜色
- `blendcolors` (行1026-1327) - 混合颜色

**行范围**: 833-1327

---

### 3. SCP-plot-dimred.R (2196行)
**功能**: 维度降维可视化 (2D和3D)

**包含函数**:
- `CellDimPlot` (行1328-2041) - 细胞维度图
- `FeatureDimPlot` (行2042-2842) - 特征维度图
- `CellDimPlot3D` (行2843-3065) - 3D细胞维度图
- `FeatureDimPlot3D` (行3066-3523) - 3D特征维度图

**行范围**: 1328-3523

---

### 4. SCP-plot-stats.R (2394行)
**功能**: 统计图表、相关图、密度图

**包含函数**:
- `FeatureStatPlot` (行3524-4462) - 特征统计图
- `CellStatPlot` (行4463-4563) - 细胞统计图
- `StatPlot` (行4564-5231) - 通用统计图
- `FeatureCorPlot` (行5232-5678) - 特征相关图
- `CellDensityPlot` (行5679-5917) - 细胞密度图

**行范围**: 3524-5917

---

### 5. SCP-plot-trajectory.R (1170行)
**功能**: 轨迹分析、网络图、速度场

**包含函数**:
- `LineagePlot` (行5918-6120) - 谱系图
- `PAGAPlot` (行6121-6286) - PAGA图
- `GraphPlot` (行6287-6665) - 图网络可视化
- `segementsDf` (行6666-6744) - 线段数据框
- `VelocityPlot` (行6745-6962) - 速度图
- `compute_velocity_on_grid` (行6963-7087) - 网格速度计算

**行范围**: 5918-7087

---

### 6. SCP-plot-volcano.R (823行)
**功能**: 火山图可视化

**包含函数**:
- `VolcanoPlot` (行7088-7910) - 火山图

**行范围**: 7088-7910

---

### 7. SCP-plot-heatmaps.R (3109行)
**功能**: 各类热图可视化

**包含函数**:
- `GroupHeatmap` (行7911-9113) - 分组热图
- `FeatureHeatmap` (行9114-10112) - 特征热图
- `CellCorHeatmap` (行10113-11019) - 细胞相关热图

**行范围**: 7911-11019

---

### 8. SCP-plot-dynamic.R (1707行)
**功能**: 动态图和投影可视化

**包含函数**:
- `DynamicHeatmap` (行11020-12101) - 动态热图
- `DynamicPlot` (行12102-12510) - 动态图
- `ProjectionPlot` (行12511-12726) - 投影图

**行范围**: 11020-12726

---

### 9. SCP-plot-enrichment.R (2068行)
**功能**: 富集分析可视化

**包含函数**:
- `EnrichmentPlot` (行12727-13562) - 富集分析图
- `GSEAPlot` (行13563-14794) - GSEA图

**行范围**: 12727-14794

---

## 拆分后统计

| 文件名 | 行数 | 函数数 | 功能描述 |
|--------|------|--------|----------|
| SCP-plot-themes.R | 832 | 6 | 主题和样式 |
| SCP-plot-helpers.R | 495 | 11 | 数据处理辅助 |
| SCP-plot-dimred.R | 2196 | 4 | 维度降维可视化 |
| SCP-plot-stats.R | 2394 | 5 | 统计和相关图 |
| SCP-plot-trajectory.R | 1170 | 6 | 轨迹和网络 |
| SCP-plot-volcano.R | 823 | 1 | 火山图 |
| SCP-plot-heatmaps.R | 3109 | 3 | 热图 |
| SCP-plot-dynamic.R | 1707 | 3 | 动态图 |
| SCP-plot-enrichment.R | 2068 | 2 | 富集分析 |
| **总计** | **14794** | **41** | **9个模块** |

---

## 实施步骤

1. ✅ 分析原文件结构
2. ✅ 识别所有函数及其行范围
3. ✅ 按功能逻辑分组
4. ⏳ 提取并创建9个新文件
5. ⏳ 验证所有函数正确迁移
6. ⏳ 删除原SCP-plot.R文件
7. ⏳ 更新NAMESPACE (roxygen2自动处理)
8. ⏳ 运行测试确保无破坏
9. ⏳ 提交更改

---

## 注意事项

1. **NAMESPACE自动生成**: 所有@export标记会被roxygen2自动处理
2. **内部函数**: 未导出的辅助函数也需要保留在相应模块
3. **依赖关系**: 确保各模块间的函数调用不会破坏
4. **测试兼容**: 现有测试应该继续正常工作
5. **文档完整**: 所有roxygen2注释必须保留

---

## 优势

- ✅ 文件大小合理 (500-3100行)
- ✅ 功能逻辑清晰
- ✅ 易于维护和扩展
- ✅ 便于代码审查
- ✅ 提高可读性
- ✅ 减少merge冲突
