# SCP-workflow.R 拆分方案

## 原文件统计
- **文件名**: SCP-workflow.R
- **大小**: 174K
- **行数**: 3,959行
- **导出函数**: 24个

## 拆分策略

将SCP-workflow.R拆分为4个功能模块文件，按工作流阶段组织。

---

## 拆分文件列表

### 1. SCP-workflow-utils.R (780行)
**功能**: 数据验证、操作和重命名工具

**包含函数**:
- `check_DataType` (行14-74) - 检查数据类型
- `check_srtList` (行75-333) - 检查Seurat对象列表
- `check_srtMerge` (行334-419) - 检查合并参数
- `RecoverCounts` (行420-499) - 恢复原始计数
- `RenameFeatures` (行500-559) - 重命名特征
- `RenameClusters` (行560-622) - 重命名聚类
- `SrtReorder` (行623-691) - 重排Seurat对象
- `SrtAppend` (行692-780) - 追加数据到Seurat对象

**行范围**: 1-780

---

### 2. SCP-workflow-reduction.R (311行)
**功能**: 降维分析

**包含函数**:
- `RunDimReduction` (行781-1033) - 运行降维分析（通用接口）
- `DefaultReduction` (行1034-1091) - 设置默认降维方法

**行范围**: 781-1091

---

### 3. SCP-workflow-integration.R (2571行)
**功能**: 批次效应校正和数据整合

**包含函数**:
- `Uncorrected_integrate` (行1092-1275) - 无校正整合
- `Seurat_integrate` (行1276-1572) - Seurat CCA整合
- `scVI_integrate` (行1573-1761) - scVI整合
- `MNN_integrate` (行1762-1967) - MNN整合
- `fastMNN_integrate` (行1968-2139) - FastMNN整合
- `Harmony_integrate` (行2140-2347) - Harmony整合
- `Scanorama_integrate` (行2348-2540) - Scanorama整合
- `BBKNN_integrate` (行2541-2764) - BBKNN整合
- `CSS_integrate` (行2765-2970) - CSS整合
- `LIGER_integrate` (行2971-3182) - LIGER整合
- `Conos_integrate` (行3183-3394) - Conos整合
- `ComBat_integrate` (行3395-3662) - ComBat整合

**行范围**: 1092-3662

---

### 4. SCP-workflow-pipelines.R (297行)
**功能**: 端到端分析管道

**包含函数**:
- `Standard_SCP` (行3663-3906) - 标准单细胞分析流程
- `Integration_SCP` (行3907-3959) - 整合分析流程

**行范围**: 3663-3959

---

## 拆分后统计

| 文件名 | 行数 | 函数数 | 功能描述 |
|--------|------|--------|----------|
| SCP-workflow-utils.R | 780 | 8 | 数据操作工具 |
| SCP-workflow-reduction.R | 311 | 2 | 降维分析 |
| SCP-workflow-integration.R | 2571 | 12 | 批次整合方法 |
| SCP-workflow-pipelines.R | 297 | 2 | 完整工作流 |
| **总计** | **3959** | **24** | **4个模块** |

---

## 功能分组逻辑

### 前处理层 (utils.R)
数据验证、质量控制、命名管理
- 检查数据类型和格式
- 验证对象列表
- 重命名和重排序
- 数据追加和合并

### 分析层 (reduction.R)
降维和特征提取
- PCA, UMAP, tSNE等降维方法
- 默认降维设置

### 整合层 (integration.R)
批次效应校正（12种方法）
- 基于锚点: Seurat CCA, MNN, fastMNN
- 基于嵌入: Harmony, Scanorama, BBKNN
- 基于因子: LIGER
- 基于图: Conos
- 基于回归: ComBat
- 深度学习: scVI
- 其他: CSS

### 管道层 (pipelines.R)
端到端自动化流程
- 标准分析管道
- 整合分析管道

---

## 特殊考虑

### 大文件模块
- **integration.R** (2,571行) - 包含12个整合方法
  - 这些方法相互独立但功能相似
  - 保持在一个文件中便于比较和选择
  - 每个方法平均约200行

### 依赖关系
- utils.R → 被所有其他模块依赖
- reduction.R → 被pipelines.R依赖
- integration.R → 被pipelines.R依赖
- pipelines.R → 使用上述所有模块

### 测试策略
- utils.R: 可以完整测试（数据操作）
- reduction.R: 可以测试（相对快速）
- integration.R: 部分跳过（计算密集）
- pipelines.R: 跳过（完整流程太慢）

---

## 实施步骤

1. ✅ 分析原文件结构
2. ✅ 识别所有函数及其行范围
3. ✅ 按功能逻辑分组
4. ⏳ 提取并创建4个新文件
5. ⏳ 验证所有函数正确迁移
6. ⏳ 删除原SCP-workflow.R文件
7. ⏳ 更新NAMESPACE (roxygen2自动处理)
8. ⏳ 运行测试确保无破坏
9. ⏳ 提交更改

---

## 优势

- ✅ 清晰的工作流阶段划分
- ✅ 整合方法集中管理
- ✅ 管道代码独立模块
- ✅ 工具函数易于重用
- ✅ 减少文件编辑冲突
- ✅ 便于添加新整合方法
