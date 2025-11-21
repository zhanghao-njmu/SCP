# SCP-analysis.R 拆分方案

## 原文件统计
- **文件名**: SCP-analysis.R
- **大小**: 289K
- **行数**: 6,228行
- **导出函数**: 20个

## 拆分策略

将SCP-analysis.R拆分为7个功能模块文件，按分析功能逻辑组织。

---

## 拆分文件列表

### 1. SCP-analysis-gene.R (661行)
**功能**: 基因ID转换和细胞周期基因

**包含函数**:
- `GeneConvert` (行70-471) - 基因ID类型转换（跨物种）
- `CC_GenePrefetch` (行472-730) - 细胞周期基因预提取

**行范围**: 1-730

---

### 2. SCP-analysis-scoring.R (498行)
**功能**: 细胞评分和基因集评分

**包含函数**:
- `CellScoring` (行731-1228) - 基因集细胞评分

**行范围**: 731-1228

---

### 3. SCP-analysis-de.R (894行)
**功能**: 差异表达分析

**包含函数**:
- `FindExpressedMarkers` (行1229-1670) - 查找表达标记基因
- `RunDEtest` (行1671-2122) - 运行差异表达测试

**行范围**: 1229-2122

---

### 4. SCP-analysis-database.R (1643行)
**功能**: 数据库准备和管理

**包含函数**:
- `ListDB` (行2123-2242) - 列出可用数据库
- `PrepareDB` (行2243-3765) - 准备注释数据库（超大函数 1,523行）

**行范围**: 2123-3765

---

### 5. SCP-analysis-enrichment.R (541行)
**功能**: 富集分析和GSEA

**包含函数**:
- `RunEnrichment` (行3766-4029) - 运行富集分析
- `RunGSEA` (行4030-4306) - 运行GSEA分析

**行范围**: 3766-4306

---

### 6. SCP-analysis-trajectory.R (911行)
**功能**: 轨迹推断和动态分析

**包含函数**:
- `RunSlingshot` (行4307-4409) - Slingshot轨迹分析
- `RunMonocle2` (行4410-4744) - Monocle2轨迹分析
- `RunMonocle3` (行4745-4947) - Monocle3轨迹分析
- `RunDynamicFeatures` (行4948-5217) - 动态特征分析
- `RunDynamicEnrichment` (行5218-5350) - 动态富集分析

**行范围**: 4307-5350

---

### 7. SCP-analysis-python.R (878行)
**功能**: Python工具集成

**包含函数**:
- `srt_to_adata` (行5351-5520) - Seurat转AnnData
- `adata_to_srt` (行5521-5784) - AnnData转Seurat
- `RunPAGA` (行5785-5894) - PAGA分析
- `RunSCVELO` (行5895-6000) - RNA速度分析
- `RunPalantir` (行6001-6085) - Palantir轨迹分析
- `RunWOT` (行6086-6228) - Waddington OT分析

**行范围**: 5351-6228

---

## 拆分后统计

| 文件名 | 行数 | 函数数 | 功能描述 |
|--------|------|--------|----------|
| SCP-analysis-gene.R | 730 | 2 | 基因转换 |
| SCP-analysis-scoring.R | 498 | 1 | 细胞评分 |
| SCP-analysis-de.R | 894 | 2 | 差异表达 |
| SCP-analysis-database.R | 1643 | 2 | 数据库管理 |
| SCP-analysis-enrichment.R | 541 | 2 | 富集分析 |
| SCP-analysis-trajectory.R | 1044 | 5 | 轨迹推断 |
| SCP-analysis-python.R | 878 | 6 | Python集成 |
| **总计** | **6228** | **20** | **7个模块** |

---

## 功能分组逻辑

### 基因层面
- **gene.R**: 基因ID转换、细胞周期基因

### 分析层面
- **scoring.R**: 基因集评分
- **de.R**: 差异表达分析
- **enrichment.R**: 功能富集

### 数据层面
- **database.R**: 数据库和注释资源

### 轨迹层面
- **trajectory.R**: 伪时序和轨迹推断（R包）
- **python.R**: Python工具集成（PAGA, scVelo等）

---

## 实施步骤

1. ✅ 分析原文件结构
2. ✅ 识别所有函数及其行范围
3. ✅ 按功能逻辑分组
4. ⏳ 提取并创建7个新文件
5. ⏳ 验证所有函数正确迁移
6. ⏳ 删除原SCP-analysis.R文件
7. ⏳ 更新NAMESPACE (roxygen2自动处理)
8. ⏳ 运行测试确保无破坏
9. ⏳ 提交更改

---

## 注意事项

### 大函数警告 ⚠️
- **PrepareDB** (1,523行) - 这个函数非常大，但由于是单一功能（数据库准备），保持在database.R中

### 依赖关系
- gene.R 和 scoring.R 可能被其他模块依赖
- python.R 依赖 reticulate 包
- trajectory.R 包含R实现的轨迹分析

### 测试考虑
- database.R 的函数需要网络访问（可能被跳过）
- python.R 的函数需要Python环境（应该跳过）
- trajectory.R 的某些函数计算密集（可能被跳过）

---

## 优势

- ✅ 功能模块化清晰
- ✅ 文件大小合理 (500-1700行)
- ✅ 易于定位特定分析功能
- ✅ 减少编辑冲突
- ✅ Python相关功能独立模块
- ✅ 数据库管理独立维护
