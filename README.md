# 酵母绝对定量蛋白质组学统一数据库

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://yeast-proteome-database.streamlit.app)
[![GitHub Deploy](https://github.com/your-username/yeast-proteome-database/actions/workflows/deploy.yml/badge.svg)](https://github.com/your-username/yeast-proteome-database/actions/workflows/deploy.yml)

## 项目简介

这是一个酵母（*Saccharomyces cerevisiae*）绝对定量蛋白质组学数据的统一数据库，包含多种实验条件下的数据集。这些数据集为酵母生理学、合成生物学和系统生物学研究提供了宝贵的资源。

## 功能特性

- 🔍 **蛋白质搜索**: 根据基因名称快速搜索蛋白质数据
- 📊 **数据可视化**: 提供多种图表展示蛋白质质量分数分布
- 🧬 **细胞器分析**: 分析不同细胞器中的蛋白质分布
- 📈 **比较分析**: 支持多条件间的比较分析
- 💾 **数据下载**: 提供完整数据集下载功能

## 在线访问

- **云端演示版**: [https://yeast-proteome-database.streamlit.app](https://yeast-proteome-database.streamlit.app)
- **GitHub仓库**: [https://github.com/your-username/yeast-proteome-database](https://github.com/your-username/yeast-proteome-database)
- **完整版本**: 请联系研究团队获取完整数据集

## 快速部署到GitHub和Streamlit Cloud

### 1. 部署到GitHub

```bash
# 1. 在GitHub上创建新仓库 yeast-proteome-database
# 2. 克隆到本地并添加文件
git clone https://github.com/your-username/yeast-proteome-database.git
cd yeast-proteome-database

# 3. 复制项目文件到仓库目录
# 4. 提交并推送
git add .
git commit -m "Initial commit: 酵母蛋白质组学数据库"
git push origin main
```

### 2. 部署到Streamlit Cloud

1. 访问 [https://streamlit.io/cloud](https://streamlit.io/cloud)
2. 使用GitHub账号登录
3. 点击 "New app"
4. 选择你的GitHub仓库 `yeast-proteome-database`
5. 设置入口文件为 `app_cloud.py`
6. 点击 "Deploy"

### 3. 自定义域名（可选）

在Streamlit Cloud应用设置中，你可以自定义域名为更友好的名称。

## 本地运行

### 环境要求

- Python 3.8+
- 依赖包见 `requirements.txt`

### 安装步骤

1. 克隆仓库:
```bash
git clone https://github.com/your-username/yeast-proteome-database.git
cd yeast-proteome-database
```

2. 安装依赖:
```bash
pip install -r requirements.txt
```

3. 运行应用:
```bash
# 云端演示版本
streamlit run app_cloud.py

# 完整版本（需要数据文件）
streamlit run new_search.py
```

## 数据说明

### 数据表结构

1. **mass_fraction_combine**: 蛋白质质量分数数据
2. **compartment_annotation_refine**: 细胞器注释数据
3. **ProMassRatio_across_compartment_combine**: 跨细胞器蛋白质质量比例数据

### 实验条件

数据库包含275个不同的实验条件（P1-P275），涵盖：
- 不同培养基条件
- 温度变化
- 压力条件
- 基因扰动

## 使用指南

### 主页搜索
在搜索框中输入基因名称（如 `YAL001C`），系统会返回匹配的蛋白质数据。

### 计算模块
1. **细胞器分析**: 分析特定细胞器中蛋白质的累积质量分数
2. **细胞器质量比例**: 比较不同条件下细胞器间的蛋白质分布
3. **蛋白质质量分布**: 展示蛋白质质量分数的整体分布

## 技术栈

- **前端**: Streamlit
- **数据处理**: Pandas, NumPy
- **可视化**: Matplotlib, Seaborn, Plotly
- **数据库**: SQLite

## 联系我们

如有任何问题或建议，请联系：

**陆宏钟教授**
- 邮箱: hongzhonglu@sjtu.edu.cn
- 主页: https://life.sjtu.edu.cn/teacher/En/luhongzhong

## 引用

如果您在研究中使用了本数据库，请引用：

```
[待发表论文引用格式]
```

## 许可证

本项目采用 MIT 许可证 - 详见 [LICENSE](LICENSE) 文件。

## 更新日志

- v1.0.0: 初始版本发布
- 持续更新中...

---

*本数据库将在新的实验数据集可用时持续更新。*
