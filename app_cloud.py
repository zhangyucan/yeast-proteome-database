import streamlit as st
import sqlite3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import base64
import os
from functools import lru_cache
from concurrent.futures import ThreadPoolExecutor
from utils import plot_cumulative_mass_fraction, plot_distribution, plot_scatter, plot_distribution_5, plot_log_scatter_5

# 设置页面标题和布局
st.set_page_config(page_title="Protein Mass Fraction Analysis", layout="wide")

executor = ThreadPoolExecutor(max_workers=2)

# 全局变量来存储文件内容
file_content = None

def prepare_file():
    """为云部署准备文件下载功能"""
    global file_content
    # 云部署时，可以使用示例数据或提供其他下载方式
    file_content = b"Demo file content - replace with actual data"

# 在后台异步准备文件
executor.submit(prepare_file)

@st.cache_data
@lru_cache(maxsize=5)
def get_rar_download_link():
    """提供数据下载链接"""
    try:
        # 云部署时提供其他下载方式
        href = '<p>数据下载功能在云版本中暂不可用。请联系 <a href="mailto:hongzhonglu@sjtu.edu.cn">hongzhonglu@sjtu.edu.cn</a> 获取完整数据集。</p>'
        return href
    except Exception as e:
        return f"Error: {str(e)}"

# 数据库连接 - 云部署版本
@st.cache_resource
def get_db_connection():
    """创建示例数据库连接"""
    # 为演示创建内存数据库
    conn = sqlite3.connect(":memory:")
    return conn

@st.cache_data
def create_sample_data():
    """创建示例数据"""
    # 创建示例的蛋白质质量分数数据
    genes = [f"YAL{str(i).zfill(3)}C" for i in range(1, 101)]
    conditions = [f"P{i}" for i in range(1, 11)]  # 简化为10个条件
    
    # 生成随机数据
    np.random.seed(42)
    data = {
        'gene': genes
    }
    for condition in conditions:
        data[condition] = np.random.exponential(0.001, len(genes))
    
    mass_fraction_df = pd.DataFrame(data)
    
    # 创建示例compartment数据
    compartments = ['cytoplasm', 'nucleus', 'mitochondria', 'ER', 'vacuole']
    compartment_data = {
        'gene': np.random.choice(genes, 200),
        'compartment': np.random.choice(compartments, 200)
    }
    compartment_df = pd.DataFrame(compartment_data)
    
    # 创建示例promass数据
    promass_data = {
        'compartment': compartments * 10,
    }
    for condition in conditions:
        promass_data[condition] = np.random.exponential(1, len(compartments) * 10)
    
    promass_df = pd.DataFrame(promass_data)
    
    # 创建mapping数据
    mapping_data = {
        'Condition_ID': conditions,
        'Description': [f"Experimental condition {i}" for i in range(1, 11)],
        'Growth_Media': np.random.choice(['YPD', 'SD', 'SC'], 10),
        'Temperature': np.random.choice([25, 30, 37], 10)
    }
    mapping_df = pd.DataFrame(mapping_data)
    
    return mass_fraction_df, compartment_df, promass_df, mapping_df

# 读取数据
mass_fraction_df, compartment_df, promass_df, mapping_df = create_sample_data()

# 顶部导航栏
st.markdown("""
<div style='background-color: #f0f2f6; padding: 1rem; border-radius: 3px; margin-bottom: 2rem;'>
    <div style='display: flex; justify-content: space-between; align-items: center;'>
        <h1 style='margin: 0; font-size: 2rem;'>A unified database of yeast absolute quantitative proteomics</h1>
    </div>
</div>
""", unsafe_allow_html=True)

# 云部署提示
st.markdown("""
<div style='background-color: #fff3cd; padding: 1rem; border-radius: 3px; margin-bottom: 1rem; border-left: 4px solid #ffc107;'>
    <strong>🌟 演示版本说明:</strong> 这是云部署的演示版本，使用模拟数据展示功能。完整数据集请联系研究团队获取。
</div>
""", unsafe_allow_html=True)

table_choice = st.sidebar.radio(
    "Please select the required function", 
    ["Main Page", "Compute", "Download", "About us"]
)

if table_choice == "Download":
    st.markdown("## Download Data")
    st.markdown("Download the database in RAR format.")
    st.markdown(get_rar_download_link(), unsafe_allow_html=True)

if table_choice == "About us":
    st.markdown("## About Us")
    st.markdown("This database contains absolute quantitative proteomic data from <i>Saccharomyces cerevisiae</i> under a variety of experimental settings. These datasets are valuable resources for yeast physiology, synthetic biology, and systems biology research. We will continue to update this database when fresh experimental datasets become available. If you have any question or suggestion, please contact with Hongzhong Lu ([hongzhonglu@sjtu.edu.cn](hongzhonglu@sjtu.edu.cn), [https://life.sjtu.edu.cn/teacher/En/luhongzhong](https://life.sjtu.edu.cn/teacher/En/luhongzhong))", unsafe_allow_html=True)

if table_choice == "Main Page":
    search_query = st.text_input("Search proteins...")

    if search_query:
        # 使用模糊匹配进行搜索
        filtered_df = mass_fraction_df[mass_fraction_df['gene'].str.contains(search_query, case=False, na=False)]
        if not filtered_df.empty:
            st.subheader(f"Search Results for '{search_query}'")
            st.dataframe(filtered_df)
            
            # 显示搜索结果的数量
            st.write(f"Found {len(filtered_df)} matching proteins")
            
            # 添加可视化搜索结果的选项
            if st.checkbox("Show visualization of search results"):
                # 选择要显示的条件
                columns_to_plot = st.multiselect(
                    "Select conditions to visualize",
                    [col for col in filtered_df.columns if col.startswith('P')],
                    default=[filtered_df.columns[1]]  # 默认选择第一个P列
                )
                
                if columns_to_plot:
                    sorted_columns = sorted(columns_to_plot, key=lambda col: filtered_df[col].sum(), reverse=True)
                    fig, ax = plt.subplots(figsize=(10, 6))
                    for col in sorted_columns:
                        plt.bar(col, filtered_df[col].sum(), label=col, alpha=0.7)
                    plt.xticks(rotation=45, ha='right')
                    plt.legend(title="Conditions")
                    plt.title(f"Mass Fraction Distribution for '{search_query}'")
                    plt.xlabel("Condition")
                    plt.ylabel("Mass Fraction")
                    st.pyplot(fig)
        else:
            st.warning(f"No proteins found matching '{search_query}'")
    else:
        # 当没有搜索时显示概览数据
        st.subheader("Mass Fraction Data Overview (Demo Data)")
        st.dataframe(mass_fraction_df.head(100))

        st.subheader("Mapping of P1-10 Overview (Demo Data)")
        st.dataframe(mapping_df)

if table_choice == "Compute":
    module = st.sidebar.selectbox(
        "Select Module",
        ["Compartment Analysis", "Compartment Mass Ratio", "Protein Mass Distribution"]
    )
    
    if module == "Compartment Analysis":
        st.subheader("Compartment Analysis")
        compartment = st.selectbox("Select Compartment", compartment_df['compartment'].unique())
        cond = st.selectbox("Select Condition", [f'P{i}' for i in range(1, 11)], key="cond1")
        
        if st.button("Generate Analysis"):
            try:
                plot_cumulative_mass_fraction(compartment, cond, mass_fraction_df, compartment_df)
            except Exception as e:
                st.error(f"分析功能在演示版本中可能不完全可用: {str(e)}")
    
    elif module == "Compartment Mass Ratio":
        st.subheader("Compartment Mass Ratio Analysis")
        analysis_type = st.radio("Select Analysis Type", ["Single Condition", "Two Conditions"])
        
        if analysis_type == "Single Condition":
            column = st.selectbox("Select Condition", [f'P{i}' for i in range(1, 11)])
            if st.button("Generate Plot"):
                try:
                    plot_distribution(promass_df, column)
                except Exception as e:
                    st.error(f"绘图功能在演示版本中可能不完全可用: {str(e)}")
        else:
            col1, col2 = st.columns(2)
            with col1:
                column1 = st.selectbox("Select First Condition", [f'P{i}' for i in range(1, 11)])
            with col2:
                column2 = st.selectbox("Select Second Condition", [f'P{i}' for i in range(1, 11)])
            if st.button("Generate Plot"):
                try:
                    plot_scatter(promass_df, column1, column2)
                except Exception as e:
                    st.error(f"绘图功能在演示版本中可能不完全可用: {str(e)}")
    
    elif module == "Protein Mass Distribution":
        st.subheader("Protein Mass Distribution Analysis")
        analysis_type = st.radio("Select Analysis Type", ["Single Condition", "Two Conditions"])
        
        if analysis_type == "Single Condition":
            column = st.selectbox("Select Condition", [f'P{i}' for i in range(1, 11)])
            if st.button("Generate Plot"):
                try:
                    plot_distribution_5(mass_fraction_df, column)
                except Exception as e:
                    st.error(f"绘图功能在演示版本中可能不完全可用: {str(e)}")
        else:
            col1, col2 = st.columns(2)
            with col1:
                column1 = st.selectbox("Select First Condition", [f'P{i}' for i in range(1, 11)])
            with col2:
                column2 = st.selectbox("Select Second Condition", [f'P{i}' for i in range(1, 11)])
            if st.button("Generate Plot"):
                try:
                    plot_log_scatter_5(mass_fraction_df, column1, column2)
                except Exception as e:
                    st.error(f"绘图功能在演示版本中可能不完全可用: {str(e)}")
