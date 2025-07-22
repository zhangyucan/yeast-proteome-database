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
    global file_content
    # 使用相对路径，适用于云部署
    file_path = "protein_database.rar"
    try:
        if os.path.exists(file_path):
            with open(file_path, 'rb') as f:
                file_content = f.read()
            print(f"Successfully loaded {len(file_content)} bytes from {file_path}")
        else:
            print(f"File {file_path} not found")
            file_content = None
    except Exception as e:
        print(f"Error loading file {file_path}: {str(e)}")
        file_content = None

# 在后台异步准备文件
executor.submit(prepare_file)

@st.cache_data
@lru_cache(maxsize=5)
def get_rar_download_link(rar_path):
    try:
        # 先检查文件是否存在
        if not os.path.exists(rar_path):
            return '<p style="color: red;">文件不存在。请联系 <a href="mailto:hongzhonglu@sjtu.edu.cn">hongzhonglu@sjtu.edu.cn</a> 获取完整数据集。</p>'
        
        # 检查file_content是否已加载
        if file_content is not None:
            file_size_mb = len(file_content) / (1024 * 1024)
            b64 = base64.b64encode(file_content).decode()
            filename = "protein_database.rar"
            href = f'<a href="data:application/x-rar-compressed;base64,{b64}" download="{filename}">📥 Download Protein Database (RAR) - {file_size_mb:.1f} MB</a>'
            return href
        else:
            # 如果file_content为空，尝试直接读取文件
            try:
                with open(rar_path, 'rb') as f:
                    content = f.read()
                file_size_mb = len(content) / (1024 * 1024)
                b64 = base64.b64encode(content).decode()
                filename = "protein_database.rar"
                href = f'<a href="data:application/x-rar-compressed;base64,{b64}" download="{filename}">📥 Download Protein Database (RAR) - {file_size_mb:.1f} MB</a>'
                return href
            except Exception as read_error:
                return f'<p style="color: orange;">读取文件时出错: {str(read_error)}。请联系 <a href="mailto:hongzhonglu@sjtu.edu.cn">hongzhonglu@sjtu.edu.cn</a> 获取完整数据集。</p>'
    except Exception as e:
        return f'<p style="color: red;">下载功能出错: {str(e)}。请联系 <a href="mailto:hongzhonglu@sjtu.edu.cn">hongzhonglu@sjtu.edu.cn</a> 获取完整数据集。</p>'

# 数据库连接
@st.cache_resource
def get_db_connection():
    # 使用相对路径，适用于云部署
    db_path = "lu_web_v3.db"
    if os.path.exists(db_path):
        return sqlite3.connect(db_path)
    else:
        st.error("Database file not found!")
        return None

# 读取数据
@st.cache_data
def load_data():
    conn = get_db_connection()
    if conn is None:
        return None, None, None
    
    try:
        mass_fraction_df = pd.read_sql('SELECT * FROM mass_fraction_combine', conn)
        compartment_df = pd.read_sql('SELECT * FROM compartment_annotation_refine', conn)
        promass_df = pd.read_sql('SELECT * FROM ProMassRatio_across_compartment_combine', conn)
        return mass_fraction_df, compartment_df, promass_df
    except Exception as e:
        st.error(f"Error loading data: {str(e)}")
        return None, None, None

# 加载数据
data_result = load_data()
if data_result[0] is None:
    st.error("Failed to load data from database!")
    st.stop()

mass_fraction_df, compartment_df, promass_df = data_result

# 顶部导航栏
st.markdown("""
<div style='background-color: #f0f2f6; padding: 1rem; border-radius: 3px; margin-bottom: 2rem;'>
    <div style='display: flex; justify-content: space-between; align-items: center;'>
        <h1 style='margin: 0; font-size: 2rem;'>A unified database of yeast absolute quantitative proteomics</h1>
    </div>
</div>
""", unsafe_allow_html=True)

table_choice = st.sidebar.radio(
        "Please select the required function", 
        ["Main Page", "Compute", "Download", "About us"]
    )

if table_choice == "Download":
    st.markdown("## Download Data")
    st.markdown("Download the database in RAR format.")
    st.markdown(get_rar_download_link("protein_database.rar"), unsafe_allow_html=True)

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
        st.subheader("Mass Fraction Data Overview")
        st.dataframe(mass_fraction_df.head(100))

        st.subheader("Mapping of P1-275 Overview")
        # 从数据库加载生理学数据集合
        try:
            conn = get_db_connection()
            if conn is not None:
                try:
                    mapping_df = pd.read_sql('SELECT * FROM physiology_collection', conn)
                    st.dataframe(mapping_df)
                except Exception as e:
                    st.error(f"Error loading mapping data: {str(e)}")
                    st.info("Mapping data could not be loaded from the database.")
            else:
                st.info("Database connection is not available.")
        except Exception as e:
            st.error(f"Error accessing database: {str(e)}")

if table_choice == "Compute":
    module = st.sidebar.selectbox(
        "Select Module",
        ["Compartment Analysis", "Compartment Mass Ratio", "Protein Mass Distribution"]
    )
    
    if module == "Compartment Analysis": # 模块三
        st.subheader("Compartment Analysis")
        compartment = st.selectbox("Select Compartment", compartment_df['compartment'].unique())
        cond = st.selectbox("Select Condition", [f'P{i}' for i in range(1, 276)], key="cond1")
        
        if st.button("Generate Analysis"):
            plot_cumulative_mass_fraction(compartment, cond, mass_fraction_df, compartment_df)
    
    elif module == "Compartment Mass Ratio": # 模块四
        st.subheader("Compartment Mass Ratio Analysis")
        analysis_type = st.radio("Select Analysis Type", ["Single Condition", "Two Conditions"])
        
        if analysis_type == "Single Condition":
            column = st.selectbox("Select Condition", [f'P{i}' for i in range(1, 276)])
            if st.button("Generate Plot"):
                plot_distribution(promass_df, column)
        else:
            col1, col2 = st.columns(2)
            with col1:
                column1 = st.selectbox("Select First Condition", [f'P{i}' for i in range(1, 276)])
            with col2:
                column2 = st.selectbox("Select Second Condition", [f'P{i}' for i in range(1, 276)])
            if st.button("Generate Plot"):
                plot_scatter(promass_df, column1, column2)  
    
    elif module == "Protein Mass Distribution":
        st.subheader("Protein Mass Distribution Analysis")
        analysis_type = st.radio("Select Analysis Type", ["Single Condition", "Two Conditions"])
        
        if analysis_type == "Single Condition":
            column = st.selectbox("Select Condition", [f'P{i}' for i in range(1, 276)])
            if st.button("Generate Plot"):
                plot_distribution_5(mass_fraction_df, column)
        else:
            col1, col2 = st.columns(2)
            with col1:
                column1 = st.selectbox("Select First Condition", [f'P{i}' for i in range(1, 276)])
            with col2:
                column2 = st.selectbox("Select Second Condition", [f'P{i}' for i in range(1, 276)])
            if st.button("Generate Plot"):
                plot_log_scatter_5(mass_fraction_df, column1, column2)  
