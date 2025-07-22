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
    file_path = "/home/yucan/yeast-proteome-database/protein_database.rar"
    with open(file_path, 'rb') as f:
        file_content = f.read()

# 在后台异步准备文件
executor.submit(prepare_file)

@st.cache_data
@lru_cache(maxsize=5)
def get_rar_download_link(rar_path):
    try:
        b64 = base64.b64encode(file_content).decode()
        filename = "protein_database.rar"
        href = f'<a href="data:application/x-rar-compressed;base64,{b64}" download="{filename}">Download Protein Database (RAR)</a>'
        return href
    except Exception as e:
        return f"Error: {str(e)}"

# 数据库连接
@st.cache_resource
def get_db_connection():
    return sqlite3.connect("/home/yucan/yeast-proteome-database/lu_web_v3.db")

# 读取数据
@st.cache_data
def load_data():
    conn = get_db_connection()
    mass_fraction_df = pd.read_sql('SELECT * FROM mass_fraction_combine', conn)
    compartment_df = pd.read_sql('SELECT * FROM compartment_annotation_refine', conn)
    promass_df = pd.read_sql('SELECT * FROM ProMassRatio_across_compartment_combine', conn)
    return mass_fraction_df, compartment_df, promass_df

# 加载数据
mass_fraction_df, compartment_df, promass_df = load_data()

# 顶部导航栏
st.markdown("""
<div style='background-color: #f0f2f6; padding: 1rem; border-radius: 3px; margin-bottom: 2rem;'>
    <div style='display: flex; justify-content: space-between; align-items: center;'>
        <h1 style='margin: 0; font-size: 2rem;'>A unified database of yeast absolute quantitative proteomics</h1>
    </div>
</div>
""", unsafe_allow_html=True)

# st.markdown("""
# <div style='background-color: #f0f2f6; padding: 1rem; border-radius: 3px; margin-bottom: 2rem;'>
#     <div style='display: flex; justify-content: space-between; align-items: center;'>
#         <h1 style='margin: 0; font-size: 1.5rem;'>Protein Mass Fraction Analysis</h1>
#     </div>
# </div>
# """, unsafe_allow_html=True)

table_choice = st.sidebar.radio(
        "Please select the required function", 
        ["Main Page", "Compute", "Download", "About us"]
    )

if table_choice == "Download":
    st.markdown("## Download Data")
    st.markdown("Download the database in RAR format.")
    st.markdown(get_rar_download_link("/home/yucan/lu_web/protein_database.rar"), unsafe_allow_html=True)

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
                        #plt.bar(filtered_df['gene'], filtered_df[col], label=col, alpha=0.5)
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
        Mapping_df = pd.read_csv("/home/yucan/lu_web/physiology_collection_v2.csv")
        st.dataframe(Mapping_df)

if table_choice == "Compute":
    #if st.sidebar.checkbox("Show Computation Modules"):
    module = st.sidebar.selectbox(
        "Select Module",
        ["Compartment Analysis", "Compartment Mass Ratio", "Protein Mass Distribution"]
    )
    
    if module == "Compartment Analysis": # 模块三
        st.subheader("Compartment Analysis")
        compartment = st.selectbox("Select Compartment", compartment_df['compartment'].unique())
        # col1, col2 = st.columns(2)
        # with col1:
        cond = st.selectbox("Select Condition", [f'P{i}' for i in range(1, 276)], key="cond1")
        # mass_fraction_df_rename = pd.read_csv("/home/yucan/lu_web/mass_fraction_combine.csv")
        # with col2:
        #     cond2 = st.selectbox("Select Condition", mass_fraction_df_rename.columns[1:], key="cond2")
        #     cond2 = mass_fraction_df.columns[mass_fraction_df_rename.columns[1:].get_loc(cond2) + 1]


        # cond = st.selectbox(f"Select Condition in {cond1} and {cond2}", [cond1, cond2],key="cond")
        
        if st.button("Generate Analysis"):
            # st.write
            plot_cumulative_mass_fraction(compartment, cond, mass_fraction_df, compartment_df)
    
    elif module == "Compartment Mass Ratio": # 模块四
        st.subheader("Compartment Mass Ratio Analysis")
        analysis_type = st.radio("Select Analysis Type", ["Single Condition", "Two Conditions"])
        
        if analysis_type == "Single Condition":
            column = st.selectbox("Select Condition", [f'P{i}' for i in range(1, 276)])
            if st.button("Generate Plot"):
                # 实现模块四的单条件分析
                plot_distribution(promass_df, column)
                pass
        else:
            col1, col2 = st.columns(2)
            with col1:
                column1 = st.selectbox("Select First Condition", [f'P{i}' for i in range(1, 276)])
            with col2:
                column2 = st.selectbox("Select Second Condition", [f'P{i}' for i in range(1, 276)])
            if st.button("Generate Plot"):
                # 实现模块四的双条件分析
                plot_scatter(promass_df, column1, column2)  
                pass
    
    elif module == "Protein Mass Distribution":
        st.subheader("Protein Mass Distribution Analysis")
        analysis_type = st.radio("Select Analysis Type", ["Single Condition", "Two Conditions"])
        
        if analysis_type == "Single Condition":
            column = st.selectbox("Select Condition", [f'P{i}' for i in range(1, 276)])
            if st.button("Generate Plot"):
                plot_distribution_5(mass_fraction_df, column)
                pass
        else:
            col1, col2 = st.columns(2)
            with col1:
                column1 = st.selectbox("Select First Condition", [f'P{i}' for i in range(1, 276)])
            with col2:
                column2 = st.selectbox("Select Second Condition", [f'P{i}' for i in range(1, 276)])
            if st.button("Generate Plot"):
                plot_log_scatter_5(mass_fraction_df, column1, column2)  
                pass
        pass


# About 页面
# if st.sidebar.checkbox("Show About"):
#     st.markdown("""
#     ## About Us
    
#     This is a placeholder text for the About section. Please replace with actual content.
#     """)




