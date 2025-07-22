import sqlite3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st  
column = 'P1'
column1 = 'P1'
column2 = 'P2'

# 从数据库中读取数据
conn = sqlite3.connect("/home/yucan/lu_web_v2.5.db")
df = pd.read_sql('SELECT * FROM mass_fraction_combine', conn)

# 选中P1列，取出所有数值，对这些数值取log，画出分布图
def plot_distribution(df, column):
    if column in df.columns:
        P1_log = np.log(df[column])
        plt.figure(figsize=(10, 6))
        sns.histplot(P1_log, kde=True)
        plt.title(column+' Log Distribution')
        plt.xlabel(f'Log({column})')
        plt.ylabel('Frequency')
        plt.savefig('./distributin.png')
        # 使其在streamlit中显示
        st.pyplot()

        # plt.show()
    else:
        print(column+"列不存在")

# 选中P1和P2列，取出所有数值，对这些数值取log
def plot_log_scatter(df, column1, column2):
    if column1 in df.columns and column2 in df.columns:
        P1_log = np.log(df[column1])
        P2_log = np.log(df[column2])

        # 画出P1和P2的散点图
        plt.figure(figsize=(10, 6))
        sns.scatterplot(x=P1_log, y=P2_log)
        plt.title(f'{column1} vs {column2} Log Scatter Plot')
        plt.xlabel(f'Log({column1})')
        plt.ylabel(f'Log({column2})')
        plt.savefig('./scatter.png')
        st.pyplot()
        # plt.show()

        # 画出P1和P2的分布图进行比较
        fig, axes = plt.subplots(1, 2, figsize=(15, 6))
        sns.histplot(P1_log, kde=True, ax=axes[0])
        axes[0].set_title(column1+' Log Distribution')
        axes[0].set_xlabel(f'Log({column1})')
        axes[0].set_ylabel('Frequency')

        sns.histplot(P2_log, kde=True, ax=axes[1])
        axes[1].set_title(column2+' Log Distribution')
        axes[1].set_xlabel(f'Log({column2})')
        axes[1].set_ylabel('Frequency')

        plt.tight_layout()
        plt.savefig('./compare.png')
        st.pyplot(fig)
        # plt.show()
    else:
        print(f"{column1}或{column2}列不存在")

plot_distribution(df, column)
plot_log_scatter(df, column1, column2)
