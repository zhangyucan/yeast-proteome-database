import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st
import seaborn as sns
import numpy as np
from adjustText import adjust_text

def plot_cumulative_mass_fraction(compartment, cond, mass_fraction_df, compartment_df):
    if compartment in compartment_df['compartment'].unique() and cond in mass_fraction_df.columns:
        nucleus_proteins = compartment_df[compartment_df['compartment'] == compartment]['gene']
        nucleus_mass_fractions = mass_fraction_df[mass_fraction_df['gene'].isin(nucleus_proteins)]
        total_mass_P = nucleus_mass_fractions[cond].sum()
        sorted_nucleus_mass_fractions = nucleus_mass_fractions[['gene', cond]].sort_values(by=cond, ascending=False)
        top_10_proteins = sorted_nucleus_mass_fractions.head(10)
        sorted_nucleus_mass_fractions['cumulative_mass'] = sorted_nucleus_mass_fractions[cond].cumsum()
        nucleus_filtered_genes = sorted_nucleus_mass_fractions[sorted_nucleus_mass_fractions['gene'].isin(nucleus_proteins)]

        # Display the total mass and top 10 proteins
        st.write(f'Total mass of proteins in the nucleus for {cond}: {total_mass_P}')
        st.write(f'Top 10 proteins by mass fraction in {cond}:')
        st.write(top_10_proteins)

        fig, axs = plt.subplots(1, 2, figsize=(15, 6))
        # First subplot: Bar chart of the top 10 proteins by mass fraction in 'P1'
        axs[0].bar(top_10_proteins['gene'], top_10_proteins[cond], color='skyblue')
        axs[0].set_title(f'Top 10 {compartment}-Related Proteins by Mass Fraction in {cond}')
        axs[0].set_xlabel(f'Gene(Top-10 proteins in {cond} mass fraction)')
        axs[0].set_ylabel(f'Mass Fraction ({cond})')
        axs[0].tick_params(axis='x', rotation=45, labelsize=8)  # Rotate gene labels for better visibility

        # Second subplot: Cumulative mass fraction plot for all nucleus-related proteins
        axs[1].plot(nucleus_filtered_genes['gene'], nucleus_filtered_genes['cumulative_mass'], marker='o', linestyle='-', color='yellow')
        axs[1].set_title(f'Cumulative Mass Fraction of Nucleus-Related Proteins in {cond}')
        axs[1].set_xlabel(f'Gene (sorted by {cond} mass fraction)')
        axs[1].set_ylabel('Cumulative Mass Fraction')
        axs[1].tick_params(axis='x', which='both', bottom=True, labelbottom=False)  # Rotate gene labels for better visibility
        axs[1].grid(True,which='major', linestyle='--', linewidth=0.5)
        plt.tight_layout()
        # plt.show()
        st.pyplot(fig)
    else:
        print(f"{compartment}或{cond}列不存在")


# def plot_distribution(data, column):
#     if column in data.columns:
#         # Select column P1 and sort it for the top 20 compartments
#         sorted_data_p = data[['compartment', column]].sort_values(by=column, ascending=False).head(20)
#         # Create a bar plot for the top 20 compartments by P1 protein mass ratio
#         plt.figure(figsize=(10, 6))
#         plt.barh(sorted_data_p['compartment'], sorted_data_p[column], color='skyblue')
#         plt.xlabel(f'Protein Mass Ratio ({column})')
#         plt.ylabel('Compartment')
#         plt.title(f'Top 20 Compartments by Protein Mass Ratio ({column})')
#         plt.gca().invert_yaxis()  # To reverse the order of compartments
#         plt.tight_layout()
#         # plt.show()
#         st.pyplot()
#     else:
#         print(column+"列不存在")

def plot_distribution(data, column):
    if column in data.columns:
        # Select column and sort it for the top 20 compartments
        sorted_data_p = data[['compartment', column]].sort_values(by=column, ascending=False).head(20)
        
        # Create a figure and axis
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Create a bar plot for the top 20 compartments by protein mass ratio
        ax.barh(sorted_data_p['compartment'], sorted_data_p[column], color='skyblue')
        ax.set_xlabel(f'Protein Mass Ratio ({column})')
        ax.set_ylabel('Compartment')
        ax.set_title(f'Top 20 Compartments by Protein Mass Ratio ({column})')
        ax.invert_yaxis()  # To reverse the order of compartments
        plt.tight_layout()
        
        # Display the plot in Streamlit
        st.pyplot(fig)
    else:
        st.error(f"The column '{column}' does not exist.")

def plot_scatter(data, column1, column2):
    if column1 in data.columns and column2 in data.columns:
        # Create a figure and axis
        fig, ax = plt.subplots(figsize=(8, 6))
        
        # Create a scatter plot between the two columns
        ax.scatter(data[column1], data[column2], color='blue', alpha=0.5)

        max_val = max(data[column1].max(), data[column2].max())
        ax.plot([0, max_val], [0, max_val], color='red', linestyle='--', label='y = x')
        ax.legend()

        # 准备标注的标签
        texts = []  # 用于存储标签
        for i, row in data.iterrows():
            if row[column1] >= 0.05 and row[column2] >= 0.05:
                texts.append(ax.text(row[column1], row[column2], row['compartment'], fontsize=8, color='green'))

        # 使用 adjust_text 自动调整标签
        adjust_text(
            texts, 
            arrowprops=dict(arrowstyle='-', color='gray', lw=0.5)  # 为调整的标签添加箭头
        )

        ax.set_xlabel(f'Protein Mass Ratio ({column1})')
        ax.set_ylabel(f'Protein Mass Ratio ({column2})')
        ax.set_title(f'Scatter Plot of Protein Mass Ratios between {column1} and {column2}')
        ax.grid(True)
        plt.tight_layout()
        
        # Display the plot in Streamlit
        st.pyplot(fig)
    else:
        st.error(f"The columns '{column1}' or '{column2}' do not exist.")

# plot_distribution(df, column)
# plot_scatter(df, column1, column2)

def plot_distribution_5(df, column):
    if column in df.columns:
        P1_log = np.log(df[column])
        
        # Create a figure and axis
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Create a histogram with KDE
        sns.histplot(P1_log, kde=True, ax=ax)
        ax.set_title(f'{column} Log Distribution')
        ax.set_xlabel(f'Log({column})')
        ax.set_ylabel('Frequency')
        
        # Display the plot in Streamlit
        st.pyplot(fig)
    else:
        st.error(f"The column '{column}' does not exist.")

# 选中P1和P2列，取出所有数值，对这些数值取log
def plot_log_scatter_5(df, column1, column2):
    if column1 in df.columns and column2 in df.columns:
        P1_log = np.log(df[column1])#.dropna()
        P2_log = np.log(df[column2])#.dropna()
        valid_mask = np.isfinite(P1_log) & np.isfinite(P2_log)    
        P1_log, P2_log = P1_log.align(P2_log, join='inner')
        correlation = np.corrcoef(P1_log[valid_mask], P2_log[valid_mask])[0, 1]

        # Create scatter plot
        fig1, ax1 = plt.subplots(figsize=(10, 6))
        sns.scatterplot(x=P1_log, y=P2_log, ax=ax1)
        ax1.set_title(f'{column1} vs {column2} Log Scatter Plot')
        ax1.set_xlabel(f'Log({column1})')
        ax1.set_ylabel(f'Log({column2})')
        ax1.text(
                x=min(P1_log) + 0.1,  # Adjust text position
                y=max(P2_log) - 0.1,
                s=f'Correlation: {correlation:.3f}',
                fontsize=12,
                color='black',
                bbox=dict(facecolor='white', alpha=0.5, edgecolor='gray')
            )
        
        # Display the scatter plot in Streamlit
        st.pyplot(fig1)

        # Create distribution plots
        fig2, axes = plt.subplots(1, 2, figsize=(15, 6))
        sns.histplot(P1_log, kde=True, ax=axes[0])
        axes[0].set_title(f'{column1} Log Distribution')
        axes[0].set_xlabel(f'Log({column1})')
        axes[0].set_ylabel('Frequency')

        sns.histplot(P2_log, kde=True, ax=axes[1])
        axes[1].set_title(f'{column2} Log Distribution')
        axes[1].set_xlabel(f'Log({column2})')
        axes[1].set_ylabel('Frequency')

        plt.tight_layout()
        
        # Display the distribution plots in Streamlit
        st.pyplot(fig2)
    else:
        st.error(f"The columns '{column1}' or '{column2}' do not exist.")