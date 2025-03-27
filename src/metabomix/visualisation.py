import pandas as pd
import matplotlib.pyplot as plt
from typing import Any,Iterable

# def plot_classes(df_in,level):
#     df_in[f'cf_{level}']
#     df_superclass = df_in[f'cf_{level}']
#     #df_superclass = df_sdf_relevant.map(lambda x: x.split("<Superclass>\n")[1].split("\n")[0]) #Parse superclass
#     plastchem_superclasscounts = df_superclass.value_counts(normalize=True)
#     d_plastchem_superclasscounts = plastchem_superclasscounts.to_dict()
#     d_temp = {}
#     for i in d_plastchem_superclasscounts:
#         d_temp[i] = d_plastchem_superclasscounts[i]
#     d_plastchem_superclasscounts = d_temp
#     df_pcc = pd.DataFrame.from_dict(d_plastchem_superclasscounts,orient="index",columns=["count"])
#     data_superclasscounts = integrated_df[f'canopus:CF_{level}'].value_counts(normalize=True)
#     d_data_superclasscounts = data_superclasscounts.to_dict()
#     df_dsc = pd.DataFrame.from_dict(d_data_superclasscounts,orient="index",columns=["count"])
#     df_joined = df_dsc.join(df_pcc,lsuffix="_data",rsuffix="_plastchem", how="outer")
#     df_joined = df_joined.fillna(0)
#     ax = df_joined[['count_data','count_plastchem']].plot(kind='bar', title ="Chemical space", figsize=(15, 10), legend=True, fontsize=12)
#     ax.set_xlabel(level, fontsize=12)
#     fig_name=f"/lustre/BIF/nobackup/hendr218/temp/Plots/plastchem_data_{level}.png"
#     plt.savefig(fig_name)


def create_counts_df(df: pd.DataFrame,col_to_plot: str,col_to_filter:str=None,normalize:bool=True) -> pd.DataFrame:
    """Creates a dataframe for """
    if col_to_filter:
        df = df[col_to_plot][(df[col_to_filter] == True) & (df["is_blank"] == False)].value_counts(normalize=normalize)
    else: 
        df = df[col_to_plot][df["is_blank"] == False].value_counts(normalize=normalize)
    df=df.to_dict()
    df = pd.DataFrame.from_dict(df,orient="index",columns=[col_to_plot])
    return df

#fig_name=f"/lustre/BIF/nobackup/hendr218/Data/cf_only_{vis_col}_by_plastic_type.png"

# def plot_count_cols(fig_name:str):
#     for vis_col in ["CF:subclass", "CF:class", "CF:superclass"]:
#         nw = onw[onw[vis_col] != "N/A"]
#         PE_PET = create_vc_df(nw,vis_col,col_to_filter="PE_PET")
#         PE_PA = create_vc_df(nw,vis_col,col_to_filter="PE_PA")
#         PE = create_vc_df(nw,vis_col,col_to_filter="PE")
#         ALL = create_vc_df(nw,vis_col)

#         df_joined = PE.join(PE_PET,lsuffix="_PE",rsuffix="_PET", how="outer")
#         df_joined2 = ALL.join(PE_PA,lsuffix="_ALL",rsuffix="_PA", how="outer")
#         df_joined = df_joined2.join(df_joined,how="outer").fillna(0)
#         df_joined = df_joined.apply(lambda x: x*100)
#         df_joined = df_joined.rename(columns=lambda x: x.split("_")[-1])
#         if vis_col == "subclass":
#             df_joined = df_joined[df_joined.apply(lambda row: any([i > 1  for i in row]),axis=1)] 


# def join_counts(*dataframes: Iterable[pd.DataFrame]):
#     if len(args) == 1:
#         return args
#     dfs = [dataframe for dataframe in dataframes]
#     df_joined = dfs[0]
#     dfs.pop[0]
#     for i,df in dfs.enumerate():
#         df_joined = df_joined.join(df,lsuffix=f"_{i-1}",rsuffix=f"_{i}",how="outer")
#     df_joined.fillna(0)


# def df_to_barplot(df:pd.DataFrame,columns:list[str],title: str,ytitle:str,fn: str) -> None:
#     ax = df[columns].plot(kind='bar', title=title, figsize=(15, 10), legend=True, fontsize=12)
#     ax.set_ylabel(ytitle, fontsize=12)
#     plt.tight_layout()
#     plt.savefig(fn)