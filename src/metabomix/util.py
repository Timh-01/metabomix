import os.path
import json
from pathlib import Path
import subprocess 
from typing import Any
from .external import batched,islice
import pandas as pd

def file_to_batches(fn_in: str, batchlength: int) -> list[str]:
    """Splits file into list of strings with every N lines as items"""
    with open(fn_in) as file:
        complete_linelist: list[str] = file.readlines()
        total_linecount: int = len(complete_linelist)
        batch_holder: list[str] = []
        batches = batched(complete_linelist,batchlength)
        for batch in batches: 
            sublist_linecount = len(batch)
            linecounter = 0
            query_input: str = ""
            for line in batch:
                if linecounter != sublist_linecount:       
                    query_input += line[:-1]+"\\n"
                    linecounter += 1
                else:
                    batch += line[-1]
            batch_holder.append(query_input)
    return batch_holder

def run_subprocess(command: str,**kwargs) -> None:
    """Run a specific command from the shell"""
    if "cwd" in kwargs:
        cwd = kwargs["cwd"]
        subprocess.run([command], cwd=cwd, shell=True, capture_output=True, text=True)
    else:
        subprocess.run([command], shell=True, capture_output=True, text=True)

def validate_dictkeys(to_check: dict,required_keys: list[str]):
    """Check if all items from a list are dictkeys"""
    return all(required_key in to_check for required_key in required_keys)

def json_file_to_dict(fn):
    """Converts a file in Json format to a dict"""
    with open(fn,'r') as file:
        content: dict = json.load(file)
    return content

def json_from_string_or_file(input_json) -> dict:
        """returns input json from file or dict as dict"""
        # self.paths,self.tools,self.other = workflowsettings.values()
        if isinstance(input_json,str):
            try:
                os.path.exists(input_json)
            except FileNotFoundError: 
                print(f"Settings file not found: {input_json}")
            with open(input_json,'r') as file:
                content: dict = json.load(file)
        elif isinstance(input_json,dict):
            content: dict = input_json
        else: 
            raise TypeError(f"json has to be file or dict, not {type(input_json)}")
        return content

def get_filelist_from_folder(folder: str) -> list:
    """Get all files from a folder as list"""
    p = Path(folder)
    filelist: list = [file.as_posix() for file in p.iterdir()]
    return filelist


# l_pos_fn, l_neg_fn = split_files_per_mode(data_folder)
#     for paramset in ["strictest","intermediate","lenient","manual"]:
#         write_mzbatch(l_files=l_pos_fn,mzbatch_fn=basis_mzbatch,metadata_fn=fn_metadata,paramset=paramset,output_folder=mzbatch_output_loc)

def query_df_value_by_other_col(query,df:pd.DataFrame,query_col,target_col) -> Any:
    """Checks if query col contains query and if so returns value from target col"""
    results = df.loc[:,[target_col]][df[query_col] == query]
    #df.loc[:,[query_col]][df  == query]
    return results.iat[0,0] if len(results) else "N/A"

def integrate_df_col_to_df(target:pd.DataFrame,t_vals,t_queries,source:pd.DataFrame,s_vals,s_queries) -> pd.DataFrame:
    """Integrate value column from source into target if query columns match"""
    target[t_vals] = target[t_queries].map(lambda t_query: 
                query_df_value_by_other_col(query=t_query,df=source,query_col=s_queries,target_col=s_vals))     
    return target

def integrate_df_cols_to_df(target_df:pd.DataFrame,source_df:pd.DataFrame,
                                target_integration_col: Any,source_integration_col: Any,
                                col_translator:dict[Any,Any]) -> pd.DataFrame:
    """Integrate values columns from source into target if query columns match for value columns in col_translator"""
    for old_name,new_name in col_translator.items():
        target_info = target_df,new_name, target_integration_col
        source_info = source_df.loc[:,[old_name,source_integration_col]],old_name,source_integration_col
        target_df = integrate_df_col_to_df(*target_info,*source_info)
    return target_df

def convert_missing(x):
    """Helper to convert missing value to the same format"""
    if (x=="nan") or (x.lower() == "n/a"):
        return "N/A"
    else:
        return x