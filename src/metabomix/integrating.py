import os
import numpy as np
import json
from typing import Any

from matchms.importing import load_from_mgf
from matplotlib import pyplot as plt
from matchms.filtering import remove_peaks_around_precursor_mz, remove_peaks_outside_top_k,select_by_mz
from matchms import calculate_scores
import matchms.Scores
from matchms.similarity import ModifiedCosine
import networkx as nx
# import matchmsextras.networking as net
from matchms.networking import SimilarityNetwork
import pandas as pd
from typing import Any,Iterable
import rdkit.Chem.PandasTools

try:
    import MS2LDA
except ImportError:
     MS2LDA = None

def network_to_edgelist_and_nodes_df(graph: nx.Graph) -> tuple[pd.DataFrame, pd.DataFrame]:
    """splits a networkx graph into an edge and nodes dataframe"""
    nodes = graph.nodes()
    for node in nodes:
      nodes[node]["placeholder"] = "placeholder"
    #    networkdf = pd.DataFrame.from_dict(dict(graph.nodes(data=True)), orient='index')  
    networkdf = pd.DataFrame.from_dict(nodes,orient='index')      

    edgelist = nx.to_pandas_edgelist(graph)
    return networkdf, edgelist


def integrate_sirius_to_graph(graph: nx.Graph,translators: tuple[dict,dict,dict],outputs: tuple[str,str,str]) -> nx.Graph:
    """Merge results from sirius tools to graph"""
    sirius_translator,fingerid_translator,canopus_translator = translators
    sirius_file,fingerid_file,canopus_file = outputs
    sirius_df = pd.read_csv(sirius_file, sep='\t')
    fingerid_df = pd.read_csv(fingerid_file, sep="\t")
    canopus_df = pd.read_csv(canopus_file, sep="\t")

    graph = df_to_graph(graph,canopus_df,canopus_translator)
    graph = df_to_graph(graph,sirius_df,sirius_translator)
    graph = df_to_graph(graph,fingerid_df,fingerid_translator)

    return graph

def df_to_graph(graph: nx.Graph,to_merge: pd.DataFrame,translator: dict[str,str]) -> nx.Graph:
    """Adds columns from dataframe to nodes, merging on "scan" columns from DF"""
    for result in to_merge.to_dict(orient="records"):
        scan = str(result["mappingFeatureId"])
        if scan in graph:
            for oldname,newname in translator.items():
                graph.nodes[scan][newname] = result[oldname]
    return graph

def integrate_ms2lda_to_df(dataframe: pd.DataFrame, dataset: str,output_folder: str,output_database: str) -> pd.DataFrame:
    """Get motifs per feature and use feature ID to merge motifs to network"""
    motifDB_query = "QUERY scaninfo(MS2DATA)"
    screening_hits = MS2LDA.screen_spectra(motifDB=output_database, dataset=dataset, motifDB_query=motifDB_query, output_folder=output_folder)
    with open(dataset) as file:
        A = file.read()
        feature_ids = [entry.split('\n')[1][11:] for entry in A.split("BEGIN IONS")[1:]]
        translation_dict = {}
        for i,j in enumerate(feature_ids):
            translation_dict[f'spec_{i}'] = j
        screening_hits["scan"] = screening_hits["hit_id"].map(lambda x: translation_dict[x])
        d_motif_per_scan = {}
        for i in screening_hits["scan"].unique():
            d_motif_per_scan[i] = tuple(screening_hits["ref_motif_id"][screening_hits["scan"] == i])
        dataframe['temp_indexcopy'] = dataframe.index
        dataframe["MS2LDA:allmotifs"] = dataframe['temp_indexcopy'].map(lambda x: map_motifs(x,d_motif_per_scan))
        dataframe = dataframe.drop('temp_indexcopy',axis=1)
    return dataframe
# integrate_df_values_if_shared_other_value(integrated_df,plastchem_df,canonical_smiles,csifingerid:smiles,plastchem:hazard_score,hazard_score)
# def integrate_df_values_if_shared_other_value(df1,df2,col_sh1,colsh2,newcol,colint2):
# #    if sh1 = sh2, col1 = col2

def create_network_from_mgf(mgf: str) -> nx.Graph:
    """process spectra from MGF, calculate similarity and create network from those scores"""
    spectrums: list[matchms.Spectrum] = load_from_mgf(mgf)
    spectrums = [peak_processing(s) for s in spectrums] 
    similarity_measure = ModifiedCosine(tolerance=0.02)
    scores: matchms.Scores = calculate_scores(spectrums, spectrums, similarity_measure,is_symmetric=True)
    scores.to_pickle(f'{mgf[:-4]}_scores.pickle')
    # scores.filter_by_range(name="ModifiedCosine_matches",low=6,below_operator="<")

    network: nx.Graph = create_network_from_scores(scores)
    return network

def peak_processing(spectrum: matchms.Spectrum) -> matchms.Spectrum:
    """Apply filters to spectrum"""
    spectrum = remove_peaks_around_precursor_mz(spectrum,mz_tolerance=17)
    spectrum = remove_peaks_outside_top_k(spectrum,k=6,mz_window=50)
    return spectrum

def create_network_from_scores(scores: matchms.Scores ) -> nx.Graph:
    """Filters scores and creates network"""
    #scores: matchms.Scores = scores_from_pickle("/lustre/BIF/nobackup/hendr218/temp/scores_001")
    newscores = scores
    newscores.filter_by_range(name="ModifiedCosine_matches",low=6,below_operator="<")
    if "scans" in newscores.queries[0].metadata:
        ms_network = SimilarityNetwork(identifier_key="scans",score_cutoff=0.7,top_n=10,link_method="mutual")
    elif "feature_id" in newscores.queries[0].metadata:
         ms_network = SimilarityNetwork(identifier_key="feature_id",score_cutoff=0.7,top_n=10,link_method="mutual")
    else:
        raise ValueError("'Scans' or 'feature_id' required in mgf to identify features in network by other tools")
    ms_network.create_network(newscores,score_name="ModifiedCosine_score")
    return ms_network


def parse_cramer_classifications(fn_cramer_csv: str) -> pd.DataFrame:
    """parses cramer classification and smiles per compound from input csv"""
    with open(fn_cramer_csv) as file:
        l_cramer = []
        for line in file:
            #skip the first line
            if line.startswith('CRAMERFLAGS,Cramer rules'):
                continue
            cramer_class,x,cramer_smiles = line.split(',')[1:4]
            l_cramer.append((cramer_class,cramer_smiles))
    df_cramer = pd.DataFrame(l_cramer,columns=("cramer_classification","SMILES"))
    return df_cramer

def get_merging_settings(config: dict,target_integration_col: str) -> tuple[str,str,dict]:
        """Gets columns to integrate on and column translations from config"""
        translations: dict[Any] = config.get("translations", {})
        merging_transl: dict[Any]  = translations.get("merging", {})
        column_transl: dict[Any] | str = translations.get("columns", "N/A")
        source_integration_col: Any | str = merging_transl.get(target_integration_col,"N/A")
        return target_integration_col,source_integration_col,column_transl
#target_int_col,source_int_col,col_translator

def map_motifs(x: Any, motifs_per_scan: dict[str,list]) -> dict | str:
    """Helper to map finding motifs per scan for MSLDA to dataframe"""
    if x in motifs_per_scan:
        return motifs_per_scan[x]
    else:
        return "N/A"

def parse_classyfire_sdf(sdf_fn: str) -> pd.DataFrame:
    """Converts SDF from classyfire to dataframe containg classifications and inchis"""
    cf_columns = ("inchi","smiles","cf_kingdom","cf_superclass","cf_class","cf_subclass")
    df_sdf = rdkit.Chem.PandasTools.LoadSDF(sdf_fn)
    df_sdf['new'] = df_sdf['InChIKey'].map(parse_cf_string)
    parsed_df = pd.DataFrame(df_sdf['new'].tolist(),columns=(cf_columns))
    return parsed_df
    
def parse_cf_string(cf_entry: str):
    """Takes classyfire SDF output entry and parses classifications and inchi"""
    split_entry = cf_entry.split('> <Intermed')[0].split('\n')
    match len(split_entry):
        case 12:
            inchi,a,smiles,b, cf_kingdom, c, cf_superclass,d,cf_class,e,cf_subclass,f = split_entry
        case 10:
            inchi,a,smiles,b, cf_kingdom, c, cf_superclass,d,cf_class,f = split_entry
            cf_subclass = "N/A"
        case 8:
            inchi,a,smiles,b, cf_kingdom, c, cf_superclass,f = split_entry
            cf_subclass = cf_class = "N/A"
        case 6:
            inchi,a,smiles,b, cf_kingdom,f = split_entry
            cf_subclass = cf_class = cf_superclass = "N/A"
        case _:
            inchi = cf_kingdom = cf_subclass = cf_class = cf_superclass = "N/A"
    return (inchi,smiles,cf_kingdom,cf_superclass,cf_class,cf_subclass)   
# plastchem_searcher = DfSearcher(datab=plastchem_db) 
# col_query = integrated_df.loc[:,["csifingerid:smiles"]]
# for new_name,plastchem_name in [
#             ("plastchem:hazard_score","Hazard_score"),
#             ("plastchem:PlastChem_lists","PlastChem_lists"),
#             ("plastchem:use","Use"),
#             ("plastchem:pcname","pubchem_name")]:
    
#     integrated_df[new_name] = col_query.map(lambda x: p.process(x,prop=plastchem_name))



def temp_metadata_adder(dataframe: pd.DataFrame,input_mgf: str,metadata_csv: str,quant_table: str) -> nx.Graph:
    """takes input mgf,df and adds metadata"""
    with open(input_mgf) as file:
        A = file.read()
        entries = A.split("BEGIN IONS")[1:]

    id_files_dict = {}
    # rt_dict = {}
    # mass_dict = {}
    for entry in entries:
        feature_id = entry.split('\n')[1][11:] 
        if "FILENAME" in entry:
            files = entry.split("FILENAME")[1][1:].split('\n')[0]
            id_files_dict[feature_id] = files.split(";")

    metadata_dict = {}
    with open(metadata_csv) as file:
        for line in file.readlines()[1:]:
            name,plastic_type,plastic_state = line.split(',')
            plastic_state=plastic_state[:-1]
            polarity = "pos" if "pos" in name else "neg"
            if plastic_state == "Blank":
                state_type = "blank"
            elif plastic_state == "procBlank":
                state_type = "proc_blank"
            else:
                state_type = f"{plastic_type}_{plastic_state}"
            if not name in metadata_dict:
                metadata_dict[name] = state_type      

    id_statetype = {}
    for id, files in id_files_dict.items():
        id_statetype[id] = ""
        for file in files: 
            statetype = metadata_dict[file]
            if not statetype in id_statetype[id]:
                id_statetype[id] += f"{statetype}, "

    #  for i in id_statetype["id"].unique():
    
    quant_table =  pd.read_csv(quant_table,low_memory=False)
    samples = [key for key in quant_table.keys() if ("mzML:area" in key) & (not "Blank" in key)]
    samples.append('id')
    blanks = [key for key in quant_table.keys() if ("mzML:area" in key) & ("Blank" in key)]
    blanks.append('id')
    s_areas = quant_table.loc[:,samples].set_index('id').transpose()
    b_areas = quant_table.loc[:,blanks].set_index('id').transpose()

    s_area_tresholded = s_areas
    is_blank_dict: dict[str,bool] = {}
    for id in s_areas.keys():
        mean = s_areas[id].mean()
        s_area_tresholded[id] = s_areas[id].map(lambda x: True if mean < 0.3 * x else False)
    for id in s_areas.keys():
        mean = quant_table.loc[:,samples].set_index('id').transpose()[id].fillna(0).mean()
        blank_max = quant_table.loc[:,blanks].set_index('id').transpose()[id].fillna(0).max()
        is_blank: bool = mean<10*blank_max
        is_blank_dict[id] = is_blank
    dataframe["is_blank"] = dataframe.index.map(lambda x: is_blank_dict[int(x)]) 
    
    # feature_metadata_dict = {}
    # for id, files in s_area_tresholded.to_dict().items():
    #     id = str(id)
    #     feature_metadata_dict[id] = ""
    #     for file,val in files.items(): 
    #         if (val) & (not file in feature_metadata_dict[id]):
    #             feature_metadata_dict[id] += f"{metadata_dict[file[9:-5]]}, "
    feature_metadata_dict = {}
    for id, files in s_area_tresholded.to_dict().items():
        id = str(id)
        feature_metadata_dict[id] = []
        for file,val in files.items(): 
            metadata: str = metadata_dict[file[9:-5]]
            if (val) & (metadata not in feature_metadata_dict[id]):
                feature_metadata_dict[id].append(metadata)

    dataframe["all_metadata"] = dataframe.index.map(lambda x: feature_metadata_dict.get(x,"N/A"))
    for metadata_type in ["PE_PA_Dec", "PE_PA_PC","PE_Dec", "PE_PC", "PE_PET_PC", "PE_PET_Dec"]:
        dataframe[metadata_type] =  dataframe["all_metadata"].map(lambda x: metadata_type in x)
    for plastic_type in ["PE","PE_PA","PE_PET"]:
        dataframe[plastic_type] = dataframe["all_metadata"].map(lambda all_metadata: True if plastic_type in ["_".join(metadata_item.split("_")[0:-1]) for metadata_item in all_metadata] else False)
    for colname,plastic_state in [("Decontaminated","Dec"),("Postconsumer","PC")]:
        dataframe[colname] = dataframe["all_metadata"].map(lambda all_metadata: True if plastic_state in [metadata_item.split("_")[-1] for metadata_item in all_metadata] else False)
    dataframe["Introduced_by_dec"] = dataframe.apply(lambda row: introduced_by_decontamination_mapper(row),axis=1)
    dataframe["Removed_by_dec"] = dataframe.apply(lambda row: removed_by_decontamination_mapper(row),axis=1)
    dataframe["Kept_by_dec"] = dataframe.apply(lambda row: kept_by_decontamination_mapper(row),axis=1)

    # dataframe["PE"] = dataframe["all_metadata"].map(lambda x: True if ("PE_Dec" in x or "PE_PC" in x) else False)
    # dataframe["plastic_types"] = dataframe["all_metadata"].apply(plastic_type_mapper)
    dataframe["plastic_types"] = dataframe.apply(plastic_type_mapper,axis=1)
    return dataframe

def introduced_by_decontamination_mapper(row)->bool:
    metadata_items = [metadata_item.split("_")[-1] for metadata_item in row["all_metadata"]]
    if "Dec" in metadata_items and not "PC" in metadata_items: 
        return True
    else:
        return False
    
def removed_by_decontamination_mapper(row)->bool:
    metadata_items = [metadata_item.split("_")[-1] for metadata_item in row["all_metadata"]]
    if "PC" in metadata_items and not "Dec" in metadata_items: 
        return True
    else:
        return False
    
def kept_by_decontamination_mapper(row)->bool:
    metadata_items = [metadata_item.split("_")[-1] for metadata_item in row["all_metadata"]]
    if "PC" in metadata_items and "Dec" in metadata_items: 
        return True
    else:
        return False

def plastic_type_mapper(row) -> str:
    if row["PE_PA"] & (row["PE_PET"]):
        return "PE-PET-PA" if row["PE"] else "PET-PA"
    if row["PE_PA"]:
        return "PE-PA" if row["PE"] else "PA"
    if row["PE_PET"]:
        return "PE-PET" if  row["PE"] else "PET"
    if row["PE"]:
        return "PE"
    else:
        return "N/A"
# def plastic_type_mapper(row) -> str:
#     if row(["PE_PA"]) & row["PE_PET"]:
#     if ("PE_PA" in all_types) & ("PE_PET" in all_types):
#         return "PE-PET-PA" if ("PE" in all_types) else "PET-PA"
#     if ("PE_PA" in all_types):
#         return "PE-PA" if ("PE" in all_types) else "PA"
#     if ("PE_PET" in all_types):
#         return "PE-PET" if  ("PE" in all_types) else "PET"
#     if "PE" in all_types:
#         return "PE"
#     else:
#         return "N/A"



    # p = DfSearcher(datab=plastchem_db) 
    # col_query = integrated_df.loc[:,["csifingerid:smiles"]]
    # for new_name,plastchem_name in [
    #         ("plastchem:hazard_score","Hazard_score"),
    #         ("plastchem:PlastChem_lists","PlastChem_lists"),
    #         ("plastchem:use","Use"),
    #         ("plastchem:pcname","pubchem_name")]:
    #     integrated_df[new_name] = col_query.map(lambda x: p.process(x,prop=plastchem_name))