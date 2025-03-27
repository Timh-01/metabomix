"""Module providings function for running various mass spectrometry tools

Included: mzmine, sirius, cramer, classyfire

"""

import json
from time import sleep
import requests
from typing import Any
import json
from time import ctime
from typing import Any
#from src.myworkflow.integrating import *
from .parsing import *
# from src.myworkflow.visualisation import *
from .util import validate_dictkeys,file_to_batches,run_subprocess,json_from_string_or_file,get_filelist_from_folder
from pathlib import Path
import xml.etree.ElementTree as ET
import os.path

import numpy as np
#import pymolnetenhancer as pme
import pandas as pd
from urllib.request import urlretrieve
import requests
import json
import networkx as nx
import subprocess
import rdkit
from time import sleep
import rdkit.Chem 
import rdkit.Chem.PandasTools
import matplotlib.pyplot as plt
import itertools
try:
    import MS2LDA
except ImportError:
    MS2LDA = None
import xml.etree.ElementTree as ET
from pathlib import Path

   
class mzbatch_writer():
        """Class to replace values from mzbatch with custom values"""
        searchterm_translator = {
        "min_samples":      (".//parameter[@name='Minimum aligned features (samples)']/abs",
                            ".//parameter[@name='Min samples filter']/parameter[@name='Min samples in all']/abs",),
        "feature_height":  (".//parameter[@name='Minimum absolute height']",),
        "RT_start":        (".//parameter[@name='Retention time']/min",),
        "RT_end":          (".//parameter[@name='Retention time']/max'",),
        "RTsearchrange":   (".//parameter[@name='Minimum search range RT/Mobility (absolute)']",),
        "top_edge_ratio":  (".//parameter[@name='Min ratio of peak top/edge']",),
        "chromtreshold":   (".//parameter[@name='Chromatographic threshold']",),
        }
        
        def __init__(self,base_mzbatch: str, data: list[str], mzmine_params: dict) -> None:
            self.base_mzbatch = base_mzbatch
            self.data = data
            self.params = mzmine_params


            tree = ET.ElementTree()
            tree.parse(self.base_mzbatch)
            self.tree = tree
            self.root: ET.Element = tree.getroot()

        def write_mzbatch(self,output_filename) -> str:
            """replace various parameters and write resulting batchfile"""
            self.replace_filenames()
            self.replace_params()
            
            # metadata=self.metadata
            # self.replace_metadata(root,metadata_fn)

            ET.indent(self.tree)
            self.tree.write(output_filename)

            return output_filename

        def replace_params(self) -> None:
            """look up mzbatch description for available params and use to set new values"""
            provided_params: list = self.get_provided_params(available_params=mzbatch_writer.searchterm_translator.keys())
            for provided_param in provided_params:
                searchterms = mzbatch_writer.searchterm_translator[provided_param] 
                value = self.params.get(provided_param)
                for searchterm in searchterms:
                    self.replace_param(param_name=searchterm,new_value=value)
               
  

        def get_provided_params(self,available_params) -> list:
            """Check which available mzmine params were selected in settings file"""
            provided_params: list = [param for param in available_params if self.params.get(param) == "True"]
            return provided_params

           
        
        def replace_param(self,param_name,new_value):
            """Replace element in mzbatch with new value"""
            for elem in self.do_mzbatch_search(param_name):
                elem.text = new_value

        def do_mzbatch_search(self,searcher_holder: list[str]) -> Any:
            "Find all instances of an element in mzbatch"
            all_findings: list[ET.Element] = []
            for searchterm in searcher_holder:
                for elem in self.root.findall(searchterm):
                    all_findings.append(elem)
                for elem in all_findings:
                    yield elem

        def replace_metadata(self,metadata_fn: str) -> None:
            """Replace metadata file with file from settings"""
            metadata_element: ET.Element | None = self.root.find(".//parameter[@name='Metadata file']/current_file") #root.findall("./batchstep/parameter/..[@name='File Names]"):
            if metadata_element is None:
                raise Settingserror("MZmine batchfile should contain 'Metadata file' parameter")
            metadata_element.text = metadata_fn   

        def replace_filenames(self) -> None:
            """Replace filenames from base mzbatch with custom ones"""
            oldfiles: list = self.root.findall(".//*[@name='File names']/") 
            fileloc: ET.Element | None = self.root.find(".//*[@name='File names']")
            if fileloc is None:
                raise Settingserror("MZmine batchfiles should contain 'File names' parameter")
            for oldfile in oldfiles:
                fileloc.remove(oldfile)
            for newfile in self.data:
                file_element = ET.SubElement(fileloc,"file")
                file_element.text=newfile



def run_mzmine(mzmine_location,mzmine_userfile_location,batchfn,temp_folder,mzmine_output_loc) -> bool:
    subprocess_command =  f"{mzmine_location} -user {mzmine_userfile_location} -b {batchfn} -o {mzmine_output_loc} -t {temp_folder}"
    run_subprocess(subprocess_command)
    return True
    
def run_sirius(input_path: str, output_path: str, sirius_path: str, **kwargs) -> bool:
    """"Run sirius from the shell

    Parameters
    ----------
    input_path
        file path to .MGF file, containing ..
    output path
        folder path to output folder
    instrument
        type of instrument, "orbitrap" or "qtof"

    Returns
    -------
    bool
        True if function sucesfully ran, False otherwise
    
    """
    #config = "--FormulaSearchSettings.applyFormulaConstraintsToBottomUp=false --IsotopeSettings.filter=true --UseHeuristic.useOnlyHeuristicAboveMz=650 --FormulaSearchDB=plastchem --Timeout.secondsPerTree=0 --FormulaSettings.enforced=H,C,N,O,F,P,I --Timeout.secondsPerInstance=0 --AlgorithmProfile=qtof --SpectralMatchingMassDeviation.allowedPeakDeviation=10.0ppm --AdductSettings.ignoreDetectedAdducts=false --AdductSettings.enforced=, --AdductSettings.prioritizeInputFileAdducts=true --UseHeuristic.useHeuristicAboveMz=300 --IsotopeMs2Settings=IGNORE --MS2MassDeviation.allowedMassDeviation=10.0ppm --SpectralMatchingMassDeviation.allowedPrecursorDeviation=10.0ppm --FormulaSearchSettings.performDeNovoBelowMz=0 --FormulaSearchSettings.applyFormulaConstraintsToDatabaseCandidates=false --EnforceElGordoFormula=true --NumberOfCandidatesPerIonization=1 --FormulaSettings.detectable=B,S,Cl,Se,Br --NumberOfCandidates=10 --AdductSettings.fallback=[[M+H]+,[M+Na]+,[M+K]+] --FormulaSearchSettings.performBottomUpAboveMz=Infinity --RecomputeResults=false spectra-search"
    instrument = kwargs["instrument"] if "instrument" in kwargs else "orbitrap"
    # subprocess_command = f"{sirius_path} --input {input_path} --output {output_path}.sirius config {config} formulas -p {instrument} fingerprints classes structures write-summaries"
    subprocess_command = f"{sirius_path} --input {input_path} --output '{output_path}.sirius' formulas fingerprints classes structures write-summaries"
    run_subprocess(subprocess_command)
    return True

def run_toxtree(input_path: str, output_path: str,toxtree_path: str, module_path: str) -> bool:
    """Runs a Toxtree module on infile from the command line

    Takes a file of SMILES as input to run toxtree from the command line. The specific toxtree module can be selected (see..). 
    Writes the results to fn_out. 

    Parameters
    ----------
    toxtree_path
        file path of toxtree jar file
    module path
        file path of toxtree module that is to be used 
    fn_in
        filename of file with SMILES strings on newlines with a "SMILES" header.
    fn_out
        output fn for toxtree classifications
        
    Returns
    -------
    bool
        True if function sucessfully ran, False otherwise
    
      """
    subprocess_command = f"java -jar {toxtree_path} -n m {module_path} -b {input_path} -i {input_path} -o {output_path}"
    cwd = Path(toxtree_path).parent
    run_subprocess(subprocess_command,cwd=cwd)
    return True
    # def run_toxtree(fn_out,fn_in):
    # toxtree_loc = "/lustre/BIF/nobackup/hendr218/Programs/Toxtree/Toxtree-v3.1.0.1851/Toxtree/Toxtree-3.1.0.1851.jar"
    # module_loc = "toxTree.tree.cramer.CramerRules"
 
    # subprocess.run([subprocess_command], cwd="/lustre/BIF/nobackup/hendr218/Programs/Toxtree/Toxtree-v3.1.0.1851/Toxtree", shell=True, capture_output=True, text=True)
    
    # df_sdf_cramer = rdkit.Chem.PandasTools.LoadSDF(fn_out)
#    # df_cramer = df_sdf_cramer.loc[:,['Cramer rules','SMILES']]
#     return fn_out

def run_classyfire(input_path: str,output_path: str) -> bool:
    """Runs classyfire for strings in infile, writes results to outfile
    
    Uses the classyfire API to classify compounds. First splits input file into batches of 999 compounds to avoid API restrictions.
    Waits between post request to avoid timeouts from the API, so can be slow. 

     Parameters
    ----------
    infile
        filename of file with SMILES strings on newlines
    outfile
        fn for output in .SDF format. 
        
    Returns
    -------
    bool
        True if function sucessfully ran, False otherwise
    
    """
    structures_batches  = file_to_batches(input_path, batchlength = 999)
    api_url: str = "http://classyfire.wishartlab.com" 
    retrieve_ids: list[str] = [post_to_classyfire(batch,api_url) for batch in structures_batches]
    classyfire_results = (get_classyfire_results(api_url,retrieve_id) for retrieve_id in retrieve_ids if check_classyfire_done(api_url, retrieve_id))
    with open (output_path,'a') as file: 
            file.writelines(classyfire_results)
    return True


def run_ms2lda(dataset,params):
    preprocessing_parameters,convergence_parameters,annotation_parameters,model_parameters,train_parameters, \
    fingerprint_parameters,dataset_parameters,n_iterations, n_motifs,motif_parameter = params
    
    MS2LDA.run(dataset, n_motifs=n_motifs, n_iterations=n_iterations,
    dataset_parameters=dataset_parameters,
    train_parameters=train_parameters,
    model_parameters=model_parameters,
    convergence_parameters=convergence_parameters,
    annotation_parameters=annotation_parameters,
    motif_parameter=motif_parameter,
    preprocessing_parameters=preprocessing_parameters,
    fingerprint_parameters=fingerprint_parameters)

    return True
  

def post_to_classyfire(structures_batch: str,base_url: str) -> str:
    """Makes a post request for classification to classyfire for structures (max 999)"""
    post_url: str = f"{base_url}/queries.json"
    data: str = '{"label": "myworkflow", "query_input": "'+structures_batch+'", "query_type": "STRUCTURE"}'
    headers: dict[str,str] = {"Content-Type": "application/json"}
    r = requests.post(post_url,data=data,headers=headers)
    retrieve_id: str = json.loads(r.text)['id']
    return retrieve_id

def check_classyfire_done(base_url: str,retrieve_id: str) -> bool:
    """Checks with API if batch is processed, waits untill it is to return True"""

#Get results from classyfire & create SDF
    status_url: str = f'{base_url}/queries/{retrieve_id}/status.json' 
    can_proceed: bool = False
    while not can_proceed:
        r = requests.get(status_url,headers={"accept": "json"})
        if r.text == "In progress":
            can_proceed = True
        elif (r.status_code != 200) | (r.text != "Done"):
            sleep(61)
        else:
            can_proceed = True
    return True
    # if requests.get(status_url,headers={"accept": "json"}).status_code == "504":
    #     sleep(60)
    # while requests.get(status_ur l,headers={"accept": "json"}).text != "Done":
    #     sleep(60) 
def get_classyfire_results(base_url: str,retrieve_id: str) -> str:
    """Obtains classyfire results in SDF for previsouly send batch via API"""
    retrieve_url: str = f'{base_url}/queries/{retrieve_id}.sdf'
    has_response: bool = False
    while not has_response:
        try:
            r = requests.get(retrieve_url,headers={"accept": "text/plain"},timeout=90)
            if r.status_code == 200:
                contents: str = r.text
                has_response = True
            else:
                sleep(60)
        except requests.exceptions.Timeout: # | (NameError):
            continue
    return contents

def set_ms2lda_params():
    preprocessing_parameters = {
        "min_mz": 0,
        "max_mz": 1000,
        "max_frags": 1000,
        "min_frags": 3,
        "min_intensity": 0.01,
        "max_intensity": 1
    }

    convergence_parameters = {
        "step_size": 50,
        "window_size": 10,
        "threshold": 0.001,
        "type": "perplexity_history"
    }
    annotation_parameters = {
        "criterium": "best", # return cluster with most compounds in it after optimization ("best" also an option)
        "cosine_similarity": 0.70, #0.8 how similar are the spectra compared to motifs in the optimization
        "n_mols_retrieved": 5 # 10 molecules retrieved from database by Spec2Vec
    }
    #Change if needed (cos score)


    model_parameters = {
        "rm_top": 0, #get rid of noisy point
        "min_cf": 0, #intensity within spectrum
        "min_df": 3, #min frequency within document
        "alpha": 0.6, #A higher alpha makes the document preferences "smoother" over topics
        "eta": 0.01, #and a higher eta makes the topic preferences "smoother" over words.
        "seed": 42,
    }

    train_parameters = {
        "parallel": 3,
        "workers": 1, 
    }

    fingerprint_parameters = {
    "fp_type": "rdkit",
    "threshold": 0.8,
    }

    dataset_parameters = {
        "acquisition_type": "DDA",
        "significant_digits": 2, #could be 3
        "charge": 1,
        "name": "pos_all",
        "output_folder": ""
    }
    return preprocessing_parameters,convergence_parameters,annotation_parameters,model_parameters,train_parameters,fingerprint_parameters,dataset_parameters
   
