Welcome to MetaboMix, a python-based workflow for metabolomics. It runs and integrates the results of various tools, finally producting a annotated molecular network in GraphML format that can be visualised in tools such as [Cytoscape](cytoscape.org). Each of these tools (with their settings) can be chosen independently for a run of the workflow and new tools can be added easily because of the modular set-up. 
Every run of the workflow is called a Mix, resulting in an annotated molecular network. Every Mix has its own Recipe, a JSON file where you choose which tools you want to use for that Mix, what settings these tools should use and where your input files can be found and output files should go. Tools can be independently ran and integrated. This means you can use a tool from Metabomix, but also add in externally obtained results. See "Quickstart" to create a simple molecular network, "how does it work" to create your own Mix or "adding new tools" to add additional tools. 

Currently supported tools include  [MZMine](https://mzio.io/mzmine-news/) (data pre-processing), [MatchMS](https://github.com/matchms/matchms) (molecular networking), the [SIRIUS/CSI:FingerID/CANOPUS pipeline](https://bio.informatik.uni-jena.de/software/sirius/) (in-silico molecular formula/compound/classification prediction), [ClassyFire](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-016-0174-y) (structure-based classification), [ToxTree](https://toxtree.sourceforge.net/) (structure based toxicity prediction) and the [PlastChem database](https://github.com/PlastChem) (on plastic related  compounds).

MetaboMix was created as a Msc thesis project at the [Van der Hooft Computational Metabolomics Group](https://vdhooftcompmet.github.io/) at the [Bioinformatics department](https://www.wur.nl/en/Research-Results/Chair-groups/Plant-Sciences/Bioinformatics.htm) of [Wageningen University & Research](https://www.wur.nl/en.htm).


# Quickstart
This will create a molecular network via MatchMS from an example MGF (A common format for mass spectrometry data).
1) create conda environment from ENV.yml 
$ conda env -f ENV.yml
2) pip install metabomix in the environment, activate environment
$ pip install metabomix
3) From the example folder: run main.py 
4) In Results > Example_run a graphML file can be found; this is the network and can be visualized in tools such as cytoscape. 
See example_notebook (in the "example" folder) to get a quick walkthrough in jupiter notebook format

# How does it work?
## Recipes
A recipe contains four sections: Paths, run_tools, integrate_tools and individual tool settings. 
Paths defines input files, output locations and tool locations. 
Tools can be independently ran and integrated, using boolean values per tool in run_tools and integrate_tools.
Individual tools settings: defined by having the name of the tools on the top level in the Recipe. Settings then dependent on the tool. 

Paths: 
    Required: 
        internal_settings: config.json for this Mix (see "config" section for more information)
        base_output_folder: path of output location.
    Optional: 
        Name: name for this Mix, otherwise the current datetime is used. 
        input_mgf: if no selected tool (such as MZMine) outputs an mgf, this path is required
        base_network: if no selected tool (such as MatchMS) outputs a molecular network, this path is required
    Tool specific:
        XXX_path: location of a tool, see "tools" section
        Other: see tool description in "tools" for other paths a tool might need
Run/integrate tools: 
    Required: 
        For each available tool:
            tool name > True/False
Individual tools settings: 
    Tool name > tools specific settings

# Config
The config file contains information on the tools that does not have to change from Mix to Mix. This includes requirements ( paths that a tool needs to be used at all), depth (a setting determining the order in which tools are processed) and translations (to unify tool-specific language into desired formats). 

Requirements: to make sure that required files are present. These can be provided in paths in the recipe if a tools is only integrated, or are dynamically generated during a Mix when a tool is ran. Besides "paths", can contain "settings" (tool specific, defined per tool in the recipe) or "optional_paths" (if one of these options is required, but not all of them). 

Depth: to determine the order of tools. 
Tools are run and integrated in the order of depth (so first all tools with depth 1, then 2 etc.). First running happensm then integration (as otherwise the result would not yet be present to integrate). There are 2 important checkpoints:
Depth 3 - graph construction. This means any tool requiring a graph (for example to add annotations to it) must be of depth 3 or more. 
Depth 4 - network construction. This means any tool requiring a network (for example to add annotations to it) must be of depth 4 or more. 
Any tool without a chosen depth in the config is ran as if it were on the final depth. 

# How to add more tools: 
New tools can be added be adding them to the config and then to See the example jupiter notebook for a tutorial. 
Briefly, adding a tool has 4 steps:
1) create code to run the tool 
2) add the tool the the list of available tools on top of the class definition in "mix.py"
3) add tool to the config, define depth and requirements if desired. 
4) add tool running and integration to the recipe in run_tools and integrate_tools
5) add block of code to get/create settings and run the tools in run_tool or integrate_tool in "mix.py"

# Tools
## MZmine
### Running
1) Installation : $ wget https://github.com/mzmine/mzmine/releases/download/v4.7.8/mzmine_Linux_portable_4.7.8.zip
                 $ unzip mzmine_Linux_portable_4.7.8.zip
2) add "mzmine_location" to recipe
    "mzmine_location":"mypath/Programs/MZMine/bin/mzmine" 
3) add "mzmine_userfile_location" to recipe
    obtain via mzmine GUI: users > open users directory 
### Integration
 Currently does not have a universal solution, on a dataset-by-dataset basis. Example for MZmine/sirius to be added soon. 

## Sirius
### Running
1) download sirius, log-in via cml (sirius login -u "username" -p)
2) recipe: add sirius program location to paths as "sirius_path" (to bin/sirius or sirius.exe?)
3) recipe: set run tools > sirius to True

### Integration
1) add sirius results files to paths (not required if running sirius): 
    formula_identifications.tsv as "sirius_tool_output"
    structure_identifications.tsv as "csi:fingerid_output"
    canopus_formula_summary.tsv as "canopus_output"
    NOTE: currently requires all of these to work
2) recipe: set integrate tools > sirius to True

## MS2LDA
### Running:
1) add MS2LDA to conda env, add mass2ql4motifs, remove massql
2) recipe: set run tools > ms2lda to True

### Integrating: 
1) add ms2lda results  to paths as "ms2lda_motifset" (not required if running ms2lda)
2) recipe: set integrate tools > ms2lda to True

## Toxtree 
### Running
1) recipe add "toxtree_path" to paths (to toxtree.jar)
2) recipe: set run tools > toxtree to True
3) recipe: add "module" to "toxtree"
    Currently only supports cramer classification - more modules can be added in translations in the config. 
### Integration
1) recipe: set integrate tools > toxtree to True

## PlastChem
### Integration
1) recipe: run

# Common questions
Q: How do I add a new database?
A: See the Plastchem integration for inspiration on how to do this!

Q: How do I add metadata?
A: Currently on a dataset-by dataset bases. You can always ask me for advice. An example for sirius and MZMine will be added in the future. 

