# RPIS_FragmentationLibraries: Molecular Fragments Useful for Design of Molecular Glues for Protein-RNA Complexes

The database lives at:
- [Zenodo](doi.org/10.5281/zenodo.17527390) (Large Data Files)
- [github](github.com/Foly93/RPIS_FragmentationLibraries) (Executables and README)
- Elsevier (Original Publication)


## Overview

This repository consists of two primary components: **chemical compound libraries in SMILES format**, a universally recognized representation in computational molecular sciences.  

The fragment libraries presented here are part of a **peer-reviewed scientific study (CURRENTLY UNDER REVIEW)**. These fragments were derived through a comprehensive *in silico* workflow designed to identify promising stabilizer candidates for protein–RNA interactions (RPIs).  

The workflow integrates several computational methods, including:

- **Binding pocket detection and evaluation**  
- **Molecular docking**  
- **Molecular dynamics simulations**  
- **Binding free energy estimations**  

This pipeline yielded a set of stabilizer candidates, which were then used to generate fragment libraries through two distinct approaches:

1. **Extended Connectivity Fingerprints (ECFP)** – to extract the most representative chemical features.  
2. **Breaking of Retrosynthetically Interesting Chemical Substructures (BRICS)** – to decompose molecules into synthetically meaningful fragments.  

Both approaches were implemented using Python’s [`rdkit.Chem`](https://www.rdkit.org/docs/source/rdkit.Chem.html) module and can be used for a large spectrum of different applications.

---

### Highlights

- **1,000 most abundant ECFP-derived fragments** for compound database filtering
- **Executable scripts** demonstrating the database filtering workflow  
- **T38DrugDB filtered database**, organized into **9 sub-databases**, containing only molecules with at least `[2, 10]` fragment matches (provided as `.tar.xz` archives for portability)  
- **All 213 BRICS fragments** decomposed from the 96 high-ranking stabilizer ligands identified in our study  
- **Executables for de novo compound generation** from the 213 BRICS fragments, with adjustable `maxDepth` parameters  
- **Comprehensive 4,000,000-compound database**, generated from all possible combinations of the 213 BRICS fragments with `maxDepth = 3`


## Intended Usage

The fragment libraries in this repository are divided into two subsets: **ECFP-derived fragments** and **BRICS-derived fragments**. Although both provide databases of SMILES strings, their intended applications differ.

- **ECFP fragments** can be used to **filter existing compound databases** for molecular patterns enriched in RPI stabilizers identified through our *in silico* workflow. These fragments are particularly useful for identifying chemical motifs associated with stabilizing protein–RNA interactions.

- **BRICS fragments**, on the other hand, also encode **combinatorial information**. In this method, chemical bonds of RPI stabilizers were broken in a retrosynthetically meaningful way using the `rdkit.Chem.BRICS` module. The resulting fragments can be **recombined** following the same retrosynthetic rules, enabling the **construction of new compounds** enriched with RPI-stabilizing features. These fragments can also be used for database filtering, but since they lack fragment importance scores, filtering must be done naively using all fragments.

While these two applications, database filtering (ECFP) and de novo compound generation (BRICS), represent the primary intended uses, they are not the only possibilities. Generative machine learning techniques, for instance, could utilize these fragment libraries to design new potential binders, particularly leveraging the ECFP dataset.

In general, the **scope of applications** for these fragment libraries is broad and limited only by the creativity of the user. The included executable scripts demonstrate the core workflows for **database filtering (ECFP)** and **de novo compound assembly (BRICS)**.

### ECFP – Database Filtering

Within the `ECFP/` directory, you will find the executable script `RefilterDatabaseWithECFP.py` and the example dataset `sampleDB.csv`. The fragment data required for execution is provided in the file `MMGBSA_ChemicalAnalysisFragments_cutoff10_ranked_datatable_3orMore.csv`. This file must be present to run `RefilterDatabaseWithECFP.py` without errors. After installing all required dependencies (listed at the beginning of `RefilterDatabaseWithECFP.py`), the script can be executed to produce nine `.csv` files named in the format: ECFP_fragmentFilteredLibrary_#N#FragsOrMore.csv. These files serve as example outputs. 
The actual results of the filtering process are available in the [Zenodo upload](doi.org/10.5281/zenodo.17527390) associated with this repository. That dataset contains the results of applying the filtering procedure to the [T38DrugDB](doi.org/10.5281/zenodo.17410437), which includes approximately **34 million drug-like compounds** published elsewhere. The filtered sub-databases are provided as compressed `.tar.xz` archives for portability. These databases can be used for **targeted virtual screening** of **RPI-stabilizing drug candidates**. This example workflow can be executed as follows:
```bash
# clone this repository
git clone git@github.com:Foly93/RPIS_FragmentationLibraries.git ./

# switch to ECFP directory
cd RPIS_FragmentationLibraries/ECFP

# INSTALL REQUIRED PYTHON PACKAGES FROM YOUR FAVOURITE PACKAGE MANAGER e.g. micromamba
micromamba activate your_env_name_goes_here

# execute the Python program
python RefilterDatabaseWithECFP.py 

# Check if the expected output was generated (assuming ubuntu, mac or power shell)
ls -rtal ECFP_fragmentFilteredLibrary_*FragsOrMore.csv
```

### BRICS – De Novo Compound Assembly

The ```BRICS/``` directory contains several executables that serve different purposes, as well as the BRICS fragment library ```RPIS_BRICS_Fragment_DB.smi```. This ```.smi``` file contains **213 BRICS fragments**, which can be combined to generate novel molecules using the ```rdkit.Chem.BRICS``` module. This functionality is demonstrated across three executables:
- ```BRICS_fragment_database_interactive.ipynb```
- ```build_one_example_Mol_from_BRICS_fragments.py```
- ```generate_N_Mols_from_BRICS_fragments.py```

The Jupyter notebook ```BRICS_fragment_database_interactive.ipynb``` requires a Python environment with ```jupyter notebook``` installed. It provides a visual and interactive overview of the fragment set, demonstrating how to:
- Display the BRICS fragments  
- Assemble random molecules from the fragment library  
- Export the resulting molecules as SMILES strings  

The assembly process uses the RDKit function ```BRICS.BRICSBuild```. Its parameters are highly sensitive—particularly the `maxDepth` option, which controls the maximum number of fragments combined into a single molecule. With 213 fragments available, the total combinatorial space is on the order of **4 × 10¹¹** possible assemblies, making exhaustive enumeration computationally infeasible. 
The script `build_one_example_Mol_from_BRICS_fragments.py` replicates the notebook’s core functionality in a standalone Python executable. It generates a single example molecule, saving the output to a file identified by its timestamp.
The final executable, `generate_N_Mols_from_BRICS_fragments.py`, creates a specified number of BRICS-assembled molecules and saves them into a timestamped `.smi` file. While it can be run without command-line arguments, optional parameters are available and can be displayed using the `-h` flag. The generated `.smi` files can be directly used for **virtual screening** or **binding affinity prediction** in **drug discovery** workflows. This showcase can be executed as follows:
```bash
# clone this repository
# not necessary if already done
git clone git@github.com:Foly93/RPIS_FragmentationLibraries.git ./

# switch to BRICS directory
cd RPIS_FragmentationLibraries/BRICS

# INSTALL REQUIRED PACKAGES FROM YOUR FAVOURITE PACKAGE MANAGER e.g. micromamba
micromamba activate your_env_name_goes_here

# open the jupyter notebook in your browser and have a good look at the functionality
jupyter notebook BRICS_fragment_database_interactive.ipynb

# execute the one-example-python-program
python build_one_example_Mol_from_BRICS_fragments.py 

# execute the batch creation python program
python generate_N_Mols_from_BRICS_fragments.py 

# display the options of the batch creation python program
python generate_N_Mols_from_BRICS_fragments.py -h 

# run the python program with custom flags
python generate_N_Mols_from_BRICS_fragments.py \
    --maxDepth 3 \
    --numMold 10 \
    --scrambleReagents True \
    --outputDirectory ../trash
```
Finally, the [Zenodo upload](doi.org/10.5281/zenodo.17527390) also contains ```BRICS_DB_BuiltMaxDepth_3.txt``` a data base that contains all possible combinations for maxDepth set to 3. This file contains 3,878,955 compound SMILES strings which is slightly less that the theoretically possible $` 213\cdot 213\cdot 213 = 9,663,597 `$ which results from incompatibilities between some BRICS fragments.

## Directory Structure 
```
RPIS_FragmentationLibraries/
├── README.md                          # THIS file
├── RPISFragmentsZenodoUpload.json     # upload json file for Zenodo (just for completeness)
│
├── ECFP/
│   ├── MMGBSA_ChemicalAnalysisFragments_cutoff10_ranked_datatable_3orMore.csv   # ECFP fragments csv file 
│   ├── RefilterDatabaseWithECFP.py      # executable for show casing data base filtering
│   ├── sampleDB.csv                     # auxiliary sample data base filtered by RefilterDatabaseWithECFP.py 
│   └── ECFP_*FragsOrMore.csv.tar.xz     # large filtered sub-databases only available from ZENODO
│
└── BRICS/
    ├── RPIS_BRICS_Fragment_DB.smi                      # BRICS fragments derived according to our publication
    ├── BRICS_fragment_database_interactive.ipynb       # jupyter notebook for exploring BRICS Assmebly interactively
    ├── build_one_example_Mol_from_BRICS_fragments.py   # program to build one example BRICS assembly
    ├── generate_N_Mols_from_BRICS_fragments.py         # program to build batches of BRICS assembled molecules, e.g. for screening 
    └── BRICS_DB_BuiltMaxDepth_3.txt                    # Data base containing SMILES of all available BRICS assemblies with maxDepth set to 3
```

## Citation
If you use these Fragment libraries, please cite:
```
[Add appropriate citation information here once available]
```

## License
This work is licensed under a Creative Commons Attribution 4.0 International License. See creativecommons.org/licenses/
by/4.0/ for further information.

## Contact
For questions, issues, or contributions:
- luis.vollmers@tum.de
- zacharias@tum.de
- %%%publication DOI once available
