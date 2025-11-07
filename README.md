# RPIS_FragmentationLibraries: Molecular Fragments Useful for Design of Molecular Glues for Protein-RNA Complexes

The database lives at:
- [Zenodo](doi.org/10.5281/zenodo.17527390) (Large Data Files)
- [github](github.com/Foly93/RPIS_FragmentationLibraries) (Executables and README)
- Elsevier (Original Publication)

## Overview
The provided fragment libraries contain molecular fragments in SMILES format for database filtering and de novo drug design. Extended connectivity fingerprints (ECFP) and Breaking of retrosynthatically interesting chemical substructures (BRICS) was used to compile the database.
- 1000 most abundant ECFP derived fragments for database filtering
- Executables showcasing the database filtering process
- T38DrugDB filtered with 43 highest ranking fragments into 9 sub-databases with only molecules containting atleast {2,10} Fragments. Saved as tar.xz for portability
- All 213 BRICS fragments that were decomposed from the 96 elevated stabilizer ligands found by our study
- Executables to create random de novo compounds from the 213 BRICS fragments with variable maxDepth values
- A 4,000,000 compounds database generated from all possible combinations of the 213 BRICS fragments with maxDepth set to 3

## Usage
The fragment libraries presented in this repository can be split into two subets. Firstly, a subset of fragments derived from ECFP and fragments derived from BRICS method. Although both approaches yield databases of SMILES strings the applications the two fragment libraries can be used for differ.
The ECFP fragment SMILES can be used to filter existing data bases for molecular patterns that are enriched in protein-RNA interaction (RPI) stabilizers found in our in silico research.
The BRICS fragments contain information about combinatorics as well. In this method, the RPI stabilizers' chemical bonds were broken in a retrosynthatically meaningful way via the ```rdkit.Chem.BRICS``` module. The bonds can be recombined in a way that can also be deemed retrosynthatically. Thus the BRICS fragments facilitate the assembly of compounds from fragments enriched with RPI stabilizing patterns.
These two applications might not be the only possibilities. Generative machine learning techniques used molecular fragments to create binders, an application the ECFP fragment library is well suited for. Generally. there is a broad spectrum of applications that these fragment libraries could be used for, limited only by the creativity of the respective user.The executables that are included in the repository showcase only the intended usage of data base filtering (ECFP) and de novo compound assembly (BRICS).

### ECFP - Data base filtering
In the ```ECFP/``` directory the executable ```RefilterDatabaseWithECFP.py``` and ```sampleDB.csv```can be found. The file that contains the fragments is called ```MMGBSA_ChemicalAnalysisFragments_cutoff10_ranked_datatable_3orMore.csv``` which is required for executing ```RefilterDatabaseWithECFP``` without error. After installing all required packages listed in the beginning of ```RefilterDatabaseWithECFP.py``` the program can be executed and should yield nine ```.csv``` files, called ```ECFP_fragmentFilteredLibrary_#N#FragsOrMore.csv```. These are only dummy files. The [Zenodo upload](doi.org/10.5281/zenodo.17527390) of this repository contains the results of the filtering applied to the [T38DrugDB](doi.org/10.5281/zenodo.17410437) which contains 34,000,000 million drug like compounds and is published elsewhere. The filtered sub-databases are saved as ```*.tar.xz``` files for portability. These database files can be used for targeted drug screening for RPI stabilizing drug candidates. This showcase can be executed in the following way:
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

### BRICS - De Novo Compound Assmebly
The ```BRICS/``` directory contains several executables serving different purposes and the BRICS fragment smiles ```RPIS_BRICS_Fragment_DB.smi```. The ```*.smi``` file contains 213 BRICS fragments that can be combined with each other to yield novel molecules via the ```rdkit.Chem.BRICS``` module. This utility is showcased in all three executables: ```BRICS_fragment_database_interactive.ipynb```, ```build_one_example_Mol_from_BRICS_fragments.py``` and ```generate_N_Mols_from_BRICS_fragments.py```. The file ```BRICS_fragment_database_interactive.ipynb``` requires a python environment with ```jupyter notebook``` installed. The notebook can then be opened and shows what the fragments look like, how to generate a random molecule assembled from the BRICS fragments and how to save it as a SMILES string. For the assembly the RDKit function ```BRICS.BRICSBuild``` is used. The options to this function are very sensitive. The option ```maxDepth``` is notworthy, since this option decides how many fragments are combined maximally. For 213 fragments the resulting combinatorics are in the order of 4e11 and should not be conducted in its entirety since this computation is likely to exceed a given timeframe.
Next ```build_one_example_Mol_from_BRICS_fragments.py``` does exactly what it says and orients its functionality closely to the jupyter notebook. The output file generated from this executable will be identified by its timestamp.
The last executable ```generate_N_Mols_from_BRICS_fragments.py``` generates a specified number of BRICS assemblies into one timestamp identified ```*.smi``` file. Techincally this program can be executed without command line arguments, but there are options that can be displayed via the ```-h``` flag. The resulting ```*.smi``` file can be used for virtual screening or binding affinity predictions in drug discovery endeavours.
The showcase can be executed in the following way:
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
If you use T38DrugDB in your research, please cite:
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
