import argparse
import datetime
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import BRICS

def main():
    parser = argparse.ArgumentParser(description='Program to fetch N randomly assembled compounds from the BRICS database. Fragments were obtained from promising RPI stabilizer candidates.\nExample usage:\n\npython generate_N_Mols_from_BRICS_fragments.py -N 10 -s True -o ../_trash -M 3', formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-N', '--numMols', default=100, type=int, help='number of molecules to generate')
    parser.add_argument('-M', '--maxDepth', default=5, type=int, help='maximum number of fragments to use for each generated molecule')
    parser.add_argument('-s', '--scrambleReagents', default=True, type=bool, help='If set to false, the order for molecule generation is always the same. Only useful for N > 1,000,000. E.g., to built libraries reliably.')
    parser.add_argument('-o', '--outputDirectory', default='./', type=str, help='Outputdirectory for molecule containing SMILES file.')
    
    args = parser.parse_args()
    
    N = args.numMols
    maxDepth = args.maxDepth
    scrambleReagents = args.scrambleReagents
    outputDirectory = args.outputDirectory

    # load the fragments from the database and show a random fragment
    with open('RPIS_BRICS_Fragment_DB.smi', 'r') as database:
        brics_fragment_smiles = [line.rstrip('\n') for line in database]
    
    brics_fragments = [Chem.MolFromSmiles(smiles) for smiles in brics_fragment_smiles]

    # Build random molecule from the BRICS fragments (up to 'maxDepth' fragments are used)
    offSwitch = 0
    singleSMIoutput = f'{outputDirectory}/brics_{datetime.datetime.now():%Y-%m-%d_%H-%M-%S}.smi'

    with open(singleSMIoutput, 'w') as outputFile:
        for mol in tqdm(BRICS.BRICSBuild(brics_fragments,
                                    maxDepth=maxDepth, 
                                    onlyCompleteMols=True,
                                    uniquify=True,
                                    scrambleReagents=scrambleReagents,
                                    seeds=None,), total=N-1):
            offSwitch += 1
            outputFile.write(Chem.MolToSmiles(mol)+'\n')

            if offSwitch >= N:
                break


if __name__ == '__main__':
    main()
