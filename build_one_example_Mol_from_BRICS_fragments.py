import argparse
import datetime
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import BRICS

def main():
    parser = argparse.ArgumentParser(description='Program to assemble a molecule randomly from the BRICS fragments. Fragments were obtained from promising RPI stabilizer candidates. Outputs a SMI and a PNG.\nExample usage:\n\npython build_one_example_Mol_from_BRICS_fragments.py -o ../_trash', formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-o', '--outputDirectory', default='./', type=str, help='Outputdirectory for molecule containing SMILES file.')
    
    args = parser.parse_args()
    outputDirectory = args.outputDirectory

    # load the fragments from the database and show a random fragment
    with open('RPIS_BRICS_Fragment_DB.smi', 'r') as database:
        brics_fragment_smiles = [line.rstrip('\n') for line in database]
    
    brics_fragments = [Chem.MolFromSmiles(smiles) for smiles in brics_fragment_smiles]

    # Build random molecule from the BRICS fragments (up to 'maxDepth' fragments are used)
    singleSMIoutput = f'{outputDirectory}/brics_{datetime.datetime.now():%Y-%m-%d_%H-%M-%S}.smi'

    for mol in BRICS.BRICSBuild(brics_fragments, maxDepth=3, onlyCompleteMols=True, scrambleReagents=True):
        break

    with open(singleSMIoutput, 'w') as outputFile:
            outputFile.write(Chem.MolToSmiles(mol)+'\n')
    Draw.MolToFile(mol, singleSMIoutput.replace('smi','png'), size=(300, 300))

if __name__ == '__main__':
    main()
