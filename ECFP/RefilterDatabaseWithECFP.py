import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem

def main():
    T38DrugDB_smiles = pd.read_csv("../../03-localCompoundDB/T38DrugDB.csv").iloc[:,0]
    fragment_smiles = list(pd.read_csv("../_output/DAT/MMGBSA_ChemicalAnalysisFragments_cutoff10_ranked_datatable_3orMore.csv",index_col=0)['SMILES'])
    fragments = [Chem.MolFromSmiles(smi,sanitize=False) for smi in fragment_smiles[:43] if Chem.MolFromSmiles(smi,sanitize=False)]
    fragmentHasFragmentsList = create_super_sub_fragment_matrix(fragments)
    
    mol_arom_super_frag_ids = []
    mol_arom_super_frag_smis = []

    for smi in tqdm(T38DrugDB_smiles):
        moleculeAndFragments = match_superfragments_with_correct_aromaticity(smi, fragments, fragmentHasFragmentsList)

        if moleculeAndFragments:
            mol_arom_super_frag_ids.append(moleculeAndFragments[0])
            mol_arom_super_frag_smis.append(moleculeAndFragments[1])

    mol_arom_super_frag_smis = np.array(mol_arom_super_frag_smis)
    molFragCounts = np.zeros(len(mol_arom_super_frag_ids))
    for i in range(len(mol_arom_super_frag_ids)):
        molFragCounts[i] = mol_arom_super_frag_ids[i].size

    for fragment_min in [2,3,4,5,6,7,8,9,10]:
        save_compound_library(mol_arom_super_frag_smis, mol_arom_super_frag_ids, molFragCounts, fragments, fragment_min=fragment_min)


def save_compound_library(mol_arom_super_frag_smis, mol_arom_super_frag_ids, molFragCounts, fragments, fragment_min=2, fileName=None):
    if not fileName:
        fileName = f'../_output/ECFP_fragmentFilteredLibrary_{fragment_min}FragsOrMore.csv'
    mol_smiles = mol_arom_super_frag_smis[np.where(molFragCounts >= fragment_min)[0]]
    frag_smiles = [[Chem.MolToSmiles(mol) for mol in np.array(fragments)[frag_id]] for frag_id in mol_arom_super_frag_ids if frag_id.size >= fragment_min]
    print(f'saving SMILES to files for Compounds with at least {fragment_min} fragments')
    pd.DataFrame(data={'SMILES_MOL': mol_smiles,
                       'SMILES_FRAGS': frag_smiles}).to_csv(fileName, index=None)


def match_superfragments_with_correct_aromaticity(smi, fragments, fragmentHasFragmentsList):
    # initialize molecule as RDKit object, and get matching fragment indices
    mol = Chem.MolFromSmiles(smi)
    hasFragmentBools = [mol.HasSubstructMatch(frag) for frag in fragments]
    hasFragmentIds = np.where(np.array(hasFragmentBools) == True)[0]
    # skip molecule if no fragments match
    if hasFragmentIds.size <= 0:
        return None
    
    correctAromaticFragIds = []
    for hasFragmentIdx in hasFragmentIds:
        mol_aromaticity = []
        frag_aromaticity = []
        # iterate over fragment atom ids and matching molecule atom ids to ensure aromaticity match
        for frag_atom_idx, mol_atom_idx in enumerate(mol.GetSubstructMatch(fragments[hasFragmentIdx])):
            # save aromaticities of molecule atoms
            mol_aromaticity.append(mol.GetAtomWithIdx(mol_atom_idx).GetIsAromatic())
            # save aromaticities of fragment atoms
            frag_aromaticity.append(fragments[hasFragmentIdx].GetAtomWithIdx(frag_atom_idx).GetIsAromatic())
        # fragment only truly matches if aromaticity matches perfectly (not automatically ensured by HasSubstructMatch(frag) s.a.)
        if all(np.array(mol_aromaticity) == np.array(frag_aromaticity)): # yes I have checcked that the Atom Ids matches on the two lists
            correctAromaticFragIds.append(hasFragmentIdx)
    hasCorrectAromFragIds = np.array(correctAromaticFragIds)
    
    # skip molecule if no fragment is left after aromaticity check
    if hasCorrectAromFragIds.size <= 0:
        return None

    # Crate Boolean Array with True for Fragment indices where applies
    hasFragmentArray = np.zeros(len(hasFragmentBools), dtype=bool)
    hasFragmentArray[hasCorrectAromFragIds] = True
    hasSuperFragmentBools = list(hasFragmentArray).copy()
    
    for CorrectAromFragIdx in hasCorrectAromFragIds:
        # Fetch all subfragments of each Fragment found in the molecule (with correct aromaticity)
        FoundSubfragmentsIds = np.where(fragmentHasFragmentsList[CorrectAromFragIdx] == True)[0]
        # Ensure that all subfragments of a matched molecule fragment are set to false
        for subfragmentId in FoundSubfragmentsIds:
            if subfragmentId != CorrectAromFragIdx:
                hasSuperFragmentBools[subfragmentId] = False
    
    # if no fragments are left after subfragment cancellation skip the molecule (I think the code is wrong, but it shouldnt happen either way)
    if np.array(hasSuperFragmentBools).size <= 0:
        return None
    mol_arom_super_frag_id = np.where(np.array(hasSuperFragmentBools) == True)[0]
    mol_arom_super_frag_smi = smi

    return mol_arom_super_frag_id, mol_arom_super_frag_smi


def create_super_sub_fragment_matrix(fragments):
    # initialize the SuperFragment / Subfragment Matrix
    fragmentHasFragmentsList = []
    for super_frag in fragments:
        preliminary_frags = np.array([super_frag.HasSubstructMatch(sub_frag) for sub_frag in fragments])
        preliminary_idcs = np.where(preliminary_frags == True)[0]
        correctAromaticFragIds = []
    
        for preliminary_idx in preliminary_idcs:
            super_frag_aromaticity = []
            sub_frag_aromaticity = []
            
            for sub_frag_atom_idx, super_frag_atom_idx in enumerate(super_frag.GetSubstructMatch(fragments[preliminary_idx])):

                super_frag_aromaticity.append(super_frag.GetAtomWithIdx(super_frag_atom_idx).GetIsAromatic())
                sub_frag_aromaticity.append(fragments[preliminary_idx].GetAtomWithIdx(sub_frag_atom_idx).GetIsAromatic())
    
            if all(np.array(super_frag_aromaticity) == np.array(sub_frag_aromaticity)):
                correctAromaticFragIds.append(preliminary_idx)
            else:
                print(f'non-matching aromaticity detected for {Chem.MolToSmiles(fragments[preliminary_idx])} and {Chem.MolToSmiles(super_frag)}')
        substructure_frag_idcs = np.array(correctAromaticFragIds)
        substructure_frags = np.zeros(preliminary_frags.size, dtype=bool)
        substructure_frags[correctAromaticFragIds] = True
        fragmentHasFragmentsList.append(substructure_frags)
    
    fragmentHasFragmentsList = np.array(fragmentHasFragmentsList)
    return fragmentHasFragmentsList

if __name__ == '__main__':
    main()
