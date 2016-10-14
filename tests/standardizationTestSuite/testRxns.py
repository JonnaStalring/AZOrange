from rdkit import Chem
from rdkit.Chem import AllChem

inMols = Chem.SDMolSupplier("diazo.sdf")
inMols = Chem.SDMolSupplier("phosphine.sdf")
inMols = Chem.SDMolSupplier("nitro.sdf")

#rxn = AllChem.ReactionFromSmarts('[C:1]=[O,N:2]>>[C:1][*:2]')
#rxn = AllChem.ReactionFromSmarts('[S+2:1]([O-:2])([O-:3])>>[S+0:1](=[O-0:2])(=[O-0:3])')
#rxn = AllChem.ReactionFromSmarts('[N:1](=[O:2])=[O:3]>>[*+1:1]([*-1:2])=[*:3]') # Transformed in RDK

# Diazo
#rxn = AllChem.ReactionFromSmarts('[C-:1][N+:2]#[N:3]>>[C+0:1]=[N+:2]=[N-:3]')
# Phosphine
#rxn = AllChem.ReactionFromSmarts('[PH3+:1]>>[PH2+0:1]')
# Nitro
rxn = AllChem.ReactionFromSmarts('[N+:1](=[O:2])([OH:3])>>[N+1:1]([O-1:2])=[O:3]')
rxn = AllChem.ReactionFromSmarts('[N:1]([OH:2])([OH:3])>>[N+1:1]([O-1:2])=[O:3]')


#smilesList = ['CC=O', 'CC=N', 'CCO']
#smilesList = ['[O-]S[O-]', 'C1CCS(=O)(=O)CC1', '[S+2:1]([O-:2])([O-:3])', 'O=N(=O)c1ccccc1', '']
#smilesList = [Chem.MolToSmiles(inMol[1])]

for mol in inMols:
    #mol = Chem.MolFromSmiles(inSmiles)
    prod = rxn.RunReactants((mol,))
    print "inSmiles ", Chem.MolToSmiles(mol)
    try:
        print "Prod ", Chem.MolToSmiles(prod[0][0]) 
    except:
        print "No product"
