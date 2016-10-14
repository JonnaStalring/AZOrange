from rdkit.Chem import Draw
from rdkit import Chem

def writeToFile(mols, names, comments, fileName):

    fid = open("CheckStandardizationTmp.txt", "a")
    idx = 0
    for mol in mols:
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        name = names[idx]
        comment = comments[idx]
        idx = idx + 1
        print smiles, name, comment
        fid.write(smiles+"\t"+name+"\t"+comment+"\n")
    fid.close()
    image = Draw.MolsToGridImage(mols, molsPerRow=6, subImgSize=(600,600), legends=names)
    image.save(fileName+".png")


fid = open("CheckStandardizationTmp.txt", "w")
fid.close()

#mols = Chem.SDMolSupplier("phosphate.sdf")
#names = ["Phosphate"]
#comments = ["Inorganic"]
#writeToFile(mols, names, comments, "phosphate")
#
#mols = Chem.SDMolSupplier("nitro.sdf")
#names = ["NitroMesomer1", "NitroMesomer2", "NitroMesomer3"]
#comments = ["Mesomer", "Mesomer", "Mesomer"]
#writeToFile(mols, names, comments, "nitro")
#
#mols = Chem.SDMolSupplier("triazine.sdf")
#names = ["TriazineTau1", "TriazineTau2"]
#comments = ["Tautomer", "Tautomer"]
#writeToFile(mols, names, comments, "triazine")
#
#mols = Chem.SDMolSupplier("propenal.sdf")
#names = ["PropenalMesomer1", "PropenalMesomer2"]
#comments = ["Mesomer", "Mesomer"]
#writeToFile(mols, names, comments, "propenal")
#
#mols = Chem.SDMolSupplier("propaneAcid.sdf")
#names = ["PropaneAcidTau1", "PropaneAcidTau2"]
#comments = ["Tautomer", "Tautomer"]
#writeToFile(mols, names, comments, "PropaneAcid")
#
#mols = Chem.SDMolSupplier("Propenol.sdf")
#names = ["PropenolMesomer1", "PropenolMesomer2"]
#comments = ["Tautomer", "Tautomer"]
#writeToFile(mols, names, comments, "Propenol")

mols = Chem.SDMolSupplier("aminopyrimidine.sdf")
names = ["AminopyrimidineTau1", "AminopyrimidineTau2"]
comments = ["Tautomer", "Tautomer"]
writeToFile(mols, names, comments, "Aminopyrimidine")
