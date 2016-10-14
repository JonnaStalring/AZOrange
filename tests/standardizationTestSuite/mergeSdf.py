from rdkit import Chem
import os

outFile = "AllTestCmpds.txt"
fid = open(outFile, "w")

files = os.listdir(".")
idx = 1
for file in files:
    fileName, ext = os.path.splitext(file)
    if ext == ".sdf":
        print fileName
        inMols = Chem.SDMolSupplier(fileName+".sdf")
        for mol in inMols:
            smiles = Chem.MolToSmiles(mol)
            fid.write(smiles+"\n")
fid.close()
