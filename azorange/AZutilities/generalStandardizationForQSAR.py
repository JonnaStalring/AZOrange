import os
import string
import random

from rdkit import Chem
from rdkit.Chem import Draw

from molvs import Standardizer
from molvs import tautomer
from molvs import charge
from molvs import fragment
from molvs import normalize
from molvs import validate_smiles

print "Need HG RDKit environment and the cheminfo module"

IMAGESDIR = "./StdzedImages"
SCRATCHDIR = "/disk1/jonsta/temp/standScratch/"


def reportDuplicates(smilesDict):
    """
    Input: Dictionary; name: smiles
    Output: List of dict with names of duplicates and all names with the canonical smiles
    The input smiles are assumed to be canonical such that the identification of duplicates can be performed by a simple 
    string comparison of smiles strings. 
    """   

    duplicates = {}
    for name, smiles in smilesDict.iteritems():
        for inname, insmiles in smilesDict.iteritems():
            if smiles == insmiles and name != inname:    # There is an identical smiles
                if duplicates.has_key(smiles):           # A duplicate with this smiles was found before
                    if name not in duplicates[smiles]:   # The name has not been added to this key
                        duplicates[smiles].append(name)
                else:                                    # The duplicate was not found before
                    duplicates[smiles] = [name, inname]
    return duplicates


def createSDF(inMol):
    sdfFile = SCRATCHDIR+'standardize'+str(random.random())
    w = Chem.SDWriter(sdfFile+".sdf")
    w.write(inMol)
    return sdfFile


def standardizeCA(inMol):

    sdfFile = createSDF(inMol)
    # Standardizer
    sdfOutFile = sdfFile+"Out.sdf"
    os.system("standardize -c StandardizerConfigForDB.xml "+sdfFile+".sdf -f sdf -o "+sdfOutFile)
    outMol = Chem.SDMolSupplier(sdfOutFile)[0]

    return outMol


def standardizeMolVS(inMol):

    t = tautomer.TautomerCanonicalizer()
    outMol = t.canonicalize(inMol)
    c = charge.Uncharger()
    outMol = c.uncharge(outMol)

    return outMol


def standardize(mol, name, texFile, verbose):
    """
    General method for standardization of smiles for QSAR modeling considering:
    1) Removal of inorganic compounds
    2) Removal of mixture
    3) Removal of counter ions
    4) Neutralization of organic ions if possible
    5) Canonical representation of mesomers (carboxylic acids, nitro groups)
    6) Canonical representation of tautomers
    7) Canonical representation of aromaticity
    8) Report duplicates
    """

    outMol = standardizeCA(mol)
    stdMol = standardizeMolVS(outMol)
    smiles = Chem.MolToSmiles(stdMol)

    if verbose:
        # For drawing structures
        image = Draw.MolsToGridImage([mol, stdMol], molsPerRow=2, subImgSize=(200,200))
        fileName = "transformation"+str(random.randint(1, 1000000000))+".png"
        image.save(os.path.join(IMAGESDIR, fileName))
        fid = open(texFile, "a")
        fid.write(string.replace(name, "_", "")+" & \\includegraphics[scale=0.8]{"+fileName+"} \\\\\n")
        fid.write("\\hline\n")
        fid.close()

    return smiles 


def createReportFiles(fileName):

    # Output file with standardized structures

    # Put output images into subdirectory of pwd
    if not os.path.isdir(IMAGESDIR):
        os.system("mkdir "+IMAGESDIR)
    else:
        os.system("rm -r "+IMAGESDIR)
        os.system("mkdir "+IMAGESDIR)

    # Latex file for compiling structure depictions
    texFile = os.path.join(IMAGESDIR, "compilation.tex")
    fid = open(texFile, "w")
    fid.write("\\documentclass[10pt, a4paper]{article}\n")
    fid.write("\\usepackage{longtable}\n")
    fid.write("\\usepackage{graphicx}\n")
    fid.write("\\begin{document}\n")
    fid.write("\\begin{longtable}{|l|l|}\n")
    fid.write("\\hline\n")
    fid.write("Name & In and Out Smiles \\\\\n")
    fid.write("\\hline\n")
    fid.close()



if __name__ == "__main__":

    # Location of input txt file
    fileName = os.sys.argv[1]
    print "Standardizing the file ", fileName

    verbose = False
    print "Verbose mode is ", verbose

    if verbose:
        createReportFiles(fileName)
    else:
        texFile = "None"

    # Read and write standardized structures
    outFile, ext = os.path.splitext(fileName)
    outFileName = outFile+"Stdzd.txt"
    outFid = open(outFileName, "w")
    outFid.write("Smiles\tName\tClass\tpIC50\n")
    outFid.close()
    fid = open(fileName)
    logFile = outFile+"log.txt"
    logFid = open(logFile, "w")
    logFid.close()
    header = fid.readline()
    smilesDict = {}
    print "Might need to change in accordance with file format "
    print "**********Please check file format*********"
    idx = 0
    for line in fid:
        idx = idx + 1
        lineList = string.split(line, "\t") 
        smiles = lineList[2]
        name = lineList[3]
        className = string.strip(lineList[4])
        pIC50 = string.strip(lineList[5])
        logFid = open(logFile, "a")
        #logFid.write("........... Processing "+smiles+" .............\n")
        info = validate_smiles(smiles)
        if info:
            logFid.write(str(info)+" "+smiles+" "+name+"\n")
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            logFid.write("Could not be transformed to rdkit object "+ smiles+" "+name+"\n")
            stdzdSmiles = None
        else:
            stdzdSmiles = standardize(mol, name, texFile, verbose)
        if stdzdSmiles:
            smilesDict[name] = stdzdSmiles
            outFid = open(outFileName, "a")
            outFid.write(str(stdzdSmiles)+"\t"+name+"\t"+className+"\t"+pIC50+"\n")
            outFid.close()
        logFid.close()
    fid.close()
    logFid.close()

    # Report duplicates
    duplicates = reportDuplicates(smilesDict)

    print "Duplicates and the number of duplicates"
    print duplicates
    print len(duplicates)
    
