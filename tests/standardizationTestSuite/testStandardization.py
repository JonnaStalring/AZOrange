import os
import time
import random
import string

from rdkit.Chem import Draw
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

from molvs import Standardizer
from molvs import tautomer
from molvs import charge
from molvs import fragment
from molvs import normalize
#import normalize


SCRATCHDIR = "/disk1/jonsta/temp/standScratch/"


def createSDF(inMol):
    sdfFile = SCRATCHDIR+'standardize'+str(random.random())
    w = Chem.SDWriter(sdfFile+".sdf")
    w.write(inMol)
    return sdfFile 

def standardizeCA(inMol):
    sdfFile = createSDF(inMol)

    # Standardizer
    sdfOutFile = sdfFile+"Out.sdf"
    #os.system("standardize -c 'removesolvents..stripsalts..removefragment:method=keeplargest..neutralize..aromatize..tautomerize..mesomerize' "+sdfFile+".sdf -f sdf -o "+sdfOutFile)
    os.system("standardize -c StandardizerConfigForDB.xml "+sdfFile+".sdf -f sdf -o "+sdfOutFile)

    # Canonical tautomer
    #sdfOutFile2 = sdfFile+"Out2.sdf"
    #time.sleep(5)
    #os.system("cxcalc canonicaltautomer -f sdf "+sdfOutFile+" > "+sdfOutFile2)
    #print sdfOutFile2
    #outMol = Chem.SDMolSupplier(sdfOutFile2)[0]

    outMol = Chem.SDMolSupplier(sdfOutFile)[0]

    #inchi = Chem.inchi.MolToInchi(outMol)
    #print inchi
    #outMol = Chem.inchi.MolFromInchi(inchi)
    return outMol

def standardizeMolVS(inMol):
    #f = fragment.LargestFragmentChooser()
    #outMol = f.choose(inMol)
    #c = charge.Uncharger()
    #outMol = c.uncharge(outMol)
    #s = Standardizer()
    #outMol = s.standardize(outMol)
    #n = normalize.Normalizer()
    #outMol = n.normalize(outMol)

    t = tautomer.TautomerCanonicalizer()
    outMol = t.canonicalize(inMol)
    c = charge.Uncharger()
    outMol = c.uncharge(outMol)

    # Transform with Inchi
    #print "inMol"
    #print Chem.MolToSmiles(inMol)
    #inchi = Chem.inchi.MolToInchi(inMol)
    #print inchi
    #print "outMol"
    #print Chem.MolToSmiles(outMol)
    #inchi = Chem.inchi.MolToInchi(outMol)
    #print inchi
    #outMol = Chem.inchi.MolFromInchi(inchi)

    return outMol


def getTanDist(outMols):
    """Get tan dist between all pairs in outMols """
    tanDists = []
    tanDistsMorgan = []
    fps = [FingerprintMols.FingerprintMol(x) for x in outMols]
    for outIdx in range(len(outMols)): 
        for inIdx in range(outIdx + 1, len(outMols)): 
            print outIdx, inIdx
            tanDist = DataStructs.FingerprintSimilarity(fps[outIdx],fps[inIdx])
            fpsM1 = AllChem.GetMorganFingerprint(outMols[outIdx],2)
            fpsM2 = AllChem.GetMorganFingerprint(outMols[inIdx],2)
            #tanDistM = DataStructs.TanimotoSimilarity(fpsM1, fpsM2)
            tanDistM = DataStructs.DiceSimilarity(fpsM1, fpsM2)
            tanDists.append(round(tanDist,2)) 
            tanDistsMorgan.append(round(tanDistM,2)) 
    return tanDists, tanDistsMorgan        


def getDescDiff(outMols):

    descDiffList = []
    nms=[x[0] for x in Descriptors._descList]
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(nms)
    for outIdx in range(len(outMols)): 
        for inIdx in range(outIdx + 1, len(outMols)): 
             descrs1 = calc.CalcDescriptors(outMols[outIdx])
             descrs2 = calc.CalcDescriptors(outMols[inIdx])
             descDiff = []
             for idx in range(len(descrs1)):
                 if descrs1[idx] != descrs2[idx]:
                     print "Differing descriptor ", nms[idx]
                     print descrs1[idx], descrs2[idx]
                     descDiff.append(string.replace(nms[idx], "_", ""))
             descDiffList.append(len(descDiff))
    return descDiffList


def process(fileName, texFile, idx):
    inMols = Chem.SDMolSupplier(fileName+".sdf")
    outMols = []
    allMols = []
    for inMol in inMols:
        if inMol:
            outMol = standardizeCA(inMol)
            #outMol = standardizeMolVS(inMol)
            outMol = standardizeMolVS(outMol)
        else:
            outMol = inMol
        allMols.append(inMol)
        allMols.append(outMol)
        outMols.append(outMol)
    try:
        tanDists, tanDistsMorgan = getTanDist(outMols)  # Between all outMols
        descDiffs = getDescDiff(outMols)
    except:
        tanDists = "None"
        tanDistsMorgan = "None"
        descDiffs = "None"
    image = Draw.MolsToGridImage(allMols, molsPerRow=2, subImgSize=(200,200))
    image.save(fileName+"MV.png")
    #image.save(fileName+"CA.png")
    fid = open(texFile, "a")
    fid.write(str(idx)+" & "+fileName+" & \\includegraphics[scale=0.6]{"+fileName+"MV.png} & & "+str(tanDists)+"& "+str(tanDistsMorgan)+" & "+str(descDiffs)+" \\\\\n")
    #fid.write(str(idx)+" & "+fileName+" & \\includegraphics[scale=0.6]{"+fileName+"CA.png} & & "+str(tanDists)+"& "+str(tanDistsMorgan)+" & "+str(descDiffs)+" \\\\\n")
    fid.write("\\hline\n")
    fid.close()




if __name__ == "__main__":

    print "Change: table name, fileName with png * 2, method call, module load cheminfo"
    #texFile = "tableCA.tex"
    texFile = "tableVS.tex"
    fid = open(texFile, "w")
    fid.write("\\begin{longtable}{|l|l|l|l|l|l|l|}\n")
    fid.write("\\hline\n")
    fid.write("Idx & Name & Structure In and Out & Comment & TanTop & DiceMorgan & No. of DescDiff\\\\\n")
    fid.write("\\hline\n")
    fid.close()
    files = os.listdir(".")
    #files = ["Viagra.sdf"]
    #files = ["101and102.sdf"]
    #files = ["diazoMVtest.sdf"]
    #files = ["sulfoneMVtest.sdf"]
    #files = ["diazo.sdf"]
    #files = ["phosphine.sdf"]
    #files = ["NNO.sdf"]
    #files = ["azide.sdf"]
    idx = 1
    for file in files:
        fileName, ext = os.path.splitext(file)
        if ext == ".sdf":
            print fileName
            process(fileName, texFile, idx)
            idx = idx + 1
    fid = open(texFile, "a")
    fid.write("\\end{longtable}")
    fid.close()
