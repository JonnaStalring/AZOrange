import scipy
import numpy
import math
from AZutilities import dataUtilities
from AZutilities import Mahalanobis
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys


def getRespVar(tanList, tanDict, train, nameAttr):
    
    # Find names of NN
    NNnames = []
    for key, value in tanDict.iteritems():
        if value in tanList:
            NNnames.append(key)

    # Find the responses of the NN
    respList = []
    for ex in train:
        if ex[nameAttr].value in NNnames:
            respList.append(ex.get_class().value)

    # Calculate the entropy amonst the responses
    if len(set(respList)) == 1:
        entropy = 0.0
    else:
        valueList = train.domain.classVar.values
        class1 = 0
        class2 = 0
        for elem in respList:
            if elem == valueList[0]:
                class1 = class1 + 1
            elif elem == valueList[1]:
                class2 = class2 + 1
        pList = [float(class1)/(class1+class2), float(class2)/(class1+class2)]
        entropy = scipy.stats.entropy(pList)

    return entropy


def getDescVect(predEx):
    descVect = []
    for attr in predEx.domain.attributes:
        descVect.append(predEx[attr.name].value)
    return descVect
 

def getXNN(trainSmilesList, train, predEx, smilesAttrName, nameAttr, X, simType):

    if simType == "Topological":
        fpsTrain = [FingerprintMols.FingerprintMol(x) for x in trainSmilesList]
        fp = FingerprintMols.FingerprintMol(Chem.MolFromSmiles(predEx[smilesAttrName].value))
    elif simType == "Morgan":
        fpsTrain = [AllChem.GetMorganFingerprint(x, 2) for x in trainSmilesList]
        fp = AllChem.GetMorganFingerprint(Chem.MolFromSmiles(predEx[smilesAttrName].value), 2)
    elif simType == "MACCS":
        fpsTrain = [MACCSkeys.GenMACCSKeys(x) for x in trainSmilesList]
        fp = MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(predEx[smilesAttrName].value))
    elif simType == "Mahalanobis":
        attrList = [smilesAttrName, nameAttr]
        predEx = dataUtilities.attributeDeselectionExample(predEx, attrList)
        fp = getDescVect(predEx)
        numTrain = dataUtilities.attributeDeselectionData(train, attrList)
        trainMat = []
        for ex in numTrain:
            descVect = getDescVect(ex) 
            trainMat.append(descVect)
        norm = Mahalanobis.create_inverse_covariance_norm(trainMat)
    else:
        print "This type of sim is not implemented ", simType

    simDict = {}
    idx = 0
    simList = []
    for ex in train:
        if simType == "Topological":
            sim = DataStructs.FingerprintSimilarity(fpsTrain[idx],fp)
        elif simType == "Morgan":
            sim = DataStructs.DiceSimilarity(fpsTrain[idx],fp)
        elif simType == "MACCS":
            sim = DataStructs.FingerprintSimilarity(fpsTrain[idx],fp)
        elif simType == "Mahalanobis":
            descVect = trainMat[idx]
            dist = Mahalanobis.compute_distance(fp, descVect, norm)
            sim = dist
        else:
            print "This type of sim is not implemented ", simType
        idx = idx + 1
        simDict[ex[nameAttr].value] = sim
        simList.append(sim)

    if simType == "Mahalanobis":  # Mahalanobis gives a distance while the other methods are similarities
        simList.sort()
    else:
        simList.sort(reverse = True)
    simList = simList[0:X]
    medSim = round(numpy.median(simList), 3)
    stdSim =  round(numpy.std(simList), 3)
    minSim =  round(min(simList), 3)
    maxSim =  round(max(simList), 3)

    entropy =  round(getRespVar(simList, simDict, train, nameAttr), 3)
    entropyClosest =  round(getRespVar(simList[0:X/2], simDict, train, nameAttr), 3)

    return medSim, stdSim, minSim, maxSim, entropy, entropyClosest


def getXNNstat(trainSmilesList, train, predEx, smilesAttrName, nameAttr, X, simTypes):

    XNNstatDict = {}

    if "Topological" in simTypes:
        medSim, stdSim, minSim, maxSim, entropy, entropyClosest = getXNN(trainSmilesList, train, predEx, smilesAttrName, nameAttr, X, simType)
        XNNstatDict["medTopSim"] = medSim
        XNNstatDict["stdTopSim"] = stdSim
        XNNstatDict["minTopSim"] = minSim
        XNNstatDict["maxTopSim"] = maxSim
        XNNstatDict["entropyTop"] = entropy
        XNNstatDict["entropyClosestTop"] = entropyClosest

    if "Morgan" in simTypes:
        medSim, stdSim, minSim, maxSim, entropy, entropyClosest = getXNN(trainSmilesList, train, predEx, smilesAttrName, nameAttr, X, simType)
        XNNstatDict["medMorSim"] = medSim
        XNNstatDict["stdMorSim"] = stdSim
        XNNstatDict["minMorSim"] = minSim
        XNNstatDict["maxMorSim"] = maxSim
        XNNstatDict["entropyMor"] = entropy
        XNNstatDict["entropyClosestMor"] = entropyClosest

    if "MACCS" in simTypes:
        medSim, stdSim, minSim, maxSim, entropy, entropyClosest = getXNN(trainSmilesList, train, predEx, smilesAttrName, nameAttr, X, simType)
        XNNstatDict["medMACCSSim"] = medSim
        XNNstatDict["stdMACCSSim"] = stdSim
        XNNstatDict["minMACCSSim"] = minSim
        XNNstatDict["maxMACCSSim"] = maxSim
        XNNstatDict["entropyMACCS"] = entropy
        XNNstatDict["entropyClosestMACCS"] = entropyClosest

    if "Mahalanobis" in simTypes:
        medSim, stdSim, minSim, maxSim, entropy, entropyClosest = getXNN(trainSmilesList, train, predEx, smilesAttrName, nameAttr, X, simType)
        XNNstatDict["medMahalDist"] = medSim
        XNNstatDict["stdMahalDist"] = stdSim
        XNNstatDict["minMahalDist"] = minSim
        XNNstatDict["maxMahalDist"] = maxSim
        XNNstatDict["entropyMahal"] = entropy
        XNNstatDict["entropyClosestMahal"] = entropyClosest
      
    return XNNstatDict


def getPredConfParam(train, predEx, smilesAttrName, nameAttr, simTypes = ["Morgan"]):

    trainSmilesList = []
    for trainEx in train:
        trainSmilesList.append(Chem.MolFromSmiles(trainEx[smilesAttrName].value))
    X = 10 # The number of NN
    XNNstatDict = getXNNstat(trainSmilesList, train, predEx, smilesAttrName, nameAttr, X, simTypes)

    return XNNstatDict


if __name__ == "__main__":

    print "******************Calculate metrics for binary data sets*******************"
    trainFile = "GlobalNoAlzheimerBulk_Train.txt"
    #testFile = "GlobalNoRSVBulk37SMARTcypSMILES_Test.txt"
    testFile = "GlobalNoAlzheimerBulk_Test.txt"
    #smilesAttrName = "SMILES"
    smilesAttrName = "Structure"
    nameAttr = "MV Number"
    outFile = "simMetricsNoAlzRandTest.txt"
     
    train = dataUtilities.DataTable(trainFile)
    test = dataUtilities.DataTable(testFile)

    fid = open(outFile, "w")
    fid.write("Name\tmedTopSim\tstdTopSim\tminTopSim\tmaxTopSim\tentropy\tentropyClosest\tmedMorSim\tstdMorSim\tminMorSim\tmaxMorSim\tentropyMor\tentropyClosestMor\tmedMACCSSim\tstdMACCSSim\tminMACCSSim\tmaxMACCSSim\tentropyMACCS\tentropyClosestMACCS\n")
    fid.close()
    for ex in test:
        XNNstatDict = getPredConfParam(train, ex, smilesAttrName, nameAttr)
        print XNNstatDict
        fid = open(outFile, "a")
        fid.write(ex[nameAttr].value+"\t"+str(XNNstatDict["medTopSim"])+"\t"+str(XNNstatDict["stdTopSim"])+"\t"+str(XNNstatDict["minTopSim"])+"\t"+str(XNNstatDict["maxTopSim"])+"\t"+str(XNNstatDict["entropyTop"])+"\t"+str(XNNstatDict["entropyClosestTop"])+"\t"+str(XNNstatDict["medMorSim"])+"\t"+str(XNNstatDict["stdMorSim"])+"\t"+str(XNNstatDict["minMorSim"])+"\t"+str(XNNstatDict["maxMorSim"])+"\t"+str(XNNstatDict["entropyMor"])+"\t"+str(XNNstatDict["entropyClosestMor"])+"\t"+str(XNNstatDict["medMACCSSim"])+"\t"+str(XNNstatDict["stdMACCSSim"])+"\t"+str(XNNstatDict["minMACCSSim"])+"\t"+str(XNNstatDict["maxMACCSSim"])+"\t"+str(XNNstatDict["entropyMACCS"])+"\t"+str(XNNstatDict["entropyClosestMACCS"])+"\n")
        fid.close()



