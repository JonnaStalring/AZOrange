import operator
import numpy
import Orange
from AZutilities import dataUtilities

# Change the resp var
#data = dataUtilities.DataTable("GlobalNoRSVdragonResp.txt")
#attrList = []
#for attr in data.domain.attributes:
#    #if attr.name == "CLint":
#    if attr.name == '"HLM_XEN025;Mean;CLint (uL/min/mg);(Num)"':
#        respAttr = attr
#    else:
#        attrList.append(attr)
#attrList.append(data.domain.classVar)
#
#print dir(Orange.data)
#newDomain = Orange.data.Domain(attrList, respAttr) 
#newData = dataUtilities.DataTable(newDomain, data)
#newData.save("GlobalNoRSVdragonResp2.txt")

def getTopTuple(pDict, pList, nTop):

    topRank = pList[0:nTop]
    topTupleList = []
    for elem in topRank:
        for key, value in pDict.iteritems():
            if  value == elem:
                topTupleList.append((key, value)) 
    return topTupleList


def getRankedAbs(data, nTop):

    # Calculate the Pearson correlation between each variable and the response. Collect in pDict
    pDict = {}
    pList = []
    fid = open("PearsonCorrResp.txt", "w")
    fid.write("Attribute\tPcorr\n")
    for attr in data.domain.attributes:
        attrName = attr.name
        valueList = []
        respList = []
        nFails = 0
        for ex in data:
            failed = False
            try:
                value = float(ex[attrName].value)
                resp = float(ex.get_class().value)
            except:
                nFails = nFails + 1
                failed = True
            if not failed:
                valueList.append(value)
                respList.append(resp)
        try:
            pCorr = abs(numpy.corrcoef(valueList, respList)[0][1])
            if not numpy.isnan(pCorr):
                pDict[attrName] = pCorr
                pList.append(pCorr)
                fid.write(attrName+"\t"+str(pCorr)+"\n")
        except:
            pass
    fid.close()

    # Return a tuple of the nTop attributes and their correlation sorted in decending order
    pList.sort(reverse = True)
    topTupleList = getTopTuple(pDict, pList, nTop)

    return topTupleList


def getRanked(data, nTop):

    # Calculate the Pearson correlation between each variable and the response. Collect in pDict
    pDict = {}
    pList = []
    fid = open("PearsonCorrResp.txt", "w")
    fid.write("Attribute\tPcorr\n")
    for attr in data.domain.attributes:
        attrName = attr.name 
        valueList = []
        respList = []
        nFails = 0
        for ex in data:
            failed = False
            try:
                value = float(ex[attrName].value)
                resp = float(ex.get_class().value)
            except:
                nFails = nFails + 1
                failed = True
            if not failed:
                valueList.append(value)
                respList.append(resp)
        try:
            pCorr = numpy.corrcoef(valueList, respList)[0][1]
            if not numpy.isnan(pCorr):
                pDict[attrName] = pCorr
                pList.append(pCorr)
                fid.write(attrName+"\t"+str(pCorr)+"\n")
        except: 
            pass
    fid.close()

    # Return a tuple of the nTop attributes and their correlation sorted in decending order
    pList.sort(reverse = True)
    topTupleList = getTopTuple(pDict, pList, nTop)
    pList.sort()
    bottomTupleList = getTopTuple(pDict, pList, nTop)

    return topTupleList, bottomTupleList
    

if __name__ == "__main__":
    #data = dataUtilities.DataTable("GlobalNoRSVDragonSMARTcypResp.txt")
    data = dataUtilities.DataTable("GlobalNoRSVmold2Resp.txt")
    nTop = 20
    #attrList = ["Project Name", "Structure",       "MV Number"]
    attrList = ["name", "Number"]
    data = dataUtilities.attributeDeselectionData(data, attrList)

    topRankTuple, bottomRankTuple = getRanked(data, nTop)

    # Select top corr and anti corr variables
    nTop = 20
    nBottom = 10
    topAttrList = []
    bAttrList = []
    for elem in topRankTuple:
        if len(topAttrList) < nTop:
            print elem
            topAttrList.append(elem[0])
    for elem in bottomRankTuple:
        if len(bAttrList) < nTop:
            print elem
            bAttrList.append(elem[0])
    #attrList = topAttrList + bAttrList+["CLint", "Score1",  "Score2",  "Score3",  "Energy1", "Energy2", "Energy3"]
    attrList = topAttrList + bAttrList+["Clint"]
    data = dataUtilities.attributeSelectionData(data, attrList)
    data.save("GlobalNoRSVmold2RespTop30.txt")





