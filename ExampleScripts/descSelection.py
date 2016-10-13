import Orange
import string
from trainingMethods import AZorngRF
from AZutilities import dataUtilities
from AZutilities import evalUtilities
import PearsonRank

def printTestSetAcc(Model, test, learners, resultsFid, projectName, isRand):

    TL = 0
    TH = 0
    FL = 0
    FH = 0
    for ex in test:
        pred = Model(ex)
        actual = ex.get_class()
        #print pred, actual
        if actual == "Low":
            if pred == "Low":
                TL = TL + 1
            elif pred == "High":
                FH = FH + 1
        elif actual == "High":
            if pred == "High":
                TH = TH + 1
            elif pred == "Low":
                FL = FL + 1
    print "TH, TL, FH, FL"
    print TH, TL, FH, FL
    CM = [[TH, FL], [FH, TL]]
    CA = round(float(TH+TL)/(TH+TL+FH+FL),3)
    print CA
    MCC = round(evalUtilities.calcMCC(CM), 3)
    print "MCC ", MCC
    if isRand:
        resultsFid.write("NO_"+projectName+"_randTest\t"+str(TH)+"\t"+str(TL)+"\t"+str(FH)+"\t"+str(FL)+"\t"+str(CA)+"\t"+str(MCC)+"\n" )
    else:
        resultsFid.write(projectName+"_test\t"+str(TH)+"\t"+str(TL)+"\t"+str(FH)+"\t"+str(FL)+"\t"+str(CA)+"\t"+str(MCC)+"\n" )
    
    return round(MCC,3)


def printCV(train, learners, resultsFid, projectName):

    cv = Orange.evaluation.testing.cross_validation(learners, train, folds = 10)
    AUC = Orange.evaluation.scoring.AUC(cv)
    CM = Orange.evaluation.scoring.confusion_matrices(cv)
    CA = Orange.evaluation.scoring.CA(cv)
    MCC = Orange.evaluation.scoring.MCC(CM)
    print CM[0].TP, CM[0].TN, CM[0].FP, CM[0].FN
    print "CA", CA[0]
    print "MCC", MCC[0]
    #print "AUC", AUC[0]
    resultsFid.write("NO_"+projectName+"_CV\t"+str(int(CM[0].TN))+"\t"+str(int(CM[0].TP))+"\t"+str(int(CM[0].FN))+"\t"+str(int(CM[0].FP))+"\t"+str(round(CA[0],3))+"\t"+str(round(MCC[0],3))+"\n" )

    return round(MCC[0],3)


def OldchangeRespVar(data, newRespName):
    """ Make the attribute newRespName the respone variable and the response variable an attribute. Do not remove any newRespName as attribute"""
    if newRespName != data.domain.classVar.name:
        attrList = []
        attrNames = []
        for attr in data.domain.attributes:
            attrList.append(attr) 
            attrNames.append(attr.name) 
            if attr.name == newRespName:
                respVar = attr
        if data.domain.classVar.name not in attrNames:
            attrList = attrList + [data.domain.classVar]
        newDomain = Orange.data.Domain(attrList, respVar)
        newData = Orange.data.Table(newDomain, data)
    else:
        newData = data
    return newData



def changeRespVar(data, newRespName):
    """ Make the attribute newRespName the respone variable and the response variable an attribute. Do not remove any newRespName as attribute"""
    if newRespName != data.domain.classVar.name:
        attrList = []
        attrNames = []
        for attr in data.domain.attributes:
            if attr.name == newRespName:
                respVar = attr
            else:
                attrList.append(attr)
                attrNames.append(attr.name)
        if data.domain.classVar.name not in attrNames:
            attrList = attrList + [data.domain.classVar]
        newDomain = Orange.data.Domain(attrList, respVar)
        newData = Orange.data.Table(newDomain, data)
    else:
        newData = data
    return newData


def tuplePrint(rankTuple, nTop):
    idx = 1
    for pair in rankTuple:
        if idx < nTop:
            print pair[0]+"\t"+str(round(pair[1], 3))
            idx = idx + 1

def getRank(scoreTuple):
    # Assume scoreTuple is soreted in descending order with absolute values
    rankTuple = []
    idx = 1
    for elem in scoreTuple:
        rankTuple.append((elem[0], int(idx)))
        idx = idx + 1 
    return rankTuple


def getRankSum(rankSumList, rankList):

    # Return the sum of ranks in rankSum
    newRankList = []
    for rankTuple in rankSumList:
        name = rankTuple[0]
        rank = rankTuple[1]
        for elem in rankList:
            if elem[0] == name:
                rank = rank + elem[1]
        newRankList.append((name, rank))
    return newRankList


def sortRankSum(rankSumTuple, nMethods):

    valueList = []
    for elem in rankSumTuple:
        valueList.append(elem[1])
    valueList.sort()
    newRankSum = []
    for rank in valueList:
        for tupElem in rankSumTuple:
            if tupElem[1] == rank:
                newRankSum.append((tupElem[0], round(float(tupElem[1])/nMethods, 3)))
    return newRankSum


def getAccStat(rankSumTuple, nDesc, train, randTest, extTest, resultsFid, projectName):

    print "Select features based on top ranked features"
    attrList = []
    for elem in rankSumTuple:
        if len(attrList) < nDesc:
            attrList.append(elem[0])
    train = dataUtilities.attributeSelectionData(train, attrList)
    train = dataUtilities.attributeDeselectionData(train, ['HLM_XEN025;Mean;CLint (uL/min/mg);(Num)'])

    print train.domain.attributes, len(train.domain.attributes), train.domain.classVar

    # Get accuracies
    learners = [AZorngRF.RFLearner(nTrees=100)]
    print "CV accuracy"
    MCC_CV = printCV(train, learners, resultsFid, projectName)
    Model = learners[0](train)
    print "Random Test set accuracy"
    MCC_rand = printTestSetAcc(Model, randTest, learners, resultsFid, projectName, True)
    print "External Test set accuracy"
    MCC_ext = printTestSetAcc(Model, extTest, learners, resultsFid, projectName, False)

    return MCC_CV, MCC_rand, MCC_ext



def Wrapper(train, randTest, extTest, resultsFid, projectName, MCCdict, descList):

    # Select RDKit desc and keep CLint for reference
    attrList = []
    fid = open("RDKitDesc.txt")
    for line in fid:
        attrList.append(string.strip(line))
    fid.close()    
    #attrList = attrList + ["CLint"] 
    attrList = attrList + ["HLM_XEN025;Mean;CLint (uL/min/mg);(Num)"]
    train = dataUtilities.attributeSelectionData(train, attrList)

    print "Ranking with contiuous response"
    train = changeRespVar(train, "HLM_XEN025;Mean;CLint (uL/min/mg);(Num)")
    print "Cont resp ", train.domain.classVar

    print "Pearson ranking"
    scoreTuple = PearsonRank.getRankedAbs(train, 178)
    tuplePrint(scoreTuple, 20)
    rankTuple = getRank(scoreTuple)
    rankSumTuple = rankTuple
    #tuplePrint(rankSumTuple, 50)
    #tuplePrint(rankTuple, 20)

    print "ReliefF ranking "
    relief = Orange.feature.scoring.Relief(k=20, m=50)
    scoreTuple = Orange.feature.scoring.score_all(train, score=relief)
    #tuplePrint(scoreTuple, 20)
    rankTuple = getRank(scoreTuple)
#    rankSumTuple = rankTuple
    rankSumTuple = getRankSum(rankSumTuple, rankTuple)
#    tuplePrint(rankSumTuple, 50)
#    tuplePrint(rankTuple, 20)
#
    print "Ranking with categorical response"
    train = changeRespVar(train, "CLint")
    print "Cat resp ", train.domain.classVar
#
    print "ReliefF ranking "
    relief = Orange.feature.scoring.Relief(k=20, m=50)
    scoreTuple = Orange.feature.scoring.score_all(train, score=relief)
    tuplePrint(scoreTuple, 20)
    rankTuple = getRank(scoreTuple)
#    print "OBS!!!!!!!!!!!! Remove if already done"
#    rankSumTuple = rankTuple
    rankSumTuple = getRankSum(rankSumTuple, rankTuple)
    #tuplePrint(rankTuple, 20)
#
    print "MDL ranking "
    score = Orange.feature.scoring.MDL()
    scoreTuple = Orange.feature.scoring.score_all(train, score)
    tuplePrint(scoreTuple, 20)
    rankTuple = getRank(scoreTuple)
#    rankSumTuple = rankTuple
    rankSumTuple = getRankSum(rankSumTuple, rankTuple)
#    #tuplePrint(rankTuple, 20)
#    #tuplePrint(rankSumTuple, 50)
#
    # Return average rank in sorted list
    rankSumTuple = sortRankSum(rankSumTuple, 1)
    tuplePrint(rankSumTuple, 50)
    #print rankSumTuple

    for nDesc in descList:
        MCC_CV, MCC_rand, MCC_ext = getAccStat(rankSumTuple, nDesc, train, randTest, extTest, resultsFid, projectName)
        MCCdict[projectName][str(nDesc)] = [MCC_CV, MCC_rand, MCC_ext]

    return MCCdict


if __name__ == "__main__":

    #projectName = "ALZHEIMER_BETA"
    projectList = ["ALZHEIMER_BETA", "CATHEPSIN_S", "HCV_POL_NUC", "MALT_1", "RENIN"]
    #projectList = ["ALZHEIMER_BETA", "CATHEPSIN_S"]
    descList = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 100, 120, 140, 160, 177]
    resultsFile = "accResults.txt"
    descFile = "descSelectionResults.txt"


    resultsFid = open(resultsFile, "w")
    resultsFid.write("Data\tTH\tTL\tFH\tFL\tCA\tMCC\n")
    resultsFid.close()
    descFid = open(descFile, "w")
    headerStr = ""
    for project in projectList:
        headerStr = headerStr + "MCC_CV_NO_"+project+"\t"+"MCC_rand_NO_"+project+"\t"+"MCC_ext"+project+"\t"
    headerStr = string.strip(headerStr)
    descFid.write("nDesc\t"+headerStr+"\tMCC_CV_AVG\tMCC_Rand_AVG\tMCC_Ext_AVG\n")
    descFid.close()
    MCCdict = {}
    for projectName in projectList:
        train = dataUtilities.DataTable("XEN025_NO_"+projectName+"Train.txt")
        randTest = dataUtilities.DataTable("XEN025_NO_"+projectName+"RandTest.txt")
        extTest = dataUtilities.DataTable("XEN025"+projectName+"Test.txt")
        resultsFid = open(resultsFile, "a")
        MCCdict[projectName] = {}
        MCCdict = Wrapper(train, randTest, extTest, resultsFid, projectName, MCCdict, descList)
        resultsFid.close()

    print MCCdict        
    descFid = open(descFile, "a")
    for nDesc in descList:
        wrtStr = ""
        descSumCV = 0
        descSumRand = 0
        descSumExt = 0
        for project in projectList:
            wrtStr = wrtStr + str(MCCdict[project][str(nDesc)][0])+"\t"+str(MCCdict[project][str(nDesc)][1])+"\t"+str(MCCdict[project][str(nDesc)][2])+"\t"
            descSumCV = descSumCV + MCCdict[project][str(nDesc)][0]
            descSumRand = descSumRand + MCCdict[project][str(nDesc)][1]
            descSumExt = descSumExt + MCCdict[project][str(nDesc)][2]
        descFid.write(str(nDesc)+"\t"+string.strip(wrtStr)+"\t"+str(descSumCV/len(projectList))+"\t"+str(descSumRand/len(projectList))+"\t"+str(descSumExt/len(projectList))+"\n")
    descFid.close()





