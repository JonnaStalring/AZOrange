import string

#smiFile = "XEN025smi.smi"
smiFile = "LiuJCIM2015smi3.smi"
csvFile = "SMARTcyp/SMARTCyp_Results_2016-09-29_10-31-50.csv"
#outFile = "XEN025SMARTcyp.txt"
outFile = "LiuJCIM2015SMARTcyp3.txt"

# Create dictionary with MVnumber and order
fid = open(smiFile)
idx = 1
mvDict = {}
for line in fid:
    lineList = string.split(line, "\t")
    mvnumber = string.strip(lineList[1])
    mvDict[mvnumber] = str(idx)
    idx = idx + 1
fid.close()


# The order of the molecules in the csv file follows the smi file
fid = open(csvFile)
header = fid.readline()
molDict = {}
for line in fid:
    lineList = string.split(line, ",")
    molIdx = lineList[0]
    rank = lineList[2] 
    score = lineList[3]
    energy = lineList[4]
    if molIdx not in molDict.keys():
        molDict[molIdx] = {}
    if rank == "1":    
        if "1" not in molDict[molIdx]:
            molDict[molIdx]["1"] = {}
            molDict[molIdx]["1"]["score"] = []
            molDict[molIdx]["1"]["energy"] = []
        molDict[molIdx]["1"]["score"].append(score)
        molDict[molIdx]["1"]["energy"].append(energy)
    elif rank == "2":    
        if "2" not in molDict[molIdx]:
            molDict[molIdx]["2"] = {}
            molDict[molIdx]["2"]["score"] = []
            molDict[molIdx]["2"]["energy"] = []
        molDict[molIdx]["2"]["score"].append(score)
        molDict[molIdx]["2"]["energy"].append(energy)
    elif rank == "3":    
        if "3" not in molDict[molIdx]:
            molDict[molIdx]["3"] = {}
            molDict[molIdx]["3"]["score"] = []
            molDict[molIdx]["3"]["energy"] = []
        molDict[molIdx]["3"]["score"].append(score)
        molDict[molIdx]["3"]["energy"].append(energy)
fid.close()

outFid = open(outFile, "w")
outFid.write("MVnumber\tMolIdx\tScore1\tScore2\tScore3\tEnergy1\tEnergy2\tEnergy3\n")
for key, value in mvDict.iteritems():
     print "Writing ", key, value, molDict[value]
     score1 = None
     score2 = None
     score3 = None
     energy1 = None
     energy2 = None
     energy3 = None
     for rank, results in molDict[value].iteritems():
         results["energy"].sort()
         if rank == "1":
             score1 = results["score"][0]
             energy1 = results["energy"][0]
             if len(results["energy"]) == 2:
                 score2 =  results["score"][0]
                 energy2 = results["energy"][1]
             elif len(results["energy"]) > 2:
                 score2 = results["score"][0]
                 energy2 =  results["energy"][1]
                 score3 = results["score"][0]
                 energy3 =  results["energy"][2]
         elif rank == "2":
             if score2 == None:
                 score2 = results["score"][0]
                 energy2 = max(results["energy"])
             if len(results["energy"]) > 1:
                 score3 =  results["score"][0]
                 energy3 = results["energy"][1]
         elif rank == "3":
             if score3 == None:
                 score3 = results["score"][0]
                 energy3 = max(results["energy"])
     print value, score1, score2, score3, energy1, energy2, energy3
     outFid.write(key+"\t"+value+"\t"+score1+"\t"+score2+"\t"+score3+"\t"+energy1+"\t"+energy2+"\t"+energy3+"\n")
outFid.close()
