import string

#smiFile = "RSV.smi"
#mold2File = "RSVmold2.txt"
#outFile = "RSVmold2Resp.txt"
smiFile = "GlobalNoRSV.txt"
mold2File = "GlobalNoRSVmold2.txt"
outFile = "GlobalNoRSVmold2Resp.txt"

smiFid = open(smiFile)
header = smiFid.readline()
respDict = {}
idx = 1
for line in smiFid:
    lineList = string.split(line, "\t")
    name = lineList[3]
    contResp = lineList[1]
    resp = string.strip(lineList[4])
    respDict[str(idx)] = [name, resp, contResp]
    idx = idx + 1
smiFid.close()
    
idx = 1
mold2Fid = open(mold2File)
header = string.strip(mold2Fid.readline())
for line in mold2Fid:
    lineStr = string.strip(line)
    respDict[str(idx)].append(lineStr)
    idx = idx + 1
mold2Fid.close()

print len(respDict)

outFid = open(outFile, "w")
outFid.write("name\t"+header+"\tClint\tCLintCont\n")
for key, value in respDict.iteritems():
    outFid.write(respDict[key][0]+"\t"+respDict[key][3]+"\t"+respDict[key][1]+"\t"+respDict[key][2]+"\n")
outFid.close()
