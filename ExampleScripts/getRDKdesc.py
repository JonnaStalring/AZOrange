import os
from AZutilities import dataUtilities
from AZutilities import getCinfonyDesc
import Orange
import orange

#fileName = "XEN025dragonNewHeaderResp.txt"
fileName = "LiuJCIM2015dragonNewHeaderResp.txt"
path, ext = os.path.splitext(fileName)
outFileName = path+"RDKbulk"+ext 

data = dataUtilities.DataTable(fileName)
descList = getCinfonyDesc.getAvailableDescs("rdkPhysChem")
newData = getCinfonyDesc.getRdkDescResult(data, descList)
#descList = getCinfonyDesc.getAvailableDescs("rdk")
#newData = getCinfonyDesc.getRdkDescResult(data, descList, radius = 3)
newData.save(outFileName)

