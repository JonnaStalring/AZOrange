from AZutilities import dataUtilities

#file1 = "XEN025dragonNewHeaderRespRDKbulk.txt"
#file2 = "XEN025SMARTcyp.txt"
file1 = "LiuJCIM2015SMARTcyp.txt"
file2 = "LiuJCIM2015dragonNewHeaderRespRDKbulk.txt"

data1 = dataUtilities.DataTable(file1)
data2 = dataUtilities.DataTable(file2)

#data = dataUtilities.horizontalMerge(data1, data2, "NAME", "MV Number")
#data = dataUtilities.horizontalMerge(data1, data2, "NAME", "MVnumber")
data = dataUtilities.horizontalMerge(data1, data2, "MVnumber", "NAME")

#data.save("XEN025dragonNewHeaderRespRDKbulkSMARTcyp.txt")
data.save("LiuJCIM2015dragonNewHeaderRespRDKbulkSMARTcyp.txt")

   


