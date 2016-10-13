from AZutilities import dataUtilities
import orange
import Orange

data = dataUtilities.DataTable("XEN025dragonNewHeaderRespRDKbulkSMARTcypSP.txt")

classVar = orange.EnumVariable("CLint", values = ["High", "Low"])
attrList = data.domain.attributes+[data.domain.classVar] 
newDomain = Orange.data.Domain(attrList, classVar)
newData = dataUtilities.DataTable(newDomain, data)

for ex in newData:
    clint = ex['HLM_XEN025;Mean;CLint (uL/min/mg);(Num)'].value
    if clint < 80:
        ex["CLint"] = "Low"
    else:
        ex["CLint"] = "High"

newData.save("XEN025AllDesc.txt")
