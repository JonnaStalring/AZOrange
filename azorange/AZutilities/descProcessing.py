from AZutilities import dataUtilities
import numpy

def filterDesc(data, zeroFracT = 0.95, LowVarT = 0.01, HighCorrT = 0.95):

    print "Initial number of descriptors ", len(data.domain.attributes)

    # Find the descriptors for which the fraction of zeros is smaller than zeroFracT - keep these
    attrList = []
    rmAttrList = []
    for attr in data.domain.attributes:
        valueList = []
        nZero = 0
        for ex in data:
            value = ex[attr.name].value
            if value == 0:
                nZero = nZero + 1
            valueList.append(value)
        zeroFrac = float(nZero)/len(valueList)
        if zeroFrac < zeroFracT:
            attrList.append(attr.name)
        else:
            rmAttrList.append(attr.name)
    print "Descriptors deselected because of a large fraction of zeros: "
    print rmAttrList
    data = dataUtilities.attributeSelectionData(data, attrList)
    print "Remaining number of descriptors ", len(data.domain.attributes)

    # Filter descriptors based on normalized variance
    rmAttrList = []
    for attr in data.domain.attributes:
        valueList = []
        for ex in data:
            value = ex[attr.name].value
            valueList.append(value)
        variance = numpy.var(valueList)
        mean = numpy.mean(valueList)
        normVar = variance/mean
        if normVar < LowVarT:
            rmAttrList.append(attr.name)
        
    print "Descriptors deselected because of low variance "
    print rmAttrList
    data = dataUtilities.attributeDeselectionData(data, rmAttrList)
    print "Remaining number of descriptors ", len(data.domain.attributes)
    
    print "Correlation filter not implemented yet"

    return data


