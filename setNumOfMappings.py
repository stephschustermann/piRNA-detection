from collections import Counter

def numOfMappingPerRead(genomeDataFrame):
    column = list(genomeDataFrame['probe_id'])
    return Counter(column)

def getNumOfMapping(numOfMappingPerReadTable, sequenceId):
    return numOfMappingPerReadTable[sequenceId]