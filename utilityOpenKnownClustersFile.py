import csv

def getKnownClusters(path, startCol=None, endCol=None):
    knownClusters = []
    with open(path, 'r') as csvFile:
        reader = csv.reader(csvFile)
        for row in reader:
            knownClusters.append(row[startCol:endCol])
    
    csvFile.close()
    return knownClusters

def removeComasFromStringNumber(table, columnsToConvert):
    for rIndex in range(1, len(table)):
       for  cIndex in range(len(columnsToConvert)):
           table[rIndex][columnsToConvert[cIndex]] = table[rIndex][columnsToConvert[cIndex]].replace(',', '')
    return table
           