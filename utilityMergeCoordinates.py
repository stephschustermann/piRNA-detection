import pandas as pd

def canMerge(coord1, coord2):
    return coord1.chrom == coord2.chrom and coord1.chromEnd == coord2.chromStart

def mergeCoordinates(coordinateTable):
    mergedTable = []
    for index, row in coordinateTable.iterrows():
        if index > 0:
            lastRow = mergedTable[-1]
            if canMerge(lastRow, row):
                lastRow.chromEnd = row.chromEnd
            else:
                mergedTable.append(row)
        else:
            mergedTable.append(row)
    mergedDataFrame = pd.DataFrame(data=mergedTable, columns=coordinateTable.columns)
    return mergedDataFrame
