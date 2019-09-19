# imports
import pandas as pd
from collections import Counter

def getElandFileContent(path):
    content = []
    header = True
    with open(path)as f:
        for line in f:
            if header:
                content.append(['probe_id', 'strand'])
                header = False
            else:
                row = line.strip().split('\t')
                strand = row[-1]
                probe_id = ''
                for cell in row:
                    if cell.startswith('HWI'):
                        probe_id = cell
                content.append([probe_id, strand])
    return content


# read TE instersection table with reads
def readIntersectionWithSequences(intersectionFilePath):
    content = getElandFileContent(intersectionFilePath)
    df = pd.DataFrame(data=content[1:], columns=content[0])
    column = list(df['probe_id'])
    return Counter(column)
# if read exists in TE intersection table return true
def mapsToTranposableElement(intersectionTEtable, sequenceId):
    if intersectionTEtable[sequenceId] > 0:
        return True
    return False