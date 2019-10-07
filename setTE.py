# imports
import pandas as pd
from collections import Counter

df = []

def getElandFileContent(path):
    content = []
    header = True
    with open(path)as f:
        for line in f:
            if header:
                content.append(['probe_id', 'strand', 'type'])
                header = False
            else:
                row = line.strip().split('\t')
                print(row)
                strand = row[-1]
                teType = row[1]
                if row[1].isdigit(): # a row with no type cell
                    teType = ''
                probe_id = ''
                for cell in row:
                    if cell.startswith('HWI'):
                        probe_id = cell
                content.append([probe_id, strand, teType])
    return content


# read TE instersection table with reads
def readIntersectionWithSequences(intersectionFilePath):
    # for debug: intersectionFilePath = './featuresData/repbase.out.eland'
    content = getElandFileContent(intersectionFilePath)
    df = pd.DataFrame(data=content[1:], columns=content[0])
    column = list(df['probe_id'])
    return Counter(column)
# if read exists in TE intersection table return true
def mapsToTranposableElement(intersectionTEtable, sequenceId):
    if intersectionTEtable[sequenceId] > 0:
        return  df.loc[df['probe_id'] == sequenceId] # TODO https://github.com/users/stephschustermann/projects/1#card-27427098
    return False

