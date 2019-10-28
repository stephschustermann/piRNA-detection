# imports
import pandas as pd
from collections import Counter

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

def getDataframeTE(intersectionFilePath):
    # intersectionFilePath = './featuresData/repbase.out.eland'
    content = getElandFileContent(intersectionFilePath)
    df = pd.DataFrame(data=content[1:], columns=content[0])
    return df
# read TE instersection table with reads
def readIntersectionWithSequences(df):
    column = list(df['probe_id'])
    return Counter(column)

def getTEtypes(df):
    types = list(df['type'])
    return list(dict.fromkeys(types)) # unique types

def getReturnObject(teList, df):
    result = {}
    result['antisense'] = False
    types = []
    for index, te in teList.iterrows():
        types.append(te['type'])
        if te[1] == '-':
            result['antisense'] = True
    types = list(dict.fromkeys(types))
    result['types'] = types
    return result

# if read exists in TE intersection table return true
def mapsToTranposableElement(intersectionTEtable, sequenceId, df):
    # if the sequence maps to more than one TE
    # return all the TE types
    # return true for mapping
    # return true if there is at least one TE mapped antisense
    if intersectionTEtable[sequenceId] > 0:
        return  getReturnObject(df.loc[df['probe_id'] == sequenceId], df)
    return False

