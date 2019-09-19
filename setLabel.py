# -*- coding: utf-8 -*-
def readLabelsFile(path):
    d = {}
    with open(path)as f:
        for line in f:
            if line.strip():
                (key, val) = line.strip().split('\t')
                d[key] = val
    return d

def findLabel(labelsTable, sequenceId):
    if sequenceId in labelsTable:
        return labelsTable[sequenceId]
    return 'unknown'

