import numpy as np

def getBedFileContent(path, ignoreHeader=False, startCol=None, endCol=None):
    content = []
    with open(path)as f:
        for line in f:
            if not ignoreHeader or not line.startswith("#"):
               content.append(np.asarray(line.strip().split())[startCol:endCol])
    return content

def getBedFileContentWithTab(path, ignoreHeader=False, startCol=None, endCol=None):
    content = []
    with open(path)as f:
        for line in f:
            if not ignoreHeader or not line.startswith("#"):
               content.append(np.asarray(line.strip().split('\t'))[startCol:endCol])
    return content
