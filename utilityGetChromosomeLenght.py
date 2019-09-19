import numpy as np

def initChromosomeLenght(path):
    content = {}
    with open(path)as f:
        for line in f:
            pairs = np.asarray(line.strip().split('\t'))
            content[pairs[0]] = pairs[1]
    return content
