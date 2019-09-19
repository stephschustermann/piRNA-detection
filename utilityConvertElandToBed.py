import numpy as np

def getElandFileContent(path):
    content = []
    with open(path)as f:
        for line in f:
            content.append(np.asarray(line.strip().split('\t')))
    return content

def setChromosomeName(content):
    for row in content:
        firstCell = row[0]
        chromosomeName = firstCell.split()[0]
        row[0] = chromosomeName
    return content

def createBedFileFromEland(newBedName, elandName):
    content = []
    content = getElandFileContent(elandName)
    content = setChromosomeName(content)
        
    f= open(newBedName,"w+")
    index = 0
    f.write("chrom\tchromStart\tchromEnd\tprobe_id\tstrand\n")
    for row in content:
        if row[0].startswith("chr"):
            coordStart = int(row[1])
            seqLength = len(row[2])
            probe_id = row[3]
            positional = row[-1]
            coordEnd = coordStart + seqLength - 1
            toPrint = [row[0], row[1], str(coordEnd), probe_id, positional]
            index = index + 1
            f.write("\t".join(toPrint)+'\n')
         
    f.close()


def intersectBedFiles():
    import pybedtools
    a = pybedtools.example_bedtool('./Users/stephanieschustermann/piRNA-classifier/piRNA-classifier/featuresData/newBedFile.bed')
    b = pybedtools.example_bedtool('./Users/stephanieschustermann/piRNA-classifier/piRNA-classifier/featuresData/GENCODE_v29_knownGene.bed')
    print(a.intersect(b))
    return a.intersect(b)

#createBedFileFromEland('../featuresData/mostCommonSequence.bed', '../featuresData/mostCommonSequence.eland')
#createBedFileFromEland('./featuresData/genome.out.bed', '../cd-hit/genome.out.eland')