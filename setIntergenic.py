# open the genes file
import pandas as pd
import utilityOpenBedFile

def initIntergenicFile(intergenicBedFilePath):
    intergenicBedFilePath = './featuresData/intersectionWithKnownGenes.bed'
    content = []
    content = utilityOpenBedFile.getBedFileContentWithTab(intergenicBedFilePath)

    # convert content to data frame
    df = pd.DataFrame(data=content, columns=['chrom', 'chromStart', 'chromEnd', 'name', 'strand'])

    # sort the heterochromatin by chromosome and locations
    df.chromStart = df.chromStart.astype(float)
    df.chromEnd = df.chromEnd.astype(float)
    df = df.sort_values(by=['chrom','chromStart'])
    df = df.drop_duplicates()
    df.index = range(len(df.index))
    genesList = df.values.tolist()
    geneDict = {}
    for gene in genesList:
        if gene[3] not in geneDict:
            geneDict[gene[3]] = []
        geneDict[gene[3]].append(gene)
    return geneDict

# check if the current coordinate is intergenic
# returns TRUE is the coordinate is intergenic otherwise FALSE
def isIntergenic(geneDict, sequence):
    if sequence['name'] in geneDict: # TODO how do we want to separate all the matches?
        return geneDict[sequence['name']] # is INTRA genic and not between genes
    
    return True

#isIntergenic({ 'chrom': 'chr22', 'chromStart': 10736285, 'chromEnd': 10736295 }, 'output.bed')
