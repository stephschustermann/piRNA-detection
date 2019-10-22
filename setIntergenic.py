# open the genes file
import pandas as pd
import utilityOpenBedFile

def initIntergenicFile(intergenicBedFilePath):
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
def isIntragenic(geneDict, sequence):
    #TODO we need to somehow get the gene info, such as type, utr3/5, intron/exon
    if sequence['name'] in geneDict:
        allSequences = geneDict[sequence['name']]
        for row in allSequences:
            if row['chromStart'] == sequence['chromStart'] and row['chromEnd'] == sequence['chromEnd'] and row['chrom'] == sequence['chrom']:
                return row  # is INTRA genic and not between genes
    
    return False
