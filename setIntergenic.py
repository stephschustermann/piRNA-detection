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
    df.index = range(len(df.index))
    return df

# check if the current coordinate is intergenic
# returns TRUE is the coordinate is intergenic otherwise FALSE
def isIntergenic(df, chromCoordinate):
    # filter locations by chromosome number
    chromName = chromCoordinate['chrom']
    chromLocations = df.loc[df['chrom'] == chromName]

    for index, row in chromLocations.iterrows():
        if row['chromStart'] == chromCoordinate['chromStart'] and row['chromEnd'] == chromCoordinate['chromEnd']:
            return False # is INTRA genic and not between genes
        
    return True

#isIntergenic({ 'chrom': 'chr22', 'chromStart': 10736285, 'chromEnd': 10736295 }, 'output.bed')
