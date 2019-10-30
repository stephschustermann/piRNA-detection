# import and init files for features
import setDistanceFromHeterochromatin
import setInKnownCluster
import setIntergenic
import setTE
import setLabel
import utilityGetChromosomeLenght

print('init files')
heterochromatinDataFrame = setDistanceFromHeterochromatin.initHeterochromatinFiles()
knownClustersDataFrame = setInKnownCluster.initKnownClustersFile('./featuresData/intersectionWithKnownClusters.bed')
transposonElementsDF = setTE.getDataframeTE('./featuresData/repbase.out.eland')
transposonElementsDictionary = setTE.readIntersectionWithSequences(transposonElementsDF)
transposonElementsTypes = setTE.getTEtypes(transposonElementsDF)
labelsTable = setLabel.readLabelsFile('./featuresData/sncRNA_labels.txt')
intergenicDataFrame = setIntergenic.initIntergenicFile('./featuresData/intersectionWithKnownGenes.bed')
chromosomeLengths = utilityGetChromosomeLenght.initChromosomeLenght('./featuresData/hg38.chrom.sizes.txt')

# convert the eland file into bed file
import utilityConvertElandToBed
smallRNAFileName = './featuresData/smallRNASequences.bed'
print('convert genome file to bed')
utilityConvertElandToBed.createBedFileFromEland(smallRNAFileName, './featuresData/genome.out.eland')
                          
# read smallRNA sequences file
import utilityOpenBedFile
sequences = []
sequences = utilityOpenBedFile.getBedFileContentWithTab(smallRNAFileName)

# convert sequences to data frame
import pandas as pd

print('get data frame')
sequencesDataFrame = pd.DataFrame(data=sequences[1:], columns=sequences[0])
print('changing some values in data frame')
# sort the sequences by chromosome and locations
sequencesDataFrame.chromStart = sequencesDataFrame.chromStart.astype(float)
sequencesDataFrame.chromEnd = sequencesDataFrame.chromEnd.astype(float)
sequencesDataFrame = sequencesDataFrame.sort_values(by=['probe_id', 'chrom','chromStart'])
sequencesDataFrame.index = range(len(sequencesDataFrame.index))

sequencesDataFrame['heterochromatinDistance'] = 0
sequencesDataFrame['isInCluster'] = False
sequencesDataFrame['isIntragenic'] = False
sequencesDataFrame['geneSense'] = 0
sequencesDataFrame['isTE'] = False
sequencesDataFrame['teAntisense'] = False

for teType in transposonElementsTypes:
    sequencesDataFrame[teType] = 0
print('adding features')
# for each sequence, add a feature
for index, sequence in sequencesDataFrame.iterrows():
    # Only known chromosomes
    if sequence['chrom'] in chromosomeLengths:
        # add heterochromatinDistance
        chromatinDistance = setDistanceFromHeterochromatin.getClosestHeterochromatinDistance(heterochromatinDataFrame, sequence)
        sequencesDataFrame.at[index, 'heterochromatinDistance'] = chromatinDistance / int(chromosomeLengths[sequence['chrom']]) # normal to chromosome length
        
        # add isInCluster
        isInClusterResult = setInKnownCluster.isSequenceInsideCluster(knownClustersDataFrame, sequence['probe_id'])
        sequencesDataFrame.at[index, 'isInCluster'] = isInClusterResult
        
        # add isIntragenic and sense
        isIntergenicResult = setIntergenic.isIntragenic(intergenicDataFrame, sequence)
        if not isIntergenicResult:
            sequencesDataFrame.at[index, 'isIntragenic'] = isIntergenicResult
            sequencesDataFrame.at[index, 'geneSense'] = -1
        else:
            sequencesDataFrame.at[index, 'isIntragenic'] = True
            if isIntergenicResult[4] == '+':
                sequencesDataFrame.at[index, 'geneSense'] = 1
            else:
                sequencesDataFrame.at[index, 'geneSense'] = 0
        
        # add maps to TE
        isTE = setTE.mapsToTranposableElement(transposonElementsDictionary, sequence['probe_id'], transposonElementsDF)
        if not isTE:
            sequencesDataFrame.at[index, 'isTE'] = isTE
        else:
            sequencesDataFrame.at[index, 'isTE'] = True
            sequencesDataFrame.at[index, 'teAntisense']= isTE['antisense']
            for teType in isTE['types']:
                sequencesDataFrame.at[index, teType] = 1
    else: #not a defined chromosome so need to drop
        sequencesDataFrame = sequencesDataFrame.drop(sequencesDataFrame.index[index])

print('features finish, creating sequences data file')
# write all sequences in data frame
sequencesDataFrame.to_csv('./featuresData/allSequencesData.txt', index=False)
print('get unique reads')
# get all unique reads
sequencesProbeId = list(set(sequencesDataFrame['probe_id']))
print('calculate heterochromatin distances')
