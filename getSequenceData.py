# import and init files for features
import heterochromatinDistance
import isInCluster
import isIntergenic
import isTransposonElement
import numOfMapping
import label
import getChromosomeLength

print('init files')
heterochromatinDataFrame = heterochromatinDistance.initHeterochromatinFiles()
knownClustersDataFrame = isInCluster.initKnownClustersFile('./featuresData/intersectionWithKnownClusters.bed')
transposonElementsDictionary = isTransposonElement.readIntersectionWithSequences('./featuresData/repbase.out.eland')
labelsTable = label.readLabelsFile('./featuresData/sncRNA_labels.txt')
intergenicDataFrame = isIntergenic.initIntergenicFile('./featuresData/intersectionWithKnownGenes.bed')
chromosomeLengths = getChromosomeLength.initChromosomeLenght('./featuresData/hg38.chrom.sizes.txt')

# convert the eland file into bed file
import convertElandToBed
smallRNAFileName = './featuresData/smallRNASequences.bed'
print('convert genome file to bed')
convertElandToBed.createBedFileFromEland(smallRNAFileName, './featuresData/genome.out.eland')
                          
# read smallRNA sequences file
import openBedFile
sequences = []
sequences = openBedFile.getBedFileContentWithTab(smallRNAFileName)

# convert sequences to data frame
import pandas as pd
import numpy as np
print('get data frame')
sequencesDataFrame = pd.DataFrame(data=sequences[1:], columns=sequences[0])
print('changing some values in data frame')
# sort the sequences by chromosome and locations
sequencesDataFrame.chromStart = sequencesDataFrame.chromStart.astype(float)
sequencesDataFrame.chromEnd = sequencesDataFrame.chromEnd.astype(float)
sequencesDataFrame = sequencesDataFrame.sort_values(by=['probe_id', 'chrom','chromStart'])
sequencesDataFrame.index = range(len(sequencesDataFrame.index))

sequencesDataFrame['heterochromatinDistance'] = np.empty((len(sequencesDataFrame), 0)).tolist()
sequencesDataFrame['isInCluster'] = np.empty((len(sequencesDataFrame), 0)).tolist()
sequencesDataFrame['isIntergenic'] = np.empty((len(sequencesDataFrame), 0)).tolist()
sequencesDataFrame['isTE'] = np.empty((len(sequencesDataFrame), 0)).tolist()
print('adding features')
# for each sequence, add a feature
for index, sequence in sequencesDataFrame.iterrows():
    # Only known chromosomes
    if sequence['chrom'] in chromosomeLengths:
        # add heterochromatinDistance
        chromatinDistance = heterochromatinDistance.getClosestHeterochromatinDistance(heterochromatinDataFrame, sequence)
        sequencesDataFrame.at[index, 'heterochromatinDistance'] = chromatinDistance / int(chromosomeLengths[sequence['chrom']]) # normal to chromosome length
        
        # add isInCluster
        isInClusterResult = isInCluster.isSequenceInsideCluster(knownClustersDataFrame, sequence['probe_id'])
        sequencesDataFrame.at[index, 'isInCluster'] = isInClusterResult
        
        # add isIntergenic
        isIntergenicResult = isIntergenic.isIntergenic(intergenicDataFrame, sequence)
        sequencesDataFrame.at[index, 'isIntergenic'] = isIntergenicResult
        
        # add maps to TE
        isTE = isTransposonElement.mapsToTranposableElement(transposonElementsDictionary, sequence['probe_id'])
        sequencesDataFrame.at[index, 'isTE'] = isTE
    else:
        sequencesDataFrame = sequencesDataFrame.drop(sequencesDataFrame.index[index])

print('features finish, creating sequences data file')
# write all sequences in data frame
sequencesDataFrame.to_csv('./featuresData/allSequencesData.txt', index=False)
print('get unique reads')
# get all unique reads
sequencesProbeId = list(set(sequencesDataFrame['probe_id']))
print('calculate heterochromatin distances')
# calculate heterochromatin distances using percentile
heterochromatinDistances = sequencesDataFrame['heterochromatinDistance']
dist25 = np.percentile(heterochromatinDistances, 25)
dist50 = np.percentile(heterochromatinDistances, 50)
dist75 = np.percentile(heterochromatinDistances, 75)
dist100 = np.percentile(heterochromatinDistances, 100)

print('read number of sequences per each read')
numOfMappingInGenome = numOfMapping.numOfMappingPerRead(sequencesDataFrame)

# import pydash
import pydash
print('read features starting')
readFeaturesList = []
for read in sequencesProbeId:
    features = []
    features.append(read)
    
    numOfMappings = numOfMapping.getNumOfMapping(numOfMappingInGenome, read)
    sequencesOfRead = sequencesDataFrame.loc[sequencesDataFrame['probe_id'] == read]

    if numOfMappings > 0:
        features.append((sum(sequencesOfRead['isInCluster'])/numOfMappings)*100)
        features.append((sum(sequencesOfRead['isIntergenic'])/numOfMappings)*100)
        features.append((sum(sequencesOfRead['isTE'])/numOfMappings)*100)
        total25 = len(pydash.collections.filter_(list(sequencesOfRead['heterochromatinDistance']), lambda x: x <= dist25))
        total50 = len(pydash.collections.filter_(list(sequencesOfRead['heterochromatinDistance']), lambda x: x <= dist50 and x > dist25))
        total75 = len(pydash.collections.filter_(list(sequencesOfRead['heterochromatinDistance']), lambda x: x <= dist75 and x > dist50))
        total100 = len(pydash.collections.filter_(list(sequencesOfRead['heterochromatinDistance']), lambda x: x <= dist100 and x > dist75))
        
        features.append((total25/numOfMappings)*100)
        features.append((total50/numOfMappings)*100)
        features.append((total75/numOfMappings)*100)
        features.append((total100/numOfMappings)*100)
    else:
        features.append(0) # inCluster
        features.append(0) # inIntergenic
        features.append(0) # isTE
        features.append(0) # 25
        features.append(0) # 50
        features.append(0) # 75
        features.append(0) # 100
        
    features.append(numOfMappings)
    labelName = label.findLabel(labelsTable, read)
    features.append(labelName)
    readFeaturesList.append(features)
        

readFeaturesDataFrame = pd.DataFrame(data=readFeaturesList, columns=['sequence', 'isInCluster', 'isIntergenic', 'isTE', 'Q25', 'Q50', 'Q75', 'Q100', 'numOfMappings', 'label'])
print('read features finish, creating read features data file')
# write all sequences in data frame
readFeaturesDataFrame.to_csv('./featuresData/allReadFeaturesData.txt', index=False)

