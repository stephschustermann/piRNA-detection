import numpy as np

# import and init files for features
import setDistanceFromHeterochromatin
import setInKnownCluster
import setIntergenic
import setTE
import utilityGetChromosomeLenght

print('init files')
heterochromatinDataFrame = setDistanceFromHeterochromatin.initHeterochromatinFiles()
knownClustersDataFrame = setInKnownCluster.initKnownClustersFile('./featuresData/intersectionWithKnownClusters.bed')
transposonElementsDictionary = setTE.readIntersectionWithSequences('./featuresData/repbase.out.eland')
intergenicDataFrame = setIntergenic.initIntergenicFile('./featuresData/intersectionWithKnownGenes.bed')
chromosomeLengths = utilityGetChromosomeLenght.initChromosomeLenght('./featuresData/hg38.chrom.sizes.txt')

def setFeaturesToElandFile(elandFileName, featuresFileName):
    featuresFile = open(featuresFileName,"w+")
    featuresFile.write("chrom,chromStart,chromEnd,probe_id,heterochromatinDistance,isInCluster,isIntergenic,isTE\n")

    with open(elandFileName)as f:
        for line in f:
            sequencesArray = np.asarray(line.strip().split('\t'))
            chrom = sequencesArray[0].split()[0]
            # Only known chromosomes
            if chrom in chromosomeLengths:
                coordStart = int(sequencesArray[1])
                seqLength = len(sequencesArray[2])
                probe_id = sequencesArray[3]
                coordEnd = coordStart + seqLength - 1
    
                sequence = {}
                sequence['chrom'] = chrom
                sequence['chromStart'] = int(sequencesArray[1])
                sequence['chromEnd'] = coordEnd
                sequence['probe_id'] = probe_id
            
                # add heterochromatinDistance
                chromatinDistance = setDistanceFromHeterochromatin.getClosestHeterochromatinDistance(heterochromatinDataFrame, sequence)
                sequence['heterochromatinDistance'] = chromatinDistance / int(chromosomeLengths[sequence['chrom']]) # normal to chromosome length
                # add isInCluster
                isInClusterResult = setInKnownCluster.isSequenceInsideCluster(knownClustersDataFrame, sequence['probe_id'])
                sequence['isInCluster'] = isInClusterResult
                
                # add isIntergenic
                isIntergenicResult = setIntergenic.isIntergenic(intergenicDataFrame, sequence)
                sequence['isIntergenic'] = isIntergenicResult
                
                # add maps to TE
                isTE = setTE.mapsToTranposableElement(transposonElementsDictionary, sequence['probe_id'])
                sequence['isTE'] = isTE
                
                sequenceToWrite = []
                for value in sequence.values():
                    sequenceToWrite.append(str(value))
                featuresFile.write(','.join(sequenceToWrite))
                featuresFile.write('\n')
                
    featuresFile.close()

setFeaturesToElandFile('./featuresData/genome.out.eland', './featuresData/allSequencesData.txt')
elandFileName = './featuresData/genome.out.eland'
featuresFileName = './featuresData/allSequencesData.txt'


