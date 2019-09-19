import pandas as pd
import numpy as np
import numOfMapping
import setLabel

sequencesDataFrame = pd.read_csv('./featuresData/allSequencesData.txt')
labelsTable = setLabel.readLabelsFile('./featuresData/sncRNA_labels.txt')

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

# readFeaturesFile = open('./featuresData/allReadFeaturesData.txt',"w+")
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
    labelName = setLabel.findLabel(labelsTable, read)
    features.append(labelName)
    readFeaturesList.append(features)
        

readFeaturesDataFrame = pd.DataFrame(data=readFeaturesList, columns=['sequence', 'isInCluster', 'isIntergenic', 'isTE', 'Q25', 'Q50', 'Q75', 'Q100', 'numOfMappings', 'label'])
print('read features finish, creating read features data file')
# write all sequences in data frame
readFeaturesDataFrame.to_csv('./featuresData/allReadFeaturesData.txt', index=False)


