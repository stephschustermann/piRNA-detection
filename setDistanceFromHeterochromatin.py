import pandas as pd
import numpy as np
import utilityOpenBedFile as openBedFile
import utilityMergeCoordinates as mergeCoordinates

def initHeterochromatinFiles():
    # open heterochromatin files and save them as a list of lists
    content = []
    content = openBedFile.getBedFileContent("./featuresData/gap.bed",
                                            startCol=1,
                                            endCol=4)
    content = content + openBedFile.getBedFileContent("./featuresData/centromeres.bed",
                                                      ignoreHeader=True,
                                                      startCol=1,
                                                      endCol=4)
    
    # convert content to data frame           
    df = pd.DataFrame(data=content[2:], columns=content[1])
    
    # sort the heterochromatin by chromosome and locations
    df.chromStart = df.chromStart.astype(float)
    df.chromEnd = df.chromEnd.astype(float)
    df = df.sort_values(by=['chrom','chromStart'])
    df.index = range(len(df.index))
    return mergeCoordinates.mergeCoordinates(df)

def closest_node(node, nodes):
    from scipy.spatial import distance
    nodes = np.asarray(nodes)
    dist_2 = np.sum(distance.euclidean(nodes, node), axis=1)
    return np.argmin(dist_2)

def getClosestHeterochromatinDistance(df, chromCoordinate):
    # disable pandas warning
    pd.options.mode.chained_assignment = None  # default='warn'

    # filter locations by chromosome number
    chromName = chromCoordinate['chrom']
    chromLocations = df.loc[df['chrom'] == chromName]
    
    if chromLocations.empty:
        return -1
    
    # convert the table to number of coordinates
    chromLocations.drop(['chrom'], axis = 1, inplace = True)
    chromLocations.chromStart = chromLocations.chromStart.astype(float)
    chromLocations.chromEnd = chromLocations.chromEnd.astype(float)
    chromLocations.index = range(len(chromLocations.index))
    
    # find the closest heterochromatin location to each given coordinate
    chromStart = chromCoordinate['chromStart']
    chromEnd = chromCoordinate['chromEnd']
    from scipy import spatial
    tree = spatial.KDTree(chromLocations)
    queryResult = tree.query([(chromStart, chromEnd)])
    
    # get the closest heterochromatin instance from data frame
    closestChromatin = chromLocations.loc[queryResult[1][0], :]

    # calculate the distance from the closest heterochromatin
    if closestChromatin[0] <= chromStart:
        if closestChromatin[1] <= chromEnd:
            return chromStart - closestChromatin[1]
        else:
            return -1
    else:
        if closestChromatin[1] <= chromStart:
            return closestChromatin[0] - chromEnd
        else:
            return -1
#df = initHeterochromatinFiles()
#getClosestHeterochromatinDistance(df,
# { 'chrom': 'chr22', 'chromStart': 1000012, 'chromEnd': 1040012 })