def initKnownClustersFile(clustersBedFilePath):
   content = []
   with open(clustersBedFilePath)as f:
       for line in f:
           row = line.strip().split('\t')
           probe_id = ''
           for cell in row:
               if cell.startswith('HWI'):
                   probe_id = cell
           content.append(probe_id)
   return content

# find if the current coordinate is inside a known cluster and find distance!
# returns zero for coordinates inside a cluster,
# otherwise the distance from the closest cluster
def isSequenceInsideCluster(df, chromCoordinate):
    if chromCoordinate in df :
        return True
    return False