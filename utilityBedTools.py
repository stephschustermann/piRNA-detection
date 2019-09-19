def callBedTools(aFile, bFile, output):
    import os
    os.system('bedtools intersect -nonamecheck -wa -a ' + aFile + ' -b ' + bFile + ' > ' + output)
    
    #bedtools intersect -nonamecheck -wa -a mostCommonSequence.bed -b GENCODE_v29_knownGene.bed intersectionWithKnownGenes.bed