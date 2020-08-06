### IMPORT
import argparse
import contextlib
import gzip
import itertools



########## ########## ########## ########## ########## ########## 
### FUNCTIONS



def fArgParse():
    '''
    Interface for collecting and organizing CLI input data
    '''    
    ## Define parser
    vParser = argparse.ArgumentParser(description="Program_Description")
    ## Parser Arguements
    vParser.add_argument("-R1",  help="Read  1 File", required=True)
    vParser.add_argument("-R2",  help="Index 1 File", required=True)
    vParser.add_argument("-R3",  help="Index 2 File", required=True)
    vParser.add_argument("-R4",  help="Read  2 File", required=True)
    vParser.add_argument("-qc",  help="Phred Score Cutoff", required=True)
    ## Store CLI input in vArgs dictionary {'option':'input'}
    vArgs = vars(vParser.parse_args())
    return(vArgs)

def fCountReads(vfInfile,vfNumLines=4):
    '''
    Count the number of illumina reads in a .gz file
    '''
    ## Count total lines in file
    with gzip.open(vfInfile, 'rb') as vfFile:
        for vfCount, vfLine in enumerate(vfFile):
            pass
    ## Divide total lines by read length
    vfCount = (vfCount+1) / int(vfNumLines)
    return(int(vfCount))

def fRevComp(vfSeq):
    '''
    Generate the reverse complement of a sequence
    '''
    ## Dictionary for determining complmentary base
    vfComp = {"A":"T","T":"A","G":"C","C":"G",}
    ## Remove newline
    vfSeq = vfSeq.rstrip()  
    ## Reverse sequence
    vfSeq = vfSeq[::-1]
    ## Complement Sequence
    vfRevComp = ''
    for vfNuc in vfSeq:
        if vfNuc in vfComp.keys():
            vfRevComp += vfComp[vfNuc]
        else:
            vfRevComp += vfNuc
    return(vfRevComp)

def fGenPermutations(vfList,vfLen):
    '''
    Outputs a list all permutations for a given input list
    (vfList)   The list used to generate permutations from
    (vfLen)     How many items to combine together
    '''
    return(itertools.product(vfList,repeat=vfLen))
    # return(list(itertools.permutations(vfList,vfLen)))

def fFormatTuple(vfTuple,vfSeperator="-"):
    '''
    Converts a 2-length tuple into a string with a seperator character
    '''
    vfOutput = []
    for vfItem in vfTuple:
        vFormat = vfItem[0] + vfSeperator + vfItem[1]
        vfOutput.append(vFormat)
    return(vfOutput)

def fWriteList(vfList,vfFile):
    '''
    Write the contents of a list to a file
    '''
    for vItem in vfList:
        vfFile.write(vItem)

def fIsSublist(vfBigList,vfSubList):
    '''
    Returns True if all elements of a list os present in another list
    (vfBigList)     List, The all-inclusive list
    (vfSubList)     List, The list to analyze
    '''
    # vReturn = True
    # for vItem in vfSubList:
    #     if vItem not in vfBigList:
    #         vReturn = False
    #         break
    # return(vReturn)
    return(all(elem in vfBigList  for elem in vfSubList))

def fProgBar(vfPlace,vfRange):
    '''
    vfPlace     Current position in total length
    vfRange     The total length to measure progress across
    '''
    vfRange = int(vfRange/100)
    if vfPlace % vfRange == 0:
        print(vfPlace/vfRange+1, "% Complete")



########## ########## ########## ########## ########## ##########
### DEFINE VARIABLES



## Extract Input data
vInput = fArgParse()

## Identify Quality Score Cutoff (Inclusive)
vPhredCut = int(vInput["qc"])
print("\nQuality Score Identified:", vPhredCut)

## Identify acceptable Phredd+33 scores based on Quality Score Cutoff
vQualityPass = list("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ")[vPhredCut:]
print(vQualityPass)

## Calculate Number of Reads
vNumReads = fCountReads(vInput["R1"])
print("\nNumber of Reads Identified:", vNumReads)

## Identify Input Files
vR1 = gzip.open(vInput["R1"],"rt")
vR2 = gzip.open(vInput["R2"],"rt")
vR3 = gzip.open(vInput["R3"],"rt")
vR4 = gzip.open(vInput["R4"],"rt")
print("\nInput Files Identified:\n-",vInput["R1"],"\n-",vInput["R2"],"\n-",vInput["R3"],"\n-",vInput["R4"])

## Barcodes
vBarcodes = [
"GTAGCGTA","CGATCGAT","GATCAAGG","AACAGCGA","TAGCCATG","CGGTAATC",
"CTCTGGAT","TACCGGAT","CTAGCTCA","CACTTCAC","GCTACTCT","ACGATCAG",
"TATGGCAC","TGTTCCGT","GTCCTAAG","TCGACAAG","TCTTCGAC","ATCATGCG",
"ATCGTGGT","TCGAGAGT","TCGGATTC","GATCTTGC","AGAGTCCA","AGGATAGC"]

## Generate Output Filenames (vFiles)
vFiles = [
"output/Low_Quality_F.txt",
"output/Low_Quality_R.txt",
"output/Index_Mixed_F.txt",
"output/Index_Mixed_R.txt",
"output/summary.txt"]
for vCode in vBarcodes:
    vFiles.append("output/Index_Match_F_" + vCode + ".txt")
    vFiles.append("output/Index_Match_R_" + vCode + ".txt")
print("\nOutput files generated", "("+str(len(vFiles)-1)+")")

## Make Dictionary for Counting Index Combinations (vCountIndexes)
vCountIndexes = {"ambiguous":0, "low_quality":0, "unknown":0}
vCombos =  fFormatTuple(fGenPermutations(vBarcodes,2))
for vItem in vCombos:
    vCountIndexes[vItem] = 0



########## ########## ########## ########## ########## ##########
### MAIN PROGRAM



## Open all output files
with contextlib.ExitStack() as vStack:
    vOpen = [vStack.enter_context(open(vName,"w")) for vName in vFiles]

    ## Store Reference of Opened Files
    vOpened = {}
    for n, vDirectory in enumerate(vFiles):
        if "Match" in vDirectory:
            vOpened[vDirectory[19:29]] = vOpen[n]

    ## Iterate through each read in parallel across R1-R4 files
    for m in range(vNumReads):

        ## Extract Single Read
        vRead1  = list(itertools.islice(vR1,0,4))
        vIndex1 = list(itertools.islice(vR2,0,4))
        vIndex2 = list(itertools.islice(vR3,0,4))
        vRead2  = list(itertools.islice(vR4,0,4))

        ## Add barcodes to read header lines
        vSeqI1 = vIndex1[1].strip()
        vSeqI2 = fRevComp(vIndex2[1].strip())
        vBarComb = vSeqI1 + "-" + vSeqI2
        vRead1[0] = vRead1[0][:-1] + "\t" + vBarComb + "\n"
        vRead2[0] = vRead2[0][:-1] + "\t" + vBarComb + "\n"



        ## FILTER 1: Remove Ambiguous Indexes
        if "N" not in vIndex1[1] and "N" not in vIndex2[1]:
            
            ## FILTER 2: Remove Low Quality Indexes
            vQual_I1 = list(vIndex1[3].strip())
            vQual_I2 = list(vIndex2[3].strip())[::-1]
            if fIsSublist(vQualityPass,vQual_I1) and fIsSublist(vQualityPass,vQual_I2):
                
                ## FILTER 3: Is Index Hopping Present?
                ## Situation: Index hopping is not present
                if vBarComb in vCountIndexes:
                    if vSeqI1 == vSeqI2:
                        vCountIndexes[vBarComb] += 1
                        fWriteList(vRead1,vOpened["F_"+vSeqI1])
                        fWriteList(vRead2,vOpened["R_"+vSeqI1])

                    ## Situation: Index hopping is present
                    else:
                        if vBarComb in vCountIndexes:
                            vCountIndexes[vBarComb] += 1
                            fWriteList(vRead1,vOpen[2])
                            fWriteList(vRead2,vOpen[3])
                
                ## Situation: Unknown index combination (write to Low Quality files)
                else:
                    vCountIndexes["unknown"] += 1
                    fWriteList(vRead1,vOpen[0])
                    fWriteList(vRead2,vOpen[1])          
            
            ## Situation: Indexes are low quality
            else:
                vCountIndexes["low_quality"] += 1
                fWriteList(vRead1,vOpen[0])
                fWriteList(vRead2,vOpen[1])

        ## Situation: Indexes are ambiguous
        else:
            vCountIndexes["ambiguous"] += 1
            fWriteList(vRead1,vOpen[0])
            fWriteList(vRead2,vOpen[1])

        ## Track Progress
        fProgBar(m,vNumReads)



    ## Write Out: Dual Matched
    vSum = 0
    vOpen[4].write("Dual Matched Indexes\nIndex\tOccurences")
    for vItem in vBarcodes:
        vKey = vItem + "-" + vItem
        vOpen[4].write("\n"+vItem+"\t"+str(vCountIndexes[vKey]))
        vSum += vCountIndexes[vKey]
    vOpen[4].write("\nTOTAL\t"+str(vSum)+"\t("+str(vSum*100/vNumReads)+"%)")
 
    ## Write Out: Index Hopping
    vSum = 0
    vOpen[4].write("\n\nIndex-Hopping\nIndex1-Index2\tOccurences")
    for vItem in vCountIndexes:
        if "-" in vItem:
            if vItem.split("-")[0] != vItem.split("-")[1]:
                vOpen[4].write("\n"+vItem+"\t"+str(vCountIndexes[vItem]))
                vSum += vCountIndexes[vItem]
    vOpen[4].write("\nTOTAL\t"+str(vSum)+"\t("+str(vSum*100/vNumReads)+"%)")

    ## Write Out: Other
    vOpen[4].write("\n\nAmbiguous Index(s)\nIndex\tOccurences")
    vOpen[4].write("\nambiguous\t"+str(vCountIndexes["ambiguous"]) + "\t("+str(vCountIndexes["ambiguous"]*100/vNumReads)+"%)")
    vOpen[4].write("\nlow_quality\t"+str(vCountIndexes["low_quality"])+"\t("+str(vCountIndexes["low_quality"]*100/vNumReads)+"%)")
    vOpen[4].write("\nunknown\t"+str(vCountIndexes["unknown"])+"\t("+str(vCountIndexes["unknown"]*100/vNumReads)+"%)")