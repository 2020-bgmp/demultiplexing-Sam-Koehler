###DEF: Import
# conda install numpy
# conda install matplotlib
import gzip
import matplotlib.pyplot as plt
import numpy as np

###DEF: Input File
# vFile = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
# vFile = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
vFile = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
# vFile = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"
# vFile = "2_test.txt.gz"


###DEF: Define Function to Convert Phred+33 Scores
def convert_phred(letter):
    """Converts a single character into a phred score"""
    return(ord(letter)-33)



###DEF: Generate List of Lists
vWidth = 8
all_qscores = np.zeros(vWidth,dtype=float)



###DEF: Populate all_qscores with all quality scores at each index
## Tracks Read Number (0-vLength)
vReadCounter = 0

with gzip.open(vFile,'rt') as vOpen:
    for n,vLine in enumerate(vOpen):
        if n % 4 == 3:
            for m, vLetter in enumerate(vLine.strip()):
                all_qscores[m] += convert_phred(vLetter)
            vReadCounter += 1
            
        ###DEF: Progress Counter
        vRange = 14600000
        if n % vRange == 0:
            print("~",n/vRange,"% complete")

###DEF: Convert sums to means
for vIndex in range(len(all_qscores)):
    all_qscores[vIndex] = all_qscores[vIndex] / vReadCounter



###DEF: Plot K-mer occurence counts
print("Initialize: Graph")
###DEF: Define Figure, specify size
vFigure = plt.figure(figsize=(8, 4))
###DEF: Defines Axes
vAxes = vFigure.add_axes([0,0,1,1])
###DEF: Populate Axes with Data Lists, Standard Deviation (yerr)
vAxes.bar([x for x in range(vWidth)],all_qscores)
###DEF: Y axis is logarithmic
###DEF: Label Axes
plt.xlabel("nt Position")
plt.ylabel("Mean Quality Score")
###DEF: Display Figure
plt.savefig("R3_output.png",bbox_inches='tight')