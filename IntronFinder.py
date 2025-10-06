#Note: '>' is used as a delimiter in multiple places. No particular reason for this choice of delimiter

#Importing packages: sys for system exit, re for regex search, argparse to setup command line input 
import sys
import re
import argparse

#User defined Functions
def isSingleSplitRead(cigar, alignNum):
    """
    Checks that a read is singly aligned and contains junctions using samfile information
    Samfile row for a singly aligned read ends in "NH:i:1"
    Cigar String for a read with junctions contains skipped regions denoted by 'N'
    :param cigar: Cigar string describing the alignment(String)
    :param alignNum: String of form NH:i:x where x is integer number of locations a read aligned(String)
    """
    if((alignNum == "NH:i:1") & (cigar.count("N") >= 1 )):
        return True
    else:
        return False

assert isSingleSplitRead("63M50N40M56N37M", "NH:i:1") == True
assert isSingleSplitRead("63M50N40M56N37M", "NH:i:2") == False #false since read aligned twice
assert isSingleSplitRead("63M", "NH:i:1") == False #false since no skipped regions

def geneLocProcessor(genomeLoc):
    """
    Returns a concatenated genome location information string after removing special characters
    and adding delimiter '>' between data points 
    :param genomeLoc: String containing chromosome, start, end positions and strand of a gene(String)
    """
    genomeLocation = genomeLoc.replace(",", "").replace(":", ">").replace("..", ">").replace("(", ">(")
    return genomeLocation
assert geneLocProcessor("TGME49_chrVIII:6,631,349..6,636,865(+)") == "TGME49_chrVIII>6631349>6636865>(+)"

def junctionCounter(rname, pos, cigar):
    """
    Returns a iterator of junctions in a given cigar string
    Traverses up the chromosome for each match or deletion in the cigar
    Records the start and end position of skipped regions
    Returns a concatenated string of Rname, start and end positions for each junction, delimited by '>'
    :param rname: Chromosome to which the read aligned(String)
    :param pos: start position of the alignment(int)
    :param cigar: A cigar string describing the alignment (String)
    """
    start = pos
    end = start
    cigParts = re.split("([A-Z])", cigar) #breaking cigar string into parts while retaining delimiters
    #Iterate over alternate parts since they are letters
    for juncCount in range(1, len(cigParts)-1, 2): 
        part = cigParts[juncCount]
        #Moves start position forward if matches or deletions are found
        if((part == "M") | (part == "D")): 
            start = start + int(cigParts[juncCount-1]) 
        #Records start and end position if junction is found
        elif(part == "N"): 
            end = start + int(cigParts[juncCount-1])
            junction = rname + ">" + str(start) + ">" + str(end) #concatenates info for returning
            start = end #moves to the end of junction and continues searching for more junctions
            yield junction

assert list(junctionCounter("TGME49_chrVIII", 6631349, "63M50N40M56N37M")) == ['TGME49_chrVIII>6631412>6631462', 'TGME49_chrVIII>6631502>6631558']

#Setting Up command line to accept the two input files using ArgParser
parser = argparse.ArgumentParser(description = "Junction Counting Assessment", formatter_class = argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("Sam", help = "Sam File") #accepts a sam file as first input
parser.add_argument("GeneLocFile", help = "Gene Locations File") #accepts a gene location table as second input
parameters = parser.parse_args()

#Assigning File Locations from argparser to input variables
Sam = parameters.Sam
GeneLocTable = parameters.GeneLocFile

Output = "Intron Junctions.txt" #Defining output table

#Extracting RNAME, POS and CIGAR from Sam file
AlignmentInfo = []
try:
    #Reading the Sam file
    with open(Sam) as SamFile: 
        for line in SamFile:
            #Skipping header rows starting with @
            if(line[0] == "@"): 
                continue
            try:
                SamInfo = line.strip().split() #split each row in the file into its constituents
                Rname = SamInfo[2]
                Pos = SamInfo[3]
                Cigar = SamInfo[5]
                AlignNum = SamInfo[len(SamInfo)-1]
                #Append a '>' delimited read information string if it is singly aligned and split
                if(isSingleSplitRead(Cigar, AlignNum)): 
                    AlignmentInfo.append(Rname + ">" + Pos + ">" + Cigar)
            #In case of missing information skip the line
            except IndexError: 
                print(line + " in samfile is lacking essential information. Line skipped")
                continue
#Exit with error message if file isnt found
except FileNotFoundError: 
    sys.exit(Sam + " not found. Please check and try again")

#Extracting Gene Locations from gene summary file
GeneInfo = {}
try:
    with open(GeneLocTable) as GeneLocFile:
        Header = GeneLocFile.readline() #store header line separately
        for line in GeneLocFile:
            GenePart = line.strip().split() #split each row in the file into its constituents
            try:
                GeneKey = geneLocProcessor(GenePart[2]) #get a > delimited string with gene info
                GeneInfo[GeneKey] = GenePart[0] #link gene info to its gene ID
            #In case of missing information skip the line
            except IndexError: 
                print(line + " in Gene Location file is lacking essential information. Line skipped")
                continue
#Exit with error message if file isnt found
except FileNotFoundError: 
    sys.exit(GeneLocTable + " not found. Please check and try again")
    
#Counting Junctions
JunctionInfo = {}
for line in AlignmentInfo:
    Rname, Pos, Cigar = line.split(">") #extract alignment info from > delimited string
    JunctionSet = junctionCounter(Rname, int(Pos), Cigar) #stores junctions in a single read
    for Junction in JunctionSet:
        #Update count of junctions found repeatedly
        try:
            JunctionInfo[Junction] += 1 
        #Add new entry if new junction is found 
        except KeyError:
            JunctionInfo[Junction] = 1 

#Mapping junctions to their genes
GeneJunctionTable = []
GeneJunctionTable.append("Gene ID\tJuncStart\tJuncEnd\tSupportReads\n") #adding a header row
for line in GeneInfo.keys():
    GeneID= GeneInfo[line] #extract GeneID from dict
    GeneJunctionCounter = 0 #counts Junctions per gene
    GeneChrom, GeneStart, GeneEnd, Strand = line.split(">") #extract gene info from > delimited string 
    for row in JunctionInfo.keys():
        JunctionCount = str(JunctionInfo[row]) #extract number of junctions at specific loc from dict
        JuncChrom, JuncStart, JuncEnd = row.split(">") #extract junction info from > delimited string
        #Append Junctions lying in gene boundries
        if((int(JuncStart) in range( int(GeneStart), int(GeneEnd) + 1, 1)) & (int(JuncEnd) in range( int(GeneStart), int(GeneEnd) + 1, 1))): 
            GeneJunctionTable.append(GeneID + "\t" + JuncStart + "\t" + JuncEnd + "\t" + JunctionCount)
            GeneJunctionCounter += 1 #update GeneJunctionCounter
    #if the previous gene had junctions then add a blankspace which will later become a newline 
    if(GeneJunctionCounter > 0):
        GeneJunctionTable.append("") 

#Writing output table and adding newline separators
with(open(Output, 'w') as out): 
     out.write('\n'.join(GeneJunctionTable))