from math import log
import sys

from ScaffoldSeq import *
from userinterface import *

##############################
###   Start of Execution   ###
##############################

'''
Silent Start Detector
Purpose: If a command line argument is provided, the program will run in silence
Notes:
'''
silent_mode = False
if len(sys.argv) == 1:
    silent_mode = False
elif len(sys.argv) == 2:  # Triggers Interface Mode
    print('Silent mode started. Interactive mode will be circumvented.')
    if os.path.isfile(sys.argv[1]):
        job_file = sys.argv[1]
        silent_mode = True
    else:
        print('Job file dictated by command line argument not located. Please check file name.')
        getch()
elif len(sys.argv) > 2:
    print('Too many input arguments. Please revise command line arguments.')
    print('Press any key to exit.')
    input()
    raise SystemExit

start = time.time()  # Starts timer


'''
Program Start
Purpose: Start the application, choosing between silent and non-silent modes
Note: You can change the default settings below if so desired
'''

if silent_mode:
    settings, job, divs = Silent_Load(job_file)

else:
    clear()
    settings = ['.8', '0.5', '10000', 'On', 'On', '10', 'Return to Main Menu    ']
    job, divs, settings = main(settings)

[thresh, damp, maxSeqCount, use_bck, pairwise, bead_ratio] = settings[0:6]
[thresh, damp, maxSeqCount, use_bck, pairwise, bead_ratio] = [float(thresh), float(damp), int(maxSeqCount), use_bck,
                                                              pairwise, float(bead_ratio)]
[filename, files], FRAMEsplit = job[0:2], []
adaptors = job[3:5]
FRAMEsplit.append(job[2])
looplength, diversityspacing = [], []
for i in range(job[5]):
    looplength.append(divs[i][1:3])
    diversityspacing.append(divs[i][3])
    FRAMEsplit.append(divs[i][0])

# Generate empty count matrices for later counting
LOOPcnt = [[[[0. * 22] for j in range(i)] for i in range(pos[0], pos[1] + 1)] for pos in looplength]

# Set up variables to be used for sequence sorting


'''
Program Start
'''

loadtime = time.time() - start
print('Job has started. \nScaffoldSeq is evaluating the data set...')
time1 = time.time()

ProteinTotal, FULLPROT, timeCnt = 0, [], []
timeCnt.append(time.time())

'''
Sequence Acquisition (Ruler Generation/Search)
Purpose: Creates a 'blank scaffold' of a certain size, uses the 5'/3' 
    boundaries to identify sequences that fit the size restrictions

FULLPROT --> Stores list of all diversified regions originating from the same 
    sequence, but only when all diversified regions are accounted for.

'''
while not os.path.isfile(files):
    files = input(
        '\nError! \nSequence data file must exist within the current directory.' +
        '\nPlease enter the sequence file name: \n')

with open(files, 'r') as Mainfile:
    LOWlmt = len(''.join(FRAMEsplit)) + 3 * (sum([i[0] for i in looplength]))
    HIGHlmt = len(''.join(FRAMEsplit)) + 3 * (sum([i[1] for i in looplength]))
    SeqCnt = 0
    LOOPseqs = [[[] for i in range(j[1] - j[0] + 1)] for j in looplength]

    for line in Mainfile:
        SeqCnt += 1
        matchFn, tempLOOPS, buildFULL = 0, [], []
        if SeqCnt >= maxSeqCount:
            break
        length, protein = Adaptor_Alignment(line, adaptors)
        if LOWlmt <= length <= HIGHlmt and length % 3 == 0:
            LOOPstart = 0
            # Process sequences that lie within an appropriate size range
            '''
            Sequence Acquisition (Region Alignment)
            Purpose: Goes through each matching sequence at tries to align sections of framework in order to identify
             diversified regions
            Note:
            '''
            for (i, loop) in zip(range(len(looplength)), looplength):
                protFrag, protFragRat, protFragaa = [], [], []

                LOOPstart += len(FRAMEsplit[i]) + loop[0] * 3

                for x in range(LOOPstart, LOOPstart + (loop[1] - loop[0] + 1) * 3, 3):
                    protFrag.append(protein[x:x + len(FRAMEsplit[i + 1])])

                for x in protFrag:
                    protFragRat.append(FamCrit(FRAMEsplit[i + 1], x))

                tag = protFragRat.index(max(protFragRat))
                LOOPseq = protein[LOOPstart - loop[0] * 3:LOOPstart + tag * 3]

                # global framework_match_threshold
                if max(protFragRat) >= framework_match_threshold:
                    matchFn += 1
                    LOOPseqs[i][tag].append(translate_dna(LOOPseq))
                    tempLOOPS.append(translate_dna(LOOPseq))

                # DRW edit on 2016-02-01 to exclude all cysteine containing clones.
                #                if max(protFragRat) >= framework_match_threshold and 'C' not in translate_dna(LOOPseq):
                #                    matchFn += 1
                #                    LOOPseqs[i][tag].append(translate_dna(LOOPseq))
                #                    tempLOOPS.append(translate_dna(LOOPseq))

                LOOPstart += 3 * tag  # Reset loopstart variable to prepare for next interation

                # Take into account loop length diversity (-1...+1)
        if len(tempLOOPS) == len(looplength):
            FULLPROT.append(tempLOOPS)

        if matchFn > 0:  # Keep track of the number of identified proteins of interest
            ProteinTotal += 1

        if SeqCnt % 100000 == 0:  # Give us some kind of meter for data processing
            print('Scanned  next 100k entries in %.1f sec' % (time.time() - timeCnt[0]))

# Remind user how long the algorithm took to identify all relavent proteins
scantime = time.time() - time1
print('Scanned full data set in %.1f sec' % (scantime))
time1b = time.time()

"""
###################################################################
##                                                               ##
##                   Operate on AA sequences                     ##
##      List, Tally  Occurrences and Sort Unique Sequences       ##
##                                                               ##
###################################################################
"""

# Enumerate and Remove Duplicates by Loop Length

UNIQUELOOPseqs = [[[[str(x), LOOPseq.count(x)] for x in list(set(LOOPseq))] for LOOPseq in LOOPset] for LOOPset in
                  LOOPseqs]
UNIQUELOOPseqs = [[sorted(LOOPseq, key=lambda x: -x[1]) for LOOPseq in LOOPset] for LOOPset in UNIQUELOOPseqs]

###################################
##           Time Point          ##
###################################

scantime1 = time.time() - time1b
print('Organized Diversified Regions in %.1f sec' % scantime1)
time2 = time.time()

#######################################
##  Remove Background (Rare Events)  ##
#######################################

FREQcntr = [[i[1] for FREQset in UNIQUELOOPseq for i in FREQset] for UNIQUELOOPseq in UNIQUELOOPseqs]
seqcnt = [[[x, cntr.count(x)] for x in list(set(cntr))] for cntr in
          FREQcntr]  # Removes repeats and stores as [sequence,frequency] instead
seqcnt = [sorted(sc, key=lambda x: -x[0]) for sc in seqcnt]  # Order based on frequency
cntset = [list(set(cntr)) for cntr in FREQcntr]
cntset = [sorted(cs) for cs in cntset]  # reverse=True   for reverse ordering
BGcap = []

'''
Sequence Analysis (magnetic bead sorting background, see Ackerman, et al. Biotechnol. Prog., 2009)
Purpose: Remove sequence background specific to bead sorting assay. 
Note: Can be turned off via the settings file
'''

# This value defines the threshold level for inclusion of sequences throughout the analysis.
# Sequences having fewer occurrences than the 'bck' threshold will be considered rare, unwanted false-positives.


if use_bck == 'On':
    bck = 0.25 / bead_ratio
    cut = [(sum(cntr) * bck) + ((sum(cntr) * bck) ** (0.5)) for cntr in FREQcntr]
    for cs, c, cntr in zip(cntset, cut, FREQcntr):
        setmark = 0.
        for i in cs:
            setmark += cntr.count(i) * i
            if setmark >= c:
                BGcap.append(i)
                break
else:
    BGcap = [bead_ratio for i in cntset]

pairBG = float(min(BGcap))

UNIQUELOOPseqs = [[[i for i in FREQset if i[1] > BG] for FREQset in UNIQUELOOPseq] for (UNIQUELOOPseq, BG) in
                  zip(UNIQUELOOPseqs, BGcap)]

print('Total Proteins', str(ProteinTotal))

looptime = time.time() - time2
print('Background Removed in %.2f sec' % (looptime))
time3 = time.time()

###################################
##       Family Clustering       ##
###################################

print('Clustering Threshold : ', thresh)

'''
Sequence Analysis (Gap alignments)
Purpose: Adds dashes to the clustered families in order to align multiple lengths of loops
Note:
'''

'''
Sequence Analysis (Clusters sequences)
Purpose: Checks for sequence homology and if above clustering threshold, stores as new cluster
Note: 
'''
if damp == 1.0:
    pass
# else:
# print 'UNIQUELOOPseqs', UNIQUELOOPseqs


CLUSTERseqs = []
fam1cnt = 1
for LOOPseqs, loop, space in zip(UNIQUELOOPseqs, looplength, diversityspacing):
    # print 'LOOPseqs', LOOPseqs[0]
    fam1cnt, FamLst1 = 1, []
    for B in range(len(LOOPseqs)):
        if len(LOOPseqs[B]) == 0:
            FamLst1.append([[]])

        else:
            FamLst1.append([])
            FamLst1[B].append([LOOPseqs[B][0]])
            LOOPseqs[B].remove(LOOPseqs[B][0])

        for ent in LOOPseqs[B]:
            skip = 0
            for zf in FamLst1[B]:
                if FamCrit(ent[0], zf[0][0]) == 1.0:
                    skip = 1
                elif skip == 0 and damp == 1 and thresh == 0:
                    zf.append(ent)
                    skip = 1
                elif skip == 0 and FamCrit(ent[0], zf[0][0]) >= thresh:
                    zf.append(ent)
                    skip = 1
                else:
                    if skip == 0 and FamLst1[B].index(zf) == len(FamLst1[B]) - 1:
                        FamLst1[B].append([ent])
                        fam1cnt += 1

    CLUSTERseqs.append(add_gap(FamLst1, loop[1], space))
    # print 'CLUSTERseqs', CLUSTERseqs
# print 'CLUSTERseqs', CLUSTERseqs
# Creates a flattened list for these clustered sequences
# C_seqs = [d for a in CLUSTERseqs for b in a for c in b for d in c]


unique_time = time.time() - time3
print('Family Clusters Identified in %.2f sec' % unique_time)
time3b = time.time()

"""
###################################################################
##           Quantify Normalized Amino Acid Frequency            ##
##      Output Includes Positional AA Frequency and Occurrence   ##
###################################################################
"""

# Counts occurrences of amino acids at each position
FREQUENCYseqs = [residue_frequency(LOOPseq, loop, damp)[0] for LOOPseq, loop in zip(CLUSTERseqs, looplength)]
#
LEADSseqs = [residue_frequency(LOOPseq, loop, damp)[1] for LOOPseq, loop in zip(CLUSTERseqs, looplength)]
# Calculates the fractional distribution of amino acids at each position
PERCENTseqs = [[[i / sum(j) for i in j] for j in loop] for loop in FREQUENCYseqs]

# print PERCENTseqs

###

sorttime = time.time() - time3b
print('Site-wise Frequency Matrix Constructed in %.2f sec' % sorttime)
time4 = time.time()
freqmatrix = [a for b in PERCENTseqs for a in b]

###

if pairwise == 'On':
    # Add place holders for loop length diversity alignment for the list of
    # intact diversified regions (FULLPROT) originating from the same sequence
    seq = full_format(FULLPROT, looplength, diversityspacing)

    # Counts and stores all unique protein sequences in single list --> [[seq0, 0#],[seq1, 1#]...]
    seq = full_freq(seq)

    print('Analyzing Pairwise Interactions Within %.i Unique Proteins ' % len(seq))

    Unique_Prot = [len(seq)]

    pairdata, PERCENTpairseqs = full_pairwise(seq, pairBG, damp)

    results, freqdata, refs = FullMatrixCompare(pairdata, PERCENTpairseqs)

    pairtime = time.time() - time4
    print('Pairwise Matrices Constructed in %.2f sec' % pairtime)
    time4 = time.time()

'''
Additional homology models can be inserted below using matching formatting.
'''

# Comparisons involving gaps are assigned scores that match the X (any amino acid) category within RBLOSUM64.
# Score associated with Stop Codons are given a score of -9.
bl_64 = {('A', 'A'): 4, ('A', 'R'): -2, ('A', 'N'): -1, ('A', 'D'): -2, ('A', 'C'): -1, ('A', 'Q'): -1, ('A', 'E'): -1,
         ('A', 'G'): 0, ('A', 'H'): -2, ('A', 'I'): -1, ('A', 'L'): -2, ('A', 'K'): -1, ('A', 'M'): -1, ('A', 'F'): -2,
         ('A', 'P'): -1, ('A', 'S'): 1, ('A', 'T'): 0, ('A', 'W'): -3, ('A', 'Y'): -2, ('A', 'V'): 0, ('A', '_'): -9,
         ('A', '-'): -1,
         ('R', 'A'): -2, ('R', 'R'): 5, ('R', 'N'): 0, ('R', 'D'): -2, ('R', 'C'): -3, ('R', 'Q'): 1, ('R', 'E'): 0,
         ('R', 'G'): -2, ('R', 'H'): 0, ('R', 'I'): -3, ('R', 'L'): -2, ('R', 'K'): 2, ('R', 'M'): -2, ('R', 'F'): -3,
         ('R', 'P'): -2, ('R', 'S'): -1, ('R', 'T'): -1, ('R', 'W'): -2, ('R', 'Y'): -2, ('R', 'V'): -3, ('R', '_'): -9,
         ('R', '-'): -1,
         ('N', 'A'): -1, ('N', 'R'): 0, ('N', 'N'): 6, ('N', 'D'): 1, ('N', 'C'): -3, ('N', 'Q'): 0, ('N', 'E'): 0,
         ('N', 'G'): -1, ('N', 'H'): 1, ('N', 'I'): -3, ('N', 'L'): -3, ('N', 'K'): 0, ('N', 'M'): -2, ('N', 'F'): -3,
         ('N', 'P'): -2, ('N', 'S'): 0, ('N', 'T'): 0, ('N', 'W'): -3, ('N', 'Y'): -2, ('N', 'V'): -3, ('N', '_'): -9,
         ('N', '-'): -1,
         ('D', 'A'): -2, ('D', 'R'): -2, ('D', 'N'): 1, ('D', 'D'): 6, ('D', 'C'): -3, ('D', 'Q'): 0, ('D', 'E'): 2,
         ('D', 'G'): -2, ('D', 'H'): -1, ('D', 'I'): -3, ('D', 'L'): -3, ('D', 'K'): -1, ('D', 'M'): -3, ('D', 'F'): -3,
         ('D', 'P'): -2, ('D', 'S'): 0, ('D', 'T'): -1, ('D', 'W'): -4, ('D', 'Y'): -3, ('D', 'V'): -3, ('D', '_'): -9,
         ('D', '-'): -1,
         ('C', 'A'): -1, ('C', 'R'): -3, ('C', 'N'): -3, ('C', 'D'): -3, ('C', 'C'): 9, ('C', 'Q'): -3, ('C', 'E'): -4,
         ('C', 'G'): -3, ('C', 'H'): -2, ('C', 'I'): -1, ('C', 'L'): -1, ('C', 'K'): -3, ('C', 'M'): -1, ('C', 'F'): -2,
         ('C', 'P'): -3, ('C', 'S'): -1, ('C', 'T'): -1, ('C', 'W'): -3, ('C', 'Y'): -2, ('C', 'V'): -1, ('C', '_'): -9,
         ('C', '-'): -2,
         ('Q', 'A'): -1, ('Q', 'R'): 1, ('Q', 'N'): 0, ('Q', 'D'): 0, ('Q', 'C'): -3, ('Q', 'Q'): 5, ('Q', 'E'): 2,
         ('Q', 'G'): -2, ('Q', 'H'): 1, ('Q', 'I'): -3, ('Q', 'L'): -2, ('Q', 'K'): 1, ('Q', 'M'): 0, ('Q', 'F'): -3,
         ('Q', 'P'): -1, ('Q', 'S'): 0, ('Q', 'T'): 0, ('Q', 'W'): -2, ('Q', 'Y'): -2, ('Q', 'V'): -2, ('Q', '_'): -9,
         ('Q', '-'): -1,
         ('E', 'A'): -1, ('E', 'R'): 0, ('E', 'N'): 0, ('E', 'D'): 2, ('E', 'C'): -4, ('E', 'Q'): 2, ('E', 'E'): 5,
         ('E', 'G'): -2, ('E', 'H'): 0, ('E', 'I'): -3, ('E', 'L'): -3, ('E', 'K'): 1, ('E', 'M'): -2, ('E', 'F'): -3,
         ('E', 'P'): -1, ('E', 'S'): 0, ('E', 'T'): -1, ('E', 'W'): -3, ('E', 'Y'): -2, ('E', 'V'): -2, ('E', '_'): -9,
         ('E', '-'): -1,
         ('G', 'A'): 0, ('G', 'R'): -2, ('G', 'N'): -1, ('G', 'D'): -2, ('G', 'C'): -3, ('G', 'Q'): -2, ('G', 'E'): -2,
         ('G', 'G'): 6, ('G', 'H'): -2, ('G', 'I'): -4, ('G', 'L'): -4, ('G', 'K'): -2, ('G', 'M'): -3, ('G', 'F'): -3,
         ('G', 'P'): -2, ('G', 'S'): -1, ('G', 'T'): -2, ('G', 'W'): -2, ('G', 'Y'): -3, ('G', 'V'): -3, ('G', '_'): -9,
         ('G', '-'): -2,
         ('H', 'A'): -2, ('H', 'R'): 0, ('H', 'N'): 1, ('H', 'D'): -1, ('H', 'C'): -2, ('H', 'Q'): 1, ('H', 'E'): 0,
         ('H', 'G'): -2, ('H', 'H'): 7, ('H', 'I'): -3, ('H', 'L'): -3, ('H', 'K'): -1, ('H', 'M'): -1, ('H', 'F'): -1,
         ('H', 'P'): -2, ('H', 'S'): -1, ('H', 'T'): -2, ('H', 'W'): -1, ('H', 'Y'): 1, ('H', 'V'): -3, ('H', '_'): -9,
         ('H', '-'): -1,
         ('I', 'A'): -1, ('I', 'R'): -3, ('I', 'N'): -3, ('I', 'D'): -3, ('I', 'C'): -1, ('I', 'Q'): -3, ('I', 'E'): -3,
         ('I', 'G'): -4, ('I', 'H'): -3, ('I', 'I'): 4, ('I', 'L'): 2, ('I', 'K'): -3, ('I', 'M'): 1, ('I', 'F'): 0,
         ('I', 'P'): -3, ('I', 'S'): -2, ('I', 'T'): -1, ('I', 'W'): -2, ('I', 'Y'): -1, ('I', 'V'): 2, ('I', '_'): -9,
         ('I', '-'): -1,
         ('L', 'A'): -2, ('L', 'R'): -2, ('L', 'N'): -3, ('L', 'D'): -3, ('L', 'C'): -1, ('L', 'Q'): -2, ('L', 'E'): -3,
         ('L', 'G'): -4, ('L', 'H'): -3, ('L', 'I'): 2, ('L', 'L'): 4, ('L', 'K'): -2, ('L', 'M'): 2, ('L', 'F'): 1,
         ('L', 'P'): -3, ('L', 'S'): -2, ('L', 'T'): -1, ('L', 'W'): -1, ('L', 'Y'): -1, ('L', 'V'): 1, ('L', '_'): -9,
         ('L', '-'): -1,
         ('K', 'A'): -1, ('K', 'R'): 2, ('K', 'N'): 0, ('K', 'D'): -1, ('K', 'C'): -3, ('K', 'Q'): 1, ('K', 'E'): 1,
         ('K', 'G'): -2, ('K', 'H'): -1, ('K', 'I'): -3, ('K', 'L'): -2, ('K', 'K'): 5, ('K', 'M'): -1, ('K', 'F'): -3,
         ('K', 'P'): -1, ('K', 'S'): 0, ('K', 'T'): -1, ('K', 'W'): -3, ('K', 'Y'): -2, ('K', 'V'): -2, ('K', '_'): -9,
         ('K', '-'): -1,
         ('M', 'A'): -1, ('M', 'R'): -2, ('M', 'N'): -2, ('M', 'D'): -3, ('M', 'C'): -1, ('M', 'Q'): 0, ('M', 'E'): -2,
         ('M', 'G'): -3, ('M', 'H'): -1, ('M', 'I'): 1, ('M', 'L'): 2, ('M', 'K'): -1, ('M', 'M'): 6, ('M', 'F'): 0,
         ('M', 'P'): -2, ('M', 'S'): -1, ('M', 'T'): -1, ('M', 'W'): -2, ('M', 'Y'): -1, ('M', 'V'): 0, ('M', '_'): -9,
         ('M', '-'): -1,
         ('F', 'A'): -2, ('F', 'R'): -3, ('F', 'N'): -3, ('F', 'D'): -3, ('F', 'C'): -2, ('F', 'Q'): -3, ('F', 'E'): -3,
         ('F', 'G'): -3, ('F', 'H'): -1, ('F', 'I'): 0, ('F', 'L'): 1, ('F', 'K'): -3, ('F', 'M'): 0, ('F', 'F'): 6,
         ('F', 'P'): -3, ('F', 'S'): -2, ('F', 'T'): -2, ('F', 'W'): 1, ('F', 'Y'): 3, ('F', 'V'): -1, ('F', '_'): -9,
         ('F', '-'): -1, ('P', 'A'): -1, ('P', 'R'): -2, ('P', 'N'): -2, ('P', 'D'): -2, ('P', 'C'): -3,
         ('P', 'Q'): -1, ('P', 'E'): -1, ('P', 'G'): -2, ('P', 'H'): -2, ('P', 'I'): -3, ('P', 'L'): -3, ('P', 'K'): -1,
         ('P', 'M'): -2, ('P', 'F'): -3, ('P', 'P'): 7, ('P', 'S'): -1, ('P', 'T'): -1, ('P', 'W'): -3, ('P', 'Y'): -3,
         ('P', 'V'): -2, ('P', '_'): -9, ('P', '-'): -2,
         ('S', 'A'): 1, ('S', 'R'): -1, ('S', 'N'): 0, ('S', 'D'): 0, ('S', 'C'): -1, ('S', 'Q'): 0, ('S', 'E'): 0,
         ('S', 'G'): -1, ('S', 'H'): -1, ('S', 'I'): -2, ('S', 'L'): -2, ('S', 'K'): 0, ('S', 'M'): -1, ('S', 'F'): -2,
         ('S', 'P'): -1, ('S', 'S'): 4, ('S', 'T'): 1, ('S', 'W'): -3, ('S', 'Y'): -2, ('S', 'V'): -2, ('S', '_'): -9,
         ('S', '-'): -1, ('T', 'A'): 0, ('T', 'R'): -1, ('T', 'N'): 0, ('T', 'D'): -1, ('T', 'C'): -1, ('T', 'Q'): 0,
         ('T', 'E'): -1,
         ('T', 'G'): -2, ('T', 'H'): -2, ('T', 'I'): -1, ('T', 'L'): -1, ('T', 'K'): -1, ('T', 'M'): -1, ('T', 'F'): -2,
         ('T', 'P'): -1, ('T', 'S'): 1, ('T', 'T'): 5, ('T', 'W'): -3, ('T', 'Y'): -2, ('T', 'V'): 0, ('T', '_'): -9,
         ('T', '-'): -1, ('W', 'A'): -3, ('W', 'R'): -2, ('W', 'N'): -3, ('W', 'D'): -4, ('W', 'C'): -3, ('W', 'Q'): -2,
         ('W', 'E'): -3,
         ('W', 'G'): -2, ('W', 'H'): -1, ('W', 'I'): -2, ('W', 'L'): -1, ('W', 'K'): -3, ('W', 'M'): -2, ('W', 'F'): 1,
         ('W', 'P'): -3, ('W', 'S'): -3, ('W', 'T'): -3, ('W', 'W'): 11, ('W', 'Y'): 2, ('W', 'V'): -3, ('W', '_'): -9,
         ('W', '-'): -2, ('Y', 'A'): -2, ('Y', 'R'): -2, ('Y', 'N'): -2, ('Y', 'D'): -3, ('Y', 'C'): -2, ('Y', 'Q'): -2,
         ('Y', 'E'): -2,
         ('Y', 'G'): -3, ('Y', 'H'): 1, ('Y', 'I'): -1, ('Y', 'L'): -1, ('Y', 'K'): -2, ('Y', 'M'): -1, ('Y', 'F'): 3,
         ('Y', 'P'): -3, ('Y', 'S'): -2, ('Y', 'T'): -2, ('Y', 'W'): 2, ('Y', 'Y'): 7, ('Y', 'V'): -1, ('Y', '_'): -9,
         ('Y', '-'): -1, ('V', 'A'): 0, ('V', 'R'): -3, ('V', 'N'): -3, ('V', 'D'): -3, ('V', 'C'): -1, ('V', 'Q'): -2,
         ('V', 'E'): -2,
         ('V', 'G'): -3, ('V', 'H'): -3, ('V', 'I'): 2, ('V', 'L'): 1, ('V', 'K'): -2, ('V', 'M'): 0, ('V', 'F'): -1,
         ('V', 'P'): -2, ('V', 'S'): -2, ('V', 'T'): 0, ('V', 'W'): -3, ('V', 'Y'): -1, ('V', 'V'): 4, ('V', '_'): -9,
         ('V', '-'): -1, ('_', 'A'): -9, ('_', 'R'): -9, ('_', 'N'): -9, ('_', 'D'): -9, ('_', 'C'): -9, ('_', 'Q'): -9,
         ('_', 'E'): -9,
         ('_', 'G'): -9, ('_', 'H'): -9, ('_', 'I'): -9, ('_', 'L'): -9, ('_', 'K'): -9, ('_', 'M'): -9, ('_', 'F'): -9,
         ('_', 'P'): -9, ('_', 'S'): -9, ('_', 'T'): -9, ('_', 'W'): -9, ('_', 'Y'): -9, ('_', 'V'): -9, ('_', '_'): -9,
         ('_', '-'): -9, ('-', 'A'): -1, ('-', 'R'): -1, ('-', 'N'): -1, ('-', 'D'): -1, ('-', 'C'): -2, ('-', 'Q'): -1,
         ('-', 'E'): -1,
         ('-', 'G'): -2, ('-', 'H'): -1, ('-', 'I'): -1, ('-', 'L'): -1, ('-', 'K'): -1, ('-', 'M'): -1, ('-', 'F'): -1,
         ('-', 'P'): -2, ('-', 'S'): -1, ('-', 'T'): -1, ('-', 'W'): -2, ('-', 'Y'): -1, ('-', 'V'): -1, ('-', '_'): -9,
         ('-', '-'): 0
         }

#   Showing similarity as large value through the addition on ("len(a) - ...")


LEADdistances, LEADdistancesBL = [], []
# Top 100 Lead Clones are used for similarity distance calculations
for LEADset in LEADSseqs:
    LEADS = [i[0][0] for i in LEADset]
    LEADdistance, LEADdistanceBL = [[''] + LEADS], [[''] + LEADS]
    for seq in LEADS[:100]:
        LEADdistance.append([seq])
        LEADdistanceBL.append([seq])
        for seq2 in LEADS[:100]:
            LEADdistance[-1].append(quick_diff(seq, seq2))

            if sum([max([v for k, v in bl_64.items() for s in (i, j) if s in k]) for i, j in zip(seq, seq2)]) == 0:
                LEADdistanceBL[-1].append(0)
                # print seq,seq2
            else:
                LEADdistanceBL[-1].append(float(sum([bl_64[(i, j)] for i, j in zip(seq, seq2)])) / (
                    sum([max([v for [k, v] in bl_64.items() for s in (i, j) if s in k]) for i, j in zip(seq, seq2)])))
    LEADdistances.append(LEADdistance)
    LEADdistancesBL.append(LEADdistanceBL)

phylo_time = time.time() - time4
print('Similarity Matrix Constructed in %.2f sec' % phylo_time)
time4b = time.time()
# print LEADdistances

# Extensive Data Set for Associations
if pairwise == 'On':
    publishdata = open(filename + '_Epistasis.csv', 'w')
    writer1 = csv.writer(publishdata)

    # Job Overview
    overview = []
    for i in range(len(aa)):
        for j in range(len(aa)):
            overview.append(aa[i] + aa[j])
    writer1.writerow(['All Pairwise interactions:'])
    writer1.writerow(['Alignment', 'Positions', 'Residues', 'Predicted Frequency', 'Pairwise Frequency', 'Difference',
                      'Epistasis (RSE)'])

    MIList = []

    for f1, p2, ref in zip(freqdata, pairdata, refs):

        MIcoord = 0.

        for a, b, c, i in zip([ii for jj in f1 for ii in jj], [iii for jjj in p2 for iii in jjj], overview,
                              range(len(overview))):
            if float(a) != 0. and float(b) != 0.: MIcoord += float(b) * log(float(b) / float(a), 2)

            # If file size is not an issue, remove this conditional before next statement: if float(a) != 0. and float(b) != 0.:
            if float(a) != 0. and float(b) != 0.: writer1.writerow(ref + ["'" + c] + [a] + [b] + [abs(a - b)] + [
                float(b) * log(float(b) / float(a), 2)])  # Space saving (this data file gets massive)

        MIList.append([ref[1], MIcoord])

    publishdata.close()

    # Build MI APC Table (Dunn, et al, Bioinformatics 2008)
    MI_APC_List = []

    for pair in MIList:
        MI_APC_List.append([int(pair[0][:pair[0].index("&")]), int(pair[0][pair[0].index("&") + 1:]), pair[1]])

        MI_APC_Set = set([x[0] for x in MI_APC_List] + [y[1] for y in MI_APC_List])

        # Use this set to guide APC calc.
        # MI_APC = [sum(x[1] for x in MIList if (int(x[0][:x[0].index("&")])== b or int(x[0][x[0].index("&")+1:]) == b)) for b in MI_APC_Set]

        MI_APC = [sum(x[2] for x in MI_APC_List if (int(x[0]) == b or int(x[1]) == b)) / (len(MI_APC_Set) - 1) for b in
                  MI_APC_Set]

        MI_avg = sum(MI_APC) / len(MI_APC)

        APC = [x[2] - (MI_APC[x[0] - 1] * MI_APC[x[1] - 1] / MI_avg) for x in MI_APC_List]

    MI_APC_MERGE = []
    for entry in enumerate(MI_APC_List):
        MI_APC_MERGE.append(entry[1])
        MI_APC_MERGE[entry[0]].append(APC[entry[0]])

    MItime = time.time() - time4b
    # print 'Mutual Information Values Calculated in %.2f sec' % (MItime)

    time4c = time.time()

    publishMIdata = open(filename + '_Mutual-Information.csv', 'w')
    writerMI = csv.writer(publishMIdata)

    # Frequency of times each residue occurs
    writerMI.writerow(['Total Proteins', str(ProteinTotal)])
    writerMI.writerow(['Unique Proteins', '  %.i  ' % (Unique_Prot[0])])
    writerMI.writerow(['Total Clusters', fam1cnt])
    # writerMI.writerow(['Organized Diversified Regions (sec)', ' %.1f ' % (scantime1)])
    # writerMI.writerow(['Scanned Full Data Set (sec)', ' %.1f ' % (scantime)])
    # writerMI.writerow(['Family Clustering (sec)', ' %.2f ' % (uniquetime)])
    # writerMI.writerow(['Total Runtime (sec)', ' %.0f ' % (time.time()-time1)])
    # writerMI.writerow(['Residue Frequency:'])
    # writerMI.writerow(list(aa))
    # for row in PERCENTpairseqs:
    # writerMI.writerow(row)
    writerMI.writerow('\n')

    writerMI.writerow(
        ['Site 1', 'Site 2', 'Mutual Information (MI) ', 'Mutual Information with Average Product Correction (MIp)'])
    for row in MI_APC_MERGE:
        writerMI.writerow(row)

    #    writerMI.writerow('\n')
    #
    #   writerMI.writerow(['Position Pair', 'Mutual Information', 'Mutual Information with Average Product Correction'])
    #    for row in zip(MIList,APC):
    #        writerMI.writerow(row)

    publishMIdata.close()

    print('Mutual Information for %i Site-pairs Calculated in %.2f sec' % (len(pairdata), MItime))

    plt.close('all')
    MI_flat_list = list(a for a in APC)
    hist_export(MI_flat_list, filename)

    print('Success!')

scriptspecs = [['Job', filename], ['files', files],
               ['thresh', str(thresh)], ['damp', str(damp)],
               ['maxSeqCount', str(maxSeqCount)], ['Background Cap', str(bead_ratio)]]

###
PublishData(FREQUENCYseqs, PERCENTseqs, LEADSseqs, CLUSTERseqs, LEADdistances, LEADdistancesBL, scriptspecs,
            looplength, filename)
# UNIQUELOOPseqs rather than CLUSTERseqs will print the full list of unique sequences with the removal of background.
print('\n')
print('ScaffoldSeq completed the requested analyses in %.0f sec' % (time.time() - time1))
print('\n')
print('Results for %s have been published to the output files.' % filename)
print('Press any key to exit...')
getch()
