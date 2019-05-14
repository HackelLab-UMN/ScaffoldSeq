from math import log
import sys

from blosum64_matrix import *
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
