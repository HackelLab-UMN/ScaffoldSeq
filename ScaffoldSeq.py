"""
ScaffoldSeq: Software for characterization of directed evolution populations
University of Minnesota
__author__ = "Daniel R. Woldring, Patrick V. Holec, and Bejnamin J. Hackel"
__version__ = 2.0
__status__ = "Research and Development"
__maintainer__ = "Adam T. McConnell"
"""
import csv
import datetime
import difflib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

# TODO: clean up documentation
# TODO: clean up docstrings
# TODO: rewrite algorithm in the form of a class?
# TODO: profile code for performance evaluation
# TODO: if this is to be scalable, will need to use generators rather than returning lists
# TODO: test everything with UnitTest
# TODO: replacing all xrange with range. Might need to use generators to yield later
# TODO: update any open() statements into with open() as f statements

'''
Define Global Variables
Purpose: Create a few standard variables to use throughout analysis
'''
# global aa, system_name, adaptor_tolerance, framework_match_threshold
aa = "ACDEFGHIKLMNPQRSTVWY_-"

# Change this number for the number of mutations you want the adaptor alignment to allow
# Note: This is typically 0, but we accept there is some experimental error with Illumina Sequencing which we estimate
# at 0.4% false negatives
adaptor_tolerance = 0

# Change this number for the fraction of framework homology that is required to identify a protein from a sequence
# Note: We recommend 0.7, but this can range based on the extent of diversification introduced outside the intentially
# diversified regions
framework_match_threshold = 0.7

curr_date = datetime.datetime.now().strftime("%y-%m-%d-%H-%M")

"""
###################################################################
##                                                               ##
##              List of Functions Used within Script             ##
##                                                               ##
###################################################################
"""


def full_format(roughseqs, big_loop, spacer):
    """
    Sequence Formatting (Compensate for loop length diversity). Formats the sequence of interest to append dashes in the
     place of the loop length diversity
    Note: The spacer exists to add dashes at a specific index so frequency analysis can be applied to loops of all sizes
    :param roughseqs:
    :param big_loop:
    :param spacer:
    :return:
    """
    allseqs = []
    for seqs in roughseqs:
        seq = ''
        for i in range(len(seqs)):
            seq += seqs[i][:(len(seqs[i]) - spacer[i])] \
                   + '-' * (big_loop[i][1] - len(seqs[i])) \
                   + seqs[i][len(seqs[i]) - spacer[i]:]
        allseqs.append(seq)
    return allseqs


def full_freq(seqs):
    """
    Sequence Analysis (Unique sequence frequency). Counts the frequency of unique sequences in data
    :param seqs:
    :return:
    """
    from collections import Counter
    counter = Counter(seqs)
    return [[i, j] for i, j in zip(counter.keys(), counter.values())]


def residue_frequency(LOOPseqs, loop, damp):
    """
    Site-wise Amino Acid Counts within clustered data. Collects residue frequencies for the clustered data
    len(LOOPcnt) = length of diversified region
    :param LOOPseqs:
    :param loop: (min region length, max region length) tuple
    :return:loop count and loop leads
    """
    memCnt, LOOPcnt, LOOP_fams, LOOP_leads = [], [], 0, []
    for a in range(loop[1]):
        LOOPcnt.append([0.] * 22)
    for i in range(loop[1] - loop[0] + 1):
        for fam in LOOPseqs[i]:
            if len(fam) > 0:
                LOOP_fams += 1
                LOOP_leads.append([fam[0], LOOP_fams])
            for k in range(loop[1]):
                memCnt.append([0.] * 22)
            for mem in fam:
                for j in range(loop[1]):
                    memCnt[j][aa.index(mem[0][j])] += mem[1]
            for m in range(loop[1]):
                for n in range(22):
                    LOOPcnt[m][n] += (memCnt[m][n]) ** damp
            memCnt = []
    return LOOPcnt, LOOP_leads


def unique_frequency(FULLseqs, damp):
    """
    Site-wise Amino Acid Counts within unclustered (pairwise) data
    Purpose: Collects residue frequencies for single list of unique sequences
    Notes:
    len(FULLcnt) = length of diversified region
    loop = (min region length, max region length)
    :param FULLseqs:
    :return:
    """
    FULLcnt = []
    for a in range(len(FULLseqs[0][0])):
        FULLcnt.append([0.] * 22)
    for mem in FULLseqs:
        for j in range(len(FULLseqs[0][0])):
            FULLcnt[j][aa.index(mem[0][j])] += (mem[1]) ** damp
    return FULLcnt


def full_pairwise(seqs, pairBG, damp):
    """
    Sequence Analysis (Pairwise frequencies)
    Purpose: Counts the occurrences/frequencies of every pairwise interaction  between two sites in a diversified region
    Background removal and dampening take place; however, clustering is less  relevant for this epistatic analysis.
    Output Structure: 22x22 matrix for each position pair, i and j
        Pair sequence (i,j) ordering: i = {0:(len(protein) - 1)}, j = {(i+1):len(protein)}
        Note: 22 = 20AA + stop + gap
    :param seqs:
    :param pairBG:
    :param damp:
    :return:
    """
    allpair = []

    FREQUENCYpairseqs = unique_frequency(seqs, damp)

    PERCENTpairseqs = [[i / sum(j) for i in j] for j in FREQUENCYpairseqs]

    for i in range(len(seqs[-1][0]) - 1):
        for j in range(i + 1, len(seqs[-1][0])):
            allpair.append([[0. for ii in range(len(aa))] for jj in range(len(aa))])
            for seq in seqs:
                try:
                    if seq[1] > pairBG: allpair[-1][aa.index(seq[0][i])][aa.index(seq[0][j])] += seq[1] ** damp
                except:
                    print('Error- skipping sequence...')
            totalentries = sum([ii for jj in allpair[-1] for ii in jj])
            allpair[-1] = [[ii / totalentries for ii in jj] for jj in allpair[-1]]
    return allpair, PERCENTpairseqs


def FullMatrixCompare(pairdata, freqmatrix):
    """
    Sequence Analysis (Output of pairwise analysis). Creates lines for later output in an excel form
    Note: Also compares experimental pairwise frequencies to what would be expected probabilistically
      what would happen if each residue operated fully independently without interaction with nearby amino acids
    Inputs:         'pairdata'
                    'freqmatrix'
    :param pairdata:
    :param freqmatrix:
    :return:
    """
    # Flatten the lists of each diversified region into one continuous list
    # freqmatrix = [a for b in freqmatrix for a in b]
    ss = len(freqmatrix)
    freqdata, refs, diff = [], [], []
    for i in range(ss - 1):
        for j in range(i + 1, ss):
            freqdata.append([[a * b for b in freqmatrix[j]] for a in freqmatrix[i]])
            refs.append(
                ["'" + '-' * i + 'X' + '-' * (j - i - 1) + 'X' + '-' * (ss - j - 1), str(i + 1) + " & " + str(j + 1)])
    # Comparison between both data sets (pairwise experimental and probability)
    for p1, p2, ref, i in zip(pairdata, freqdata, refs, range(len(pairdata))):
        diff.append([ref[0], ref[1], sum(
            [abs(a - b) for a, b in zip([ii for jj in p1 for ii in jj], [iii for jjj in p2 for iii in jjj])])])
    return diff, freqdata, refs


def FamCrit(seq1, seq2):  # Two input strings and finding ratio of difference
    """
    Purpose: Computes sitewise similarity between the input strings, seq1 and seq2
    :param seq1:
    :param seq2:
    :return:
    """
    return difflib.SequenceMatcher(a=seq1.lower(), b=seq2.lower()).ratio()


def find(lst, predicate):  #
    return (i for i, j in enumerate(lst) if predicate(j)).next()


def PublishData(FREQUENCYseqs, PERCENTseqs, LEADSseqs, UNIQUEseqs, LEADdistances, LEADdistancesBL, scriptspecs,
                looplength, filename):
    """
    Data Output (Basic CSV Display). Bare minimum data output using only default Python libraries
    :param FREQUENCYseqs:
    :param PERCENTseqs:
    :param LEADSseqs:
    :param UNIQUEseqs:
    :param LEADdistances:
    :param LEADdistancesBL:
    :param scriptspecs:
    :param looplength:
    :param filename:
    :return:
    """
    for FREQUENCY, PERCENT, LEADS, UNIQUE, LEADD, LEADBL, LOOP, i in zip(FREQUENCYseqs, PERCENTseqs, LEADSseqs,
                                                                         UNIQUEseqs, LEADdistances, LEADdistancesBL,
                                                                         looplength, range(len(looplength))):
        publishdata = open(filename + '_Region_' + str(i + 1) + '.csv', 'w')
        writer1 = csv.writer(publishdata)

        # Job Overview
        writer1.writerow(['JOB SUMMARY:'])
        for specs in scriptspecs:
            writer1.writerow(specs)
        writer1.writerow(['Diversification region: ', '#' + str(i + 1)])
        writer1.writerow(['Length diversity range:', LOOP[0], LOOP[1]])
        writer1.writerow('\n')

        # Number of times each residue occurs
        writer1.writerow(['Residue Occurrences:'])
        writer1.writerow(list(aa))
        for row in FREQUENCY:
            writer1.writerow(row)
        writer1.writerow('\n')

        # Frequency of times each residue occurs
        writer1.writerow(['Residue Frequency:'])
        writer1.writerow(list(aa))
        for row in PERCENT:
            writer1.writerow(row)
        writer1.writerow('\n')

        '''
        Take PERCENTseqs (site-wise diveristy matrix) as Heat map figure maker.
        Upper and lower bounds (vmin=, vmax=) for color bar can be adjusted to highlight poorly
        represented residues. 
        '''
        vmin = 0.0
        vmax = 0.6

        plt.close('all')

        freq_norm = PERCENT[:]
        norm_mat = [np.sum(a[:20]) for a in PERCENT]

        for pi in range(len(PERCENT)):
            for pj in range(len(PERCENT[pi]) - 2):
                freq_norm[pi][pj] = PERCENT[pi][pj] / norm_mat[pi]

        data = zip(*np.array(freq_norm))
        df = pd.DataFrame(data)

        Cols = [str(s + 1) for s in range(len(PERCENT))]

        Index = list('ACDEFGHIKLMNPQRSTVWY*-')

        df = pd.DataFrame(data, index=Index, columns=Cols)
        df_alt = df.copy()

        df_max = df_alt.drop(df_alt.index[[20, 21]]).max().max()
        df_min = df_alt.drop(df_alt.index[[20, 21]]).min().min()
        df = df_alt.drop(df_alt.index[[20, 21]])

        fig, ax = plt.subplots(figsize=(15, 10))

        ax.xaxis.tick_top()
        ax.yaxis.tick_left()

        plt.text(0.5, 1.08, 'Sitewise Frequency Analysis - Region ' + str(i + 1), fontsize=25,
                 horizontalalignment='center', transform=ax.transAxes)

        plt.gca().invert_yaxis()
        plt.pcolor(df, cmap='Blues', vmin=vmin, vmax=vmax)

        plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
        plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns)

        plt.ylabel('Residue', fontsize=20)
        plt.xlabel('Site', fontsize=20)

        plt.colorbar()
        plt.show()
        plt.savefig(filename + '_Region-' + str(i + 1) + '_' + curr_date + '.png')

        if i + 1 == len(looplength):
            print('Summarizing figures have been stored within {0} for future reference.'.format(os.getcwd()))

        # Ranking of major protein clusters
        writer1.writerow(['Familial Interconnectivity (Hamming distance scoring):'])
        for row in LEADD:
            writer1.writerow(row)
        writer1.writerow('\n')

        # Ranking of major protein clusters revised with Blosum 64
        writer1.writerow(['Familial Interconnectivity (revised BLOSUM64 scoring):'])
        for row in LEADBL:
            writer1.writerow(row)
        writer1.writerow('\n')

        # Ranking of major protein clusters
        writer1.writerow(['Sequence Leads:'])
        writer1.writerow(['Sequence', 'Frequency', 'Index'])
        for row in LEADS:
            writer1.writerow(row[0] + [row[1]])
        writer1.writerow('\n')

        # All unique sequences listed
        writer1.writerow(['All Unique Sequences (BG removed):'])
        writer1.writerow(['Sequence', 'Frequency'])
        for size in UNIQUE:
            for fam in size[:100]:
                for row in fam:
                    if row: writer1.writerow(row)
        writer1.writerow('\n')

        publishdata.close()


def hist_export(prob, sample_id):
    """
        prob = list of values to be plotted as histogram
        score_type is specific to the method used to calculate plotted values
        sample_id refers to the source of data
    :param prob:
    :param sample_id:
    :return:
    """
    fig, ax = plt.subplots(1)

    x = np.array(prob)
    mu = x.mean()
    median = np.median(x)
    sigma = x.std()

    #   http://matplotlib.org/users/mathtext.html
    # textstr = '$\mu=%.3f$\n$\mathrm{median}=%.3f$\n$\sigma=%.3f$\n$\delta=%.0e$'%(mu, median, sigma, ep_prob)
    textstr = '$\mu=%.3f$\n$\mathrm{median}=%.3f$\n$\sigma=%.3f$' % (mu, median, sigma)

    ax.hist(x, bins=50, align='mid', range=(x.min(), x.max()), color='gray')

    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='cornflowerblue', alpha=0.5)

    # place a text box in upper left in axes coords
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=16,
            verticalalignment='top', bbox=props)
    header = 'Mutual Information (MIp) Histogram'
    ax.set_title(sample_id + ' ' + header)
    plt.savefig(sample_id + '_' + curr_date + '.png', dpi=300)
    return


def Adaptor_Alignment(seq, adaptors):
    """
    Sequence Acquisition (Adaptor Alignment)
    Purpose: This function allows you to check for where the adaptors match the sequence, and output the identified sequence/length
    Note: You can change tolerance by modifying global variable "adaptor_tolerance"
    :param seq:
    :param adaptors:
    :return:
    """
    # global adaptor_tolerance
    if adaptor_tolerance == 0:
        length = len(seq[seq.find(adaptors[0]):seq.find(adaptors[1])]) - len(adaptors[0])
        protein = seq[seq.find(adaptors[0]) + len(adaptors[0]):seq.find(adaptors[1])]
        return length, protein
    else:
        adaptor_index = []
        for i in range(2):
            temp = []
            for j in range(len(seq) - len(adaptors[i]) + 1):
                temp.append(FamCrit(seq[j:j + len(adaptors[i])], adaptors[i]) * (len(adaptors[i])))
                print('Seq:', seq[j:j + len(adaptors[i])], adaptors[i])
                print(temp[-1])
            print(len(adaptors[i]) - max(temp))
            print(adaptors)
            if len(adaptors[i]) - max(temp) <= adaptor_tolerance:
                print(temp.index(max(temp)))
                adaptor_index.append(temp.index(max(temp)))
            else:
                return 0, ''
        length = len(seq[adaptor_index[0]:adaptor_index[1]]) - len(adaptors[0])
        protein = seq[adaptor_index[0] + len(adaptors[0]):adaptor_index[1]]
        return length, protein


# DRW edit 1-25-2016
##def add_gap(loop_list,big_loop,depth):
##    for i in xrange(len(loop_list)):
##        for fam in xrange(len(loop_list[i])):
##            for seq in loop_list[i][fam]:
##                a = 0
##                try:
##                    a += 1
##                    seq[0] = seq[0][:depth + 1] \
##                    + '-' * (big_loop - len(seq[0])) \
##                    + seq[0][depth + 1:]
##                except IndexError:
##                    pass
##    return loop_list

def add_gap(loop_list, big_loop, depth):
    for i in range(len(loop_list)):
        for fam in range(len(loop_list[i])):
            for seq in loop_list[i][fam]:
                try:
                    seq[0] = seq[0][:(len(seq[0]) - depth)] \
                             + '-' * (big_loop - len(seq[0])) \
                             + seq[0][len(seq[0]) - depth:]
                except IndexError:
                    pass
    return loop_list


def quick_diff(a, b):
    """
    Sequence Analysis (Number of differences)
    Purpose: Quickly calculates character difference count between two strings or lists
    :param a: some string or list
    :param b: some string or list
    :return:
    """
    return len(a) - sum(a[i] != b[i] for i in range(len(a)))