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
from math import log
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
from platform import system
import sys
import textwrap
import time

# TODO: clean up documentation
# TODO: clean up docstrings
# TODO: separate out GUI from algorithm
# TODO: rewrite algorithm in the form of a class?
# TODO: profile code for performance evaluation
# TODO: if this is to be scalable, will need to use generators rather than returning lists
# TODO: test everything with UnitTest
# TODO: replacing all xrange with range. Might need to use generators to yield later
# TODO: update any open() statements into with open() as f statements


if system() != 'Windows':
    #  termios is Unix specific. It does not exist for Windows platforms
    import sys
    import tty
    import termios

'''
Define Global Variables
Purpose: Create a few standard variables to use throughout analysis
'''
# global aa, system_name, adaptor_tolerance, framework_match_threshold
aa = "ACDEFGHIKLMNPQRSTVWY_-"
system_name = system()

# Change this number for the number of mutations you want the adaptor alignment to allow
# Note: This is typically 0, but we accept there is some experimental error with Illumina Sequencing which we estimate
# at 0.4% false negatives
adaptor_tolerance = 0

# Change this number for the fraction of framework homology that is required to identify a protein from a sequence
# Note: We recommend 0.7, but this can range based on the extent of diversification introduced outside the intentially
# diversified regions
framework_match_threshold = 0.7


###################################################################
##                                                               ##
##              List of Functions Used within Script             ##
##                                                               ##
###################################################################

# clear = lambda: os.system('cls')


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


def residue_frequency(LOOPseqs, loop):
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
                    LOOPcnt[m][n] += (memCnt[m][n]) ** (damp)
            memCnt = []
    return LOOPcnt, LOOP_leads


'''
Site-wise Amino Acid Counts within unclustered (pairwise) data
Purpose: Collects residue frequencies for single list of unique sequences
Notes:
len(FULLcnt) = length of diversified region
loop = (min region length, max region length)
'''


def unique_frequency(FULLseqs):
    FULLcnt = []
    for a in range(len(FULLseqs[0][0])):
        FULLcnt.append([0.] * 22)
    for mem in FULLseqs:
        for j in range(len(FULLseqs[0][0])):
            FULLcnt[j][aa.index(mem[0][j])] += (mem[1]) ** damp
    return FULLcnt


'''
Sequence Analysis (Pairwise frequencies)
Purpose: Counts the occurrences/frequencies of every pairwise interaction 
    between two sites in a diversified region
Background removal and dampening take place; however, clustering is less 
    relevant for this epistatic analysis.
Output Structure: 22x22 matrix for each position pair, i and j
    Pair sequence (i,j) ordering: i = {0:(len(protein) - 1)}, j = {(i+1):len(protein)}
Note: 22 = 20AA + stop + gap
'''


def full_pairwise(seqs, pairBG, damp):
    allpair = []

    FREQUENCYpairseqs = unique_frequency(seqs)

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


'''
Sequence Analysis (Output of pairwise analysis)
Purpose: Creates lines for later output in an excel form
Note: Also compares experimental pairwise frequencies to what would be expected probabilistically
      what would happen if each residue operated fully independently without interaction with nearby amino acids
Inputs:         'pairdata'
                'freqmatrix'

'''


def FullMatrixCompare(pairdata, freqmatrix):
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


###################################################################
###     The below section is devoted to the user interface      ###
###################################################################


def clear():
    """
    Interface Modification (Clear characters)
    Purpose: Clears all text currently in interface
    :return:
    """
    os.system(['clear', 'cls'][os.name == 'nt'])


def translate_dna(sequence):
    """
    Sequence Analysis (DNA to AAs)
    Purpose: Translates any input DNA sequence into amino acids
    Note: Notably, used in interface to provide translated sequence while entering DNA information real-time
    :param sequence:
    :return:
    """
    codontable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
        '---': '-', 'NNN': 'X', '***': '-', '&&&': '+'
    }
    proteinsequence = ''
    for n in range(0, len(sequence), 3):
        if sequence[n:(n + 3)] in codontable:
            proteinsequence += codontable[sequence[n:n + 3]]
        else:
            proteinsequence += '!'
    if not proteinsequence: return ' ? '
    return proteinsequence


def LoadJobs():
    """
    Interface Menu (Load Job)
    Purpose: Attempts to access pickled 'SavedJobs.p' folder for previously worked files, with corresponding menu
    Note: If no file is detected, it will generate some default protein scaffold settings
    :return:
    """
    try:
        saved = open("SavedJobs.p", "rb")
        DefaultJobs = pickle.load(saved)
        saved.close()
    except:
        DefaultJobs = [['Affibody_ABY025', '', 'GCCGAAGCGAAATAC', 'GCTAGC', 'GGATCC', 2,
                        [['CTGCCGAACCTGACC', 13, 13, 2],
                         ['GACCCGTCCCAGAGCTCTGAACTCCTGTCTGAGGCGAAGAAACTGAACGATTCCCAAGCACCAAAA', 13, 13, 2],
                         ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2]]]]
        DefaultJobs.append(['DARPin', '', 'GACGTTAACGCT', 'GGATCC', 'AAGCTT', 2, [['ACTCCGCTGCACCTGGCTGCT', 6, 6, 0], [
            'GGTCACCTGGAAATCGTTGAAGTTCTGCTGAAGTACGGTGCT', 2, 2, 0], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2],
                                                                                  ['', 4, 6, 2], ['', 4, 6, 2],
                                                                                  ['', 4, 6, 2]]])
        DefaultJobs.append(['Fibronectin_Fn3HP', 'High_affinity.fasta',
                            'TCCTCCGACTCTCCGCGTAACCTGGAGGTTACCAACGCAACTCCGAACTCTCTGACTATTTCTTGG', 'GCTAGC', 'GGATCC', 3,
                            [['TACCGTATCACCTACGGCGAAACTGGTGGTAACTCCCCGAGCCAGGAATTCACTGTTCCG', 6, 10, 3],
                             ['GCGACCATCAGCGGTCTGAAACCGGGCCAGGATTATACCATTACCGTGTACGCTGTA', 3, 7, 1],
                             ['CCAATCAGCATCAATTATCGCACCGAAATCGACAAACCGTCTCAG', 6, 12, 3]] + [['', 4, 6, 2] for j in
                                                                                             range(5)]])
        DefaultJobs.append(
            ['Gene-2-Protein_Gp2', 'Gp2_evolved_binders.fasta', 'AAATTTTGGGCGACTGTA', 'GCTAGC', 'GGATCC', 2,
             [['TTCGAGGTTCCGGTTTATGCTGAAACCCTGGACGAAGCACTGGAACTGGCCGAATGGCAGTAC', 6, 8, 6],
              ['GTGACCCGCGTGCGTCCG', 6, 8, 6], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2],
              ['', 4, 6, 2], ['', 4, 6, 2]]])
        DefaultJobs.append(
            ['Knottin', '', 'GGCCAGTCTGGCCAGGGCACCTGCAACACCCCGGGCTGCACCTGCAGCTGGCCGGTGTGC', 'TGACTAGCAATGCTGACTGA',
             'TCTGGTGACTACAACAAAAAC', 1,
             [['TGCGGCGAAACCTGCGTGGGCGGAGGGCAGTCTGGGCAG', 7, 7, 0], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2],
              ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2]]])

        saved = open("SavedJobs.p", "wb")
        pickle.dump(DefaultJobs, saved)
        saved.close()
    print('\n' * 40 + 'Loaded Jobs:')
    for job in DefaultJobs:
        print(' -', job[0])
    print('\n' * (20 - len(DefaultJobs)))
    __getch()

    JobMenu = ['Job Name: ', 'FASTA/FASTQ File: ', 'Gene Start:', "5' Anchor: ", "3' Anchor: ",
               '# of Diversified Regions: ']
    OptionsMenu = ['Select', 'Delete', 'Return to Main Menu', 'Exit']
    ss, cc, ii = 0, 79, 0
    while True:
        clear()
        print(' ' * 70 + 'Set', str(ii + 1), 'of', str(len(DefaultJobs)))
        print('\n' * 10 + 'Saved Files:'.center(cc))
        for i in range(len(JobMenu)):
            print((JobMenu[i] + '       ' + str(DefaultJobs[ii][i]).center(30) + '   ').center(cc))
        print('')
        for i in range(len(OptionsMenu)):
            if i == ss:
                print(('>  ' + OptionsMenu[i].center(10) + '  <').center(cc))
            else:
                print(('   ' + OptionsMenu[i].center(10) + '   ').center(cc))
        print('\n' * 10)
        key = __getch()
        if key == 'Up':
            ss = (ss - 1) % len(OptionsMenu)
        elif key == 'Down':
            ss = (ss + 1) % len(OptionsMenu)
        elif key == 'Left':
            ii = (ii - 1) % len(DefaultJobs)
        elif key == 'Right':
            ii = (ii + 1) % len(DefaultJobs)
        elif key == 'Enter':
            if ss == 0:
                return DefaultJobs[ii][:-1], DefaultJobs[ii][-1]
            elif ss == 1:
                print('Are you sure?'.center(cc) + '\n' + 'Press escape to abort...'.center(cc))
                __getch()
                del DefaultJobs[ii]
                saved = open("SavedJobs.p", "wb")
                pickle.dump(DefaultJobs, saved)
                saved.close()
                print('Deleted job template.')
                time.sleep(1.0)
                ii = 0
            elif ss == 3:
                raise SystemExit
            elif ss == 2:
                return False, False


def __getch():
    """
    Interface Menu (Keypress Recording)
    Purpose: Logs keypress without using input functions for quicker, smoother experience
    Note: Interprets arrowkeys as well! Should be able to handle Windows, Mac, and Linux system
    :return:
    """
    if system_name == 'Windows':
        from msvcrt import getch
        while True:
            key = ord(getch())
            if key == 27:  # ESC
                raise SystemExit
            elif key == 13:  # Enter
                return 'Enter'
            elif key == 8:
                return 'Delete'
            elif key == 224:  # Special keys (arrows, f keys, ins, del, etc.)
                key = ord(getch())
                if key == 80:  # Down arrow
                    return 'Down'
                elif key == 72:  # Up arrow
                    return 'Up'
                elif key == 77:  # Right arrow
                    return 'Right'
                elif key == 75:  # Left arrow
                    return 'Left'
            else:
                return chr(key)
    else:
        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setraw(sys.stdin.fileno())
            ch = sys.stdin.read(1)
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
        if ch == '\r':
            return 'Enter'
        elif ch == '`':
            exit()
        elif ord(ch) == 127:
            return 'Delete'
        elif ord(ch) == 27:
            ch = sys.stdin.read(2)
            if ch == '[A':
                return 'Up'
            elif ch == '[B':
                return 'Down'
            elif ch == '[C':
                return 'Right'
            elif ch == '[D':
                return 'Left'
        else:
            return ch


def main(settings):
    """
    Interface Menu (Main Menu)
    Purpose: First menu to operate the algorithm, branches to Load Jobs, Settings, and Start Job
    :param settings:
    :return:
    """
    main_options = ['Start Job', 'Load  Job', 'Settings', 'Information', 'Exit']
    clear()
    cc = 79
    print('\n' * 4)
    print('##################################################################'.center(cc))
    print('##                                                              ##'.center(cc))
    print('##                          ScaffoldSeq                         ##'.center(cc))
    print('##                                                              ##'.center(cc))
    print('##                          Hackel Lab                          ##'.center(cc))
    print('##                    University of Minnesota                   ##'.center(cc))
    print('##                                                              ##'.center(cc))
    print('##                          Version 2.0                         ##'.center(cc))
    print('##                    Last Updated: 2/16/2016                   ##'.center(cc))
    print('##                   Woldring - Holec - Hackel                  ##'.center(cc))
    print('##                                                              ##'.center(cc))
    print('##                  Press any key to continue...                ##'.center(cc))
    print('##                                                              ##'.center(cc))
    print('##################################################################'.center(cc))

    __getch()
    clear()

    cc, ss = 79, 0
    while True:
        clear()
        print('\n' * 5 + '-- Main Menu --'.center(cc) + '\n')
        for i in range(len(main_options)):
            if i == ss:
                print('>  ' + main_options[i] + '  <'.center(cc))
            else:
                print(main_options[i].center(cc))
        print('\n' * 9)
        key = __getch()
        if key == 'Up':
            ss = (ss - 1) % len(main_options)
        elif key == 'Down':
            ss = (ss + 1) % len(main_options)
        elif key == 'Enter':
            if ss == 0:
                while True:
                    job, divs = ModifyJob()
                    if job:
                        break
                if divs:
                    return job, divs, settings
            if ss == 1:
                job, divs = LoadJobs()
                if job:
                    while True:
                        job, divs = ModifyJob(job, divs)
                        if job:
                            break
                    if divs:
                        return job, divs, settings
            elif ss == 2:
                settings = MainSettings(settings)
            elif ss == 3:
                clear()
                for i in textwrap.wrap('ScaffoldSeq was designed for quickly analyzing diversified regions within a ' +
                                       'protein or scaffold using high-throughput sequence data. The algorithm uniquely ' +
                                       'clusters similar clones, dampens dominant sequence motifs, and eliminates ' +
                                       'background signals to accurately highlight both heterogeneities and trends ' +
                                       'within very large datasets. The software produces sitewise and pairwise analyses ' +
                                       'for the user-defined regions of interest.'): print(i)
                key = __getch()
            elif ss == 4:
                raise SystemExit


def ModifyJob(job='', divs=''):
    """
    Interface Menu (Start Job)
    Purpose: Allow a user to modify job, whether it is generated from scratch or an existing template

    Variable Structure Notes:
    job = ["" for x in xrange(len(JobMenu)-1)]+[2]    --> ['', '', '', '', '', 2]
    divs = [['',4,6,2] for j in xrange(8)]
    :param job:
    :param divs:
    :return:
    """
    cc, ss, rs = 79, 0, 0
    JobMenu = ['Job Name', 'FASTA/FASTQ File', 'Gene Start', "5' Anchor", "3' Anchor", '# of Diversified Regions']
    RegionMenu = ['DNA After Region', 'Minimum Region Length', 'Maximum Region Length', 'Insert after # Position']
    OptionsMenu = ['Accept', 'Save', 'Cancel']
    if not job:
        job = ["" for x in range(len(JobMenu) - 1)] + [2]
        divs = [['', 4, 6, 2] for j in range(8)]
    limit = len(JobMenu) + len(RegionMenu) + len(OptionsMenu) + 1
    while True:
        clear()
        construct = translate_dna(job[2])
        for i in range(job[5]):
            construct += translate_dna(
                divs[i][3] * '***' + (divs[i][2] - divs[i][1]) * '&&&' + (divs[i][1] - divs[i][3]) * '***')
            construct += translate_dna(divs[i][0])
        print('-- Job Settings --'.center(cc))

        for i in range(len(JobMenu)):
            if i == ss:
                print(JobMenu[i] + ' >' + str(job[i]).center(30) + '< '.center(cc))
            else:
                print(JobMenu[i] + '  ' + str(job[i]).center(30) + '  '.center(cc))
        print('')
        if ss == len(JobMenu):
            print((' >' + '<Region ' + str(rs + 1) + '>' + '< ').center(20).center(cc))
        else:
            print(('  ' + '<Region ' + str(rs + 1) + '>' + '  ').center(20).center(cc))

        for i in range(len(RegionMenu)):
            if i + len(JobMenu) + 1 == ss:
                print(RegionMenu[i] + ' >' + str(divs[rs][i]).center(30) + '< '.center(cc))
            else:
                print(RegionMenu[i] + '  ' + str(divs[rs][i]).center(30) + '  '.center(cc))
        print('')
        for i in range(len(OptionsMenu)):
            if i + len(JobMenu) + len(RegionMenu) + 1 == ss:
                print('>  ' + OptionsMenu[i].center(30) + '  <'.center(cc))
            else:
                print('   ' + OptionsMenu[i].center(30) + '   '.center(cc))
        print('\n' + 'Translated Gene of Interest:'.center(cc))
        print(construct.center(cc) + '\n')
        print('? = undeclared  - = diversified  + = loop length  ! = translation error'.center(cc))

        key = __getch()
        if key == 'Up':
            ss = (ss - 1) % limit
        elif key == 'Down':
            ss = (ss + 1) % limit
        elif key == 'Enter':
            if ss == limit - 3:
                try:
                    if construct.count('!') + construct.count('?') != 0: raise SystemError
                    if job[3].count('A') + job[3].count('T') + job[3].count('C') + job[3].count('G') != len(job[3]):
                        raise SystemError
                    if job[4].count('A') + job[4].count('T') + job[4].count('C') + job[4].count('G') != len(job[4]):
                        raise SystemError
                    print('Press any key to start job...'.center(cc) + '\n' + 'Esc to abort'.center(cc))
                    __getch()
                    return job, divs
                except:
                    print('\n' * 30)
                    print('Invalid entry in settings'.center(cc) + '\n' + 'Be sure anchors are properly defined'.center(
                        cc) + 'There are no undeclared or translation errors'.center(cc))
                    print('\n' * 19)
                    __getch()
            elif ss == limit - 2:
                try:
                    saved = open("SavedJobs.p", "rb")
                    DefaultJobs = pickle.load(saved)
                    saved.close()
                except:
                    DefaultJobs = [['Affibody_ABY025', '', 'GCCGAAGCGAAATAC', 'GCTAGC', 'GGATCC', 2,
                                    [['CTGCCGAACCTGACC', 13, 13, 2],
                                     ['GACCCGTCCCAGAGCTCTGAACTCCTGTCTGAGGCGAAGAAACTGAACGATTCCCAAGCACCAAAA', 13, 13, 2],
                                     ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2],
                                     ['', 4, 6, 2]]]]
                    DefaultJobs.append(['DARPin', '', 'GACGTTAACGCT', 'GGATCC', 'AAGCTT', 2,
                                        [['ACTCCGCTGCACCTGGCTGCT', 6, 6, 0],
                                         ['GGTCACCTGGAAATCGTTGAAGTTCTGCTGAAGTACGGTGCT', 2, 2, 0], ['', 4, 6, 2],
                                         ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2]]])
                    DefaultJobs.append(['Fibronectin_Fn3HP', 'High_affinity.fasta',
                                        'TCCTCCGACTCTCCGCGTAACCTGGAGGTTACCAACGCAACTCCGAACTCTCTGACTATTTCTTGG', 'GCTAGC',
                                        'GGATCC', 3,
                                        [['TACCGTATCACCTACGGCGAAACTGGTGGTAACTCCCCGAGCCAGGAATTCACTGTTCCG', 6, 10, 3],
                                         ['GCGACCATCAGCGGTCTGAAACCGGGCCAGGATTATACCATTACCGTGTACGCTGTA', 3, 7, 1],
                                         ['CCAATCAGCATCAATTATCGCACCGAAATCGACAAACCGTCTCAG', 6, 12, 3]] + [['', 4, 6, 2]
                                                                                                         for j in
                                                                                                         range(5)]])
                    DefaultJobs.append(
                        ['Gene-2-Protein_Gp2', 'Gp2_evolved_binders.fasta', 'AAATTTTGGGCGACTGTA', 'GCTAGC', 'GGATCC', 2,
                         [['TTCGAGGTTCCGGTTTATGCTGAAACCCTGGACGAAGCACTGGAACTGGCCGAATGGCAGTAC', 6, 8, 6],
                          ['GTGACCCGCGTGCGTCCG', 6, 8, 6], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2],
                          ['', 4, 6, 2], ['', 4, 6, 2]]])
                    DefaultJobs.append(['Knottin', '', 'GGCCAGTCTGGCCAGGGCACCTGCAACACCCCGGGCTGCACCTGCAGCTGGCCGGTGTGC',
                                        'TGACTAGCAATGCTGACTGA', 'TCTGGTGACTACAACAAAAAC', 1,
                                        [['TGCGGCGAAACCTGCGTGGGCGGAGGGCAGTCTGGGCAG', 7, 7, 0], ['', 4, 6, 2],
                                         ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2],
                                         ['', 4, 6, 2]]])
                    saved = open("SavedJobs.p", "wb")
                    pickle.dump(DefaultJobs, saved)
                    saved.close()
                DefaultJobs.append(job + [divs])
                saved = open("SavedJobs.p", "wb")
                pickle.dump(DefaultJobs, saved)
                saved.close()
                print('\n' * 10 + 'Saved!'.center(cc) + '\n' * 18)
                __getch()
            elif ss == limit - 1:
                return True, False
        elif ss == 5:
            clear()
            if key == 'Left':
                job[5] = (job[5] - 2) % 8 + 1
            if key == 'Right':
                job[5] = (job[5]) % 8 + 1
        elif ss == 6:
            if key == 'Left':
                rs = (rs - 1) % job[5]
            if key == 'Right':
                rs = (rs + 1) % job[5]
        elif ss == 8:
            if key == 'Left' and divs[rs][1] > 1: divs[rs][1] -= 1
            if key == 'Right' and divs[rs][1] < divs[rs][2]: divs[rs][1] += 1
        elif ss == 9:
            if key == 'Left' and divs[rs][2] > divs[rs][1]: divs[rs][2] -= 1
            if key == 'Right': divs[rs][2] += 1
        elif ss == 10:
            if key == 'Left' and divs[rs][3] > 0: divs[rs][3] -= 1
            if key == 'Right' and divs[rs][3] < divs[rs][1]: divs[rs][3] += 1
        elif key == 'Left' or key == 'Right' or key == 'Enter':
            pass
        elif ss <= len(JobMenu):
            job[ss] = job[ss][:-1] if key == "Delete" else job[ss] + key
            # if key == 'Delete': job[ss] = job[ss][:-1]
            # else: job[ss] += key
        elif ss == 7:
            if key == 'Delete':
                divs[rs][ss - len(JobMenu) - 1] = divs[rs][ss - len(JobMenu) - 1][:-1]
            else:
                divs[rs][ss - len(JobMenu) - 1] += key


def MainSettings(settings):
    """
    Interface Menu (Settings Menu)
    Purpose: Allows user to modify the behind the scenes settings for the analysis
    :param settings:
    :return:
    """
    settings_menu = ["Sequence Similarity Threshold", "Frequency Dampening Power", "Maximum Sequence Count",
                     "Assay Background Filter", "Pairwise Analysis", ["Filter Coefficient", "Hard Cap Filter"], ""]
    cc, ss, limit, switch, pairwise_analysis = 79, 0, len(settings_menu), ['On', 'Off'], ['On', 'Off']
    while True:
        clear()
        print('\n' * 2)
        print('-- System Settings --'.center(cc) + '\n')
        for i in range(len(settings_menu)):
            if i == 5:
                display = settings_menu[i][switch.index(settings[3])]
            else:
                display = settings_menu[i]
            if i == ss:
                print(display + '    >  ' + settings[i] + '  <'.center(cc))
            else:
                print(display + '       ' + settings[i] + '   '.center(cc))
        print('\n' * 10)
        key = __getch()
        if key == 'Up':
            ss = (ss - 1) % len(settings_menu)
        elif key == 'Down':
            ss = (ss + 1) % len(settings_menu)
        elif key == 'Left':
            if ss == 3:
                settings[3] = switch[1 - switch.index(settings[3])]
            if ss == 4:
                settings[4] = pairwise_analysis[1 - pairwise_analysis.index(settings[4])]
        elif key == 'Right':
            if ss == 3:
                settings[3] = switch[1 - switch.index(settings[3])]
            if ss == 4:
                settings[4] = pairwise_analysis[1 - pairwise_analysis.index(settings[4])]
        elif ss != limit - 3 and ss != limit - 1:
            if key == 'Delete':
                settings[ss] = settings[ss][:-1]
            else:
                settings[ss] += key
        elif key == 'Enter' and ss == len(settings_menu) - 1:
            try:
                [float(i) for i in settings[0:3] + [settings[5]]]
                return settings
            except:
                print('\n' * 30)
                print('Invalid entry in settings'.center(cc) + '\n' + 'Be sure all values are numbers'.center(cc))
                print('\n' * 19)
                __getch()


def Silent_Load(file_name):
    text = open(file_name, 'r')
    lines = text.readlines()
    JobMenu = ['Job Name', 'FASTA/FASTQ File', 'Gene Start', "5' Anchor", "3' Anchor", '# of Diversified Regions']
    RegionMenu = ['DNA After Region', 'Minimum Region Length', 'Maximum Region Length', 'Insert after # Position']
    SettingMenu = ["Sequence Similarity Threshold", "Frequency Dampening Power", "Maximum Sequence Count",
                   "Assay Background Filter", "Pairwise Analysis", "Filter Coefficient", "Hard Cap Filter"]
    jobs, regions, settings = ['' for _ in range(len(JobMenu))], \
                              ['' for _ in range(len(RegionMenu))], \
                              ['' for _ in range(len(SettingMenu))]
    for i in range(len(JobMenu)):
        for line in lines:
            if JobMenu[i] in line:
                jobs[i] = line[line.index(JobMenu[i]) + len(JobMenu[i]):].replace(' ', '').replace('\n', '').replace(
                    ':', '')
    for i in range(len(RegionMenu)):
        for line in lines:
            if RegionMenu[i] in line:
                regions[i] = line[line.index(RegionMenu[i]) + len(RegionMenu[i]):].replace(' ', '').replace('\n',
                                                                                                            '').replace(
                    ':', '')
    for i in range(len(SettingMenu)):
        for line in lines:
            if SettingMenu[i] in line:
                settings[i] = line[line.index(SettingMenu[i]) + len(SettingMenu[i]):].replace(' ', '').replace('\n',
                                                                                                               '').replace(
                    ':', '')
    temp1 = jobs[:-1] + [int(jobs[-1])]
    temp2 = [[a, int(b), int(c), int(d)] for a, b, c, d in
             zip(regions[0].split(','), regions[1].split(','), regions[2].split(','), regions[3].split(','))]
    temp3 = [float(settings[0]), float(settings[1]), int(settings[2])] + settings[3:5] + [int(settings[5])]
    return temp3, temp1, temp2


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
        __getch()
elif len(sys.argv) > 2:
    print('Too many input arguments. Please revise command line arguments.')
    print('Press any key to exit.')
    input()
    raise SystemExit

start = time.time()  # Starts timer
curr_date = datetime.datetime.now().strftime("%y-%m-%d-%H-%M")

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

###################################################################
###     The above section is devoted to the user interface      ###
###################################################################


###################################################################
###       The below section is devoted to the algorithm         ###
###################################################################

# Generate empty count matrices for later counting
LOOPcnt = [[[[0. * 22] for j in range(i)] for i in range(pos[0], pos[1] + 1)] for pos in looplength]


###      Set up variables to be used for sequence sorting


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

                #global framework_match_threshold
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

                ###         Take into account loop length diversity (-1...+1)
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

###################################################################
##                                                               ##
##                   Operate on AA sequences                     ##
##      List, Tally  Occurrences and Sort Unique Sequences       ##
##                                                               ##
###################################################################


###         Enumerate and Remove Duplicates by Loop Length

UNIQUELOOPseqs = [[[[str(x), LOOPseq.count(x)] for x in list(set(LOOPseq))] for LOOPseq in LOOPset] for LOOPset in
                  LOOPseqs]
UNIQUELOOPseqs = [[sorted(LOOPseq, key=lambda x: -x[1]) for LOOPseq in LOOPset] for LOOPset in UNIQUELOOPseqs]

###################################
##           Time Point          ##
###################################

scantime1 = time.time() - time1b
print('Organized Diversified Regions in %.1f sec' % (scantime1))
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

###################################################################
##           Quantify Normalized Amino Acid Frequency            ##
##      Output Includes Positional AA Frequency and Occurrence   ##
###################################################################

# Counts occurrences of amino acids at each position
FREQUENCYseqs = [residue_frequency(LOOPseq, loop)[0] for LOOPseq, loop in zip(CLUSTERseqs, looplength)]
#
LEADSseqs = [residue_frequency(LOOPseq, loop)[1] for LOOPseq, loop in zip(CLUSTERseqs, looplength)]
# Calculates the fractional distribution of amino acids at each position
PERCENTseqs = [[[i / sum(j) for i in j] for j in loop] for loop in FREQUENCYseqs]

# print PERCENTseqs

###

sorttime = time.time() - time3b
print('Site-wise Frequency Matrix Constructed in %.2f sec' % (sorttime))
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
         ('A', '-'): -1
    , ('R', 'A'): -2, ('R', 'R'): 5, ('R', 'N'): 0, ('R', 'D'): -2, ('R', 'C'): -3, ('R', 'Q'): 1, ('R', 'E'): 0,
         ('R', 'G'): -2, ('R', 'H'): 0, ('R', 'I'): -3, ('R', 'L'): -2, ('R', 'K'): 2, ('R', 'M'): -2, ('R', 'F'): -3,
         ('R', 'P'): -2, ('R', 'S'): -1, ('R', 'T'): -1, ('R', 'W'): -2, ('R', 'Y'): -2, ('R', 'V'): -3, ('R', '_'): -9,
         ('R', '-'): -1
    , ('N', 'A'): -1, ('N', 'R'): 0, ('N', 'N'): 6, ('N', 'D'): 1, ('N', 'C'): -3, ('N', 'Q'): 0, ('N', 'E'): 0,
         ('N', 'G'): -1, ('N', 'H'): 1, ('N', 'I'): -3, ('N', 'L'): -3, ('N', 'K'): 0, ('N', 'M'): -2, ('N', 'F'): -3,
         ('N', 'P'): -2, ('N', 'S'): 0, ('N', 'T'): 0, ('N', 'W'): -3, ('N', 'Y'): -2, ('N', 'V'): -3, ('N', '_'): -9,
         ('N', '-'): -1
    , ('D', 'A'): -2, ('D', 'R'): -2, ('D', 'N'): 1, ('D', 'D'): 6, ('D', 'C'): -3, ('D', 'Q'): 0, ('D', 'E'): 2,
         ('D', 'G'): -2, ('D', 'H'): -1, ('D', 'I'): -3, ('D', 'L'): -3, ('D', 'K'): -1, ('D', 'M'): -3, ('D', 'F'): -3,
         ('D', 'P'): -2, ('D', 'S'): 0, ('D', 'T'): -1, ('D', 'W'): -4, ('D', 'Y'): -3, ('D', 'V'): -3, ('D', '_'): -9,
         ('D', '-'): -1
    , ('C', 'A'): -1, ('C', 'R'): -3, ('C', 'N'): -3, ('C', 'D'): -3, ('C', 'C'): 9, ('C', 'Q'): -3, ('C', 'E'): -4,
         ('C', 'G'): -3, ('C', 'H'): -2, ('C', 'I'): -1, ('C', 'L'): -1, ('C', 'K'): -3, ('C', 'M'): -1, ('C', 'F'): -2,
         ('C', 'P'): -3, ('C', 'S'): -1, ('C', 'T'): -1, ('C', 'W'): -3, ('C', 'Y'): -2, ('C', 'V'): -1, ('C', '_'): -9,
         ('C', '-'): -2
    , ('Q', 'A'): -1, ('Q', 'R'): 1, ('Q', 'N'): 0, ('Q', 'D'): 0, ('Q', 'C'): -3, ('Q', 'Q'): 5, ('Q', 'E'): 2,
         ('Q', 'G'): -2, ('Q', 'H'): 1, ('Q', 'I'): -3, ('Q', 'L'): -2, ('Q', 'K'): 1, ('Q', 'M'): 0, ('Q', 'F'): -3,
         ('Q', 'P'): -1, ('Q', 'S'): 0, ('Q', 'T'): 0, ('Q', 'W'): -2, ('Q', 'Y'): -2, ('Q', 'V'): -2, ('Q', '_'): -9,
         ('Q', '-'): -1
    , ('E', 'A'): -1, ('E', 'R'): 0, ('E', 'N'): 0, ('E', 'D'): 2, ('E', 'C'): -4, ('E', 'Q'): 2, ('E', 'E'): 5,
         ('E', 'G'): -2, ('E', 'H'): 0, ('E', 'I'): -3, ('E', 'L'): -3, ('E', 'K'): 1, ('E', 'M'): -2, ('E', 'F'): -3,
         ('E', 'P'): -1, ('E', 'S'): 0, ('E', 'T'): -1, ('E', 'W'): -3, ('E', 'Y'): -2, ('E', 'V'): -2, ('E', '_'): -9,
         ('E', '-'): -1
    , ('G', 'A'): 0, ('G', 'R'): -2, ('G', 'N'): -1, ('G', 'D'): -2, ('G', 'C'): -3, ('G', 'Q'): -2, ('G', 'E'): -2,
         ('G', 'G'): 6, ('G', 'H'): -2, ('G', 'I'): -4, ('G', 'L'): -4, ('G', 'K'): -2, ('G', 'M'): -3, ('G', 'F'): -3,
         ('G', 'P'): -2, ('G', 'S'): -1, ('G', 'T'): -2, ('G', 'W'): -2, ('G', 'Y'): -3, ('G', 'V'): -3, ('G', '_'): -9,
         ('G', '-'): -2
    , ('H', 'A'): -2, ('H', 'R'): 0, ('H', 'N'): 1, ('H', 'D'): -1, ('H', 'C'): -2, ('H', 'Q'): 1, ('H', 'E'): 0,
         ('H', 'G'): -2, ('H', 'H'): 7, ('H', 'I'): -3, ('H', 'L'): -3, ('H', 'K'): -1, ('H', 'M'): -1, ('H', 'F'): -1,
         ('H', 'P'): -2, ('H', 'S'): -1, ('H', 'T'): -2, ('H', 'W'): -1, ('H', 'Y'): 1, ('H', 'V'): -3, ('H', '_'): -9,
         ('H', '-'): -1
    , ('I', 'A'): -1, ('I', 'R'): -3, ('I', 'N'): -3, ('I', 'D'): -3, ('I', 'C'): -1, ('I', 'Q'): -3, ('I', 'E'): -3,
         ('I', 'G'): -4, ('I', 'H'): -3, ('I', 'I'): 4, ('I', 'L'): 2, ('I', 'K'): -3, ('I', 'M'): 1, ('I', 'F'): 0,
         ('I', 'P'): -3, ('I', 'S'): -2, ('I', 'T'): -1, ('I', 'W'): -2, ('I', 'Y'): -1, ('I', 'V'): 2, ('I', '_'): -9,
         ('I', '-'): -1
    , ('L', 'A'): -2, ('L', 'R'): -2, ('L', 'N'): -3, ('L', 'D'): -3, ('L', 'C'): -1, ('L', 'Q'): -2, ('L', 'E'): -3,
         ('L', 'G'): -4, ('L', 'H'): -3, ('L', 'I'): 2, ('L', 'L'): 4, ('L', 'K'): -2, ('L', 'M'): 2, ('L', 'F'): 1,
         ('L', 'P'): -3, ('L', 'S'): -2, ('L', 'T'): -1, ('L', 'W'): -1, ('L', 'Y'): -1, ('L', 'V'): 1, ('L', '_'): -9,
         ('L', '-'): -1
    , ('K', 'A'): -1, ('K', 'R'): 2, ('K', 'N'): 0, ('K', 'D'): -1, ('K', 'C'): -3, ('K', 'Q'): 1, ('K', 'E'): 1,
         ('K', 'G'): -2, ('K', 'H'): -1, ('K', 'I'): -3, ('K', 'L'): -2, ('K', 'K'): 5, ('K', 'M'): -1, ('K', 'F'): -3,
         ('K', 'P'): -1, ('K', 'S'): 0, ('K', 'T'): -1, ('K', 'W'): -3, ('K', 'Y'): -2, ('K', 'V'): -2, ('K', '_'): -9,
         ('K', '-'): -1
    , ('M', 'A'): -1, ('M', 'R'): -2, ('M', 'N'): -2, ('M', 'D'): -3, ('M', 'C'): -1, ('M', 'Q'): 0, ('M', 'E'): -2,
         ('M', 'G'): -3, ('M', 'H'): -1, ('M', 'I'): 1, ('M', 'L'): 2, ('M', 'K'): -1, ('M', 'M'): 6, ('M', 'F'): 0,
         ('M', 'P'): -2, ('M', 'S'): -1, ('M', 'T'): -1, ('M', 'W'): -2, ('M', 'Y'): -1, ('M', 'V'): 0, ('M', '_'): -9,
         ('M', '-'): -1
    , ('F', 'A'): -2, ('F', 'R'): -3, ('F', 'N'): -3, ('F', 'D'): -3, ('F', 'C'): -2, ('F', 'Q'): -3, ('F', 'E'): -3,
         ('F', 'G'): -3, ('F', 'H'): -1, ('F', 'I'): 0, ('F', 'L'): 1, ('F', 'K'): -3, ('F', 'M'): 0, ('F', 'F'): 6,
         ('F', 'P'): -3, ('F', 'S'): -2, ('F', 'T'): -2, ('F', 'W'): 1, ('F', 'Y'): 3, ('F', 'V'): -1, ('F', '_'): -9,
         ('F', '-'): -1
    , ('P', 'A'): -1, ('P', 'R'): -2, ('P', 'N'): -2, ('P', 'D'): -2, ('P', 'C'): -3, ('P', 'Q'): -1, ('P', 'E'): -1,
         ('P', 'G'): -2, ('P', 'H'): -2, ('P', 'I'): -3, ('P', 'L'): -3, ('P', 'K'): -1, ('P', 'M'): -2, ('P', 'F'): -3,
         ('P', 'P'): 7, ('P', 'S'): -1, ('P', 'T'): -1, ('P', 'W'): -3, ('P', 'Y'): -3, ('P', 'V'): -2, ('P', '_'): -9,
         ('P', '-'): -2
    , ('S', 'A'): 1, ('S', 'R'): -1, ('S', 'N'): 0, ('S', 'D'): 0, ('S', 'C'): -1, ('S', 'Q'): 0, ('S', 'E'): 0,
         ('S', 'G'): -1, ('S', 'H'): -1, ('S', 'I'): -2, ('S', 'L'): -2, ('S', 'K'): 0, ('S', 'M'): -1, ('S', 'F'): -2,
         ('S', 'P'): -1, ('S', 'S'): 4, ('S', 'T'): 1, ('S', 'W'): -3, ('S', 'Y'): -2, ('S', 'V'): -2, ('S', '_'): -9,
         ('S', '-'): -1
    , ('T', 'A'): 0, ('T', 'R'): -1, ('T', 'N'): 0, ('T', 'D'): -1, ('T', 'C'): -1, ('T', 'Q'): 0, ('T', 'E'): -1,
         ('T', 'G'): -2, ('T', 'H'): -2, ('T', 'I'): -1, ('T', 'L'): -1, ('T', 'K'): -1, ('T', 'M'): -1, ('T', 'F'): -2,
         ('T', 'P'): -1, ('T', 'S'): 1, ('T', 'T'): 5, ('T', 'W'): -3, ('T', 'Y'): -2, ('T', 'V'): 0, ('T', '_'): -9,
         ('T', '-'): -1
    , ('W', 'A'): -3, ('W', 'R'): -2, ('W', 'N'): -3, ('W', 'D'): -4, ('W', 'C'): -3, ('W', 'Q'): -2, ('W', 'E'): -3,
         ('W', 'G'): -2, ('W', 'H'): -1, ('W', 'I'): -2, ('W', 'L'): -1, ('W', 'K'): -3, ('W', 'M'): -2, ('W', 'F'): 1,
         ('W', 'P'): -3, ('W', 'S'): -3, ('W', 'T'): -3, ('W', 'W'): 11, ('W', 'Y'): 2, ('W', 'V'): -3, ('W', '_'): -9,
         ('W', '-'): -2
    , ('Y', 'A'): -2, ('Y', 'R'): -2, ('Y', 'N'): -2, ('Y', 'D'): -3, ('Y', 'C'): -2, ('Y', 'Q'): -2, ('Y', 'E'): -2,
         ('Y', 'G'): -3, ('Y', 'H'): 1, ('Y', 'I'): -1, ('Y', 'L'): -1, ('Y', 'K'): -2, ('Y', 'M'): -1, ('Y', 'F'): 3,
         ('Y', 'P'): -3, ('Y', 'S'): -2, ('Y', 'T'): -2, ('Y', 'W'): 2, ('Y', 'Y'): 7, ('Y', 'V'): -1, ('Y', '_'): -9,
         ('Y', '-'): -1
    , ('V', 'A'): 0, ('V', 'R'): -3, ('V', 'N'): -3, ('V', 'D'): -3, ('V', 'C'): -1, ('V', 'Q'): -2, ('V', 'E'): -2,
         ('V', 'G'): -3, ('V', 'H'): -3, ('V', 'I'): 2, ('V', 'L'): 1, ('V', 'K'): -2, ('V', 'M'): 0, ('V', 'F'): -1,
         ('V', 'P'): -2, ('V', 'S'): -2, ('V', 'T'): 0, ('V', 'W'): -3, ('V', 'Y'): -1, ('V', 'V'): 4, ('V', '_'): -9,
         ('V', '-'): -1
    , ('_', 'A'): -9, ('_', 'R'): -9, ('_', 'N'): -9, ('_', 'D'): -9, ('_', 'C'): -9, ('_', 'Q'): -9, ('_', 'E'): -9,
         ('_', 'G'): -9, ('_', 'H'): -9, ('_', 'I'): -9, ('_', 'L'): -9, ('_', 'K'): -9, ('_', 'M'): -9, ('_', 'F'): -9,
         ('_', 'P'): -9, ('_', 'S'): -9, ('_', 'T'): -9, ('_', 'W'): -9, ('_', 'Y'): -9, ('_', 'V'): -9, ('_', '_'): -9,
         ('_', '-'): -9
    , ('-', 'A'): -1, ('-', 'R'): -1, ('-', 'N'): -1, ('-', 'D'): -1, ('-', 'C'): -2, ('-', 'Q'): -1, ('-', 'E'): -1,
         ('-', 'G'): -2, ('-', 'H'): -1, ('-', 'I'): -1, ('-', 'L'): -1, ('-', 'K'): -1, ('-', 'M'): -1, ('-', 'F'): -1,
         ('-', 'P'): -2, ('-', 'S'): -1, ('-', 'T'): -1, ('-', 'W'): -2, ('-', 'Y'): -1, ('-', 'V'): -1, ('-', '_'): -9,
         ('-', '-'): 0
         }


#   Showing similarity as large value through the addition on ("len(a) - ...")


def quick_diff(a, b):
    """
    Sequence Analysis (Number of differences)
    Purpose: Quickly calculates character difference count between two strings or lists
    :param a: some string or list
    :param b: some string or list
    :return:
    """
    return len(a) - sum(a[i] != b[i] for i in range(len(a)))


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
__getch()
