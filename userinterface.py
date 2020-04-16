import os
import pickle
from platform import system
import textwrap
import time

system_name = system()

if system() != 'Windows':
    #  termios is Unix specific. It does not exist for Windows platforms
    import sys
    import tty
    import termios


"""
###################################################################
###     The below section is devoted to the user interface      ###
###################################################################
"""


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
                         ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2]]],
                       ['DARPin', '', 'GACGTTAACGCT', 'GGATCC', 'AAGCTT', 2, [['ACTCCGCTGCACCTGGCTGCT', 6, 6, 0], [
                           'GGTCACCTGGAAATCGTTGAAGTTCTGCTGAAGTACGGTGCT', 2, 2, 0], ['', 4, 6, 2], ['', 4, 6, 2],
                                                                              ['', 4, 6, 2],
                                                                              ['', 4, 6, 2], ['', 4, 6, 2],
                                                                              ['', 4, 6, 2]]],
                       ['Fibronectin_Fn3HP', 'High_affinity.fasta',
                        'TCCTCCGACTCTCCGCGTAACCTGGAGGTTACCAACGCAACTCCGAACTCTCTGACTATTTCTTGG', 'GCTAGC', 'GGATCC', 3,
                        [['TACCGTATCACCTACGGCGAAACTGGTGGTAACTCCCCGAGCCAGGAATTCACTGTTCCG', 6, 10, 3],
                         ['GCGACCATCAGCGGTCTGAAACCGGGCCAGGATTATACCATTACCGTGTACGCTGTA', 3, 7, 1],
                         ['CCAATCAGCATCAATTATCGCACCGAAATCGACAAACCGTCTCAG', 6, 12, 3]] + [['', 4, 6, 2] for j in
                                                                                         range(5)]],
                       ['Gene-2-Protein_Gp2', 'Gp2_evolved_binders.fasta', 'AAATTTTGGGCGACTGTA', 'GCTAGC', 'GGATCC', 2,
                        [['TTCGAGGTTCCGGTTTATGCTGAAACCCTGGACGAAGCACTGGAACTGGCCGAATGGCAGTAC', 6, 8, 6],
                         ['GTGACCCGCGTGCGTCCG', 6, 8, 6], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2],
                         ['', 4, 6, 2], ['', 4, 6, 2]]],
                       ['Knottin', '', 'GGCCAGTCTGGCCAGGGCACCTGCAACACCCCGGGCTGCACCTGCAGCTGGCCGGTGTGC',
                        'TGACTAGCAATGCTGACTGA',
                        'TCTGGTGACTACAACAAAAAC', 1,
                        [['TGCGGCGAAACCTGCGTGGGCGGAGGGCAGTCTGGGCAG', 7, 7, 0], ['', 4, 6, 2], ['', 4, 6, 2],
                         ['', 4, 6, 2],
                         ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2]]]]

        saved = open("SavedJobs.p", "wb")
        pickle.dump(DefaultJobs, saved)
        saved.close()
    print('\n' * 40 + 'Loaded Jobs:')
    for job in DefaultJobs:
        print(' -', job[0])
    print('\n' * (20 - len(DefaultJobs)))
    getch()

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
        key = getch()
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
                getch()
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


def getch():
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

    getch()
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
        key = getch()
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
                key = getch()
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

        key = getch()
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
                    getch()
                    return job, divs
                except:
                    print('\n' * 30)
                    print('Invalid entry in settings'.center(cc) + '\n' + 'Be sure anchors are properly defined'.center(
                        cc) + 'There are no undeclared or translation errors'.center(cc))
                    print('\n' * 19)
                    getch()
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
                                     ['', 4, 6, 2]]], ['DARPin', '', 'GACGTTAACGCT', 'GGATCC', 'AAGCTT', 2,
                                                       [['ACTCCGCTGCACCTGGCTGCT', 6, 6, 0],
                                                        ['GGTCACCTGGAAATCGTTGAAGTTCTGCTGAAGTACGGTGCT', 2, 2, 0],
                                                        ['', 4, 6, 2],
                                                        ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2],
                                                        ['', 4, 6, 2]]], ['Fibronectin_Fn3HP', 'High_affinity.fasta',
                                                                          'TCCTCCGACTCTCCGCGTAACCTGGAGGTTACCAACGCAACTCCGAACTCTCTGACTATTTCTTGG',
                                                                          'GCTAGC',
                                                                          'GGATCC', 3,
                                                                          [[
                                                                               'TACCGTATCACCTACGGCGAAACTGGTGGTAACTCCCCGAGCCAGGAATTCACTGTTCCG',
                                                                               6, 10, 3],
                                                                           [
                                                                               'GCGACCATCAGCGGTCTGAAACCGGGCCAGGATTATACCATTACCGTGTACGCTGTA',
                                                                               3, 7, 1],
                                                                           [
                                                                               'CCAATCAGCATCAATTATCGCACCGAAATCGACAAACCGTCTCAG',
                                                                               6, 12, 3]] + [['', 4, 6, 2]
                                                                                             for j in
                                                                                             range(5)]],
                                   ['Gene-2-Protein_Gp2', 'Gp2_evolved_binders.fasta', 'AAATTTTGGGCGACTGTA', 'GCTAGC',
                                    'GGATCC', 2,
                                    [['TTCGAGGTTCCGGTTTATGCTGAAACCCTGGACGAAGCACTGGAACTGGCCGAATGGCAGTAC', 6, 8, 6],
                                     ['GTGACCCGCGTGCGTCCG', 6, 8, 6], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2],
                                     ['', 4, 6, 2],
                                     ['', 4, 6, 2], ['', 4, 6, 2]]],
                                   ['Knottin', '', 'GGCCAGTCTGGCCAGGGCACCTGCAACACCCCGGGCTGCACCTGCAGCTGGCCGGTGTGC',
                                    'TGACTAGCAATGCTGACTGA', 'TCTGGTGACTACAACAAAAAC', 1,
                                    [['TGCGGCGAAACCTGCGTGGGCGGAGGGCAGTCTGGGCAG', 7, 7, 0], ['', 4, 6, 2],
                                     ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2], ['', 4, 6, 2],
                                     ['', 4, 6, 2]]]]
                    saved = open("SavedJobs.p", "wb")
                    pickle.dump(DefaultJobs, saved)
                    saved.close()
                DefaultJobs.append(job + [divs])
                saved = open("SavedJobs.p", "wb")
                pickle.dump(DefaultJobs, saved)
                saved.close()
                print('\n' * 10 + 'Saved!'.center(cc) + '\n' * 18)
                getch()
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
        key = getch()
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
                getch()


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
