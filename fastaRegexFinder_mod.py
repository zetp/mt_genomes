#!/usr/bin/env python

# This is a mofdified version of script from Dario Beraldi
# https://github.com/dariober/bioinformatics-cafe/tree/master/fastaRegexFinder
# under MIT license

import re
import sys
import string
import argparse
import operator

VERSION='0.1.1'

parser = argparse.ArgumentParser(description="""
DESCRIPTION

    Search a multi-fasta file for matches to a regex and return results.
    This scirpt is modified to be used with accompanied script that generates random sequences,
    which is source of multifasta as we can also provide the data about numer of REGEX occurneces in
    the original genome for comparison.
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--fasta', '-f',
                   type= str,
                   help='''Input fasta file to search. Use '-' to read the file from stdin.

                   ''',
                   required= True)

parser.add_argument('--regex', '-r',
                   type= str,
                   help='''Regex to be searched in the fasta input.
Matches to the reverse complement will have - strand.
The default regex is '([gG]{3,}\w{1,7}){3,}[gG]{3,}' which searches
for G-quadruplexes.
                   ''',
                   default= '([gG]{3,}\w{1,7}){3,}[gG]{3,}')

parser.add_argument('--matchcase', '-m',
                   action= 'store_true',
                   help='''Match case while searching for matches. Default is
to ignore case (I.e. 'ACTG' will match 'actg').
                   ''')

parser.add_argument('--noreverse',
                   action= 'store_true',
                   help='''Do not search the reverse complement of the input fasta.
Use this flag to search protein sequences.
                   ''')

parser.add_argument('--maxstr',
                   type= int,
                   required= False,
                   default= 10000,
                   help='''Maximum length of the match to report in the 7th column of the output.
Default is to report up to 10000nt.
Truncated matches are reported as <ACTG...ACTG>[<maxstr>,<tot length>]
                   ''')

parser.add_argument('--zzz', '-z',
                   type= int,
                   required= False,
                   default= 5,
                   help='''How many times G4 was found in original genome
                   ''')

parser.add_argument('--seqnames', '-s',
                   type= str,
                   nargs= '+',
                   default= [None],
                   required= False,
		   help='''List of fasta sequences in --fasta to
search. E.g. use --seqnames chr1 chr2 chrM to search only these crhomosomes.
Default is to search all the sequences in input.
                   ''')
parser.add_argument('--quiet', '-q',
                   action= 'store_true',
                   help='''Do not print progress report (i.e. sequence names as they are scanned).
                   ''')

parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + VERSION)

args = parser.parse_args()

" --------------------------[ Check and parse arguments ]---------------------- "

if args.matchcase:
    flag= 0
else:
    flag= re.IGNORECASE

" ------------------------------[  Functions ]--------------------------------- "

def sort_table(table, cols):
    """ Code to sort list of lists
    see http://www.saltycrane.com/blog/2007/12/how-to-sort-table-by-columns-in-python/
    sort a table by multiple columns
        table: a list of lists (or tuple of tuples) where each inner list
               represents a row
        cols:  a list (or tuple) specifying the column numbers to sort by
               e.g. (1,0) would sort by column 1, then by column 0
    """
    for col in reversed(cols):
        table = sorted(table, key=operator.itemgetter(col))
    return(table)

def trimMatch(x, n):
    """ Trim the string x to be at most length n. Trimmed matches will be reported
    with the syntax ACTG[a,b] where Ns are the beginning of x, a is the length of
    the trimmed strng (e.g 4 here) and b is the full length of the match
    EXAMPLE:
        trimMatch('ACTGNNNN', 4)
        >>>'ACTG[4,8]'
        trimMatch('ACTGNNNN', 8)
        >>>'ACTGNNNN'
    """
    if len(x) > n and n is not None:
        m= x[0:n] + '[' + str(n) + ',' + str(len(x)) + ']'
    else:
        m= x
    return(m)

def revcomp(x):
    """Reverse complement string x. Ambiguity codes are handled and case conserved.

    Test
    x= 'ACGTRYSWKMBDHVNacgtryswkmbdhvn'
    revcomp(x)
    """
    compdict=  {'A':'T',
                'C':'G',
                'G':'C',
                'T':'A',
                'R':'Y',
                'Y':'R',
                'S':'W',
                'W':'S',
                'K':'M',
                'M':'K',
                'B':'V',
                'D':'H',
                'H':'D',
                'V':'B',
                'N':'N',
                'a':'t',
                'c':'g',
                'g':'c',
                't':'a',
                'r':'y',
                'y':'r',
                's':'w',
                'w':'s',
                'k':'m',
                'm':'k',
                'b':'v',
                'd':'h',
                'h':'d',
                'v':'b',
                'n':'n'}
    xrc= []
    for n in x:
        xrc.append(compdict[n])
    xrc= ''.join(xrc)[::-1]
    return(xrc)
# -----------------------------------------------------------------------------
LISTA = []
LP=0

psq_re_f= re.compile(args.regex, flags= flag)
## psq_re_r= re.compile(regexrev)

if args.fasta != '-':
    ref_seq_fh= open(args.fasta)
else:
    ref_seq_fh= sys.stdin

ref_seq=[]
line= (ref_seq_fh.readline()).strip()
chr= re.sub('^>', '', line)
line= (ref_seq_fh.readline()).strip()
gquad_list= []
while True:
#    if not args.quiet:
#        sys.stderr.write('Processing %s\n' %(chr))
    while line.startswith('>') is False:
        ref_seq.append(line)
        line= (ref_seq_fh.readline()).strip()
        if line == '':
            break
    ref_seq= ''.join(ref_seq)
    if args.seqnames == [None] or chr in args.seqnames:
        for m in re.finditer(psq_re_f, ref_seq):
            matchstr= trimMatch(m.group(0), args.maxstr)
            quad_id= str(chr) + '_' + str(m.start()) + '_' + str(m.end()) + '_for'
            gquad_list.append([chr, m.start(), m.end(), quad_id, len(m.group(0)), '+', matchstr])
        if args.noreverse is False:
            ref_seq= revcomp(ref_seq)
            seqlen= len(ref_seq)
            for m in re.finditer(psq_re_f, ref_seq):
                matchstr= trimMatch(revcomp(m.group(0)), args.maxstr)
                mstart= seqlen - m.end()
                mend= seqlen - m.start()
                quad_id= str(chr) + '_' + str(mstart) + '_' + str(mend) + '_rev'
                gquad_list.append([chr, mstart, mend, quad_id, len(m.group(0)), '-', matchstr])
        gquad_sorted= sort_table(gquad_list, (1,2,3))
        gquad_list= []
###        print "Processing %s found:\t%d" %(chr[:40], len(gquad_sorted))
        LISTA.append(len(gquad_sorted))
        LP = LP+1
        for xline in gquad_sorted:
            xline= '\t'.join([str(x) for x in xline])
#            print(xline)
    chr= re.sub('^>', '', line)
    ref_seq= []
    line= (ref_seq_fh.readline()).strip()
    if line == '':
        break

print "regex: %s" % (args.regex)
print "analysed sequences: %d" % (LP)
print "regex found: %d" % (sum(LISTA))
print "average: %f" % (sum(LISTA)/float(LP))
print "G4 was found at least %d times in %d sequences out of %d total sequences" % (args.zzz, sum(i >= args.zzz for i in LISTA), len(LISTA))
print "found G4 (No per sequence):"
print LISTA


#gquad_sorted= sort_table(gquad_list, (0,1,2,3))
#
#for line in gquad_sorted:
#    line= '\t'.join([str(x) for x in line])
#    print(line)
sys.exit()
