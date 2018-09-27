#!/usr/bin/env python

# Random nucleotide sequence generator
# with given lenght, AT% and GC skew

from sys import argv
from sympy.solvers import solve
from sympy import Symbol
import numpy as np

display_help = 0
user_input = 0

LEN = 100
AT = 0.5 # AT content in %
REP = 1
ACID = "DNA"
_nucleotides = ['A', 'T', 'C', 'G']

SKEW = 0

# COMMAND LINE ARGUMENTS
if (len(argv) > 1):
    argv.pop(0) # remove first item from list as it is path to itself
    for arg in argv: # loop thru all arguments
		if "-n=RNA" in arg:
			try:
				_nucleotides = ['A', 'U', 'C', 'G']
				ACID = "RNA"
			except: continue

		elif "-at=" in arg:
			arg = arg.replace("-at=", "")
			try:
				_at = float(arg)
				AT = _at
			except: continue

		elif "-l=" in arg:
			arg = arg.replace("-l=", "")
			try:
				_L = int(arg)
				LEN = _L
			except: continue

		elif "-r=" in arg:
			arg = arg.replace("-r=", "")
			try:
				_R = int(arg)
				REP = _R
			except: continue

		elif "-s=" in arg:
			arg = arg.replace("-s=", "")
			try:
				_S = float(arg)
				SKEW = _S
			except: continue

		elif ("-h" in arg or "-help" in arg or "--h" in arg or "--help" in arg) is True:
			display_help = 1

### HELP
if display_help == 1:
    print "\n Random nucleotide sequence generator\n\
 Zbyszek Pietras 2016 zp@cantab.net\n\
 Use of command line argments:\n\
 -n=[RNA/DNA] type of nucleic acid (default value: %s)\n\
 -at=[float] %% of AT (default value: %.2f)\n\
 -r=[integer] No. of sequences to generate (default value: %d)\n\
 -l=[integer] sequence lenght (default value: %d)\n\
 -s=[float] GC Skew = (G - C)/(G + C) (default value: %.2f)\n" % (ACID, AT, REP, LEN, SKEW)
    exit()
### /HELP

lista__prcntAT = []
lista_GCskew = []

# here we will calculate multipliers for GC skew
if SKEW == 0:
	_c = 0.5
	_g = 0.5
else:
	C_ = Symbol('C_')
	G_ = Symbol('G_')
	_c = solve(((1 - C_) - C_)/((1 - C_) + C_) - SKEW)
	_g = solve((G_ - (1 - G_))/(G_ + (1 - G_)) - SKEW)
	_c = _c[0] # solve function returns a list we need to change it into number
	_g = _g[0]

# remember! 'A' 'T/U' 'C' 'G'
_pro = [(AT/2), (AT/2), ((1-AT)*_c), ((1-AT)*_g)]
# print _pro

while REP > 0:
	_acid = ''

### below is commendet out old and slow generator
#	while len(_acid) < LEN:
#		_nuc = random.choice(_nucleotides)
#		_razor = random.random()
#		if _razor < _pro[_nuc]:
#			_acid += _nuc

### here is new and fast generator using numpy
	_acid = ''.join(np.random.choice(_nucleotides, size=LEN, replace=True, p=_pro))

    # here we calculate AT%, GCskew i C/G ratio for each sequence and it will be added to name of fasta sequence
	_prcntAT = float(((float(_acid.count("A"))+float(_acid.count("T")))/float(LEN))*100)
	_gs = _acid.count("G")
	_cs = _acid.count("C")
	if (_gs + _cs) == 0:
		_GCskew = "noGsCs"
	else:
		_GCskew = float((_gs - _cs))/(_gs + _cs)

	if _gs == 0:
		_CoverG = "noGs"
	else:
		_CoverG = float(_cs)/_gs
	print ">%d_pAT=%.1f_GCskew=%.3f_CtoG=%.2f" % (REP, _prcntAT, _GCskew, _CoverG)
	REP = REP-1
	print _acid

### code below is to check the accuracy of sequence generation
# 	lista__prcntAT.append(_prcntAT)
# 	try:
# 		lista_GCskew.append(float(_GCskew))
# 	except: continue
#
# print "requested AT:%.1f obtained: %.1f delta: %.1f"  % ((AT*100), np.mean(lista__prcntAT), ((AT*100)-(np.mean(lista__prcntAT))))
# print "requested GC Skew:%.1f obtained:%.1f delta: %.1f"  % (SKEW, np.mean(lista_GCskew), (SKEW-(np.mean(lista_GCskew))))
