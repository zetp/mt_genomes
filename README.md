# Collection of scripts to analyse G-quadruplex sequences in mitochondrial genomes.

Some of those were used to generate analyses published in Nature Communications article
"[Dedicated surveillance mechanism controls G-quadruplex forming non-coding RNAs in human mitochondria](https://doi.org/10.1038/s41467-018-05007-9)"

both scripts can be linked using command line so one can avoid writing big multifasta
files onto disk:
'''
RNSG.py -s=0.6 -l=10000 -r=100 | fastaRegexFinder_mod.py -f -
'''
