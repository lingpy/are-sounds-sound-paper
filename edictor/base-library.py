from lingreg.util import prep_wordlist, prep_alignments
from lingpy import *

def run(wordlist):
    wordlist.renumber("cog")
    wordlist = prep_wordlist(wordlist)
    alms = Alignments(wordlist, ref="cogid", transcription="form")
    alms = prep_alignments(alms)
    alms.align(method="library")
    return alms
