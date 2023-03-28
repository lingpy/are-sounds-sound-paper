from lingreg.util import prep_wordlist, prep_alignments
from lingpy import *
from collections import defaultdict

def run(wordlist):
    wordlist.renumber("cog")
    # get only those languages which have a glottocode
    lng2gl = defaultdict(set)
    for idx, language, glottocode in wordlist.iter_rows("doculect",
                                                        "glottocode"):
        lng2gl[language] = glottocode

    selected = [lng for lng in lng2gl if lng2gl[lng]]
    print(selected)
    D = {0: wordlist.columns}
    for idx in wordlist:
        if wordlist[idx, "doculect"] in selected:
            D[idx] = wordlist[idx]
    wl = Wordlist(D)
    wordlist = prep_wordlist(wl)
    alms = Alignments(wordlist, ref="cogid", transcription="form")
    alms = prep_alignments(alms)
    alms.align()
    return alms
