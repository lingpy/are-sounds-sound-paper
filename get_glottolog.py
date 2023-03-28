from lingpy import *
from glob import glob

files = glob("datasets/*.tsv")
for f in files:

    wl = Wordlist(f)
    l2g = {wl[idx, "doculect"]: wl[idx, "glottocode"] for idx in wl}

    with open(f.replace("datasets/", "languages/"), "w") as f:
        f.write("Language\tGlottocode\n")
        for k, v in sorted(l2g.items()):
            f.write(k+'\t'+v+'\n')
