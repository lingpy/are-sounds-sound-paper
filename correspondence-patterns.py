from lingpy import *
from glob import glob
from lingrex.copar import CoPaR
from lingreg.trim import trim_by_gap, trim_by_core, get_skeleton, apply_trim, trim_random


def trim_alignments(
        wordlist,
        trim_fun,
        pattern_threshold=3,
        word_threshold=0.75,
        ref="cogid",
        gap_threshold=0.5,
        skeletons=("CV", "VC")
        ):
    """
    Function takes a wordlist and then applies the trimming procedure.
    """
    new_alms = {}
    for cogid, msa in wordlist.msa[ref].items():
        trimmed = apply_trim(
                msa["alignment"],
                trim_fun(
                    msa["alignment"],
                    threshold=gap_threshold,
                    skeletons=skeletons
                    )
                )
        for idx, row in zip(msa["ID"], trimmed):
            new_alms[idx] = row
    wordlist.add_entries("original_alignment", "alignment", lambda x: x)
    wordlist.add_entries("original_tokens", "tokens", lambda x: x)
    wordlist.add_entries("original_structure", "structure", lambda x: x)
    wordlist.add_entries("alignment", new_alms, lambda x: x, override=True)
    wordlist.add_entries(
            "tokens",
            new_alms, lambda x: [c for c in x if c != "-"],
            override=True)
    wordlist.add_entries(
            "structure",
            new_alms, lambda x: get_skeleton([[c for c in x if c != "-"]]),
            override=True
            )


files = glob("datasets/*.tsv")
for f in files:

    alms = Alignments(f, ref="cogid", transcription="form")
    trim_alignments(alms, trim_by_gap)
    alms.output("tsv", filename="trimmed/"+f.split("/")[-1][:-4],
                prettify=False)
    cop = CoPaR("trimmed/"+f.split("/")[-1], alignment="alignment", ref="cogid", transcription="form")
    cop.get_sites()
    cop.cluster_sites()
    cop.sites_to_pattern()
    cop.add_patterns()
    cop.write_patterns("correspondences/"+f.split("/")[-1])

