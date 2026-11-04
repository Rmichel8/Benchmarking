import pandas as pd
import re

def canon(s: str) -> str:
    s = str(s).strip().lower()
    s = re.sub(r"[^a-z0-9\-\(\)\. ]+", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s

# alias 
TRUTH_ALIAS = {
    "escherichia coli amph beta lactamase": "amph",
    "escherichia coli uhpt with mutation conferring resistance to fosfomycin": "uhpt",
    "morganella morganii gyrb conferring resistance to fluoroquinolone": "gyrb",
    "haemophilus influenzae pbp3 conferring resistance to beta lactam antibiotics": "pbp3",
    "escherichia coli ef-tu mutants conferring resistance to pulvomycin": "ef-tu",
    "escherichia coli ef tu mutants conferring resistance to pulvomycin": "ef-tu",
    "klebsiella pneumoniae kpne": "kpne",
    "klebsiella pneumoniae kpnf": "kpnf",
    "klebsiella pneumoniae kpnh": "kpnh",
    "acinetobacter baumannii abaq": "abaq",
    "pseudomonas aeruginosa soxr": "soxr",

}

PRED_ALIAS = {
    # variants fosB
    "fosb2": "fosb",
    "fosbx1": "fosb",

    # organism + gene →  gene
    "klebsiellapneumoniaekpne": "kpne",
    "klebsiella pneumoniae kpne": "kpne",

    "klebsiellapneumoniaekpnf": "kpnf",
    "klebsiella pneumoniae kpnf": "kpnf",

    "klebsiellapneumoniaekpnh": "kpnh",
    "klebsiella pneumoniae kpnh": "kpnh",
    # beta-lactamase
    "amphbeta-lactamase": "amph",
    "amph beta-lactamase": "amph",
    "escherichiacoliamph": "amph",

    # other organism prefix
    "pseudomonasaeruginosacpxr": "cpxr",
    "pseudomonas aeruginosa cpxr": "cpxr",
    "pseudomonasaeruginosasoxr": "soxr",
    "gyrbconferringresistancetofluoroquinolones": "gyrb",
    "pbp3conferringresistancetobeta-lactamantibiotics": "pbp3",
    "staphylococcusaureusfosb": "fosb",
    "enterobactercloacaeacra": "acra",
    "escherichiacoliacra": "acra",
    "klebsiellapneumoniaeacra": "acra",
    "shigellaflexneriacra": "acra",
    "haemophilusparainfluenzaegyra": "gyra",
    "capnocytophagagingivalisgyra": "gyra",
    "burkholderiadolosagyra": "gyra",
    #NanoBsi RGI
    "elfamycinresistantef-tu": "ef-tu",
    "fosfomycisthioltransferase": "fosA8",
}


def truth_to_token(gene_symbol: str) -> str:
    s = canon(gene_symbol)
    if s in TRUTH_ALIAS:
        return TRUTH_ALIAS[s]
    if len(s.split()) == 1:
        return s

    # take the last “token gene-like” 
    generic_tail = {"betalactamase","beta-lactamase", "mutants", "resistance", "conferring", "antibiotics", "antibiotic"}
    parts = s.split()
    for tok in reversed(parts):
        if tok in generic_tail:
            continue
        if re.search(r"[a-z]", tok):
            return tok

    # fallback
    return parts[-1]

def pred_to_token(gene_symbol: str) -> str:
    s = str(gene_symbol).strip().lower()
    s = re.sub(r"[^a-z0-9\-\(\)\.]+", "", s)

    if s in PRED_ALIAS:
        return PRED_ALIAS[s]

    return s



