"""figure_notes_late.py — plain-English figure descriptions from S9 on.

`figure_notes.py` reached the 500-line budget; it appends this module's
`FIGURES` to its own, so README and dashboard read one list.
"""

FIGURES = [
    {
        "path": "results/filter_atlas/figures/filter_atlas.png",
        "task": "S9", "script": "scripts/s9_figures.py",
        "title": "Do selectivity filters follow the family tree?",
        "shows": "The selectivity filter, the few residues of the pore that decide "
                 "which ion passes, read in every P-loop channel and laid on the "
                 "family trees. A — the calcium-channel (Cav) tree with every "
                 "sequence coloured by its filter. B — which filters the "
                 "four-repeat channels (sodium, calcium, NALCN, CatSper, "
                 "two-pore) carry. C — for every family, how closely each filter "
                 "position follows the tree compared with every other position "
                 "in the same alignment.",
        "how": "Every pore region was aligned onto one shared reference alignment "
               "(MAFFT), and the filter read at the columns occupied by the "
               "potassium channel KcsA's TVGYG motif and by the four filter "
               "residues of the human heart sodium channel Nav1.5 (DEKA). The "
               "reads were checked against the classifier's independent "
               "method (98 % agreement on 363 four-repeat channels). Each "
               "filter was then scored on the maximum-likelihood trees by "
               "parsimony: how many times it must have changed, compared with "
               "every other alignment column (the retention index, where 1 "
               "means each variant arose once).",
        "read": "A: each horizontal line ending at the right is one sequence; "
                "horizontal distance is substitutions per site, and the three "
                "light-grey lines at the bottom are the bacterial potassium "
                "channels used as the outgroup. Blue = EEEE (the high-voltage "
                "calcium channels), orange = EEDD (the low-voltage T-type "
                "channels), green = DDDD (all from the ciliate Paramecium), "
                "violet = any other filter, light grey = a repeat that could "
                "not be read. Forty-seven of the 49 EEDD sequences form the one "
                "orange block, a branch with bootstrap support 100: the T-type "
                "filter arose once. B: each bar is one family (number of "
                "sequences in brackets), split into its most common filters "
                "from dark to light blue, then 'other' (off-white) and unread "
                "(grey); the letters are the residues of repeats I to IV "
                "(CatSper and the two-pore channels have one and two repeats). "
                "The unlabelled second segments are DEEA (Nav), EKEE (NALCN). "
                "C: each dot is one filter position in one family, placed by "
                "its percentile among all the alignment's columns; right of "
                "the 50 line = follows the tree more closely than a typical "
                "position. Orange = a four-repeat filter position; blue = the "
                "variable middle residue of the potassium motif TxGYG; green = "
                "its Y/F; grey = other positions of the window. Four-repeat "
                "positions and the Y/F sit mostly far right; the middle "
                "residue x is often far left, i.e. it has changed many times "
                "independently.",
    },
]
