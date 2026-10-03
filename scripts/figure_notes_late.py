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
{
        "path": "results/repertoire/figures/repertoire.png",
        "task": "S10a", "script": "scripts/s10_figures.py",
        "title": "Which ion channels did each ancestor carry?",
        "shows": "For each of the 75 census channel families, whether it was present "
                 "in fifteen ancestors on the tree of eukaryotes, from the last "
                 "eukaryotic common ancestor (LECA) to the ancestor of mammals, and "
                 "how often it was gained and lost. A — the reconstructed state at "
                 "each ancestor. B — the number of gains and losses placed on the "
                 "tree for each family. C — every inferred loss sorted by how well "
                 "it is backed by the genome search.",
        "how": "Each family was scored present or absent in one reference proteome "
               "for each of 437 eukaryotic orders (present = a high-confidence "
               "profile match; a medium-confidence match, or an absence in a "
               "proteome less than 70 % complete by BUSCO, counts as unknown). The "
               "states were laid on the NCBI Taxonomy tree, unresolved branchings "
               "kept as they are, and the fewest gains and losses that explain "
               "them were found by parsimony, with a gain costing two losses. The "
               "same was repeated with a gain costing one loss and with only one "
               "gain allowed, to see which conclusions depend on that choice. "
               "Losses were then compared with the genome search of the 37 "
               "panel species that sit in these orders.",
        "read": "A: rows are families, grouped by superfamily (black lines); "
                "columns are ancestors. Dark blue = present, pale blue = absent, "
                "grey = the data fit present and absent equally well. A dot marks "
                "an ancestor whose state changes under one of the other two "
                "costs; cells without a dot hold under all three. B: bars to the "
                "left are gains, to the right losses (blue and orange); the solid "
                "part is placed identically under all three costs, the pale part "
                "only under the main one. A family can be present at LECA (A) and "
                "still show no gain in B, because its origin before LECA is not a "
                "branch of this tree. C: one bar of all 673 losses inferred under "
                "the main cost: grey = seen only as missing from proteomes; green "
                "= confirmed by a controlled genome absence in a species below "
                "the loss; orange = contradicted, because the genome search found "
                "an intact gene the proteome lacks. Most losses are grey: they "
                "are annotation-level absences, not proven gene losses.",
    },
    {
        "path": "results/repertoire/figures/s10b_checks.png",
        "task": "S10b", "script": "scripts/s10b_figures.py",
        "title": "Do animals carry bacterial-type channels, and which absences hold?",
        "shows": "Four checks on the repertoire reconstruction. A — every animal gene "
                 "or protein that our profiles call a bacterial-type channel "
                 "(the mechanosensitive channels MscS and MscL, the potassium "
                 "channel KcsA, the glutamate receptor GluR0), sorted by what the "
                 "other genes on the same stretch of assembled DNA look like. "
                 "B — the water flea Daphnia pulex searched on its 2011 and its "
                 "2021 genome assembly. C — the genome search's absences of three "
                 "channels (ZAC, PACC1, CLCC1) compared with two public ortholog "
                 "databases. D — how many orders in each kingdom carry a two-pore "
                 "potassium (K2P) channel call.",
        "how": "A: each sequence's genomic contig (a continuous piece of the "
               "assembly) was compared by translated BLAST (blastx) with the "
               "proteomes of the 50 panel species, its own species removed; each "
               "other gene on the contig was assigned to animals, other eukaryotes "
               "or bacteria by its best match, only when that match beat the best "
               "match from any other group by at least 10 %. A gene lying among "
               "animal genes is likely a real animal gene; one among bacterial "
               "genes is likely contamination from a symbiont. B: the project's "
               "genome search (protein-to-genome alignment of known channels, "
               "then profile scoring) was rerun unchanged on the newest Daphnia "
               "assembly. C: NCBI Gene and Ensembl were asked for orthologs of the "
               "human gene; any listed protein was scored with our profiles and "
               "located on our genome. D: counts from the order-level proteome "
               "panel.",
        "read": "A: one bar per channel type; MscS is split into high- and "
                "medium-confidence calls. Blue = at least one neighbouring gene is "
                "animal (embedded); orange = the neighbours are non-animal "
                "(foreign); light grey = no other gene on the contig; numbers are "
                "counts. No high-confidence animal MscS is embedded. B: one row "
                "per channel family whose verdict was 'absent' or changed; left "
                "square = 2011 assembly, right = 2021. Orange = absent (a "
                "controlled absence); dark grey = partial (a weak alignment no "
                "profile scores); blue = found; pale blue = found at medium "
                "confidence; light grey = nothing found and nothing to compare. "
                "C: one row per absence; green = neither database lists an "
                "ortholog (agrees), orange = a database lists one (disputed), "
                "grey = neither database covers that species; the text gives each "
                "database's answer. D: bars are orders; blue = K2P called at high "
                "confidence, orange = no call in a complete proteome, grey = "
                "unknown (medium call or incomplete proteome). 'Protists' are "
                "eukaryotes outside animals, fungi and plants.",
    },
    {
        "path": "results/duplication/figures/duplications.png",
        "task": "S11a", "script": "scripts/s11_figures.py",
        "title": "When were the ion-channel genes duplicated?",
        "shows": "Gene duplications in 63 ion-channel family trees, placed on the tree of "
                 "the 52 panel species. A — duplications whose two copies are found in "
                 "more than one species, by the group of species they date to. B — "
                 "duplications inside a single species. C — human gene pairs that an "
                 "independent database (OHNOLOGS v2) attributes to the two whole-genome "
                 "duplications at the origin of vertebrates ('2R'), and the age our gene "
                 "trees give them. D — how each family tree was rooted.",
        "how": "Each family's maximum-likelihood tree was compared with the NCBI "
               "taxonomy tree of the panel species (reconciliation): a node is a "
               "duplication when the species below its two branches overlap, and it is "
               "dated to the smallest group of species containing all its descendants. "
               "Trees without a trusted outgroup root were rooted where the fewest "
               "duplications (then losses) are needed; a duplication is counted only if "
               "it holds under every such root and both branches below it have "
               "bootstrap support of at least 95 %. Human pairs were matched to OHNOLOGS "
               "through HGNC gene identifiers; the database never changed a call.",
        "read": "A and B: bar length = number of supported duplications, with the count "
                "at the end. In A, dark blue = the two vertebrate whole-genome "
                "duplications' window (Vertebrata, Gnathostomata = jawed vertebrates), "
                "green = the teleost-fish genome duplication's window (Clupeocephala = "
                "the group holding zebrafish and pufferfish), light blue = any other "
                "group. In B, blue = duplications among proteome genes, grey = ones "
                "involving a locus found only in the genome (Cornu, the garden snail, "
                "has no proteome). C: one bar per family with at least 4 such pairs; "
                "blue = our tree dates the pair to the vertebrate window (agreement), "
                "orange = our tree dates it older than vertebrates (usually because a "
                "non-vertebrate sequence sits among the vertebrate copies), grey = "
                "other. D: number of trees; 'declared root' = rooted on the outgroup "
                "named in the catalogue, blue where that root is also a "
                "fewest-duplication root, orange where it is not; light blue / grey = "
                "rooted by fewest duplications, at one edge or tied between several.",
    },
]
