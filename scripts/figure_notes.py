"""figure_notes.py — plain-English descriptions of every headline figure.

The one place a figure is explained. `scripts/dashboard.py` shows these under
each figure and `scripts/readme_figures.py` writes them into README.md's
"Results in figures" section, so the two cannot drift apart. Each entry says
what the figure shows, how it was made, and how to read its colours and
marks. Written after looking at the figure (D11); when a figure changes,
look at it again and update its entry here.

Fields: `path`, `task`, `title` (one line), `shows`, `how`, `read` (each a
short paragraph, plain English, no unexplained abbreviations).
"""

from __future__ import annotations

FIGURES: list[dict] = [
    {
        "path": "results/s0_baseline/figures/catalogue_scope.png",
        "task": "S0", "script": "scripts/s0_figures.py",
        "title": "What the project counts as an ion channel, and which protein "
                 "domains are misleading",
        "shows": "A — how many human pore-forming (channel) genes fall into each "
                 "superfamily, a superfamily being a group of channel families "
                 "that share a common ancestor. One superfamily, the P-loop "
                 "channels (voltage-gated potassium, sodium and calcium channels "
                 "and their relatives), holds 143 of the 329 genes. B — every "
                 "protein domain (identified by its Pfam code, e.g. PF00520) "
                 "that occurs in more than one catalogued family, and how many "
                 "families carry it.",
        "how": "Counted from the project's hand-built catalogue of channel "
               "families, after every entry was checked against the live "
               "UniProt and InterPro databases.",
        "read": "A: bar length = number of human genes; bar colour = superfamily "
                "(blue P-loop, orange Cys-loop, green iGluR, violet CLC, greys "
                "for the rest). The numbers at the bar ends read 'census + "
                "control', e.g. '23+9' = 23 channel families plus 9 look-alike "
                "families kept only so they can be recognised and excluded. "
                "B: red bars are domains also found in proteins that are not "
                "channels, so finding that domain does not prove a protein is a "
                "channel; blue bars occur only in channels. PF00520 is found in "
                "20 families, one of them an enzyme.",
    },
    {
        "path": "results/benchmark_controls/figures/benchmark.png",
        "task": "S1", "script": "scripts/s1_figures.py",
        "title": "How well the classifier names known channels, and which test "
                 "made each decision",
        "shows": "A — 71 well-known channel proteins of known family, run through "
                 "the project's classifier, grouped by superfamily. B — the "
                 "'hazards': recorded ways two different families can be "
                 "confused, and how the test proteins that touch each hazard "
                 "were called.",
        "how": "Each protein was classified with itself removed from the "
               "reference set (so it cannot simply match itself). The "
               "classifier has three kinds of test: domain-architecture rules, "
               "a selectivity-filter sequence motif, and percent identity to "
               "reference proteins.",
        "read": "A: each bar is one superfamily's test proteins, split by "
                "outcome — dark blue = named correctly by a domain rule, green = "
                "by the filter motif, light blue = by similarity to a reference "
                "protein, red = named as the wrong family, grey = no family "
                "named. 50 of 71 were named correctly; only 22 of those needed "
                "the reference comparison, so the classifier is not just "
                "finding the nearest known protein. B: blue = protein called "
                "right, red = called wrong. H2 (red) was later rewritten (S2b). "
                "H17–H20 have no bars (marked on the figure) because those hazards "
                "were added after this benchmark; they are tested in the S2d "
                "figure.",
    },
    {
        "path": "results/census_v2/figures/census_v2.png",
        "task": "S2", "script": "scripts/s2_figures.py",
        "title": "The first full census: every database record carrying a channel "
                 "domain",
        "shows": "A — the 1.28 million UniProt protein records that carry at least "
                 "one channel-pore domain, per superfamily, split by whether the "
                 "domain rules could name the family. B — the 16 largest channel "
                 "families by number of records named.",
        "how": "Every protein in UniProt with any of the listed pore domains was "
               "downloaded (the count checked exactly against UniProt's own), "
               "then classified using its domains only.",
        "read": "A: dark blue = named to a family; light blue = placed in the "
                "superfamily but the family could not be told apart. Cys-loop, "
                "iGluR, CLC, DEG/ENaC and P2X are all light blue because every "
                "family in them has the same domains. Records the rules could "
                "not place at all (524,177) are not drawn. B: bar colour = "
                "superfamily (blue P-loop; grey others).",
    },
    {
        "path": "results/census_v3/figures/profiles.png",
        "task": "S3a", "script": "scripts/s3_figures.py",
        "title": "Profile models: checked against the domain rules, then used "
                 "where those rules gave up",
        "shows": "A — for each family, how often the new profile models agree "
                 "with the domain rules where both make a call. B — the records "
                 "that the domain rules could place only in a superfamily, and "
                 "what the profile models then made of them.",
        "how": "One statistical profile (a hidden Markov model built from a "
               "family's aligned sequences) per family; each record is given to "
               "the best-scoring profile only if it clearly beats the "
               "runner-up. Sequences used to build a profile were left out of "
               "the comparison in A.",
        "read": "A: one dot per family; x = how many records the domain rules "
                "named (log scale), y = fraction on which the profile agrees "
                "(note the y-axis runs only from 0.993 to 1: agreement is 99.9 % "
                "overall). Dot colour = superfamily. B: each bar is 100 % of a "
                "superfamily's unresolved records (total at the right); dark "
                "blue = now named to a family, light blue = still superfamily "
                "only, grey = other or unassigned.",
    },
    {
        "path": "results/proteome_scope/figures/panel.png",
        "task": "S4", "script": "scripts/s4_figures.py",
        "title": "The 52 species the census is measured in, and how good their "
                 "data are",
        "shows": "The 50 species whose complete protein sets (reference "
                 "proteomes) form the census's search space, plus two species "
                 "with a genome but no protein set (a land snail, Cornu, and an "
                 "electric ray, Torpedo; not plotted).",
        "how": "Each species' proteome was chosen by a fixed rule and its "
               "download verified by checksum. Quality numbers come from "
               "UniProt and NCBI.",
        "read": "One dot per species. x = BUSCO completeness: the percentage of a "
                "standard set of near-universal genes found in the protein set "
                "(higher = more complete). y = scaffold N50 of the underlying "
                "genome assembly (log scale; higher = longer continuous "
                "stretches of DNA, a less fragmented genome). Colour = lineage "
                "group. The dashed line marks 90 % completeness. Labelled species "
                "are the weak spots: those left of the line may lack genes for "
                "technical reasons, and the sponge's genome is very fragmented, "
                "so an absence in any of them is weaker evidence.",
    },
    {
        "path": "results/panel_sweep/figures/domain_search_missed.png",
        "task": "S3b", "script": "scripts/s3b_figures.py",
        "title": "Channels that a domain-based search would never have found",
        "shows": "Channel proteins in the 50 species that the profile models "
                 "identify confidently but that carry none of the pore domains "
                 "the first census searched for. A — by family. B — by lineage.",
        "how": "Every protein of the 50 proteomes was scored against all family "
               "profiles; high-confidence channel calls absent from the "
               "domain-based census were counted.",
        "read": "A: bar length = number of missed members; the label gives the "
                "share of that family's members in these species (e.g. '76 % of "
                "33' for the proton channel Hv1, MITOK 100 % because it has no "
                "catalogued domain at all). Colour = superfamily. B: the miss "
                "rate per lineage group (missed / all members) — about 1 % in "
                "vertebrates, 12 % in plants: domain databases are built mostly "
                "from well-studied animals.",
    },
    {
        "path": "results/genome_sweep/figures/presence_matrix.png",
        "task": "S5b", "script": "scripts/s5_figures.py",
        "title": "Which channel family is in which species — checked in the "
                 "genome, not just the protein list",
        "shows": "Every channel family (columns, grouped by superfamily with "
                 "vertical rules) in every one of the 52 species (rows, "
                 "bacteria at the top, then fungi, plants, protists, simple "
                 "animals, invertebrates, vertebrates, human, and three viruses "
                 "at the bottom).",
        "how": "Known channel proteins from related species were aligned to "
               "each genome (miniprot, backed up by tblastn) to find genes the "
               "protein list missed. An absence is accepted only if the same "
               "method found the family's relatives reliably in that genome "
               "and the genome is contiguous enough to hold the gene.",
        "read": "Dark blue = found in the species' protein set. Bright blue = "
                "found only in the genome (the protein set missed a real gene). "
                "Pale blue = species with a genome but no protein set (Cornu, "
                "Torpedo). Mid grey = genome match too weak to call; light grey "
                "= partial gene or an assembly gap. Pink = only a fragment "
                "(trace). Red = a controlled absence: the gene is genuinely not "
                "there by every test. Off-white = not informative: no test "
                "could decide (e.g. no related species to search with).",
    },
    {
        "path": "results/alignments/figures/alignments_modules.png",
        "task": "S6", "script": "scripts/s6_figures.py",
        "title": "Building the alignments, and cutting out each channel's pore",
        "shows": "A — which census sequences went into each family's alignment. "
                 "B — how much of each alignment survives trimming. C — whether "
                 "the automatically cut-out pore regions ('modules') match the "
                 "pore region UniProt annotates. D — how accurate the method is "
                 "for families with no annotated pore.",
        "how": "Only confidently assigned, intact sequences were aligned (MAFFT "
               "L-INS-i). The pore module was cut from each sequence by "
               "aligning it to its family profile and taking a fixed stretch "
               "of that profile.",
        "read": "A: dark blue = included; light blue = left out because the "
                "family call was only medium confidence; grey, pink, pale = "
                "left out for other reasons (no profile call, broken genome "
                "gene, identical duplicate). B: one dot per family; x = "
                "alignment length, y = columns kept after trimming (log "
                "scales; dashed line = nothing removed); dot size = number of "
                "sequences; colour = superfamily. C: one dot per checked "
                "protein; 1.0 = the cut-out pore exactly matches UniProt's "
                "annotation (overlap score, Jaccard); 436 of 441 score above "
                "0.8 (dashed line). D: bar = error in placing the pore when "
                "the family's own annotation is hidden; under the dashed line "
                "(12 positions) counts as accurate. Innexins fail (80), so the "
                "method is used only for the six families listed in violet.",
    },
    {
        "path": "results/phylogeny/figures/tier1_trim.png",
        "task": "S7a", "script": "scripts/s7_figures.py",
        "title": "Choosing how to trim alignments before building trees",
        "shows": "For each family alignment, how many informative positions "
                 "(columns that can distinguish between branches of a tree) two "
                 "trimming settings keep.",
        "how": "Alignments contain gappy, unreliable columns that are usually "
               "removed before tree-building. Two settings of the trimAl tool "
               "were compared on the alignments alone — before any tree was "
               "built, so the choice could not be steered by the trees.",
        "read": "One dot per family (colour = superfamily). x = informative "
                "positions kept by the automatic setting used in S6; y = kept "
                "by the 'keep any column at least half filled' setting (both "
                "log scales). Dots above the dashed diagonal mean the second "
                "setting keeps more; it does in 62 of 63 families, so it was "
                "adopted.",
    },
    {
        "path": "results/phylogeny/figures/tier1_trees.png",
        "task": "S7b/S7d", "script": "scripts/s7_figures.py",
        "title": "Family trees: how well supported they are, and whether they "
                 "can be rooted",
        "shows": "A — the trimming choice (as in the S7a figure). B — how "
                 "confident each of the 63 family trees is. C — whether each "
                 "tree's root (its oldest split) is reliable. D — one example "
                 "tree, the ryanodine receptors.",
        "how": "One maximum-likelihood tree per family (IQ-TREE), with 1000 "
               "ultrafast bootstrap replicates (UFBoot): a support score from "
               "0–100 for each branch, where ≥ 95 is conventionally strong. "
               "Trees are rooted by adding a related outgroup family named in "
               "the catalogue in advance; the root sits where the outgroup "
               "joins.",
        "read": "B: one dot per family; x = sequences in the tree, y = fraction "
                "of branches with support ≥ 95; dot size = alignment length; "
                "colour = superfamily. C: bar = support for the root branch "
                "(dashed line = 95); 'occ.' = how much of the alignment the "
                "outgroup sequences fill. Hatched bars = the outgroup did not "
                "stay together as one group, so the root is undefined. The six "
                "potassium-type families re-tested with two new outgroups each "
                "show three bars (grey = original root, then the two new ones) "
                "and are all marked 'unresolved'. Text rows = outgroup of one "
                "sequence, which gives no support value. D: the tree itself; "
                "branch length = amount of sequence change; tip colour = "
                "lineage; black dots = branches with support ≥ 95; rooted on "
                "the three human IP3 receptors (bottom left).",
    },
    {
        "path": "results/auxiliary/figures/auxiliary.png",
        "task": "S20", "script": "scripts/s20_figures.py",
        "title": "Why published counts of human ion channels disagree "
                 "(240–400)",
        "shows": "A — three curated database lists of human ion channels, split "
                 "into what each actually contains. B — real pore-forming genes "
                 "each list leaves out. C — 'auxiliary' subunits (proteins that "
                 "sit on channels but do not form the pore), grouped by true "
                 "relatedness. D — how many such auxiliary proteins the profile "
                 "models find in the 50 species, and how many are genuine.",
        "how": "The three lists (GtoPdb, HGNC, UniProt keyword 'ion channel') "
               "were matched gene by gene to the catalogue. Relatedness in C "
               "comes from all-against-all sequence comparison; in D a hit "
               "counts as genuine if its best match in the human proteome is "
               "that auxiliary group.",
        "read": "A: dark blue = genuine pore-forming channel genes; orange = "
                "auxiliary subunits; other colours = aquaporins, transporters, "
                "enzymes, claudins, pseudogenes and similar non-channels. Totals "
                "at the right. B: one bar per list (dark blue GtoPdb, light "
                "blue HGNC, orange UniProt); e.g. UniProt omits all 21 "
                "connexins. C: each bar is one auxiliary family; separate "
                "coloured segments are groups of proteins unrelated to each "
                "other — 6 of 11 families mix unrelated proteins. D: log scale; "
                "dark blue = genuine; pink = hits whose best human match is an "
                "unrelated protein (shared repeat domains).",
    },
    {
        "path": "results/method_contribution/figures/method_contribution.png",
        "task": "S15", "script": "scripts/s15_figures.py",
        "title": "Which search method finds each channel",
        "shows": "A — every member of the final census, per superfamily, by the "
                 "first method that finds it. B — the curated human channel "
                 "genes, and how many each method finds and names correctly.",
        "how": "Methods applied in order of cost: domain search, then profile "
               "models, then the genome sweep. A member is credited to the "
               "first method that finds it.",
        "read": "A: each bar = 100 % of a superfamily's members (count at the "
                "right). Dark blue = domain search found it and named the right "
                "family; light blue = domain search found it but could not name "
                "the family; orange = only the profile models found it; green "
                "= only the genome sweep found it. Overall domain search finds "
                "94 % but names only 44 %. B: three bars per superfamily — "
                "light blue = found by domain search, dark blue = named "
                "correctly by domain rules, orange = named correctly by "
                "profiles (count of human genes at the right).",
    },
    {
        "path": "results/census_v3/figures/census_r4.png",
        "task": "S2d", "script": "scripts/s3r4_figures.py",
        "title": "Eight newly added channels brought into the census, and "
                 "whether their look-alikes can be told apart",
        "shows": "A — the 26,783 database records added for the newly "
                 "catalogued channel families and their non-channel "
                 "look-alikes (decoys). B — for well-studied (reviewed) "
                 "records, how clearly the correct profile beats the next best "
                 "one.",
        "how": "The new families' domains were searched as an addition to the "
               "existing census, profiles built for each, and every record "
               "assigned to its best profile.",
        "read": "A: bar = records by their final call; 'decoy:' = a "
                "non-channel look-alike; a name in [brackets] = placed in that "
                "superfamily but not named to a family. Colour = confidence of "
                "the profile call (dark blue high, light blue medium, greys "
                "low or none). B: x = margin over the runner-up profile (0 = "
                "tie, 1 = no contest); dashed line = 0.30, the bar for a high "
                "confidence call. Open circles = proteins used to build the "
                "profile; filled = proteins held out as a fair test. Shaded "
                "bands pair each channel with its look-alike (hazards "
                "H17–H19); all separate cleanly.",
    },
    {
        "path": "results/genome_sweep/figures/genome_r4.png",
        "task": "S5c", "script": "scripts/s5r4_figures.py",
        "title": "The seven newly added channel families across the 52 genomes",
        "shows": "Each new family (columns) in each species (rows, bacteria at "
                 "the top to human, viruses at the bottom), checked in the "
                 "genome as in the S5b figure.",
        "how": "Same genome search and absence tests as S5b, run for the new "
               "families only; all earlier results were confirmed unchanged.",
        "read": "Same colours as the S5b matrix: dark blue = in the protein "
                "set; bright blue = only in the genome (protein set missed it); "
                "pale blue = genome-only species; grey = weak or partial; pink "
                "= trace; red = controlled absence (truly not there); "
                "off-white = cannot be decided. TMEM87 is present in every "
                "animal; its presence in fungi, plants and some protists is "
                "real (see the GOST tree figure).",
    },
    {
        "path": "results/phylogeny/gost/figures/gost_tree.png",
        "task": "S2f", "script": "scripts/s7_gost_tree.py",
        "title": "TMEM87 is an ancient lineage found across eukaryotes",
        "shows": "A tree of 115 proteins from the GOST protein superfamily: the "
                 "proposed channel TMEM87 and its non-channel relatives "
                 "(GPR107/108 and the fungal, plant and protist GOST proteins).",
        "how": "Maximum-likelihood tree (IQ-TREE) with 1000 bootstrap "
               "replicates. The tree is unrooted; it is drawn hanging from an "
               "arbitrary point, so left-to-right order does not mean older "
               "to younger.",
        "read": "Tip colour = what the protein is: dark blue = animal TMEM87 "
                "(carrying TMEM87's animal-specific GOLD domain); light blue = "
                "other animal GOST proteins; green = plants and algae; orange "
                "= fungi; violet = single-celled relatives of animals; grey = "
                "other protists. Squares = non-animal proteins the profile "
                "models call TMEM87. Small black dots = branches with support "
                "≥ 95. Branch length = amount of sequence change. The animal "
                "TMEM87s and the squares fall in one strongly supported group "
                "(support 100), apart from GPR107/108.",
    },
    {
        "path": "results/panel_density/figures/order_matrix.png",
        "task": "S4b", "script": "scripts/s4b_report.py",
        "title": "Every channel family across 439 orders of eukaryotes",
        "shows": "Presence of each channel family (rows) in one representative "
                 "protein set per eukaryotic order (columns: 439 orders, "
                 "grouped into animals, fungi, plants and protists).",
        "how": "One reference proteome per order, chosen by a fixed rule, "
               "scanned with all family profiles.",
        "read": "A blue mark = the family was confidently found in that order's "
                "protein set; blank = not found. Long unbroken runs show "
                "families present throughout a kingdom (e.g. VDAC, OSCA, GPHR "
                "in nearly all eukaryotes); runs only in the middle of the "
                "animal block are vertebrate-specific (glycine and 5-HT3 "
                "receptors, connexins, pannexins, CFTR). Blank cells are "
                "weaker evidence than in the S5b figure: a protein set can "
                "simply miss a gene, and these were not checked in the genome.",
    },
    {
        "path": "results/phylogeny/figures/fold_network.png",
        "task": "S8a", "script": "scripts/s8_figures.py",
        "title": "Which channel superfamilies share a 3-D shape",
        "shows": "Most channel superfamilies have no sequence similarity, so "
                 "they cannot be placed in one tree. Instead this compares "
                 "their predicted 3-D structures. A — how similar in shape "
                 "every pair of superfamilies is. B — five relationships the "
                 "literature proposes, tested.",
        "how": "One AlphaFold-predicted structure per family (only the "
               "confidently predicted parts), cut to the part each comparison "
               "is about (the pore, or the voltage sensor), compared pairwise "
               "with TM-align. TM-score runs from 0 (unrelated shapes) to 1 "
               "(identical); above about 0.5 usually means the same fold. The "
               "pass rule was fixed before measuring: median score ≥ 0.5 and "
               "each side's closest match among all other superfamilies is the "
               "other side.",
        "read": "A: each cell = the median similarity between two "
                "superfamilies' families (darker blue = more similar; scale "
                "0–1). Diagonal cells compare different families within the "
                "same superfamily (a structure is never compared with itself); "
                "hatched = only one family, so nothing to compare. Rows are "
                "ordered so similar superfamilies sit together. Boxed cells = "
                "the literature's proposed relationships (violet box = "
                "supported, black = not distinguished). B: each dot = one "
                "family pair behind a proposed relationship (blue = supported, "
                "grey = not); black line = median; violet dashed line = the "
                "best score either side reaches with anything else; grey "
                "dashed line = the 0.5 bar. Supported: glutamate-receptor pore "
                "vs potassium-channel pore, Hv1 vs the Kv voltage sensor, "
                "innexins vs connexins.",
    },
]


def by_path() -> dict[str, dict]:
    return {f["path"]: f for f in FIGURES}


def caption_text(f: dict) -> str:
    """One plain-text block (dashboard)."""
    return (f"{f['task']} — {f['title']}. What it shows: {f['shows']} "
            f"How it was made: {f['how']} How to read it: {f['read']}")


def caption_markdown(f: dict) -> str:
    """README block: heading, image, three labelled paragraphs."""
    return (f"**{f['task']} — {f['title']}.**\n\n"
            f"![{f['title']}]({f['path']})\n\n"
            f"*What it shows.* {f['shows']}\n\n"
            f"*How it was made.* {f['how']}\n\n"
            f"*How to read it.* {f['read']}\n\n"
            f"<sub>Drawn by `{f['script']}` from the task's committed tables.</sub>\n")
