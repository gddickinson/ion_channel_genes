"""The review's source list — titles to resolve, not a bibliography.

Every entry is `(key, title)` or `(key, title, constraint)`, where the
constraint is extra Europe PMC query syntax — `AUTH:"Catterall WA" AND
PUB_YEAR:2011`. Short generic titles ("TRP channels", "Store-Operated
Calcium Channels") match hundreds of papers and the search ranks the wrong
one; the constraint is how a real paper with an unhelpfully short title gets
resolved without loosening the title match that keeps the bibliography
honest. Nothing here is a citation yet: the titles are
*queries*, and `scripts/s0_review_refs.py` resolves each one against Europe
PMC and writes the verified metadata to
`results/s0_baseline/references.tsv`. A title that does not resolve, or that
resolves to something whose title does not match, is reported and **not**
written — so a misremembered paper becomes a build failure rather than a
fabricated reference.

This matters more than usual here. A review is the one artefact in a
publication project where an invented citation is both easy to produce and
almost impossible for a reader to catch, so the rule is absolute: the
bibliography is machine-generated from live records, and prose may cite only
keys that survived resolution (`scripts/s0_review_build.py --check`).

Keys are stable and arbitrary. The build renumbers them into order of first
appearance, so adding a paper mid-section does not renumber the file.
"""

from __future__ import annotations

SOURCES: list[tuple[str, str]] = [
    # ---------------------------------------------------- foundations
    ("hodgkin1952", "A quantitative description of membrane current and its application to conduction and excitation in nerve"),
    ("neher1976", "Single-channel currents recorded from membrane of denervated frog muscle fibres"),
    ("catterall2000", "From ionic currents to molecular mechanisms: the structure and function of voltage-gated sodium channels"),
    ("yu2004", "The VGL-chanome: a protein superfamily specialized for electrical signaling and ionic homeostasis"),
    ("yellen2002", "The voltage-gated potassium channels and their relatives"),
    # -------------------------------------------- potassium and the filter
    ("doyle1998", "The structure of the potassium channel: molecular basis of K+ conduction and selectivity"),
    ("heginbotham1994", "Mutations in the K+ channel signature sequence"),
    ("zhou2001", "Chemistry of ion coordination and hydration revealed by a K+ channel-Fab complex at 2.0 A resolution"),
    ("mackinnon2003", "Potassium channels and the atomic basis of selective ion conduction (Nobel Lecture)"),
    ("long2005", "Crystal structure of a mammalian voltage-dependent Shaker family K+ channel"),
    ("papazian1987", "Cloning of genomic and complementary DNA from Shaker, a putative potassium channel gene from Drosophila"),
    ("tempel1987", "Sequence of a probable potassium channel component encoded at Shaker locus of Drosophila"),
    ("kuo2003", "Crystal structure of the potassium channel KirBac1.1 in the closed state"),
    ("whorton2011", "Crystal structure of the mammalian GIRK2 K+ channel and gating regulation by G proteins, PIP2, and sodium"),
    ("brohawn2012", "Crystal structure of the human K2P TRAAK, a lipid- and mechano-sensitive K+ ion channel"),
    ("enyedi2010", "Molecular background of leak K+ currents: two-pore domain potassium channels"),
    ("hibino2010", "Inwardly rectifying potassium channels: their structure, function, and physiological roles"),
    ("jiang2002", "Crystal structure and mechanism of a calcium-gated potassium channel"),
    # --------------------------------------------- sodium, calcium, leak
    ("noda1984", "Primary structure of Electrophorus electricus sodium channel deduced from cDNA sequence"),
    ("payandeh2011", "The crystal structure of a voltage-gated sodium channel"),
    ("heinemann1992", "Calcium channel characteristics conferred on the sodium channel by single mutations"),
    ("wu2015", "Structure of the voltage-gated calcium channel Cav1.1 complex"),
    ("catterall2011", "Voltage-gated calcium channels", 'AUTH:"Catterall" AND PUB_YEAR:2011'),
    ("lu2007", "The neuronal channel NALCN contributes resting sodium permeability and is required for normal respiratory rhythm"),
    ("ren2001", "A sperm ion channel required for sperm motility and male fertility"),
    ("kirichok2006", "Whole-cell patch-clamp measurements of spermatozoa reveal an alkaline-activated Ca2+ channel"),
    ("calcraft2009", "NAADP mobilizes calcium from acidic organelles through two-pore channels"),
    # ------------------------------------------------------- CNG and HCN
    ("kaupp2002", "Cyclic nucleotide-gated ion channels", 'AUTH:"Kaupp" AND PUB_YEAR:2002'),
    ("robinson2003", "Hyperpolarization-activated cation currents: from molecules to physiological function"),
    ("lee2017", "Structures of the human HCN1 hyperpolarization-activated channel"),
    # ------------------------------------------------------------- TRP
    ("montell1989", "Molecular characterization of the Drosophila trp locus: a putative integral membrane protein required for phototransduction"),
    ("caterina1997", "The capsaicin receptor: a heat-activated ion channel in the pain pathway"),
    ("liao2013", "Structure of the TRPV1 ion channel determined by electron cryo-microscopy"),
    ("clapham2003", "TRP channels as cellular sensors"),
    ("venkatachalam2007", "TRP channels", 'AUTH:"Venkatachalam" AND PUB_YEAR:2007'),
    ("nilius2011", "The transient receptor potential family of ion channels"),
    ("mccleskey2005", "Ion channels of nociception"),
    ("liman1999", "TRP2: a candidate transduction channel for mammalian pheromone sensory signaling"),
    # --------------------------------------------------- Cys-loop / pLGIC
    ("noda1982", "Primary structure of alpha-subunit precursor of Torpedo californica acetylcholine receptor deduced from cDNA sequence"),
    ("unwin2005", "Refined structure of the nicotinic acetylcholine receptor at 4A resolution"),
    ("brejc2001", "Crystal structure of an ACh-binding protein reveals the ligand-binding domain of nicotinic receptors"),
    ("hilf2008", "X-ray structure of a prokaryotic pentameric ligand-gated ion channel"),
    ("bocquet2009", "X-ray structure of a pentameric ligand-gated ion channel in an apparently open conformation"),
    ("hibbs2011", "Principles of activation and permeation in an anion-selective Cys-loop receptor"),
    ("miller2014", "Crystal structure of a human GABAA receptor"),
    ("du2015", "Glycine receptor mechanism elucidated by electron cryo-microscopy"),
    ("galzi1992", "Mutations in the channel domain of a neuronal nicotinic receptor convert ion selectivity from cationic to anionic"),
    ("thompson2010", "The structural basis of function in Cys-loop receptors"),
    # ------------------------------------------------------------ iGluR
    ("hollmann1994", "Cloned glutamate receptors"),
    ("sobolevsky2009", "X-ray structure, symmetry and mechanism of an AMPA-subtype glutamate receptor"),
    ("chen1999", "Functional characterization of a potassium-selective prokaryotic glutamate receptor"),
    ("traynelis2010", "Glutamate receptor ion channels: structure, regulation, and function"),
    ("karakas2014", "Crystal structure of a heterotetrameric NMDA receptor ion channel"),
    ("lee2014", "NMDA receptor structures reveal subunit arrangement and pore architecture"),
    ("sommer1991", "RNA editing in brain controls a determinant of ion flow in glutamate-gated channels"),
    ("lam1998", "Glutamate-receptor genes in plants"),
    ("benton2009", "Variant ionotropic glutamate receptors as chemosensory receptors in Drosophila"),
    ("kohda2000", "Mutation of a glutamate receptor motif reveals its role in gating and delta2 receptor channel properties"),
    # -------------------------------------------------------- P2X, ENaC
    ("kawate2009", "Crystal structure of the ATP-gated P2X(4) ion channel in the closed state"),
    ("north2002", "Molecular physiology of P2X receptors"),
    ("fountain2007", "An intracellular P2X receptor required for osmoregulation in Dictyostelium discoideum"),
    ("canessa1994", "Amiloride-sensitive epithelial Na+ channel is made of three homologous subunits"),
    ("jasti2007", "Structure of acid-sensing ion channel 1 at 1.9 A resolution and low pH"),
    ("chalfie1981", "Developmental genetics of the mechanosensory neurons of Caenorhabditis elegans"),
    ("kellenberger2002", "Epithelial sodium channel/degenerin family of ion channels: a variety of functions for a shared structure"),
    # ------------------------------------------------------ anion channels
    ("jentsch2002", "Molecular structure and physiological function of chloride channels"),
    ("dutzler2002", "X-ray structure of a ClC chloride channel at 3.0 A reveals the molecular basis of anion selectivity"),
    ("accardi2004", "Secondary active transport mediated by a prokaryotic homologue of ClC Cl- channels"),
    ("miller2006", "ClC chloride channels viewed through a transporter lens"),
    ("riordan1989", "Identification of the cystic fibrosis gene: cloning and characterization of complementary DNA"),
    ("zhang2017", "Conformational Change of the Extracellular Parts of the CFTR Protein during Channel Gating"),
    ("gadsby2006", "The ABC protein turned chloride channel whose failure causes cystic fibrosis"),
    ("yang2008", "TMEM16A confers receptor-activated calcium-dependent chloride conductance"),
    ("caputo2008", "TMEM16A, a membrane protein associated with calcium-dependent chloride channel activity"),
    ("schroeder2008", "Expression cloning of TMEM16A as a calcium-activated chloride channel subunit"),
    ("brunner2014", "X-ray structure of a calcium-activated TMEM16 lipid scramblase"),
    ("suzuki2010", "Calcium-dependent phospholipid scrambling by TMEM16F"),
    ("pedemonte2014", "Structure and function of TMEM16 proteins (anoctamins)"),
    ("sun2002", "The vitelliform macular dystrophy protein defines a new family of chloride channels"),
    ("dickson2014", "Structure and insights into the function of a Ca(2+)-activated Cl(-) channel"),
    ("voss2014", "Identification of LRRC8 heteromers as an essential component of the volume-regulated anion channel VRAC"),
    ("qiu2014", "SWELL1, a plasma membrane protein, is an essential component of volume-regulated anion channel"),
    ("deneka2018", "Structure of a volume-regulated anion channel of the LRRC8 family"),
    # ------------------------------------------------------- large pores
    ("maeda2009", "Structure of the connexin 26 gap junction channel at 3.5 A resolution"),
    ("panchin2000", "A ubiquitous family of putative gap junction molecules"),
    ("sosinsky2005", "Structural organization of gap junction channels"),
    ("ma2012", "Calcium homeostasis modulator 1 (CALHM1) is the pore-forming subunit of an ion channel that mediates extracellular Ca2+ regulation of neuronal excitability"),
    ("taruno2013", "CALHM1 ion channel mediates purinergic neurotransmission of sweet, bitter and umami tastes"),
    ("phelan2001", "Innexins: a family of invertebrate gap-junction proteins"),
    # ------------------------------------------------------ mechanosensation
    ("coste2010", "Piezo1 and Piezo2 are essential components of distinct mechanically activated cation channels"),
    ("ge2015", "Architecture of the mammalian mechanosensitive Piezo1 channel"),
    ("ranade2015", "Mechanically activated ion channels", 'AUTH:"Ranade" AND PUB_YEAR:2015'),
    ("sukharev1994", "A large-conductance mechanosensitive channel in E. coli encoded by mscL alone"),
    ("chang1998", "Structure of the MscL homolog from Mycobacterium tuberculosis: a gated mechanosensitive ion channel"),
    ("bass2002", "Crystal structure of Escherichia coli MscS, a voltage-modulated and mechanosensitive channel"),
    ("yuan2014", "OSCA1 mediates osmotic-stress-evoked Ca2+ increases vital for osmosensing in Arabidopsis"),
    ("murthy2018", "OSCA/TMEM63 are an evolutionarily conserved family of mechanically activated ion channels"),
    ("kawashima2011", "Mechanotransduction in mouse inner ear hair cells requires transmembrane channel-like genes"),
    ("pan2018", "TMC1 Forms the Pore of Mechanosensory Transduction Channels in Vertebrate Inner Ear Hair Cells"),
    ("walker2000", "A Drosophila mechanosensory transduction channel"),
    # -------------------------------------------------- intracellular channels
    ("furuichi1989", "Primary structure and functional expression of the inositol 1,4,5-trisphosphate-binding protein P400"),
    ("takeshima1989", "Primary structure and expression from complementary DNA of skeletal muscle ryanodine receptor"),
    ("desgeorges2016", "Structural Basis for Gating and Activation of RyR1"),
    ("fan2015", "Gating machinery of InsP3R channels revealed by electron cryomicroscopy"),
    ("berridge2016", "The Inositol Trisphosphate/Calcium Signaling Pathway in Health and Disease"),
    ("yazawa2007", "TRIC channels are essential for Ca2+ handling in intracellular stores"),
    ("baughman2011", "Integrative genomics identifies MCU as an essential component of the mitochondrial calcium uniporter"),
    ("destefani2011", "A forty-kilodalton protein of the inner membrane is the mitochondrial calcium uniporter"),
    ("colombini1979", "A candidate for the permeability pathway of the outer mitochondrial membrane"),
    ("cang2015", "TMEM175 Is an Organelle K(+) Channel Regulating Lysosomal Function"),
    ("lee2017tmem175", "The lysosomal potassium channel TMEM175 adopts a novel tetrameric architecture"),
    # ------------------------------------------------ proton, CRAC, misc
    ("ramsey2006", "A voltage-gated proton-selective channel lacking the pore domain"),
    ("sasaki2006", "A voltage sensor-domain protein is a voltage-gated proton channel"),
    ("murata2005", "Phosphoinositide phosphatase activity coupled to an intrinsic voltage sensor"),
    ("decoursey2013", "Voltage-gated proton channels: molecular biology, physiology, and pathophysiology of the H(V) family"),
    ("tu2018", "An evolutionarily conserved gene family encodes proton-selective ion channels"),
    ("feske2006", "A mutation in Orai1 causes immune deficiency by abrogating CRAC channel function"),
    ("vig2006", "CRACM1 is a plasma membrane protein essential for store-operated Ca2+ entry"),
    ("hou2012", "Crystal structure of the calcium release-activated calcium channel Orai"),
    ("prakriya2015", "Store-Operated Calcium Channels", 'AUTH:"Prakriya" AND PUB_YEAR:2015'),
    ("littleton2000", "Ion channels and synaptic organization: analysis of the Drosophila genome"),
    # ------------------------------------------------------- evolution
    ("anderson2001", "Phylogeny of ion channels: clues to structure and function"),
    ("jegla2009", "Evolution of the human ion channel set"),
    ("liebeskind2011", "Evolution of sodium channels predates the origin of nervous systems in animals"),
    ("liebeskind2015", "Convergence of ion channel genome content in early animal evolution"),
    ("moran2015", "Evolution of voltage-gated ion channels at the emergence of Metazoa"),
    ("hedrich2012", "Ion channels in plants"),
    ("martinac2008", "Ion channels in microbes"),
    ("derst1998", "Evolutionary link between prokaryotic and eukaryotic K+ channels"),
    # ----------------------------------------------- disease, pharmacology
    ("ashcroft2006", "From molecule to malady"),
    ("kullmann2010", "Neurological channelopathies"),
    ("abriel2015", "Ion channel macromolecular complexes in cardiomyocytes: roles in sudden cardiac death"),
    ("bagal2013", "Ion channels as therapeutic targets: a drug discovery perspective"),
    ("santos2017", "A comprehensive map of molecular drug targets"),
    ("overington2006", "How many drug targets are there?"),
    ("wulff2009", "Voltage-gated potassium channels as therapeutic targets"),
    ("alexander2023", "THE CONCISE GUIDE TO PHARMACOLOGY 2023/24: Ion channels"),
    # ------------------------------------------------ databases and tools
    ("paysanlafosse2023", "InterPro in 2022"),
    ("mistry2021", "Pfam: The protein families database in 2021"),
    ("uniprot2023", "UniProt: the Universal Protein Knowledgebase in 2023"),
    ("katoh2013", "MAFFT multiple sequence alignment software version 7: improvements in performance and usability"),
    ("minh2020", "IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era"),
    ("kalyaanamoorthy2017", "ModelFinder: fast model selection for accurate phylogenetic estimates"),
    ("capellagutierrez2009", "trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses"),
    ("vankempen2024", "Fast and accurate protein structure search with Foldseek"),
    ("jumper2021", "Highly accurate protein structure prediction with AlphaFold"),
    ("varadi2022", "AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models"),
    ("eddy2011", "Accelerated Profile HMM Searches"),
]
