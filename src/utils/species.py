"""Species table — taxon ids, common names and Ensembl slugs.

Ion channels are not a vertebrate story. Potassium channels are in every
domain of life, the mechanosensitive MscS family is in bacteria, plants and
fungi and absent from animals, the OSCA channels were found in plants first,
P2X receptors are in amoebae and gone from *Drosophila*, and TRPN was lost
in mammals. A species panel drawn from the usual model organisms would make
the census a description of vertebrate biology with a few footnotes.

So the panel spans the tree: bacteria and archaea (where the structural
prototypes live), the non-metazoan eukaryotes that decide the origin
questions, the early-branching animals, the invertebrate models, and the
vertebrate series. `PANELS` groups them so a task can ask for the slice it
needs — a tier-2 rooting run wants the prokaryotes, a duplication-timing run
wants the vertebrate series, and neither wants the other.

Every `taxon_id` here is a claim, checked against the UniProt taxonomy API
by `scripts/s0_catalogue_verify.py`. Ensembl slugs are best-effort: most of
these species are not in Ensembl Vertebrates at all, and the client falls
back to a name query.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SpeciesInfo:
    scientific: str
    common: str
    taxon_id: int
    ensembl_slug: str
    group: str = ""       # the panel this species anchors
    note: str = ""


_SPECIES: list[SpeciesInfo] = [
    # --- bacteria & archaea: the structural prototypes -------------------
    SpeciesInfo("Escherichia coli", "E. coli", 83333, "escherichia_coli",
                "prokaryote", "MscL, MscS, ClC-ec1, KefC"),
    SpeciesInfo("Bacillus subtilis", "B. subtilis", 1423, "bacillus_subtilis",
                "prokaryote", "YugO, MscL"),
    SpeciesInfo("Streptomyces lividans", "S. lividans", 1916, "",
                "prokaryote", "KcsA — the first K+ channel structure"),
    SpeciesInfo("Synechocystis sp. PCC 6803", "Synechocystis", 1148, "",
                "prokaryote", "GluR0, the K+-selective glutamate receptor"),
    SpeciesInfo("Gloeobacter violaceus", "Gloeobacter", 33072, "",
                "prokaryote", "GLIC, the proton-gated pentameric channel"),
    SpeciesInfo("Aliarcobacter butzleri", "Aliarcobacter", 367737, "",
                "prokaryote", "NavAb, the prototype voltage-gated Na+ channel; "
                "the genus was renamed from Arcobacter"),
    SpeciesInfo("Methanothermobacter thermautotrophicus", "M. thermautotrophicus",
                187420, "", "prokaryote", "MthK, the Ca2+-gated K+ channel"),
    SpeciesInfo("Thermus thermophilus", "T. thermophilus", 274, "",
                "prokaryote", "NaK and thermostable channel structures"),
    # --- non-metazoan eukaryotes: where the origins are decided ----------
    SpeciesInfo("Saccharomyces cerevisiae", "Baker's yeast", 4932,
                "saccharomyces_cerevisiae", "fungi", "TRPY1/Yvc1, TOK1"),
    SpeciesInfo("Schizosaccharomyces pombe", "Fission yeast", 4896,
                "schizosaccharomyces_pombe", "fungi", ""),
    SpeciesInfo("Arabidopsis thaliana", "Thale cress", 3702,
                "arabidopsis_thaliana", "plant", "OSCA, MSL, GLR, TPC1, SLAC1"),
    SpeciesInfo("Oryza sativa", "Rice", 4530, "oryza_sativa", "plant", ""),
    SpeciesInfo("Physcomitrium patens", "Moss", 3218, "", "plant", ""),
    SpeciesInfo("Chlamydomonas reinhardtii", "Chlamydomonas", 3055,
                "chlamydomonas_reinhardtii", "algae",
                "channelrhodopsins; the light-gated outgroup"),
    SpeciesInfo("Dictyostelium discoideum", "Slime mould", 44689,
                "dictyostelium_discoideum", "amoebozoa",
                "P2X receptors, IP3R-like iplA"),
    SpeciesInfo("Paramecium tetraurelia", "Paramecium", 5888, "", "ciliate",
                "the ciliate voltage-gated Ca2+ channels"),
    SpeciesInfo("Tetrahymena thermophila", "Tetrahymena", 5911, "", "ciliate", ""),
    SpeciesInfo("Trypanosoma brucei", "Trypanosome", 5691, "", "excavate", ""),
    SpeciesInfo("Plasmodium falciparum", "Malaria parasite", 5833, "",
                "apicomplexa", ""),
    SpeciesInfo("Monosiga brevicollis", "Choanoflagellate", 81824, "",
                "holozoa", "the animals' closest unicellular relatives"),
    SpeciesInfo("Capsaspora owczarzaki", "Capsaspora", 595528, "", "holozoa", ""),
    # --- early-branching animals ----------------------------------------
    SpeciesInfo("Amphimedon queenslandica", "Sponge", 400682, "", "basal_metazoan",
                "a sponge with channel genes and no neurons"),
    SpeciesInfo("Trichoplax adhaerens", "Placozoan", 10228, "", "basal_metazoan", ""),
    SpeciesInfo("Nematostella vectensis", "Starlet sea anemone", 45351,
                "", "cnidarian", "a nearly complete channel repertoire"),
    SpeciesInfo("Hydra vulgaris", "Hydra", 6087, "", "cnidarian", ""),
    # --- protostome invertebrates ---------------------------------------
    SpeciesInfo("Caenorhabditis elegans", "C. elegans", 6239,
                "caenorhabditis_elegans", "invertebrate",
                "the degenerins; ~90 K+ channel genes"),
    SpeciesInfo("Drosophila melanogaster", "Fruit fly", 7227,
                "drosophila_melanogaster", "invertebrate",
                "Shaker, trp, nompC, ppk, the ionotropic receptors"),
    SpeciesInfo("Daphnia pulex", "Water flea", 6669, "", "invertebrate", ""),
    SpeciesInfo("Lottia gigantea", "Owl limpet", 225164, "", "invertebrate", ""),
    SpeciesInfo("Aplysia californica", "Sea hare", 6500, "", "invertebrate",
                "AChBP; the classic electrophysiology preparation"),
    SpeciesInfo("Lymnaea stagnalis", "Pond snail", 6523, "", "invertebrate",
                "the acetylcholine-binding protein (hazard H2)"),
    SpeciesInfo("Cornu aspersum", "Garden snail", 6535, "", "invertebrate",
                "FaNaC, the peptide-gated channel; formerly Helix aspersa"),
    # --- deuterostomes to vertebrates -----------------------------------
    SpeciesInfo("Strongylocentrotus purpuratus", "Purple sea urchin", 7668,
                "", "deuterostome", "CatSper and sperm channels"),
    SpeciesInfo("Branchiostoma floridae", "Amphioxus", 7739, "", "deuterostome",
                "pre-2R chordate outgroup"),
    SpeciesInfo("Ciona intestinalis", "Sea squirt", 7719, "ciona_intestinalis",
                "deuterostome", "Ci-VSP and Ci-Hv1 (hazard H9)"),
    SpeciesInfo("Petromyzon marinus", "Sea lamprey", 7757, "petromyzon_marinus",
                "vertebrate", "cyclostome; brackets the 2R duplications"),
    SpeciesInfo("Callorhinchus milii", "Elephant shark", 7868,
                "callorhinchus_milii", "vertebrate", "slow-evolving chondrichthyan"),
    SpeciesInfo("Torpedo marmorata", "Marbled electric ray", 7788, "",
                "vertebrate", "the nicotinic receptor's source tissue"),
    SpeciesInfo("Danio rerio", "Zebrafish", 7955, "danio_rerio", "vertebrate",
                "teleost 3R duplicates"),
    SpeciesInfo("Takifugu rubripes", "Fugu", 31033, "takifugu_rubripes",
                "vertebrate", "compact teleost genome"),
    SpeciesInfo("Latimeria chalumnae", "Coelacanth", 7897, "latimeria_chalumnae",
                "vertebrate", "pre-tetrapod sarcopterygian"),
    SpeciesInfo("Xenopus tropicalis", "Western clawed frog", 8364,
                "xenopus_tropicalis", "vertebrate",
                "the expression system most channel physiology was done in"),
    SpeciesInfo("Anolis carolinensis", "Green anole", 28377, "anolis_carolinensis",
                "vertebrate", ""),
    SpeciesInfo("Gallus gallus", "Chicken", 9031, "gallus_gallus", "vertebrate", ""),
    SpeciesInfo("Ornithorhynchus anatinus", "Platypus", 9258,
                "ornithorhynchus_anatinus", "vertebrate", "monotreme"),
    SpeciesInfo("Monodelphis domestica", "Opossum", 13616, "monodelphis_domestica",
                "vertebrate", "marsupial"),
    SpeciesInfo("Mus musculus", "Mouse", 10090, "mus_musculus", "vertebrate", ""),
    SpeciesInfo("Rattus norvegicus", "Rat", 10116, "rattus_norvegicus",
                "vertebrate", ""),
    SpeciesInfo("Homo sapiens", "Human", 9606, "homo_sapiens", "vertebrate",
                "the census reference"),
    # --- viruses (the viroporin division) -------------------------------
    SpeciesInfo("Influenza A virus", "Influenza A", 11320, "", "virus",
                "M2, the amantadine target"),
    SpeciesInfo("Human immunodeficiency virus type 1", "HIV-1", 11676, "",
                "virus", "Vpu"),
    SpeciesInfo("Severe acute respiratory syndrome coronavirus 2", "SARS-CoV-2",
                2697049, "", "virus", "E and ORF3a"),
]


SPECIES_LOOKUP: dict[str, SpeciesInfo] = {}
for _sp in _SPECIES:
    SPECIES_LOOKUP[_sp.scientific.lower()] = _sp
    SPECIES_LOOKUP[_sp.common.lower()] = _sp
    if _sp.ensembl_slug:
        SPECIES_LOOKUP[_sp.ensembl_slug] = _sp


#: Named slices of the table. A task asks for the panel its question needs.
PANELS: dict[str, list[str]] = {
    "human": ["Homo sapiens"],
    "vertebrate": [s.scientific for s in _SPECIES if s.group == "vertebrate"],
    "metazoan": [s.scientific for s in _SPECIES
                 if s.group in ("vertebrate", "deuterostome", "invertebrate",
                                "cnidarian", "basal_metazoan")],
    "eukaryote": [s.scientific for s in _SPECIES if s.group not in
                  ("prokaryote", "virus")],
    "prokaryote": [s.scientific for s in _SPECIES if s.group == "prokaryote"],
    "outgroup": [s.scientific for s in _SPECIES
                 if s.group in ("prokaryote", "fungi", "plant", "algae",
                                "amoebozoa", "ciliate", "excavate",
                                "apicomplexa", "holozoa")],
    "virus": [s.scientific for s in _SPECIES if s.group == "virus"],
    "tree_of_life": [
        "Escherichia coli", "Synechocystis sp. PCC 6803",
        "Saccharomyces cerevisiae", "Arabidopsis thaliana",
        "Dictyostelium discoideum", "Monosiga brevicollis",
        "Amphimedon queenslandica", "Nematostella vectensis",
        "Caenorhabditis elegans", "Drosophila melanogaster",
        "Branchiostoma floridae", "Danio rerio", "Xenopus tropicalis",
        "Gallus gallus", "Mus musculus", "Homo sapiens",
    ],
}

#: Default for an unscoped search: broad enough to see a family's range,
#: small enough to finish.
DEFAULT_SPECIES_PANEL: list[str] = PANELS["tree_of_life"]


def resolve_species(name: str) -> SpeciesInfo | None:
    if not name:
        return None
    return SPECIES_LOOKUP.get(name.strip().lower())


def ensembl_species_slug(name: str) -> str:
    info = resolve_species(name)
    if info and info.ensembl_slug:
        return info.ensembl_slug
    return name.strip().lower().replace(" ", "_")


def panel(name: str) -> list[str]:
    """A named species panel, or a single species name passed straight through."""
    if name in PANELS:
        return list(PANELS[name])
    info = resolve_species(name)
    return [info.scientific] if info else ([name] if name else [])


def all_species() -> list[SpeciesInfo]:
    return list(_SPECIES)


def groups() -> dict[str, int]:
    out: dict[str, int] = {}
    for s in _SPECIES:
        out[s.group] = out.get(s.group, 0) + 1
    return out
