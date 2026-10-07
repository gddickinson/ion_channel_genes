# S12 — Structures

AlphaFold DB coverage of the final census, experimental references and a model check, and S8a's fold network re-read with experimental units, TM-region units and a dense set of representatives. Rules: **D54**, fixed before any structure was fetched; S8a's verdicts (D48) stay primary.

## 1. AlphaFold DB and PDB coverage of the census (D54 (1))

Frame: S15's final census — **7,061 proteome members** in 75 census families (plus 336 genome loci, which have no UniProt accession and so no AlphaFold DB entry).

* **AlphaFold DB model of the exact accession: 5,655 (80.1 %)**; model sequence identical to the census sequence 5,512 (78.1 %).
* **Usable** (exact, global pLDDT ≥ 70): **4,103 (58.1 %)**; median global pLDDT 76.2.
* **Any experimental structure** (UniProt PDB cross-reference): **390 (5.5 %)**.

**Families with no usable model**: ryr (48 members, 0 models, median pLDDT –), clcc1 (22 members, 18 models, median pLDDT 56.4), tmem109 (12 members, 11 models, median pLDDT 57.5).

**Lowest usable share** (families with ≥ 1 usable model):

| family | superfamily | members | model | usable | frac_usable | median_plddt |
|---|---|---|---|---|---|---|
| cav | ploop | 197 | 47 | 5 | 0.0254 | 61.9 |
| kv_kcnq | ploop | 89 | 72 | 3 | 0.0337 | 58.4 |
| itpr | ca_release | 86 | 5 | 5 | 0.0581 | 74.6 |
| nav | ploop | 134 | 42 | 13 | 0.097 | 69.1 |
| piezo | piezo | 44 | 7 | 6 | 0.1364 | 72.1 |
| mitok | mitok | 14 | 14 | 2 | 0.1429 | 68.0 |
| trpn | ploop | 13 | 2 | 2 | 0.1538 | 75.3 |
| hv1 | hv | 23 | 19 | 4 | 0.1739 | 58.6 |

**Why members lack a model.** Of the 1,406 members with no AlphaFold DB entry, only 128 are longer than 2,700 residues (AlphaFold DB's length limit); the rest are simply not in the database. Lowest model share (≥ 20 members):

| family | members | model | frac_model | no_model_over_2700aa |
|---|---|---|---|---|
| ryr | 48 | 0 | 0.0 | 48 |
| itpr | 86 | 5 | 0.0581 | 47 |
| piezo | 44 | 7 | 0.1591 | 15 |
| cav | 197 | 47 | 0.2386 | 6 |
| nav | 134 | 42 | 0.3134 | 1 |

Per family: `coverage_families.tsv`; per member: `coverage_members.tsv`.

## 2. The experimental reference and the model check (D54 (2)–(3))

**56 of 75 census families have an experimental reference** among their catalogue exemplars (48 Electron Microscopy, 5 X-ray diffraction, 3 Solution NMR).
Chains are cut by SIFTS per-residue UniProt numbering; residues SIFTS maps to an isoform are renumbered to the canonical sequence by global alignment (addendum).

**Model check**: 50 of 52 comparable AFDB units are consistent with experiment (TM-score ≥ 0.5, normalised by the experimental unit), 41 at ≥ 0.8; median 0.917. 3 families' reference is a different exemplar from S8a's node, 1 has no model (RyR).

Inconsistent:

| node | pdb | chain | method | observed | tm_by_exp | rmsd |
|---|---|---|---|---|---|---|
| hv1 | 5oqk | A | Solution NMR | 118 | 0.4975 | 3.54 |
| viroporin | 2n29 | A | Solution NMR | 54 | 0.3129 | 2.12 |

**RyR** (no AFDB model) enters through Hs_RYR2 7u9x chain A (Electron Microscopy, 2.58 Å), module residues 4769–4870.

No experimental reference (19): deg_invertebrate, gphr, ampa, iglur_prok, kainate, innexin, clcc1, mitok, orai, otop, p2x_nonmetazoan, catsper, kv_modifier, trpn, tmem109, ano_channel, ano_scramblase, tmc, tric.

## 3. The literature edges read four ways (D54 (4)–(6))

S8a's AFDB reading is the primary verdict (D48); the other three are sensitivity readings. Each cell: median TM-score, rank b-for-a / a-for-b, verdict.

| edge | s8a_afdb | experimental | tm_region | dense |
|---|---|---|---|---|
| tmem16_like – tmem16_like | 0.4839 (1/1) **not distinguished** | 0.4886 (1/1) **not distinguished** | 0.5334 (1/1) **supported** | 0.5114 (1/1) **supported** |
| iglur – ploop | 0.6364 (1/1) **supported** | 0.6147 (1/1) **supported** | 0.6364 (1/1) **supported** | 0.6528 (1/1) **supported** |
| innexin_like – connexin | 0.5648 (1/1) **supported** | 0.5592 (1/1) **supported** | 0.6029 (1/1) **supported** | 0.5639 (1/1) **supported** |
| ca_release – ploop | 0.556 (2/2) **not distinguished** | 0.5888 (1/2) **not distinguished** | 0.556 (2/2) **not distinguished** | 0.6248 (1/2) **not distinguished** |
| hv – ploop:VSD | 0.5922 (1/1) **supported** | 0.4148 (1/1) **not distinguished** | 0.5922 (1/1) **supported** | unmeasured |

Experimental variant: 57 nodes on an experimental unit, 19 on S8a's AFDB unit.

## 4. The dense set: recovery and detectability (D54 (6))

**268 representatives** of 72 census families (one per family × S4 group; 9 dropped, reasons in `dense_reps.tsv`), all-vs-all TM-align and Foldseek.

**Superfamily recovery** (superfamilies with ≥ 2 measured families, 161 representatives): the closest structure outside a representative's own family is in its own superfamily for **156 (96.9 %) by TM-align** and **156 (96.9 %) by Foldseek E-value**.

| superfamily | families | reps | tm_recovered | foldseek_recovered |
|---|---|---|---|---|
| cysloop | 7 | 14 | 14 | 14 |
| deg_enac | 3 | 7 | 7 | 7 |
| iglur | 6 | 12 | 11 | 11 |
| innexin_like | 3 | 5 | 5 | 5 |
| msc | 2 | 4 | 0 | 0 |
| p2x | 2 | 10 | 10 | 10 |
| ploop | 23 | 85 | 85 | 85 |
| tmem16_like | 4 | 24 | 24 | 24 |

Not recovered by TM-align — the family the closest partner belongs to:

| superfamily | closest_family | reps |
|---|---|---|
| msc | calhm | 2 |
| iglur | kca_sk | 1 |
| msc | viroporin | 1 |
| msc | orai | 1 |

**Foldseek detectability** between superfamilies (share of representative pairs with E ≤ 1e-3): **4 of 435 superfamily pairs** have any detected pair:

| a | b | rep_pairs | foldseek_detected | frac_detected | median_tm |
|---|---|---|---|---|---|
| p2x | pac | 10 | 10 | 1.0 | 0.3441 |
| deg_enac | pac | 7 | 4 | 0.5714 | 0.4267 |
| connexin | innexin_like | 10 | 2 | 0.2 | 0.5565 |
| iglur | ploop | 1020 | 9 | 0.0088 | 0.6475 |

Figure: `results/structures/figures/structures.png`.
