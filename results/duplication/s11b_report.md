# S11b — repeat-duplication order and outgroup-free roots

Rendered by `scripts/s11b_report.py` from the tables beside it (D13). Rules: D53. Bulk under `<data root>/trees/s11b/`.

## 1. The order of the internal repeat duplications

**Rules (D53 (1)–(4)), fixed before any input was built.** Whole S1–S6 repeats (14 spans in profile states, from the S6 references' TRANSMEM helices) of one chain per family × S4 group — the most central chain whose every repeat covers ≥ 90 % of its profile states — aligned together (L-INS-i, trimAl `-gt 0.5`); the pore module of the same chains is the sensitivity unit. Three pairings of the four 4×6TM repeat classes are each an ML search constrained by their bipartition (TPC's two repeats free) and compared with the unconstrained ML tree by the AU test (rejected iff p < 0.05).

**Representatives**: 25 chains nav 5, cav 7, nalcn 4, tpc 9; cells with no eligible chain: 1.

| unit | tips | aligned_cols | cols | informative | gap_frac |
|---|---|---|---|---|---|
| repeat | 82 | 755 | 231 | 229 | 0.055 |
| pore | 82 | 323 | 106 | 106 | 0.0831 |

**AU test** (ΔlnL from the best tree):

| unit | tree | model | logL | deltaL | p_SH | p_AU | satisfies | rejected |
|---|---|---|---|---|---|---|---|---|
| repeat | ML | LG+F+R5 | -32310.54322 | 4.0318 | 0.542 | 0.419 |  | False |
| repeat | H13 | LG+F+R5 | -32306.51145 | 0.0 | 1.0 | 0.712 | True | False |
| repeat | H12 | LG+F+R5 | -32327.77949 | 21.268 | 0.171 | 0.0993 | True | False |
| repeat | H14 | LG+F+R5 | -32329.08792 | 22.576 | 0.148 | 0.0619 | True | False |
| pore | ML | Q.pfam+G4 | -13640.24866 | 0.0 | 1.0 | 0.683 |  | False |
| pore | H13 | Q.pfam+G4 | -13643.47308 | 3.2244 | 0.607 | 0.471 | True | False |
| pore | H12 | Q.pfam+G4 | -13644.43633 | 4.1877 | 0.558 | 0.354 | True | False |
| pore | H14 | Q.pfam+G4 | -13646.31294 | 6.0643 | 0.461 | 0.258 | True | False |

* **repeat** unit: **unresolved:H13,H12,H14**
* **pore** unit: **unresolved:H13,H12,H14**

Pairings: **H13** {I,III}|{II,IV} (a tandem duplication of a two-repeat, TPC-like ancestor — stated in advance as the only pairing that fits it), **H12** {I,II}|{III,IV}, **H14** {I,IV}|{II,III}.

**Repeat classes in the unconstrained ML tree** (descriptive): is each class one clade across Nav, Cav and NALCN, how many other tips stand in the way, and the class make-up of the smallest edge side that holds the whole class plus anything else.

| unit | class | tips | one_clade | ufboot | intruders | intruder_classes | nearest_side |
|---|---|---|---|---|---|---|---|
| repeat | I | 16 | False |  | 11 | III | III:11 |
| repeat | II | 16 | False |  | 25 | IV,tII | IV:16;tII:9 |
| repeat | III | 16 | False |  | 16 | I | I:16 |
| repeat | IV | 16 | True | 76.0 | 0 |  | II:12 |
| repeat | tI | 9 | True | 100.0 | 0 |  | I:16;III:16 |
| repeat | tII | 9 | False |  | 4 | II | II:4 |
| pore | I | 16 | False |  | 50 | II,III,tI,tII | II:16;III:16;tI:9;tII:9 |
| pore | II | 16 | False |  | 20 | I,IV | I:4;IV:16 |
| pore | III | 16 | True | 97.0 | 0 |  | I:12;tI:9;tII:9 |
| pore | IV | 16 | True | 90.0 | 0 |  | I:4;II:4 |
| pore | tI | 9 | True | 97.0 | 0 |  | I:12 |
| pore | tII | 9 | True | 75.0 | 0 |  | I:12;tI:9 |

## 2. Outgroup-free roots (non-reversible model)

**Rules (D53 (5)–(7)).** S7b's ingroup alignment for each family (outgroup rows removed, D41 mask unchanged), NQ.pfam with S7b's rate-heterogeneity terms, 1,000 UFBoot (→ rootstrap) and an AU test over every root branch. A root is **resolved** iff rootstrap ≥ 95 and the non-reversible model fits better than Q.pfam on the same tree (ΔlnL > 0). It **agrees** with S11a iff its split is one of S11a's optimal reconciliation roots on S7b's tree. Controls: families whose declared outgroup root S11a accepted and also found optimal.

| family | role | n_tips | cols | s7b_model | nq_model |
|---|---|---|---|---|---|
| cng | root | 313 | 768 | Q.pfam+F+R10 | NQ.pfam+R10 |
| k2p | root | 380 | 302 | Q.pfam+R9 | NQ.pfam+R9 |
| kv_kcnq | root | 97 | 667 | JTT+I+R5 | NQ.pfam+I+R5 |
| kv_shaker | root | 312 | 458 | Q.pfam+I+R9 | NQ.pfam+I+R9 |
| kca_slo | root | 106 | 1033 | LG+F+I+R7 | NQ.pfam+I+R7 |
| kv_eag | root | 158 | 909 | JTT+F+I+R7 | NQ.pfam+I+R7 |
| enac | control | 41 | 654 | JTT+R4 | NQ.pfam+R4 |
| glyr | control | 65 | 454 | JTT+R4 | NQ.pfam+R4 |
| ht3 | control | 52 | 457 | JTT+R5 | NQ.pfam+R5 |
| iglur_nonvertebrate | control | 38 | 916 | Q.pfam+R5 | NQ.pfam+R5 |
| nalcn | control | 24 | 1738 | LG+F+R4 | NQ.pfam+R4 |

**Controls** (declared root known):

| family | n_tips | dlnl_nonrev | root_rootstrap | root_split | ml_root_tip | au_set | branches_tested | au_set_frac | max_rootstrap | control_recovered | declared_in_au_set |
|---|---|---|---|---|---|---|---|---|---|---|---|
| enac | 41 | -2.81 | 30.2 | 1|40 | A0AAJ7X142__Petmar | 46 | 79 | 0.58 | 30.2 | False | False |
| glyr | 65 | 8.88 | 0.0 | 1|64 | GCA_977017645.1_OZ360232.1_52680552-52824349-__Tormar | 78 | 127 | 0.61 | 43.2 | False | True |
| ht3 | 52 | 5.86 | 0.0 | 1|51 | A0A674P1E2__Takrub | 68 | 101 | 0.67 | 30.2 | False | False |
| iglur_nonvertebrate | 38 | 9.01 | 0.0 | 1|37 | O04660__Aratha | 46 | 73 | 0.63 | 35.6 | False | True |
| nalcn | 24 | 22.31 | 32.0 | 7|17 |  | 32 | 45 | 0.71 | 32.0 | False | True |

**0 / 5 controls recover the declared root as the ML root; the declared root is inside the AU confidence set in 3 / 5.**

**S7d's six:**

| family | n_tips | dlnl_nonrev | root_rootstrap | root_split | small_side | au_set | au_set_frac | max_rootstrap | resolved | s11a_split | agrees_s11a | recon_optimum_nq_topology | s11a_root_in_au_set |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| cng | 313 | 144.6 | 0.0 | 1|312 | Monosiga brevicollis | 409 | 0.66 | 72.8 | False | 142|171 | False | False | True |
| k2p | 380 | 128.2 | 0.4 | 1|379 | Daphnia pulex | 456 | 0.6 | 38.1 | False | 98|282 | False | False | False |
| kv_kcnq | 97 | 3.41 | 43.5 | 1|96 | Hydra vulgaris | 132 | 0.69 | 43.5 | False | 29|68 | False | False | True |
| kv_shaker | 311 | 49.91 | 0.0 | 1|310 | Homo sapiens | 320 | 0.52 | 60.6 | False |  | False | False | True |
| kca_slo | 106 | 36.2 | 21.7 | 36|70 | Amphimedon queenslandica; Anolis carolinensis; Aplysia californica; Branchiostoma floridae … | 148 | 0.71 | 21.7 | False | 37|69 | False | False | True |
| kv_eag | 158 | -11.59 | 95.0 | 1|157 | Cornu aspersum | 161 | 0.51 | 95.0 | False | 7|151 | False | False | True |

**0 / 6 roots resolved under D53 (6); 0 / 6 agree with S11a's reconciliation root**; S11a's root lies inside the AU confidence set in 5 / 6. Disagreements are reported, not resolved.

**What the roots look like** (descriptive): in 9 of 11 families the non-reversible ML root falls on a single terminal branch, and the AU confidence set holds 51%–71% of all branches in every family — the non-reversible signal barely discriminates among roots at these depths. The controls, whose roots are known, measure it directly.

