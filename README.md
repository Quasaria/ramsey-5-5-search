# ramsey-5-5-search

Heuristic search for graphs on 43 vertices avoiding monochromatic K₅ —
an experimental attack on the Ramsey number bound **R(5,5) ≥ 44**.

## TL;DR

Using simulated annealing seeded from exhaustively-enumerated Cayley
graphs on ℤ₄₃ (the **circulant family**), we descend from random
graphs (K₅_total ≈ 2700) all the way down to **K₅_total = 2** on 43
vertices — a hard wall reproducible across 13 structurally distinct
witnesses. Exhaustive Hamming-3 analysis (over critical edges and an
extended neighborhood) confirms that this K₅ = 2 plateau is a strict
local minimum: no combination of 1, 2, or 3 edge flips reduces it.

We also compare with a structurally different family — Cayley graphs
on **ℤ₄₂ + a central vertex** — and observe a much higher plateau
(K₅_total = 840 minimum after full enumeration of ~2M configurations,
falling to ~105 after SA), illustrating concretely why the primality
of 43 matters for Ramsey-style problems.

**No witness of R(5,5) ≥ 44 was found.** This is consistent with
the long-standing community suspicion that R(5,5) = 43 exactly, but
this work does **not** prove that — only SAT-based exhaustive methods
could.

## Status of this work

This is a personal experimental project, not a research contribution.
R(5,5) has been an open problem since the late 1980s and is currently
known to lie in [43, 48]; existence of K₅-free 2-colorings on 42
vertices has been established by McKay–Radziszowski. The methods used
here (simulated annealing, tabu search, iterated local search,
warmstarting from algebraic graphs) are standard in the literature.

The contributions of this repository are:

1. A **clean, reproducible Python toolkit** for K₅-counting and
   heuristic search on 43-vertex graphs, with bitset-based exact
   counting at ~5000–8000 iterations/sec.
2. A **systematic comparison** of two Cayley families (ℤ₄₃ vs
   ℤ₄₂ + center), with concrete numbers showing the gap.
3. An **exhaustive Hamming-1/2/3 analysis** of the K₅ = 2 plateau,
   verifying it is a strict local minimum across 13 distinct
   witnesses.
4. A **multi-process parallel warmstart pipeline** for running 21
   independent simulated-annealing trajectories from the best
   circulants.

None of these are mathematical novelties — but the toolkit and
numbers are honest, reproducible, and might be useful to anyone
attacking Ramsey-style combinatorial problems on a personal machine.

## Background: the R(5,5) problem

The Ramsey number R(s, t) is the smallest n such that any 2-coloring
of the edges of K_n contains either a red K_s or a blue K_t. R(5, 5)
is famously hard:

- **Known bounds**: 43 ≤ R(5,5) ≤ 48 (Angeltveit–McKay, 2017).
- **Witness for R(5,5) ≥ 43**: explicit graphs on 42 vertices with
  no monochromatic K₅ are known.
- **Witness for R(5,5) ≥ 44**: would require a graph on 43 vertices
  with no monochromatic K₅. **None has ever been found.**

A graph G on n vertices avoids a monochromatic K₅ iff
`K₅(G) + K₅(Ḡ) = 0`. We call this quantity `K₅_total`. The search
problem is: minimize `K₅_total` over all graphs on 43 vertices, and
ideally reach 0.

## Method

### Cost function

```
cost(G) = w₅ · K₅_total(G) + w₄ · K₄_total(G)
```

where `K₄_total` provides a smoother secondary signal when `K₅_total`
saturates. Defaults: `w₅ = 100`, `w₄ = 5` (full search) or `10`
(warmstart).

### Bitset representation

The adjacency matrix is stored as 43 Python integers (bitmasks). All
counts use bit-tricks: K₅-enumeration runs as a 5-level nested loop
over `bits_b & higher_mask(a)` etc., with `int.bit_count()` for the
innermost summation. Performance: ≈ 5000–8000 SA iterations per
second on a single core (Python 3.13, Windows 11, Intel i7).

### Simulated annealing (v5)

Located in `ramsey_55_search.py`. Implements:

- **Move types**: 1-flip (default), 2-flip with lookahead-K sampling
  (default K=4), vertex-block moves (flip 3–5 edges incident to one
  vertex).
- **Tabu list** of length 30 over recently flipped edges.
- **Iterated Local Search**: on stagnation, perturb from `best`
  (not `current`) by flipping `perturb_size` random edges, then
  resume.
- **Periodic integrity checks** (every 500k iterations): recount
  from scratch and assert against incremental counters.
- **Robust saves** with atomic-rename + retry, designed to survive
  Windows antivirus interference.
- **Resume support** (`--resume` or `--start_from FILE`).

### Circulant enumeration (warmstart source)

`ramsey_55_circulant.py` exhaustively enumerates the 2²¹ ≈ 2M
circulant graphs Circ(43, S) where S ⊂ {1, ..., 21} is symmetric.
Multiprocessing over chunks; top-K saved.

For circulants, the **vertex-0 trick** applies (orbits under ℤ₄₃
have trivial stabilizers because 43 is prime and 5 ∤ 43, so every
K₅ has orbit size 43):

```
K₅(Circ(43, S)) = 43 · K₄(N(0)) / 5
```

This is exact, and lets us count K₅ in O(|S|³) instead of O(43³).

### ℤ₄₂ + center family (alternative)

`ramsey_55_z42_center.py` enumerates a structurally different family:
Cayley graphs on ℤ₄₂ (42 ring vertices) plus one central vertex
connected to all or none. The full enumeration covers 2²² ≈ 4M
configurations.

### Parallel warmstart

`run_warmstarts.ps1` (PowerShell) launches N independent SA jobs in
parallel, each seeded from a different top circulant. Used with
`-NumStarts 21 -Parallel 8 -Hours 4` to attack the top-21 circulants
in three waves on an 8-core machine (84 CPU-hours total).

### Hamming-1/2/3 analysis

`ramsey_55_analyze.py` enumerates the explicit K₅ instances in a
saved graph, identifies "critical edges" (edges belonging to at
least one K₅), and computes the exact `cost` delta for:

- All single edge flips (`--pairs` enables also paired flips)
- All pairs of critical-edge flips (`--pairs`)
- All triples of edge flips (`--triples`, restricted to critical
  edges; `--triples_extended` widens to all edges incident to
  any vertex appearing in some K₅)

## Numerical results

### Trajectory of descent (ℤ₄₃ circulants)

| Stage | K₅_total | How |
|---|---|---|
| Random init | ~2700 | Random graph, p = 0.5 |
| Full-space SA, no warmstart | 151–187 | Hits wall around K₅ ≈ 170 |
| 2-flip lookahead + warmstart from best circulant | 9 | After ~340k iterations |
| v5 with tabu + ILS + warmstart | **2** | After ~3k iterations, then stuck |

The best circulant has S = {1, 4, 5, 6, 7, 8, 9, 12, 14, 17}, |S| =
10, and `K₅_total(Circ) = 43`, `K₄_total = 2623`. From this seed, SA
descends rapidly: 43 → 9 → 2.

### Top-21 circulants

Exhaustive enumeration of the 2²¹ circulants found **21 distinct S**
achieving `K₅_total = 43` (all with |S| ∈ {10, 11}). The next tier
sits at `K₅_total = 172`.

### Parallel warmstart from the top 21

Running v5 for 4 hours from each of these 21 circulants
(84 CPU-hours total) yielded final K₅_total values:

| K₅_final | Count |
|---|---|
| 2 | 13 ranks |
| 3 | 7 ranks |
| 4 | 1 rank |
| 0 or 1 | 0 ranks |

The 13 distinct K₅=2 witnesses split into **two K₄-signatures**:
`K₄_total = 2566` (11 ranks) and `K₄_total = 2567` (2 ranks).
Hash-distinct, possibly grouped into ~3 isomorphism classes.

### Hamming-1/2/3 verdict on the K₅ = 2 plateau

On all 13 witnesses:

- Critical edges: 14 per witness.
- **Hamming-1**: no improving flip exists. Best Δ = 0 (neutral
  moves only).
- **Hamming-2** (C(14, 2) = 91 pairs per witness, 1183 total):
  no pair reduces K₅. Best Δ = 0.
- **Hamming-3** restricted (C(14, 3) = 364 triples per witness,
  4732 total): no triple reduces K₅. Best Δ = 0.
- **Hamming-3 extended** (≈210 candidate edges × C(·, 3) ≈ 1.5M
  triples per witness, tested on 3 witnesses spanning both
  K₄-signatures): no triple reduces K₅. Best Δ = 0.

**The K₅ = 2 wall is structurally robust under all Hamming-3 moves
tested.**

### ℤ₄₂ + center family

Exhaustive enumeration (2²² ≈ 4.2M configurations, ~4 hours on
8 cores) found minimum `K₅_total = 840` (achieved by 12 ranks split
6/6 between central-ON and central-OFF). Starting SA from rank 1
descends to K₅ ≈ 105 in 5 minutes and then stagnates, with reheats
not breaking below.

**Conclusion**: the ℤ₄₂ + center family is structurally
~50× worse than ℤ₄₃ for this problem. The arithmetic property
that 43 is prime (no proper subgroups → no forced K₃/K₄ structures
from cosets) appears to matter substantially.

## Repository structure

```
ramsey-5-5-search/
├── README.md                       (this file)
├── LICENSE                          (MIT)
├── .gitignore
├── scripts/
│   ├── ramsey_55_search.py          main SA driver (v5)
│   ├── ramsey_55_circulant.py       circulant enumeration + verify
│   ├── ramsey_55_z42_center.py      ℤ₄₂+center enumeration + verify
│   ├── ramsey_55_analyze.py         Hamming-1/2/3 analysis
│   └── run_warmstarts.ps1           parallel warmstart (PowerShell)
└── results/
    ├── ramsey55_K2_backup.json      one of the 13 K₅=2 witnesses
    ├── ramsey55_K2_backup.matrix.txt
    ├── ramsey55_K9_backup.json      intermediate K₅=9
    ├── ramsey55_K6_rank4_backup.json   K₅=6 from interrupted run
    ├── ramsey55_circulant_best.json    Circ(43, {1,4,5,6,7,8,9,12,14,17})
    └── ramsey55_checkpoint.json     example checkpoint format
```

## Reproducibility

### Requirements

- Python 3.10+ (uses `int.bit_count()`)
- No external Python dependencies (stdlib only)
- PowerShell 5+ for the parallel runner (Windows) — bash equivalents
  trivial to write

### Workflow

```bash
# 1. Enumerate top circulants on ℤ₄₃
python scripts/ramsey_55_circulant.py enumerate --top_k 50
# → creates circulants_top/rank_01.json .. rank_50.json + _index.json

# 2. Run 21 parallel warmstarts (4h each, on 8 cores)
.\scripts\run_warmstarts.ps1 -NumStarts 21 -Parallel 8 -Hours 4
# → creates warmstart_runs/rank_01/.. with ramsey55_best.json

# 3. Verify a result
python scripts/ramsey_55_circulant.py verify warmstart_runs/rank_01/ramsey55_best.json

# 4. Hamming analysis on a K₅=2 witness
python scripts/ramsey_55_analyze.py results/ramsey55_K2_backup.json --pairs --triples
python scripts/ramsey_55_analyze.py results/ramsey55_K2_backup.json --triples_extended
```

### Single SA run (fastest demo)

```bash
# 5-minute SA from a circulant (should reach K₅ ≈ 5–15)
python scripts/ramsey_55_search.py \
    --start_from results/ramsey55_circulant_best.json \
    --time 300 --t_start 8.0
```

### Alternative family (ℤ₄₂ + center)

```bash
# Test a small subset first
python scripts/ramsey_55_z42_center.py enumerate --top_k 5 --max_codes 5000
# Full enumeration (~4 hours on 8 cores)
python scripts/ramsey_55_z42_center.py enumerate --top_k 50
```

## File formats

All result files are JSON with:

```json
{
  "N": 43,
  "adj_bitmasks": ["5483129885682", "2170168846309", ...],
  "info": {
    "k5": 2,
    "k4": 2566,
    "cost": 38690.0,
    "iter": 2699,
    "tag": "improvement",
    ...
  }
}
```

`adj_bitmasks[v]` is the bitmask of v's neighbors, stored as a
string to avoid JSON's 53-bit int limit. The `.matrix.txt`
companions hold the same data as a 43×43 0/1 matrix.

## Known limitations

- All methods here are **local heuristics**. Reaching K₅ = 0 (or
  proving it impossible) requires either better non-local search
  (SAT solving, integer programming, or non-Cayley structured
  searches) or proving R(5,5) = 43 — neither attempted here.
- The Cayley restriction (whether ℤ₄₃ or ℤ₄₂+center) means we only
  explore a small algebraic slice of all 2⁹⁰³ ≈ 10²⁷² graphs.
  Witnesses outside this slice could exist.
- Analysis is restricted to Hamming-3 in critical-edge neighborhoods.
  Deeper moves were not attempted, but the universality of the K₅=2
  wall across 13 distinct witnesses makes them unlikely to help.

## Possible directions for future work

For anyone who wants to push further:

- **SAT encoding** of "no K₅ in G or Ḡ on 43 vertices" → ~2M
  clauses, attempt with Kissat or CaDiCaL plus aggressive symmetry
  breaking. This is the standard route to a rigorous proof.
- **Non-abelian Cayley graphs** of order near 43 (D_{21}, semidirect
  products, etc.).
- **Construction by fusion**: stitch together two known smaller
  K₅-free graphs along a bipartite layer.
- **Machine learning** (GNN or RL) trained on SA trajectories to
  propose multi-edge moves outside the Hamming-3 neighborhood.

## Citation

If you find this useful (e.g. as a starting toolkit for similar
combinatorial Ramsey problems), please cite as:

```
Faivre, E. (2026). ramsey-5-5-search: heuristic exploration
of K₅-free graphs on 43 vertices. GitHub repository.
https://github.com/Quasaria/ramsey-5-5-search
```

## License

MIT — see `LICENSE`.

---

## Résumé en français

Ce dépôt contient une **recherche heuristique** de graphes à
43 sommets sans K₅ monochromatique — c'est-à-dire une tentative
expérimentale de montrer **R(5,5) ≥ 44**, où R(5,5) est le nombre
de Ramsey à deux couleurs pour K₅.

### Ce que le code fait

Le programme principal `ramsey_55_search.py` implémente un recuit
simulé sur représentation bitset (43 sommets ↔ 43 entiers de 43 bits
chacun) avec moves multi-types (1-flip, 2-flip avec lookahead,
vertex-block), tabou, ILS, et reprise depuis fichier. Performances :
≈ 5000–8000 itérations/sec sur 1 cœur.

`ramsey_55_circulant.py` énumère exhaustivement les 2²¹ ≈ 2 millions
de graphes circulants Cayley sur ℤ₄₃ (groupe de symétrie cyclique
d'ordre 43, premier). Les meilleurs servent de points de départ
pour le SA.

`ramsey_55_z42_center.py` fait la même chose pour une famille
algébriquement différente : graphes Cayley sur ℤ₄₂ avec un sommet
central isolé ou universel.

`ramsey_55_analyze.py` analyse finement un graphe trouvé : énumère
les K₅ explicites, identifie les arêtes critiques, calcule les
deltas exacts pour tous les flips Hamming-1, Hamming-2, et
Hamming-3 (au choix).

### Résultats numériques

- Sur la famille ℤ₄₃ : descente jusqu'à **K₅_total = 2** sur 13
  graphes témoins distincts. Mur dur sous tous les flips
  Hamming-1/2/3 testés.
- Sur la famille ℤ₄₂ + centre : minimum **K₅_total = 840** après
  énumération exhaustive de 4 millions de configurations.
  Comparaison numérique frappante : la primalité de 43 compte.
- **Aucun témoin de R(5,5) ≥ 44 trouvé.** Compatible avec la
  conjecture R(5,5) = 43, mais ce travail ne prouve rien — seule
  une approche SAT exhaustive pourrait trancher.

### Statut

Projet personnel d'optimisation combinatoire. R(5,5) est un problème
ouvert depuis les années 1990, attaqué par de nombreuses équipes
académiques avec des moyens bien supérieurs. Ce dépôt ne prétend pas
à la nouveauté mathématique — sa valeur est dans la **toolkit
reproductible** et dans la comparaison numérique entre les deux
familles algébriques.

### Utilisation rapide

```bash
# Trouver les meilleurs circulants ℤ₄₃ (~5 minutes sur 8 cœurs)
python scripts/ramsey_55_circulant.py enumerate --top_k 50

# Recuit simulé de 5 minutes depuis le meilleur
python scripts/ramsey_55_search.py \
    --start_from results/ramsey55_circulant_best.json \
    --time 300

# Vérifier un résultat
python scripts/ramsey_55_circulant.py verify results/ramsey55_K2_backup.json
```
