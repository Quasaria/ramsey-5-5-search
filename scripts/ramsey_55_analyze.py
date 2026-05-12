#!/usr/bin/env python3
"""
ramsey_55_analyze.py — Analyse fine d'un graphe à K5 résiduels faibles
=======================================================================

Étant donné un graphe sauvegardé (au format ramsey_55_search.py), liste :
  1. Les K5 explicites de G et du complément (5-uplets de sommets).
  2. La distribution des arêtes dans ces K5 (combien de K5 partagent
     chaque arête).
  3. Pour chaque arête appartenant à au moins un K5, le delta exact
     du flip : combien de K5 retirés / créés dans G et complément.
  4. Idem pour les paires d'arêtes (lookahead 2-flip exact, restreint
     aux paires utiles).

Usage :
    python ramsey_55_analyze.py ramsey55_K2_backup.json
    python ramsey_55_analyze.py ramsey55_K2_backup.json --pairs
"""
import sys
import json
import argparse
import itertools

# Force UTF-8 sur stdout/stderr pour eviter UnicodeEncodeError sous Windows
# quand on redirige vers un fichier via PowerShell (qui passe en cp1252).
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except AttributeError:
    pass

N = 43
MASK = (1 << N) - 1

if sys.version_info >= (3, 10):
    popcount = int.bit_count
else:
    def popcount(x):  # noqa: E301
        return bin(x).count("1")

def higher_mask(w):
    return MASK & ~((1 << (w + 1)) - 1)

def complement_of(adj):
    return [(MASK ^ adj[v]) & ~(1 << v) for v in range(N)]

# -------------------------------------------------------------
# Énumération explicite des K5 dans G
# -------------------------------------------------------------
def enumerate_k5(adj):
    """Retourne la liste des 5-uplets (a,b,c,d,e), a<b<c<d<e, qui forment
    un K5 dans le graphe d'adjacence `adj`."""
    out = []
    for a in range(N):
        Na = adj[a] & higher_mask(a)
        bits_b = Na
        while bits_b:
            low_b = bits_b & -bits_b
            b = low_b.bit_length() - 1
            bits_b ^= low_b
            Nab = Na & adj[b]
            bits_c = Nab & higher_mask(b)
            while bits_c:
                low_c = bits_c & -bits_c
                c = low_c.bit_length() - 1
                bits_c ^= low_c
                Nabc = Nab & adj[c]
                bits_d = Nabc & higher_mask(c)
                while bits_d:
                    low_d = bits_d & -bits_d
                    d = low_d.bit_length() - 1
                    bits_d ^= low_d
                    Nabcd = Nabc & adj[d]
                    bits_e = Nabcd & higher_mask(d)
                    while bits_e:
                        low_e = bits_e & -bits_e
                        e = low_e.bit_length() - 1
                        bits_e ^= low_e
                        out.append((a, b, c, d, e))
    return out

def edges_in_in(S, adj):
    """Nombre d'arêtes dans S induit (helper)."""
    cnt = 0
    bits = S
    while bits:
        low = bits & -bits
        w = low.bit_length() - 1
        bits ^= low
        cnt += popcount(adj[w] & S & higher_mask(w))
    return cnt

def triangles_in(S, adj):
    cnt = 0
    bits = S
    while bits:
        low = bits & -bits
        w = low.bit_length() - 1
        bits ^= low
        Nw = adj[w] & S & higher_mask(w)
        nb = Nw
        while nb:
            low2 = nb & -nb
            x = low2.bit_length() - 1
            nb ^= low2
            cnt += popcount(adj[x] & Nw & higher_mask(x))
    return cnt

def delta_counts(u, v, is_edge, adj, comp):
    common_in  = adj[u]  & adj[v]
    common_out = comp[u] & comp[v]
    T_G = triangles_in(common_in,  adj)
    T_C = triangles_in(common_out, comp)
    if is_edge:
        return (-T_G + T_C)  # delta K5 total
    else:
        return ( T_G - T_C)

def flip_edge(u, v, adj, comp):
    bv, bu = 1 << v, 1 << u
    adj[u]  ^= bv;  adj[v]  ^= bu
    comp[u] ^= bv;  comp[v] ^= bu

# -------------------------------------------------------------
# Programme principal
# -------------------------------------------------------------
def main():
    p = argparse.ArgumentParser()
    p.add_argument("file")
    p.add_argument("--pairs", action="store_true",
                   help="énumère aussi les paires (2-flip) prometteuses")
    p.add_argument("--triples", action="store_true",
                   help="énumère aussi les triplets (3-flip) — peut être lent")
    p.add_argument("--triples_extended", action="store_true",
                   help="--triples mais étend la liste des arêtes candidates "
                        "en y ajoutant celles qui touchent un sommet critique "
                        "(plus lent, mais plus de chances de trouver K5=0)")
    args = p.parse_args()

    with open(args.file, "r") as f:
        data = json.load(f)
    if data.get("N") != N:
        print(f"N mismatch: got {data.get('N')}, expected {N}")
        return
    adj = [int(x) for x in data["adj_bitmasks"]]
    comp = complement_of(adj)

    # === Étape 1 : énumération explicite des K5 ===
    k5_g = enumerate_k5(adj)
    k5_c = enumerate_k5(comp)
    print(f"=== ÉNUMÉRATION DES K5 ===")
    print(f"  K5 dans G        : {len(k5_g)}")
    for clique in k5_g:
        print(f"    {clique}")
    print(f"  K5 dans complément : {len(k5_c)}")
    for clique in k5_c:
        print(f"    {clique}")

    # === Étape 2 : arêtes "critiques" (qui appartiennent à au moins un K5) ===
    print(f"\n=== ARÊTES PRÉSENTES DANS DES K5 ===")
    edge_to_k5g = {}  # (u,v) -> liste de K5 (G) qui contiennent (u,v)
    edge_to_k5c = {}  # idem complément
    for clique in k5_g:
        for u, v in itertools.combinations(clique, 2):
            edge_to_k5g.setdefault((u, v), []).append(clique)
    for clique in k5_c:
        for u, v in itertools.combinations(clique, 2):
            edge_to_k5c.setdefault((u, v), []).append(clique)

    # arêtes de G présentes dans des K5 de G : leur retrait casse ces K5
    print(f"\n  Arêtes de G dans des K5 de G  ({len(edge_to_k5g)} arêtes) :")
    for (u, v), cliques in sorted(edge_to_k5g.items()):
        print(f"    ({u:2d},{v:2d}) ∈ {len(cliques)} K5(G)")

    # non-arêtes de G présentes dans des K5 du complément : leur ajout casserait ces K5(C)
    print(f"\n  Non-arêtes de G dans des K5 du complément  ({len(edge_to_k5c)} arêtes) :")
    for (u, v), cliques in sorted(edge_to_k5c.items()):
        print(f"    ({u:2d},{v:2d}) ∈ {len(cliques)} K5(complément)")

    # === Étape 3 : delta exact pour chaque flip d'arête critique ===
    print(f"\n=== DELTA POUR CHAQUE FLIP HAMMING-1 ===")
    print(f"(format: (u,v)  is_edge  delta_K5_total  →  K5 final si flippé)")

    # On regroupe : flips qui pourraient potentiellement réduire K5
    candidates = set(edge_to_k5g.keys()) | set(edge_to_k5c.keys())

    current_k5 = len(k5_g) + len(k5_c)
    print(f"K5 total courant = {current_k5}")
    best_delta = 0
    best_moves = []
    for (u, v) in sorted(candidates):
        is_edge = (adj[u] >> v) & 1
        d = delta_counts(u, v, is_edge, adj, comp)
        new_k5 = current_k5 + d
        marker = ""
        if new_k5 < current_k5:
            marker = "  ← AMÉLIORE"
            if d < best_delta:
                best_delta = d
                best_moves = [(u, v, is_edge, d)]
            elif d == best_delta:
                best_moves.append((u, v, is_edge, d))
        if new_k5 == 0:
            marker = "  ★★★ K5=0 EN UN FLIP ! ★★★"
        print(f"  ({u:2d},{v:2d})  is_edge={is_edge}  d={d:+3d}  →  K5={new_k5}{marker}")

    if best_delta < 0:
        print(f"\nMeilleur 1-flip : delta = {best_delta} (K5 → {current_k5 + best_delta})")
        for (u, v, ie, d) in best_moves:
            print(f"  flip ({u},{v}) [is_edge={ie}]")
    else:
        print(f"\n!! Aucun 1-flip n'améliore — minimum local Hamming-1 confirmé.")

    # === Étape 4 : 2-flip exhaustif sur les arêtes critiques ===
    if args.pairs:
        print(f"\n=== 2-FLIP EXHAUSTIF (paires d'arêtes critiques) ===")
        cand_list = sorted(candidates)
        n_cand = len(cand_list)
        print(f"({n_cand} arêtes critiques → {n_cand*(n_cand-1)//2} paires à tester)")

        best_pair_delta = 0
        best_pairs = []
        zero_pairs = []
        for i in range(n_cand):
            u1, v1 = cand_list[i]
            is_edge1 = (adj[u1] >> v1) & 1
            d1 = delta_counts(u1, v1, is_edge1, adj, comp)
            flip_edge(u1, v1, adj, comp)
            for j in range(i+1, n_cand):
                u2, v2 = cand_list[j]
                is_edge2 = (adj[u2] >> v2) & 1
                d2 = delta_counts(u2, v2, is_edge2, adj, comp)
                tot_d = d1 + d2
                new_k5 = current_k5 + tot_d
                if new_k5 == 0:
                    zero_pairs.append(((u1, v1), (u2, v2), is_edge1, is_edge2))
                if tot_d < best_pair_delta:
                    best_pair_delta = tot_d
                    best_pairs = [((u1, v1), (u2, v2), tot_d)]
                elif tot_d == best_pair_delta and len(best_pairs) < 10:
                    best_pairs.append(((u1, v1), (u2, v2), tot_d))
            flip_edge(u1, v1, adj, comp)  # undo

        if zero_pairs:
            print(f"\n★★★ {len(zero_pairs)} PAIRE(S) MENANT À K5=0 :")
            for (e1, e2, ie1, ie2) in zero_pairs[:20]:
                print(f"  flip {e1} (is_edge={ie1}) + flip {e2} (is_edge={ie2})")
        else:
            print(f"\nAucune paire ne donne K5=0.")
            print(f"Meilleur 2-flip : delta = {best_pair_delta} "
                  f"(K5 → {current_k5 + best_pair_delta})")
            for (e1, e2, d) in best_pairs[:5]:
                print(f"  flip {e1} + flip {e2}")

    # === Étape 5 : 3-flip exhaustif (Hamming-3) ===
    if args.triples or args.triples_extended:
        print(f"\n=== 3-FLIP EXHAUSTIF (triplets d'arêtes) ===")

        # Construire la liste des arêtes candidates
        cand_set = set(edge_to_k5g.keys()) | set(edge_to_k5c.keys())

        if args.triples_extended:
            # Élargir : on ajoute toutes les arêtes (et non-arêtes) incidentes
            # à un sommet apparaissant dans au moins un K5
            critical_vertices = set()
            for clique in k5_g:
                critical_vertices.update(clique)
            for clique in k5_c:
                critical_vertices.update(clique)
            print(f"(extended : {len(critical_vertices)} sommets critiques)")
            for u in critical_vertices:
                for v in range(N):
                    if v == u:
                        continue
                    a, b = (u, v) if u < v else (v, u)
                    cand_set.add((a, b))

        cand_list = sorted(cand_set)
        n_cand = len(cand_list)
        n_triples = n_cand * (n_cand - 1) * (n_cand - 2) // 6
        print(f"({n_cand} arêtes candidates → {n_triples} triplets à tester)")

        if n_triples > 200_000:
            print(f"!! Beaucoup de triplets ({n_triples}). Peut prendre du temps.")

        best_triple_delta = 0
        best_triples = []
        zero_triples = []

        for i in range(n_cand):
            u1, v1 = cand_list[i]
            is_edge1 = (adj[u1] >> v1) & 1
            d1 = delta_counts(u1, v1, is_edge1, adj, comp)
            flip_edge(u1, v1, adj, comp)
            for j in range(i + 1, n_cand):
                u2, v2 = cand_list[j]
                is_edge2 = (adj[u2] >> v2) & 1
                d2 = delta_counts(u2, v2, is_edge2, adj, comp)
                flip_edge(u2, v2, adj, comp)
                for k in range(j + 1, n_cand):
                    u3, v3 = cand_list[k]
                    is_edge3 = (adj[u3] >> v3) & 1
                    d3 = delta_counts(u3, v3, is_edge3, adj, comp)
                    tot_d = d1 + d2 + d3
                    new_k5 = current_k5 + tot_d
                    if new_k5 == 0:
                        zero_triples.append((
                            (u1, v1), (u2, v2), (u3, v3),
                            is_edge1, is_edge2, is_edge3,
                        ))
                    if tot_d < best_triple_delta:
                        best_triple_delta = tot_d
                        best_triples = [((u1, v1), (u2, v2), (u3, v3), tot_d)]
                    elif tot_d == best_triple_delta and len(best_triples) < 10:
                        best_triples.append(((u1, v1), (u2, v2), (u3, v3), tot_d))
                flip_edge(u2, v2, adj, comp)  # undo e_j
            flip_edge(u1, v1, adj, comp)  # undo e_i

        if zero_triples:
            print(f"\n★★★ {len(zero_triples)} TRIPLET(S) MENANT À K5=0 :")
            for (e1, e2, e3, ie1, ie2, ie3) in zero_triples[:20]:
                print(f"  flip {e1} (is_edge={ie1}) + "
                      f"flip {e2} (is_edge={ie2}) + "
                      f"flip {e3} (is_edge={ie3})")
        else:
            print(f"\nAucun triplet ne donne K5=0.")
            print(f"Meilleur 3-flip : delta = {best_triple_delta} "
                  f"(K5 → {current_k5 + best_triple_delta})")
            for (e1, e2, e3, d) in best_triples[:5]:
                print(f"  flip {e1} + flip {e2} + flip {e3}")

if __name__ == "__main__":
    main()
