#!/usr/bin/env python3
"""
ramsey_55_z42_center.py - Enumeration des graphes "Cayley Z_42 + sommet central"
================================================================================

Construction :
  - Sommets 0..41 organises en cycle Z_42 avec structure Cayley : pour S ⊂ Z_42\\{0}
    symetrique (S = -S mod 42), on a u ~ v ssi (v-u mod 42) ∈ S.
  - Sommet 42 = "centre", connecte a AUCUN ou a TOUS les sommets de l'anneau.

C'est une famille differente des circulants Z_43, parce que Z_42 a un sous-groupe
de symetrie d'ordre 42 (au lieu de 43), et le sommet central distingue topologi-
quement le centre. Les K_5 doivent respecter cette double structure.

Espace de recherche : 2^21 configurations pour S × 2 choix pour central
                      = 4,194,304 graphes au total.

Modes :
  enumerate : enumeration exhaustive multi-coeurs, top-K sauve
  verify    : verification rigoureuse d'un fichier .json
  list      : afficher le top-K deja sauve
"""
import sys
import os
import time
import json
import argparse
import multiprocessing

# Force UTF-8 sur stdout/stderr (probleme PowerShell Windows)
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except AttributeError:
    pass

# -----------------------------------------------------------------------------
# Constantes
# -----------------------------------------------------------------------------
N_RING = 42         # sommets 0..41 forment l'anneau Z_42
N_TOTAL = 43        # sommets 0..42, le 42 est le centre
CENTER = 42
MASK = (1 << N_TOTAL) - 1

if sys.version_info >= (3, 10):
    popcount = int.bit_count
else:
    def popcount(x):
        return bin(x).count("1")


def higher_mask(w):
    return MASK & ~((1 << (w + 1)) - 1)


# -----------------------------------------------------------------------------
# Comptage K_4 et K_5 par bitset
# -----------------------------------------------------------------------------
def count_k4(adj):
    cnt = 0
    for a in range(N_TOTAL):
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
                cnt += popcount(Nabc & higher_mask(c))
    return cnt


def count_k5(adj):
    cnt = 0
    for a in range(N_TOTAL):
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
                    cnt += popcount(Nabcd & higher_mask(d))
    return cnt


def complement_of(adj):
    return [(MASK ^ adj[v]) & ~(1 << v) for v in range(N_TOTAL)]


# -----------------------------------------------------------------------------
# Encodage / decodage des configurations
# -----------------------------------------------------------------------------
def build_representatives():
    """Liste des orbites de Z_42\\{0} sous l'involution a -> -a mod 42.

    Resultat : 20 paires (a, -a) + 1 auto-inverse (21,) = 21 elements.
    Chaque entree est un tuple : soit (a, -a), soit (a,) si auto-inverse.
    """
    seen = set()
    reps = []
    for a in range(1, N_RING):
        if a in seen:
            continue
        inv = (-a) % N_RING
        if inv == a:
            reps.append((a,))
            seen.add(a)
        else:
            reps.append((a, inv))
            seen.add(a)
            seen.add(inv)
    return reps


REPS = build_representatives()
assert len(REPS) == 21, f"Erreur de construction des representants : {len(REPS)} au lieu de 21"


def code_to_S(code):
    """Retourne le set S correspondant au code 21 bits."""
    S = set()
    for bit_idx, rep in enumerate(REPS):
        if (code >> bit_idx) & 1:
            for x in rep:
                S.add(x)
    return S


def build_adj(code, central_on):
    """Construit la liste adj (43 bitmasks) pour la configuration donnee."""
    S = code_to_S(code)
    adj = [0] * N_TOTAL
    # Aretes dans l'anneau Z_42
    for u in range(N_RING):
        for d in S:
            v = (u + d) % N_RING
            if u < v:
                adj[u] |= (1 << v)
                adj[v] |= (1 << u)
    # Sommet central (42) connecte a TOUS les sommets de l'anneau si central_on
    if central_on:
        for v in range(N_RING):
            adj[CENTER] |= (1 << v)
            adj[v] |= (1 << CENTER)
    return adj


# -----------------------------------------------------------------------------
# Evaluation d'un graphe
# -----------------------------------------------------------------------------
def evaluate(code, central_on, w_k5=100, w_k4=10):
    adj = build_adj(code, central_on)
    comp = complement_of(adj)
    k5g = count_k5(adj)
    k5c = count_k5(comp)
    k4g = count_k4(adj)
    k4c = count_k4(comp)
    k5_total = k5g + k5c
    k4_total = k4g + k4c
    cost = w_k5 * k5_total + w_k4 * k4_total
    return cost, k5_total, k4_total


# -----------------------------------------------------------------------------
# Worker pour multiprocessing
# -----------------------------------------------------------------------------
def worker_chunk(args):
    start, end, top_k, w_k5, w_k4 = args
    local_top = []
    for code in range(start, end):
        for central_on in (False, True):
            cost, k5, k4 = evaluate(code, central_on, w_k5, w_k4)
            local_top.append((cost, k5, k4, code, central_on))
        if len(local_top) > 4 * top_k:
            local_top.sort()
            local_top = local_top[:top_k]
    local_top.sort()
    return local_top[:top_k]


def enumerate_all(top_k=50, num_workers=None, w_k5=100, w_k4=10,
                  output_dir="z42_top", max_codes=None):
    total = 1 << 21
    if max_codes is not None and max_codes < total:
        total = max_codes
        print(f"[enumerate] LIMITE A {max_codes} codes (mode test)")
    if num_workers is None:
        num_workers = multiprocessing.cpu_count()

    n_chunks = max(num_workers * 16, 64)
    chunk_size = max(1, total // n_chunks + 1)
    chunks = []
    for i in range(n_chunks):
        s = i * chunk_size
        e = min((i + 1) * chunk_size, total)
        if s < total:
            chunks.append((s, e, top_k, w_k5, w_k4))

    print(f"[enumerate] N=43 (Z_42 + centre)")
    print(f"            espace=2^21 * 2 = {2 * total}")
    print(f"            workers={num_workers}  chunks={len(chunks)} de ~{chunk_size}")
    print(f"            top_k={top_k}", flush=True)

    start_time = time.time()
    global_top = []

    with multiprocessing.Pool(num_workers) as pool:
        for chunk_idx, result in enumerate(pool.imap_unordered(worker_chunk, chunks), 1):
            global_top.extend(result)
            global_top.sort()
            global_top = global_top[:top_k]

            elapsed = time.time() - start_time
            done_frac = chunk_idx / len(chunks)
            eta = elapsed / done_frac * (1 - done_frac) if done_frac > 0 else 0
            best = global_top[0] if global_top else None
            if best is not None:
                cost, k5, k4, code, central_on = best
                ctr = "ON " if central_on else "OFF"
                rank50 = global_top[-1][0] if len(global_top) == top_k else "-"
                print(f"[{chunk_idx:4d}/{len(chunks)}]  "
                      f"{done_frac*100:5.1f}%  #1 K5={k5} K4={k4} cost={cost} "
                      f"ctr={ctr} (rank{top_k}: cost={rank50})  ETA={int(eta)}s",
                      flush=True)

    os.makedirs(output_dir, exist_ok=True)
    index_entries = []
    for rank, (cost, k5, k4, code, central_on) in enumerate(global_top, 1):
        adj = build_adj(code, central_on)
        filename = f"rank_{rank:02d}.json"
        S_list = sorted(code_to_S(code))
        data = {
            "N": N_TOTAL,
            "adj_bitmasks": [str(x) for x in adj],
            "info": {
                "tag": "z42_center_cayley",
                "rank": rank,
                "k5": k5,
                "k4": k4,
                "cost": cost,
                "code": code,
                "central_on": central_on,
                "S": S_list,
            }
        }
        with open(os.path.join(output_dir, filename), "w") as f:
            json.dump(data, f)
        index_entries.append({
            "rank": rank,
            "filename": filename,
            "k5_total": k5,
            "k4_total": k4,
            "cost": cost,
            "code": code,
            "central_on": central_on,
            "S": S_list,
        })

    with open(os.path.join(output_dir, "_index.json"), "w") as f:
        json.dump({"entries": index_entries}, f, indent=2)

    elapsed = time.time() - start_time
    print(f"\n[done] {int(elapsed)}s - top-{top_k} sauvegarde dans {output_dir}/")
    for entry in index_entries[:20]:
        ctr = "ON " if entry["central_on"] else "OFF"
        print(f"  rank {entry['rank']:2d}: K5={entry['k5_total']:5d}  "
              f"K4={entry['k4_total']:5d}  cost={entry['cost']:8d}  "
              f"central={ctr}  |S|={len(entry['S']):2d}  S={entry['S']}")
    print(f"\nIndex sauvegarde : {output_dir}/_index.json")


def verify_file(path):
    with open(path, "r") as f:
        data = json.load(f)
    if data.get("N") != N_TOTAL:
        print(f"N mismatch: got {data.get('N')}, expected {N_TOTAL}")
        return
    adj = [int(x) for x in data["adj_bitmasks"]]
    comp = complement_of(adj)
    k5g = count_k5(adj)
    k5c = count_k5(comp)
    k4g = count_k4(adj)
    k4c = count_k4(comp)
    print(f"[verify] type=Z_42 + centre (brute-force)")
    print(f"  K5(G)={k5g}  K5(complement)={k5c}  K5_total={k5g+k5c}")
    print(f"  K4(G)={k4g}  K4(complement)={k4c}  K4_total={k4g+k4c}")


def list_top(output_dir="z42_top"):
    index_path = os.path.join(output_dir, "_index.json")
    if not os.path.exists(index_path):
        print(f"Pas d'index trouve : {index_path}")
        return
    with open(index_path) as f:
        index = json.load(f)
    print(f"Top du dossier {output_dir}/ :")
    print(f"  rank    K5     K4     cost  central |S|  S")
    print(f"  ----  -----  -----  -------  ------- ---  --------------------------------")
    for entry in index["entries"]:
        ctr = "ON " if entry["central_on"] else "OFF"
        print(f"  {entry['rank']:4d}  {entry['k5_total']:5d}  {entry['k4_total']:5d}  "
              f"{entry['cost']:7d}  {ctr:7s} {len(entry['S']):3d}  {entry['S']}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd")

    pe = sub.add_parser("enumerate")
    pe.add_argument("--top_k", type=int, default=50)
    pe.add_argument("--workers", type=int, default=0)
    pe.add_argument("--w_k5", type=int, default=100)
    pe.add_argument("--w_k4", type=int, default=10)
    pe.add_argument("--max_codes", type=int, default=0,
                    help="Si >0, limite l'enumeration aux premiers N codes "
                         "(utile pour test rapide). 0 = espace complet.")

    pv = sub.add_parser("verify")
    pv.add_argument("file")

    pl = sub.add_parser("list")

    args = p.parse_args()

    if args.cmd == "enumerate":
        nw = args.workers if args.workers > 0 else None
        mc = args.max_codes if args.max_codes > 0 else None
        enumerate_all(top_k=args.top_k, num_workers=nw,
                      w_k5=args.w_k5, w_k4=args.w_k4, max_codes=mc)
    elif args.cmd == "verify":
        verify_file(args.file)
    elif args.cmd == "list":
        list_top()
    else:
        p.print_help()
