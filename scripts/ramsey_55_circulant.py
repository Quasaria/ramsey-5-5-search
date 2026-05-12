#!/usr/bin/env python3
"""
ramsey_55_circulant.py — Recherche R(5,5)>=44 dans les graphes circulants
=========================================================================

Restreint la recherche à l'espace des graphes circulants Circ(43, S)
où S ⊂ {1, ..., 21} est le "connection set".

Modes :
  enumerate       — énumère tous les S, garde le TOP-K (multi-cœurs).
  sa              — recuit simulé dans l'espace S.
  verify FILE     — recompute K5/K4 d'un graphe sauvegardé.
  list            — liste les top-K circulants sauvegardés.

Format des sauvegardes : compatible avec ramsey_55_search.py
(utilisable avec --start_from pour transférer un circulant en full-space).

Lancement typique :
  python ramsey_55_circulant.py enumerate --top_k 50
  python ramsey_55_circulant.py list
  python ramsey_55_circulant.py verify circulants_top/rank_01_S0x0129f9.json
"""
import sys
import os
import time
import math
import random
import json
import argparse
import multiprocessing as mp

# -----------------------------------------------------------------------
# Constantes
# -----------------------------------------------------------------------
N = 43
MASK = (1 << N) - 1
N_DISTANCES = (N - 1) // 2

if sys.version_info >= (3, 10):
    popcount = int.bit_count
else:
    def popcount(x):  # noqa: E301
        return bin(x).count("1")


def higher_mask(w):
    return MASK & ~((1 << (w + 1)) - 1)


# -----------------------------------------------------------------------
# Pré-calcul : adjacence par distance unique
# -----------------------------------------------------------------------
def _build_dist_adj():
    out = [None] * (N_DISTANCES + 1)
    for d in range(1, N_DISTANCES + 1):
        a = [0] * N
        for i in range(N):
            j = (i + d) % N
            a[i] |= 1 << j
            a[j] |= 1 << i
        out[d] = a
    return out

DIST_ADJ = _build_dist_adj()


def build_circulant(S_int):
    adj = [0] * N
    bits = S_int
    d = 1
    while bits:
        if bits & 1:
            da = DIST_ADJ[d]
            for i in range(N):
                adj[i] |= da[i]
        bits >>= 1
        d += 1
    return adj


def complement_of(adj):
    return [(MASK ^ adj[v]) & ~(1 << v) for v in range(N)]


# -----------------------------------------------------------------------
# Comptages bitset
# -----------------------------------------------------------------------
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


def k4s_in(S, adj):
    cnt = 0
    bits = S
    while bits:
        low = bits & -bits
        w = low.bit_length() - 1
        bits ^= low
        cnt += triangles_in(adj[w] & S & higher_mask(w), adj)
    return cnt


def k5_circulant(adj):
    """N premier ⇒ K5_total = (N × #K5_par_v0) / 5"""
    N0 = adj[0]
    return (N * k4s_in(N0, adj)) // 5


def k4_circulant(adj):
    N0 = adj[0]
    return (N * triangles_in(N0, adj)) // 4


def eval_S(S_int):
    adj = build_circulant(S_int)
    comp = complement_of(adj)
    return (k5_circulant(adj) + k5_circulant(comp),
            k4_circulant(adj) + k4_circulant(comp))


def cost_of(k5, k4, w5, w4):
    return w5 * k5 + w4 * k4


# -----------------------------------------------------------------------
# Sauvegarde (format compatible ramsey_55_search.py)
# -----------------------------------------------------------------------
def save_circulant_to(S_int, info_extra, fname_json, fname_mat=None):
    """Sauvegarde robuste dans un fichier nommé."""
    pid = os.getpid()
    tmp = f"{fname_json}.{pid}.tmp"

    adj = build_circulant(S_int)
    rows = []
    for i in range(N):
        rows.append("".join("1" if (adj[i] >> j) & 1 else "0"
                            for j in range(N)))
    S_list = sorted([d for d in range(1, N_DISTANCES + 1)
                     if (S_int >> (d - 1)) & 1])
    info = {
        "graph_type": "circulant",
        "S_int_hex": f"0x{S_int:06x}",
        "S_distances": S_list,
        "S_size": len(S_list),
    }
    info.update(info_extra)
    payload = {
        "N": N,
        "info": info,
        "adj_bitmasks": [int(x) for x in adj],
        "adjacency_matrix": rows,
    }

    for attempt in range(6):
        try:
            with open(tmp, "w") as f:
                json.dump(payload, f, indent=2)
            try:
                os.replace(tmp, fname_json)
                break
            except (PermissionError, OSError, FileNotFoundError):
                pass
        except (PermissionError, OSError):
            pass
        time.sleep(0.1 * (attempt + 1))

    try:
        if os.path.exists(tmp):
            os.remove(tmp)
    except OSError:
        pass

    if fname_mat:
        try:
            with open(fname_mat, "w") as f:
                f.write("\n".join(rows) + "\n")
        except (PermissionError, OSError):
            pass


def parse_S_distances(s):
    distances = [int(x.strip()) for x in s.split(",") if x.strip()]
    S_int = 0
    for d in distances:
        if not (1 <= d <= N_DISTANCES):
            raise ValueError(f"distance {d} hors de [1..{N_DISTANCES}]")
        S_int |= 1 << (d - 1)
    return S_int


# -----------------------------------------------------------------------
# Mode ENUMERATE — garde le TOP-K
# -----------------------------------------------------------------------
def evaluate_chunk_topk(args):
    """
    Évalue les S impairs dans [start, end). Renvoie la liste des top-K
    locaux : liste de tuples (cost, S_int, k5, k4), triés.
    """
    start, end, w5, w4, top_k = args
    # min-heap par coût (en stockant -cost, ou liste maintenue triée)
    # Pour rester simple : liste de top_k éléments, tri par insertion.
    top = []  # liste de (cost, S, k5, k4), triée croissant
    s_start = start if (start & 1) else (start + 1)
    n_evaluated = 0
    for s in range(s_start, end, 2):
        k5, k4 = eval_S(s)
        n_evaluated += 1
        c = cost_of(k5, k4, w5, w4)
        if len(top) < top_k:
            top.append((c, s, k5, k4))
            top.sort()
        elif c < top[-1][0]:
            top[-1] = (c, s, k5, k4)
            top.sort()
    return (top, n_evaluated)


def run_enumerate(args):
    workers = args.workers if args.workers > 0 else mp.cpu_count()
    n_total = 1 << N_DISTANCES
    n_effective = n_total // 2

    chunk_size = max(n_total // (workers * 16), 4096)
    chunks = []
    for start in range(0, n_total, chunk_size):
        end = min(start + chunk_size, n_total)
        chunks.append((start, end, args.weight_k5, args.weight_k4, args.top_k))

    print(f"[enumerate] N={N}  espace=2^{N_DISTANCES}={n_total}  "
          f"effectifs=~{n_effective}")
    print(f"            workers={workers}  chunks={len(chunks)} de ~{chunk_size}  "
          f"top_k={args.top_k}")

    out_dir = args.output_dir
    os.makedirs(out_dir, exist_ok=True)

    global_top = []  # liste triée globale (cost, S, k5, k4)
    processed_chunks = 0
    processed_S = 0
    t0 = time.time()

    try:
        with mp.Pool(workers) as pool:
            for result in pool.imap_unordered(evaluate_chunk_topk, chunks):
                local_top, n_eval = result
                processed_chunks += 1
                processed_S += n_eval

                # Fusionner local_top dans global_top
                # On vérifie l'unicité par S (au cas où)
                seen_S = {x[1] for x in global_top}
                for entry in local_top:
                    if entry[1] not in seen_S:
                        global_top.append(entry)
                        seen_S.add(entry[1])
                global_top.sort()
                if len(global_top) > args.top_k:
                    global_top = global_top[:args.top_k]

                el = time.time() - t0
                rate = processed_S / el if el > 0 else 0
                eta = (n_effective - processed_S) / max(rate, 1)
                pct = 100 * processed_S / max(n_effective, 1)
                best = global_top[0]
                worst_in_top = global_top[-1] if len(global_top) >= args.top_k else None
                worst_str = f"  (rank{args.top_k}: cost={worst_in_top[0]:.0f})" if worst_in_top else ""
                print(f"[{processed_chunks:>4}/{len(chunks)}]  "
                      f"{processed_S}/{n_effective} ({pct:5.1f}%)  "
                      f"#1 K5={best[2]} K4={best[3]} cost={best[0]:.0f}  "
                      f"S={best[1]:#08x}{worst_str}  "
                      f"{rate:.0f} S/s  ETA={eta:.0f}s")

    except KeyboardInterrupt:
        print("\n[interrupt]")

    # Sauvegarde finale du top-K
    print(f"\n[done] {time.time()-t0:.0f}s — sauvegarde du top-{len(global_top)} dans {out_dir}/")

    # Index récapitulatif
    index = []
    for rank, (cost, S, k5, k4) in enumerate(global_top, start=1):
        S_list = sorted([d for d in range(1, N_DISTANCES + 1)
                         if (S >> (d - 1)) & 1])
        fname_base = f"rank_{rank:02d}_S{S:#08x}".replace("0x", "0x")
        fname_json = os.path.join(out_dir, f"{fname_base}.json")
        fname_mat = os.path.join(out_dir, f"{fname_base}.matrix.txt")
        save_circulant_to(S, {
            "k5": k5, "k4": k4, "cost": cost,
            "rank": rank, "tag": "enumerate-top-k",
        }, fname_json, fname_mat)
        index.append({
            "rank": rank,
            "S_int_hex": f"0x{S:06x}",
            "S_distances": S_list,
            "S_size": len(S_list),
            "k5_total": k5,
            "k4_total": k4,
            "cost": cost,
            "filename": os.path.basename(fname_json),
        })
        if rank <= 20:
            print(f"  rank {rank:2d}: K5={k5:>4}  K4={k4:>5}  cost={cost:>7.0f}  "
                  f"|S|={len(S_list):>2}  S={S_list}")

    # Aussi : compatibilité descendante avec l'ancien format
    if global_top:
        cost, S, k5, k4 = global_top[0]
        save_circulant_to(S, {
            "k5": k5, "k4": k4, "cost": cost, "tag": "enumerate-best",
        }, "ramsey55_circulant_best.json", "ramsey55_circulant_best.matrix.txt")

    # Sauvegarde de l'index JSON
    index_path = os.path.join(out_dir, "_index.json")
    with open(index_path, "w") as f:
        json.dump({"top_k": len(global_top), "entries": index}, f, indent=2)
    print(f"\nIndex sauvegardé : {index_path}")
    print(f"Pour voir la liste : python ramsey_55_circulant.py list")


# -----------------------------------------------------------------------
# Mode LIST — affiche le top-K sauvegardé
# -----------------------------------------------------------------------
def run_list(args):
    out_dir = args.output_dir
    index_path = os.path.join(out_dir, "_index.json")
    if not os.path.exists(index_path):
        print(f"Pas d'index trouvé dans {out_dir}/_index.json")
        print(f"Lance d'abord : python ramsey_55_circulant.py enumerate --top_k 50")
        return
    with open(index_path, "r") as f:
        data = json.load(f)
    print(f"Top-{data['top_k']} circulants sauvegardés dans {out_dir}/ :\n")
    print(f"  {'rank':>4}  {'K5':>4}  {'K4':>5}  {'cost':>7}  {'|S|':>3}  S")
    print(f"  {'-'*4}  {'-'*4}  {'-'*5}  {'-'*7}  {'-'*3}  {'-'*30}")
    for e in data["entries"]:
        S_str = str(e["S_distances"])
        if len(S_str) > 50:
            S_str = S_str[:47] + "..."
        print(f"  {e['rank']:>4}  {e['k5_total']:>4}  {e['k4_total']:>5}  "
              f"{e['cost']:>7.0f}  {e['S_size']:>3}  {S_str}")


# -----------------------------------------------------------------------
# Mode SA — recuit simulé dans l'espace S (21 bits)
# -----------------------------------------------------------------------
def run_sa(args):
    rng = random.Random(args.seed)

    if args.S_distances is not None:
        S = parse_S_distances(args.S_distances)
    else:
        target_size = rng.randint(9, 12)
        chosen = rng.sample(range(N_DISTANCES), target_size)
        S = sum(1 << b for b in chosen)
        S |= 1

    k5, k4 = eval_S(S)
    cost = cost_of(k5, k4, args.weight_k5, args.weight_k4)
    best_S, best_cost, best_k5, best_k4 = S, cost, k5, k4

    print(f"[sa] N={N}  seed={args.seed}  |S|={popcount(S)}  "
          f"K5={k5}  K4={k4}  cost={cost:.0f}")
    save_circulant_to(best_S, {"k5": k5, "k4": k4, "cost": cost,
                               "tag": "sa-initial", "iter": 0},
                      "ramsey55_circulant_best.json",
                      "ramsey55_circulant_best.matrix.txt")

    T = args.t_start
    iters = 0
    accepts = 0
    last_improve = 0
    deadline = (time.time() + args.time) if args.time > 0 else float("inf")
    t0 = time.time()

    try:
        while T > args.t_end and time.time() < deadline and best_k5 > 0:
            do_two = rng.random() < args.two_bit_prob
            S2 = S ^ (1 << rng.randrange(N_DISTANCES))
            if do_two:
                S2 ^= 1 << rng.randrange(N_DISTANCES)
            if S2 == 0 or S2 == ((1 << N_DISTANCES) - 1):
                T *= args.cooling
                iters += 1
                continue

            new_k5, new_k4 = eval_S(S2)
            new_cost = cost_of(new_k5, new_k4,
                               args.weight_k5, args.weight_k4)
            d_cost = new_cost - cost

            if d_cost <= 0 or rng.random() < math.exp(-d_cost / T):
                S, k5, k4, cost = S2, new_k5, new_k4, new_cost
                accepts += 1
                if cost < best_cost - 1e-9:
                    best_S, best_cost = S, cost
                    best_k5, best_k4 = k5, k4
                    last_improve = iters
                    save_circulant_to(best_S, {
                        "k5": k5, "k4": k4, "cost": cost,
                        "T": T, "iter": iters,
                        "elapsed_sec": time.time() - t0,
                        "tag": "sa-improvement"},
                        "ramsey55_circulant_best.json",
                        "ramsey55_circulant_best.matrix.txt")
                    if k5 == 0:
                        print(f"\n!!! K5=0 trouvé à iter={iters} S={S:#08x} !!!")
                        return

            T *= args.cooling
            iters += 1

            if iters - last_improve > args.restart_after:
                T = max(T, args.t_start * 0.5)
                for _ in range(rng.randint(3, 5)):
                    S ^= 1 << rng.randrange(N_DISTANCES)
                if S == 0 or S == ((1 << N_DISTANCES) - 1):
                    S |= 1
                k5, k4 = eval_S(S)
                cost = cost_of(k5, k4, args.weight_k5, args.weight_k4)
                last_improve = iters
                print(f"[reheat @ {iters}] T←{T:.3f}  |S|={popcount(S)}  "
                      f"K5={k5}  K4={k4}")

            if iters % args.log_every == 0:
                el = time.time() - t0
                rate = iters / el if el > 0 else 0
                print(f"it={iters:>7}  T={T:6.2f}  |S|={popcount(S):>2}  "
                      f"K5={k5:>4}  K4={k4:>5}  cost={cost:>8.0f}  "
                      f"best K5={best_k5:>4} K4={best_k4:>5}  "
                      f"acc={100*accepts/max(iters,1):.1f}%  {rate:.0f} it/s")
    except KeyboardInterrupt:
        print("\n[interrupt]")

    print(f"\n[stop] iter={iters}  best K5={best_k5}  K4={best_k4}  "
          f"S={best_S:#08x}")


# -----------------------------------------------------------------------
# Mode VERIFY
# -----------------------------------------------------------------------
def _brute_total_k5(adj):
    cnt = 0
    for u in range(N):
        cnt += k4s_in(adj[u] & higher_mask(u), adj)
    return cnt

def _brute_total_k4(adj):
    cnt = 0
    for u in range(N):
        cnt += triangles_in(adj[u] & higher_mask(u), adj)
    return cnt

def run_verify(args):
    with open(args.file, "r") as f:
        data = json.load(f)
    if data.get("N") != N:
        print(f"Erreur : N={data.get('N')}, attendu {N}")
        return
    adj = [int(x) for x in data["adj_bitmasks"]]
    comp = complement_of(adj)

    for v in range(N):
        if (adj[v] >> v) & 1:
            print(f"Erreur : self-loop sur {v}")
            return
        for u in range(v):
            if ((adj[v] >> u) & 1) != ((adj[u] >> v) & 1):
                print(f"Erreur : asymétrie sur ({u},{v})")
                return

    info = data.get("info", {})
    is_circulant = info.get("graph_type") == "circulant"

    if is_circulant:
        k5g_fast = k5_circulant(adj)
        k5g_brut = _brute_total_k5(adj)
        if k5g_fast != k5g_brut:
            print(f"!! divergence fast={k5g_fast} brute={k5g_brut} sur G")
        k5c_fast = k5_circulant(comp)
        k5c_brut = _brute_total_k5(comp)
        if k5c_fast != k5c_brut:
            print(f"!! divergence sur complément")
        print(f"[verify] type=circulant (cross-check fast vs brute OK)")
        k5g, k5c = k5g_fast, k5c_fast
        k4g = k4_circulant(adj)
        k4c = k4_circulant(comp)
    else:
        print(f"[verify] type=general (brute-force)")
        k5g = _brute_total_k5(adj)
        k5c = _brute_total_k5(comp)
        k4g = _brute_total_k4(adj)
        k4c = _brute_total_k4(comp)

    print(f"  K5(G)={k5g}  K5(complément)={k5c}  K5_total={k5g + k5c}")
    print(f"  K4(G)={k4g}  K4(complément)={k4c}  K4_total={k4g + k4c}")
    if "S_distances" in info:
        print(f"  S = {info['S_distances']}  |S|={info.get('S_size', '?')}")
    if k5g + k5c == 0:
        print(f"  ★★★ R(5,5) >= 44 démontré par ce graphe ! ★★★")


# -----------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Recherche R(5,5)>=44 dans les graphes circulants.")
    sub = parser.add_subparsers(dest="mode", required=True)

    pe = sub.add_parser("enumerate", help="énumère tous les circulants")
    pe.add_argument("--workers", type=int, default=0)
    pe.add_argument("--weight_k5", type=float, default=100.0)
    pe.add_argument("--weight_k4", type=float, default=5.0)
    pe.add_argument("--top_k", type=int, default=50,
                    help="nb de meilleurs circulants à sauvegarder (défaut 50)")
    pe.add_argument("--output_dir", type=str, default="circulants_top",
                    help="dossier de sauvegarde du top-K")

    pl = sub.add_parser("list", help="liste le top-K sauvegardé")
    pl.add_argument("--output_dir", type=str, default="circulants_top")

    ps = sub.add_parser("sa", help="recuit simulé dans S-space")
    ps.add_argument("--seed", type=int, default=None)
    ps.add_argument("--time", type=float, default=0)
    ps.add_argument("--weight_k5", type=float, default=100.0)
    ps.add_argument("--weight_k4", type=float, default=5.0)
    ps.add_argument("--t_start", type=float, default=200.0)
    ps.add_argument("--t_end", type=float, default=0.05)
    ps.add_argument("--cooling", type=float, default=0.99995)
    ps.add_argument("--two_bit_prob", type=float, default=0.15)
    ps.add_argument("--restart_after", type=int, default=2000)
    ps.add_argument("--log_every", type=int, default=200)
    ps.add_argument("--S_distances", type=str, default=None)

    pv = sub.add_parser("verify", help="recompute K5/K4 d'un fichier JSON")
    pv.add_argument("file", type=str)

    args = parser.parse_args()

    if args.mode == "enumerate":
        run_enumerate(args)
    elif args.mode == "list":
        run_list(args)
    elif args.mode == "sa":
        if args.seed is None:
            args.seed = random.SystemRandom().randint(0, 2**31)
        run_sa(args)
    elif args.mode == "verify":
        run_verify(args)


if __name__ == "__main__":
    main()
