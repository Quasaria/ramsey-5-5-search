#!/usr/bin/env python3
"""
ramsey_55_search.py — Recherche heuristique pour R(5,5) >= 44 (v5)
==================================================================

Cherche un graphe de N=43 sommets sans K5 monochrome (ni clique de 5,
ni stable de 5). Trouver un tel graphe prouverait R(5,5) >= 44.

Méthode : recuit simulé sur représentation bitset.
- adj[v] : bitmask des voisins de v.
- coût  : w5 * #K5(G+complément) + w4 * #K4(G+complément).
- évaluation locale : seuls les K5/K4 contenant l'arête flippée changent.

Évolutions :
- v2 : mutation ciblée K4.
- v3 : moves à 2 arêtes (second flip aléatoire).
- v4 : moves à 2 arêtes avec lookahead intelligent + sauvegarde robuste
       Windows + --resume + weight_k4 défaut 5.
- v5 :
  * Tabu list (fenêtre glissante des dernières arêtes flippées).
  * Vertex-block moves (flip de k arêtes incidentes à un même sommet).
  * ILS restart : sur stagnation, on perturbe le best, pas le current.
  * Vérification d'intégrité périodique (compteurs vs recompute).
  * --start_from FILE : démarrage depuis un graphe arbitraire.
  * Compteurs d'améliorations par type de move (1f / 2f / vertex).

Sauvegardes :
- ramsey55_best.{json,matrix.txt}       : meilleur graphe trouvé.
- ramsey55_checkpoint.{json,matrix.txt} : état courant (périodique).

Lancements typiques :
    # depuis zéro
    python ramsey_55_search.py --seed 42 --time 3600
    # reprise du meilleur, exploration plus agressive
    python ramsey_55_search.py --resume --time 7200 --t_start 30 \
        --vertex_move_prob 0.10 --two_flip_prob 0.35
"""
import sys
import os
import time
import math
import random
import json
import argparse
import collections

# Force UTF-8 sur stdout/stderr pour eviter UnicodeEncodeError sous Windows
# (notamment quand on redirige vers un fichier via PowerShell qui passe en cp1252)
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except AttributeError:
    pass  # Python < 3.7, peu probable

# -------------------------------------------------------------
# Constantes du problème
# -------------------------------------------------------------
N    = 43
MASK = (1 << N) - 1

if sys.version_info >= (3, 10):
    popcount = int.bit_count
else:
    def popcount(x):  # noqa: E301
        return bin(x).count("1")

K4_EDGE_PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))

# -------------------------------------------------------------
# Hyperparamètres par défaut
# -------------------------------------------------------------
DEFAULTS = dict(
    weight_k5            = 100.0,
    weight_k4            =   5.0,
    t_start              =  80.0,
    t_end                =   0.05,
    cooling              =   0.9999985,
    restart_after        = 250_000,
    perturb_size         =      60,
    targeted_prob        =     0.5,
    targeted_tries       =       8,
    two_flip_prob        =     0.3,
    two_flip_lookahead   =       4,
    vertex_move_prob     =    0.05,   # NEW v5
    vertex_move_kmin     =       3,
    vertex_move_kmax     =       5,
    tabu_size            =      30,   # NEW v5
    integrity_every      = 500_000,   # NEW v5
    log_every            =   5_000,
    checkpoint_every     = 100_000,
)

# -------------------------------------------------------------
# Helpers bitset
# -------------------------------------------------------------
def higher_mask(w):
    return MASK & ~((1 << (w + 1)) - 1)

def random_graph(p, rng):
    adj = [0] * N
    for i in range(N):
        for j in range(i + 1, N):
            if rng.random() < p:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
    return adj

def complement_of(adj):
    return [(MASK ^ adj[v]) & ~(1 << v) for v in range(N)]

def random_bit(mask, rng):
    """Index uniformément aléatoire d'un bit à 1, ou -1 si mask vide."""
    k = popcount(mask)
    if k == 0:
        return -1
    r = rng.randrange(k)
    while r > 0:
        mask &= mask - 1
        r -= 1
    return (mask & -mask).bit_length() - 1

def flip_edge(u, v, adj, comp):
    """Bascule en place l'arête (u,v) dans G et son complément."""
    bv, bu = 1 << v, 1 << u
    adj[u]  ^= bv;  adj[v]  ^= bu
    comp[u] ^= bv;  comp[v] ^= bu

# -------------------------------------------------------------
# Comptages dans un sous-graphe induit
# -------------------------------------------------------------
def edges_in(S, adj):
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

def k4s_in(S, adj):
    cnt = 0
    bits = S
    while bits:
        low = bits & -bits
        w = low.bit_length() - 1
        bits ^= low
        cnt += triangles_in(adj[w] & S & higher_mask(w), adj)
    return cnt

def total_k5(adj):
    cnt = 0
    for u in range(N):
        cnt += k4s_in(adj[u] & higher_mask(u), adj)
    return cnt

def total_k4(adj):
    cnt = 0
    for u in range(N):
        cnt += triangles_in(adj[u] & higher_mask(u), adj)
    return cnt

# -------------------------------------------------------------
# Évaluation delta lors d'une bascule d'arête
# -------------------------------------------------------------
def delta_counts(u, v, is_edge, adj, comp):
    common_in  = adj[u]  & adj[v]
    common_out = comp[u] & comp[v]
    T_G = triangles_in(common_in,  adj)
    T_C = triangles_in(common_out, comp)
    E_G = edges_in    (common_in,  adj)
    E_C = edges_in    (common_out, comp)
    if is_edge:
        return (-T_G + T_C, -E_G + E_C)
    else:
        return ( T_G - T_C,  E_G - E_C)

# -------------------------------------------------------------
# Tabu list (set + deque pour O(1) lookup et FIFO eviction)
# -------------------------------------------------------------
class TabuList:
    def __init__(self, capacity):
        self.capacity = capacity
        self.queue = collections.deque()
        self.set = set()

    def __contains__(self, edge):
        return edge in self.set

    def push(self, edge):
        if edge in self.set:
            return
        self.queue.append(edge)
        self.set.add(edge)
        while len(self.queue) > self.capacity:
            old = self.queue.popleft()
            self.set.discard(old)

    def clear(self):
        self.queue.clear()
        self.set.clear()

# -------------------------------------------------------------
# Sélection d'arête (ciblée K4 ou uniforme), avec filtre tabu
# -------------------------------------------------------------
def sample_k4_edge(target_adj, rng, max_tries):
    for _ in range(max_tries):
        u = rng.randrange(N)
        Nu = target_adj[u]
        if popcount(Nu) < 3:
            continue
        v = random_bit(Nu, rng)
        Nuv = Nu & target_adj[v]
        if popcount(Nuv) < 2:
            continue
        w = random_bit(Nuv, rng)
        Nuvw = Nuv & target_adj[w]
        if popcount(Nuvw) < 1:
            continue
        x = random_bit(Nuvw, rng)
        verts = (u, v, w, x)
        i, j = K4_EDGE_PAIRS[rng.randrange(6)]
        a, b = verts[i], verts[j]
        if a > b:
            a, b = b, a
        return (a, b)
    return None

def select_edge(adj, comp, targeted_prob, max_tries, rng, tabu=None):
    """
    Tire (u, v, was_targeted_hit). Si tabu fourni, jusqu'à 4 retries
    pour éviter une arête tabu, puis fallback non-filtré pour garantir
    un retour.
    """
    for _ in range(4):
        if targeted_prob > 0.0 and rng.random() < targeted_prob:
            target = adj if rng.random() < 0.5 else comp
            picked = sample_k4_edge(target, rng, max_tries)
            if picked is not None:
                if tabu is None or picked not in tabu:
                    return picked[0], picked[1], True
                continue
        u = rng.randrange(N)
        v = rng.randrange(N - 1)
        if v >= u:
            v += 1
        if u > v:
            u, v = v, u
        if tabu is None or (u, v) not in tabu:
            return u, v, False
    # fallback : on rend sans filtrer
    u = rng.randrange(N)
    v = rng.randrange(N - 1)
    if v >= u:
        v += 1
    if u > v:
        u, v = v, u
    return u, v, False

# -------------------------------------------------------------
# Vertex-block move : flip de k arêtes incidentes à un sommet
# -------------------------------------------------------------
def vertex_move(adj, comp, weight_k5, weight_k4, T, rng, k_min, k_max, tabu):
    """
    Pick u ∈ V, pick k voisins distincts v1..vk parmi V\\{u}, flip les
    k arêtes (u,vi), Metropolis sur la somme des deltas. Retourne
    (delta_k5, delta_k4, accepted, applied_edges).
    """
    u = rng.randrange(N)
    others = list(range(N))
    others.pop(u)
    k = rng.randint(k_min, min(k_max, len(others)))
    chosen = rng.sample(others, k)

    applied = []
    tot_k5 = 0
    tot_k4 = 0
    for v in chosen:
        a, b = (u, v) if u < v else (v, u)
        if tabu is not None and (a, b) in tabu:
            continue  # skip silencieusement
        is_edge = (adj[a] >> b) & 1
        d_k5, d_k4 = delta_counts(a, b, is_edge, adj, comp)
        flip_edge(a, b, adj, comp)
        tot_k5 += d_k5
        tot_k4 += d_k4
        applied.append((a, b))

    if not applied:
        return 0, 0, False, []

    d_cost = weight_k5 * tot_k5 + weight_k4 * tot_k4

    if d_cost <= 0 or rng.random() < math.exp(-d_cost / T):
        return tot_k5, tot_k4, True, applied
    else:
        # rejet : undo en ordre inverse
        for a, b in reversed(applied):
            flip_edge(a, b, adj, comp)
        return 0, 0, False, []

# -------------------------------------------------------------
# Sauvegarde robuste (anti-AV / anti-indexation Windows)
# -------------------------------------------------------------
_save_counter = [0]

def save_graph(adj, info, prefix):
    """
    Écriture qui ne jette JAMAIS d'exception au-dehors :
      - .tmp unique (PID + compteur) pour éviter les conflits AV.
      - 6 retries avec backoff exponentiel.
      - cleanup systématique.
    """
    _save_counter[0] += 1
    counter = _save_counter[0]
    pid = os.getpid()

    fname_json = prefix + ".json"
    fname_mat  = prefix + ".matrix.txt"
    tmp = f"{fname_json}.{pid}.{counter}.tmp"

    rows = []
    payload = None
    try:
        for i in range(N):
            rows.append("".join("1" if (adj[i] >> j) & 1 else "0"
                                for j in range(N)))
        payload = {
            "N": N,
            "info": info,
            "adj_bitmasks": [int(x) for x in adj],
            "adjacency_matrix": rows,
        }
    except Exception:
        return

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

    try:
        with open(fname_mat, "w") as f:
            f.write("\n".join(rows) + "\n")
    except (PermissionError, OSError):
        pass

def load_graph(filename):
    with open(filename, "r") as f:
        data = json.load(f)
    if data.get("N") != N:
        raise ValueError(f"N mismatch in {filename}: got {data.get('N')}, expected {N}")
    return [int(x) for x in data["adj_bitmasks"]]

# -------------------------------------------------------------
# Boucle principale
# -------------------------------------------------------------
def run(args):
    rng = random.Random(args.seed)

    # Initialisation
    adj = None
    if args.start_from:
        try:
            adj = load_graph(args.start_from)
            print(f"[start_from] graphe chargé depuis {args.start_from}")
        except (OSError, ValueError, json.JSONDecodeError) as e:
            print(f"[start_from] échec ({e})")
            adj = None
    if adj is None and args.resume:
        try:
            adj = load_graph("ramsey55_best.json")
            print(f"[resume] graphe chargé depuis ramsey55_best.json")
        except (OSError, ValueError, json.JSONDecodeError) as e:
            print(f"[resume] échec ({e}), fallback aléatoire")
            adj = None
    if adj is None:
        adj = random_graph(0.5, rng)

    comp = complement_of(adj)
    k5 = total_k5(adj) + total_k5(comp)
    k4 = total_k4(adj) + total_k4(comp)
    cost = args.weight_k5 * k5 + args.weight_k4 * k4

    best_cost = cost
    best_k5   = k5
    best_adj  = list(adj)

    tabu = TabuList(args.tabu_size) if args.tabu_size > 0 else None

    print(f"[init v5] N={N}  seed={args.seed}  K5={k5}  K4={k4}  cost={cost:.1f}")
    print(f"          w4={args.weight_k4}  tgt={args.targeted_prob}  "
          f"2f={args.two_flip_prob}(la={args.two_flip_lookahead})  "
          f"vm={args.vertex_move_prob}(k={args.vertex_move_kmin}-{args.vertex_move_kmax})  "
          f"tabu={args.tabu_size}  ils_from_best={args.ils_from_best}")

    save_graph(best_adj,
               {"iter": 0, "k5": k5, "k4": k4, "cost": cost, "tag": "initial"},
               prefix="ramsey55_best")

    T = args.t_start
    iters         = 0
    accepts       = 0

    one_flip_used = 0
    one_flip_acc  = 0
    two_flip_used = 0
    two_flip_acc  = 0
    vertex_used   = 0
    vertex_acc    = 0
    targeted_used = 0
    targeted_acc  = 0

    improv_by_one = 0
    improv_by_two = 0
    improv_by_vtx = 0
    last_improve  = 0

    deadline = (time.time() + args.time) if args.time > 0 else float("inf")
    t0 = time.time()

    try:
        while T > args.t_end and time.time() < deadline and best_k5 > 0:

            roll = rng.random()
            improved = False

            # ============= VERTEX-BLOCK MOVE =============
            if roll < args.vertex_move_prob:
                vertex_used += 1
                d_k5, d_k4, accepted, applied = vertex_move(
                    adj, comp, args.weight_k5, args.weight_k4,
                    T, rng, args.vertex_move_kmin, args.vertex_move_kmax, tabu)
                if accepted:
                    k5 += d_k5
                    k4 += d_k4
                    cost += args.weight_k5 * d_k5 + args.weight_k4 * d_k4
                    accepts += 1
                    vertex_acc += 1
                    if tabu is not None:
                        for ab in applied:
                            tabu.push(ab)
                    if cost < best_cost - 1e-9:
                        best_cost = cost
                        best_k5 = k5
                        best_adj = list(adj)
                        last_improve = iters
                        improv_by_vtx += 1
                        improved = True
                        save_graph(best_adj, {
                            "iter": iters, "k5": k5, "k4": k4, "cost": cost,
                            "T": T, "elapsed_sec": time.time() - t0,
                            "tag": "improvement", "move": "vertex",
                            "k_flips": len(applied),
                        }, prefix="ramsey55_best")
                        if k5 == 0:
                            print(f"\n!!! GRAPHE TROUVÉ (vertex) à it={iters} !!!")
                            return best_adj

            # ============= 2-FLIP AVEC LOOKAHEAD =============
            elif roll < args.vertex_move_prob + args.two_flip_prob:
                two_flip_used += 1
                u1, v1, hit1 = select_edge(adj, comp, args.targeted_prob,
                                           args.targeted_tries, rng, tabu)
                is_edge1 = (adj[u1] >> v1) & 1
                d1_k5, d1_k4 = delta_counts(u1, v1, is_edge1, adj, comp)
                flip_edge(u1, v1, adj, comp)

                best_d2_cost = math.inf
                best_d2_k5 = 0
                best_d2_k4 = 0
                best_uv2 = None
                hit2 = False
                for _ in range(args.two_flip_lookahead):
                    cu2, cv2, chit2 = select_edge(adj, comp, args.targeted_prob,
                                                  args.targeted_tries, rng, tabu)
                    if (cu2, cv2) == (u1, v1):
                        continue
                    cis_edge2 = (adj[cu2] >> cv2) & 1
                    cd2_k5, cd2_k4 = delta_counts(cu2, cv2, cis_edge2, adj, comp)
                    cd_cost = args.weight_k5 * cd2_k5 + args.weight_k4 * cd2_k4
                    if cd_cost < best_d2_cost:
                        best_d2_cost = cd_cost
                        best_d2_k5 = cd2_k5
                        best_d2_k4 = cd2_k4
                        best_uv2 = (cu2, cv2)
                        hit2 = chit2

                hit_used = hit1 or hit2
                if hit_used:
                    targeted_used += 1

                if best_uv2 is None:
                    flip_edge(u1, v1, adj, comp)  # défaire le flip 1
                else:
                    tot_k5 = d1_k5 + best_d2_k5
                    tot_k4 = d1_k4 + best_d2_k4
                    d_cost = args.weight_k5 * tot_k5 + args.weight_k4 * tot_k4

                    if d_cost <= 0 or rng.random() < math.exp(-d_cost / T):
                        flip_edge(*best_uv2, adj, comp)
                        k5 += tot_k5
                        k4 += tot_k4
                        cost += d_cost
                        accepts += 1
                        two_flip_acc += 1
                        if hit_used:
                            targeted_acc += 1
                        if tabu is not None:
                            tabu.push((u1, v1))
                            tabu.push(best_uv2)
                        if cost < best_cost - 1e-9:
                            best_cost = cost
                            best_k5 = k5
                            best_adj = list(adj)
                            last_improve = iters
                            improv_by_two += 1
                            improved = True
                            save_graph(best_adj, {
                                "iter": iters, "k5": k5, "k4": k4, "cost": cost,
                                "T": T, "elapsed_sec": time.time() - t0,
                                "tag": "improvement", "move": "two_flip",
                            }, prefix="ramsey55_best")
                            if k5 == 0:
                                print(f"\n!!! GRAPHE TROUVÉ (2-flip) à it={iters} !!!")
                                return best_adj
                    else:
                        flip_edge(u1, v1, adj, comp)  # rejet → undo

            # ============= 1-FLIP =============
            else:
                one_flip_used += 1
                u, v, hit = select_edge(adj, comp, args.targeted_prob,
                                        args.targeted_tries, rng, tabu)
                if hit:
                    targeted_used += 1
                is_edge = (adj[u] >> v) & 1
                d_k5, d_k4 = delta_counts(u, v, is_edge, adj, comp)
                d_cost = args.weight_k5 * d_k5 + args.weight_k4 * d_k4

                if d_cost <= 0 or rng.random() < math.exp(-d_cost / T):
                    flip_edge(u, v, adj, comp)
                    k5 += d_k5
                    k4 += d_k4
                    cost += d_cost
                    accepts += 1
                    one_flip_acc += 1
                    if hit:
                        targeted_acc += 1
                    if tabu is not None:
                        tabu.push((u, v))

                    if cost < best_cost - 1e-9:
                        best_cost = cost
                        best_k5 = k5
                        best_adj = list(adj)
                        last_improve = iters
                        improv_by_one += 1
                        improved = True
                        save_graph(best_adj, {
                            "iter": iters, "k5": k5, "k4": k4, "cost": cost,
                            "T": T, "elapsed_sec": time.time() - t0,
                            "tag": "improvement", "move": "one_flip",
                        }, prefix="ramsey55_best")
                        if k5 == 0:
                            print(f"\n!!! GRAPHE TROUVÉ (1-flip) à it={iters} !!!")
                            return best_adj

            T *= args.cooling
            iters += 1

            # ===================== ILS RESTART =====================
            if iters - last_improve > args.restart_after:
                if args.ils_from_best:
                    adj = list(best_adj)
                    comp = complement_of(adj)
                if tabu is not None:
                    tabu.clear()
                T = max(T, args.t_start * 0.5)
                for _ in range(args.perturb_size):
                    a = rng.randrange(N)
                    b = rng.randrange(N)
                    if a == b:
                        continue
                    if a > b:
                        a, b = b, a
                    flip_edge(a, b, adj, comp)
                k5 = total_k5(adj) + total_k5(comp)
                k4 = total_k4(adj) + total_k4(comp)
                cost = args.weight_k5 * k5 + args.weight_k4 * k4
                last_improve = iters
                source = "best" if args.ils_from_best else "current"
                print(f"[reheat @ {iters} from {source}] T←{T:.3f}  "
                      f"K5={k5}  K4={k4}  cost={cost:.1f}")

            # ===================== INTEGRITY CHECK =====================
            if iters % args.integrity_every == 0 and iters > 0:
                real_k5 = total_k5(adj) + total_k5(comp)
                real_k4 = total_k4(adj) + total_k4(comp)
                if real_k5 != k5 or real_k4 != k4:
                    print(f"[DRIFT @ {iters}] running k5={k5}/k4={k4}, "
                          f"real k5={real_k5}/k4={real_k4} — corrigé.")
                    k5 = real_k5
                    k4 = real_k4
                    cost = args.weight_k5 * k5 + args.weight_k4 * k4
                else:
                    print(f"[integrity @ {iters}] OK (k5={k5} k4={k4})")

            # ===================== LOG =====================
            if iters % args.log_every == 0:
                el = time.time() - t0
                rate = iters / el if el > 0 else 0
                two_pct = 100 * two_flip_used / max(iters, 1)
                two_acc = 100 * two_flip_acc  / max(two_flip_used, 1)
                vm_pct  = 100 * vertex_used   / max(iters, 1)
                vm_acc  = 100 * vertex_acc    / max(vertex_used, 1)
                tgt_pct = 100 * targeted_used / max(iters, 1)
                tgt_acc = 100 * targeted_acc  / max(targeted_used, 1)
                print(f"it={iters:>9}  T={T:6.2f}  K5={k5:>4}  K4={k4:>5}  "
                      f"cost={cost:>8.0f}  best={best_k5:>4}  "
                      f"acc={100*accepts/max(iters,1):4.1f}%  "
                      f"2f={two_pct:4.1f}%({two_acc:4.1f}%) "
                      f"vm={vm_pct:4.1f}%({vm_acc:4.1f}%) "
                      f"tgt={tgt_pct:4.1f}%({tgt_acc:4.1f}%)  "
                      f"impr=1f:{improv_by_one}/2f:{improv_by_two}/vm:{improv_by_vtx}  "
                      f"{rate:>5.0f} it/s")

            # ===================== CHECKPOINT =====================
            if iters % args.checkpoint_every == 0 and iters > 0:
                save_graph(adj, {"iter": iters, "k5": k5, "k4": k4,
                                 "cost": cost, "T": T, "tag": "checkpoint"},
                           prefix="ramsey55_checkpoint")

    except KeyboardInterrupt:
        print("\n[interrupt] sauvegarde finale...")
        save_graph(adj, {"iter": iters, "k5": k5, "k4": k4, "cost": cost,
                         "T": T, "tag": "interrupted"},
                   prefix="ramsey55_checkpoint")

    el = time.time() - t0
    print(f"\n[stop] iter={iters}  best_K5={best_k5}  best_cost={best_cost:.1f}  "
          f"({el:.0f}s)")
    print(f"       improvements found by:  1-flip={improv_by_one}  "
          f"2-flip={improv_by_two}  vertex={improv_by_vtx}")
    print(f"       graphe : ramsey55_best.{{json,matrix.txt}}")
    return best_adj

# -------------------------------------------------------------
# CLI
# -------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Recherche heuristique R(5,5)>=44 — recuit simulé v5 "
                    "(tabu + vertex-move + ILS-restart + integrity check).")
    for k, v in DEFAULTS.items():
        if isinstance(v, int):
            p.add_argument(f"--{k}", type=int, default=v)
        else:
            p.add_argument(f"--{k}", type=float, default=v)
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--time", type=float, default=0,
                   help="durée maximale en secondes (0 = illimité)")
    p.add_argument("--resume", action="store_true",
                   help="reprendre depuis ramsey55_best.json si dispo")
    p.add_argument("--start_from", type=str, default=None,
                   help="charger un graphe depuis un fichier JSON donné")
    p.add_argument("--ils_from_best", action="store_true", default=True,
                   help="ILS : sur reheat, repartir du meilleur (défaut on)")
    p.add_argument("--no_ils_from_best", action="store_false",
                   dest="ils_from_best",
                   help="désactive ILS : reheat depuis l'état courant")
    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()
    if args.seed is None:
        args.seed = random.SystemRandom().randint(0, 2**31)
    run(args)