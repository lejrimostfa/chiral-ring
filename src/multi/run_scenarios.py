"""Multi-ring scenarios and figures (Paper 2, work-order step 4).

M=2, N=64, R=10 disorder realizations, collision readout, n_traj=20k.
One figure = one function, fixed seed, regenerable in one command:

    python3 run_scenarios.py a     # control: disjoint baths, exact additivity
    python3 run_scenarios.py b     # DFS coherence + objectivity (C1)
    python3 run_scenarios.py c3    # concurrence C(t) vs D_rel(t)  (C3 figure)
    python3 run_scenarios.py c     # generic common bath: deficit vs multi-info (C2)
    python3 run_scenarios.py d     # sweep shared fraction f
    python3 run_scenarios.py e     # frustration: forward + mirror (C5)
    python3 run_scenarios.py f     # sigma statistics + fluctuation theorem (C4)
    python3 run_scenarios.py all   # everything

Numbers printed to stdout are the ones quoted in results_multi.md.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import model_multi as mm

FIG = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "paper2", "")
N = 64
R = 10
NTRAJ = 20_000
GMIN, GMAX = 0.1, 0.4


# ----------------------------------------------------------------------
# coupling builders
# ----------------------------------------------------------------------
def couplings_shared(Nq, f, rng, gmin=GMIN, gmax=GMAX):
    """Ring0/ring1 couplings with a shared fraction f. f=0 disjoint, f=1 all
    shared. Returns (g_mat, mirror_mask) where mirror_mask marks ring-1-private
    qubits (used only by the frustration scenario)."""
    n_sh = int(round(f * Nq))
    n_pr = Nq - n_sh
    n1 = n_pr // 2
    g = np.zeros((2, Nq))
    g[0, :n_sh] = rng.uniform(gmin, gmax, n_sh)         # shared: both couple
    g[1, :n_sh] = rng.uniform(gmin, gmax, n_sh)
    g[0, n_sh:n_sh + n1] = rng.uniform(gmin, gmax, n1)  # ring-0 private
    g[1, n_sh + n1:] = rng.uniform(gmin, gmax, n_pr - n1)  # ring-1 private
    mirror_mask = np.zeros(Nq, bool)
    mirror_mask[n_sh + n1:] = True
    return g, mirror_mask


def couplings_identical(Nq, rng, gmin=GMIN, gmax=GMAX):
    """DFS case: g_{0,n}=g_{1,n} for all n (all qubits shared, identical)."""
    g0 = rng.uniform(gmin, gmax, Nq)
    return np.vstack([g0, g0])


def couplings_random_shared(Nq, rng, gmin=GMIN, gmax=GMAX):
    """Generic case: both rings couple to all qubits, independent couplings."""
    return rng.uniform(gmin, gmax, (2, Nq))


def _avg_over_realizations(fn, seeds):
    """fn(rng)->1d array; return (mean, std) over seeds."""
    vals = np.array([fn(np.random.default_rng(s)) for s in seeds])
    return vals.mean(axis=0), vals.std(axis=0)


# ----------------------------------------------------------------------
# (a) control: disjoint baths -> exact additivity, zero correlation
# ----------------------------------------------------------------------
def scenario_a():
    ts = np.linspace(0.1, 2.0, 20)
    seeds = [1000 + r for r in range(R)]

    def one(rng):
        g, _ = couplings_shared(N, 0.0, rng)
        out = []
        for t in ts:
            d = mm.sample_joint(g, t, NTRAJ, np.random.default_rng(7))
            out.append([d["deficit"], d["cov"], d["tc"]])
        return np.array(out).ravel()

    m, s = _avg_over_realizations(one, seeds)
    m = m.reshape(len(ts), 3); s = s.reshape(len(ts), 3)
    deficit, cov, tc = m[:, 0], m[:, 1], m[:, 2]
    print("=== (a) CONTROL disjoint baths (exactness test) ===")
    print(f"  max|deficit| over t = {np.abs(deficit).max():.3e}  (should be ~0)")
    print(f"  max|cov|     over t = {np.abs(cov).max():.3e}  (should be ~0)")
    print(f"  max|tc|      over t = {np.abs(tc).max():.3e}  (should be ~0)")

    plt.figure(figsize=(6, 4))
    plt.axhline(0, color="k", lw=0.5)
    plt.plot(ts, deficit, "o-", label=r"deficit $\sigma_S-\sum_i\sigma_i$")
    plt.plot(ts, cov, "s-", label=r"posterior cov $s_1 s_2$")
    plt.plot(ts, tc, "^-", label="record multi-info (tc)")
    plt.xlabel("t"); plt.ylabel("nats")
    plt.title("(a) Disjoint baths: exact additivity, zero correlation")
    plt.legend(); plt.tight_layout()
    plt.savefig(FIG + "fig_multi_a_control.png", dpi=150); plt.close()


# ----------------------------------------------------------------------
# (b) DFS coherence + objectivity of the total current  (C1)
# ----------------------------------------------------------------------
def scenario_b():
    ts = np.linspace(0.0, 3.0, 60)
    rng = np.random.default_rng(2001)
    g_id = couplings_identical(N, rng)
    g_rd = couplings_random_shared(N, np.random.default_rng(2002))
    S_pm, S_mp = np.array([1, -1]), np.array([-1, 1])
    S_pp, S_mm = np.array([1, 1]), np.array([-1, -1])

    def curves(g):
        pm = np.array([abs(mm.coherence(g, S_pm, S_mp, t)) for t in ts])
        pp = np.array([abs(mm.coherence(g, S_pp, S_mm, t)) for t in ts])
        return pm, pp

    pm_id, pp_id = curves(g_id)
    pm_rd, pp_rd = curves(g_rd)
    print("=== (b) DFS coherence (C1) ===")
    print(f"  identical: |Gamma_(+-,-+)| min={pm_id.min():.4f} (protected~1), "
          f"|Gamma_(++,--)| min={pp_id.min():.4f} (decays)")
    print(f"  random:    |Gamma_(+-,-+)| min={pm_rd.min():.4f}, "
          f"|Gamma_(++,--)| min={pp_rd.min():.4f}")

    fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for a, (pm, pp, ttl) in zip(ax, [(pm_id, pp_id, "identical couplings (DFS)"),
                                     (pm_rd, pp_rd, "random couplings")]):
        a.plot(ts, pm, label=r"$|\Gamma_{(+-),(-+)}|$", lw=2)
        a.plot(ts, pp, label=r"$|\Gamma_{(++),(--)}|$", lw=2)
        a.set_xlabel("t"); a.set_title(ttl); a.legend()
    ax[0].set_ylabel(r"$|\Gamma|$")
    fig.suptitle("(b) DFS: zero-total-chirality coherence is protected")
    fig.tight_layout(); fig.savefig(FIG + "fig_multi_b_dfs_coherence.png", dpi=150)
    plt.close(fig)

    # objectivity (C1): contrast identical vs random couplings.
    # identical (DFS): individual P_i stays at chance ~0.5 (chirality hidden),
    #   but the TOTAL current objectifies where it is recordable (nonzero sector).
    # random: individual P_i objectifies normally.
    P1_id, Ptot_id, P1_rd = [], [], []
    for t in ts:
        di = mm.sample_joint(g_id, t, NTRAJ, np.random.default_rng(55))
        dr = mm.sample_joint(g_rd, t, NTRAJ, np.random.default_rng(56))
        P1_id.append(di["P_i"][0]); Ptot_id.append(di["pol_total_cond"])
        P1_rd.append(dr["P_i"][0])
    print(f"  objectivity at t={ts[-1]:.1f}: identical P_1={P1_id[-1]:.3f} "
          f"(chance), P_total|nonzero={Ptot_id[-1]:.3f} (objectifies); "
          f"random P_1={P1_rd[-1]:.3f}")

    plt.figure(figsize=(6, 4))
    plt.plot(ts, Ptot_id, lw=2, color="C2",
             label=r"$P_{\rm total}$ identical (recordable sector)")
    plt.plot(ts, P1_rd, "-", color="C1", label=r"$P_1$ random couplings")
    plt.plot(ts, P1_id, "--", color="C0", label=r"$P_1$ identical (DFS)")
    plt.axhline(0.5, color="k", lw=0.5, ls=":")
    plt.xlabel("t"); plt.ylabel("polarization")
    plt.title("(b) DFS objectivity: total current records; individuals don't")
    plt.legend(); plt.tight_layout()
    plt.savefig(FIG + "fig_multi_b_objectivity.png", dpi=150); plt.close()


# ----------------------------------------------------------------------
# (b/C3) concurrence C(t) vs relative distinguishability D_rel(t)
# ----------------------------------------------------------------------
def scenario_c3():
    ts = np.linspace(0.0, 3.0, 80)
    g_id = couplings_identical(N, np.random.default_rng(3001))
    g_rd = couplings_random_shared(N, np.random.default_rng(3002))
    S_pm, S_mp = np.array([1, -1]), np.array([-1, 1])

    def CD(g):
        C = np.array([mm.concurrence_pm(g, t) for t in ts])
        Dr = np.array([mm.relative_distinguishability(g, S_pm, S_mp, t) for t in ts])
        return C, Dr

    C_id, Dr_id = CD(g_id)
    C_rd, Dr_rd = CD(g_rd)
    print("=== (C3) concurrence vs relative distinguishability ===")
    print(f"  identical (DFS): D_rel max={Dr_id.max():.2e} (~0), "
          f"C min={C_id.min():.4f} (protected ~1)")
    print(f"  random: C at D_rel~1 nat, C(final)={C_rd[-1]:.4f}, "
          f"monotone decreasing in D_rel = "
          f"{np.all(np.diff(C_rd[np.argsort(Dr_rd)]) <= 1e-9)}")

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(ts, C_id, lw=2, label="C(t) identical (DFS)")
    ax[0].plot(ts, C_rd, lw=2, label="C(t) random")
    ax[0].set_xlabel("t"); ax[0].set_ylabel("concurrence C(t)")
    ax[0].legend(); ax[0].set_title("concurrence vs time")
    ax[1].plot(Dr_rd, C_rd, "o-", ms=3, label="random")
    ax[1].scatter(Dr_id, C_id, s=12, color="C0",
                  label="identical (D_rel$\\equiv$0)")
    ax[1].set_xlabel(r"relative distinguishability $D_{\rm rel}(t)$ [nats]")
    ax[1].set_ylabel("concurrence C(t)")
    ax[1].legend(); ax[1].set_title("C decreases as records form")
    fig.suptitle("(C3) Entanglement survives exactly where no records form")
    fig.tight_layout(); fig.savefig(FIG + "fig_multi_c3_concurrence.png", dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------------
# (c) generic common bath: deficit vs record multi-information  (C2)
# ----------------------------------------------------------------------
def scenario_c():
    ts = np.linspace(0.1, 2.0, 20)
    seeds = [4000 + r for r in range(R)]

    def one(rng):
        g = couplings_random_shared(N, rng)
        out = []
        for t in ts:
            d = mm.sample_joint(g, t, NTRAJ, np.random.default_rng(77))
            out.append([d["deficit"], d["tc"], d["cov"]])
        return np.array(out).ravel()

    m, s = _avg_over_realizations(one, seeds)
    m = m.reshape(len(ts), 3); s = s.reshape(len(ts), 3)
    deficit, tc, cov = m[:, 0], m[:, 1], m[:, 2]
    ed, et = s[:, 0], s[:, 1]
    M = 2
    ceiling = M * np.log(2)   # hard bound on any multi-information of M bits
    print("=== (c) generic common bath: additivity deficit (C2) ===")
    print(f"  deficit >= 0 everywhere: {np.all(deficit > -1e-6)}  (firm part of C2)")
    print(f"  deficit(final)={deficit[-1]:.3f}+/-{ed[-1]:.3f} nats  "
          f"vs multi-info ceiling M*ln2={ceiling:.3f} nats")
    print(f"  deficit/tc (median) = {np.nanmedian(deficit/np.where(tc>1e-6,tc,np.nan)):.1f}"
          f"  -> deficit is EXTENSIVE, tc is bounded: literal identity FAILS")
    print(f"  at t={ts[-1]:.1f}: tc={tc[-1]:.3f}+/-{et[-1]:.3f}, cov={cov[-1]:.3f}")

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].errorbar(ts, deficit, yerr=ed, fmt="o-", capsize=3,
                   label=r"deficit $\Delta=\sigma_S-\sum_i\sigma_i$")
    ax[0].plot(ts, tc, "s-", label="record multi-info (posterior TC)")
    ax[0].axhline(ceiling, color="C3", ls="--", label=r"$M\ln 2$ multi-info bound")
    ax[0].set_xlabel("t"); ax[0].set_ylabel("nats"); ax[0].legend()
    ax[0].set_title("deficit is extensive; multi-info is bounded")
    ax[1].plot(ts, cov, "^-", color="C3")
    ax[1].axhline(0, color="k", lw=0.5)
    ax[1].set_xlabel("t"); ax[1].set_ylabel(r"posterior cov $\langle s_1s_2\rangle$")
    ax[1].set_title("record-induced ring-ring correlation")
    fig.suptitle(r"(c) C2: $\Delta\geq0$ holds, but $\Delta$ is not a "
                 r"(bounded) multi-information")
    fig.tight_layout(); fig.savefig(FIG + "fig_multi_c_deficit.png", dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------------
# (d) sweep shared fraction f
# ----------------------------------------------------------------------
def scenario_d():
    fs = np.linspace(0.0, 1.0, 9)
    t = 1.0
    seeds = [5000 + r for r in range(R)]

    def one(rng):
        out = []
        for f in fs:
            g, _ = couplings_shared(N, f, rng)
            d = mm.sample_joint(g, t, NTRAJ, np.random.default_rng(88))
            out.append([d["deficit"], d["cov"], d["tc"]])
        return np.array(out).ravel()

    m, s = _avg_over_realizations(one, seeds)
    m = m.reshape(len(fs), 3); s = s.reshape(len(fs), 3)
    deficit, cov, tc = m[:, 0], m[:, 1], m[:, 2]
    ed, ec = s[:, 0], s[:, 1]
    print("=== (d) sweep shared fraction f (t=1.0) ===")
    for f, D, C, T in zip(fs, deficit, cov, tc):
        print(f"  f={f:.2f}: deficit={D:7.3f}  cov={C:+.4f}  tc={T:.3f}")

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].errorbar(fs, deficit, yerr=ed, fmt="o-", capsize=3, label="deficit")
    ax[0].errorbar(fs, tc, yerr=s[:, 2], fmt="s-", capsize=3, label="multi-info tc")
    ax[0].set_xlabel("shared fraction f"); ax[0].set_ylabel("nats"); ax[0].legend()
    ax[0].set_title("deficit / multi-info vs f")
    ax[1].errorbar(fs, cov, yerr=ec, fmt="^-", capsize=3, color="C3")
    ax[1].axhline(0, color="k", lw=0.5)
    ax[1].set_xlabel("shared fraction f")
    ax[1].set_ylabel(r"posterior cov $\langle s_1 s_2\rangle$")
    ax[1].set_title("record-induced ring-ring correlation vs f")
    fig.suptitle("(d) Bath sharing turns on correlation and the deficit")
    fig.tight_layout(); fig.savefig(FIG + "fig_multi_d_sweep_f.png", dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------------
# (e) frustration: ring0 forward + ring1 mirror, sharing fraction f  (C5)
# ----------------------------------------------------------------------
def scenario_e():
    ts = np.linspace(0.05, 1.5, 30)
    T = 1.5
    fvals = [0.0, 0.25, 0.5, 1.0]
    seeds = [6000 + r for r in range(R)]

    # sigma_1(t), sigma_2(t) averaged over realizations, per f
    data = {}
    for f in fvals:
        def one(rng, f=f):
            g, mask = couplings_shared(N, f, rng)
            rec = []
            for t in ts:
                d = mm.sample_frustration(g, t, T, mask, NTRAJ,
                                          np.random.default_rng(99))
                rec.append([d["sigma"][0], d["sigma"][1],
                            d["P_i"][0], d["P_i"][1]])
            return np.array(rec).ravel()
        m, s = _avg_over_realizations(one, seeds)
        data[f] = m.reshape(len(ts), 4)

    print("=== (e) FRUSTRATION forward(ring1) + mirror(ring2) (C5) ===")
    for f in fvals:
        s1, s2 = data[f][:, 0], data[f][:, 1]
        # crossing of ring2's arrow (sign change of sigma_2 in time)
        sign_change = np.where(np.diff(np.sign(s2)) != 0)[0]
        tc_txt = (f"sigma_2 crosses 0 at t~{ts[sign_change[0]]:.2f}"
                  if len(sign_change) else "sigma_2 no zero-crossing in t")
        print(f"  f={f:.2f}: sigma_1(final)={s1[-1]:+.2f} "
              f"sigma_2(final)={s2[-1]:+.2f}  [{tc_txt}]")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for f in fvals:
        ax[0].plot(ts, data[f][:, 0], label=f"f={f:.2f}")
        ax[1].plot(ts, data[f][:, 1], label=f"f={f:.2f}")
    for a, ttl in zip(ax, [r"$\sigma_1$ (ring 1, forward)",
                           r"$\sigma_2$ (ring 2, mirror)"]):
        a.axhline(0, color="k", lw=0.5); a.set_xlabel("t")
        a.set_ylabel("nats"); a.set_title(ttl); a.legend()
    fig.suptitle("(e) C5 arrow frustration: shared qubits pull ring 2's arrow up")
    fig.tight_layout(); fig.savefig(FIG + "fig_multi_e_frustration.png", dpi=150)
    plt.close(fig)

    # sigma_2 vs f at an INTERMEDIATE time (mirror still has records to consume;
    # at t=T the mirror is fully consumed and sigma_2->0 for every f, which would
    # make the crossing degenerate).
    t_star = 0.5
    f_dense = np.linspace(0, 1, 11)
    s2_f = []
    for f in f_dense:
        vals = []
        for r in range(R):
            g, mask = couplings_shared(N, f, np.random.default_rng(6000 + r))
            d = mm.sample_frustration(g, t_star, T, mask, NTRAJ,
                                      np.random.default_rng(99))
            vals.append(d["sigma"][1])
        s2_f.append(np.mean(vals))
    s2_f = np.array(s2_f)
    cross = np.where(np.diff(np.sign(s2_f)) != 0)[0]
    fstar = (f_dense[cross[0]] + (0 - s2_f[cross[0]]) *
             (f_dense[cross[0] + 1] - f_dense[cross[0]]) /
             (s2_f[cross[0] + 1] - s2_f[cross[0]])) if len(cross) else np.nan
    print(f"  ring-2 arrow annulation at f* = {fstar:.3f} (t={t_star:.2f})")
    plt.figure(figsize=(6, 4))
    plt.axhline(0, color="k", lw=0.5)
    plt.plot(f_dense, s2_f, "o-")
    if not np.isnan(fstar):
        plt.axvline(fstar, color="C3", ls="--", label=f"f* = {fstar:.2f}")
        plt.legend()
    plt.xlabel("shared fraction f")
    plt.ylabel(r"$\sigma_2$ (mirror ring) [nats]")
    plt.title(f"(e) Mirror arrow flips sign as sharing increases (t={t_star})")
    plt.tight_layout(); plt.savefig(FIG + "fig_multi_e_crossing.png", dpi=150)
    plt.close()


# ----------------------------------------------------------------------
# (f) sigma statistics + detailed fluctuation theorem  (C4)
# ----------------------------------------------------------------------
def scenario_f():
    g = couplings_random_shared(N, np.random.default_rng(7001))
    # moderate-D times: the detailed FT needs both sigma and -sigma sampled, and
    # the sigma<0 tail is exponentially suppressed, so very large D is untestable
    # (that suppression is itself the C4 content, shown by frac(sigma<0)).
    times = [0.1, 0.2, 0.35]
    print("=== (f) sigma statistics + fluctuation theorem (C4) ===")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    NT = 400_000
    for t in times:
        d = mm.sample_joint(g, t, NT, np.random.default_rng(321))
        sig = d["sigma_S_traj"]
        frac_neg = float(np.mean(sig < 0))
        # symmetric bins around 0 so bin k and its mirror -k align exactly
        lim = np.percentile(np.abs(sig), 99)
        edges = np.linspace(-lim, lim, 61)
        counts, _ = np.histogram(sig, bins=edges)
        centers = 0.5 * (edges[:-1] + edges[1:])
        dens = counts / (counts.sum() * (edges[1] - edges[0]))
        ax[0].plot(centers, dens, label=f"t={t}")
        # detailed FT: ln[p(sigma)/p(-sigma)] vs sigma, expect slope 1.
        # Restrict to bins with adequate counts on BOTH sides (the sigma<0 tail
        # is sparse; including under-sampled bins biases the slope low).
        pos = centers > 0
        cp = counts[pos]; cn = counts[::-1][pos]
        good = (cp >= 25) & (cn >= 25)
        x = centers[pos][good]; y = np.log(cp[good] / cn[good])
        slope = np.polyfit(x, y, 1)[0] if len(x) > 3 else np.nan
        ax[1].plot(x, y, "o", ms=3, label=f"t={t} (slope={slope:.2f})")
        print(f"  t={t}: <sigma>={sig.mean():.3f}  D={d['D']:.3f}  "
              f"frac(sigma<0)={frac_neg:.4f}  FT slope(well-sampled)={slope:.3f}")
    xx = np.linspace(0, np.percentile(np.abs(sig), 95), 10)
    ax[1].plot(xx, xx, "k:", label="slope 1 (FT)")
    ax[0].set_xlabel(r"$\sigma_S$"); ax[0].set_ylabel("p($\\sigma$)")
    ax[0].set_yscale("log"); ax[0].legend(); ax[0].set_title("p($\\sigma$) at 3 times")
    ax[1].set_xlabel(r"$\sigma$"); ax[1].set_ylabel(r"$\ln[p(\sigma)/p(-\sigma)]$")
    ax[1].legend(); ax[1].set_title("detailed FT: slope 1")
    fig.suptitle("(f) C4: p($\\sigma$)/p($-\\sigma$) = e^$\\sigma$")
    fig.tight_layout(); fig.savefig(FIG + "fig_multi_f_fluctuation.png", dpi=150)
    plt.close(fig)


SCENARIOS = {"a": scenario_a, "b": scenario_b, "c3": scenario_c3,
             "c": scenario_c, "d": scenario_d, "e": scenario_e, "f": scenario_f}

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    keys = list(SCENARIOS) if which == "all" else [which]
    for k in keys:
        SCENARIOS[k]()
    print("\nDone:", ", ".join(keys))
