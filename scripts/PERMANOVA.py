#!/usr/bin/env python3
"""
PERMANOVA-style partial explained variance with a block factor ("study"),
excluding "response" from the tested variables, across three omics tables
(taxonomy, BGC, pathway). Outputs:
  - CSV of partial R² (% of total)
  - CSV of PERMANOVA p-values (restricted permutations within 'study')
  - Grouped bar plot with "*" for p<0.05, figsize=(7.2, 4), dpi=300

Distances: Bray–Curtis on row-normalized tables
Partial R²: adonis-like via Gower-centered distances and projection matrices
"""

import argparse, math
from typing import Dict, List, Tuple, Optional
import numpy as np, pandas as pd, matplotlib.pyplot as plt
try:
    from scipy.spatial.distance import pdist, squareform
    SCIPY_AVAILABLE = True
except Exception:
    SCIPY_AVAILABLE = False

# ---- Plasmid theme (colors & typography) ------------------------------------
PLASMID_COLORS = {
    "taxa":    "#1F77B4",  # blue
    "bgc":     "#D62728",  # dark red
    "pathway": "#87CEEB",  # light blue
}
PLASMID_FONT_FAMILY = "DejaVu Sans"   # portable sans-serif family
PLASMID_FONT_WEIGHT = "bold"          # label/tick/star emphasis


def read_metadata(meta_path: str, id_candidates: Optional[List[str]] = None) -> pd.DataFrame:
    if id_candidates is None:
        id_candidates = ["sample","sample_id","Sample","SampleID","sampleID",
                         "subject","Subject","id","ID"]
    md = pd.read_csv(meta_path, low_memory=False)
    use_index = next((c for c in md.columns if c in id_candidates), md.columns[0])
    md = md.set_index(use_index)
    md.index = md.index.astype(str).str.strip()
    return md


def read_omics_table(path: str, metadata_index: pd.Index) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=0, index_col=0)
    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()

    # If samples are columns, transpose
    overlap_rows = len(set(df.index) & set(metadata_index))
    overlap_cols = len(set(df.columns) & set(metadata_index))
    if overlap_cols > overlap_rows:
        df = df.transpose()
        df.index = df.index.astype(str).str.strip()
        df.columns = df.columns.astype(str).str.strip()

    # Numeric + drop all-zero features
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    nz = df.columns[(df.sum(axis=0) != 0)]
    if len(nz) > 0:
        df = df.loc[:, nz]
    return df


def align_samples(md: pd.DataFrame, tables: Dict[str, pd.DataFrame]):
    common = set(md.index)
    for _, tab in tables.items():
        common &= set(tab.index)
    common = sorted(list(common))
    return md.loc[common].copy(), {k: v.loc[common].copy() for k, v in tables.items()}


def compute_bray_curtis_matrix(data: np.ndarray) -> np.ndarray:
    if SCIPY_AVAILABLE:
        return squareform(pdist(data, metric="braycurtis"))
    n = data.shape[0]
    D = np.zeros((n, n), float)
    for i in range(n):
        diffs = np.abs(data[i+1:] - data[i])
        sums  = (data[i+1:] + data[i])
        num = diffs.sum(axis=1)
        den = sums.sum(axis=1)
        vals = np.divide(num, den, out=np.zeros_like(num), where=den != 0)
        D[i, i] = 0.0
        D[i, i+1:] = vals
        D[i+1:, i] = vals
    return D


def gower_center(D: np.ndarray) -> np.ndarray:
    A = -0.5 * (D**2)
    n = A.shape[0]
    H = np.eye(n) - np.ones((n, n)) / n
    return H @ A @ H


def design_matrix_for_variable(s: pd.Series) -> Optional[np.ndarray]:
    sn = pd.to_numeric(s, errors="coerce")
    if sn.notna().mean() >= 0.7:
        snc = sn.dropna()
        if snc.nunique() <= 1:
            return None
        x = (sn - sn.mean()) / sn.std(ddof=0)
        X = x.to_numpy().reshape(-1, 1)
        return None if not np.isfinite(X).all() else X
    sc = s.astype(str).fillna("unknown")
    if pd.Series(sc).nunique() <= 1:
        return None
    d = pd.get_dummies(sc, drop_first=True)
    return None if d.shape[1] == 0 else d.to_numpy(float)


def matrix_rank(X: np.ndarray, tol: float = 1e-10) -> int:
    return int(np.linalg.matrix_rank(X, tol=tol))


def ss_effect(G: np.ndarray, X: np.ndarray) -> float:
    XtX = X.T @ X
    M = X.T @ G @ X
    return float(np.trace(M @ np.linalg.pinv(XtX)))


def partial_stats(G: np.ndarray, Z: Optional[np.ndarray], X: np.ndarray):
    n = G.shape[0]
    SS_total = float(np.trace(G))
    if SS_total <= 0:
        return (float("nan"), float("nan"), SS_total, 0, 0, n)
    df_total = n - 1

    if Z is None or Z.size == 0:
        df_X  = matrix_rank(X)
        df_res= df_total - df_X
        SS_X  = ss_effect(G, X)
        SS_res= SS_total - SS_X
        return SS_X, SS_res, SS_total, df_X, df_res, n

    df_Z  = matrix_rank(Z)
    XZ    = np.concatenate([Z, X], axis=1)
    df_XZ = matrix_rank(XZ)
    df_X  = max(0, df_XZ - df_Z)
    df_res= max(0, df_total - df_XZ)

    SS_Z  = ss_effect(G, Z)
    SS_XZ = ss_effect(G, XZ)
    SS_X_given_Z = SS_XZ - SS_Z
    if SS_X_given_Z < 0 and SS_X_given_Z > -1e-10:
        SS_X_given_Z = 0.0
    SS_res = SS_total - SS_XZ
    return SS_X_given_Z, SS_res, SS_total, df_X, df_res, n


def pseudo_F(SS_X: float, df_X: int, SS_res: float, df_res: int) -> float:
    if df_X <= 0 or df_res <= 0:
        return float("nan")
    MS_X = SS_X / df_X
    MS_res = SS_res / df_res
    return float("nan") if MS_res <= 0 else MS_X / MS_res


def permute_within_groups(groups: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    n = len(groups)
    idx = np.arange(n)
    _, inv = np.unique(groups, return_inverse=True)
    for g in np.unique(inv):
        gidx = np.where(inv == g)[0]
        if len(gidx) > 1:
            idx[gidx] = rng.permutation(gidx)
    return idx


def permanova_pvalue(G: np.ndarray, s_var: pd.Series, s_block: pd.Series,
                     n_perm: int = 999, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    X = design_matrix_for_variable(s_var)
    if X is None:
        return float("nan"), float("nan")
    Z = design_matrix_for_variable(s_block)
    SS_X, SS_res, SS_total, df_X, df_res, _ = partial_stats(G, Z, X)
    R2 = SS_X / SS_total if SS_total > 0 else float("nan")
    F_obs = pseudo_F(SS_X, df_X, SS_res, df_res)
    if not np.isfinite(F_obs):
        return R2, float("nan")

    groups = s_block.to_numpy().astype(str)
    ge = 0
    for _ in range(n_perm):
        pidx = permute_within_groups(groups, rng)
        Xp = X[pidx, :]
        SSXp, SSresp, SS_totp, dfXp, dfresp, _ = partial_stats(G, Z, Xp)
        Fp = pseudo_F(SSXp, dfXp, SSresp, dfresp)
        if np.isfinite(Fp) and Fp >= F_obs - 1e-12:
            ge += 1
    p = (ge + 1.0) / (n_perm + 1.0)
    return R2, p


def compute_distance_mats(omics: Dict[str, pd.DataFrame]) -> Dict[str, np.ndarray]:
    Dmats = {}
    for name, tab in omics.items():
        arr = tab.to_numpy(float)
        arr[arr < 0] = 0.0
        rs = arr.sum(axis=1, keepdims=True)
        rs[rs == 0] = 1.0
        Dmats[name] = compute_bray_curtis_matrix(arr / rs)
    return Dmats


def compute_partial_r2_and_pvals(Dmats, md, vars_to_test, block_var,
                                 n_perm=999, seed=None):
    rng = np.random.default_rng(seed)
    r2_rows = []
    p_rows  = []
    s_block_full = md[block_var].astype(str)

    for var in vars_to_test:
        s_var_full = md[var]
        sn = pd.to_numeric(s_var_full, errors="coerce")
        mask = sn.notna().values if sn.notna().mean() >= 0.7 else np.ones(len(s_var_full), bool)
        if mask.sum() < 3:
            continue
        s_var   = s_var_full[mask]
        s_block = s_block_full[mask]

        r2_row = {"metadata": var}
        p_row  = {"metadata": var}
        for name, D in Dmats.items():
            G = gower_center(D[np.ix_(mask, mask)])
            r2p, pval = permanova_pvalue(G, s_var, s_block, n_perm=n_perm, rng=rng)
            r2_row[name] = (r2p * 100.0) if (r2p == r2p) else float("nan")
            p_row[name]  = pval

        r2_rows.append(r2_row)
        p_rows.append(p_row)

    if not r2_rows:
        empty = pd.DataFrame(columns=["metadata"] + list(Dmats.keys())).set_index("metadata")
        return empty, empty

    r2_df = pd.DataFrame(r2_rows).set_index("metadata")
    p_df  = pd.DataFrame(p_rows).set_index("metadata")
    r2_df = r2_df.loc[r2_df.max(axis=1).sort_values(ascending=False).index]
    p_df  = p_df.loc[r2_df.index]
    return r2_df, p_df


def plot_grouped_bars_with_stars(df_r2, df_p, out_png,
                                 label_fs=13, tick_fs=11, legend_fs=11,
                                 title_fs=None, star_fs=16, legend_ncol=1,
                                 font_family=PLASMID_FONT_FAMILY,
                                 font_weight=PLASMID_FONT_WEIGHT,
                                 # NEW: negative values pull the star closer to the bar
                                 star_offset_pts=-2):
    """
    Grouped bars for partial R² with "*" for PERMANOVA p<0.05.

    Customizations:
      - Legend lives INSIDE the axes at top-right (loc='upper right').
      - Stars are placed closer to each bar via a small point offset.
      - X-axis title ("Metadata") removed.
      - Axis tick labels prettified:
          * 'gender' shown as 'Sex'
          * 'age_years' shown as 'Age'
          * underscores replaced with spaces
          * 'anti_pd_1_anti_ctla_4' shown as two rows:
              "Anti-PD-1" on top, "and Anti-CTLA-4" below
    """
    if df_r2.empty:
        print("No results to plot.")
        return

    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import numpy as np

    # --- Helper to pretty-print metadata names on the x-axis ---
    def _format_metadata_label(name: str) -> str:
        # 1 & 3: rename gender -> sex; age_years -> age (for display)
        if name == "gender":
            name = "sex"
        if name == "age_years":
            name = "age"

        # 5: special-case long combo therapy label, make it two rows
        if name == "anti_pd_1_anti_ctla_4":
            return "Anti-PD-1\nand Anti-CTLA-4"

        # 4: underscores -> spaces (for all other labels)
        label = name.replace("_", " ")

        lower = label.lower()
        # Nice capitalization for common fields and treatments
        special = {
            "sex": "Sex",
            "age": "Age",
            "bmi": "BMI",
            "area": "Area",
            "instrument": "Instrument",
            "year": "Year",
            "anti pd 1": "Anti-PD-1",
            "anti ctla 4": "Anti-CTLA-4",
        }
        if lower in special:
            return special[lower]

        # Fallback: title case (e.g. "some var" -> "Some Var")
        return label.title()

    fig, ax = plt.subplots(figsize=(8, 4), dpi=300)

    x = np.arange(len(df_r2.index))
    width = 0.25

    # Bars with plasmid colors
    bars1 = ax.bar(x - width, df_r2["taxa"].values,    width, label="Taxa",
                   color=PLASMID_COLORS.get("taxa", None))
    bars2 = ax.bar(x,          df_r2["bgc"].values,     width, label="BGC",
                   color=PLASMID_COLORS.get("bgc", None))
    bars3 = ax.bar(x + width,  df_r2["pathway"].values, width, label="Pathway",
                   color=PLASMID_COLORS.get("pathway", None))

    # Keep the light y-grid behind bars
    ax.set_axisbelow(True)
    ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.45)

    # Y-axis label only (2: remove x-axis title)
    ax.set_ylabel("% explained variance\n(block = study)",
                  fontsize=label_fs, fontfamily=font_family, fontweight=font_weight)
    # No x-axis label:
    # ax.set_xlabel("Metadata", ...)  # <- removed on purpose

    # X ticks / labels with pretty formatting
    ax.set_xticks(x)
    pretty_labels = [_format_metadata_label(m) for m in df_r2.index]
    ax.set_xticklabels(pretty_labels, rotation=30, ha="right",
                       fontsize=tick_fs, fontfamily=font_family, fontweight=font_weight)

    # Y tick styling
    ax.tick_params(axis="y", labelsize=tick_fs)
    for t in ax.get_yticklabels():
        t.set_fontfamily(font_family)
        t.set_fontweight(font_weight)

    # (Title is still optional; currently disabled as in your version)
    # if title_fs is not None:
    #     fig.suptitle("Partial explained variance by metadata",
    #                  fontsize=title_fs, fontfamily=font_family, fontweight=font_weight)

    # --- Stars: closer to bars via small offset in *points* ---
    # Negative offset pulls the star down (closer); positive pushes up (farther).
    def annotate_stars(bars, col):
        for i, rect in enumerate(bars):
            meta = df_r2.index[i]
            y    = rect.get_height()
            p    = df_p.loc[meta, col]
            if np.isfinite(y) and p == p and p < 0.05:
                x_c = rect.get_x() + rect.get_width() / 2.0
                ax.annotate(
                    "*", xy=(x_c, y), xytext=(0, star_offset_pts),
                    textcoords="offset points",
                    ha="center", va="bottom",
                    fontsize=star_fs, fontfamily=font_family, fontweight=font_weight
                )

    annotate_stars(bars1, "taxa")
    annotate_stars(bars2, "bgc")
    annotate_stars(bars3, "pathway")

    # --- Legend INSIDE the axes at the top-right ---
    star = Line2D([0], [0], marker="*", linestyle="None", label="PERMANOVA p < 0.05")
    h, l = ax.get_legend_handles_labels()
    h.append(star); l.append("PERMANOVA p < 0.05")
    ax.legend(
        h, l,
        fontsize=legend_fs,
        ncol=legend_ncol,
        frameon=False,
        columnspacing=1.2,
        loc="upper right"      # inside, top-right
    )

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)



def run(metadata_path, taxonomy_path, bgc_path, pathway_path,
        block_var="study", vars_to_try=None,
        out_prefix="./permanova_partial_explained_variance_braycurtis_BLOCK_study.no.response",
        n_perm=999, seed=0):
    if vars_to_try is None:
        vars_to_try = [
            "area","age_years","bmi","gender","anti_pd_1","anti_ctla_4",
            "anti_pd_1_anti_ctla_4","instrument","year"
        ]

    md = read_metadata(metadata_path)
    omics_raw = {
        "taxa":    read_omics_table(taxonomy_path, md.index),
        "bgc":     read_omics_table(bgc_path,      md.index),
        "pathway": read_omics_table(pathway_path,  md.index)
    }
    md, omics = align_samples(md, omics_raw)
    if len(md) < 3:
        raise RuntimeError("Not enough overlapping samples.")
    if block_var not in md.columns:
        raise RuntimeError(f"Block variable '{block_var}' not found.")

    Dmats = compute_distance_mats(omics)
    vars_present = [v for v in vars_to_try if v in md.columns and v != block_var]

    r2_df, p_df = compute_partial_r2_and_pvals(
        Dmats, md, vars_present, block_var, n_perm=n_perm, seed=seed
    )

    out_csv   = f"{out_prefix}.csv"
    out_png   = f"{out_prefix}.png"
    out_pcsv  = f"{out_prefix}.pvals.csv"
    out_block = f"{out_prefix}.block.csv"

    r2_df.round(2).to_csv(out_csv)
    p_df.to_csv(out_pcsv)

    # Block-alone R² for 'study'
    block_rows = []
    s_block = md[block_var]
    for name, D in Dmats.items():
        G = gower_center(D)
        Z = design_matrix_for_variable(s_block)
        r2b = float("nan") if Z is None else (ss_effect(G, Z) / float(np.trace(G)))
        block_rows.append({
            "omics": name,
            "R2_percent": round((r2b * 100.0) if r2b == r2b else float("nan"), 2)
        })
    pd.DataFrame(block_rows).set_index("omics").to_csv(out_block)

    # Plot (sizes match your previous call; change here if you want different sizes)
    plot_grouped_bars_with_stars(
        r2_df, p_df, out_png,
        label_fs=14,    # axis labels
        tick_fs=12,     # tick labels
        legend_fs=12,   # legend text
        title_fs=15,    # figure title (None to omit)
        star_fs=18,     # asterisk size
        legend_ncol=1   # stacked legend under plot
    )

    print("Top rows:\n", r2_df.head(20))
    print("\nP-values:\n", p_df.head(20))
    return out_csv, out_png, out_pcsv, out_block


def parse_args():
    p = argparse.ArgumentParser(
        description="PERMANOVA-style partial R² with BLOCK=study; excludes 'response'; with p-values and stars"
    )
    p.add_argument("--metadata", default="/mnt/data/combined_metadata_unknown_everywhere.csv")
    p.add_argument("--taxonomy", default="/mnt/data/taxonomy_fit_adjust_batch2_adj.tsv")
    p.add_argument("--bgc", default="/mnt/data/bgc_fit_adjust_batch2_adj.tsv")
    p.add_argument("--pathway", default="/mnt/data/pathway_fit_adjust_batch2_adj.tsv")
    p.add_argument("--block", default="study")
    p.add_argument("--vars", default="area,age_years,bmi,gender,anti_pd_1,anti_ctla_4,anti_pd_1_anti_ctla_4,instrument,year")
    p.add_argument("--out-prefix", default="/mnt/data/permanova_partial_explained_variance_braycurtis_BLOCK_study.no.response")
    p.add_argument("--perms", type=int, default=999)
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def main():
    a = parse_args()
    vars_to_try = [v.strip() for v in a.vars.split(",") if v.strip()]
    run(a.metadata, a.taxonomy, a.bgc, a.pathway, a.block,
        vars_to_try, a.out_prefix, a.perms, a.seed)


if __name__ == "__main__":
    main()
