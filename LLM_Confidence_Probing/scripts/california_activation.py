"""Experiment B — California's contextual embedding across the 4 prompt-framing runs.

Builds on scripts/prompt_compare.py. In `p3_completion_chat` the answer-slot probability
for "California" collapsed to 0.55 (vs 0.98 for the direct question), with most of the
missing mass going to hedge words. Question: when the model hedges, does it still
*represent* California the same way internally, or has its internal "idea" of California
changed? I.e. is the low P(California) a surface/output effect (hedging) or a deeper
representational difference?

For each run we rebuild the EXACT probe context from prompt_compare.run_probe, append the
single " California" token, run one forward with output_hidden_states=True, and read:
  - cal_emb  = hidden_states[-1][0, -1]  -> last-layer contextual embedding of "California"
  - decision = hidden_states[-1][0, -2]  -> state at the token before it ("is"), whose
               unembedding produces the California logit (where hedging shows up)
  - all-layer California vectors          -> to see at which depth runs diverge
We compare cal_emb / decision across runs (cosine, L2), calibrated against three baselines
(same-context California-vs-Texas, mean-centering, static-embedding displacement), and run a
2x2 ablation testing whether the p3-chat anomaly is the "double-tap" (echo of a dangling
user sentence) or the dangling fragment itself.

Run:  python scripts/california_activation.py
"""
from __future__ import annotations

import os

os.environ.setdefault("HF_HOME", "/workspace/.hf_home")  # cache on workspace volume

import json
from pathlib import Path

import torch
from transformers import AutoModelForCausalLM, AutoTokenizer

MODEL = "Qwen/Qwen2.5-1.5B-Instruct"
READOUT = "The most populous state is"          # experimenter-authored assistant lead-in
NEUTRAL_LEADIN = "The answer is"                # non-echoing lead-in for the ablation
OUT_DIR = Path("/workspace/outputs/california_activation")
PROMPT_COMPARE_RESULTS = Path("/workspace/outputs/prompt_compare/results.json")

# Same 4 runs as prompt_compare.py (id, prompt, mode).
RUNS = [
    ("p1_question_chat",   "What is the most populous state in the United States?",      "chat"),
    ("p2_rephrase_chat",   "Among the 50 states, which one has the largest population?", "chat"),
    ("p3_completion_chat", "Everyone knows the most populous state is",                  "chat"),
    ("p3_completion_raw",  "Everyone knows the most populous state is",                  "raw"),
]


# ----------------------------------------------------------------- context building
def probe_context(tok, prompt: str, mode: str, readout: str = READOUT) -> str:
    """Reproduce prompt_compare.run_probe's `rendered` string exactly."""
    if mode == "chat":
        return tok.apply_chat_template(
            [{"role": "user", "content": prompt}],
            add_generation_prompt=True, tokenize=False) + readout
    return prompt


def ctx_ids_for(tok, rendered: str, device):
    return tok(rendered, add_special_tokens=False, return_tensors="pt").input_ids[0].to(device)


# ------------------------------------------------------------- activation extraction
@torch.no_grad()
def state_vectors(tok, model, device, ctx_ids, state_id: int) -> dict:
    """Append `state_id` to ctx, forward once, and pull the activations of interest.

    Returns float32 CPU tensors.
      emb        : last-layer activation AT the state token (contextual embedding)
      decision   : last-layer activation at the token BEFORE it (predicts the state)
      all_layers : (num_layers+1, hidden) activation at the state token, every layer
      static     : layer-0 (embedding) activation at the state token (context-independent)
      pred_prob  : P(state_id | context) recovered from this same forward's logits
    """
    seq = torch.cat([ctx_ids, torch.tensor([state_id], device=device)])
    out = model(seq.unsqueeze(0), output_hidden_states=True)
    hs = out.hidden_states                                  # tuple len num_layers+1
    pred_logits = out.logits[0, -2].float()                 # logits predicting the state token
    pred_prob = torch.log_softmax(pred_logits, dim=-1)[state_id].exp().item()
    return {
        "emb": hs[-1][0, -1].float().cpu(),
        "decision": hs[-1][0, -2].float().cpu(),
        "all_layers": torch.stack([h[0, -1].float().cpu() for h in hs]),
        "static": hs[0][0, -1].float().cpu(),
        "pred_prob": pred_prob,
    }


# --------------------------------------------------------------------- comparisons
def cos(a: torch.Tensor, b: torch.Tensor) -> float:
    return torch.nn.functional.cosine_similarity(a, b, dim=0).item()


def l2(a: torch.Tensor, b: torch.Tensor) -> float:
    return (a - b).norm().item()


def pairwise(vectors: dict, fn) -> dict:
    ids = list(vectors)
    return {i: {j: fn(vectors[i], vectors[j]) for j in ids} for i in ids}


def main() -> None:
    device = "cuda" if torch.cuda.is_available() else "cpu"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Loading {MODEL} on {device} ...")
    tok = AutoTokenizer.from_pretrained(MODEL)
    model = AutoModelForCausalLM.from_pretrained(MODEL, dtype=torch.bfloat16).to(device).eval()

    cal_id = tok.encode(" California", add_special_tokens=False)
    tex_id = tok.encode(" Texas", add_special_tokens=False)
    assert len(cal_id) == 1 and len(tex_id) == 1, "expected single-token California/Texas"
    cal_id, tex_id = cal_id[0], tex_id[0]

    # Reference probabilities from the original experiment, for the sanity check.
    ref_prob = {}
    if PROMPT_COMPARE_RESULTS.exists():
        ref = json.loads(PROMPT_COMPARE_RESULTS.read_text())
        for r in ref["runs"]:
            ca = next(c for c in r["probe"]["candidates"] if c["state"] == "California")
            ref_prob[r["id"]] = ca["first_token_prob"]

    # ---- main 4 runs: California (+ Texas baseline) at the measured slot --------
    cal, tex, dec, layers, static = {}, {}, {}, {}, {}
    sanity = []
    print("\n=== Extracting California activations per run ===")
    for run_id, prompt, mode in RUNS:
        ctx = ctx_ids_for(tok, probe_context(tok, prompt, mode), device)
        c = state_vectors(tok, model, device, ctx, cal_id)
        t = state_vectors(tok, model, device, ctx, tex_id)
        cal[run_id], tex[run_id] = c["emb"], t["emb"]
        dec[run_id], layers[run_id], static[run_id] = c["decision"], c["all_layers"], c["static"]
        ref_p = ref_prob.get(run_id)
        sanity.append({"run": run_id, "recovered_P(Cal)": c["pred_prob"],
                       "results.json_P(Cal)": ref_p})
        tag = f"(ref {ref_p:.4f})" if ref_p is not None else ""
        print(f"  {run_id:<20} recovered P(California)={c['pred_prob']:.4f} {tag}")

    run_ids = [r[0] for r in RUNS]

    # ---- comparisons -----------------------------------------------------------
    cal_cos, cal_l2 = pairwise(cal, cos), pairwise(cal, l2)
    dec_cos, dec_l2 = pairwise(dec, cos), pairwise(dec, l2)

    # Baseline (a): same-context California-vs-Texas = "one state's worth of difference".
    ca_tx_same_ctx = {r: cos(cal[r], tex[r]) for r in run_ids}

    # Baseline (b): mean-center the 4 California embeddings, inspect residuals.
    cal_stack = torch.stack([cal[r] for r in run_ids])
    cal_mean = cal_stack.mean(0)
    resid = {r: cal[r] - cal_mean for r in run_ids}
    resid_norm = {r: resid[r].norm().item() for r in run_ids}
    resid_cos = pairwise(resid, cos)

    # Baseline (c): how far context moved California from its static embedding,
    # vs the static California<->Texas gap (same number for all runs; static is
    # context-independent so we take it from the first run).
    static_ca, static_tx = static[run_ids[0]], state_vectors(
        tok, model, device, ctx_ids_for(tok, probe_context(tok, *RUNS[0][1:]), device), tex_id)["static"]
    static_ca_tx_gap = l2(static_ca, static_tx)
    ctx_displacement = {r: l2(cal[r], static_ca) for r in run_ids}

    # Layer-wise cosine of California's vector, each run vs p1 (the confident baseline).
    n_layers = layers[run_ids[0]].shape[0]
    layerwise = {r: [cos(layers[r][L], layers["p1_question_chat"][L]) for L in range(n_layers)]
                 for r in run_ids}

    # ---- 2x2 ablation: dangling fragment vs echo/double-tap --------------------
    QUESTION_USER = "What is the most populous state in the United States?"
    DANGLING_USER = "Everyone knows the most populous state is"
    ABLATION = [
        ("question_echo",   QUESTION_USER, READOUT),         # == p1 family
        ("question_neutral", QUESTION_USER, NEUTRAL_LEADIN),
        ("dangling_echo",   DANGLING_USER, READOUT),         # == p3-chat
        ("dangling_neutral", DANGLING_USER, NEUTRAL_LEADIN),
    ]
    print("\n=== 2x2 ablation (P(California) at the slot) ===")
    abl = {}
    abl_emb = {}
    for cell_id, user, leadin in ABLATION:
        ctx = ctx_ids_for(tok, probe_context(tok, user, "chat", leadin), device)
        v = state_vectors(tok, model, device, ctx, cal_id)
        abl[cell_id] = v["pred_prob"]
        abl_emb[cell_id] = v["emb"]
        print(f"  {cell_id:<18} P(California)={v['pred_prob']:.4f}")
    abl_emb_cos = pairwise(abl_emb, cos)

    # ---- save ------------------------------------------------------------------
    analysis = {
        "model": MODEL,
        "sanity_check": sanity,
        "cal_emb": {"cosine": cal_cos, "l2": cal_l2},
        "decision": {"cosine": dec_cos, "l2": dec_l2},
        "baseline_a_california_vs_texas_same_context": ca_tx_same_ctx,
        "baseline_b_mean_centered": {"residual_norm": resid_norm, "residual_cosine": resid_cos},
        "baseline_c_static_displacement": {
            "static_california_texas_gap_L2": static_ca_tx_gap,
            "context_displacement_from_static_L2": ctx_displacement,
        },
        "layerwise_cosine_vs_p1": layerwise,
        "ablation_2x2": {"P(California)": abl, "cal_emb_cosine": abl_emb_cos},
    }
    (OUT_DIR / "analysis.json").write_text(json.dumps(analysis, indent=2))

    torch.save({"cal_emb": cal, "decision": dec, "texas_emb": tex,
                "cal_all_layers": layers, "cal_static": static,
                "ablation_cal_emb": abl_emb}, OUT_DIR / "activations.pt")

    _plot(run_ids, cal_cos, dec_cos, ca_tx_same_ctx, layerwise, abl)
    _print_summary(run_ids, cal_cos, ca_tx_same_ctx, resid_norm, ctx_displacement,
                   static_ca_tx_gap, abl)
    print(f"\nSaved {OUT_DIR/'analysis.json'}, {OUT_DIR/'activations.pt'}, "
          f"{OUT_DIR/'similarity.png'}")


def _print_summary(run_ids, cal_cos, ca_tx, resid_norm, displacement, static_gap, abl) -> None:
    p3 = "p3_completion_chat"
    print("\n" + "=" * 78 + "\nSUMMARY\n" + "=" * 78)
    print("California-embedding cosine, p3_completion_chat vs the others "
          "(yardstick = CA-vs-TX in same context):")
    for r in run_ids:
        if r == p3:
            continue
        print(f"  cos(CA_p3, CA_{r:<18}) = {cal_cos[p3][r]:.4f}   "
              f"| CA-vs-TX in {r} = {ca_tx[r]:.4f}")
    print("\nMean-centered residual norms (larger = more run-specific deviation):")
    for r in run_ids:
        print(f"  {r:<20} {resid_norm[r]:.3f}")
    print(f"\nContext displacement from static California embedding "
          f"(static CA-TX gap = {static_gap:.2f}):")
    for r in run_ids:
        print(f"  {r:<20} {displacement[r]:.2f}")
    print("\n2x2 ablation P(California):")
    print(f"  {'':<12}{'echo':>10}{'neutral':>10}")
    print(f"  {'question':<12}{abl['question_echo']:>10.4f}{abl['question_neutral']:>10.4f}")
    print(f"  {'dangling':<12}{abl['dangling_echo']:>10.4f}{abl['dangling_neutral']:>10.4f}")


def _plot(run_ids, cal_cos, dec_cos, ca_tx, layerwise, abl) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    short = [r.replace("_completion", "").replace("_question", "").replace("_chat", "-chat")
             .replace("_raw", "-raw") for r in run_ids]
    fig, axes = plt.subplots(2, 2, figsize=(13, 11))

    def heat(ax, mat, title):
        m = np.array([[mat[i][j] for j in run_ids] for i in run_ids])
        im = ax.imshow(m, vmin=m[~np.eye(len(m), dtype=bool)].min(), vmax=1.0, cmap="viridis")
        ax.set_xticks(range(len(short))); ax.set_xticklabels(short, rotation=30, ha="right", fontsize=8)
        ax.set_yticks(range(len(short))); ax.set_yticklabels(short, fontsize=8)
        for i in range(len(m)):
            for j in range(len(m)):
                ax.text(j, i, f"{m[i, j]:.3f}", ha="center", va="center",
                        color="w" if m[i, j] < 0.5 * (m.max() + m.min()) else "k", fontsize=8)
        ax.set_title(title, fontsize=10); fig.colorbar(im, ax=ax, fraction=0.046)

    heat(axes[0, 0], cal_cos, "cos(cal_emb) across runs")
    heat(axes[0, 1], dec_cos, "cos(decision state) across runs")

    ax = axes[1, 0]
    for r in run_ids:
        ax.plot(range(len(layerwise[r])), layerwise[r], marker=".", label=r, linewidth=1)
    ax.set_xlabel("layer (0 = static embedding)"); ax.set_ylabel("cosine vs p1_question_chat")
    ax.set_title("California's representation by depth (vs p1)", fontsize=10)
    ax.legend(fontsize=7); ax.grid(alpha=0.3)

    ax = axes[1, 1]
    grid = np.array([[abl["question_echo"], abl["question_neutral"]],
                     [abl["dangling_echo"], abl["dangling_neutral"]]])
    im = ax.imshow(grid, vmin=0, vmax=1, cmap="RdYlGn")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["echo lead-in", "neutral lead-in"])
    ax.set_yticks([0, 1]); ax.set_yticklabels(["question user", "dangling user"])
    for i in range(2):
        for j in range(2):
            ax.text(j, i, f"{grid[i, j]:.3f}", ha="center", va="center", fontsize=12)
    ax.set_title("2x2 ablation: P(California)", fontsize=10); fig.colorbar(im, ax=ax, fraction=0.046)

    fig.suptitle("Experiment B — California's contextual embedding across framings", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(OUT_DIR / "similarity.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
