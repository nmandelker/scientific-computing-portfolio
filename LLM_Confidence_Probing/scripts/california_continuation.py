"""Experiment B, continuation probe — what does each run predict AFTER "...is California"?

The last-layer activation at the California token (cal_emb) is a next-token *planning*
state, so the four runs that all confidently answer "California" can still differ in what
they intend to say next. This script makes that concrete: for each run it appends the
" California" token to the exact probe context and reads the next-token distribution,
saving the top-5 continuations as a table + bar plot.

Outputs: outputs/california_activation/continuation.json and continuation_top5.png
Run:  python scripts/california_continuation.py
"""
from __future__ import annotations

import json
from pathlib import Path

import torch
from transformers import AutoModelForCausalLM, AutoTokenizer

from scripts.california_activation import (MODEL, RUNS, OUT_DIR, probe_context, ctx_ids_for)

TOP_K = 5


@torch.no_grad()
def next_token_dist(tok, model, device, ctx_ids, cal_id: int):
    seq = torch.cat([ctx_ids, torch.tensor([cal_id], device=device)])
    logits = model(seq.unsqueeze(0)).logits[0, -1].float()      # predicts token AFTER California
    return torch.softmax(logits, dim=-1)


def main() -> None:
    device = "cuda" if torch.cuda.is_available() else "cpu"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Loading {MODEL} on {device} ...")
    tok = AutoTokenizer.from_pretrained(MODEL)
    model = AutoModelForCausalLM.from_pretrained(MODEL, dtype=torch.bfloat16).to(device).eval()
    cal_id = tok.encode(" California", add_special_tokens=False)[0]

    results = {}
    for run_id, prompt, mode in RUNS:
        ctx = ctx_ids_for(tok, probe_context(tok, prompt, mode), device)
        probs = next_token_dist(tok, model, device, ctx, cal_id)
        v, i = probs.topk(TOP_K)
        results[run_id] = [{"token": tok.decode([t]), "repr": repr(tok.decode([t])),
                            "token_id": int(t), "prob": float(p)} for p, t in zip(v, i)]

    (OUT_DIR / "continuation.json").write_text(json.dumps(results, indent=2))

    # ---- console table --------------------------------------------------------
    print('\nTop-5 tokens predicted AFTER "...the most populous state is California":\n')
    header = "rank  " + "".join(f"{rid:<26}" for rid in results)
    print(header)
    for k in range(TOP_K):
        row = f"{k+1:<6}"
        for rid in results:
            c = results[rid][k]
            row += f"{c['repr']+' '+format(c['prob'],'.3f'):<26}"
        print(row)

    _plot(results)
    print(f"\nSaved {OUT_DIR/'continuation.json'} and {OUT_DIR/'continuation_top5.png'}")


def _plot(results: dict) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for ax, (rid, top) in zip(axes.flat, results.items()):
        labels = [c["repr"] for c in top][::-1]
        probs = [c["prob"] for c in top][::-1]
        ax.barh(labels, probs, color="#4c72b0")
        ax.set_xlim(0, 1)
        ax.set_title(rid, fontsize=10)
        ax.set_xlabel("P(next token | ...is California)")
        for y, p in enumerate(probs):
            ax.text(p + 0.01, y, f"{p:.2f}", va="center", fontsize=8)
    fig.suptitle('What each run plans to say AFTER "California" '
                 '(top-5 next tokens)', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(OUT_DIR / "continuation_top5.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
