"""Compare Qwen's response & confidence across three same-meaning prompts.

Three prompts mean the same thing but are framed differently:
  1. direct question         2. rephrased question      3. leading completion stub
Prompt #3 is run BOTH as a raw completion and through the chat template, so we
also see the chat-vs-raw format effect. Four runs total.

Each run gets two analyses:
  A. Free generation  — greedy decode + per-token logprobs (what it says, how sure).
  B. Answer-slot probe — at a fixed slot whose next token is the state name, compute
     P("<state>") for all 50 states and compare (esp. California vs Texas).

Probe efficiency: a forward pass gives next-token conditionals only along the tokens
actually in its input, so we use exactly TWO forwards per prompt:
  - one over the context  -> P(first token) for all 50 states at once
    (the 40 single-token names, incl. California/Texas/Florida, are fully done here);
  - one BATCHED over the 5 distinct multi-token prefixes (New/North/South/Rhode/West)
    -> P(second token | prefix) for all 10 two-token names at once.

Outputs are printed and saved under outputs/prompt_compare/.
Run:  python scripts/prompt_compare.py
"""
from __future__ import annotations

import os

os.environ.setdefault("HF_HOME", "/workspace/.hf_home")  # cache on workspace volume

import json
from pathlib import Path

import torch
from transformers import AutoModelForCausalLM, AutoTokenizer

MODEL = "Qwen/Qwen2.5-1.5B-Instruct"
READOUT = "The most populous state is"      # assistant lead-in for chat runs
MAX_NEW_TOKENS = 100
TOP_N = 10
OUT_DIR = Path("/workspace/outputs/prompt_compare")

STATES = [
    "Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado",
    "Connecticut", "Delaware", "Florida", "Georgia", "Hawaii", "Idaho", "Illinois",
    "Indiana", "Iowa", "Kansas", "Kentucky", "Louisiana", "Maine", "Maryland",
    "Massachusetts", "Michigan", "Minnesota", "Mississippi", "Missouri", "Montana",
    "Nebraska", "Nevada", "New Hampshire", "New Jersey", "New Mexico", "New York",
    "North Carolina", "North Dakota", "Ohio", "Oklahoma", "Oregon", "Pennsylvania",
    "Rhode Island", "South Carolina", "South Dakota", "Tennessee", "Texas", "Utah",
    "Vermont", "Virginia", "Washington", "West Virginia", "Wisconsin", "Wyoming",
]

# (id, prompt text, mode) — mode is "chat" (chat template) or "raw" (plain completion).
RUNS = [
    ("p1_question_chat",   "What is the most populous state in the United States?",      "chat"),
    ("p2_rephrase_chat",   "Among the 50 states, which one has the largest population?", "chat"),
    ("p3_completion_chat", "Everyone knows the most populous state is",                  "chat"),
    ("p3_completion_raw",  "Everyone knows the most populous state is",                  "raw"),
]


def main() -> None:
    device = "cuda" if torch.cuda.is_available() else "cpu"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Loading {MODEL} on {device} ...")
    tok = AutoTokenizer.from_pretrained(MODEL)
    model = AutoModelForCausalLM.from_pretrained(MODEL, dtype=torch.bfloat16).to(device).eval()
    state_ids = {s: tok.encode(" " + s, add_special_tokens=False) for s in STATES}
    assert all(len(v) <= 2 for v in state_ids.values()), "probe assumes <=2-token state names"

    results = []
    for run_id, prompt, mode in RUNS:
        print("\n" + "=" * 78)
        print(f"RUN {run_id}  (mode={mode})\n  prompt: {prompt!r}")
        gen_block = run_generation(tok, model, device, prompt, mode)
        probe_block = run_probe(tok, model, device, prompt, mode, state_ids)
        results.append({"id": run_id, "mode": mode, "prompt": prompt,
                        "generation": gen_block, "probe": probe_block})

    print_summary(results)
    save_artifacts(results)


# ---------------------------------------------------------------- generation (A)
def build_inputs(tok, prompt: str, mode: str, device):
    if mode == "chat":
        enc = tok.apply_chat_template(
            [{"role": "user", "content": prompt}],
            add_generation_prompt=True, return_tensors="pt", return_dict=True,
        )
    else:
        enc = tok(prompt, return_tensors="pt")
    return {k: v.to(device) for k, v in enc.items()}


def run_generation(tok, model, device, prompt: str, mode: str) -> dict:
    inputs = build_inputs(tok, prompt, mode, device)
    prompt_len = inputs["input_ids"].shape[1]
    with torch.no_grad():
        gen = model.generate(
            **inputs, max_new_tokens=MAX_NEW_TOKENS, do_sample=False,
            return_dict_in_generate=True, output_logits=True,
            pad_token_id=tok.eos_token_id,
        )
    new_ids = gen.sequences[0, prompt_len:]
    text = tok.decode(new_ids, skip_special_tokens=True).strip()
    per_token = []
    for step, logits in enumerate(gen.logits):
        tid = new_ids[step].item()
        lp = torch.log_softmax(logits[0].float(), dim=-1)[tid].item()
        per_token.append({"token": tok.decode([tid]), "logprob": lp})
    print(f"  generated ({len(per_token)} tok): {text!r}")
    return {"text": text, "tokens": per_token}


# --------------------------------------------------------------------- probe (B)
def run_probe(tok, model, device, prompt: str, mode: str, state_ids: dict) -> dict:
    if mode == "chat":
        rendered = tok.apply_chat_template(
            [{"role": "user", "content": prompt}],
            add_generation_prompt=True, tokenize=False) + READOUT
    else:
        rendered = prompt
    ctx_ids = tok(rendered, add_special_tokens=False, return_tensors="pt").input_ids[0].to(device)
    ctx_len = ctx_ids.shape[0]

    # Forward #1: distribution over the first answer token (shared by all 50 states).
    with torch.no_grad():
        first_logdist = torch.log_softmax(
            model(ctx_ids.unsqueeze(0)).logits[0, -1].float(), dim=-1)

    # Forward #2 (batched): one row per distinct first token of a two-token state,
    # giving P(second token | context, prefix) for every state in that group.
    prefixes = sorted({ids[0] for ids in state_ids.values() if len(ids) == 2})
    prefix_logdist = {}
    if prefixes:
        batch = torch.stack([torch.cat([ctx_ids, torch.tensor([p], device=device)])
                             for p in prefixes])           # (n_prefix, ctx_len+1)
        with torch.no_grad():
            logits = model(batch).logits[:, -1].float()    # (n_prefix, V)
        for row, p in enumerate(prefixes):
            prefix_logdist[p] = torch.log_softmax(logits[row], dim=-1)

    cands = []
    for state, ids in state_ids.items():
        lp = first_logdist[ids[0]].item()
        if len(ids) == 2:
            lp += prefix_logdist[ids[0]][ids[1]].item()
        cands.append({"state": state, "n_tokens": len(ids), "logprob": lp,
                      "prob": float(torch.tensor(lp).exp()),
                      "first_token_prob": first_logdist[ids[0]].exp().item()})

    tot = sum(c["prob"] for c in cands)
    for c in cands:
        c["relative_share"] = c["prob"] / tot if tot > 0 else 0.0
    cands.sort(key=lambda c: -c["prob"])
    ranks = {c["state"]: i + 1 for i, c in enumerate(cands)}

    top = torch.topk(first_logdist, 5)
    top5 = [{"token": tok.decode([i]), "prob": lp.exp().item()}
            for lp, i in zip(top.values, top.indices)]

    pc = next(c["prob"] for c in cands if c["state"] == "California")
    pt = next(c["prob"] for c in cands if c["state"] == "Texas")
    print("  top states: " + ", ".join(f"{c['state']}={c['prob']:.3f}" for c in cands[:5]))
    print(f"  California rank {ranks['California']} (p={pc:.3f}) | "
          f"Texas rank {ranks['Texas']} (p={pt:.3f}) | "
          f"argmax {top5[0]['token']!r} ({top5[0]['prob']:.3f})")
    return {"readout_context": rendered, "candidates": cands, "ranks": ranks, "top5": top5}


# ----------------------------------------------------------------- summary + save
def _state(probe, name, key="prob"):
    return next(c[key] for c in probe["candidates"] if c["state"] == name)


def print_summary(results: list[dict]) -> None:
    print("\n" + "=" * 92 + "\nSUMMARY\n" + "=" * 92)
    print(f"{'run':<20} {'mode':<5} {'P(Cal)':>8} {'P(Tex)':>8} {'Cal:Tex':>8} "
          f"{'argmax':<14} {'CalRk':>5} {'TexRk':>5}")
    for r in results:
        p = r["probe"]
        pc, pt = _state(p, "California"), _state(p, "Texas")
        ratio = pc / pt if pt > 0 else float("inf")
        print(f"{r['id']:<20} {r['mode']:<5} {pc:>8.4f} {pt:>8.4f} {ratio:>8.1f} "
              f"{p['top5'][0]['token']!r:<14} {p['ranks']['California']:>5} {p['ranks']['Texas']:>5}")


def save_artifacts(results: list[dict]) -> None:
    with open(OUT_DIR / "results.json", "w") as f:
        json.dump({"model": MODEL, "states": STATES, "runs": results}, f, indent=2)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    for ax, r in zip(axes.flat, results):
        top = r["probe"]["candidates"][:TOP_N][::-1]
        names = [c["state"] for c in top]
        probs = [c["prob"] for c in top]
        colors = ["#d62728" if n == "California" else "#ff7f0e" if n == "Texas"
                  else "#7f7f7f" for n in names]
        ax.barh(names, probs, color=colors)
        ax.set_xlim(0, 1)                       # same scale across all panels
        ax.set_title(f"{r['id']} ({r['mode']})", fontsize=10)
        ax.set_xlabel("P(model answers with this state)")
    fig.suptitle("Top-10 state probabilities by prompt framing "
                 "(California=red, Texas=orange)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(OUT_DIR / "candidate_probs.png", dpi=150)
    plt.close(fig)
    print(f"\nSaved {OUT_DIR/'results.json'} and {OUT_DIR/'candidate_probs.png'}")


if __name__ == "__main__":
    main()
