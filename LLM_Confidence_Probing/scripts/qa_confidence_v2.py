"""Experiment A v2 — de-confounded, de-saturated rerun.

v1 found only weak support, and Exp B showed the *original* motivating observation (p3 "confident"
phrasing -> lowest P(California)) was a double-tap artifact, not a confidence effect. v1 also had
confounds: mixed prompt endings ('?' vs 'is'/'be'), first token = answer-entity vs "The", saturated
easy facts, and a DEAD hedge measure (hedge words were in the prompt, never the response). v2 fixes:

  - All prompts are open-ended STUBS ending in "is" (3 framings, one per tier) -> d is built on the
    same final "is" token, isolating the preamble's confidence effect.
  - 11 facts = 6 easy + 5 hard/uncertain (de-saturate belief).
  - Hedge measure = fraction of *generated* tokens that are hedge/uncertainty words (whole-response
    scan, expanded list) -> actually fires on hard questions.

Reuses the structure of scripts/qa_confidence.py.
Outputs: outputs/qa_confidence_v2/{data.json, d.pt, scatter_qa.png}
Run:  python scripts/qa_confidence_v2.py
"""
from __future__ import annotations

import os

os.environ.setdefault("HF_HOME", "/workspace/.hf_home")

import json
from pathlib import Path

import numpy as np
import torch
from transformers import AutoModelForCausalLM, AutoTokenizer

MODEL = "Qwen/Qwen2.5-1.5B-Instruct"
OUT_DIR = Path("/workspace/outputs/qa_confidence_v2")
MAX_NEW_TOKENS = 50

# (id, subject phrase [lowercase 'the ...'], reference answer for eyeballing — NOT auto-graded)
FACTS = [
    ("populous_state",   "the most populous US state",                       "California"),
    ("tallest_mountain", "the tallest mountain in the world",                "Mount Everest"),
    ("largest_ocean",    "the largest ocean on Earth",                       "the Pacific"),
    ("capital_france",   "the capital of France",                            "Paris"),
    ("largest_planet",   "the largest planet in the solar system",           "Jupiter"),
    ("longest_river",    "the longest river in the world",                   "the Nile (or Amazon)"),
    ("bbn_element",      "the heaviest element created during Big Bang nucleosynthesis",
                         "Lithium (~Li-7)"),
    ("us_tax_2025",      "the total tax revenue collected by the US federal government in 2025, "
                         "rounded to the nearest trillion dollars",          "~$5 trillion"),
    ("piano_tuners",     "the number of piano tuners living in Chicago",     "~60-130 (Fermi)"),
    ("asteroid_belt",    "the largest object in the main asteroid belt",     "Ceres"),
    ("smallest_animal",  "the smallest animal in the world",                 "myxozoan spores (debated)"),
]

# (framing id, tier, builder).  confident/hedged build d; neutral is a held-out midpoint.
FRAMINGS = [
    ("confident", "confident", lambda s: f"Everyone knows that {s} is"),
    ("neutral",   "neutral",   lambda s: f"{s[0].upper()}{s[1:]} is"),
    ("hedged",    "hedged",    lambda s: f"I'm not sure, but I think that {s} is"),
]

HEDGE_WORDS = ["currently", "usually", "typically", "generally", "maybe", "perhaps", "probably",
               "possibly", "likely", "might", "sometimes", "often", "arguably", "supposedly",
               "apparently", "presumably", "approximately", "around", "roughly", "about",
               "estimated", "estimate", "difficult", "hard", "unknown", "uncertain", "unclear",
               "debated", "disputed", "somewhere", "thought", "considered", "believed", "could",
               "potentially", "speculative", "depends", "varies", "may"]
HEDGE_SET = {w.lower() for w in HEDGE_WORDS}

TIER_COLOR = {"confident": "#2ca02c", "neutral": "#7f7f7f", "hedged": "#d62728"}


def render(tok, prompt: str) -> str:
    return tok.apply_chat_template([{"role": "user", "content": prompt}],
                                   add_generation_prompt=True, tokenize=False)


def final_question_idx(tok, ids: list[int]) -> int:
    im_end = tok.convert_tokens_to_ids("<|im_end|>")
    ends = [i for i, t in enumerate(ids) if t == im_end]
    return ends[-1] - 1                                     # the final "is" token


def prev_content_idx(tok, ids: list[int], start: int) -> int:
    """Walk back from `start` over tokens with no alphanumeric char (-> subject's last word)."""
    i = start
    while i > 0 and not any(c.isalnum() for c in tok.decode([ids[i]])):
        i -= 1
    return i


@torch.no_grad()
def process(tok, model, device, prompt: str, hedge_ids: list[int]) -> dict:
    enc = tok(render(tok, prompt), add_special_tokens=False, return_tensors="pt").to(device)
    prompt_ids = enc.input_ids[0]
    plen = prompt_ids.shape[0]
    ids = prompt_ids.tolist()
    is_idx = final_question_idx(tok, ids)
    word_idx = prev_content_idx(tok, ids, is_idx - 1)

    gen = model.generate(**enc, max_new_tokens=MAX_NEW_TOKENS, do_sample=False,
                         return_dict_in_generate=True, output_logits=True,
                         pad_token_id=tok.eos_token_id)
    new_ids = gen.sequences[0, plen:]
    # truncate at eos for the response scan
    eos = tok.eos_token_id
    resp_ids = []
    for t in new_ids.tolist():
        if t == eos:
            break
        resp_ids.append(t)
    text = tok.decode(resp_ids, skip_special_tokens=True).strip()

    first_logprobs = torch.log_softmax(gen.logits[0][0].float(), dim=-1)
    first_id = new_ids[0].item()
    first_token_logprob = first_logprobs[first_id].item()
    probs = first_logprobs.exp()
    hedge_mass_first = float(sum(probs[h].item() for h in hedge_ids))

    # whole-response hedge scan (the fixed measure)
    n_hedge = sum(1 for t in resp_ids if tok.decode([t]).strip().lower() in HEDGE_SET)
    hedge_frac = n_hedge / max(len(resp_ids), 1)

    out = model(torch.cat([prompt_ids, new_ids[:1]]).unsqueeze(0), output_hidden_states=True)
    last = out.hidden_states[-1][0].float().cpu()
    return {
        "text": text, "first_token": tok.decode([first_id]),
        "first_token_logprob": first_token_logprob,
        "hedge_mass_first": hedge_mass_first, "hedge_frac": hedge_frac,
        "q_is": last[is_idx], "q_word": last[word_idx], "resp": last[plen],
        "is_tok": tok.decode([ids[is_idx]]), "word_tok": tok.decode([ids[word_idx]]),
    }


def unit(v):
    return v / v.norm()


def build_d(Q, tiers, facts, exclude=None):
    conf = [Q[i] for i in range(len(Q)) if tiers[i] == "confident" and facts[i] != exclude]
    hedg = [Q[i] for i in range(len(Q)) if tiers[i] == "hedged" and facts[i] != exclude]
    return unit(torch.stack(conf).mean(0) - torch.stack(hedg).mean(0))


def pearson(x, y):
    return float(np.corrcoef(x, y)[0, 1])


def spearman(x, y):
    return float(np.corrcoef(np.argsort(np.argsort(x)), np.argsort(np.argsort(y)))[0, 1])


def main() -> None:
    device = "cuda" if torch.cuda.is_available() else "cpu"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Loading {MODEL} on {device} ...")
    tok = AutoTokenizer.from_pretrained(MODEL)
    model = AutoModelForCausalLM.from_pretrained(MODEL, dtype=torch.bfloat16).to(device).eval()
    hedge_ids = sorted({tok.encode(" " + w, add_special_tokens=False)[0] for w in HEDGE_WORDS})

    rows = []
    print("\n=== Generating (33 prompts) ===")
    for fid, subject, answer in FACTS:
        print(f"\n--- {fid}  (ref: {answer}) ---")
        for frid, tier, build in FRAMINGS:
            prompt = build(subject)
            r = process(tok, model, device, prompt, hedge_ids)
            r.update({"fact": fid, "framing": frid, "tier": tier, "prompt": prompt, "answer": answer})
            rows.append(r)
            print(f"  [{tier:<9}] lp={r['first_token_logprob']:>6.3f} hedge_frac={r['hedge_frac']:.2f}"
                  f"  first={r['first_token']!r:<14} | {r['text'][:64]!r}")

    tiers = [r["tier"] for r in rows]
    facts = [r["fact"] for r in rows]

    proj, d_all = {}, {}
    for key in ("q_is", "q_word"):
        Q = [r[key] for r in rows]
        R = [r["resp"] for r in rows]
        d_all[key] = build_d(Q, tiers, facts)
        qp, rp = [], []
        for i in range(len(rows)):
            d = build_d(Q, tiers, facts, exclude=facts[i])
            qp.append(float(torch.dot(Q[i], d)))
            rp.append(float(torch.dot(R[i], d)))
        proj[key] = {"q_proj": qp, "resp_proj": rp}

    flp = [r["first_token_logprob"] for r in rows]
    hfrac = [r["hedge_frac"] for r in rows]

    corr = {}
    for key in ("q_is", "q_word"):
        qp = proj[key]["q_proj"]
        corr[key] = {name: {"pearson": pearson(qp, y), "spearman": spearman(qp, y)}
                     for name, y in (("first_token_logprob", flp), ("hedge_frac", hfrac),
                                     ("resp_proj", proj[key]["resp_proj"]))}

    _save(rows, proj, corr, d_all)
    _plot(rows, proj, flp, hfrac, corr)
    _print_corr(corr)
    print(f"\nSaved {OUT_DIR/'data.json'}, {OUT_DIR/'d.pt'}, {OUT_DIR/'scatter_qa.png'}")


def _print_corr(corr) -> None:
    print("\n" + "=" * 72 + "\nCORRELATIONS  (question-confidence Q·d  vs  answer measure)\n" + "=" * 72)
    print(f"{'Q definition':<12}{'answer measure':<22}{'Pearson':>9}{'Spearman':>10}")
    for key in corr:
        for name, c in corr[key].items():
            print(f"{key:<12}{name:<22}{c['pearson']:>9.3f}{c['spearman']:>10.3f}")


def _save(rows, proj, corr, d_all) -> None:
    data = {"model": MODEL, "facts": FACTS, "hedge_words": HEDGE_WORDS,
            "correlations": corr, "prompts": []}
    for i, r in enumerate(rows):
        data["prompts"].append({
            "fact": r["fact"], "framing": r["framing"], "tier": r["tier"], "prompt": r["prompt"],
            "reference_answer": r["answer"], "generated": r["text"], "first_token": r["first_token"],
            "first_token_logprob": r["first_token_logprob"],
            "hedge_mass_first": r["hedge_mass_first"], "hedge_frac": r["hedge_frac"],
            "is_tok": r["is_tok"], "word_tok": r["word_tok"],
            "q_proj_is": proj["q_is"]["q_proj"][i], "q_proj_word": proj["q_word"]["q_proj"][i],
            "resp_proj_is": proj["q_is"]["resp_proj"][i],
            "resp_proj_word": proj["q_word"]["resp_proj"][i],
        })
    (OUT_DIR / "data.json").write_text(json.dumps(data, indent=2))
    torch.save({"d_all_is": d_all["q_is"], "d_all_word": d_all["q_word"]}, OUT_DIR / "d.pt")


def _plot(rows, proj, flp, hfrac, corr) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = [TIER_COLOR[r["tier"]] for r in rows]
    measures = [("first_token_logprob", flp), ("hedge_frac", hfrac), ("resp_proj", None)]
    defs = [("q_is", "final token: 'is'"), ("q_word", "subject last word")]
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    for row, (key, dlabel) in enumerate(defs):
        qp = proj[key]["q_proj"]
        for col, (name, y) in enumerate(measures):
            yy = proj[key]["resp_proj"] if name == "resp_proj" else y
            ax = axes[row, col]
            ax.scatter(qp, yy, c=colors, s=45, edgecolor="k", linewidth=0.3)
            c = corr[key][name]
            ax.set_title(f"{name}  (r={c['pearson']:.2f}, ρ={c['spearman']:.2f})", fontsize=10)
            ax.set_xlabel(f"question confidence  Q·d  [{dlabel}]")
            ax.set_ylabel(name)
            ax.grid(alpha=0.3)
    handles = [plt.Line2D([0], [0], marker="o", ls="", mec="k", mfc=c, label=t)
               for t, c in TIER_COLOR.items()]
    fig.legend(handles=handles, loc="upper right", fontsize=9)
    fig.suptitle("Experiment A v2 — question confidence (Q·d) vs answer confidence "
                 "(stubs ending in 'is'; 11 facts incl. hard)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(OUT_DIR / "scatter_qa.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
