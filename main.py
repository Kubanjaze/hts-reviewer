import sys
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

import argparse, os, json, re, warnings
warnings.filterwarnings("ignore")
import pandas as pd
from dotenv import load_dotenv
import anthropic
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")
load_dotenv()
os.environ.setdefault("ANTHROPIC_API_KEY", os.getenv("ANTHROPIC_API_KEY", ""))

TRIAGE_SYSTEM = """You are an HTS hit triage expert. Given a compound's activity data and computed properties, assess its development priority.

Respond with ONLY valid JSON:
{
  "priority": "high" | "medium" | "low",
  "liability_flags": ["flag1", "flag2"],
  "rationale": "1-2 sentence explanation",
  "next_step": "recommended next action"
}

Priority criteria:
- HIGH: pIC50 >= 8.0, good drug-likeness (RO5 pass), no major flags
- MEDIUM: pIC50 7.0-8.0, acceptable properties, minor flags
- LOW: pIC50 < 7.0, poor properties, or major liability flags

Be specific about liability flags (e.g., "high_logp", "low_solubility_risk", "reactive_group")."""


def compute_properties(smiles: str) -> dict:
    """Compute RDKit drug-likeness properties."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"error": "Invalid SMILES"}
    return {
        "mw": round(Descriptors.MolWt(mol), 1),
        "logp": round(Descriptors.MolLogP(mol), 2),
        "hba": Descriptors.NumHAcceptors(mol),
        "hbd": Descriptors.NumHDonors(mol),
        "tpsa": round(Descriptors.TPSA(mol), 1),
        "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
        "ro5_violations": sum([
            Descriptors.MolWt(mol) > 500,
            Descriptors.MolLogP(mol) > 5,
            Descriptors.NumHAcceptors(mol) > 10,
            Descriptors.NumHDonors(mol) > 5,
        ]),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Phase 77 — Claude reviewer for HTS top-N hits",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", required=True, help="Compounds CSV with pic50")
    parser.add_argument("--n", type=int, default=3, help="Number of top hits to review")
    parser.add_argument("--model", default="claude-haiku-4-5-20251001", help="Model ID")
    parser.add_argument("--output-dir", default="output", help="Output directory")
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    df = pd.read_csv(args.input).sort_values("pic50", ascending=False).head(args.n)
    client = anthropic.Anthropic()

    print(f"\nPhase 77 — Claude HTS Hit Triage")
    print(f"Model: {args.model} | Top-{args.n} hits by pIC50\n")

    reviews = []
    total_input = 0
    total_output = 0

    for rank, (_, row) in enumerate(df.iterrows(), 1):
        props = compute_properties(row["smiles"])
        props_str = ", ".join(f"{k}={v}" for k, v in props.items())

        user_msg = (
            f"HTS Hit #{rank}\n"
            f"Name: {row['compound_name']}\n"
            f"SMILES: {row['smiles']}\n"
            f"pIC50: {row['pic50']:.2f}\n"
            f"Computed properties: {props_str}\n\n"
            f"Triage this hit. Respond with JSON only."
        )

        response = client.messages.create(
            model=args.model,
            max_tokens=256,
            system=TRIAGE_SYSTEM,
            messages=[{"role": "user", "content": user_msg}],
        )
        text = "".join(b.text for b in response.content if hasattr(b, "text"))
        total_input += response.usage.input_tokens
        total_output += response.usage.output_tokens

        match = re.search(r'\{.*\}', text, re.DOTALL)
        if match:
            review = json.loads(match.group())
        else:
            review = {"priority": "parse_error", "raw": text}

        review["rank"] = rank
        review["compound"] = row["compound_name"]
        review["pic50"] = row["pic50"]
        review["properties"] = props
        reviews.append(review)

        print(f"  #{rank} {row['compound_name']:20s} pIC50={row['pic50']:.2f} -> {review.get('priority', '?').upper()}")
        if review.get("liability_flags"):
            print(f"     Flags: {', '.join(review['liability_flags'])}")
        print(f"     Next: {review.get('next_step', 'N/A')[:70]}")
        print()

    cost = (total_input / 1e6 * 0.80) + (total_output / 1e6 * 4.0)

    # Priority distribution
    priority_counts = {}
    for r in reviews:
        p = r.get("priority", "unknown")
        priority_counts[p] = priority_counts.get(p, 0) + 1

    report = (
        f"Phase 77 — HTS Hit Triage Report\n{'='*50}\n"
        f"Model: {args.model} | Hits reviewed: {len(reviews)}\n\n"
        f"Priority distribution: {priority_counts}\n\n"
    )
    for r in reviews:
        report += f"#{r['rank']} {r['compound']}: pIC50={r['pic50']:.2f} -> {r.get('priority', '?')}\n"
        report += f"  Flags: {r.get('liability_flags', [])}\n"
        report += f"  Rationale: {r.get('rationale', 'N/A')}\n"
        report += f"  Next step: {r.get('next_step', 'N/A')}\n\n"
    report += f"Tokens: in={total_input} out={total_output}\nCost: ${cost:.4f}\n"
    print(report)

    with open(os.path.join(args.output_dir, "hts_reviews.json"), "w") as f:
        json.dump(reviews, f, indent=2)
    with open(os.path.join(args.output_dir, "hts_review_report.txt"), "w") as f:
        f.write(report)
    print("Done.")


if __name__ == "__main__":
    main()
