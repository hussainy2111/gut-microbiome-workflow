
#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate genus-level composition and phylotype-level diversity results."
    )
    parser.add_argument("--shared", required=True)
    parser.add_argument("--taxonomy", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def read_shared_file(path):
    shared = pd.read_csv(path, sep=r"\s+")

    required = {"label", "Group"}
    if not required.issubset(shared.columns):
        raise ValueError("Shared file must contain label and Group columns.")

    return shared


def read_taxonomy_file(path):
    taxonomy = pd.read_csv(path, sep="\t")

    required = {"OTU", "Taxonomy"}
    if not required.issubset(taxonomy.columns):
        raise ValueError("Taxonomy file must contain OTU and Taxonomy columns.")

    return taxonomy


def genus_from_taxonomy(taxonomy_string):
    if pd.isna(taxonomy_string):
        return "Unclassified"

    ranks = str(taxonomy_string).split(";")

    # The sixth rank represents genus in this Mothur taxonomy hierarchy.
    if len(ranks) < 6:
        return "Unclassified"

    genus = ranks[5].split("(")[0].replace('"', "").strip()

    if (
        not genus
        or genus.lower() in {"unclassified", "unknown"}
        or genus.lower().endswith("_unclassified")
    ):
        return "Unclassified"

    return genus


def calculate_shannon(counts):
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()

    if total == 0:
        return 0.0

    proportions = counts[counts > 0] / total
    return float(-np.sum(proportions * np.log(proportions)))


def rarefied_shannon(counts, depth, rng, iterations=100):
    counts = np.asarray(counts, dtype=np.int64)

    if counts.sum() < depth:
        raise ValueError("A sample has fewer reads than the rarefaction depth.")

    estimates = []

    for _ in range(iterations):
        rarefied = rng.multivariate_hypergeometric(counts, depth)
        estimates.append(calculate_shannon(rarefied))

    return float(np.mean(estimates)), float(np.std(estimates, ddof=1))


def plot_diversity(diversity, output_path):
    fig, ax = plt.subplots(figsize=(8, 6))

    periods = ["early", "late"]

    for position, period in enumerate(periods):
        values = diversity.loc[
            diversity["period"] == period, "shannon_mean"
        ]

        if values.empty:
            continue

        x = np.full(len(values), position, dtype=float)
        x += np.linspace(-0.08, 0.08, len(values))

        ax.scatter(x, values, label=period, s=55)

        ax.hlines(
            values.mean(),
            position - 0.2,
            position + 0.2,
            linewidth=2
        )

    ax.set_xticks(range(len(periods)))
    ax.set_xticklabels(["Early", "Late"])
    ax.set_ylabel("Shannon diversity (phylotype level)")
    ax.set_xlabel("Sampling period")
    ax.set_title("Phylotype-level Shannon diversity")
    ax.legend(title="Sampling period")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    args = parse_arguments()

    output_dir = Path(args.outdir)
    output_dir.mkdir(parents=True, exist_ok=True)

    shared = read_shared_file(args.shared)
    taxonomy = read_taxonomy_file(args.taxonomy)
    metadata = pd.read_csv(args.metadata)

    required_metadata = {"sample_id", "day", "period"}
    if not required_metadata.issubset(metadata.columns):
        raise ValueError("Metadata must contain sample_id, day, and period.")

    excluded_columns = {"label", "Group", "numOtus", "numPhylos"}
    otu_columns = [
        column for column in shared.columns
        if column not in excluded_columns
    ]

    if not otu_columns:
        raise ValueError("No abundance feature columns were found.")

    if not set(otu_columns).issubset(set(taxonomy["OTU"])):
        missing = sorted(set(otu_columns) - set(taxonomy["OTU"]))
        raise ValueError(f"Shared-file features missing from taxonomy: {missing}")

    # Original sample-by-phylotype count table.
    phylotype_table = (
        shared.set_index("Group")[otu_columns]
        .apply(pd.to_numeric, errors="raise")
        .astype(int)
    )

    if (phylotype_table < 0).any().any():
        raise ValueError("Abundance counts cannot be negative.")

    phylotype_table.to_csv(
        output_dir / "phylotype_abundance.tsv",
        sep="\t"
    )

    # Map each phylotype to its genus for taxonomic composition analysis.
    taxonomy["Genus"] = taxonomy["Taxonomy"].apply(genus_from_taxonomy)
    genus_lookup = taxonomy.set_index("OTU")["Genus"].to_dict()

    records = []

    for sample_id, row in phylotype_table.iterrows():
        for otu, count in row.items():
            if count > 0:
                records.append({
                    "sample_id": sample_id,
                    "Genus": genus_lookup.get(otu, "Unclassified"),
                    "count": int(count),
                })

    abundance = pd.DataFrame(records)

    if abundance.empty:
        raise ValueError("No abundance records were found.")

    genus_table = (
        abundance.groupby(["sample_id", "Genus"])["count"]
        .sum()
        .unstack(fill_value=0)
    )

    genus_table.to_csv(
        output_dir / "genus_abundance.tsv",
        sep="\t"
    )

    # Calculate genus-level relative abundance.
    totals = genus_table.sum(axis=1)
    relative_abundance = genus_table.div(totals, axis=0)

    relative_abundance.to_csv(
        output_dir / "genus_relative_abundance.tsv",
        sep="\t"
    )

    # Plot the ten most abundant genera.
    top_genera = (
        genus_table.sum(axis=0)
        .sort_values(ascending=False)
        .head(10)
        .index
    )

    plot_data = relative_abundance[top_genera].copy()
    plot_data["Other"] = 1 - plot_data.sum(axis=1)

    ax = plot_data.plot(
        kind="bar",
        stacked=True,
        figsize=(12, 6),
        width=0.85
    )

    ax.set_ylabel("Relative abundance")
    ax.set_xlabel("Sample")
    ax.set_title("Genus-level community composition")
    ax.legend(
        title="Genus",
        bbox_to_anchor=(1.02, 1),
        loc="upper left"
    )

    plt.tight_layout()
    plt.savefig(
        output_dir / "genus_composition.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

    # Rarefy original phylotype counts to a common sequencing depth.
    sample_depths = phylotype_table.sum(axis=1)
    rarefaction_depth = int(sample_depths.min())

    if rarefaction_depth <= 0:
        raise ValueError("At least one sample has zero reads.")

    rng = np.random.default_rng(20261120)
    diversity_records = []

    for sample_id, counts in phylotype_table.iterrows():
        mean_shannon, sd_shannon = rarefied_shannon(
            counts.values,
            rarefaction_depth,
            rng
        )

        diversity_records.append({
            "sample_id": sample_id,
            "total_reads": int(counts.sum()),
            "rarefaction_depth": rarefaction_depth,
            "shannon_mean": mean_shannon,
            "shannon_sd": sd_shannon,
        })

    diversity = pd.DataFrame(diversity_records)

    diversity = diversity.merge(
        metadata,
        on="sample_id",
        how="left",
        validate="one_to_one"
    )

    if diversity["period"].isna().any():
        missing = diversity.loc[
            diversity["period"].isna(), "sample_id"
        ].tolist()
        raise ValueError(f"Missing metadata for samples: {missing}")

    diversity.to_csv(
        output_dir / "shannon_diversity.tsv",
        sep="\t",
        index=False
    )

    plot_diversity(
        diversity,
        output_dir / "shannon_comparison.png"
    )

    print(f"Samples analysed: {len(diversity)}")
    print(f"Rarefaction depth: {rarefaction_depth}")
    print("\nPhylotype-level Shannon diversity by period:")
    print(
        diversity.groupby("period")["shannon_mean"]
        .agg(["count", "mean", "std"])
    )

    print(f"\nResults saved to: {output_dir}")


if __name__ == "__main__":
    main()

