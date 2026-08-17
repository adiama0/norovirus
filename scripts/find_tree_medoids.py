#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
from Bio import Phylo, SeqIO


def find_tree_medoid(tree, selected_names, candidate_names):

    selected_names = set(str(name) for name in selected_names)
    candidate_names = set(str(name) for name in candidate_names)

    root = tree.root
    edge_length = {}
    order = []

    stack = [root]

    while stack:

        node = stack.pop()
        order.append(node)

        for child in node.clades:

            edge_length[child] = child.branch_length or 0.0
            stack.append(child)

    selected_count = {}
    distance_below = {}

    for node in reversed(order):

        if node.is_terminal():

            selected_count[node] = int(node.name in selected_names)
            distance_below[node] = 0.0

        else:

            selected_count[node] = 0
            distance_below[node] = 0.0

            for child in node.clades:

                selected_count[node] += selected_count[child]

                distance_below[node] += (
                    distance_below[child]
                    + selected_count[child] * edge_length[child]
                )

    number_selected = selected_count[root]

    total_distance = {
        root: distance_below[root]
    }

    for node in order:

        for child in node.clades:

            total_distance[child] = (
                total_distance[node]
                + (
                    number_selected
                    - 2 * selected_count[child]
                )
                * edge_length[child]
            )

    terminals = {
        terminal.name: terminal
        for terminal in tree.get_terminals()
    }

    candidates = [
        terminals[name]
        for name in sorted(candidate_names)
        if name in terminals
        and name in selected_names
    ]

    if not candidates:

        raise ValueError(
            "No full-genome candidates matched terminal tips in the tree."
        )

    medoid = min(
        candidates,
        key=lambda terminal: (
            total_distance[terminal],
            terminal.name
        )
    )

    if number_selected > 1:

        mean_distance = (
            total_distance[medoid]
            / (number_selected - 1)
        )

    else:
        mean_distance = 0.0

    return (
        medoid.name,
        number_selected,
        len(candidates),
        mean_distance,
    )


def main():

    parser = argparse.ArgumentParser(
        description="Find one full-genome VP1-tree medoid for each norovirus type."
    )

    parser.add_argument("--tree", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--sequences", required=True)
    parser.add_argument("--min-genome-length", type=int, default=7000)
    parser.add_argument("--group-dir", required=True)

    parser.add_argument(
        "--groups",
        nargs="+",
        required=True
    )

    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    tree = Phylo.read(args.tree, "newick")

    all_metadata = pd.read_csv(
        args.metadata,
        sep="\t",
        dtype=str
    )

    required_columns = {
        "name",
        "VP1_type"
    }

    if not required_columns.issubset(all_metadata.columns):

        raise ValueError(
            "Filtered metadata must contain name and VP1_type columns."
        )

    sequence_lengths = {
        record.id: len(record.seq)
        for record in SeqIO.parse(
            args.sequences,
            "fasta"
        )
    }

    results = []

    for group in args.groups:

        group_metadata_path = (
            Path(args.group_dir)
            / f"metadata_group_{group}.tsv"
        )

        group_metadata = pd.read_csv(
            group_metadata_path,
            sep="\t",
            dtype=str
        )

        group_names = (
            group_metadata["name"]
            .dropna()
            .astype(str)
            .tolist()
        )

        full_genome_names = [
            name
            for name in group_names
            if sequence_lengths.get(name, 0) >= args.min_genome_length
        ]

        if not full_genome_names:

            raise ValueError(
                f"No sequences >= {args.min_genome_length} nt "
                f"were found for {group}."
            )

        medoid_name, number_selected, number_candidates, mean_distance = (
            find_tree_medoid(
                tree,
                group_names,
                full_genome_names
            )
        )

        medoid_row = all_metadata[
            all_metadata["name"].astype(str) == str(medoid_name)
        ]

        if medoid_row.empty:

            vp1_type = "UNKNOWN"

        else:

            vp1_type = medoid_row.iloc[0]["VP1_type"]

        results.append(
            {
                "group": group,
                "medoid_accession": medoid_name,
                "VP1_type": vp1_type,
                "group_tree_tips": number_selected,
                "full_genome_candidates": number_candidates,
                "medoid_genome_length": sequence_lengths[medoid_name],
                "mean_patristic_distance": mean_distance,
            }
        )

    results = pd.DataFrame(results)

    results.to_csv(
        args.output,
        sep="\t",
        index=False
    )

    print()
    print(results.to_string(index=False))
    print()


if __name__ == "__main__":
    main()