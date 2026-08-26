import argparse
import csv
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--proportions", required=True)
parser.add_argument("--bootstraps", required=True)
parser.add_argument("--sample", required=True)
parser.add_argument("--output", required=True)

parser.add_argument(
    "--detection-threshold",
    type=float,
    default=0.0
)

args = parser.parse_args()

# Read original competitive-mapping estimate
original_row = None

with open(args.proportions) as handle:
    reader = csv.DictReader(
        handle,
        delimiter="\t"
    )

    for row in reader:
        if row["sample"] == args.sample:
            original_row = row
            break


if original_row is None:
    raise ValueError(
        f"Could not find {args.sample} "
        f"in {args.proportions}"
    )

# Read bootstrap proportions
with open(args.bootstraps) as handle:
    reader = csv.DictReader(
        handle,
        delimiter="\t"
    )

    vp1_types = [
        column
        for column in reader.fieldnames
        if column != "replicate"
    ]

    bootstrap_rows = list(reader)


if len(bootstrap_rows) == 0:
    raise ValueError(
        "No bootstrap replicates were found."
    )

# Calculate confidence statistics
with open(args.output, "w", newline="") as handle:
    writer = csv.writer(
        handle,
        delimiter="\t"
    )

    writer.writerow(
        [
            "VP1_type",
            "Estimate",
            "95% CI",
            "Bootstrap support"
        ]
    )

    for vp1_type in vp1_types:
        estimate = float(
            original_row[vp1_type]
        )

        bootstrap_values = np.array(
            [
                float(row[vp1_type])
                for row in bootstrap_rows
            ],
            dtype=float
        )

        lower_ci, upper_ci = np.percentile(
            bootstrap_values,
            [2.5, 97.5]
        )

        support = (
            np.mean(
                bootstrap_values
                >= args.detection_threshold
            )
            * 100
        )

        # Skip VP1 types never detected
        if estimate == 0 and support == 0:
            continue

        confidence_interval = (
            f"{lower_ci:.6f}-{upper_ci:.6f}"
        )

        bootstrap_support = (
            f"{support:.0f}%"
        )

        writer.writerow(
            [
                vp1_type,
                f"{estimate:.6f}",
                confidence_interval,
                bootstrap_support
            ]
        )