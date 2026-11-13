"""Select explicitly listed tissues or samples without dataset-specific branches."""

import csv
from pathlib import Path


def add_selection_args(parser):
    parser.add_argument(
        "--tissues", type=Path, help="TSV with a tissue column; omit for all tissues"
    )
    parser.add_argument(
        "--samples", type=Path, help="TSV with an ont_sample column; omit for all samples"
    )


def listed(path, column):
    with path.open(encoding="utf8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if column not in (reader.fieldnames or []):
            raise ValueError(f"Manifest requires {column}")
        values = [row[column] for row in reader]
    if not values or any(not x for x in values) or len(values) != len(set(values)):
        raise ValueError("Selection must be nonempty and unique")
    return values


def select_frame(frame, args):
    for option, column in [("tissues", "tissue"), ("samples", "ont_sample")]:
        path = getattr(args, option, None)
        if path is None:
            continue
        values = listed(path, column)
        absent = set(values) - set(frame[column])
        if absent:
            raise ValueError(f"Unknown {column} values: {sorted(absent)}")
        frame = frame[frame[column].isin(values)].copy()
    if frame.empty:
        raise ValueError("Selection contains no samples")
    return frame
