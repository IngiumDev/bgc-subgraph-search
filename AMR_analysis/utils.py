from __future__ import annotations

from collections import Counter

from typing import Tuple, Optional

import sys

sys.path.insert(0, "../../operon-project/code/graph_api")
sys.path.insert(0, "../../operon-project/code/search_pipeline")

from graph_classes import Subgraph, ReferenceGraph, MergedGraph, Edge, Node
from graph_classes import nameddict, count_fg_occurence, read_meta_data, get_graph_id, get_subgraph_id, \
    has_existing_graph_builds, delete_file, load_modified_from_pickle_if_exists, reset_subgraph_and_delete_modified
from graph_classes import ABUNDANCE_DETECTION_LIMIT, TXT_POSTFIX, JSON_POSTFIX, PICKLE_POSTFIX, MODIFIED_POSTFIX
import pandas as pd
import seaborn as sns



import numpy as np
from pathlib import Path
import pickle
import matplotlib.pyplot as plt


def _load_single_subgraph_pkl(graph_folder: Path):
    pkls = sorted(graph_folder.glob("subgraph_*.pkl"))
    if len(pkls) != 1:
        raise RuntimeError(
            f"{graph_folder.name}: expected exactly 1 subgraph_*.pkl, found {len(pkls)} -> {[p.name for p in pkls]}")
    with pkls[0].open("rb") as f:
        return pickle.load(f)


def graph_taxonomy_summary(
        graph,  # expects something with .nodes: dict[id -> node] and node.taxonomy: dict
        tax_key: str = "gtdbtk_ann",
) -> Tuple[float, Optional[str]]:
    """
    Returns (percent_assigned_to_majority, majority_taxonomy_value or None).

    Rules:
    - Only count nodes where node.taxonomy has a non-empty value for tax_key.
    - Majority vote is computed over annotated nodes.
    - If there's a tie for top count -> error.
    - If no annotated nodes -> (0.0, None).
    - Percent is (top_count / total_nodes) * 100 (not / annotated_nodes).
    """
    nodes = list(getattr(graph, "nodes").values())
    total = len(nodes)

    if total == 0:
        return (0.0, None)

    vals = [
        (n.taxonomy.get(tax_key) if getattr(n, "taxonomy", None) else None)
        for n in nodes
    ]
    vals = [v for v in vals if v not in (None, "", {})]

    if not vals:
        return (0.0, None)

    counts = Counter(vals)
    (top_val, top_count) = counts.most_common(1)[0]

    # error if no strict majority among annotated nodes (tie for top)
    top_count_2 = counts.most_common(2)[1][1] if len(counts) > 1 else 0
    if top_count == top_count_2:
        raise RuntimeError(
            f"No unique majority for {tax_key}: {counts} in graph {graph.filepath}"
        )

    percent = (top_count / total) * 100.0
    return (percent, top_val)


def load_functional_groups_file(functional_groups_path):
    protein_id_to_element = {}
    with open(functional_groups_path, 'r') as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                protein_id, element_symbol = parts
                protein_id_to_element[protein_id] = element_symbol
    return protein_id_to_element
