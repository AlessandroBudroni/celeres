#!/usr/bin/env python3
"""Print final average signature-size estimates in KB, with KB = bytes / 1024.

The formulas mirror the compact serialized layouts implemented in C. Averages
use the exact expected number of published seed-tree nodes for a uniform
fixed-weight challenge.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from decimal import Decimal, ROUND_CEILING


Q = 127


@dataclass(frozen=True)
class Level:
    name: str
    category: int
    target: int
    n: int
    k: int
    seed_bytes: int
    rounds: int
    weight: int
    max_published_seeds: int
    perm_bits: int | None
    tree_offsets: tuple[int, ...]
    tree_nodes_per_level: tuple[int, ...]
    tree_leaves_start_indices: tuple[int, ...]
    tree_consecutive_leaves: tuple[int, ...]


CUSTOM_LEVELS = (
    Level(
        "l1", 272, 192, 272, 136, 16, 192, 36, 87, 1819,
        (0, 0, 0, 0, 0, 0, 0, 0, 128),
        (1, 2, 4, 8, 16, 32, 64, 128, 128),
        (255, 191),
        (128, 64),
    ),
    Level(
        "l3", 416, 220, 416, 208, 24, 220, 68, 119, 3035,
        (0, 0, 0, 0, 0, 0, 0, 8, 56),
        (1, 2, 4, 8, 16, 32, 64, 120, 192),
        (247, 223, 123),
        (192, 24, 4),
    ),
    Level(
        "l5", 564, 345, 564, 282, 32, 345, 75, 169, 4364,
        (0, 0, 0, 0, 0, 2, 2, 2, 50, 178),
        (1, 2, 4, 8, 16, 30, 60, 120, 192, 256),
        (433, 369, 217, 30),
        (256, 64, 24, 1),
    ),
)


def ceil_div(x: int, y: int) -> int:
    return (x + y - 1) // y


def bits_to_represent(x: int) -> int:
    return max(1, x.bit_length())


def ceil_log2(x: int) -> int:
    if x < 1:
        raise ValueError("ring sizes must be positive")
    return (x - 1).bit_length()


def parent(i: int) -> int:
    return (i - 1) // 2 if i % 2 else (i - 2) // 2


def sibling(i: int) -> int:
    return i + 1 if i % 2 else i - 1


def rref_mat_packed_bytes(level: Level) -> int:
    data_bits = bits_to_represent(Q) * (level.n - level.k) * level.k
    return ceil_div(data_bits, 8) + ceil_div(level.n, 8)


def celeres_self_orthogonal_tag_bytes(level: Level) -> int:
    pivot_bytes = ceil_div(level.n, 8)
    order_bits = level.k * ceil_log2(level.k)
    matrix_elements = (level.k * (level.k - 1)) // 2
    root_bits = level.k
    antiorthogonal_bits = bits_to_represent(Q) * matrix_elements
    return (pivot_bytes + ceil_div(order_bits, 8) +
            ceil_div(root_bits, 8) + ceil_div(antiorthogonal_bits, 8))


def celeres_tag_bytes(level: Level, tag_mode: str) -> int:
    if tag_mode == "rref":
        return rref_mat_packed_bytes(level)
    if tag_mode == "self_orthogonal":
        return celeres_self_orthogonal_tag_bytes(level)
    raise ValueError(f"unsupported CELERES tag mode: {tag_mode}")


def seed_count_for_challenge(level: Level, challenge: list[int]) -> int:
    """Mirror GGMPath(): challenge bit 0 means release this leaf seed."""

    not_to_publish = 1
    to_publish = 0
    flags = [not_to_publish] * ((2 * level.rounds) - 1)

    challenge_index = 0
    for start, count in zip(level.tree_leaves_start_indices, level.tree_consecutive_leaves):
        for j in range(count):
            flags[start + j] = challenge[challenge_index]
            challenge_index += 1

    start_node = level.tree_leaves_start_indices[0]
    for level_idx in range(ceil_log2(level.rounds), 0, -1):
        for i in range(level.tree_nodes_per_level[level_idx] - 2, -1, -2):
            current_node = start_node + i
            parent_node = parent(current_node) + (level.tree_offsets[level_idx - 1] >> 1)
            if flags[current_node] == to_publish and flags[sibling(current_node)] == to_publish:
                flags[parent_node] = to_publish
            else:
                flags[parent_node] = not_to_publish
        start_node -= level.tree_nodes_per_level[level_idx - 1]

    start_node = 1
    seed_count = 0
    for level_idx in range(1, ceil_log2(level.rounds) + 1):
        for node_in_level in range(level.tree_nodes_per_level[level_idx]):
            current_node = start_node + node_in_level
            father_node = parent(current_node) + (level.tree_offsets[level_idx - 1] >> 1)
            if flags[current_node] == to_publish and flags[father_node] == not_to_publish:
                seed_count += 1
        start_node += level.tree_nodes_per_level[level_idx]
    return seed_count


def comb_ratio_no_opened_leaves(total_leaves: int, subtree_leaves: int, weight: int) -> float:
    """Probability that a fixed subtree contains no challenge-1 leaves."""

    if subtree_leaves <= 0:
        return 1.0
    if total_leaves - subtree_leaves < weight:
        return 0.0
    probability = 1.0
    for i in range(weight):
        probability *= (total_leaves - subtree_leaves - i) / (total_leaves - i)
    return probability


def exact_expected_seed_count(level: Level) -> float:
    """Expected GGMPath() seed count for a uniform fixed-weight challenge."""

    node_count = sum(level.tree_nodes_per_level)
    leaf_counts = [0] * node_count
    for start, count in zip(level.tree_leaves_start_indices,
                            level.tree_consecutive_leaves):
        for j in range(count):
            leaf_counts[start + j] = 1

    start_node = level.tree_leaves_start_indices[0]
    for level_idx in range(ceil_log2(level.rounds), 0, -1):
        for i in range(level.tree_nodes_per_level[level_idx] - 2, -1, -2):
            current_node = start_node + i
            sibling_node = sibling(current_node)
            parent_node = parent(current_node) + (
                level.tree_offsets[level_idx - 1] >> 1
            )
            leaf_counts[parent_node] = (
                leaf_counts[current_node] + leaf_counts[sibling_node]
            )
        start_node -= level.tree_nodes_per_level[level_idx - 1]

    expected = 0.0
    start_node = 1
    for level_idx in range(1, ceil_log2(level.rounds) + 1):
        for node_in_level in range(level.tree_nodes_per_level[level_idx]):
            current_node = start_node + node_in_level
            father_node = parent(current_node) + (
                level.tree_offsets[level_idx - 1] >> 1
            )
            node_leaves = leaf_counts[current_node]
            parent_leaves = leaf_counts[father_node]
            expected += (
                comb_ratio_no_opened_leaves(level.rounds, node_leaves,
                                            level.weight)
                - comb_ratio_no_opened_leaves(level.rounds, parent_leaves,
                                              level.weight)
            )
        start_node += level.tree_nodes_per_level[level_idx]
    return expected


def celeres_size(level: Level,
                 ring_size: int,
                 tag_mode: str) -> tuple[int, dict[str, int]]:
    if level.perm_bits is None:
        raise ValueError(f"no CELERES permutation-compression parameters for {level.name}")
    hash_bytes = 2 * level.seed_bytes
    tag_bytes = celeres_tag_bytes(level, tag_mode)
    perm_bytes = ceil_div(level.perm_bits, 8)
    path_bytes = ceil_log2(ring_size) * hash_bytes
    base_bytes = tag_bytes + (2 * hash_bytes) + 1
    seed_tree_bytes = level.max_published_seeds * level.seed_bytes
    response_bytes = perm_bytes + level.seed_bytes + path_bytes
    total = base_bytes + seed_tree_bytes + (level.weight * response_bytes)
    return total, {
        "base_bytes": base_bytes,
        "tag_bytes": tag_bytes,
        "seed_tree_bytes": seed_tree_bytes,
        "response_bytes": response_bytes,
        "perm_bytes": perm_bytes,
        "path_bytes": path_bytes,
    }


def parse_ring_sizes(value: str) -> list[int]:
    ring_sizes = []
    for item in value.split(","):
        item = item.strip()
        if not item:
            continue
        ring_size = int(item, 0)
        if ring_size <= 0:
            raise argparse.ArgumentTypeError("ring sizes must be positive")
        ring_sizes.append(ring_size)
    if not ring_sizes:
        raise argparse.ArgumentTypeError("at least one ring size is required")
    return ring_sizes


def expected_signature_bytes(base_bytes: int,
                             response_bytes: int,
                             level: Level,
                             expected_seeds: float) -> float:
    return (
        base_bytes +
        (level.weight * response_bytes) +
        (expected_seeds * level.seed_bytes)
    )


def celeres_expected_bytes(level: Level,
                           ring_size: int,
                           tag_mode: str,
                           expected_seeds: float) -> float:
    _, parts = celeres_size(level, ring_size, tag_mode)
    return expected_signature_bytes(
        parts["base_bytes"], parts["response_bytes"], level, expected_seeds
    )


def kb(value_bytes: float) -> str:
    value_kb = Decimal(str(value_bytes)) / Decimal(1024)
    return str(value_kb.quantize(Decimal("0.01"), rounding=ROUND_CEILING))


def table_rows(ring_sizes: list[int]) -> list[dict[str, int | str]]:
    rows: list[dict[str, int | str]] = []
    for level in CUSTOM_LEVELS:
        expected_seeds = exact_expected_seed_count(level)
        for ring_size in ring_sizes:
            rows.append({
                "level": level.name,
                "n": level.n,
                "k": level.k,
                "t": level.rounds,
                "w": level.weight,
                "ring_size": ring_size,
                "celeres_no_so_avg_kb": kb(celeres_expected_bytes(
                    level, ring_size, "rref", expected_seeds
                )),
                "celeres_so_avg_kb": kb(celeres_expected_bytes(
                    level, ring_size, "self_orthogonal", expected_seeds
                )),
            })
    return rows


def print_markdown(rows: list[dict[str, int | str]]) -> None:
    headers = (
        "level",
        "n",
        "k",
        "t",
        "w",
        "ring_size",
        "celeres_no_so_avg_kb",
        "celeres_so_avg_kb",
    )
    print("| " + " | ".join(headers) + " |")
    print("| " + " | ".join("---" for _ in headers) + " |")
    for row in rows:
        print("| " + " | ".join(str(row[h]) for h in headers) + " |")


def print_csv(rows: list[dict[str, int | str]]) -> None:
    fieldnames = list(rows[0].keys()) if rows else []
    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Print final average signature sizes in KB, with KB = bytes / 1024."
    )
    parser.add_argument(
        "-r",
        "--ring-sizes",
        type=parse_ring_sizes,
        default=parse_ring_sizes("2,4,8,16,32,64,128"),
        help="comma-separated ring sizes; each value is r, not log2(r)",
    )
    parser.add_argument(
        "--format",
        choices=("markdown", "csv"),
        default="markdown",
        help="output format",
    )
    args = parser.parse_args()

    rows = table_rows(args.ring_sizes)
    if args.format == "csv":
        print_csv(rows)
    else:
        print_markdown(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
