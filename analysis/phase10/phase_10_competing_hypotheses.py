#!/usr/bin/env python3
"""Render the Phase 10 competing-hypotheses schematic with Graphviz DOT."""

from __future__ import annotations

import argparse
import html
import subprocess
from pathlib import Path


PALETTE = {
    "ink": "#263238",
    "muted": "#5F6B72",
    "main_line": "#7A858B",
    "ai_accent": "#238B45",
    "gov_accent": "#238B45",
    "hub_accent": "#B8860B",
    "observed": "#2F6FB0",
    "matrix": "#C81D25",
    "myeloid": "#7B3294",
    "challenge": "#D97706",
    "neutral": "#667085",
    "band1": "#EFF6FF",
    "band2": "#FFF7ED",
    "band3": "#F0FDF4",
    "band4": "#F8FAFC",
}


def card(heading: str, lines: list[str], accent: str, width: int = 238) -> str:
    text = "<BR/>".join(
        f'<FONT FACE="Arial" POINT-SIZE="10" COLOR="{PALETTE["muted"]}">{html.escape(line)}</FONT>'
        for line in lines
    )
    return (
        "<<TABLE BORDER=\"0\" CELLBORDER=\"0\" CELLSPACING=\"0\" CELLPADDING=\"0\">"
        f'<TR><TD BGCOLOR="{accent}" HEIGHT="6" WIDTH="{width}"></TD></TR>'
        f'<TR><TD BGCOLOR="white" FIXEDSIZE="true" WIDTH="{width}" HEIGHT="74">'
        f'<FONT FACE="Arial" POINT-SIZE="15"><B>{html.escape(heading)}</B></FONT><BR/>{text}'
        "</TD></TR></TABLE>>"
    )


def layer_label(number: str, title: str, color: str) -> str:
    return card(f"{number}  {title}", ["conceptual layer"], color, width=238)


def build_dot() -> str:
    graph_attrs = {
        "rankdir": "TB",
        "splines": "ortho",
        "outputorder": "edgesfirst",
        "compound": "true",
        "newrank": "true",
        "size": "16,9!",
        "ratio": "fill",
        "pad": "0.25",
        "margin": "0.04",
        "nodesep": "0.42",
        "ranksep": "0.62",
        "bgcolor": "white",
        "fontname": "Arial",
        "labelloc": "t",
        "label": "Phase 10 competing hypotheses and falsification map (schematic)",
        "fontsize": "22",
        "fontcolor": PALETTE["ink"],
    }
    lines = ["digraph Phase10 {", "  graph [" + ", ".join(f'{k}=\"{v}\"' for k, v in graph_attrs.items()) + "];", "  node [shape=plain, fontname=\"Arial\"];", "  edge [fontname=\"Arial\", fontsize=\"9\", color=\"#7A858B\", fontcolor=\"#5F6B72\", penwidth=\"1.5\", arrowsize=\"0.75\"];" ]

    # Nodes by conceptual layer.
    bands = [
        ("observations", PALETTE["band1"], [
            ("l1", layer_label("01", "Observed evidence", PALETTE["observed"])),
            ("axis", card("Day7 spatial continuum", ["3/3 slices", "fixed Phase 9 axis • Grade A"], PALETTE["observed"])),
            ("genes", card("Non-circular signature", ["ECM/myofibroblast ↔ contractile", "200 displayed genes • Grade A"], PALETTE["matrix"])),
            ("composition", card("Axis-extreme composition", ["Fibroblast + Trem2/Spp1 rise", "Cardiomyocyte falls • Grade A"], PALETTE["myeloid"])),
        ]),
        ("hypotheses", PALETTE["band2"], [
            ("l2", layer_label("02", "Competing explanations", PALETTE["challenge"])),
            ("h1", card("H1 Matrix-remodeling state", ["activated fibroblast program", "replaces viable myocardium • A/B"], PALETTE["matrix"])),
            ("h2", card("H2 Myeloid upstream driver", ["Trem2/Spp1 macrophage niche", "direction remains unproven • B"], PALETTE["myeloid"])),
            ("h3", card("H3 Technical/injury confound", ["RNA quality or necrotic burden", "not fully excluded • B/D"], PALETTE["challenge"])),
            ("h4", card("H4 Generic wound response", ["shared across injury models", "not MI-specific • B"], PALETTE["neutral"])),
        ]),
        ("tests", PALETTE["band3"], [
            ("l3", layer_label("03", "Discriminating evidence", PALETTE["ai_accent"])),
            ("pb", card("snRNA pseudobulk", ["17 sample_ids • 11 classes", "sample-level effects / LOSO"], PALETTE["observed"])),
            ("external", card("GSE280373 projection", ["Postn/Thbs4 fibroblasts high", "descriptive; no animal IDs"], PALETTE["ai_accent"])),
            ("histology", card("Multiplex spatial validation", ["Cthrc1/Serpinh1 + Spp1/Trem2", "distance to viable myocardium"], PALETTE["hub_accent"])),
            ("perturb", card("Priority perturbation", ["block candidate myeloid→fibroblast axis", "measure scar + function"], PALETTE["hub_accent"])),
        ]),
        ("conclusion", PALETTE["band4"], [
            ("l4", layer_label("04", "Current boundary", PALETTE["neutral"])),
            ("model", card("Supported working model", ["fibrotic-to-contractile continuum", "multi-cell spatial association"], PALETTE["matrix"])),
            ("boundary", card("Not yet established", ["no causal cell direction", "no MI specificity or intervention proof"], PALETTE["challenge"])),
            ("next", card("Fastest decisive experiment", ["multiplex RNAscope/IF first", "then targeted perturbation"], PALETTE["ai_accent"])),
        ]),
    ]
    for band_name, fill, nodes in bands:
        lines.extend([
            f"  subgraph cluster_{band_name} {{",
            f'    label=""; style="rounded,filled"; color="{fill}"; fillcolor="{fill}"; penwidth="1.0"; margin="14";',
        ])
        for node_id, label in nodes:
            lines.append(f"    {node_id} [label={label}];")
        lines.append("  }")

    # Every conceptual row has an explicit rank=same constraint.
    row_nodes = [
        ["l1", "axis", "genes", "composition"],
        ["l2", "h1", "h2", "h3", "h4"],
        ["l3", "pb", "external", "histology", "perturb"],
        ["l4", "model", "boundary", "next"],
    ]
    for idx, row in enumerate(row_nodes, start=1):
        lines.append(f"  subgraph rank_row_{idx} {{ rank=same; " + "; ".join(row) + "; }")

    # Invisible ordering and vertical alignment edges.
    for row in row_nodes:
        for left, right in zip(row, row[1:]):
            lines.append(f'  {left} -> {right} [style="invis", weight="100"];')
    for top, bottom in [("l1", "l2"), ("l2", "l3"), ("l3", "l4"), ("axis", "h1"), ("genes", "h2"), ("composition", "h3"), ("h1", "pb"), ("h2", "external"), ("h3", "histology"), ("pb", "model"), ("external", "boundary"), ("histology", "next")]:
        lines.append(f'  {top} -> {bottom} [style="invis", weight="45"];')

    # Primary association flow edges.
    main = PALETTE["main_line"]
    for source, target, label in [
        ("axis", "h1", "supports"), ("genes", "h1", "supports"),
        ("composition", "h1", "supports"), ("composition", "h2", "supports"),
        ("axis", "h3", "challenge"), ("genes", "h4", "alternative"),
        ("h1", "model", "best-supported"), ("h2", "boundary", "direction unknown"),
        ("h3", "boundary", "residual risk"), ("h4", "boundary", "specificity unknown"),
    ]:
        lines.append(f'  {source} -> {target} [color="{main}", xlabel="{label}"];')

    # Evidence-integration arrows.
    green = PALETTE["ai_accent"]
    for source, target, label in [
        ("pb", "h1", "cell-source test"), ("external", "h1", "external support"),
        ("external", "h4", "model challenge"), ("histology", "next", "localization"),
        ("perturb", "next", "causality"), ("model", "next", "prioritize"),
    ]:
        lines.append(f'  {source} -> {target} [color="{green}", penwidth="1.7", xlabel="{label}"];')

    # Feedback/falsification loops are dashed, green, and never constrain layout.
    for source, target, label in [
        ("histology", "h1", "falsify localization"),
        ("histology", "h2", "falsify adjacency"),
        ("perturb", "h2", "falsify direction"),
        ("next", "boundary", "recalibrate claim"),
    ]:
        lines.append(f'  {source} -> {target} [color="{PALETTE["gov_accent"]}", style="dashed", constraint="false", penwidth="1.4", xlabel="{label}"];')

    lines.append("}")
    return "\n".join(lines) + "\n"


def render(output_stem: Path, dot_bin: str) -> None:
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    dot_path = output_stem.with_suffix(".dot")
    dot_text = build_dot()
    dot_path.write_text(dot_text, encoding="utf-8")
    for fmt in ("svg", "pdf", "png"):
        target = output_stem.with_suffix(f".{fmt}")
        subprocess.run([dot_bin, f"-T{fmt}", str(dot_path), "-o", str(target)], check=True)
        if not target.exists() or target.stat().st_size < 10_000:
            raise RuntimeError(f"Graphviz output is absent or too small: {target}")

    required = ["rankdir=\"TB\"", "splines=\"ortho\"", "outputorder=\"edgesfirst\"", "compound=\"true\"", "newrank=\"true\""]
    for token in required:
        if token not in dot_text:
            raise RuntimeError(f"Missing Graphviz invariant: {token}")
    if dot_text.count("rank=same") != 4:
        raise RuntimeError("Every conceptual row must have rank=same")
    for line in dot_text.splitlines():
        if "style=\"dashed\"" in line and "constraint=\"false\"" not in line:
            raise RuntimeError("A feedback edge constrains layout")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_stem", type=Path)
    parser.add_argument("--dot-bin", default="/opt/homebrew/bin/dot")
    args = parser.parse_args()
    render(args.output_stem, args.dot_bin)


if __name__ == "__main__":
    main()
