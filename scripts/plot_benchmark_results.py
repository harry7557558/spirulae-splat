#!/usr/bin/env python3
"""
Fetch a raw markdown benchmark report and turn it into a multi-page PDF.

Changes requested:
- Maximum one page per section.
- Any data points with the same metric are plotted on the same subplot.
- Each subplot uses a cartesian layout:
    x-axis = datasets (categorical)
    y-axis = metric value
    scatter points = presets for that dataset/metric
- Parse tables separated by multiple spaces instead of tabs.
- Fetch markdown using urllib.

Usage:
    python plot_benchmark_markdown.py "https://example.com/report.md" -o report.pdf
"""

from __future__ import annotations

import argparse
import math
import re
import sys
import urllib.request
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


HR_RE = re.compile(r"^\s*(?:-{3,}|\*{3,}|_{3,})\s*$")
H2_RE = re.compile(r"^\s*##\s+(.*)$")
FENCE_RE = re.compile(r"^\s*```(.*)$")
NUM_RE = re.compile(r"[-+]?(?:\d[\d,]*\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")

DATASET_KEYS = {
    "dataset", "scene", "benchmark", "bench", "task", "split", "domain",
    "dataset_name", "name", "testset", "evalset"
}
PRESET_KEYS = {
    "preset", "method", "model", "variant", "setting", "approach",
    "config", "policy", "recipe", "baseline"
}
META_KEYS = {
    "id", "idx", "index", "epoch", "step", "iter", "iteration",
    "seed", "trial", "rank", "run"
}

# Hard-coded metric direction hints. Extend as needed.
LOWER_IS_BETTER = {
    "lpips", "l1", "l2", "mae", "mse", "rmse", "error", "loss", "time",
    "latency", "duration", "runtime", "wer", "cer", "edit", "perplexity",
    "kl", "kld", "ade", "fde", "ate", "rte", "aae", "mota", "cost", "distance"
}
HIGHER_IS_BETTER = {
    "psnr", "ssim", "accuracy", "acc", "auc", "iou", "f1", "precision",
    "recall", "ap", "map", "bleu", "rouge", "chrf", "mrr", "ndcg",
    "dice", "fps", "throughput", "success", "sr", "pass@1"
}


@dataclass
class Observation:
    section_index: int
    section_title: str
    dataset: str
    preset: str
    metric: str
    value: float
    source_heading: str
    source_block_index: int


def fetch_markdown(url: str) -> str:
    with urllib.request.urlopen(url) as resp:
        raw = resp.read()
    return raw.decode("utf-8", errors="replace")


def split_sections(md: str) -> List[str]:
    sections: List[str] = []
    buf: List[str] = []
    for line in md.splitlines():
        if HR_RE.match(line):
            chunk = "\n".join(buf).strip("\n")
            if chunk.strip():
                sections.append(chunk)
            buf = []
        else:
            buf.append(line)
    chunk = "\n".join(buf).strip("\n")
    if chunk.strip():
        sections.append(chunk)
    return sections


def extract_section_title(section_text: str, fallback: str) -> str:
    for line in section_text.splitlines():
        m = H2_RE.match(line)
        if m:
            return m.group(1).strip()
    return fallback


def parse_h2_hints(heading: str) -> Tuple[Optional[str], Optional[str]]:
    tokens = re.findall(r"`([^`]+)`", heading)
    dataset = tokens[0].strip() if len(tokens) >= 1 else None
    preset = tokens[1].strip() if len(tokens) >= 2 else None
    return dataset, preset


def extract_code_blocks_with_context(section_text: str) -> List[Tuple[str, str, str]]:
    """
    Returns (preceding_h2_heading, fence_info, code_text)
    """
    blocks: List[Tuple[str, str, str]] = []
    last_h2 = ""
    in_code = False
    fence_info = ""
    code_lines: List[str] = []

    for line in section_text.splitlines():
        h2 = H2_RE.match(line)
        if h2 and not in_code:
            last_h2 = h2.group(1).strip()

        fence = FENCE_RE.match(line)
        if fence:
            if not in_code:
                in_code = True
                fence_info = fence.group(1).strip()
                code_lines = []
            else:
                in_code = False
                blocks.append((last_h2, fence_info, "\n".join(code_lines)))
            continue

        if in_code:
            code_lines.append(line)

    return blocks


def split_row_by_pipe(line: str) -> List[str]:
    return [cell.strip() for cell in line.strip().strip("|").split("|")]


def split_row_by_tabs(line: str) -> List[str]:
    return [cell.strip() for cell in line.rstrip("\n").split("\t")]


def split_row_by_spaces(line: str) -> List[str]:
    # Two or more spaces are treated as the delimiter.
    # This is what the user requested for the common "aligned columns" format.
    return [cell.strip() for cell in re.split(r"\s{2,}", line.strip()) if cell.strip()]


def is_separator_row(cells: Sequence[str]) -> bool:
    if not cells:
        return False
    for c in cells:
        c = c.strip()
        if not c:
            continue
        if not re.fullmatch(r"[:\- ]{3,}", c):
            return False
    return True


def parse_number(text: object) -> Optional[float]:
    if text is None:
        return None
    s = str(text).strip()
    if not s:
        return None
    if s.lower() in {"-", "—", "–", "n/a", "na", "nan", "none", "null"}:
        return None
    s = s.replace(",", "")
    m = NUM_RE.search(s)
    if not m:
        return None
    try:
        return float(m.group(0))
    except ValueError:
        return None


def normalize_table_lines(block_text: str) -> List[str]:
    return [ln.rstrip("\n") for ln in block_text.splitlines() if ln.strip()]


def parse_table_block(block_text: str) -> Optional[Tuple[List[str], List[List[str]]]]:
    """
    Try to parse a markdown-ish table from a fenced code block.

    Supports:
    - pipe-separated markdown tables
    - tab-separated tables
    - space-separated tables using 2+ spaces as delimiters
    """
    lines = normalize_table_lines(block_text)
    if not lines:
        return None

    # 1) Pipe table
    pipe_lines = [ln for ln in lines if "|" in ln]
    if len(pipe_lines) >= 2:
        rows = [split_row_by_pipe(ln) for ln in pipe_lines]
        if len(rows) >= 2 and is_separator_row(rows[1]):
            header = rows[0]
            data_rows = rows[2:]
        else:
            header = rows[0]
            data_rows = rows[1:]
        if len(header) >= 2 and data_rows:
            return header, data_rows

    # 2) Tab table
    tab_lines = [ln for ln in lines if "\t" in ln]
    if len(tab_lines) >= 2:
        rows = [split_row_by_tabs(ln) for ln in tab_lines]
        header = rows[0]
        data_rows = rows[1:]
        if len(header) >= 2 and data_rows:
            return header, data_rows

    # 3) Multiple-space table
    # Accept if most lines split into at least 2 columns.
    space_rows = [split_row_by_spaces(ln) for ln in lines]
    valid_rows = [r for r in space_rows if len(r) >= 2]
    if len(valid_rows) >= 2 and len(valid_rows) >= max(2, len(lines) // 2):
        header = valid_rows[0]
        data_rows = valid_rows[1:]
        if len(header) >= 2 and data_rows:
            return header, data_rows

    return None


def header_matches_any(header: str, keywords: Iterable[str]) -> bool:
    h = header.lower().strip()
    return any(k in h for k in keywords)


def metric_is_higher_better(metric_name: str) -> bool:
    name = metric_name.lower()
    if any(tok in name for tok in LOWER_IS_BETTER):
        return False
    if any(tok in name for tok in HIGHER_IS_BETTER):
        return True
    return None


def infer_row_roles(headers: List[str], rows: List[List[str]]) -> Dict[str, object]:
    """
    Infer dataset/preset columns and metric columns.
    """
    headers = [(h.strip() if h is not None else "") for h in headers]
    ncols = len(headers)

    padded_rows = [r + [""] * max(0, ncols - len(r)) for r in rows]

    numeric_count = [0] * ncols
    nonempty_count = [0] * ncols

    for r in padded_rows:
        for i in range(ncols):
            cell = r[i].strip() if i < len(r) else ""
            if cell:
                nonempty_count[i] += 1
                if parse_number(cell) is not None:
                    numeric_count[i] += 1

    dataset_col = None
    preset_col = None
    categorical_cols: List[int] = []

    for i, h in enumerate(headers):
        hl = h.lower()
        if header_matches_any(hl, DATASET_KEYS):
            dataset_col = i
            categorical_cols.append(i)
        elif header_matches_any(hl, PRESET_KEYS):
            preset_col = i
            categorical_cols.append(i)
        elif header_matches_any(hl, META_KEYS):
            categorical_cols.append(i)

    for i in range(ncols):
        if i in categorical_cols:
            continue
        if nonempty_count[i] == 0:
            continue
        ratio = numeric_count[i] / nonempty_count[i]
        if ratio < 0.5:
            categorical_cols.append(i)

    categorical_cols = sorted(set(categorical_cols))

    metric_cols = [i for i in range(ncols) if i not in categorical_cols and numeric_count[i] > 0]

    return {
        "dataset_col": dataset_col,
        "preset_col": preset_col,
        "metric_cols": metric_cols,
        "categorical_cols": categorical_cols,
    }


def choose_dataset_and_preset(
    row: List[str],
    headers: List[str],
    roles: Dict[str, object],
    section_dataset_hint: Optional[str],
    preset_hint: Optional[str],
) -> Tuple[str, str]:
    headers = [(h.strip() if h is not None else "") for h in headers]
    row = row + [""] * max(0, len(headers) - len(row))

    dataset_col = roles["dataset_col"]
    preset_col = roles["preset_col"]
    cat_cols = roles["categorical_cols"] or []

    dataset = None
    preset = None

    if dataset_col is not None and dataset_col < len(row):
        dataset = row[dataset_col].strip() or None
    if preset_col is not None and preset_col < len(row):
        preset = row[preset_col].strip() or None

    if dataset is None and preset is None and cat_cols:
        if len(cat_cols) == 1:
            value = row[cat_cols[0]].strip() if cat_cols[0] < len(row) else ""
            if section_dataset_hint and preset_hint:
                dataset = value
                preset = preset_hint
            elif section_dataset_hint:
                dataset = section_dataset_hint
                preset = value
            elif preset_hint:
                dataset = value
                preset = preset_hint
            else:
                dataset = value
        else:
            dataset = row[cat_cols[0]].strip() if cat_cols[0] < len(row) else None
            preset = row[cat_cols[1]].strip() if cat_cols[1] < len(row) else None

    if dataset is None:
        dataset = section_dataset_hint or "unknown_dataset"
    if preset is None:
        preset = preset_hint or "unknown_preset"
    assert len(set(list(preset))) > 1

    return (dataset.strip() if dataset else "unknown_dataset",
            preset.strip() if preset else "unknown_preset")


def parse_table_as_observations(
    headers: List[str],
    rows: List[List[str]],
    section_index: int,
    section_title: str,
    heading_text: str,
    block_index: int,
) -> List[Observation]:
    roles = infer_row_roles(headers, rows)
    metric_cols: List[int] = roles["metric_cols"] or []

    if not metric_cols:
        return []

    section_dataset_hint, preset_hint = parse_h2_hints(heading_text)
    observations: List[Observation] = []

    for row in rows:
        row = row + [""] * max(0, len(headers) - len(row))
        dataset, preset = choose_dataset_and_preset(
            row=row,
            headers=headers,
            roles=roles,
            section_dataset_hint=section_dataset_hint,
            preset_hint=preset_hint,
        )

        for idx in metric_cols:
            if idx >= len(row):
                continue
            metric_name = headers[idx].strip() or f"metric_{idx}"
            value = parse_number(row[idx])
            if value is None:
                continue

            observations.append(
                Observation(
                    section_index=section_index,
                    section_title=section_title,
                    dataset=dataset,
                    preset=preset,
                    metric=metric_name,
                    value=value,
                    source_heading=heading_text,
                    source_block_index=block_index,
                )
            )

    return observations


def transpose_table(headers: List[str], rows: List[List[str]]) -> Tuple[List[str], List[List[str]]]:
    """
    Transpose a ragged table safely.
    """
    width = max([len(headers)] + [len(r) for r in rows]) if rows else len(headers)
    padded = [headers + [""] * max(0, width - len(headers))]
    for r in rows:
        padded.append(r + [""] * max(0, width - len(r)))

    transposed = list(map(list, zip(*padded)))
    new_headers = transposed[0]
    new_rows = transposed[1:]
    return new_headers, new_rows


def table_block_to_observations(
    block_text: str,
    section_index: int,
    section_title: str,
    heading_text: str,
    block_index: int,
) -> List[Observation]:
    parsed = parse_table_block(block_text)
    if parsed is None:
        return []

    headers, rows = parsed
    rows = [r for r in rows if set(list(''.join(r))) != {'-'}]
    if len(headers) != len(rows[0]):
        headers = [''] + headers

    t_headers, t_rows = transpose_table(headers, rows)
    return parse_table_as_observations(
        headers=t_headers,
        rows=t_rows,
        section_index=section_index,
        section_title=section_title,
        heading_text=heading_text,
        block_index=block_index,
    )


def collect_observations(md: str) -> List[Observation]:
    all_obs: List[Observation] = []
    sections = split_sections(md)

    for sec_idx, section_text in enumerate(sections, start=1):
        section_title = extract_section_title(section_text, fallback=f"Section {sec_idx}")
        blocks = extract_code_blocks_with_context(section_text)

        for block_idx, (heading_text, fence_info, code_text) in enumerate(blocks, start=1):
            # Skip obvious non-table code blocks.
            if "|" not in code_text and "\t" not in code_text and re.search(r"\s{2,}", code_text) is None:
                continue

            obs = table_block_to_observations(
                block_text=code_text,
                section_index=sec_idx,
                section_title=section_title,
                heading_text=heading_text,
                block_index=block_idx,
            )
            all_obs.extend(obs)

    return all_obs


def make_color_map(presets: List[str]) -> Dict[str, object]:
    cmap = plt.get_cmap("tab20")
    return {preset: cmap(i % cmap.N) for i, preset in enumerate(presets)}
    # return {preset: f"C{i%10}" for i, preset in enumerate(presets)}


def sort_datasets_for_metric(records: List[Observation], metric: str) -> List[str]:
    by_dataset: Dict[str, List[float]] = defaultdict(list)
    for r in records:
        if r.metric == metric:
            by_dataset[r.dataset].append(r.value)

    items = []
    higher = metric_is_higher_better(metric)

    for dataset, vals in by_dataset.items():
        if not vals:
            continue
        avg = sum(vals) / len(vals)
        # Worst-to-best left-to-right.
        # higher-is-better => ascending by average
        # lower-is-better  => descending by average
        sort_score = avg if higher else -avg
        items.append((sort_score, dataset))

    items.sort(key=lambda t: (t[0], t[1]))
    return [dataset for _, dataset in items]


def plot_section_to_pdf(section_title: str, records: List[Observation], pdf: PdfPages) -> None:
    if not records:
        return

    metrics = list(dict.fromkeys(r.metric for r in records))
    presets = list(dict.fromkeys(r.preset for r in records))
    preset_colors = make_color_map(presets)

    n_metrics = len(metrics)
    ncols = min(n_metrics, 3)
    nrows = math.ceil(n_metrics / ncols)

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(5 * ncols, 3 * nrows),
        squeeze=False,
    )
    # fig.suptitle(section_title, fontsize=14, y=0.98)

    legend_handles = [
        plt.Line2D([], [], linestyle="", marker="o", color=preset_colors[preset], label=preset.replace('?', '\n').replace('&', '\n'))
        for preset in presets
    ]

    for idx, metric in enumerate(metrics):
        ax = axes[idx // ncols][idx % ncols]
        metric_records = [r for r in records if r.metric == metric]

        dataset_order = sort_datasets_for_metric(metric_records, metric)
        if not dataset_order:
            ax.set_title(metric)
            ax.axis("off")
            continue

        x_map = {dataset: i for i, dataset in enumerate(dataset_order)}

        # Scatter points at the same x location for a dataset; y encodes the metric value.
        # This is the cartesian version of a radar-style comparison across presets.
        for preset in presets:
            sub = [r for r in metric_records if r.preset == preset and r.dataset in x_map]
            if not sub:
                continue

            xs = [x_map[r.dataset] for r in sub]
            ys = [r.value for r in sub]
            # print(preset, xs, ys)

            ax.scatter(
                xs,
                ys,
                marker='x',
                s=32,
                alpha=0.9,
                color=preset_colors[preset],
                label=preset,
                zorder=3,
            )

        is_higher_better = metric_is_higher_better(metric)
        if is_higher_better != None:
            metric = metric + (" ↑" if is_higher_better else " ↓")
        ax.set_title(metric)
        ax.set_xticks(range(len(dataset_order)))
        ax.set_xticklabels(dataset_order, rotation=35, ha="right")
        # ax.set_xlabel("dataset")
        # ax.set_ylabel("metric value")
        ax.grid(True, axis="y", linestyle="--", alpha=0.35)
        ax.margins(x=0.05)

    # Turn off unused axes.
    total_axes = nrows * ncols
    for idx in range(n_metrics, total_axes):
        axes[idx // ncols][idx % ncols].axis("off")

    if presets:
        fig.legend(
            handles=legend_handles,
            loc="lower center",
            ncol=min(len(presets), 6),
            frameon=False,
            bbox_to_anchor=(0.5, 0.01),
        )

    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    pdf.savefig(fig)
    plt.close(fig)


def group_records_by_section(records: List[Observation]) -> List[Tuple[str, List[Observation]]]:
    """
    One page per section, as requested.
    """
    grouped: Dict[int, List[Observation]] = defaultdict(list)
    titles: Dict[int, str] = {}

    for r in records:
        grouped[r.section_index].append(r)
        titles[r.section_index] = r.section_title

    pages: List[Tuple[str, List[Observation]]] = []
    for sec_idx in sorted(grouped):
        pages.append((titles.get(sec_idx, f"Section {sec_idx}"), grouped[sec_idx]))
    return pages


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot benchmark metrics from raw markdown into a PDF.")
    parser.add_argument("url", help="URL to the raw markdown document")
    parser.add_argument("-o", "--output", default="benchmark_plots.pdf", help="Output PDF path")
    args = parser.parse_args()

    md = fetch_markdown(args.url)
    records = collect_observations(md)

    if not records:
        print("No metric tables were found.", file=sys.stderr)
        return 1

    pages = group_records_by_section(records)
    out_path = Path(args.output).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(out_path) as pdf:
        for title, recs in pages:
            plot_section_to_pdf(title, recs, pdf)

    print(f"Saved {len(pages)} page(s) to: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
