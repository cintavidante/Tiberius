from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional, Sequence
import pandas as pd
import html
import json


DEFAULT_SHOW_COLS = [
    "created_at",
    "run_id",
    "spock",
    "spock_name",
    "config_hash",
    "input_hash",
]


# External CDN assets (works when you view via http.server + port forward).
# If you MUST be fully offline, see note below.
DATATABLES_CSS = "https://cdn.datatables.net/1.13.8/css/jquery.dataTables.min.css"
JQUERY_JS = "https://code.jquery.com/jquery-3.7.1.min.js"
DATATABLES_JS = "https://cdn.datatables.net/1.13.8/js/jquery.dataTables.min.js"


def _make_link(rel_href: str, label: str) -> str:
    # rel_href should be relative to the HTML location (stage root)
    rel_href = html.escape(rel_href, quote=True)
    label = html.escape(label)
    return f'<a href="{rel_href}">{label}</a>'


def _safe_relpath(target: Path, base: Path) -> str:
    # produce a relative path from base -> target, using forward slashes for browsers
    return target.relative_to(base).as_posix()


def write_registry_html(
    stage_root: str | Path,
    spock: int,
    *,
    registry_parquet: str | Path = "_registry.parquet",
    out_html: str | Path = "_registry.html",
    title: Optional[str] = None,
    subtitle: Optional[str] = None,
    show_first_cols: Sequence[str] = DEFAULT_SHOW_COLS,
    add_links: bool = True,
    add_reports_links: bool = True,
    max_rows: Optional[int] = None,
) -> Path:
    """
    Reads stage_root/_registry.parquet and writes a readable, scrollable HTML table.

    - Sticky header + horizontal scroll
    - Clickable links to meta.json / qa.parquet / raw_index.parquet / product file / open folder
    - Optional report links (if reports/ exists inside product_dir)
    - Uses DataTables for search/sort/filter in the browser

    Assumes 'product_dir' column exists and points to a folder under stage_root (recommended).
    """
    stage_root = Path(stage_root)
    reg_path = stage_root / registry_parquet
    out_path = stage_root / out_html

    if not reg_path.exists():
        raise FileNotFoundError(f"Registry not found: {reg_path}")

    df = pd.read_parquet(reg_path)

    if df.empty:
        # still write a minimal HTML
        html_page = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>{html.escape(title or stage_root.name)}</title></head>
<body><h2>{html.escape(title or stage_root.name)} — Registry</h2><p>No runs found.</p></body></html>"""
        out_path.write_text(html_page, encoding="utf-8")
        return out_path

    if max_rows is not None:
        df = df.head(max_rows).copy()

    # Ensure stable ordering: most recent first if created_at exists
    if "created_at" in df.columns:
        df = df.sort_values("created_at", ascending=False).reset_index(drop=True)
    else:
        df = df.reset_index(drop=True)

    # Build link columns
    if add_links and "product_dir" in df.columns:
        # Make product_dir relative (recommended) or keep as-is if already relative
        # We'll create link columns, and optionally keep product_dir itself too.
        def mk_links(prod_dir_str: str) -> dict:
            # product_dir may be relative (best) or absolute.
            prod_dir = Path(prod_dir_str)
            # If it's absolute, try to relativize to stage_root; if not possible, keep absolute href.
            try:
                rel_prod = prod_dir.relative_to(stage_root)
                base_for_links = stage_root
                prod_href = rel_prod.as_posix()
                meta_href = f"{prod_href}/_info/meta.json"
                qa_href = f"{prod_href}/_info/qa.parquet"
                qas_href = f"{prod_href}/_info/qa_summary.parquet"
                raw_href = f"{prod_href}/_info/raw_index.parquet"
                # dataset file depends on storage; try common names
                prod_zarr = f"{prod_href}/product.zarr"
                prod_h5 = f"{prod_href}/product.h5"
            except Exception:
                # product_dir not under stage_root; fall back to direct href
                base_for_links = None
                prod_href = prod_dir.as_posix()
                meta_href = f"{prod_href}/_info/meta.json"
                qa_href = f"{prod_href}/_info/qa.parquet"
                qas_href = f"{prod_href}/_info/qa_summary.parquet"
                raw_href = f"{prod_href}/_info/raw_index.parquet"
                prod_zarr = f"{prod_href}/product.zarr"
                prod_h5 = f"{prod_href}/product.h5"

            links = {
                "open": _make_link(prod_href, "open"),
                "meta": _make_link(meta_href, "meta"),
            }

            # Only show links that are likely to exist; in HTML we can still link even if missing,
            # but it’s nicer to hide missing.
            # We can check existence only if prod_dir is under stage_root.
            if base_for_links is not None:
                abs_prod = stage_root / prod_href
                if (abs_prod / "_info/qa.parquet").exists():
                    links["qa"] = _make_link(f"{prod_href}/_info/qa_preview.html", "qa")
                if (abs_prod / "_info/qa_summary.parquet").exists():
                    links["qa_summary"] = _make_link(qas_href, "qa_summary")
                if (abs_prod / "_info/raw_index.parquet").exists():
                    links["raw_index"] = _make_link(f"{prod_href}/_info/raw_index_preview.html", "raw_index")
                if (abs_prod / "product.zarr").exists():
                    links["product"] = _make_link(prod_zarr, "product.zarr")
                elif (abs_prod / "product.h5").exists():
                    links["product"] = _make_link(prod_h5, "product.h5")
            else:
                # cannot check remotely, include basic ones
                links["qa"] = _make_link(qa_href, "qa")
                links["raw_index"] = _make_link(raw_href, "raw_index")
                links["product"] = _make_link(prod_zarr, "product")

            # Reports (PNG/HTML) links
            if add_reports_links and base_for_links is not None:
                abs_prod = stage_root / prod_href
                reports_dir = abs_prod / "reports"
                if reports_dir.exists() and reports_dir.is_dir():
                    # link a few common files if present
                    candidates = [
                        "master_bias.png",
                        "master_flat_stack.png",
                        "master_flat_model.png",
                        "ratio_to_stack.png",
                        "profile_1d.png",
                        "report.html",
                    ]
                    for name in candidates:
                        p = reports_dir / name
                        if p.exists():
                            rel = _safe_relpath(p, stage_root)
                            links[f"rep:{name}"] = _make_link(rel, name)

            return links

        link_rows = df["product_dir"].astype(str).apply(mk_links)
        link_df = pd.DataFrame(list(link_rows))
        # put links first
        df = pd.concat([link_df, df], axis=1)

    # Choose column order: show_first_cols + rest
    def order_cols(cols: Iterable[str], spock) -> list[str]:
        cols = list(cols)
        ordered = []

        def add(c):
            if c in cols and c not in ordered:
                ordered.append(c)
        
        add("run_id")
        if spock != 0:
            add("created_at")

        add("spock")
        add("spock_name")

        add("raw_index")
        add("qa")
        add("qa_summary")

        for c in sorted([c for c in cols if c.startswith("spock") and "__" in c]):
            ordered.append(c)

        add("config_hash")
        add("input_hash")

        if spock != 0:
            add("open")
    
        add("meta")

        for c in cols:
            if c not in ordered and c not in ["open", "created_at", "config_hash", "run_id", "input_hash", "product_dir", "product", "context__instrument_id", "context__project_name", "context__storage_format", "context__project_date"]:
                ordered.append(c)
        
        return ordered

    df = df[order_cols(df.columns, spock)]

    stage_title = title or f"{stage_root.name} — Registry"
    n_runs = len(df)

    subtitle = subtitle

    # Convert dataframe to HTML (escape=False so our <a> links work)
    table_html = df.to_html(index=False, escape=False, table_id="registryTable")

    # Full page
    html_page = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>{html.escape(stage_title)}</title>
  <link rel="stylesheet" href="{DATATABLES_CSS}">
  <script src="{JQUERY_JS}"></script>
  <script src="{DATATABLES_JS}"></script>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 16px; }}
    .subtitle {{font-size: 18px; color: #666; margin-bottom:15px; }}
    .meta {{ color: #444; margin-bottom: 10px; }}
    .table-wrap {{ overflow-x: auto; }}
    table.dataTable thead th {{
      background: #f2f2f2;
      white-space: nowrap;
    }}
    td {{ white-space: nowrap; }}
    a {{ text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    code {{ background: #f6f6f6; padding: 2px 4px; border-radius: 4px; }}
  </style>
</head>
<body>
  <h2>{html.escape(stage_title)}</h2>
  <div class="subtitle">{html.escape(subtitle)}</div>
  <div class="meta">
    <div>Runs: <b>{n_runs}</b></div>
    <div>Registry: <code>{html.escape(reg_path.name)}</code></div>
    <div>Tip: use the search box + click column headers to sort.</div>
  </div>

  <div class="table-wrap">
    {table_html}
  </div>

  <script>
    $(document).ready(function() {{
      var table = $('#registryTable').DataTable({{
        pageLength: 25,
        lengthMenu: [10, 25, 50, 100, 250],
        scrollX: true,
        autoWidth: false
      }});

      table.columns.adjust();

      $(window).on('resize', function() {{
          table.columns.adjust();
      }});
    }});
  </script>
</body>
</html>
"""
    out_path.write_text(html_page, encoding="utf-8")
    return out_path



def _write_table_preview(df, path: Path, title: str):
    html_table = df.to_html(index=False, escape=True)
    page = f"""<!doctype html>
<html><head>
<meta charset="utf-8">
<title>{title}</title>
<style>
body {{ font-family: Arial; margin: 16px; }}
.table-wrap {{ overflow-x: auto; }}
td, th {{ white-space: nowrap; padding: 6px; border: 1px solid #ddd; }}
th {{ position: sticky; top: 0; background: #f2f2f2; }}
table {{ border-collapse: collapse; font-size: 14px; }}
</style>
</head><body>
<h2>{title}</h2>
<div class="table-wrap">
{html_table}
</div>
</body></html>
"""
    path.write_text(page, encoding="utf-8")