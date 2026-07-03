from __future__ import annotations

import csv
import html
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple, Union


MOLSTAR_VERSION = "5.4.2"
MOLSTAR_JS_URL = f"https://cdn.jsdelivr.net/npm/molstar@{MOLSTAR_VERSION}/build/viewer/molstar.js"
MOLSTAR_CSS_URL = f"https://cdn.jsdelivr.net/npm/molstar@{MOLSTAR_VERSION}/build/viewer/molstar.css"
MANIFEST_SCHEMA_VERSION = "flipbook-molstar-viewer/v1"
MASK_OPACITY = 0.30
DEFAULT_RADIUS_MIN = 0.63
DEFAULT_RADIUS_MAX = 3.18
DEFAULT_THICKNESS_SCALE = 1.0
DEFAULT_TILE_PADDING_FACTOR = 1.55
DEFAULT_SPACING_STEP = 0.025
DEFAULT_HTML_NAME = "rmsx_molstar_flipbook.html"
DEFAULT_MANIFEST_NAME = "rmsx_molstar_manifest.json"
SLICE_RE = re.compile(r"^slice_(\d+)_first_frame\.pdb$")
ASSET_DIR = Path(__file__).resolve().parent / "molstar_static"
CAMERA_MODES = {"orthographic", "perspective"}

DEFAULT_COLOR_PALETTES = {
    "magma": [
        "#000004", "#120D32", "#331068", "#5A167E", "#7D2482",
        "#A3307E", "#C83E73", "#E95562", "#F97C5D", "#FEA873",
        "#FED395", "#FCFDBF",
    ],
    "inferno": [
        "#000004", "#140B35", "#3A0963", "#60136E", "#85216B",
        "#A92E5E", "#CB4149", "#E65D2F", "#F78311", "#FCAD12",
        "#F5DB4B", "#FCFFA4",
    ],
    "plasma": [
        "#0D0887", "#3E049C", "#6300A7", "#8707A6", "#A62098",
        "#C03A83", "#D5546E", "#E76F5A", "#F58C46", "#FDAD32",
        "#FCD225", "#F0F921",
    ],
    "viridis": [
        "#440154", "#482173", "#433E85", "#38598C", "#2D708E",
        "#25858E", "#1E9B8A", "#2BB07F", "#51C56A", "#85D54A",
        "#C2DF23", "#FDE725",
    ],
    "cividis": [
        "#00204D", "#00306F", "#2A406C", "#48526B", "#5E626E",
        "#727374", "#878479", "#9E9677", "#B6A971", "#D0BE67",
        "#EAD357", "#FFEA46",
    ],
    "rocket": [
        "#03051A", "#221331", "#451C47", "#6A1F56", "#921C5B",
        "#B91657", "#D92847", "#ED513E", "#F47C56", "#F6A47B",
        "#F7C9AA", "#FAEBDD",
    ],
    "mako": [
        "#0B0405", "#231526", "#35264C", "#403A75", "#3D526D",
        "#366DA0", "#3487A6", "#35A1AB", "#43BBAD", "#6CD3AD",
        "#ADE3C0", "#DEF5E5",
    ],
    "turbo": [
        "#30123B", "#4454C4", "#4490FE", "#1FC8DE", "#29EFA2",
        "#7DFF56", "#C1F334", "#F1CA3A", "#FE922A", "#EA4F0D",
        "#BE2102", "#7A0403",
    ],
}


@dataclass
class MolstarFlipbookResult:
    """Paths and notebook HTML for a generated Molstar Flipbook."""

    html_path: Path
    manifest_path: Path
    manifest: Dict[str, Any]
    html: str
    iframe_html: str

    def _repr_html_(self) -> str:
        return self.iframe_html

    def display(self) -> "MolstarFlipbookResult":
        try:
            from IPython.display import HTML, display  # type: ignore
        except Exception as exc:  # pragma: no cover - depends on notebook runtime
            raise RuntimeError("IPython is required to display the Molstar Flipbook inline.") from exc
        display(HTML(self.iframe_html))
        return self

    def __str__(self) -> str:
        return str(self.html_path)


def _asset_text(relative_path: str) -> str:
    return (ASSET_DIR / relative_path).read_text(encoding="utf-8")


def _asset_uri(relative_path: str) -> str:
    return (ASSET_DIR / relative_path).resolve().as_uri()


def _slice_sort_key(path: Path) -> int:
    match = SLICE_RE.match(path.name)
    return int(match.group(1)) if match else 0


def _slice_paths(directory: Path) -> List[Path]:
    paths = [path for path in directory.glob("slice_*_first_frame.pdb") if SLICE_RE.match(path.name)]
    return sorted(paths, key=_slice_sort_key)


def _parse_float(text: str) -> Optional[float]:
    try:
        value = float(text.strip())
    except (TypeError, ValueError):
        return None
    if not math.isfinite(value):
        return None
    return value


def _residue_key(chain_id: str, residue_id: str) -> str:
    return f"{chain_id}:{residue_id}" if chain_id else residue_id


def _format_residue_label(residue_id: str, chain_id: str, segid: str) -> str:
    if chain_id and segid and chain_id != segid:
        return f"{residue_id} / chain {chain_id} (segid {segid})"
    if chain_id:
        return f"{residue_id} / chain {chain_id}"
    if segid:
        return f"{residue_id} / segid {segid}"
    return residue_id


def _parse_residue_bfactors(pdb_text: str) -> Tuple[List[Dict[str, Any]], Dict[str, str]]:
    residues: Dict[str, Dict[str, Any]] = {}
    chain_aliases: Dict[str, str] = {}
    order = 0

    for line in pdb_text.splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        padded = line.ljust(80)
        residue_id = (padded[22:26].strip() + padded[26:27].strip()).strip()
        if not residue_id:
            continue
        chain_id = padded[21:22].strip()
        segid = padded[72:76].strip()
        key_chain_id = chain_id or (segid if len(segid) == 1 else "")
        key = _residue_key(key_chain_id, residue_id)
        bfactor = _parse_float(padded[60:66])
        if bfactor is None:
            continue

        if key_chain_id:
            if chain_id:
                chain_aliases.setdefault(chain_id, key_chain_id)
            if segid:
                chain_aliases.setdefault(segid, key_chain_id)

        atom_name = padded[12:16].strip()
        existing = residues.get(key)
        use_record = existing is None or (atom_name == "CA" and not existing.get("_is_ca"))
        if use_record:
            residues[key] = {
                "id": residue_id,
                "chain": key_chain_id,
                "segid": segid,
                "key": key,
                "label": _format_residue_label(residue_id, key_chain_id, segid),
                "bfactor": bfactor,
                "_order": existing.get("_order", order) if existing else order,
                "_is_ca": atom_name == "CA",
            }
        if existing is None:
            order += 1

    sorted_residues = sorted(residues.values(), key=lambda item: item["_order"])
    for residue in sorted_residues:
        residue.pop("_order", None)
        residue.pop("_is_ca", None)
    return sorted_residues, chain_aliases


def _parse_mask_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def read_mask_summary(mask_path: Path, chain_aliases: Optional[Mapping[str, str]] = None) -> Dict[str, Any]:
    if not mask_path.is_file():
        return {"maskedResidues": 0, "totalResidues": 0, "masked": [], "maskedKeys": []}

    aliases = dict(chain_aliases or {})
    masked: List[Dict[str, str]] = []
    masked_keys: List[str] = []
    seen_keys = set()
    total = 0

    with mask_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            total += 1
            if not _parse_mask_bool(row.get("Masked", "")):
                continue
            residue_id = (row.get("ResidueID") or "").strip()
            if not residue_id:
                continue
            original_chain = (row.get("ChainID") or "").strip()
            mapped_chain = aliases.get(original_chain, original_chain)
            key = _residue_key(mapped_chain, residue_id)
            label = (
                f"{residue_id} / chain {mapped_chain}"
                if mapped_chain and mapped_chain == original_chain
                else f"{residue_id} / chain {mapped_chain} (source {original_chain})"
                if mapped_chain
                else residue_id
            )
            masked.append(
                {
                    "id": residue_id,
                    "chain": mapped_chain,
                    "sourceChain": original_chain,
                    "key": key,
                    "label": label,
                }
            )
            if key not in seen_keys:
                masked_keys.append(key)
                seen_keys.add(key)

    return {
        "maskedResidues": len(masked),
        "totalResidues": total,
        "masked": masked,
        "maskedKeys": masked_keys,
    }


def _read_slices_and_residues(directory: Path) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], Dict[str, float], Dict[str, Any], Dict[str, str]]:
    paths = _slice_paths(directory)
    if not paths:
        raise FileNotFoundError(
            f"No files found matching pattern 'slice_<number>_first_frame.pdb' in directory '{directory}'."
        )

    slices: List[Dict[str, Any]] = []
    residues_by_key: Dict[str, Dict[str, Any]] = {}
    summaries: Dict[str, Any] = {}
    all_values: List[float] = []
    chain_aliases: Dict[str, str] = {}

    for path in paths:
        slice_index = _slice_sort_key(path)
        rmsx_column = f"slice_{slice_index}.dcd"
        pdb_text = path.read_text(encoding="utf-8")
        parsed_residues, aliases = _parse_residue_bfactors(pdb_text)
        chain_aliases.update(aliases)
        values: List[float] = []
        max_value: Optional[float] = None
        max_residue = ""

        for parsed in parsed_residues:
            key = parsed["key"]
            value = float(parsed["bfactor"])
            values.append(value)
            all_values.append(value)
            if max_value is None or value > max_value:
                max_value = value
                max_residue = str(parsed["id"])
            residue = residues_by_key.setdefault(
                key,
                {
                    "id": str(parsed["id"]),
                    "chain": str(parsed["chain"]),
                    "segid": str(parsed.get("segid", "")),
                    "key": key,
                    "label": str(parsed["label"]),
                    "values": {},
                },
            )
            residue["values"][rmsx_column] = value

        if values:
            summaries[rmsx_column] = {
                "min": min(values),
                "max": max(values),
                "mean": sum(values) / len(values),
                "maxResidue": max_residue,
                "residueCount": len(values),
            }

        slices.append(
            {
                "index": slice_index,
                "id": f"slice_{slice_index}",
                "label": f"Slice {slice_index}",
                "filename": path.name,
                "rmsxColumn": rmsx_column,
                "pdb": pdb_text,
            }
        )

    if not residues_by_key:
        raise ValueError(f"No ATOM/HETATM residue B-factors found in {directory}.")
    if not all_values:
        raise ValueError(f"No numeric B-factors found in {directory}.")

    residues = list(residues_by_key.values())
    domain = {"min": min(all_values), "max": max(all_values)}
    return slices, residues, domain, summaries, chain_aliases


def _safe_float(value: Any, default: float) -> float:
    parsed = _parse_float(str(value))
    return default if parsed is None else parsed


def _normalize_camera_mode(camera_mode: str) -> str:
    normalized = str(camera_mode).strip().lower()
    if normalized not in CAMERA_MODES:
        raise ValueError("camera_mode must be one of: 'orthographic' or 'perspective'.")
    return normalized


def _color_stops(domain: Mapping[str, float], palette_colors: Iterable[str]) -> List[Dict[str, Any]]:
    colors = list(palette_colors)
    if len(colors) < 2:
        raise ValueError("At least two colors are required to build a Molstar palette.")
    interval = (float(domain["max"]) - float(domain["min"])) / (len(colors) - 1)
    return [
        {"bfactor": round(float(domain["min"]) + (index * interval), 2), "color": color}
        for index, color in enumerate(colors)
    ]


def _build_flipbook_reference(
    slices: List[Dict[str, Any]],
    domain: Mapping[str, float],
    palette_name: str,
    palette_colors: Iterable[str],
    spacing_factor: Any,
) -> Dict[str, Any]:
    default_spacing = _safe_float(spacing_factor, 1.0)
    color_stops = _color_stops(domain, palette_colors)
    color_mapping = ":".join(f"{stop['bfactor']},{stop['color']}" for stop in color_stops)
    num_models = len(slices)
    return {
        "viewer": "molstar",
        "defaultColumns": num_models,
        "minimumSpacingFactor": 0.0,
        "maximumSpacingFactor": 2.5,
        "defaultSpacingFactor": default_spacing,
        "tilePaddingFactor": DEFAULT_TILE_PADDING_FACTOR,
        "palette": palette_name,
        "colorStops": color_stops,
        "commands": [
            "view",
            "define axis",
            "turn x 90",
            "color byattribute bfactor",
            "worm bfactor",
            "lighting soft",
            "graphics silhouettes true",
            "set bgColor white",
            f"color byattribute a:bfactor #1-{num_models} target absc palette {color_mapping}",
            f"tile all columns {num_models} spacingFactor {default_spacing:g}",
        ],
    }


def build_molstar_manifest(
    directory: Union[str, Path],
    palette: str = "viridis",
    min_bfactor: Optional[float] = None,
    max_bfactor: Optional[float] = None,
    spacing_factor: Any = 1,
    camera_mode: str = "orthographic",
    title: str = "RMSX Molstar Flipbook",
    mask_filename: str = "masked_residues.csv",
    color_palettes: Optional[Mapping[str, List[str]]] = None,
) -> Dict[str, Any]:
    """Build a Molstar Flipbook manifest from RMSX slice PDB outputs."""
    directory_path = Path(directory)
    palettes = {name: list(colors) for name, colors in (color_palettes or DEFAULT_COLOR_PALETTES).items()}
    if palette not in palettes:
        raise ValueError(f"Palette '{palette}' is not defined.")
    normalized_camera_mode = _normalize_camera_mode(camera_mode)

    slices, residues, observed_domain, summaries, chain_aliases = _read_slices_and_residues(directory_path)
    domain_min = observed_domain["min"] if min_bfactor is None else min(observed_domain["min"], float(min_bfactor))
    domain_max = observed_domain["max"] if max_bfactor is None else max(observed_domain["max"], float(max_bfactor))
    domain = {"min": domain_min, "max": domain_max}
    default_color_min = observed_domain["min"] if min_bfactor is None else float(min_bfactor)
    default_color_max = observed_domain["max"] if max_bfactor is None else float(max_bfactor)
    domain_span = max(0.000001, domain["max"] - domain["min"])
    color_domain_step = round(max(0.1, domain_span / 50), 3)

    mask_summary = read_mask_summary(directory_path / mask_filename, chain_aliases=chain_aliases)
    if not mask_summary.get("totalResidues"):
        mask_summary["totalResidues"] = len(residues)

    palette_colors = palettes[palette]
    pdb_byte_sizes = [len(slice_entry["pdb"].encode("utf-8")) for slice_entry in slices]
    manifest: Dict[str, Any] = {
        "schemaVersion": MANIFEST_SCHEMA_VERSION,
        "title": title,
        "citation": {
            "text": "RMSX/Flipbook paper, Scientific Reports (2026)",
            "doi": "10.1038/s41598-026-39869-7",
            "url": "https://doi.org/10.1038/s41598-026-39869-7",
        },
        "slices": slices,
        "summaries": summaries,
        "domain": domain,
        "maskSummary": mask_summary,
        "residues": residues,
        "palette": {"name": palette, "colors": palette_colors},
        "availablePalettes": {name: palettes[name] for name in sorted(palettes)},
        "maskOpacity": MASK_OPACITY,
        "presentation": {
            "defaultLayout": "tiled",
            "availableLayouts": ["tiled"],
            "layoutUrlParam": "layout",
        },
        "playback": {
            "defaultDelayMs": 900,
            "minDelayMs": 150,
            "maxDelayMs": 5000,
            "delayStepMs": 50,
            "delayUrlParam": "delayMs",
        },
        "rotationModel": {
            "mode": "per-slice local coordinate transform",
            "pivot": "geometric center of each full slice before mask splitting",
            "layoutOrder": "rotate around local center, then place that center on a shared Flipbook slot anchor plus tile offset",
            "defaultRotation": {"x": 90, "y": 0, "z": 0},
            "defaultRotationSource": "ChimeraX Flipbook command: turn x 90",
        },
        "molstarRenderStyle": {
            "preset": "clean-interactive",
            "backgroundColor": "#ffffff",
            "cameraMode": normalized_camera_mode,
            "outline": True,
            "ambientOcclusion": False,
            "illumination": False,
            "multiSample": False,
        },
        "visualMapping": {
            "mode": "read RMSX values from the manifest and normalize them into the PDB B-factor column for Molstar uncertainty rendering",
            "bfactorDomain": [0, 1],
            "defaultRadiusMin": DEFAULT_RADIUS_MIN,
            "defaultRadiusMax": DEFAULT_RADIUS_MAX,
            "defaultThicknessScale": DEFAULT_THICKNESS_SCALE,
            "radiusStep": 0.05,
            "defaultColorMin": default_color_min,
            "defaultColorMax": default_color_max,
            "colorDomainStep": color_domain_step,
            "colorTheme": "uncertainty",
            "molstarUncertaintyReversesColorList": True,
            "paletteOrder": "Flipbook low-to-high palette is reversed before passing to Molstar because Molstar uncertainty coloring reverses its color list internally",
            "radiusUrlParams": ["radiusMin", "radiusMax"],
            "thicknessUrlParam": "thickness",
            "colorDomainUrlParams": ["colorMin", "colorMax"],
            "paletteUrlParam": "palette",
        },
        "selectedResidueMarker": {
            "enabledDefault": False,
            "type": "spacefill",
            "color": "#111827",
            "alpha": 0.86,
            "sizeFactor": 0.36,
            "toggleUrlParam": "marker",
            "residueUrlParam": "residue",
        },
        "keyboardShortcuts": {
            "enabled": True,
            "source": "VMD Flipbook hotkey parity subset",
            "rotationStepDegrees": 5,
            "spacingStep": DEFAULT_SPACING_STEP,
            "thicknessStep": 0.05,
        },
        "reportPayload": {
            "embeddedPdbCount": len(slices),
            "embeddedPdbBytes": sum(pdb_byte_sizes),
            "largestEmbeddedPdbBytes": max(pdb_byte_sizes) if pdb_byte_sizes else 0,
            "residueCount": len(residues),
            "largeReportWarningBytes": 10_000_000,
            "strategy": "notebook-friendly HTML with embedded PDB slices",
        },
        "flipbookReference": _build_flipbook_reference(
            slices,
            domain,
            palette,
            palette_colors,
            spacing_factor,
        ),
        "molstar": {
            "version": MOLSTAR_VERSION,
            "jsUrl": MOLSTAR_JS_URL,
            "cssUrl": MOLSTAR_CSS_URL,
            "localJsUrl": _asset_uri("vendor/molstar/5.4.2/molstar.js"),
            "localCssUrl": _asset_uri("vendor/molstar/5.4.2/molstar.css"),
        },
    }
    manifest["reportPayload"]["estimatedJsonBytes"] = len(json.dumps(manifest).encode("utf-8"))
    return manifest


def _script_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=True).replace("</", "<\\/")


def build_molstar_html(
    manifest: Mapping[str, Any],
    asset_mode: str = "cdn",
    title: Optional[str] = None,
) -> str:
    """Build standalone HTML that renders a Molstar Flipbook manifest."""
    normalized_mode = str(asset_mode).strip().lower()
    if normalized_mode not in {"cdn", "local", "inline"}:
        raise ValueError("asset_mode must be one of: 'cdn', 'local', or 'inline'.")

    escaped_title = html.escape(title or str(manifest.get("title") or "RMSX Molstar Flipbook"))
    viewer_script = _asset_text("script.js")
    manifest_json = _script_json({"manifest": manifest})

    if normalized_mode == "inline":
        molstar_css = _asset_text("vendor/molstar/5.4.2/molstar.css")
        molstar_js = _asset_text("vendor/molstar/5.4.2/molstar.js")
        asset_tags = f"<style>\n{molstar_css}\n</style>\n<script>\n{molstar_js}\n</script>"
    elif normalized_mode == "local":
        css_url = _asset_uri("vendor/molstar/5.4.2/molstar.css")
        js_url = _asset_uri("vendor/molstar/5.4.2/molstar.js")
        asset_tags = (
            f'<link rel="stylesheet" type="text/css" href="{html.escape(css_url, quote=True)}">\n'
            f'<script src="{html.escape(js_url, quote=True)}"></script>'
        )
    else:
        asset_tags = (
            f'<link rel="stylesheet" type="text/css" href="{MOLSTAR_CSS_URL}">\n'
            f'<script src="{MOLSTAR_JS_URL}"></script>'
        )

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escaped_title}</title>
  {asset_tags}
</head>
<body>
  <div id="app"></div>
  <script>
    document.getElementById("app").dataset.incoming = JSON.stringify({manifest_json});
  </script>
  <script>
{viewer_script}
  </script>
</body>
</html>
"""


def build_notebook_iframe_html(html_document: str, height: int = 720, width: str = "100%") -> str:
    height_value = max(320, int(height))
    escaped_srcdoc = html.escape(html_document, quote=True)
    escaped_width = html.escape(str(width), quote=True)
    return (
        f'<iframe title="RMSX Molstar Flipbook" srcdoc="{escaped_srcdoc}" '
        f'width="{escaped_width}" height="{height_value}" '
        'style="border:1px solid #d7dce2; border-radius:6px; width:100%;" '
        'sandbox="allow-scripts allow-same-origin allow-downloads" '
        'referrerpolicy="no-referrer"></iframe>'
    )


def write_molstar_flipbook(
    directory: Union[str, Path],
    palette: str = "viridis",
    min_bfactor: Optional[float] = None,
    max_bfactor: Optional[float] = None,
    spacing_factor: Any = 1,
    camera_mode: str = "orthographic",
    output_html: Optional[Union[str, Path]] = None,
    output_manifest: Optional[Union[str, Path]] = None,
    asset_mode: str = "cdn",
    iframe_height: int = 720,
    title: str = "RMSX Molstar Flipbook",
    mask_filename: str = "masked_residues.csv",
    color_palettes: Optional[Mapping[str, List[str]]] = None,
) -> MolstarFlipbookResult:
    """Write a Molstar Flipbook HTML file and JSON manifest next to RMSX PDB slices."""
    directory_path = Path(directory)
    manifest = build_molstar_manifest(
        directory_path,
        palette=palette,
        min_bfactor=min_bfactor,
        max_bfactor=max_bfactor,
        spacing_factor=spacing_factor,
        camera_mode=camera_mode,
        title=title,
        mask_filename=mask_filename,
        color_palettes=color_palettes,
    )
    html_document = build_molstar_html(manifest, asset_mode=asset_mode, title=title)

    html_path = Path(output_html) if output_html is not None else directory_path / DEFAULT_HTML_NAME
    manifest_path = (
        Path(output_manifest)
        if output_manifest is not None
        else directory_path / DEFAULT_MANIFEST_NAME
    )
    html_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    html_path.write_text(html_document, encoding="utf-8")
    iframe_html = build_notebook_iframe_html(html_document, height=iframe_height)
    return MolstarFlipbookResult(
        html_path=html_path,
        manifest_path=manifest_path,
        manifest=manifest,
        html=html_document,
        iframe_html=iframe_html,
    )
