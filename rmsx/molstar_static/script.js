(function () {
  "use strict";

  const SCHEMA_VERSION = "flipbook-molstar-viewer/v1";
  const STATIC_BASE = "/static/plugins/visualizations/flipbook_molstar/static/";
  const scriptUrl = document.currentScript?.src || "";
  function resolveStaticAsset(path) {
    const candidates = [
      scriptUrl,
      `${window.location.origin}${STATIC_BASE}`
    ];
    for (const candidate of candidates) {
      try {
        return new URL(path, candidate).href;
      } catch (_error) {
        // Galaxy may inject visualization scripts with a non-URL placeholder
        // during local Planemo runs; fall through to the fixed static route.
      }
    }
    return `${STATIC_BASE}${path}`;
  }
  const LOCAL_MOLSTAR_JS = resolveStaticAsset("vendor/molstar/5.4.2/molstar.js");
  const LOCAL_MOLSTAR_CSS = resolveStaticAsset("vendor/molstar/5.4.2/molstar.css");
  const CDN_MOLSTAR_JS = "https://cdn.jsdelivr.net/npm/molstar@5.4.2/build/viewer/molstar.js";
  const CDN_MOLSTAR_CSS = "https://cdn.jsdelivr.net/npm/molstar@5.4.2/build/viewer/molstar.css";
  const VISUAL_MIN = 0;
  const VISUAL_MAX = 1;
  const DEFAULT_MAX_SPACING_FACTOR = 1.5;
  const DEFAULT_SPACING_STEP = 0.025;
  const incoming = JSON.parse(document.getElementById("app")?.dataset?.incoming || "{}");
  const visualizationConfig = incoming.visualization_config || {};
  const MANAGED_URL_PARAMS = [
    "layout", "slice", "slices", "sliceA", "sliceB", "sliceC",
    "thickness", "radiusMin", "radiusMax", "colorMin", "colorMax", "palette",
    "spacing", "columns", "rotX", "rotY", "rotZ", "residue", "marker", "delayMs",
    "localDrag", "rotateSensitivity", "render", "outline", "cameraMode"
  ];
  const URL_PARAMS = initialUrlParams();
  const LAYOUTS = new Set(["tiled"]);
  const CONTROL_PANEL_KEYS = ["view", "style", "rotation", "metrics"];
  const RENDER_PRESETS = new Set(["clean-interactive", "soft"]);
  const CAMERA_MODES = new Set(["orthographic", "perspective"]);

  let REPORT = null;
  let viewer = null;
  let dragState = null;
  let renderToken = 0;
  let queuedSceneUpdate = null;
  let interactiveFrame = null;
  let resizeObserver = null;
  let resizeResetTimer = null;

  const state = {
    layout: "tiled",
    currentIndex: 0,
    visible: new Set(),
    paletteName: "viridis",
    thickness: 1,
    spacing: 1,
    columns: 1,
    rotation: { x: 90, y: 0, z: 0 },
    rotationMatrix: null,
    colorMin: 0,
    colorMax: 1,
    radiusMin: 0.63,
    radiusMax: 3.18,
    marker: false,
    localDrag: true,
    rotationSensitivity: 0.35,
    renderMode: "clean-interactive",
    outline: true,
    cameraMode: "orthographic",
    activePanel: "view",
    selectedResidueKey: "",
    representationMode: "-",
    records: [],
    loaded: false,
    liveTransforms: false,
    forceCoordinateFallback: false,
    molstarAssetSource: "-",
    manifestSource: "-",
    keyboardShortcutCount: 0,
    lastKeyboardAction: "-",
    lastRotationGesture: null
  };

  document.body.innerHTML = `
    <main class="rmsx-app">
      <aside class="rmsx-controls" data-testid="molstar-controls-sidebar">
        <div class="controls-heading">
          <h1>RMSX Flipbook</h1>
          <div class="control-row primary-row">
            <button id="resetViewButton" type="button" data-testid="molstar-reset">Reset View</button>
          </div>
        </div>
        <div id="status" class="status sidebar-status">Loading RMSX manifest...</div>
        <p class="citation-note">Please cite: RMSX/Flipbook paper, Scientific Reports (2026), doi:<a href="https://doi.org/10.1038/s41598-026-39869-7" target="_blank" rel="noopener noreferrer">10.1038/s41598-026-39869-7</a>.</p>
        <div class="control-panels" data-testid="molstar-control-panels">
          <details class="control-panel active" open data-panel="view" data-testid="molstar-panel-layout">
            <summary>View</summary>
            <div class="panel-grid">
              <label>Spacing <input id="spacingRange" type="range" min="0" max="1.5" value="1" step="0.025" data-testid="molstar-spacing-range"><input id="spacingNumber" type="number" min="0" max="1.5" value="1" step="0.025" data-testid="molstar-spacing-number"></label>
              <label>Cols <input id="columnsNumber" type="number" min="1" value="1" step="1" data-testid="molstar-columns-number"></label>
              <div class="slice-visibility">
                <div class="field-label">Slices</div>
                <div id="sliceChips" class="chips sidebar-chips" data-testid="molstar-slice-chips"></div>
              </div>
            </div>
          </details>
          <details class="control-panel" data-panel="style" data-testid="molstar-panel-scale">
            <summary>Style</summary>
            <div class="panel-grid">
              <label>Palette <select id="paletteSelect" data-testid="molstar-palette-select"></select></label>
              <div class="legend" data-testid="molstar-rmsx-legend">
                <div id="legendColorBar" class="bar" aria-hidden="true"></div>
                <div class="legend-values" aria-label="RMSX color domain">
                  <span class="legend-stop"><span id="legendLowSwatch" class="legend-swatch" aria-hidden="true"></span><span id="domainMin">-</span></span>
                  <span class="legend-stop"><span id="legendMidSwatch" class="legend-swatch" aria-hidden="true"></span><span id="domainMid">-</span></span>
                  <span class="legend-stop"><span id="legendHighSwatch" class="legend-swatch" aria-hidden="true"></span><span id="domainMax">-</span></span>
                </div>
                <div class="radius-legend" aria-label="RMSX radius domain" data-testid="molstar-radius-legend">
                  <span class="radius-stop"><span id="legendLowRadius" class="radius-dot" aria-hidden="true"></span><span id="legendLowRadiusLabel">-</span></span>
                  <span class="radius-stop"><span id="legendMidRadius" class="radius-dot" aria-hidden="true"></span><span id="legendMidRadiusLabel">-</span></span>
                  <span class="radius-stop"><span id="legendHighRadius" class="radius-dot" aria-hidden="true"></span><span id="legendHighRadiusLabel">-</span></span>
                </div>
              </div>
              <label>Thickness <input id="thicknessRange" type="range" min="0.25" max="2.5" value="1" step="0.05" data-testid="molstar-thickness-range"><input id="thicknessNumber" type="number" min="0.25" max="2.5" value="1" step="0.05" data-testid="molstar-thickness-number"></label>
              <label>Color min <input id="colorMinNumber" type="number" value="0" step="0.1" data-testid="molstar-color-min-number"></label>
              <label>Color max <input id="colorMaxNumber" type="number" value="1" step="0.1" data-testid="molstar-color-max-number"></label>
              <label>Radius min <input id="radiusMinNumber" type="number" min="0.05" max="5" value="0.63" step="0.05" data-testid="molstar-radius-min-number"></label>
              <label>Radius max <input id="radiusMaxNumber" type="number" min="0.1" max="8" value="3.18" step="0.05" data-testid="molstar-radius-max-number"></label>
              <label class="check-row">Outline <input id="outlineCheckbox" type="checkbox" checked data-testid="molstar-outline-checkbox"></label>
              <button id="resetScaleButton" type="button" data-testid="molstar-reset-scale">Reset Scale</button>
            </div>
          </details>
          <details class="control-panel" data-panel="rotation" data-testid="molstar-panel-rotation">
            <summary>Rotation</summary>
            <div class="panel-grid">
              <label>Rot X <input id="rotationXRange" type="range" min="-180" max="180" value="90" step="1" data-testid="molstar-rotation-x-range"><input id="rotationXNumber" type="number" min="-180" max="180" value="90" step="1" data-testid="molstar-rotation-x-number"></label>
              <label>Rot Y <input id="rotationYRange" type="range" min="-180" max="180" value="0" step="1" data-testid="molstar-rotation-y-range"><input id="rotationYNumber" type="number" min="-180" max="180" value="0" step="1" data-testid="molstar-rotation-y-number"></label>
              <label>Rot Z <input id="rotationZRange" type="range" min="-180" max="180" value="0" step="1" data-testid="molstar-rotation-z-range"><input id="rotationZNumber" type="number" min="-180" max="180" value="0" step="1" data-testid="molstar-rotation-z-number"></label>
              <label>Drag speed <input id="rotateSensitivityRange" type="range" min="0.1" max="3" value="0.35" step="0.05" data-testid="molstar-rotate-sensitivity-range"><input id="rotateSensitivityNumber" type="number" min="0.1" max="3" value="0.35" step="0.05" data-testid="molstar-rotate-sensitivity-number"></label>
              <div class="button-group">
                <button id="rotateXButton" type="button" data-testid="molstar-rotate-x">X +15</button>
                <button id="rotateYButton" type="button" data-testid="molstar-rotate-y">Y +15</button>
                <button id="rotateZButton" type="button" data-testid="molstar-rotate-z">Z +15</button>
                <button id="resetRotationButton" type="button" data-testid="molstar-reset-rotation">Reset</button>
              </div>
            </div>
          </details>
          <details class="control-panel" data-panel="metrics" data-testid="molstar-panel-diagnostics">
            <summary>Metrics</summary>
            <dl class="metrics">
              <div><dt>Slices</dt><dd id="currentSliceMetric">-</dd></div>
              <div><dt>RMSX range</dt><dd id="peakMetric">-</dd></div>
              <div><dt>Mean RMSX</dt><dd id="meanMetric">-</dd></div>
              <div><dt>Peak residue</dt><dd id="peakResidueMetric">-</dd></div>
              <div><dt>Residues</dt><dd id="residueCountMetric">-</dd></div>
              <div><dt>Masked</dt><dd id="maskedMetric">-</dd></div>
            </dl>
          </details>
        </div>
      </aside>
      <section class="rmsx-viewer" data-testid="molstar-report">
        <div id="molstarViewport" class="viewport" data-testid="molstar-viewport"></div>
        <pre id="diagnostics" data-testid="molstar-diagnostics" hidden>{}</pre>
      </section>
    </main>
  `;

  const elements = {
    status: document.getElementById("status"),
    resetViewButton: document.getElementById("resetViewButton"),
    controlPanels: [...document.querySelectorAll("[data-panel]")],
    outlineCheckbox: document.getElementById("outlineCheckbox"),
    paletteSelect: document.getElementById("paletteSelect"),
    thicknessRange: document.getElementById("thicknessRange"),
    thicknessNumber: document.getElementById("thicknessNumber"),
    spacingRange: document.getElementById("spacingRange"),
    spacingNumber: document.getElementById("spacingNumber"),
    columnsNumber: document.getElementById("columnsNumber"),
    rotationXRange: document.getElementById("rotationXRange"),
    rotationXNumber: document.getElementById("rotationXNumber"),
    rotationYRange: document.getElementById("rotationYRange"),
    rotationYNumber: document.getElementById("rotationYNumber"),
    rotationZRange: document.getElementById("rotationZRange"),
    rotationZNumber: document.getElementById("rotationZNumber"),
    rotateSensitivityRange: document.getElementById("rotateSensitivityRange"),
    rotateSensitivityNumber: document.getElementById("rotateSensitivityNumber"),
    rotateXButton: document.getElementById("rotateXButton"),
    rotateYButton: document.getElementById("rotateYButton"),
    rotateZButton: document.getElementById("rotateZButton"),
    resetRotationButton: document.getElementById("resetRotationButton"),
    colorMinNumber: document.getElementById("colorMinNumber"),
    colorMaxNumber: document.getElementById("colorMaxNumber"),
    radiusMinNumber: document.getElementById("radiusMinNumber"),
    radiusMaxNumber: document.getElementById("radiusMaxNumber"),
    resetScaleButton: document.getElementById("resetScaleButton"),
    sliceChips: document.getElementById("sliceChips"),
    viewport: document.getElementById("molstarViewport"),
    diagnostics: document.getElementById("diagnostics"),
    legendColorBar: document.getElementById("legendColorBar"),
    domainMin: document.getElementById("domainMin"),
    domainMid: document.getElementById("domainMid"),
    domainMax: document.getElementById("domainMax"),
    legendLowSwatch: document.getElementById("legendLowSwatch"),
    legendMidSwatch: document.getElementById("legendMidSwatch"),
    legendHighSwatch: document.getElementById("legendHighSwatch"),
    legendLowRadius: document.getElementById("legendLowRadius"),
    legendMidRadius: document.getElementById("legendMidRadius"),
    legendHighRadius: document.getElementById("legendHighRadius"),
    legendLowRadiusLabel: document.getElementById("legendLowRadiusLabel"),
    legendMidRadiusLabel: document.getElementById("legendMidRadiusLabel"),
    legendHighRadiusLabel: document.getElementById("legendHighRadiusLabel"),
    currentSliceMetric: document.getElementById("currentSliceMetric"),
    meanMetric: document.getElementById("meanMetric"),
    peakMetric: document.getElementById("peakMetric"),
    peakResidueMetric: document.getElementById("peakResidueMetric"),
    residueCountMetric: document.getElementById("residueCountMetric"),
    maskedMetric: document.getElementById("maskedMetric"),
    selectedResidueMetric: null,
    selectedRadiusMetric: null,
    selectedColorSwatch: null,
    selectedColorMetric: null,
    styleMetric: null,
    assetMetric: null
  };
  addStyles();

  function addStyles() {
    const style = document.createElement("style");
    style.textContent = `
      * { box-sizing: border-box; }
      body { margin: 0; background: #f7f8fa; color: #1d2630; font: 14px -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; letter-spacing: 0; }
      .rmsx-app { display: grid; grid-template-columns: minmax(236px, 286px) minmax(0, 1fr); grid-template-rows: minmax(0, 1fr); height: 100vh; min-height: 560px; overflow: hidden; }
      .rmsx-controls { position: relative; z-index: 5; min-height: 0; background: #fff; border-right: 1px solid #d7dce2; padding: 9px; overflow-y: auto; }
      .controls-heading { display: grid; grid-template-columns: 1fr; align-items: start; gap: 7px; }
      h1 { margin: 0; font-size: 15px; font-weight: 650; }
      .control-panels { display: grid; grid-template-columns: 1fr; gap: 6px; margin-top: 7px; }
      .control-panel { border: 1px solid #d7dce2; border-radius: 6px; background: #fff; min-width: 0; }
      .control-panel.active { border-color: #b8c7dc; box-shadow: inset 3px 0 0 #1f6feb; }
      .control-panel summary { min-height: 29px; padding: 5px 8px; cursor: pointer; font-weight: 650; color: #344054; user-select: none; }
      .control-panel[open] summary { border-bottom: 1px solid #e8ebef; }
      .panel-grid { display: grid; grid-template-columns: 1fr; gap: 6px 8px; padding: 7px; }
      label { display: grid; grid-template-columns: 68px minmax(58px, 1fr) auto; gap: 6px; align-items: center; min-width: 0; margin: 0; color: #5f6b7a; }
      label.check-row { grid-template-columns: 1fr auto; min-height: 28px; padding: 0 7px; border: 1px solid #d7dce2; border-radius: 6px; background: #fff; color: #1d2630; }
      label.inline-check { flex: 1 1 110px; min-width: 110px; }
      select, input, button { min-height: 28px; border: 1px solid #d7dce2; border-radius: 6px; background: #fff; color: #1d2630; font: inherit; }
      input[type="checkbox"] { min-height: 0; width: 16px; height: 16px; accent-color: #1f6feb; }
      input[type="range"] { width: 100%; }
      input[type="number"] { width: 58px; padding: 0 5px; }
      select { padding: 0 8px; }
      button { padding: 0 8px; cursor: pointer; }
      button.primary, button.active { border-color: #1f6feb; background: #1f6feb; color: #fff; }
      .control-row, .button-group { display: flex; flex-wrap: wrap; gap: 6px; }
      .primary-row { margin-left: 0; width: 100%; }
      .primary-row button { width: 100%; }
      .button-group button { flex: 1 1 auto; min-width: 62px; }
      .metrics { margin: 0; }
      .metrics div { display: flex; justify-content: space-between; gap: 8px; border-bottom: 1px solid #e8ebef; padding: 4px 0; }
      .metrics div:last-child { border-bottom: 0; }
      .metrics dt { color: #5f6b7a; }
      .metrics dd { margin: 0; font-weight: 650; }
      .color-metric { display: inline-flex; align-items: center; gap: 6px; }
      .color-swatch { width: 14px; height: 14px; }
      .rmsx-viewer { position: relative; z-index: 1; display: grid; grid-template-rows: minmax(0, 1fr); min-width: 0; min-height: 0; height: 100vh; }
      .status { margin: 7px 0 0; padding: 6px 8px; border: 1px solid #d7dce2; border-radius: 6px; background: #fbfcfd; color: #13795b; font-size: 12px; line-height: 1.3; }
      .status.error { color: #b42318; border-color: #f1b4ad; background: #fff4f2; }
      .citation-note { margin: 6px 0 0; padding: 0 2px; color: #5f6b7a; font-size: 11px; line-height: 1.35; }
      .citation-note a { color: #1f6feb; text-decoration: none; }
      .citation-note a:hover { text-decoration: underline; }
      .slice-visibility { display: grid; grid-column: 1 / -1; gap: 6px; }
      .field-label { color: #5f6b7a; }
      .chips { display: flex; flex-wrap: wrap; gap: 5px; padding: 0; }
      .chip { min-width: 30px; border-radius: 999px; padding: 3px 8px; background: #fff; border: 1px solid #d7dce2; }
      .chip.active { border-color: #1f6feb; background: #eef4ff; color: #1f6feb; }
      .legend { display: grid; grid-column: 1 / -1; gap: 7px; color: #5f6b7a; font-size: 12px; }
      .bar { height: 10px; border: 1px solid #d7dce2; border-radius: 999px; background: linear-gradient(90deg, #440154, #21918c, #fde725); }
      .legend-values, .radius-legend { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 8px; align-items: center; }
      .legend-stop, .radius-stop { display: inline-flex; align-items: center; min-width: 0; gap: 6px; font-variant-numeric: tabular-nums; white-space: nowrap; }
      .legend-stop:nth-child(2), .radius-stop:nth-child(2) { justify-content: center; }
      .legend-stop:nth-child(3), .radius-stop:nth-child(3) { justify-content: flex-end; }
      .legend-swatch { width: 12px; height: 12px; flex: 0 0 auto; border-radius: 999px; border: 1px solid rgba(29, 38, 48, 0.25); background: #5f6b7a; }
      .radius-dot { width: 8px; height: 8px; flex: 0 0 auto; border-radius: 999px; border: 1px solid rgba(29, 38, 48, 0.25); background: #5f6b7a; }
      .viewport { min-height: 0; height: 100%; background: #fff; border-top: 1px solid #d7dce2; cursor: grab; }
      .viewport.local-drag-disabled { cursor: default; }
      .viewport.dragging { cursor: grabbing; }
      @media (min-width: 1280px) {
        .rmsx-app { grid-template-columns: minmax(250px, 300px) minmax(0, 1fr); }
      }
      @media (max-width: 900px) {
        .rmsx-app { grid-template-columns: 1fr; grid-template-rows: auto minmax(420px, 1fr); height: auto; min-height: 100vh; overflow: visible; }
        .rmsx-controls { min-height: auto; border-right: 0; border-bottom: 1px solid #d7dce2; overflow: visible; }
        .control-panels { grid-template-columns: 1fr; }
        .panel-grid { grid-template-columns: repeat(auto-fit, minmax(235px, 1fr)); }
        .rmsx-viewer { height: auto; min-height: calc(100vh - 240px); }
        .viewport { min-height: 360px; }
      }
      @media (max-width: 720px) {
        .controls-heading { align-items: stretch; }
        .primary-row { margin-left: 0; width: 100%; }
        .primary-row button { flex: 1 1 0; }
      }
    `;
    document.head.appendChild(style);
  }

  function setStatus(message, isError) {
    elements.status.textContent = message;
    elements.status.classList.toggle("error", Boolean(isError));
  }

  function uniqueStrings(values) {
    return [...new Set(values.filter(Boolean).map(String))];
  }

  function datasetIdCandidates() {
    return uniqueStrings([
      visualizationConfig.dataset_id,
      visualizationConfig.id,
      visualizationConfig.hda_id,
      visualizationConfig.dataset?.id,
      incoming.dataset_id,
      incoming.id
    ]);
  }

  function historyIdCandidates() {
    return uniqueStrings([
      visualizationConfig.history_id,
      visualizationConfig.historyId,
      visualizationConfig.history?.id,
      visualizationConfig.dataset?.history_id,
      visualizationConfig.dataset?.historyId,
      incoming.history_id,
      incoming.historyId,
      incoming.history?.id,
      incoming.dataset?.history_id,
      incoming.dataset?.historyId,
      URL_PARAMS.get("history_id")
    ]);
  }

  function parseManifestResponseText(text) {
    const parsed = JSON.parse(text);
    if (parsed?.schemaVersion) {
      return parsed;
    }
    const envelopeText = parsed?.item_data || parsed?.data || parsed?.contents || parsed?.content;
    if (typeof envelopeText === "string") {
      return JSON.parse(envelopeText);
    }
    return parsed;
  }

  async function fetchManifest() {
    const inlineManifest = incoming.manifest || visualizationConfig.manifest;
    if (inlineManifest?.schemaVersion) {
      state.manifestSource = "inline-harness-manifest";
      return inlineManifest;
    }
    const manifestText = document.getElementById("app")?.dataset?.manifest;
    if (manifestText) {
      state.manifestSource = "inline-harness-dataset";
      return JSON.parse(manifestText);
    }
    const datasetCandidates = datasetIdCandidates();
    const historyCandidates = historyIdCandidates();
    if (!datasetCandidates.length) {
      throw new Error("No Galaxy dataset id was provided to the RMSX Flipbook visualization.");
    }
    const urls = [];
    for (const historyId of historyCandidates) {
      const encodedHistoryId = encodeURIComponent(historyId);
      for (const datasetId of datasetCandidates) {
        if (/^https?:\/\//.test(datasetId)) {
          continue;
        }
        const encodedDatasetId = encodeURIComponent(datasetId);
        urls.push(`/api/histories/${encodedHistoryId}/contents/${encodedDatasetId}/display?to_ext=json`);
        urls.push(`/api/histories/${encodedHistoryId}/contents/${encodedDatasetId}/display`);
      }
    }
    for (const datasetId of datasetCandidates) {
      if (/^https?:\/\//.test(datasetId)) {
        urls.push(datasetId);
        continue;
      }
      const encoded = encodeURIComponent(datasetId);
      urls.push(`/api/datasets/${encoded}/get_content_as_text`);
      urls.push(`/api/datasets/${encoded}/display?to_ext=json`);
      urls.push(`/api/datasets/${encoded}/content/data`);
      urls.push(`/datasets/${encoded}/display?to_ext=json`);
      urls.push(`/datasets/${encoded}/display/?preview=True`);
      urls.push(`/api/datasets/${encoded}/display`);
    }
    const errors = [];
    for (const url of [...new Set(urls)]) {
      try {
        const response = await fetch(url, { credentials: "same-origin" });
        if (!response.ok) {
          errors.push(`${url}: ${response.status}`);
          continue;
        }
        const text = await response.text();
        state.manifestSource = url;
        return parseManifestResponseText(text);
      } catch (error) {
        errors.push(`${url}: ${error.message}`);
      }
    }
    throw new Error(`Could not load RMSX manifest from Galaxy dataset. Tried ${errors.slice(0, 4).join("; ")}`);
  }

  function validateManifest(manifest) {
    if (!manifest || manifest.schemaVersion !== SCHEMA_VERSION) {
      throw new Error(`This dataset is not an RMSX Flipbook manifest. Expected schemaVersion ${SCHEMA_VERSION}.`);
    }
    const required = ["title", "slices", "residues", "summaries", "domain", "maskSummary", "palette", "availablePalettes", "presentation", "visualMapping", "rotationModel", "molstarRenderStyle"];
    const missing = required.filter((key) => manifest[key] === undefined);
    if (missing.length) {
      throw new Error(`RMSX manifest is missing required field(s): ${missing.join(", ")}.`);
    }
    if (!Array.isArray(manifest.slices) || !manifest.slices.length) {
      throw new Error("RMSX manifest does not contain any embedded PDB slices.");
    }
    const badSlice = manifest.slices.find((slice) => !["index", "id", "label", "filename", "rmsxColumn", "pdb"].every((key) => slice[key] !== undefined));
    if (badSlice) {
      throw new Error("RMSX manifest slice entries must include index, id, label, filename, rmsxColumn, and embedded pdb text.");
    }
    if (!Array.isArray(manifest.residues) || !manifest.residues.length) {
      throw new Error("RMSX manifest does not contain residue-level RMSX values.");
    }
  }

  function loadCss(url) {
    if (!url || [...document.styleSheets].some((sheet) => sheet.href === url)) {
      return null;
    }
    const link = document.createElement("link");
    link.rel = "stylesheet";
    link.href = url;
    document.head.appendChild(link);
    return link;
  }

  function loadScript(url) {
    return new Promise((resolve, reject) => {
      if (window.molstar?.Viewer) {
        resolve(url);
        return;
      }
      const script = document.createElement("script");
      script.src = url;
      script.onload = () => resolve(url);
      script.onerror = () => reject(new Error(`Molstar failed to load from ${url}.`));
      document.head.appendChild(script);
    });
  }

  async function loadMolstarAssets() {
    const cssCandidates = [
      LOCAL_MOLSTAR_CSS,
      REPORT.molstar?.localCssUrl,
      REPORT.molstar?.cssUrl,
      CDN_MOLSTAR_CSS
    ].filter(Boolean);
    const scriptCandidates = [
      LOCAL_MOLSTAR_JS,
      REPORT.molstar?.localJsUrl,
      REPORT.molstar?.jsUrl,
      CDN_MOLSTAR_JS
    ].filter(Boolean);
    loadCss([...new Set(cssCandidates)][0]);
    const failures = [];
    for (const candidate of [...new Set(scriptCandidates)]) {
      try {
        const loaded = await loadScript(candidate);
        state.molstarAssetSource = loaded === LOCAL_MOLSTAR_JS ? "local molstar@5.4.2" : loaded;
        if (loaded !== LOCAL_MOLSTAR_JS) {
          loadCss(loaded === REPORT.molstar?.jsUrl ? REPORT.molstar?.cssUrl : CDN_MOLSTAR_CSS);
        }
        return;
      } catch (error) {
        failures.push(error.message);
      }
    }
    throw new Error(`Molstar could not be loaded from local or fallback assets. ${failures.join(" ")}`);
  }

  function clamp(value, min, max) {
    return Math.min(max, Math.max(min, value));
  }

  function formatNumber(value) {
    return Number.isFinite(value) ? value.toFixed(3) : "-";
  }

  function numericParam(name, min, max, fallback) {
    const raw = URL_PARAMS.get(name);
    if (raw === null || raw.trim() === "") {
      return fallback;
    }
    const value = Number(raw);
    return Number.isFinite(value) ? clamp(value, min, max) : fallback;
  }

  function integerParam(name, min, max, fallback) {
    const raw = URL_PARAMS.get(name);
    if (raw === null || raw.trim() === "") {
      return fallback;
    }
    const value = Math.round(Number(raw));
    return Number.isFinite(value) ? clamp(value, min, max) : fallback;
  }

  function booleanParam(name, fallback = false) {
    const raw = URL_PARAMS.get(name);
    if (raw === null || raw.trim() === "") {
      return fallback;
    }
    return ["1", "true", "yes", "on"].includes(raw.toLowerCase());
  }

  function compactNumber(value, precision = 3) {
    const numeric = Number(value);
    if (!Number.isFinite(numeric)) {
      return "";
    }
    const rounded = Number(numeric.toFixed(precision));
    return Object.is(rounded, -0) ? "0" : String(rounded);
  }

  function numbersClose(left, right, tolerance = 0.0005) {
    return Math.abs(Number(left) - Number(right)) <= tolerance;
  }

  function parentWindowIfSameOrigin() {
    const candidates = [window.top, window.parent];
    for (const candidate of candidates) {
      try {
        if (!candidate || candidate === window || !candidate.location?.href) {
          continue;
        }
        const href = candidate.location.href;
        const documentUrl = candidate.document?.URL || href;
        if (href !== "about:blank" && documentUrl !== "about:blank") {
          return candidate;
        }
      } catch (_error) {
        // Cross-origin embedding is not expected in local Galaxy, but iframe-local
        // URL persistence still works as a fallback when parent access is blocked.
      }
    }
    return null;
  }

  function initialUrlParams() {
    const params = paramsFromVisualizationConfig();
    new URLSearchParams(window.location.search).forEach((value, name) => {
      params.set(name, value);
    });
    const parentWindow = parentWindowIfSameOrigin();
    if (parentWindow) {
      new URLSearchParams(parentWindow.location.search).forEach((value, name) => {
        params.set(name, value);
      });
    }
    return params;
  }

  function paramsFromVisualizationConfig() {
    const params = new URLSearchParams();
    for (const source of visualizationConfigStateSources()) {
      for (const name of MANAGED_URL_PARAMS) {
        const value = source.value[name];
        if (value === undefined || value === null || value === "") {
          continue;
        }
        params.set(name, Array.isArray(value) ? value.join(",") : String(value));
      }
    }
    return params;
  }

  function visualizationConfigStateSources() {
    return [
      ["visualization_config.flipbookMolstarState", visualizationConfig.flipbookMolstarState],
      ["visualization_config.flipbook_molstar_state", visualizationConfig.flipbook_molstar_state],
      ["visualization_config.viewerState", visualizationConfig.viewerState],
      ["visualization_config.state", visualizationConfig.state],
      ["visualization_config.settings", visualizationConfig.settings],
      ["visualization_config.config", visualizationConfig.config],
      ["incoming.flipbookMolstarState", incoming.flipbookMolstarState],
      ["incoming.flipbook_molstar_state", incoming.flipbook_molstar_state]
    ].filter(([, value]) => value && typeof value === "object" && !Array.isArray(value))
      .map(([name, value]) => ({ name, value }));
  }

  function visualizationConfigStateSummary() {
    return visualizationConfigStateSources().map((source) => ({
      source: source.name,
      managedKeys: MANAGED_URL_PARAMS.filter((key) => source.value[key] !== undefined)
    })).filter((source) => source.managedKeys.length > 0);
  }

  function urlStateContext() {
    const parentWindow = parentWindowIfSameOrigin();
    const targetWindow = parentWindow || window;
    let href = "about:blank";
    let fallbackOrigin = "http://localhost";
    try {
      href = targetWindow.location.href;
      const scriptOrigin = scriptUrl ? new URL(scriptUrl).origin : "";
      const windowOrigin = targetWindow.location?.origin || window.location?.origin || "";
      fallbackOrigin = [scriptOrigin, windowOrigin].find((origin) => origin && origin !== "null") || fallbackOrigin;
    } catch (_error) {
      href = "about:blank";
    }
    const writable = href !== "about:blank" && /^https?:\/\//.test(href);
    return {
      targetWindow,
      url: writable ? new URL(href) : new URL("/visualizations/blank", fallbackOrigin),
      scope: writable ? (parentWindow ? "parent-visualization-url" : "iframe-url") : "about-blank-iframe-url-unavailable",
      writable
    };
  }

  function defaultLayoutName() {
    return LAYOUTS.has(REPORT.presentation?.defaultLayout) ? REPORT.presentation.defaultLayout : "tiled";
  }

  function availablePalettes() {
    const palettes = REPORT.availablePalettes || {};
    if (Object.keys(palettes).length) {
      return palettes;
    }
    return { [REPORT.palette?.name || "viridis"]: REPORT.palette?.colors || [] };
  }

  function defaultPaletteName() {
    const requested = REPORT.palette?.name || "viridis";
    return availablePalettes()[requested] ? requested : paletteNames()[0] || requested;
  }

  function defaultColorMin() {
    return Number(REPORT.visualMapping?.defaultColorMin ?? REPORT.domain.min);
  }

  function defaultColorMax() {
    return Number(REPORT.visualMapping?.defaultColorMax ?? REPORT.domain.max);
  }

  function defaultRadiusMin() {
    return Number(REPORT.visualMapping?.defaultRadiusMin ?? 0.63);
  }

  function defaultRadiusMax() {
    return Number(REPORT.visualMapping?.defaultRadiusMax ?? 3.18);
  }

  function defaultThickness() {
    return Number(REPORT.visualMapping?.defaultThicknessScale ?? 1);
  }

  function minSpacing() {
    return Number(REPORT.flipbookReference?.minimumSpacingFactor ?? 0.1);
  }

  function maxSpacing() {
    return Number(REPORT.flipbookReference?.maximumSpacingFactor ?? DEFAULT_MAX_SPACING_FACTOR);
  }

  function defaultSpacing() {
    if (REPORT.flipbookReference?.defaultSpacingMode === "auto") {
      return automaticSpacingFactor();
    }
    return Number(REPORT.flipbookReference?.defaultSpacingFactor ?? 1);
  }

  function spacingStep() {
    return Number(REPORT.keyboardShortcuts?.spacingStep ?? DEFAULT_SPACING_STEP);
  }

  function defaultTileColumns() {
    return clamp(Math.round(Number(REPORT.flipbookReference?.defaultColumns ?? REPORT.slices.length)), 1, Math.max(1, REPORT.slices.length));
  }

  function defaultRenderMode() {
    const preset = String(REPORT.molstarRenderStyle?.preset || "clean-interactive").toLowerCase();
    return RENDER_PRESETS.has(preset) ? preset : "clean-interactive";
  }

  function defaultOutline() {
    return REPORT.molstarRenderStyle?.outline !== false;
  }

  function defaultCameraMode() {
    const mode = String(REPORT.molstarRenderStyle?.cameraMode || "orthographic").toLowerCase();
    return CAMERA_MODES.has(mode) ? mode : "orthographic";
  }

  function renderModeFromParam() {
    const requested = (URL_PARAMS.get("render") || "").toLowerCase();
    if (requested === "clean") {
      return "clean-interactive";
    }
    return RENDER_PRESETS.has(requested) ? requested : defaultRenderMode();
  }

  function cameraModeFromParam() {
    const requested = (URL_PARAMS.get("cameraMode") || "").toLowerCase();
    return CAMERA_MODES.has(requested) ? requested : defaultCameraMode();
  }

  function outlineFromParam() {
    const raw = URL_PARAMS.get("outline");
    if (raw === null || raw.trim() === "") {
      return defaultOutline();
    }
    return ["1", "true", "yes", "on"].includes(raw.toLowerCase());
  }

  function initialLayout() {
    const requested = (URL_PARAMS.get("layout") || "").toLowerCase();
    return LAYOUTS.has(requested) ? requested : defaultLayoutName();
  }

  function sliceIndexFromParam(name) {
    const raw = URL_PARAMS.get(name);
    if (!raw) {
      return null;
    }
    const value = Math.round(Number(raw));
    return Number.isInteger(value) && value >= 1 && value <= REPORT.slices.length ? value - 1 : null;
  }

  function sliceIndexesFromListParam(name) {
    const raw = URL_PARAMS.get(name);
    if (!raw) {
      return [];
    }
    const indexes = [];
    const seen = new Set();
    raw.split(/[,\s]+/).forEach((entry) => {
      const requested = Number(entry);
      if (!Number.isInteger(requested) || requested < 1 || requested > REPORT.slices.length) {
        return;
      }
      const index = requested - 1;
      if (!seen.has(index)) {
        seen.add(index);
        indexes.push(index);
      }
    });
    return indexes;
  }

  function initialVisibleSliceIndexes() {
    const listed = sliceIndexesFromListParam("slices");
    if (listed.length) {
      return listed;
    }
    const explicit = ["sliceA", "sliceB", "sliceC"].map(sliceIndexFromParam).filter((index) => index !== null);
    return explicit.length ? [...new Set(explicit)] : REPORT.slices.map((_, index) => index);
  }

  function initialSliceIndex() {
    const requested = sliceIndexFromParam("slice");
    return requested === null ? 0 : requested;
  }

  function defaultResidueKey() {
    return REPORT.residues[0]?.key || "";
  }

  function initialResidueKey() {
    const requested = URL_PARAMS.get("residue");
    if (requested && REPORT.residues.some((residue) => residue.key === requested || residue.id === requested)) {
      return requested;
    }
    return defaultResidueKey();
  }

  function paletteNames() {
    return Object.keys(availablePalettes()).sort((a, b) => a.localeCompare(b));
  }

  function currentPaletteColors() {
    const palettes = availablePalettes();
    return (palettes[state.paletteName] || palettes[defaultPaletteName()] || REPORT.palette.colors || []).map((hex) => String(hex).toUpperCase());
  }

  function currentMolstarUncertaintyColors() {
    return [...currentPaletteColors()].reverse().map((hex) => Number.parseInt(hex.slice(1), 16));
  }

  function uncertaintyColorParams() {
    return {
      domain: [VISUAL_MIN, VISUAL_MAX],
      list: { kind: "interpolate", colors: currentMolstarUncertaintyColors() }
    };
  }

  function colorDomainMin() {
    return Math.min(state.colorMin, state.colorMax - 0.000001);
  }

  function colorDomainMax() {
    return Math.max(state.colorMax, state.colorMin + 0.000001);
  }

  function normalizedRmsx(value) {
    if (!Number.isFinite(value)) {
      return 0;
    }
    return clamp((value - colorDomainMin()) / Math.max(0.000001, colorDomainMax() - colorDomainMin()), 0, 1);
  }

  function wormRadiusMin() {
    return Math.max(0.01, Math.min(state.radiusMin, state.radiusMax - 0.01) * state.thickness);
  }

  function wormRadiusMax() {
    return Math.max(wormRadiusMin() + 0.01, Math.max(state.radiusMax, state.radiusMin + 0.01) * state.thickness);
  }

  function wormRadiusSpan() {
    return Math.max(0.01, wormRadiusMax() - wormRadiusMin());
  }

  function maskedKeys() {
    return new Set(REPORT.maskSummary?.maskedKeys || []);
  }

  function residueKey(chainId, residueId) {
    return chainId ? `${chainId}:${residueId}` : residueId;
  }

  function isMasked(chainId, residueId) {
    const keys = maskedKeys();
    return keys.has(residueKey(chainId, residueId)) || keys.has(residueId);
  }

  function residueByKey() {
    const map = new Map();
    REPORT.residues.forEach((residue) => {
      map.set(residue.key, residue);
      map.set(residue.id, residue);
    });
    return map;
  }

  function selectedResidue() {
    return REPORT?.residues?.find((residue) => residue.key === state.selectedResidueKey || residue.id === state.selectedResidueKey) || REPORT?.residues?.[0] || null;
  }

  function selectedResidueRmsx() {
    const slice = REPORT?.slices?.[state.currentIndex];
    const residue = selectedResidue();
    const value = slice && residue ? Number(residue.values?.[slice.rmsxColumn]) : NaN;
    return Number.isFinite(value) ? value : NaN;
  }

  function visualRadiusForRmsx(value) {
    if (!Number.isFinite(value)) {
      return NaN;
    }
    return wormRadiusMin() + (wormRadiusSpan() * normalizedRmsx(value));
  }

  function selectedResidueColor() {
    const value = selectedResidueRmsx();
    if (!Number.isFinite(value)) {
      return "-";
    }
    return expectedColorForNormalizedRmsx(normalizedRmsx(value));
  }

  function structureStats(pdb) {
    const stats = { minX: Infinity, maxX: -Infinity, minY: Infinity, maxY: -Infinity, minZ: Infinity, maxZ: -Infinity, sumX: 0, sumY: 0, sumZ: 0, count: 0 };
    pdb.split(/\r?\n/).forEach((line) => {
      if (!line.startsWith("ATOM") && !line.startsWith("HETATM")) {
        return;
      }
      const x = Number(line.slice(30, 38));
      const y = Number(line.slice(38, 46));
      const z = Number(line.slice(46, 54));
      if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
        return;
      }
      stats.minX = Math.min(stats.minX, x);
      stats.maxX = Math.max(stats.maxX, x);
      stats.minY = Math.min(stats.minY, y);
      stats.maxY = Math.max(stats.maxY, y);
      stats.minZ = Math.min(stats.minZ, z);
      stats.maxZ = Math.max(stats.maxZ, z);
      stats.sumX += x;
      stats.sumY += y;
      stats.sumZ += z;
      stats.count += 1;
    });
    if (!stats.count) {
      return { ...stats, width: 30, height: 30, depth: 30, center: { x: 0, y: 0, z: 0 } };
    }
    return {
      ...stats,
      width: Math.max(1, stats.maxX - stats.minX),
      height: Math.max(1, stats.maxY - stats.minY),
      depth: Math.max(1, stats.maxZ - stats.minZ),
      center: { x: stats.sumX / stats.count, y: stats.sumY / stats.count, z: stats.sumZ / stats.count }
    };
  }

  function degreesToRadians(value) {
    return value * Math.PI / 180;
  }

  function radiansToDegrees(value) {
    return value * 180 / Math.PI;
  }

  function wrapAngle(value) {
    if (!Number.isFinite(Number(value))) {
      return 0;
    }
    let angle = Number(value);
    while (angle > 180) {
      angle -= 360;
    }
    while (angle < -180) {
      angle += 360;
    }
    return angle;
  }

  function multiplyMatrices(left, right) {
    return left.map((row) => right[0].map((_, columnIndex) => row.reduce((sum, value, index) => sum + value * right[index][columnIndex], 0)));
  }

  function rotationMatrixFor(rotation) {
    const x = degreesToRadians(Number(rotation.x ?? 90));
    const y = degreesToRadians(Number(rotation.y ?? 0));
    const z = degreesToRadians(Number(rotation.z ?? 0));
    const sx = Math.sin(x);
    const cx = Math.cos(x);
    const sy = Math.sin(y);
    const cy = Math.cos(y);
    const sz = Math.sin(z);
    const cz = Math.cos(z);
    const rx = [[1, 0, 0], [0, cx, -sx], [0, sx, cx]];
    const ry = [[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]];
    const rz = [[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]];
    return multiplyMatrices(rz, multiplyMatrices(ry, rx));
  }

  function rotationMatrix() {
    if (!state.rotationMatrix) {
      syncRotationMatrixFromEuler();
    }
    return state.rotationMatrix;
  }

  function defaultRotationMatrix() {
    return rotationMatrixFor(REPORT.rotationModel?.defaultRotation || { x: 90, y: 0, z: 0 });
  }

  function visualEnvelope() {
    const layoutMatrix = defaultRotationMatrix();
    const projectedWidths = REPORT.slices.map((slice) => {
      const stats = structureStats(slice.pdb);
      return rotatedExtentX(stats, layoutMatrix);
    });
    return (Math.max(30, ...projectedWidths) + visualRadiusPadding()) * tilePaddingFactor();
  }

  function automaticSpacingFactor() {
    const layoutMatrix = defaultRotationMatrix();
    const widestSlice = Math.max(30, ...REPORT.slices.map((slice) =>
      rotatedExtentX(structureStats(slice.pdb), layoutMatrix)
    ));
    const gapFraction = Number(REPORT.flipbookReference?.autoSpacingGapFraction ?? 0.10);
    const minimumGap = Number(REPORT.flipbookReference?.autoSpacingMinimumGap ?? 8);
    const targetSlot = widestSlice + Math.max(minimumGap, widestSlice * gapFraction);
    return clamp(targetSlot / visualEnvelope(), minSpacing(), maxSpacing());
  }

  function rotatedExtentX(stats, matrix) {
    const corners = [
      [stats.minX, stats.minY, stats.minZ],
      [stats.minX, stats.minY, stats.maxZ],
      [stats.minX, stats.maxY, stats.minZ],
      [stats.minX, stats.maxY, stats.maxZ],
      [stats.maxX, stats.minY, stats.minZ],
      [stats.maxX, stats.minY, stats.maxZ],
      [stats.maxX, stats.maxY, stats.minZ],
      [stats.maxX, stats.maxY, stats.maxZ]
    ].map(([x, y, z]) => transformPoint(matrix, stats.center, stats.center, x, y, z).x);
    return Math.max(1, Math.max(...corners) - Math.min(...corners));
  }

  function tilePaddingFactor() {
    return Math.max(1, Number(REPORT.flipbookReference?.tilePaddingFactor ?? 1.55));
  }

  function visualRadiusPadding() {
    return Math.max(24, wormRadiusMax() * 8 + 12);
  }

  function tileOffset(index) {
    if (state.layout !== "tiled") {
      return { x: 0, y: 0, z: 0 };
    }
    const columns = clamp(Math.round(state.columns), 1, REPORT.slices.length);
    const rows = Math.ceil(REPORT.slices.length / columns);
    const row = Math.floor(index / columns);
    const column = index % columns;
    const rowLength = row === rows - 1 ? REPORT.slices.length - row * columns : columns;
    const slot = visualEnvelope() * state.spacing;
    return {
      x: (column - (rowLength - 1) / 2) * slot,
      y: ((rows - 1) / 2 - row) * slot * 0.82,
      z: 0
    };
  }

  function transformPoint(matrix, center, target, x, y, z) {
    const lx = x - center.x;
    const ly = y - center.y;
    const lz = z - center.z;
    return {
      x: target.x + matrix[0][0] * lx + matrix[0][1] * ly + matrix[0][2] * lz,
      y: target.y + matrix[1][0] * lx + matrix[1][1] * ly + matrix[1][2] * lz,
      z: target.z + matrix[2][0] * lx + matrix[2][1] * ly + matrix[2][2] * lz
    };
  }

  function identityRotationMatrix() {
    return [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  }

  function cloneRotationMatrix(matrix) {
    return (matrix || identityRotationMatrix()).map((row) => row.slice());
  }

  function eulerFromRotationMatrix(matrix) {
    const m = matrix || identityRotationMatrix();
    const y = Math.asin(clamp(-m[2][0], -1, 1));
    const cy = Math.cos(y);
    let x = 0;
    let z = 0;
    if (Math.abs(cy) > 0.000001) {
      x = Math.atan2(m[2][1], m[2][2]);
      z = Math.atan2(m[1][0], m[0][0]);
    } else {
      z = Math.atan2(-m[0][1], m[1][1]);
    }
    return {
      x: radiansToDegrees(x),
      y: radiansToDegrees(y),
      z: radiansToDegrees(z)
    };
  }

  function setRotationMatrix(matrix, options = {}) {
    state.rotationMatrix = cloneRotationMatrix(matrix);
    if (options.updateEuler !== false) {
      const euler = eulerFromRotationMatrix(state.rotationMatrix);
      state.rotation = {
        x: wrapAngle(euler.x),
        y: wrapAngle(euler.y),
        z: wrapAngle(euler.z)
      };
    }
  }

  function syncRotationMatrixFromEuler() {
    state.rotationMatrix = rotationMatrixFor(state.rotation);
  }

  function vectorFromArray(value, fallback) {
    const vector = {
      x: Number(value?.[0]),
      y: Number(value?.[1]),
      z: Number(value?.[2])
    };
    if (Number.isFinite(vector.x) && Number.isFinite(vector.y) && Number.isFinite(vector.z)) {
      return vector;
    }
    return { ...fallback };
  }

  function subtractVectors(left, right) {
    return {
      x: left.x - right.x,
      y: left.y - right.y,
      z: left.z - right.z
    };
  }

  function crossVectors(left, right) {
    return {
      x: (left.y * right.z) - (left.z * right.y),
      y: (left.z * right.x) - (left.x * right.z),
      z: (left.x * right.y) - (left.y * right.x)
    };
  }

  function normalizeVector(vector, fallback) {
    const length = Math.sqrt((vector.x * vector.x) + (vector.y * vector.y) + (vector.z * vector.z));
    if (!Number.isFinite(length) || length < 0.000001) {
      return { ...fallback };
    }
    return {
      x: vector.x / length,
      y: vector.y / length,
      z: vector.z / length
    };
  }

  function cameraSnapshot() {
    const camera = viewer?.plugin?.canvas3d?.camera;
    try {
      return camera?.getSnapshot?.() || camera?.state || null;
    } catch (error) {
      console.debug("Molstar camera snapshot unavailable for local rotation.", error);
      return null;
    }
  }

  function currentScreenRotationAxes() {
    const snapshot = cameraSnapshot();
    const position = vectorFromArray(snapshot?.position, { x: 0, y: 0, z: 100 });
    const target = vectorFromArray(snapshot?.target, { x: 0, y: 0, z: 0 });
    const view = normalizeVector(subtractVectors(target, position), { x: 0, y: 0, z: -1 });
    let up = normalizeVector(vectorFromArray(snapshot?.up, { x: 0, y: 1, z: 0 }), { x: 0, y: 1, z: 0 });
    const right = normalizeVector(crossVectors(view, up), { x: 1, y: 0, z: 0 });
    up = normalizeVector(crossVectors(right, view), up);
    return { right, up, view };
  }

  function axisAngleRotationMatrix(axis, degrees) {
    const unit = normalizeVector(axis, { x: 0, y: 1, z: 0 });
    const angle = degreesToRadians(degrees);
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    const t = 1 - c;
    const x = unit.x;
    const y = unit.y;
    const z = unit.z;
    return [
      [(t * x * x) + c, (t * x * y) - (s * z), (t * x * z) + (s * y)],
      [(t * y * x) + (s * z), (t * y * y) + c, (t * y * z) - (s * x)],
      [(t * z * x) - (s * y), (t * z * y) + (s * x), (t * z * z) + c]
    ];
  }

  function rotationDeltaMatrixForScreenDrag(dx, dy, axes) {
    const horizontal = axisAngleRotationMatrix(axes.up, dx * state.rotationSensitivity);
    const vertical = axisAngleRotationMatrix(axes.right, dy * state.rotationSensitivity);
    return multiplyMatrices(vertical, horizontal);
  }

  function roundedVector(vector) {
    return {
      x: Number(vector.x.toFixed(4)),
      y: Number(vector.y.toFixed(4)),
      z: Number(vector.z.toFixed(4))
    };
  }

  function applyScreenRotationDrag(dx, dy, axes) {
    const delta = rotationDeltaMatrixForScreenDrag(dx, dy, axes);
    setRotationMatrix(multiplyMatrices(delta, rotationMatrix()));
    state.lastRotationGesture = {
      mode: "screen-axis delta, per-slice local pivot",
      dx: Number(dx.toFixed(2)),
      dy: Number(dy.toFixed(2)),
      sensitivity: state.rotationSensitivity,
      axes: {
        right: roundedVector(axes.right),
        up: roundedVector(axes.up),
        view: roundedVector(axes.view)
      }
    };
  }

  function transformedPdb(slice, index, mode, applySceneTransform = true) {
    const stats = structureStats(slice.pdb);
    const anchor = structureStats(REPORT.slices[0].pdb).center;
    const offset = tileOffset(index);
    const target = { x: anchor.x + offset.x, y: anchor.y + offset.y, z: anchor.z + offset.z };
    const matrix = rotationMatrix();
    const residues = residueByKey();
    let atomCount = 0;
    const lines = slice.pdb.split(/\r?\n/).map((line) => {
      if (!line.startsWith("ATOM") && !line.startsWith("HETATM")) {
        return line;
      }
      const padded = line.padEnd(80, " ");
      const x = Number(padded.slice(30, 38));
      const y = Number(padded.slice(38, 46));
      const z = Number(padded.slice(46, 54));
      const chainId = padded.slice(21, 22).trim();
      const residueId = padded.slice(22, 26).trim();
      const masked = isMasked(chainId, residueId);
      if (mode === "unmasked" && masked) {
        return null;
      }
      if (mode === "masked" && !masked) {
        return null;
      }
      const residue = residues.get(residueKey(chainId, residueId)) || residues.get(residueId);
      const rmsx = residue?.values?.[slice.rmsxColumn];
      const bfactor = normalizedRmsx(Number(rmsx));
      const point = applySceneTransform && Number.isFinite(x) && Number.isFinite(y) && Number.isFinite(z)
        ? transformPoint(matrix, stats.center, target, x, y, z)
        : null;
      atomCount += 1;
      const xText = point ? point.x.toFixed(3).padStart(8) : padded.slice(30, 38);
      const yText = point ? point.y.toFixed(3).padStart(8) : padded.slice(38, 46);
      const zText = point ? point.z.toFixed(3).padStart(8) : padded.slice(46, 54);
      const bText = bfactor.toFixed(2).padStart(6);
      return `${padded.slice(0, 30)}${xText}${yText}${zText}${padded.slice(54, 60)}${bText}${padded.slice(66)}`.trimEnd();
    }).filter((line) => line !== null);
    return { pdb: lines.join("\n"), atomCount };
  }

  function sceneTransformForSlice(slice, index) {
    const stats = structureStats(slice.pdb);
    const anchor = structureStats(REPORT.slices[0].pdb).center;
    const offset = tileOffset(index);
    const target = { x: anchor.x + offset.x, y: anchor.y + offset.y, z: anchor.z + offset.z };
    const matrix = rotationMatrix();
    return { matrix, center: stats.center, target };
  }

  function molstarTransformForSlice(slice, index) {
    const transform = sceneTransformForSlice(slice, index);
    const r = transform.matrix;
    const c = transform.center;
    const t = transform.target;
    const tx = t.x - ((r[0][0] * c.x) + (r[0][1] * c.y) + (r[0][2] * c.z));
    const ty = t.y - ((r[1][0] * c.x) + (r[1][1] * c.y) + (r[1][2] * c.z));
    const tz = t.z - ((r[2][0] * c.x) + (r[2][1] * c.y) + (r[2][2] * c.z));
    return [
      r[0][0], r[1][0], r[2][0], 0,
      r[0][1], r[1][1], r[2][1], 0,
      r[0][2], r[1][2], r[2][2], 0,
      tx, ty, tz, 1
    ];
  }

  function selectedResiduePdb(slice, index, applySceneTransform = true) {
    if (!state.marker || !state.selectedResidueKey) {
      return { pdb: "", atomCount: 0 };
    }
    const selected = REPORT.residues.find((residue) => residue.key === state.selectedResidueKey) || REPORT.residues[0];
    const stats = structureStats(slice.pdb);
    const anchor = structureStats(REPORT.slices[0].pdb).center;
    const offset = tileOffset(index);
    const target = { x: anchor.x + offset.x, y: anchor.y + offset.y, z: anchor.z + offset.z };
    const matrix = rotationMatrix();
    let atomCount = 0;
    const lines = slice.pdb.split(/\r?\n/).map((line) => {
      if (!line.startsWith("ATOM") && !line.startsWith("HETATM")) {
        return null;
      }
      const padded = line.padEnd(80, " ");
      const chainId = padded.slice(21, 22).trim();
      const residueId = padded.slice(22, 26).trim();
      if (selected.chain && chainId !== selected.chain) {
        return null;
      }
      if (residueId !== selected.id && residueId !== selected.key) {
        return null;
      }
      const x = Number(padded.slice(30, 38));
      const y = Number(padded.slice(38, 46));
      const z = Number(padded.slice(46, 54));
      const point = applySceneTransform && Number.isFinite(x) && Number.isFinite(y) && Number.isFinite(z)
        ? transformPoint(matrix, stats.center, target, x, y, z)
        : null;
      atomCount += 1;
      const xText = point ? point.x.toFixed(3).padStart(8) : padded.slice(30, 38);
      const yText = point ? point.y.toFixed(3).padStart(8) : padded.slice(38, 46);
      const zText = point ? point.z.toFixed(3).padStart(8) : padded.slice(46, 54);
      return `${padded.slice(0, 30)}${xText}${yText}${zText}${padded.slice(54)}`.trimEnd();
    }).filter(Boolean);
    return { pdb: lines.join("\n"), atomCount };
  }

  function estimatedVisualRadius() {
    return Math.max(1, wormRadiusMax());
  }

  function sceneBoundsForEntry(entry) {
    const stats = structureStats(entry.slice.pdb);
    if (!stats.count) {
      return null;
    }
    const transform = sceneTransformForSlice(entry.slice, entry.index);
    const radius = estimatedVisualRadius();
    const corners = [
      [stats.minX, stats.minY, stats.minZ],
      [stats.minX, stats.minY, stats.maxZ],
      [stats.minX, stats.maxY, stats.minZ],
      [stats.minX, stats.maxY, stats.maxZ],
      [stats.maxX, stats.minY, stats.minZ],
      [stats.maxX, stats.minY, stats.maxZ],
      [stats.maxX, stats.maxY, stats.minZ],
      [stats.maxX, stats.maxY, stats.maxZ]
    ].map(([x, y, z]) => transformPoint(transform.matrix, transform.center, transform.target, x, y, z));
    const bounds = corners.reduce((acc, point) => ({
      minX: Math.min(acc.minX, point.x),
      maxX: Math.max(acc.maxX, point.x),
      minY: Math.min(acc.minY, point.y),
      maxY: Math.max(acc.maxY, point.y),
      minZ: Math.min(acc.minZ, point.z),
      maxZ: Math.max(acc.maxZ, point.z)
    }), {
      minX: Infinity,
      maxX: -Infinity,
      minY: Infinity,
      maxY: -Infinity,
      minZ: Infinity,
      maxZ: -Infinity
    });
    return {
      minX: bounds.minX - radius,
      maxX: bounds.maxX + radius,
      minY: bounds.minY - radius,
      maxY: bounds.maxY + radius,
      minZ: bounds.minZ - radius,
      maxZ: bounds.maxZ + radius
    };
  }

  function focusEntries() {
    return activeEntries();
  }

  function sceneFocusSphere(entries = focusEntries()) {
    const sceneBounds = entries.map(sceneBoundsForEntry).filter(Boolean);
    if (!sceneBounds.length) {
      return null;
    }
    const bounds = sceneBounds.reduce((acc, current) => ({
      minX: Math.min(acc.minX, current.minX),
      maxX: Math.max(acc.maxX, current.maxX),
      minY: Math.min(acc.minY, current.minY),
      maxY: Math.max(acc.maxY, current.maxY),
      minZ: Math.min(acc.minZ, current.minZ),
      maxZ: Math.max(acc.maxZ, current.maxZ)
    }), {
      minX: Infinity,
      maxX: -Infinity,
      minY: Infinity,
      maxY: -Infinity,
      minZ: Infinity,
      maxZ: -Infinity
    });
    if (!Number.isFinite(bounds.minX) || !Number.isFinite(bounds.maxX)) {
      return null;
    }
    const center = [
      (bounds.minX + bounds.maxX) / 2,
      (bounds.minY + bounds.maxY) / 2,
      (bounds.minZ + bounds.maxZ) / 2
    ];
    const dx = bounds.maxX - bounds.minX;
    const dy = bounds.maxY - bounds.minY;
    const dz = bounds.maxZ - bounds.minZ;
    return {
      center,
      radius: Math.max(1, Math.sqrt((dx * dx) + (dy * dy) + (dz * dz)) / 2)
    };
  }

  function cameraFocusExtraRadius(sphere) {
    return Math.max(4, sphere.radius * 0.9);
  }

  function tiledPlacementSummary(focusSphere = null) {
    const envelope = visualEnvelope();
    const slot = envelope * state.spacing;
    return {
      envelope: Number(envelope.toFixed(3)),
      slot: Number(slot.toFixed(3)),
      spacing: state.spacing,
      columns: state.columns,
      tilePaddingFactor: tilePaddingFactor(),
      visualRadiusPadding: Number(visualRadiusPadding().toFixed(3)),
      cameraExtraRadius: focusSphere ? Number(cameraFocusExtraRadius(focusSphere).toFixed(3)) : null
    };
  }

  function hexColorToMolstarNumber(hex, fallback = 0xffffff) {
    const normalized = String(hex || "").replace("#", "");
    if (!/^[0-9a-fA-F]{6}$/.test(normalized)) {
      return fallback;
    }
    return Number.parseInt(normalized, 16);
  }

  function molstarCanvasProps(options = {}) {
    const style = REPORT.molstarRenderStyle || {};
    const softRender = ["soft", "studio", "cinematic"].includes(state.renderMode);
    const occlusionEnabled = !options.interactive && (softRender ? style.ambientOcclusion !== "never" : style.ambientOcclusion === true);
    const illuminationEnabled = !options.interactive && (softRender ? style.illumination !== "never" : style.illumination === true);
    return {
      transparentBackground: false,
      dpoitIterations: 1,
      userInteractionReleaseMs: 0,
      multiSample: {
        mode: "off",
        sampleLevel: 0,
        reduceFlicker: false,
        reuseOcclusion: false
      },
      camera: {
        mode: state.cameraMode,
        fov: 45,
        manualReset: false
      },
      cameraFog: { name: "off", params: {} },
      renderer: {
        backgroundColor: hexColorToMolstarNumber(style.backgroundColor, 0xffffff),
        ambientIntensity: 0.78
      },
      postprocessing: {
        enabled: state.outline || occlusionEnabled,
        outline: state.outline
          ? { name: "on", params: { scale: 0.55, threshold: 0.22, color: 0x1f2937, includeTransparent: true } }
          : { name: "off", params: {} },
        occlusion: occlusionEnabled
          ? {
              name: "on",
              params: {
                samples: 8,
                multiScale: { name: "off", params: {} },
                radius: 3.2,
                bias: 0.85,
                blurKernelSize: 11,
                blurDepthBias: 0.5,
                resolutionScale: 0.5,
                color: 0x000000,
                transparentThreshold: 0.4
              }
            }
          : { name: "off", params: {} },
        antialiasing: { name: "smaa", params: {} },
        shadow: { name: "off", params: {} },
        dof: { name: "off", params: {} },
        sharpening: { name: "off", params: {} },
        bloom: { name: "off", params: {} },
        background: { variant: { name: "off", params: {} } }
      },
      marking: {
        enabled: false,
        highlightEdgeColor: 0x000000,
        selectEdgeColor: 0x000000,
        ghostEdgeStrength: 0,
        innerEdgeFactor: 1
      },
      illumination: {
        enabled: illuminationEnabled,
        maxIterations: 4,
        denoise: true
      }
    };
  }

  function applyMolstarRenderStyle(options = {}) {
    const plugin = viewer?.plugin;
    if (!plugin?.canvas3d?.setProps) {
      return false;
    }
    try {
      plugin.canvas3d.setProps(molstarCanvasProps(options));
      plugin.canvas3d.requestDraw?.();
      return true;
    } catch (error) {
      console.warn("RMSX Flipbook render style could not be applied.", error);
      return false;
    }
  }

  function canvasRenderStyleSummary() {
    const props = viewer?.plugin?.canvas3d?.props || {};
    const post = props.postprocessing || {};
    return {
      preset: state.renderMode,
      outline: post.outline?.name || (state.outline ? "on" : "off"),
      cameraMode: props.camera?.mode || viewer?.plugin?.canvas3d?.camera?.state?.mode || state.cameraMode,
      occlusion: post.occlusion?.name || "off",
      illumination: Boolean(props.illumination?.enabled),
      cameraFog: props.cameraFog?.name || "off",
      multiSampleMode: props.multiSample?.mode || "off",
      ambientIntensity: props.renderer?.ambientIntensity ?? null,
      backgroundColor: props.renderer?.backgroundColor ?? null
    };
  }

  async function createViewer() {
    elements.viewport.replaceChildren();
    viewer = await window.molstar.Viewer.create("molstarViewport", {
      layoutIsExpanded: false,
      layoutShowControls: false,
      layoutShowRemoteState: false,
      layoutShowSequence: false,
      layoutShowLog: false,
      layoutShowLeftPanel: false,
      viewportShowExpand: true,
      viewportShowSelectionMode: false,
      viewportShowAnimation: false
    });
    applyMolstarRenderStyle();
    setupViewportResizeObserver();
  }

  function requestMolstarDraw() {
    try {
      viewer?.handleResize?.();
      viewer?.plugin?.layout?.events?.updated?.next?.(void 0);
      window.dispatchEvent(new Event("resize"));
      viewer?.plugin?.canvas3d?.requestDraw?.();
    } catch (error) {
      console.debug("Molstar resize/draw refresh failed.", error);
    }
  }

  function schedulePostLayoutReset() {
    if (!state.loaded) {
      return;
    }
    const resetAfterLayout = () => {
      requestMolstarDraw();
      resetView();
      updateDiagnostics();
    };
    const scheduleFrame = window.requestAnimationFrame || ((callback) => window.setTimeout(callback, 16));
    scheduleFrame(() => {
      scheduleFrame(() => {
        resetAfterLayout();
        window.setTimeout(resetAfterLayout, 180);
        window.setTimeout(resetAfterLayout, 700);
      });
    });
  }

  function setupViewportResizeObserver() {
    if (resizeObserver || typeof ResizeObserver === "undefined" || !elements.viewport) {
      return;
    }
    let lastWidth = 0;
    let lastHeight = 0;
    resizeObserver = new ResizeObserver((entries) => {
      const rect = entries?.[0]?.contentRect;
      if (!rect) {
        return;
      }
      const widthChanged = Math.abs(rect.width - lastWidth) > 2;
      const heightChanged = Math.abs(rect.height - lastHeight) > 2;
      lastWidth = rect.width;
      lastHeight = rect.height;
      if (!state.loaded || (!widthChanged && !heightChanged)) {
        return;
      }
      window.clearTimeout(resizeResetTimer);
      resizeResetTimer = window.setTimeout(() => {
        requestMolstarDraw();
        resetView();
        updateDiagnostics();
      }, 120);
    });
    resizeObserver.observe(elements.viewport);
  }

  async function addStructure(plugin, pdb, label, alpha, marker) {
    if (!pdb.trim()) {
      return null;
    }
    const data = await plugin.builders.data.rawData({ data: pdb, label });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, "pdb");
    const model = await plugin.builders.structure.createModel(trajectory);
    const structure = await plugin.builders.structure.createStructure(model);
    if (marker) {
      return plugin.builders.structure.representation.addRepresentation(structure, {
        type: "spacefill",
        typeParams: { sizeFactor: 0.36, alpha: 0.86, quality: "high" },
        color: "uniform",
        colorParams: { value: 0x111827 }
      });
    }
    const rep = {
      type: "putty",
      typeParams: { sizeFactor: 1, quality: "high", alpha },
      color: "uncertainty",
      colorParams: uncertaintyColorParams(),
      size: "uncertainty",
      sizeParams: { bfactorFactor: wormRadiusSpan(), rmsfFactor: 0, baseSize: wormRadiusMin() }
    };
    try {
      const representation = await plugin.builders.structure.representation.addRepresentation(structure, rep);
      state.representationMode = "putty";
      return representation;
    } catch (error) {
      const representation = await plugin.builders.structure.representation.addRepresentation(structure, {
        ...rep,
        type: "cartoon",
        typeParams: { aspectRatio: 1.2, sizeFactor: Math.max(0.22, wormRadiusMax() / 2.7), quality: "high", alpha }
      });
      state.representationMode = "cartoon";
      return representation;
    }
  }

  function representationObject(representation) {
    return representation?.cell?.obj?.data?.repr
      || representation?.obj?.data?.repr
      || representation?.data?.repr
      || null;
  }

  function recordRepresentation(record) {
    return representationObject(record?.representation);
  }

  function activeEntries() {
    return REPORT.slices.map((slice, index) => ({ slice, index })).filter((entry) => state.visible.has(entry.index));
  }

  function allEntries() {
    return REPORT.slices.map((slice, index) => ({ slice, index }));
  }

  function isSliceVisible(index) {
    return state.visible.has(index);
  }

  function disposeViewer() {
    if (viewer?.dispose) {
      viewer.dispose();
    } else if (viewer?.plugin?.dispose) {
      viewer.plugin.dispose();
    }
    viewer = null;
    state.records = [];
    state.loaded = false;
    state.liveTransforms = false;
    elements.viewport.replaceChildren();
  }

  function flushMolstarDraw(representations, fast = false) {
    const canvas = viewer?.plugin?.canvas3d;
    if (!canvas) {
      return;
    }
    for (const repr of representations) {
      try {
        if (repr && typeof canvas.update === "function") {
          canvas.update(repr, fast);
        }
      } catch (error) {
        console.debug("Molstar representation update failed.", error);
      }
    }
    try {
      canvas.commit?.(fast);
      canvas.requestDraw?.();
    } catch (error) {
      console.debug("Molstar draw flush failed.", error);
    }
  }

  function applyRecordTransform(record) {
    const repr = recordRepresentation(record);
    if (!repr?.setState) {
      return null;
    }
    const visible = isSliceVisible(record.index) && (record.kind !== "marker" || state.marker);
    repr.setState({
      transform: molstarTransformForSlice(record.slice, record.index),
      visible,
      pickable: visible
    });
    return repr;
  }

  function applyLiveTransforms(autoView = false, fast = false) {
    const updated = [];
    for (const record of state.records) {
      const repr = applyRecordTransform(record);
      if (repr) {
        updated.push(repr);
      }
    }
    state.liveTransforms = updated.length > 0;
    if (updated.length) {
      flushMolstarDraw(updated, fast);
    }
    if (autoView) {
      resetView();
    }
    setLoadedSceneStatus();
    updateMetrics();
    return updated.length;
  }

  async function loadLiveScene(autoView) {
    disposeViewer();
    await createViewer();
    setStatus(`Loading ${REPORT.slices.length} RMSX slices in one native Molstar scene...`);
    const plugin = viewer.plugin;
    const hasMask = (REPORT.maskSummary?.maskedKeys || []).length > 0;
    for (const entry of allEntries()) {
      if (hasMask) {
        const unmasked = transformedPdb(entry.slice, entry.index, "unmasked", false);
        const unmaskedRep = await addStructure(plugin, unmasked.pdb, `${entry.slice.label} unmasked`, 1, false);
        state.records.push({ ...entry, kind: "unmasked", representation: unmaskedRep });
        const masked = transformedPdb(entry.slice, entry.index, "masked", false);
        const maskedRep = await addStructure(plugin, masked.pdb, `${entry.slice.label} masked`, REPORT.maskOpacity || 0.3, false);
        state.records.push({ ...entry, kind: "masked", representation: maskedRep });
      } else {
        const all = transformedPdb(entry.slice, entry.index, "all", false);
        const representation = await addStructure(plugin, all.pdb, entry.slice.label, 1, false);
        state.records.push({ ...entry, kind: "all", representation });
      }
      const marker = selectedResiduePdb(entry.slice, entry.index, false);
      const markerRep = await addStructure(plugin, marker.pdb, `${entry.slice.label} selected residue`, 0.86, true);
      if (markerRep) {
        state.records.push({ ...entry, kind: "marker", representation: markerRep });
      }
    }
    state.loaded = true;
    const updated = applyLiveTransforms(autoView);
    if (!updated) {
      state.forceCoordinateFallback = true;
      await renderCoordinateScene(autoView);
    } else if (autoView !== false) {
      schedulePostLayoutReset();
    }
  }

  async function renderCoordinateScene(autoView) {
    const currentToken = ++renderToken;
    setStatus("Rendering RMSX slices in Molstar...");
    disposeViewer();
    await createViewer();
    if (currentToken !== renderToken) {
      return;
    }
    const plugin = viewer.plugin;
    const hasMask = (REPORT.maskSummary?.maskedKeys || []).length > 0;
    for (const entry of activeEntries()) {
      if (hasMask) {
        const unmasked = transformedPdb(entry.slice, entry.index, "unmasked", true);
        await addStructure(plugin, unmasked.pdb, `${entry.slice.label} unmasked`, 1, false);
        const masked = transformedPdb(entry.slice, entry.index, "masked", true);
        await addStructure(plugin, masked.pdb, `${entry.slice.label} masked`, REPORT.maskOpacity || 0.3, false);
      } else {
        const all = transformedPdb(entry.slice, entry.index, "all", true);
        await addStructure(plugin, all.pdb, entry.slice.label, 1, false);
      }
      const marker = selectedResiduePdb(entry.slice, entry.index, true);
      await addStructure(plugin, marker.pdb, `${entry.slice.label} selected residue`, 0.86, true);
    }
    if (autoView !== false) {
      resetView();
      schedulePostLayoutReset();
    }
    setLoadedSceneStatus();
    updateMetrics();
  }

  async function renderScene(autoView) {
    if (!REPORT) {
      return;
    }
    if (!state.forceCoordinateFallback && state.loaded && state.liveTransforms) {
      applyLiveTransforms(autoView !== false);
      return;
    }
    if (!state.forceCoordinateFallback) {
      await loadLiveScene(autoView !== false);
      return;
    }
    await renderCoordinateScene(autoView !== false);
  }

  function resetView() {
    const plugin = viewer?.plugin;
    const sphere = sceneFocusSphere();
    if (sphere && plugin?.managers?.camera?.focusSphere) {
      plugin.managers.camera.focusSphere(sphere, { durationMs: 0, extraRadius: cameraFocusExtraRadius(sphere) });
    } else if (plugin?.managers?.camera?.reset) {
      plugin.managers.camera.reset();
    } else if (plugin?.canvas3d?.requestCameraReset) {
      plugin.canvas3d.requestCameraReset();
    }
  }

  function hasMaskedResidues() {
    return (REPORT.maskSummary?.maskedKeys || []).length > 0;
  }

  function setLoadedSceneStatus() {
    if (!REPORT) {
      return;
    }
    const visibleCount = visibleSliceIndexes().length;
    const maskText = hasMaskedResidues() ? `; ${Number(REPORT.maskSummary?.maskedResidues ?? REPORT.maskSummary?.maskedKeys?.length ?? 0)} masked` : "";
    setStatus(`${visibleCount}/${REPORT.slices.length} slices visible; ${state.paletteName}; ${state.representationMode}${maskText}.`);
  }

  function sequenceRmsxStats() {
    const stats = { min: Infinity, max: -Infinity, sum: 0, count: 0, peakResidue: "-" };
    for (const residue of REPORT?.residues || []) {
      for (const value of Object.values(residue.values || {})) {
        const numeric = Number(value);
        if (!Number.isFinite(numeric)) {
          continue;
        }
        stats.min = Math.min(stats.min, numeric);
        stats.sum += numeric;
        stats.count += 1;
        if (numeric > stats.max) {
          stats.max = numeric;
          stats.peakResidue = residue.label || residue.key || "-";
        }
      }
    }
    stats.mean = stats.count ? stats.sum / stats.count : NaN;
    if (!stats.count) {
      stats.min = NaN;
      stats.max = NaN;
    }
    return stats;
  }

  function updateMetrics() {
    const stats = sequenceRmsxStats();
    const visibleCount = visibleSliceIndexes().length;
    elements.currentSliceMetric.textContent = REPORT ? `${visibleCount}/${REPORT.slices.length}` : "-";
    elements.meanMetric.textContent = formatNumber(stats.mean);
    elements.peakMetric.textContent = `${formatNumber(stats.min)} - ${formatNumber(stats.max)}`;
    elements.peakResidueMetric.textContent = stats.peakResidue;
    elements.residueCountMetric.textContent = String(REPORT?.residues?.length || "-");
    elements.maskedMetric.textContent = `${Number(REPORT.maskSummary?.maskedResidues ?? REPORT.maskSummary?.maskedKeys?.length ?? 0)} / ${Number(REPORT.maskSummary?.totalResidues ?? REPORT.residues?.length ?? 0)}`;
    elements.outlineCheckbox.checked = state.outline;
    elements.viewport.classList.toggle("local-drag-disabled", !state.localDrag);
    elements.thicknessRange.value = String(state.thickness);
    elements.thicknessNumber.value = String(state.thickness);
    elements.spacingRange.value = String(state.spacing);
    elements.spacingNumber.value = String(state.spacing);
    elements.colorMinNumber.value = String(Number(state.colorMin.toFixed(3)));
    elements.colorMaxNumber.value = String(Number(state.colorMax.toFixed(3)));
    elements.radiusMinNumber.value = String(Number(state.radiusMin.toFixed(3)));
    elements.radiusMaxNumber.value = String(Number(state.radiusMax.toFixed(3)));
    elements.rotateSensitivityRange.value = String(Number(state.rotationSensitivity.toFixed(3)));
    elements.rotateSensitivityNumber.value = String(Number(state.rotationSensitivity.toFixed(3)));
    updateLegend();
    updateDiagnostics();
  }

  function hexToRgb(hex) {
    const normalized = String(hex || "").replace("#", "");
    if (!/^[0-9a-fA-F]{6}$/.test(normalized)) {
      return { r: 0, g: 0, b: 0 };
    }
    return {
      r: Number.parseInt(normalized.slice(0, 2), 16),
      g: Number.parseInt(normalized.slice(2, 4), 16),
      b: Number.parseInt(normalized.slice(4, 6), 16)
    };
  }

  function rgbToHex(rgb) {
    return `#${[rgb.r, rgb.g, rgb.b].map((value) => {
      return clamp(Math.round(value), 0, 255).toString(16).padStart(2, "0").toUpperCase();
    }).join("")}`;
  }

  function interpolateHexColor(leftHex, rightHex, fraction) {
    const left = hexToRgb(leftHex);
    const right = hexToRgb(rightHex);
    const t = clamp(fraction, 0, 1);
    return rgbToHex({
      r: left.r + ((right.r - left.r) * t),
      g: left.g + ((right.g - left.g) * t),
      b: left.b + ((right.b - left.b) * t)
    });
  }

  function expectedColorForNormalizedRmsx(normalized) {
    const colors = currentPaletteColors();
    if (!colors.length) {
      return "#000000";
    }
    if (colors.length === 1) {
      return colors[0];
    }
    const scaled = clamp(normalized, 0, 1) * (colors.length - 1);
    const lower = Math.floor(scaled);
    const upper = Math.ceil(scaled);
    if (lower === upper) {
      return colors[lower];
    }
    return interpolateHexColor(colors[lower], colors[upper], scaled - lower);
  }

  function currentPaletteGradient() {
    const colors = currentPaletteColors();
    if (!colors.length) {
      return "#000000, #000000";
    }
    if (colors.length === 1) {
      return `${colors[0]}, ${colors[0]}`;
    }
    return colors.map((color, index) => {
      const pct = colors.length === 1 ? 0 : (index / (colors.length - 1)) * 100;
      return `${color} ${pct.toFixed(2)}%`;
    }).join(", ");
  }

  function mappingLegendStops() {
    const min = colorDomainMin();
    const max = colorDomainMax();
    const mid = min + ((max - min) / 2);
    return [
      { key: "Low", rmsx: min, normalized: 0, radius: wormRadiusMin() },
      { key: "Mid", rmsx: mid, normalized: 0.5, radius: wormRadiusMin() + (wormRadiusSpan() / 2) },
      { key: "High", rmsx: max, normalized: 1, radius: wormRadiusMax() }
    ].map((stop) => ({
      ...stop,
      color: expectedColorForNormalizedRmsx(stop.normalized)
    }));
  }

  function radiusDotSize(radius) {
    return clamp(radius * 5, 7, 24);
  }

  function updateLegend() {
    if (elements.legendColorBar) {
      elements.legendColorBar.style.background = `linear-gradient(90deg, ${currentPaletteGradient()})`;
    }
    const textElements = {
      Low: elements.domainMin,
      Mid: elements.domainMid,
      High: elements.domainMax
    };
    const swatches = {
      Low: elements.legendLowSwatch,
      Mid: elements.legendMidSwatch,
      High: elements.legendHighSwatch
    };
    const radiusDots = {
      Low: elements.legendLowRadius,
      Mid: elements.legendMidRadius,
      High: elements.legendHighRadius
    };
    const radiusLabels = {
      Low: elements.legendLowRadiusLabel,
      Mid: elements.legendMidRadiusLabel,
      High: elements.legendHighRadiusLabel
    };
    mappingLegendStops().forEach((stop) => {
      const value = textElements[stop.key];
      const swatch = swatches[stop.key];
      const dot = radiusDots[stop.key];
      const label = radiusLabels[stop.key];
      if (value) {
        value.textContent = formatNumber(stop.rmsx);
      }
      if (swatch) {
        swatch.style.background = stop.color;
      }
      if (dot) {
        const size = radiusDotSize(stop.radius);
        dot.style.width = `${size.toFixed(1)}px`;
        dot.style.height = `${size.toFixed(1)}px`;
        dot.style.background = stop.color;
      }
      if (label) {
        label.textContent = stop.radius.toFixed(2);
      }
    });
  }

  function legendSummary() {
    return {
      activeRmsxColorDomain: {
        min: Number(colorDomainMin().toFixed(4)),
        max: Number(colorDomainMax().toFixed(4))
      },
      radiusRange: {
        min: Number(wormRadiusMin().toFixed(4)),
        max: Number(wormRadiusMax().toFixed(4))
      },
      stops: mappingLegendStops().map((stop) => ({
        key: stop.key.toLowerCase(),
        rmsx: Number(stop.rmsx.toFixed(4)),
        normalized: stop.normalized,
        radius: Number(stop.radius.toFixed(4)),
        color: stop.color
      })),
      elements: {
        legend: Boolean(document.querySelector('[data-testid="molstar-rmsx-legend"]')),
        colorBar: Boolean(elements.legendColorBar),
        radiusLegend: Boolean(document.querySelector('[data-testid="molstar-radius-legend"]')),
        lowValue: elements.domainMin?.textContent || "",
        midValue: elements.domainMid?.textContent || "",
        highValue: elements.domainMax?.textContent || "",
        lowRadius: elements.legendLowRadiusLabel?.textContent || "",
        midRadius: elements.legendMidRadiusLabel?.textContent || "",
        highRadius: elements.legendHighRadiusLabel?.textContent || ""
      }
    };
  }

  function visibleSliceIndexes() {
    return REPORT.slices.map((_, index) => index).filter((index) => state.visible.has(index));
  }

  function hiddenSliceIndexes() {
    const visible = new Set(visibleSliceIndexes());
    return REPORT.slices.map((_, index) => index).filter((index) => !visible.has(index));
  }

  function firstVisibleSliceIndex() {
    return visibleSliceIndexes()[0] ?? 0;
  }

  function expectedManagedUrlParams() {
    const allVisible = visibleSliceIndexes().length === REPORT.slices.length;
    const residue = selectedResidue();
    return {
      layout: state.layout === defaultLayoutName() ? null : state.layout,
      slice: state.currentIndex === 0 ? null : String(state.currentIndex + 1),
      slices: allVisible ? null : visibleSliceIndexes().map((index) => index + 1).join(","),
      sliceA: null,
      sliceB: null,
      sliceC: null,
      thickness: numbersClose(state.thickness, defaultThickness()) ? null : compactNumber(state.thickness),
      radiusMin: numbersClose(state.radiusMin, defaultRadiusMin()) ? null : compactNumber(state.radiusMin),
      radiusMax: numbersClose(state.radiusMax, defaultRadiusMax()) ? null : compactNumber(state.radiusMax),
      colorMin: numbersClose(colorDomainMin(), defaultColorMin()) ? null : compactNumber(colorDomainMin()),
      colorMax: numbersClose(colorDomainMax(), defaultColorMax()) ? null : compactNumber(colorDomainMax()),
      palette: state.paletteName === defaultPaletteName() ? null : state.paletteName,
      spacing: numbersClose(state.spacing, defaultSpacing()) ? null : compactNumber(state.spacing),
      columns: state.columns === defaultTileColumns() ? null : String(state.columns),
      rotX: numbersClose(state.rotation.x, REPORT.rotationModel?.defaultRotation?.x ?? 90) ? null : compactNumber(state.rotation.x, 0),
      rotY: numbersClose(state.rotation.y, REPORT.rotationModel?.defaultRotation?.y ?? 0) ? null : compactNumber(state.rotation.y, 0),
      rotZ: numbersClose(state.rotation.z, REPORT.rotationModel?.defaultRotation?.z ?? 0) ? null : compactNumber(state.rotation.z, 0),
      residue: residue?.key === defaultResidueKey() ? null : residue?.key,
      marker: state.marker === (REPORT.selectedResidueMarker?.enabledDefault !== false) ? null : state.marker ? "1" : "0",
      delayMs: null,
      localDrag: state.localDrag === true ? null : "0",
      rotateSensitivity: numbersClose(state.rotationSensitivity, 0.35) ? null : compactNumber(state.rotationSensitivity),
      render: state.renderMode === defaultRenderMode() ? null : state.renderMode,
      outline: state.outline === defaultOutline() ? null : state.outline ? "1" : "0",
      cameraMode: state.cameraMode === defaultCameraMode() ? null : state.cameraMode
    };
  }

  function syncUrlState() {
    if (!window.history?.replaceState || !REPORT) {
      return;
    }
    const context = urlStateContext();
    if (!context.writable) {
      return;
    }
    const url = context.url;
    const expected = expectedManagedUrlParams();
    Object.entries(expected).forEach(([name, value]) => {
      if (value === null || value === "") {
        url.searchParams.delete(name);
      } else {
        url.searchParams.set(name, value);
      }
    });
    const current = `${context.targetWindow.location.pathname}${context.targetWindow.location.search}${context.targetWindow.location.hash}`;
    const next = `${url.pathname}${url.search}${url.hash}`;
    if (next !== current) {
      try {
        context.targetWindow.history.replaceState({ flipbookMolstarState: true, scope: context.scope }, "", next);
      } catch (error) {
        console.debug("RMSX Flipbook URL state could not be written; controls will continue without shareable URL sync.", error);
      }
    }
  }

  function urlStateSummary() {
    const context = urlStateContext();
    if (!context.writable) {
      return {
        scope: context.scope,
        writable: false,
        synced: true,
        expectedManagedParams: expectedManagedUrlParams(),
        matches: {},
        note: "Galaxy did not expose a writable iframe or parent visualization URL; viewer state remains live in controls and diagnostics."
      };
    }
    const params = context.url.searchParams;
    const expected = expectedManagedUrlParams();
    const matches = {};
    MANAGED_URL_PARAMS.forEach((name) => {
      matches[name] = expected[name] === null ? !params.has(name) : params.get(name) === expected[name];
    });
    return {
      scope: context.scope,
      writable: true,
      synced: Object.values(matches).every(Boolean),
      expectedManagedParams: expected,
      matches
    };
  }

  function updateDiagnostics() {
    const focusSphere = sceneFocusSphere();
    const diagnostics = {
      schemaVersion: REPORT?.schemaVersion,
      manifestSource: state.manifestSource,
      sliceCount: REPORT?.slices?.length || 0,
      layout: state.layout,
      loadedAllSlices: state.loaded && state.records.length >= (REPORT?.slices?.length || 0),
      liveTransforms: state.liveTransforms,
      coordinateFallback: state.forceCoordinateFallback,
      statusText: elements.status.textContent,
      geometryMode: state.liveTransforms ? "molstar-representation-transform" : state.forceCoordinateFallback ? "browser-side-pdb-copies" : "native-transform-pending",
      stateInitialization: {
        visualizationConfigSources: visualizationConfigStateSummary(),
        urlParamsOverrideVisualizationConfig: true
      },
      mask: {
        maskedResidues: Number(REPORT.maskSummary?.maskedResidues ?? REPORT.maskSummary?.maskedKeys?.length ?? 0),
        totalResidues: Number(REPORT.maskSummary?.totalResidues ?? REPORT.residues?.length ?? 0),
        opacity: Number(REPORT.maskOpacity ?? 0.3)
      },
      palette: state.paletteName,
      visibleSlices: visibleSliceIndexes().map((index) => index + 1),
      focusSphere,
      visibility: {
        loadedSliceIndexes: state.loaded ? REPORT.slices.map((_, index) => index + 1) : [],
        visibleSliceIndexes: visibleSliceIndexes().map((index) => index + 1),
        hiddenSliceIndexes: hiddenSliceIndexes().map((index) => index + 1),
        chipCount: elements.sliceChips.querySelectorAll("[data-testid='molstar-slice-chip']").length
      },
      presentation: {
        layout: state.layout,
        currentIndex: state.currentIndex,
        currentSlice: state.currentIndex + 1,
        spacing: state.spacing,
        columns: state.columns,
        tiledPlacement: tiledPlacementSummary(focusSphere),
        rotation: {
          x: Number(state.rotation.x.toFixed(3)),
          y: Number(state.rotation.y.toFixed(3)),
          z: Number(state.rotation.z.toFixed(3))
        },
        marker: state.marker,
        localDrag: state.localDrag,
        rotationSensitivity: state.rotationSensitivity
      },
      renderStyle: canvasRenderStyleSummary(),
      visualMapping: {
        colorMin: colorDomainMin(),
        colorMax: colorDomainMax(),
        configuredRadiusMin: state.radiusMin,
        configuredRadiusMax: state.radiusMax,
        radiusMin: wormRadiusMin(),
        radiusMax: wormRadiusMax(),
        radiusSpan: wormRadiusSpan(),
        thickness: state.thickness,
        colorPalette: state.paletteName,
        defaultColorPalette: defaultPaletteName(),
        availableColorPalettes: paletteNames(),
        colorTheme: {
          requestedTheme: REPORT.visualMapping?.colorTheme || "uncertainty",
          molstarUncertaintyReversesColorList: REPORT.visualMapping?.molstarUncertaintyReversesColorList === true,
          flipbookLowToHighColors: currentPaletteColors(),
          sentToMolstarUncertaintyColors: [...currentPaletteColors()].reverse().map((hex) => hex.toUpperCase()),
          effectiveOrder: REPORT.visualMapping?.paletteOrder || "low-to-high"
        },
        selectedResidue: selectedResidue()?.key || null,
        selectedRmsx: selectedResidueRmsx(),
        selectedVisualRadius: visualRadiusForRmsx(selectedResidueRmsx()),
        selectedVisualColor: selectedResidueColor(),
        legend: legendSummary()
      },
      controls: {
        sidebarLayout: true,
        compactControlPanels: true,
        compactTabs: false,
        accordionControls: true,
        activeControlPanel: state.activePanel,
        playback: false,
        layoutModes: ["tiled"],
        paletteSwitching: true,
        colorDomain: true,
        radiusRange: true,
        resetScale: true,
        renderStyle: false,
        outlineToggle: true,
        spacing: true,
        thickness: true,
        columns: true,
        rotateSensitivity: true,
        localRotation: true,
        selectedResidueMarker: false,
        maskedResidueStyling: true,
        urlStatePersistence: true
      },
      keyboard: {
        enabled: true,
        rotationStepDegrees: 5,
        spacingStep: spacingStep(),
        thicknessStep: 0.05,
        bindings: ["u/i", "n/m", "j/k", "[/]", "-/=", ",/."],
        shortcutCount: state.keyboardShortcutCount,
        lastAction: state.lastKeyboardAction
      },
      lastRotationGesture: state.lastRotationGesture,
      urlState: urlStateSummary(),
      molstarAssetSource: state.molstarAssetSource
    };
    elements.diagnostics.textContent = JSON.stringify(diagnostics, null, 2);
  }

  function populateControls() {
    const defaultRotation = REPORT.rotationModel?.defaultRotation || { x: 90, y: 0, z: 0 };
    state.layout = initialLayout();
    state.paletteName = availablePalettes()[(URL_PARAMS.get("palette") || "").toLowerCase()]
      ? (URL_PARAMS.get("palette") || "").toLowerCase()
      : defaultPaletteName();
    state.colorMin = numericParam("colorMin", REPORT.domain.min, REPORT.domain.max, defaultColorMin());
    state.colorMax = numericParam("colorMax", REPORT.domain.min, REPORT.domain.max, defaultColorMax());
    state.radiusMin = numericParam("radiusMin", 0.05, 5, defaultRadiusMin());
    state.radiusMax = numericParam("radiusMax", 0.1, 8, defaultRadiusMax());
    state.thickness = numericParam("thickness", 0.25, 2.5, defaultThickness());
    state.spacing = numericParam("spacing", minSpacing(), maxSpacing(), defaultSpacing());
    state.renderMode = renderModeFromParam();
    state.outline = outlineFromParam();
    state.cameraMode = cameraModeFromParam();
    state.rotation = {
      x: numericParam("rotX", -180, 180, Number(defaultRotation.x ?? 90)),
      y: numericParam("rotY", -180, 180, Number(defaultRotation.y ?? 0)),
      z: numericParam("rotZ", -180, 180, Number(defaultRotation.z ?? 0))
    };
    syncRotationMatrixFromEuler();
    state.rotationSensitivity = numericParam("rotateSensitivity", 0.1, 3, 0.35);
    state.columns = integerParam("columns", 1, Math.max(1, REPORT.slices.length), defaultTileColumns());
    state.visible = new Set(initialVisibleSliceIndexes());
    state.currentIndex = initialSliceIndex();
    if (!state.visible.has(state.currentIndex)) {
      state.currentIndex = firstVisibleSliceIndex();
    }
    state.marker = booleanParam("marker", false);
    state.localDrag = booleanParam("localDrag", true);
    state.selectedResidueKey = initialResidueKey();

    elements.columnsNumber.max = String(Math.max(1, REPORT.slices.length));
    elements.columnsNumber.value = String(state.columns);
    elements.spacingRange.min = String(minSpacing());
    elements.spacingNumber.min = String(minSpacing());
    elements.spacingRange.max = String(maxSpacing());
    elements.spacingNumber.max = String(maxSpacing());
    elements.spacingRange.step = String(spacingStep());
    elements.spacingNumber.step = String(spacingStep());
    elements.colorMinNumber.min = String(REPORT.domain.min);
    elements.colorMinNumber.max = String(REPORT.domain.max);
    elements.colorMaxNumber.min = String(REPORT.domain.min);
    elements.colorMaxNumber.max = String(REPORT.domain.max);
    elements.colorMinNumber.step = String(REPORT.visualMapping?.colorDomainStep ?? 0.1);
    elements.colorMaxNumber.step = String(REPORT.visualMapping?.colorDomainStep ?? 0.1);
    elements.radiusMinNumber.step = String(REPORT.visualMapping?.radiusStep ?? 0.05);
    elements.radiusMaxNumber.step = String(REPORT.visualMapping?.radiusStep ?? 0.05);
    elements.paletteSelect.replaceChildren(...paletteNames().map((name) => new Option(name.replace(/[-_]+/g, " "), name)));
    elements.paletteSelect.value = state.paletteName;
    elements.outlineCheckbox.checked = state.outline;
    setActiveControlPanel(state.activePanel);
    updateMetrics();
    renderChips();
    syncUrlState();
  }

  function setActiveControlPanel(panel) {
    const next = CONTROL_PANEL_KEYS.includes(panel) ? panel : "view";
    state.activePanel = next;
    elements.controlPanels.forEach((panelElement) => {
      const active = panelElement.dataset.panel === next;
      panelElement.classList.toggle("active", active);
      if (active) {
        panelElement.open = true;
      }
    });
    updateDiagnostics();
  }

  function renderChips() {
    elements.sliceChips.replaceChildren(...REPORT.slices.map((slice, index) => {
      const button = document.createElement("button");
      button.type = "button";
      button.className = `chip${isSliceVisible(index) ? " active" : ""}`;
      button.dataset.testid = "molstar-slice-chip";
      button.dataset.sliceIndex = String(index + 1);
      button.setAttribute("aria-pressed", state.visible.has(index) ? "true" : "false");
      button.setAttribute("aria-label", state.visible.has(index) ? `Hide item ${index + 1}` : `Show item ${index + 1}`);
      button.title = state.visible.has(index) ? `Hide item ${index + 1}` : `Show item ${index + 1}`;
      button.textContent = String(slice.index ?? index + 1);
      button.addEventListener("click", () => {
        if (state.visible.has(index) && state.visible.size > 1) {
          state.visible.delete(index);
        } else {
          state.visible.add(index);
        }
        state.currentIndex = index;
        renderChips();
        syncUrlState();
        renderScene(true);
      });
      return button;
    }));
  }

  function setLayout(layout) {
    if (!LAYOUTS.has(layout) || state.layout === layout) {
      return;
    }
    state.layout = layout;
    if (!state.visible.has(state.currentIndex)) {
      state.currentIndex = firstVisibleSliceIndex();
    }
    syncUrlState();
    renderChips();
    renderScene(true);
  }

  function reloadScene(autoView = false) {
    state.loaded = false;
    state.liveTransforms = false;
    renderScene(autoView);
  }

  function queueSceneReload(autoView = false, delay = 120) {
    window.clearTimeout(queuedSceneUpdate);
    queuedSceneUpdate = window.setTimeout(() => reloadScene(autoView), delay);
  }

  function queueGeometryUpdate(autoView = false, delay = 60) {
    window.clearTimeout(queuedSceneUpdate);
    queuedSceneUpdate = window.setTimeout(() => {
      if (state.loaded && state.liveTransforms) {
        applyLiveTransforms(autoView, true);
      } else {
        reloadScene(autoView);
      }
    }, delay);
  }

  function queueInteractiveGeometryUpdate(autoView = false) {
    if (interactiveFrame !== null) {
      return;
    }
    const scheduleFrame = window.requestAnimationFrame || ((callback) => window.setTimeout(callback, 16));
    interactiveFrame = scheduleFrame(() => {
      interactiveFrame = null;
      if (state.loaded && state.liveTransforms) {
        applyLiveTransforms(autoView, true);
      } else {
        queueSceneReload(autoView, 60);
      }
    });
  }

  function syncRotationControls() {
    for (const axis of ["x", "y", "z"]) {
      elements[`rotation${axis.toUpperCase()}Range`].value = String(Math.round(state.rotation[axis]));
      elements[`rotation${axis.toUpperCase()}Number`].value = String(Math.round(state.rotation[axis]));
    }
  }

  function addRotation(axis, degrees) {
    state.rotation[axis] = ((Number(state.rotation[axis]) || 0) + degrees + 540) % 360 - 180;
    syncRotationMatrixFromEuler();
    syncRotationControls();
    syncUrlState();
    queueGeometryUpdate(false);
  }

  function updateThickness(value) {
    const next = clamp(Number(value), 0.25, 2.5);
    if (!Number.isFinite(next)) {
      return;
    }
    state.thickness = next;
    elements.thicknessRange.value = next.toFixed(3);
    elements.thicknessNumber.value = next.toFixed(3);
    updateMetrics();
    syncUrlState();
    queueSceneReload(false);
  }

  function updateSpacing(value) {
    const next = clamp(Number(value), minSpacing(), maxSpacing());
    if (!Number.isFinite(next)) {
      return;
    }
    state.spacing = next;
    elements.spacingRange.value = next.toFixed(3);
    elements.spacingNumber.value = next.toFixed(3);
    updateMetrics();
    syncUrlState();
    if (state.layout === "tiled") {
      queueGeometryUpdate(false);
    }
  }

  function updateTileColumns(value) {
    const next = clamp(Math.round(Number(value)), 1, REPORT.slices.length);
    if (!Number.isFinite(next)) {
      return;
    }
    state.columns = next;
    elements.columnsNumber.value = String(next);
    updateMetrics();
    syncUrlState();
    if (state.layout === "tiled") {
      queueGeometryUpdate(false);
    }
  }

  function updatePalette(value) {
    const requested = String(value || "").toLowerCase();
    if (!availablePalettes()[requested]) {
      return;
    }
    state.paletteName = requested;
    elements.paletteSelect.value = requested;
    updateMetrics();
    syncUrlState();
    queueSceneReload(false);
  }

  function updateRenderMode(value) {
    const requested = String(value || "").toLowerCase();
    const next = requested === "clean" ? "clean-interactive" : requested;
    if (!RENDER_PRESETS.has(next)) {
      return;
    }
    state.renderMode = next;
    if (elements.renderSelect) {
      elements.renderSelect.value = next;
    }
    applyMolstarRenderStyle();
    updateMetrics();
    syncUrlState();
  }

  function setOutline(enabled) {
    state.outline = Boolean(enabled);
    elements.outlineCheckbox.checked = state.outline;
    applyMolstarRenderStyle();
    updateMetrics();
    syncUrlState();
  }

  function updateColorDomain(bound, value) {
    const next = clamp(Number(value), REPORT.domain.min, REPORT.domain.max);
    if (!Number.isFinite(next)) {
      return;
    }
    if (bound === "min") {
      state.colorMin = Math.min(next, colorDomainMax() - 0.000001);
    } else {
      state.colorMax = Math.max(next, colorDomainMin() + 0.000001);
    }
    elements.colorMinNumber.value = colorDomainMin().toFixed(3);
    elements.colorMaxNumber.value = colorDomainMax().toFixed(3);
    updateMetrics();
    syncUrlState();
    queueSceneReload(false);
  }

  function updateRadiusRange(bound, value) {
    const next = clamp(Number(value), 0.05, 8);
    if (!Number.isFinite(next)) {
      return;
    }
    if (bound === "min") {
      state.radiusMin = Math.min(next, state.radiusMax - 0.01);
    } else {
      state.radiusMax = Math.max(next, state.radiusMin + 0.01);
    }
    elements.radiusMinNumber.value = state.radiusMin.toFixed(3);
    elements.radiusMaxNumber.value = state.radiusMax.toFixed(3);
    updateMetrics();
    syncUrlState();
    queueSceneReload(false);
  }

  function resetScale() {
    state.colorMin = defaultColorMin();
    state.colorMax = defaultColorMax();
    state.radiusMin = defaultRadiusMin();
    state.radiusMax = defaultRadiusMax();
    state.thickness = defaultThickness();
    updateMetrics();
    syncUrlState();
    queueSceneReload(false);
  }

  function setLocalDrag(enabled) {
    state.localDrag = Boolean(enabled);
    syncUrlState();
    updateMetrics();
  }

  function updateRotateSensitivity(value) {
    const next = clamp(Number(value), 0.1, 3);
    if (!Number.isFinite(next)) {
      return;
    }
    state.rotationSensitivity = next;
    elements.rotateSensitivityRange.value = next.toFixed(3);
    elements.rotateSensitivityNumber.value = next.toFixed(3);
    syncUrlState();
    updateMetrics();
  }

  function setMarker(enabled) {
    state.marker = Boolean(enabled);
    syncUrlState();
    renderScene(false);
  }

  function wireEvents() {
    elements.controlPanels.forEach((panelElement) => {
      panelElement.addEventListener("toggle", () => {
        if (panelElement.open) {
          setActiveControlPanel(panelElement.dataset.panel);
        }
      });
    });
    elements.resetViewButton.addEventListener("click", resetView);
    elements.paletteSelect.addEventListener("change", (event) => updatePalette(event.target.value));
    elements.outlineCheckbox.addEventListener("change", (event) => setOutline(event.target.checked));
    elements.thicknessRange.addEventListener("input", (event) => updateThickness(event.target.value));
    elements.thicknessRange.addEventListener("change", (event) => updateThickness(event.target.value));
    elements.thicknessNumber.addEventListener("input", (event) => updateThickness(event.target.value));
    elements.thicknessNumber.addEventListener("change", (event) => updateThickness(event.target.value));
    elements.spacingRange.addEventListener("input", (event) => updateSpacing(event.target.value));
    elements.spacingRange.addEventListener("change", (event) => updateSpacing(event.target.value));
    elements.spacingNumber.addEventListener("input", (event) => updateSpacing(event.target.value));
    elements.spacingNumber.addEventListener("change", (event) => updateSpacing(event.target.value));
    elements.columnsNumber.addEventListener("input", (event) => updateTileColumns(event.target.value));
    elements.columnsNumber.addEventListener("change", (event) => updateTileColumns(event.target.value));
    for (const axis of ["x", "y", "z"]) {
      const range = elements[`rotation${axis.toUpperCase()}Range`];
      const number = elements[`rotation${axis.toUpperCase()}Number`];
      const handler = (event) => {
        state.rotation[axis] = Number(event.target.value);
        syncRotationMatrixFromEuler();
        syncRotationControls();
        syncUrlState();
        queueGeometryUpdate(false);
      };
      range.addEventListener("input", handler);
      range.addEventListener("change", handler);
      number.addEventListener("input", handler);
      number.addEventListener("change", handler);
    }
    elements.rotateXButton.addEventListener("click", () => addRotation("x", 15));
    elements.rotateYButton.addEventListener("click", () => addRotation("y", 15));
    elements.rotateZButton.addEventListener("click", () => addRotation("z", 15));
    elements.rotateSensitivityRange.addEventListener("input", (event) => updateRotateSensitivity(event.target.value));
    elements.rotateSensitivityRange.addEventListener("change", (event) => updateRotateSensitivity(event.target.value));
    elements.rotateSensitivityNumber.addEventListener("input", (event) => updateRotateSensitivity(event.target.value));
    elements.rotateSensitivityNumber.addEventListener("change", (event) => updateRotateSensitivity(event.target.value));
    elements.resetRotationButton.addEventListener("click", () => {
      state.rotation = { ...(REPORT.rotationModel?.defaultRotation || { x: 90, y: 0, z: 0 }) };
      syncRotationMatrixFromEuler();
      syncRotationControls();
      syncUrlState();
      queueGeometryUpdate(false);
    });
    elements.colorMinNumber.addEventListener("input", (event) => updateColorDomain("min", event.target.value));
    elements.colorMinNumber.addEventListener("change", (event) => updateColorDomain("min", event.target.value));
    elements.colorMaxNumber.addEventListener("input", (event) => updateColorDomain("max", event.target.value));
    elements.colorMaxNumber.addEventListener("change", (event) => updateColorDomain("max", event.target.value));
    elements.radiusMinNumber.addEventListener("input", (event) => updateRadiusRange("min", event.target.value));
    elements.radiusMinNumber.addEventListener("change", (event) => updateRadiusRange("min", event.target.value));
    elements.radiusMaxNumber.addEventListener("input", (event) => updateRadiusRange("max", event.target.value));
    elements.radiusMaxNumber.addEventListener("change", (event) => updateRadiusRange("max", event.target.value));
    elements.resetScaleButton.addEventListener("click", resetScale);
    elements.viewport.addEventListener("pointerdown", (event) => {
      if (!state.localDrag || event.button !== 0) {
        return;
      }
      event.preventDefault();
      event.stopPropagation();
      dragState = { pointerId: event.pointerId, x: event.clientX, y: event.clientY, axes: currentScreenRotationAxes() };
      elements.viewport.classList.add("dragging");
      elements.viewport.setPointerCapture?.(event.pointerId);
    });
    elements.viewport.addEventListener("pointermove", (event) => {
      if (!dragState || dragState.pointerId !== event.pointerId) {
        return;
      }
      event.preventDefault();
      event.stopPropagation();
      const coalesced = typeof event.getCoalescedEvents === "function" ? event.getCoalescedEvents() : [];
      const samples = coalesced.length ? coalesced : [event];
      const latest = samples[samples.length - 1] || event;
      const dx = latest.clientX - dragState.x;
      const dy = latest.clientY - dragState.y;
      const axes = dragState.axes || currentScreenRotationAxes();
      dragState = { pointerId: event.pointerId, x: latest.clientX, y: latest.clientY, axes };
      if (dx === 0 && dy === 0) {
        return;
      }
      applyScreenRotationDrag(dx, dy, axes);
      syncRotationControls();
      queueInteractiveGeometryUpdate(false);
    });
    const endDrag = (event) => {
      if (dragState && dragState.pointerId !== event.pointerId) {
        return;
      }
      dragState = null;
      elements.viewport.classList.remove("dragging");
      elements.viewport.releasePointerCapture?.(event.pointerId);
      syncUrlState();
      queueGeometryUpdate(false, 20);
    };
    elements.viewport.addEventListener("pointerup", endDrag);
    elements.viewport.addEventListener("pointercancel", endDrag);
    document.addEventListener("keydown", (event) => {
      if (event.target && ["INPUT", "SELECT", "TEXTAREA"].includes(event.target.tagName)) {
        return;
      }
      const colorStep = Number(REPORT.visualMapping?.colorDomainStep ?? 0.5);
      const actions = {
        u: ["rotate-x-positive", () => addRotation("x", 5)],
        i: ["rotate-x-negative", () => addRotation("x", -5)],
        n: ["rotate-y-positive", () => addRotation("y", 5)],
        m: ["rotate-y-negative", () => addRotation("y", -5)],
        j: ["rotate-z-positive", () => addRotation("z", 5)],
        k: ["rotate-z-negative", () => addRotation("z", -5)],
        "[": ["thickness-increase", () => updateThickness(state.thickness + 0.05)],
        "]": ["thickness-decrease", () => updateThickness(state.thickness - 0.05)],
        "-": ["spacing-decrease", () => updateSpacing(state.spacing - spacingStep())],
        "=": ["spacing-increase", () => updateSpacing(state.spacing + spacingStep())],
        "+": ["spacing-increase", () => updateSpacing(state.spacing + spacingStep())],
        ",": ["color-domain-low-increase", () => updateColorDomain("min", state.colorMin + colorStep)],
        ".": ["color-domain-high-decrease", () => updateColorDomain("max", state.colorMax - colorStep)]
      };
      const normalizedKey = event.key.length === 1 ? event.key.toLowerCase() : event.key;
      const action = actions[normalizedKey];
      if (action) {
        event.preventDefault();
        event.stopPropagation();
        action[1]();
        state.keyboardShortcutCount += 1;
        state.lastKeyboardAction = action[0];
        updateDiagnostics();
      }
    });
  }

  async function init() {
    try {
      REPORT = await fetchManifest();
      validateManifest(REPORT);
      document.title = REPORT.title || "RMSX Flipbook";
      await loadMolstarAssets();
      populateControls();
      wireEvents();
      syncRotationControls();
      await renderScene();
    } catch (error) {
      setStatus(error.message, true);
      console.error(error);
    }
  }

  init();
}());
