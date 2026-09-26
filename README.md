# qPCR Primer Validation

> Open-source standard curve engine and MIQE quality control for real-time PCR.  
> If Bio-Rad CFX gives you raw numbers, this tool tells you whether your assay actually works.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Streamlit](https://img.shields.io/badge/UI-Streamlit-FF4B4B.svg)](https://streamlit.io)
[![Plotly](https://img.shields.io/badge/Plots-Plotly-3F4F75.svg)](https://plotly.com)

---

## The 30-Second Workflow

| Step | Action | What Happens Under the Hood |
|---|---|---|
| **01** | **Drop the folder or ZIP** | Recursively pairs Bio-Rad CSV exports with your 8×12 `plate.xlsx` layout. |
| **02** | **Spot the limit (LOQ)** | Detects plateaus where low concentrations hit the noise floor ($C_q$ tailing). |
| **03** | **Export clean reports** | Real-time reactive curves, MIQE validation metrics, and multi-sheet Excel workbooks. |

---

## Problems It Solves

| In the Lab Without This Tool | With This Tool |
|---|---|
| ❌ **Invisible detection limits:** Your dilution curve hits 116% efficiency because the 1:10,000 point flattened out. You spend hours wondering why your slope is wrong. | ✅ **Smart LOQ advisor:** Automatically flags when your lowest dilution hits the limit of quantification and simulates the curve without the plateau. |
| ❌ **Excel formula breakage:** You delete an outlier replicate to see how the curve changes; `=AVERAGE()` throws `#DIV/0!`, infecting `=SLOPE()` and killing the whole sheet. | ✅ **Fault-tolerant calculation:** Native reactive state handles missing wells, zero replicates, or trimmed concentrations without formula errors. |
| ❌ **Copy-paste naming traps:** Lab technicians write `10x_4` or `10x_16` in 10-fold dilution series. Naive scripts calculate slope against 4-fold factors, yielding bizarre –5.51 slopes. | ✅ **Heuristic dilution parser:** Reconciles biological intent. Handles `1k`, `10k`, metric multipliers, and legacy copy-paste labels automatically. |
| ❌ **Blind well lookups:** You have to cross-reference row A05 in three different text files to see if the melt peak matches the amplification curve. | ✅ **Unified 8×12 spatial heatmap:** High-contrast 96-well grid view displays actual $C_q$ values inside the wells alongside sample tags. |

---

## The Four Pillars

```
                     RAW RUN FILES (*.csv) + plate.xlsx
                                    │
                                    ▼
                        [ AUTO-INGEST ENGINE ]
             Regex Well Sanitizer · Heuristic Series Resolver
                                    │
               ┌────────────────────┴────────────────────┐
               ▼                                         ▼
      [ SPATIAL MATRIX ]                        [ LOQ ADVISOR ]
  8×12 Contrast Cq Heatmap                 Limit of Quantification
               │                                         │
               └────────────────────┬────────────────────┘
                                    ▼
                        [ INTERACTIVE DASHBOARD ]
                   Plotly Triple Panel · Live Replicates
                                    │
                                    ▼
                          [ ARTIFACT EXPORT ]
             Summary CSV · Multi-Sheet Formatted Excel (96-Well)
```

1. **Ingest Engine** — Zero-config discovery. Drop a raw run directory or ZIP archive. It maps `A1` to `A01`, strips whitespace, and identifies negative controls (NTC) without manual column mapping.
2. **LOQ Advisor** — Standard curves fail most often because the highest dilution falls below the assay's dynamic range. The advisor calculates both baseline and trimmed models so you see the impact immediately.
3. **Reactive Replicate Filter** — Toggle individual wells or entire dilution tiers with checkboxes to evaluate pipetting artifacts in real time.
4. **Spatial Verification** — View the whole plate at a glance. Text contrast auto-adjusts against background fluorescence intensity so $C_q$ numbers remain readable.

---

## What This Tool Is NOT

- **Not a $\Delta\Delta C_T$ expression differential calculator.** It evaluates primer pairs and assay efficiency *before* you run experimental biological samples.
- **Not a cloud SaaS.** 100% local execution. No biological sequence metadata, primer names, or thermocycler files leave your machine.
- **Not an instrument-locked vendor utility.** Reads open tabular outputs, allowing cross-platform evaluation regardless of software licensing.

---

## Quickstart

### 1. Installation

```bash
git clone https://github.com/pavrostiva/qpcr_primer_validation.git
cd qpcr_primer_validation

pip install streamlit plotly pandas numpy scipy openpyxl
```

### 2. Launch

```bash
streamlit run app.py
```

The browser UI will open at `http://localhost:8501`.

### 3. Try with Synthetic Demo Data
Check **`🧪 Load demo synthetic data`** in the sidebar.  
The tool will synthesize an in-memory 96-well run with amplification kinetics, Gaussian melt peaks, and an intentional LOQ plateau to demonstrate live outlier handling.

---

## Input Expectations

The app looks for files matching standard Bio-Rad CFX export patterns:

```
my_run_folder/
├── my_run - Quantification Cq Results.csv
├── my_run - Quantification Amplification Results_SYBR.csv
├── my_run - Melt Curve Derivative Results_SYBR.csv
└── plate.xlsx
```

### Layout Spreadsheet (`plate.xlsx`)
An 8×12 grid with rows labeled **A–H** and columns **1–12**. Cells use the standard notation:  
`{Primer}_{Series}_{Dilution}`

| | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| **A** | `Actin_4x_1` | `Actin_4x_16` | `GeneX_10x_1` | `GeneX_10x_100` |
| **B** | `Actin_4x_1` | `Actin_4x_64` | `GeneX_10x_1` | `GeneX_10x_1k` |
| ... | ... | ... | ... | ... |
| **H** | `Actin_4x_16` | `Actin_4x_NTC` | `GeneX_10x_100` | `GeneX_10x_NTC` |

- `NTC` tags are routed to negative control evaluation (alerted if $C_q < 35$).
- Metric abbreviations (`1k` $\to$ 1,000, `10k` $\to$ 10,000) parse automatically.

---

## Generated Artifacts

Exports automatically inherit the directory or archive name (e.g., `2026-09-24_qPCR_...`):

- **`*_summary.csv`** — Executive table pairing raw metrics against LOQ-adjusted values (Slope, Efficiency %, $R^2$, NTC pass/fail).
- **`*_report.xlsx`** — Complete Excel package:
  - `Executive_Summary` — Compact summary matrix.
  - `Plate_96_Map` — Visual 8×12 grids containing both user sample tags and instrument $C_q$ values.
  - `Raw_Wells` — Flat auditing table linking every well to its kinetic values.

---

## License

MIT