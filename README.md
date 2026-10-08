<div align="center">

# qPCR Primer Validation

**Check whether your primers work, in a couple of clicks.**
Upload a Bio-Rad CFX run, get PCR efficiency, R², no-template-control status and an interactive plate view.

<br>

<a href="https://qpcr-primer-validation.streamlit.app/">
  <img src="https://img.shields.io/badge/%E2%96%B6%20%20OPEN%20THE%20APP-qpcr--primer--validation.streamlit.app-1F4E79?style=for-the-badge" alt="Open the app" height="56">
</a>

<br><br>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Streamlit](https://img.shields.io/badge/UI-Streamlit-FF4B4B.svg)](https://streamlit.io)

</div>

---

## What is this for?

Before you use a primer pair on real samples, you run a **dilution series** (a *standard curve*) and ask:

- Does the product double every cycle? (**PCR efficiency**, ideally 90–110 %)
- Is the response linear across the dilutions? (**R²**, ideally ≥ 0.98)
- Is the water control clean? (**NTC**)
- Is there one specific product? (**melt curve**, one sharp peak)

The Bio-Rad software exports the raw numbers but does not answer these questions. This app does, and lets you see what happens to the curve when you remove a suspicious well.

It is meant for primer **validation**, not for gene-expression analysis (no ΔΔCt).

---

## How to use it

1. **Export from CFX Maestro** the CSV files of your run (see [Input files](#input-files)).
2. **Fill in a plate layout** in Excel and save it as `plate.xlsx` (see [Plate layout](#plate-layout)).
3. **Put everything in one folder, zip it** and drop the ZIP into the app.

That's it. You get:

| You see | What it tells you |
|---|---|
| **Summary table** | Efficiency, R², NTC status and a green / red verdict for every primer |
| **96-well heatmap** | Cq of every well, so you can spot a missing or odd well at a glance |
| **Well manager** | Every primer is shown on one page; untick a well and its efficiency and R² update instantly |
| **Three plots** | Amplification curves, melt peaks and the standard curve, per primer |
| **Downloads** | A PDF report (one page per primer) and an Excel report with the plate map |

> **Try it first:** tick *Load demo synthetic data* in the sidebar to explore the app without any files.

---

## How to read the result

| Value | Good | What a bad value usually means |
|---|---|---|
| **Efficiency** | 90 – 110 % | < 90 %: inhibition or poorly designed primers. > 110 %: pipetting error, or the most dilute point is beyond the detection limit |
| **R²** | ≥ 0.98 | Noisy replicates or a point that does not fit |
| **Slope** | about −3.32 (= 100 %) | Efficiency is calculated as `10^(−1/slope) − 1` |
| **NTC** | no Cq, or Cq ≥ 35 | Contamination or primer dimers. The app uses the *lowest* NTC Cq, so one dirty replicate is not averaged away |

### Detection-limit advisor

If a curve fails, the app checks whether it would pass without the **most dilute point(s)**. A very dilute sample often drifts to the plateau of the assay's detection limit, which pushes efficiency above 100 %. If dropping that point brings the curve into range, the app says so, shows the improved numbers and unticks those wells for you. You can always tick them back.

The advisor tries removing the one or two most dilute points (at least 3 dilutions always remain, and it warns when fewer than the 5 recommended by MIQE are left). It is a hint, not a verdict: **dropping points must be reported in your methods.**

---

## Input files

Export these from the Bio-Rad CFX software as CSV. The app finds them by name, anywhere inside the folder or ZIP:

| File (name contains) | Required | Used for |
|---|---|---|
| `Quantification Cq Results` | yes | Cq of each well |
| `Quantification Amplification Results` | no | Amplification curves |
| `Melt Curve Derivative Results` | no | Melt peaks |
| `plate.xlsx` | yes | Which primer and dilution is in which well |

A run folder looks like this:

```
my_run/
├── my_run - Quantification Cq Results.csv
├── my_run - Quantification Amplification Results_SYBR.csv
├── my_run - Melt Curve Derivative Results_SYBR.csv
└── plate.xlsx
```

Other files in the folder (Run Information, End Point, etc.) are ignored.

### Plate layout

An Excel sheet with rows **A–H** in the first column and columns **1–12** across. Each cell is one tag:

```
Primer_Series_Dilution
```

| | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| **A** | `H1_10x_1` | `H1_10x_100` | `M1_10x_1` | `M1_10x_100` |
| **B** | `H1_10x_1` | `H1_10x_1k` | `M1_10x_1` | `M1_10x_1k` |
| **C** | `H1_10x_10` | `H1_10x_10k` | `M1_10x_10` | `M1_10x_10k` |
| **D** | `H1_10x_100` | `H1_10x_NTC` | `M1_10x_100` | `M1_10x_NTC` |

- **Primer**: any name; it may contain underscores (`GAPDH_human_10x_100`).
- **Series**: a label for the dilution series, for example `10x` or `4x`. One primer can have several series, and each is evaluated separately.
- **Dilution**: the dilution factor as a number (`1`, `10`, `100`) or with a suffix (`1k` = 1000, `10k` = 10 000).
- **`NTC`** in the dilution position marks a no-template control.
- **Leave unused wells empty.** A partly filled plate is fine, and replicates can be anywhere on the plate; wells are grouped by tag, not by position.
- If a tag cannot be read, the app lists it in a warning instead of silently skipping it.

**A note on copy-pasted labels.** If a series is called `10x` but the dilutions are written `1, 4, 16, 64, 256` (copied from a 4-fold series), the app reads them as `1, 10, 100, 1000, 10000`. Without this the slope would be calculated against the wrong dilution factors (−5.5 instead of −3.3). Please still check your layout.

---

## Run it on your own computer

The hosted app is convenient, but it runs on Streamlit Community Cloud. If your data is confidential, run the app locally; then nothing leaves your machine.

```bash
git clone https://github.com/pavrostiva/qpcr_primer_validation.git
cd qpcr_primer_validation
pip install -r requirements.txt
streamlit run app.py
```

The app opens at `http://localhost:8501`. Locally you can also point it at a folder instead of uploading a ZIP.

---

## Output files

Names start with the name of your ZIP or folder. Both reports follow the wells you have ticked in the well manager.

- **`…_report.pdf`**: click *Build PDF report*, then *Download PDF*. Page 1 has the summary table and the 96-well plate. Then there is one landscape page per primer with efficiency, R², slope, NTC, amplification curves, melt peaks, the standard curve and a table of every well (excluded wells are struck through).
- **`…_report.xlsx`**: `Executive_Summary`, `Plate_96_Map` (your layout and the Cq values as 8×12 grids) and `Raw_Wells` (one row per well, with an `Included_in_fit` column).

If you change a tick after building the PDF, build it again.

---

## Limitations

- Bio-Rad CFX exports only. Other instruments need a different file parser.
- SYBR-style single-channel runs. Multi-channel (e.g. FAM + HEX) runs are not supported.
- Efficiency and R² are calculated on the mean Cq of each dilution, with at least 3 dilutions required. MIQE recommends 5 points over at least 3 orders of magnitude, which is worth keeping in mind when you design the series.
- The melt curve is for you to inspect; the app does not judge it automatically.

## Validation

On the same data, the slopes and efficiencies match the published tool [Auto-qPCR](https://github.com/neuroeddu/Auto-qPCR) to two decimals (for example slope −3.32, efficiency 100.15 % for one primer pair).

## License

MIT
