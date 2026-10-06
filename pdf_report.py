"""PDF report builder: one summary page, then one landscape page per primer."""
import io
from datetime import date
from xml.sax.saxutils import escape

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import (Image, KeepTogether, PageBreak, Paragraph,
                                SimpleDocTemplate, Spacer, Table, TableStyle)

NAVY = colors.HexColor("#1F4E79")
GREY = colors.HexColor("#6B7280")
LIGHT = colors.HexColor("#F3F4F6")
STATUS_COLORS = {
    "OPTIMAL": colors.HexColor("#2E7D32"),
    "SUBOPTIMAL": colors.HexColor("#C62828"),
    "NO DATA": colors.HexColor("#757575"),
}

_ss = getSampleStyleSheet()
H1 = ParagraphStyle("H1", parent=_ss["Title"], fontSize=22, leading=26, textColor=NAVY, alignment=0, spaceAfter=2)
H2 = ParagraphStyle("H2", parent=_ss["Heading2"], fontSize=13, leading=16, textColor=NAVY, spaceBefore=6, spaceAfter=4)
BODY = ParagraphStyle("Body", parent=_ss["BodyText"], fontSize=9, leading=12)
SMALL = ParagraphStyle("Small", parent=BODY, fontSize=8, leading=10, textColor=GREY)
CELL = ParagraphStyle("Cell", parent=BODY, fontSize=8.5, leading=10.5)
CELL_C = ParagraphStyle("CellC", parent=CELL, alignment=TA_CENTER)
HEAD = ParagraphStyle("Head", parent=CELL_C, textColor=colors.white, fontName="Helvetica-Bold")

PAGE = landscape(A4)
MARGIN = 1.4 * cm
CONTENT_W = PAGE[0] - 2 * MARGIN


def _fmt(v, spec):
    return "n/a" if v is None or (isinstance(v, float) and np.isnan(v)) else format(v, spec)


def _dil(v):
    if v is None or pd.isna(v):
        return "Unknown"
    return f"1/{int(v)}" if float(v).is_integer() else f"1/{v:.1f}"


def _png(fig, width):
    """Render a matplotlib figure into a reportlab Image of the given width."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    w_px, h_px = plt.imread(io.BytesIO(buf.getvalue())).shape[1::-1]
    buf.seek(0)
    return Image(buf, width=width, height=width * h_px / w_px)


def _plate_heatmap(merged_df):
    rows, cols = list("ABCDEFGH"), list(range(1, 13))
    grid = merged_df.pivot(index="Row", columns="Col", values="Cq").reindex(index=rows, columns=cols)
    tags = merged_df.pivot(index="Row", columns="Col", values="Tag").reindex(index=rows, columns=cols)
    data = grid.to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(10.5, 3.9))
    cmap = plt.get_cmap("Blues_r").copy()
    cmap.set_bad("#F3F4F6")
    vmin = np.nanmin(data) if np.isfinite(data).any() else 0
    vmax = np.nanmax(data) if np.isfinite(data).any() else 40
    ax.imshow(np.ma.masked_invalid(data), cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    mid = (vmin + vmax) / 2
    for i in range(8):
        for j in range(12):
            v, tag = data[i, j], tags.iloc[i, j]
            if np.isnan(v):
                txt, col = ("-" if tag in ("Empty", None) or pd.isna(tag) else "no Cq"), "#888888"
            else:
                txt, col = f"{v:.1f}", ("white" if v < mid else "#111111")
            ax.text(j, i, txt, ha="center", va="center", fontsize=8, color=col)
    ax.set_xticks(range(12), cols)
    ax.set_yticks(range(8), rows)
    ax.set_xticks(np.arange(-.5, 12, 1), minor=True)
    ax.set_yticks(np.arange(-.5, 8, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="both", length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.xaxis.tick_top()
    return _png(fig, CONTENT_W * 0.86)


def _primer_figure(res, amp_df, melt_df):
    """Amplification, melt peaks and standard curve side by side."""
    active = set(res["std"].loc[res["std"]["Active"], "Well"])
    wells = pd.concat([res["std"][["Well"]].assign(ntc=False), res["ntc"][["Well"]].assign(ntc=True)])

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.3))
    for ax, df, xcol, color, title, xl, yl in (
        (axes[0], amp_df, "Cycle", "#1f77b4", "Amplification", "Cycle", "RFU"),
        (axes[1], melt_df, "Temperature", "#2ca02c", "Melt peaks", "Temperature (°C)", "-d(RFU)/dT"),
    ):
        if df is not None and xcol in df.columns:
            x = df[xcol].values
            for w, is_ntc in zip(wells["Well"], wells["ntc"]):
                if w not in df.columns:
                    continue
                if is_ntc:
                    ax.plot(x, df[w].values, color="#d62728", ls="--", lw=1.1, zorder=3)
                elif w in active:
                    ax.plot(x, df[w].values, color=color, lw=1.1, zorder=2)
                else:
                    ax.plot(x, df[w].values, color="#cccccc", lw=0.8, zorder=1)
        else:
            ax.text(0.5, 0.5, "no data in export", ha="center", va="center", transform=ax.transAxes, color="#888")
        ax.set(title=title, xlabel=xl, ylabel=yl)

    ax = axes[2]
    std = res["std"]
    if not std.empty:
        x_all = -np.log10(std["Dilution_Factor"].values.astype(float))
        ax.scatter(x_all[std["Active"].values], std["Cq"].values[std["Active"].values], s=10, color="#9db7d5", zorder=2)
        ax.scatter(x_all[~std["Active"].values], std["Cq"].values[~std["Active"].values], s=22, marker="x", color="#999999", zorder=2)
    g = res["grp_sel"]
    if len(g) >= 3 and not np.isnan(res["slope"]):
        xs = -np.log10(g["Dilution_Factor"].values.astype(float))
        ax.errorbar(xs, g["mean"].values, yerr=g["std"].fillna(0).values, fmt="o", color="#1F4E79", ms=6, capsize=3, zorder=4)
        xl_ = np.linspace(xs.min() - 0.2, xs.max() + 0.2, 20)
        ax.plot(xl_, res["slope"] * xl_ + res["intercept"], color="#C00000", lw=1.6, zorder=3)
        ax.text(0.04, 0.05, f"Eff {res['eff']:.1f}%\nR² {res['r2']:.4f}\nslope {res['slope']:.3f}",
                transform=ax.transAxes, fontsize=8, va="bottom",
                bbox=dict(boxstyle="round", fc="white", ec="#cccccc"))
    else:
        ax.text(0.5, 0.5, "fewer than 3 dilutions", ha="center", va="center", transform=ax.transAxes, color="#888")
    ax.set(title="Standard curve", xlabel="Log10 relative concentration", ylabel="Cq")

    for a in axes:
        a.grid(alpha=0.25)
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _png(fig, CONTENT_W)


def _dilution_table(res):
    """One row per dilution; excluded wells are struck through."""
    std = res["std"]
    header = [Paragraph(h, HEAD) for h in ("Dilution", "Wells (Cq)", "Used", "Mean Cq", "SD")]
    rows = [header]
    for dil, grp in std.groupby("Dilution_Factor"):
        parts = []
        for _, r in grp.iterrows():
            txt = f"{escape(r['Well'])} ({r['Cq']:.2f})"
            parts.append(txt if r["Active"] else f"<strike>{txt}</strike>")
        act = grp[grp["Active"]]["Cq"]
        rows.append([
            Paragraph(_dil(dil), CELL_C), Paragraph(",&nbsp; ".join(parts), CELL),
            Paragraph(f"{len(act)}/{len(grp)}", CELL_C),
            Paragraph(_fmt(float(act.mean()), ".2f") if len(act) else "excluded", CELL_C),
            Paragraph(_fmt(float(act.std()), ".2f") if len(act) > 1 else "-", CELL_C),
        ])
    if not res["ntc"].empty:
        ntc_txt = ",&nbsp; ".join(
            f"{escape(r['Well'])} ({'no Cq' if pd.isna(r['Cq']) else format(r['Cq'], '.2f')})" for _, r in res["ntc"].iterrows())
        rows.append([Paragraph("NTC", CELL_C), Paragraph(ntc_txt, CELL), Paragraph("", CELL_C),
                     Paragraph(escape(res["ntc_label"]), CELL_C), Paragraph("", CELL_C)])
    t = Table(rows, colWidths=[2.4 * cm, CONTENT_W - 12.4 * cm, 2.2 * cm, 3.2 * cm, 2.2 * cm], repeatRows=1)
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), NAVY),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, LIGHT]),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#D1D5DB")),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
    ]))
    return t


def _metric_strip(res):
    def cell(label, value):
        return [Paragraph(f'<font size="7.5" color="#6B7280">{label}</font>', CELL_C),
                Paragraph(f'<font size="14"><b>{escape(value)}</b></font>', CELL_C)]
    items = [
        ("PCR EFFICIENCY", _fmt(res["eff"], ".1f") + ("%" if not np.isnan(res["eff"]) else "")),
        ("R²", _fmt(res["r2"], ".4f")),
        ("SLOPE", _fmt(res["slope"], ".3f")),
        ("WELLS USED", f"{res['n_active']}/{res['n_total']}"),
        ("NTC", res["ntc_label"]),
    ]
    cols = [cell(l, v) for l, v in items]
    t = Table([[c[0] for c in cols], [c[1] for c in cols]], colWidths=[CONTENT_W / 5] * 5)
    t.setStyle(TableStyle([
        ("BOX", (0, 0), (-1, -1), 0.5, colors.HexColor("#D1D5DB")),
        ("INNERGRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#E5E7EB")),
        ("BACKGROUND", (0, 0), (-1, -1), LIGHT),
        ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    return t


def _primer_page(res, amp_df, melt_df):
    color = STATUS_COLORS[res["status"]]
    head = Table([[Paragraph(f"{escape(res['primer'])} <font color='#6B7280' size='12'>({escape(res['series'])})</font>", H1),
                   Paragraph(f"<font color='white'><b>{res['status']}</b></font>", CELL_C)]],
                 colWidths=[CONTENT_W - 3.6 * cm, 3.6 * cm])
    head.setStyle(TableStyle([("BACKGROUND", (1, 0), (1, 0), color), ("VALIGN", (0, 0), (-1, -1), "MIDDLE")]))

    flow = [head, Spacer(1, 4), _metric_strip(res), Spacer(1, 5)]
    if res["advice"]:
        flow.append(Paragraph(escape(res["advice"]), BODY))
    if res["excluded_wells"]:
        flow.append(Paragraph("Wells excluded from the fit: " + escape(", ".join(res["excluded_wells"])), BODY))
    if res["n_no_cq"]:
        flow.append(Paragraph(f"{res['n_no_cq']} sample well(s) have no Cq (no amplification).", BODY))
    flow += [Spacer(1, 4)]
    flow.append(_primer_figure(res, amp_df, melt_df))
    flow += [Spacer(1, 6), _dilution_table(res)]
    return flow


def _footer(run_name):
    def draw(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 7.5)
        canvas.setFillColor(GREY)
        canvas.drawString(MARGIN, 0.8 * cm, f"qPCR primer validation  |  {run_name}")
        canvas.drawRightString(PAGE[0] - MARGIN, 0.8 * cm, f"Page {doc.page}")
        canvas.restoreState()
    return draw


def build_pdf_report(run_name, layout_name, merged_df, results, amp_df, melt_df):
    """Return the PDF as bytes. `results` is the list of per-primer dicts built in app.py."""
    buf = io.BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=PAGE, leftMargin=MARGIN, rightMargin=MARGIN,
                            topMargin=MARGIN, bottomMargin=1.4 * cm,
                            title=f"qPCR primer validation - {run_name}")

    story = [Paragraph("qPCR primer validation report", H1),
             Paragraph(f"Run: <b>{escape(run_name)}</b> &nbsp;|&nbsp; Layout: {escape(layout_name)} "
                       f"&nbsp;|&nbsp; Generated: {date.today().isoformat()}", SMALL),
             Spacer(1, 8), Paragraph("Summary", H2)]

    cols = ["Primer", "Series", "Raw eff.", "Raw R²", "Eff.", "R²", "Slope", "Wells", "Recommendation", "Status", "NTC"]
    widths = [2.9, 1.7, 1.9, 1.9, 1.9, 1.9, 1.7, 1.5, 5.2, 2.4, 2.6]
    scale = CONTENT_W / (sum(widths) * cm)
    rows = [[Paragraph(c, HEAD) for c in cols]]
    status_cells = []
    for i, r in enumerate(results, start=1):
        rows.append([
            Paragraph(escape(r["primer"]), CELL), Paragraph(escape(r["series"]), CELL_C),
            Paragraph(_fmt(r["eff_raw"], ".1f") + ("%" if not np.isnan(r["eff_raw"]) else ""), CELL_C),
            Paragraph(_fmt(r["r2_raw"], ".4f"), CELL_C),
            Paragraph(_fmt(r["eff"], ".1f") + ("%" if not np.isnan(r["eff"]) else ""), CELL_C),
            Paragraph(_fmt(r["r2"], ".4f"), CELL_C), Paragraph(_fmt(r["slope"], ".3f"), CELL_C),
            Paragraph(f"{r['n_active']}/{r['n_total']}", CELL_C),
            Paragraph(escape(r["recommendation"]), CELL),
            Paragraph(f"<font color='white'><b>{r['status']}</b></font>", CELL_C),
            Paragraph(escape(r["ntc_label"]), CELL_C),
        ])
        status_cells.append(("BACKGROUND", (9, i), (9, i), STATUS_COLORS[r["status"]]))
    summary = Table(rows, colWidths=[w * cm * scale for w in widths], repeatRows=1)
    summary.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), NAVY),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, LIGHT]),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#D1D5DB")),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ] + status_cells))
    story += [summary, Spacer(1, 4),
              Paragraph("Raw = all wells included. Eff. / R² / Slope = the wells ticked in the app at the time of export. "
                        "Optimal = efficiency 90-110% and R² >= 0.98. NTC passes with no Cq or Cq >= 35.", SMALL),
              Spacer(1, 8),
              KeepTogether([Paragraph("96-well plate, Cq values", H2), _plate_heatmap(merged_df)])]

    for res in results:
        story.append(PageBreak())
        story += _primer_page(res, amp_df, melt_df)

    footer = _footer(run_name)
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    return buf.getvalue()
