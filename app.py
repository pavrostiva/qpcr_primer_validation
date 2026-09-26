import os
import glob
import re
import tempfile
import zipfile
import numpy as np
import pandas as pd
from scipy import stats
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

st.set_page_config(
    page_title="qPCR primer validation",
    page_icon="🧬",
    layout="wide"
)

# ----------------- Helper & Parser Functions -----------------

def normalize_well(val: str) -> str:
    """Normalize well identifiers to standard A01-H12 format."""
    m = re.search(r"([A-Ha-h])\s*0?(\d{1,2})", str(val))
    return f"{m.group(1).upper()}{int(m.group(2)):02d}" if m else str(val).strip()


def parse_dilution_value(val_str: str) -> float:
    """Parse biological dilutions, supporting '1k' (1000) and '10k' (10000)."""
    s = str(val_str).lower().strip()
    if s.endswith("k"):
        return float(s[:-1]) * 1000.0
    return float(s)


def find_files_recursively(base_dir: str):
    """
    Search recursively for Bio-Rad files and user plate layouts.
    Strictly blacklists thermocycler system files from layout matching.
    """
    files_map = {"cq": None, "amp": None, "melt": None, "layout": None}
    layout_candidates = []

    BIORAD_SYSTEM_KEYWORDS = [
        "quantification", "melt curve", "end point", "results_sybr", 
        "run information", "summary", "standard curve results"
    ]

    for root, _, files in os.walk(base_dir):
        for f in files:
            full_path = os.path.join(root, f)
            f_lower = f.lower()

            if "quantification cq results" in f_lower and f.endswith(".csv"):
                files_map["cq"] = full_path
            elif "quantification amplification results" in f_lower and f.endswith(".csv"):
                files_map["amp"] = full_path
            elif "melt curve derivative results" in f_lower and f.endswith(".csv"):
                files_map["melt"] = full_path

            if any(k in f_lower for k in BIORAD_SYSTEM_KEYWORDS):
                continue

            if f_lower.endswith((".xlsx", ".xls", ".csv")):
                score = 0
                name_no_ext = os.path.splitext(f_lower)[0]
                if name_no_ext in ["plate", "layout", "plate_layout", "platemap", "plate_map", "plan", "grid"]:
                    score = 100 if f_lower.endswith((".xlsx", ".xls")) else 80
                elif any(k in name_no_ext for k in ["plate", "layout", "plan", "grid"]):
                    score = 50 if f_lower.endswith((".xlsx", ".xls")) else 30

                if score > 0:
                    layout_candidates.append((score, full_path))

    if layout_candidates:
        layout_candidates.sort(key=lambda x: x[0], reverse=True)
        files_map["layout"] = layout_candidates[0][1]

    return files_map


def parse_plate_file(layout_filepath: str) -> pd.DataFrame:
    """
    Parse 8x12 plate layout grid from Excel or CSV.
    Cleans leading/trailing whitespace and accommodates copy-pasted 10x tags.
    """
    if layout_filepath.endswith((".xlsx", ".xls")):
        df_raw = pd.read_excel(layout_filepath, header=None)
    else:
        df_raw = pd.read_csv(layout_filepath, header=None)

    df_raw = df_raw.fillna("").astype(str)
    records = []
    rows_letters = ["A", "B", "C", "D", "E", "F", "G", "H"]
    COPY_PASTE_10X_MAP = {"1": 1.0, "4": 10.0, "16": 100.0, "64": 1000.0, "256": 10000.0}

    for r_letter in rows_letters:
        row_match = df_raw[df_raw.iloc[:, 0].str.strip().str.upper() == r_letter]
        if row_match.empty and df_raw.shape[1] > 1:
            row_match = df_raw[df_raw.iloc[:, 1].str.strip().str.upper() == r_letter]

        if not row_match.empty:
            row_data = row_match.iloc[0].tolist()
            cells = [c.strip() for c in row_data if c.strip() and c.strip().upper() != r_letter]
            for col_idx, tag in enumerate(cells[:12], start=1):
                well_id = f"{r_letter}{col_idx:02d}"
                parts = [p.strip() for p in tag.split("_") if p.strip()]

                if len(parts) >= 3:
                    primer = parts[0]
                    series = parts[1]
                    d_str = parts[2]
                elif len(parts) == 2:
                    primer = parts[0]
                    series = "std"
                    d_str = parts[1]
                else:
                    primer = parts[0] if parts else tag.strip()
                    series = "std"
                    d_str = "1"

                is_ntc = d_str.upper() == "NTC"
                factor = np.nan

                if not is_ntc:
                    clean_s = series.lower().strip()
                    clean_d = d_str.lower().strip()
                    if "10" in clean_s and clean_d in COPY_PASTE_10X_MAP:
                        factor = COPY_PASTE_10X_MAP[clean_d]
                    else:
                        try:
                            factor = parse_dilution_value(clean_d)
                        except ValueError:
                            factor = np.nan

                records.append({
                    "Well": well_id,
                    "Row": r_letter,
                    "Col": col_idx,
                    "Primer": primer,
                    "Series": series,
                    "Tag": tag.strip(),
                    "Dilution_Factor": factor,
                    "Is_NTC": is_ntc
                })
    return pd.DataFrame(records)


def clean_matrix_file(filepath: str, axis_name: str) -> pd.DataFrame:
    """Clean curve data matrices and drop duplicate columns."""
    if not filepath or not os.path.exists(filepath):
        return None
    df = pd.read_csv(filepath)
    unnamed = [c for c in df.columns if "unnamed" in str(c).lower()]
    if unnamed and len(df.columns) > 1:
        df = df.drop(columns=unnamed)

    new_cols = {}
    for c in df.columns:
        if axis_name.lower() in str(c).lower():
            new_cols[c] = axis_name
        else:
            new_cols[c] = normalize_well(c)

    df = df.rename(columns=new_cols)
    df = df.loc[:, ~df.columns.duplicated(keep="last")]
    return df


def calc_stats(x_vals, y_vals):
    if len(x_vals) < 3:
        return np.nan, np.nan, np.nan, np.nan
    slope, intercept, r_val, p_val, se = stats.linregress(x_vals, y_vals)
    eff = (10 ** (-1.0 / slope) - 1.0) * 100.0
    return slope, intercept, r_val ** 2, eff


# ----------------- Synthetic Demo Generator -----------------

def generate_synthetic_demo_data():
    """Generate artificial synthetic qPCR run for demonstration."""
    np.random.seed(42)
    rows_letters = ["A", "B", "C", "D", "E", "F", "G", "H"]
    records = []
    
    for r in rows_letters:
        for c in range(1, 13):
            well = f"{r}{c:02d}"
            primer = "TargetX" if c <= 6 else "TargetY"
            dil_idx = (c - 1) % 6
            if dil_idx < 5:
                dil_val = 10 ** dil_idx
                tag = f"{primer}_10x_{dil_val}"
                is_ntc = False
                cq = 22.0 + 3.32 * dil_idx + np.random.normal(0, 0.12)
                if primer == "TargetX" and dil_idx == 4:
                    cq = 34.2  # Plateau LOQ demonstration
            else:
                dil_val = np.nan
                tag = f"{primer}_10x_NTC"
                is_ntc = True
                cq = np.nan

            records.append({
                "Well": well, "Row": r, "Col": c, "Primer": primer, "Series": "10x",
                "Tag": tag, "Dilution_Factor": dil_val, "Is_NTC": is_ntc, "Cq": cq
            })

    merged_df = pd.DataFrame(records)

    cycles = np.arange(1, 41)
    amp_dict = {"Cycle": cycles}
    for _, row in merged_df.iterrows():
        w = row["Well"]
        if row["Is_NTC"]:
            amp_dict[w] = np.random.normal(5, 2, 40)
        else:
            cq_val = row["Cq"]
            amp_dict[w] = 3000.0 / (1.0 + np.exp(-(cycles - cq_val) / 1.4)) + np.random.normal(0, 8, 40)
    amp_df = pd.DataFrame(amp_dict)

    temps = np.linspace(60, 95, 71)
    melt_dict = {"Temperature": temps}
    for _, row in merged_df.iterrows():
        w = row["Well"]
        if row["Is_NTC"]:
            melt_dict[w] = np.random.normal(2, 1, len(temps))
        else:
            peak_t = 83.5 if row["Primer"] == "TargetX" else 81.8
            melt_dict[w] = 450.0 * np.exp(-((temps - peak_t) ** 2) / (2 * 1.5 ** 2)) + np.random.normal(0, 4, len(temps))
    melt_df = pd.DataFrame(melt_dict)

    return merged_df, amp_df, melt_df


# ----------------- Sidebar: Controls -----------------

with st.sidebar:
    st.markdown("## 🧬 qPCR primer validation")
    st.caption("Standard curve & MIQE quality control.")

    st.markdown("---")
    use_demo = st.checkbox("🧪 Load demo synthetic data", value=False, help="Try the tool with simulated data without uploading files.")

    active_data_dir = None
    run_name = "qPCR_Run"

    if not use_demo:
        input_mode = st.radio("Input mode:", ["Local folder", "Upload ZIP archive"])
        
        if input_mode == "Upload ZIP archive":
            uploaded_zip = st.file_uploader("Upload .ZIP archive:", type=["zip"])
            if uploaded_zip is not None:
                temp_dir_obj = tempfile.TemporaryDirectory()
                zip_path = os.path.join(temp_dir_obj.name, uploaded_zip.name)
                with open(zip_path, "wb") as f:
                    f.write(uploaded_zip.getbuffer())
                with zipfile.ZipFile(zip_path, "r") as zf:
                    zf.extractall(temp_dir_obj.name)
                active_data_dir = temp_dir_obj.name
                run_name = os.path.splitext(uploaded_zip.name)[0]
        else:
            # Smart Folder Path with detected subfolders
            current_exec_dir = os.getcwd()
            detected_dirs = []
            for root, dirs, _ in os.walk(current_exec_dir):
                for d in dirs:
                    rel_p = os.path.relpath(os.path.join(root, d), current_exec_dir)
                    if not rel_p.startswith(".") and "env" not in rel_p and "git" not in rel_p:
                        detected_dirs.append(rel_p)

            # Quick helper dropdown
            selected_sub = st.selectbox(
                "Quick folder pick (optional):", 
                ["Current directory"] + sorted(detected_dirs)[:15]
            )
            
            initial_path = current_exec_dir if selected_sub == "Current directory" else os.path.join(current_exec_dir, selected_sub)
            local_path = st.text_input("Folder path:", value=initial_path)
            
            if os.path.exists(local_path):
                active_data_dir = local_path
                run_name = os.path.basename(os.path.abspath(local_path))
            else:
                st.error("Folder not found.")

        st.markdown("---")
        st.info("📋 **Expected files:**\n- `*Quantification Cq Results.csv`\n- `*Amplification Results.csv`\n- `*Melt Curve Derivative.csv`\n- `plate.xlsx` (or layout file)")


# ----------------- Data Loading -----------------

if use_demo:
    merged_df, amp_df, melt_df = generate_synthetic_demo_data()
    run_name = "Synthetic_Demo_Experiment"
    layout_name = "Simulated_Plate_96.xlsx"
else:
    if not active_data_dir:
        st.info("👋 Upload a ZIP archive or specify a folder path in the sidebar to begin analysis.")
        st.stop()

    found = find_files_recursively(active_data_dir)

    if not found["cq"]:
        st.error("❌ Could not find `*Quantification Cq Results.csv` in this folder.")
        st.stop()

    if not found["layout"]:
        st.error("❌ Plate layout file not found! Please place an Excel file named `plate.xlsx` or `layout.xlsx` in the folder.")
        st.stop()

    with st.spinner("Processing thermocycler dataset..."):
        cq_df = pd.read_csv(found["cq"])
        cq_df["Well"] = cq_df["Well"].apply(normalize_well)
        cq_df["Cq"] = pd.to_numeric(cq_df["Cq"], errors="coerce")

        amp_df = clean_matrix_file(found["amp"], "Cycle")
        melt_df = clean_matrix_file(found["melt"], "Temperature")
        layout_df = parse_plate_file(found["layout"])

        if layout_df.empty:
            st.error("Failed to parse 8x12 plate layout. Please ensure rows are labeled A through H.")
            st.stop()

        merged_df = pd.merge(layout_df, cq_df[["Well", "Cq"]], on="Well", how="left")
        layout_name = os.path.basename(found["layout"])


# ==============================================================================
#                      SINGLE SCROLLABLE DASHBOARD
# ==============================================================================

st.title(f"🧬 qPCR primer validation: {run_name}")
st.caption(f"Layout source: `{layout_name}` | Rows: 8 (A-H) | Columns: 12 (1-12)")

conditions = merged_df[["Primer", "Series"]].drop_duplicates().sort_values(by=["Primer", "Series"])

# ----------------- SECTION 1: EXECUTIVE SUMMARY -----------------
st.markdown("### 📊 Executive summary & recommendations")

summary_records = []
for _, c_r in conditions.iterrows():
    p, s = c_r["Primer"], c_r["Series"]
    sub_c = merged_df[(merged_df["Primer"] == p) & (merged_df["Series"] == s)]
    
    ntc_val = sub_c[sub_c["Is_NTC"]]["Cq"].mean()
    ntc_stat = "Clean" if pd.isna(ntc_val) else (f"Pass ({ntc_val:.1f})" if ntc_val >= 35 else f"High ({ntc_val:.1f})")

    valid_c = sub_c[~sub_c["Is_NTC"]].dropna(subset=["Cq"])
    g_c = valid_c.groupby("Dilution_Factor")["Cq"].mean().reset_index()

    s_raw, _, r2_raw, eff_raw = calc_stats(-np.log10(g_c["Dilution_Factor"].values), g_c["Cq"].values) if len(g_c) >= 3 else (np.nan, np.nan, np.nan, np.nan)

    # Smart LOQ check
    action = "Keep as-is"
    s_opt, r2_opt, eff_opt = s_raw, r2_raw, eff_raw
    if not np.isnan(eff_raw) and not (90 <= eff_raw <= 110 and r2_raw >= 0.98):
        g_trim = g_c.iloc[:-1]
        if len(g_trim) >= 3:
            s_t, _, r2_t, eff_t = calc_stats(-np.log10(g_trim["Dilution_Factor"].values), g_trim["Cq"].values)
            if 90 <= eff_t <= 110 and r2_t >= 0.98:
                action = f"Exclude 1/{int(g_c.iloc[-1]['Dilution_Factor'])} (plateau/LOQ)"
                s_opt, r2_opt, eff_opt = s_t, r2_t, eff_t

    summary_records.append({
        "Primer": p,
        "Series": s,
        "Raw Eff (%)": f"{eff_raw:.1f}%" if not np.isnan(eff_raw) else "No Amp",
        "Raw R²": f"{r2_raw:.4f}" if not np.isnan(r2_raw) else "N/A",
        "Optimized Eff (%)": f"{eff_opt:.1f}%" if not np.isnan(eff_opt) else "No Amp",
        "Optimized R²": f"{r2_opt:.4f}" if not np.isnan(r2_opt) else "N/A",
        "Smart recommendation": action,
        "Status": "🟢 Optimal" if (90 <= eff_opt <= 110 and r2_opt >= 0.98) else "🔴 Suboptimal",
        "NTC QC": ntc_stat
    })

sum_table = pd.DataFrame(summary_records)
st.dataframe(sum_table, use_container_width=True)

# Export Buttons with Complete Excel Format (both table & 96-well grid)
col_dl1, col_dl2 = st.columns(2)
export_name = f"{run_name}_qPCR"
with col_dl1:
    st.download_button(
        f"📥 Download summary CSV ({export_name}_summary.csv)",
        data=sum_table.to_csv(index=False),
        file_name=f"{export_name}_summary.csv",
        mime="text/csv"
    )
with col_dl2:
    excel_tmp = tempfile.NamedTemporaryFile(suffix=".xlsx", delete=False)
    with pd.ExcelWriter(excel_tmp.name, engine="openpyxl") as writer:
        # Sheet 1: Executive Summary
        sum_table.to_excel(writer, sheet_name="Executive_Summary", index=False)

        # Sheet 2: 96-Well Plate Map (Layout + Cq)
        ws_plate = writer.book.create_sheet(title="Plate_96_Map")
        plate_cq_grid = merged_df.pivot(index="Row", columns="Col", values="Cq")
        plate_tag_grid = merged_df.pivot(index="Row", columns="Col", values="Tag")
        rows_letters = ["A", "B", "C", "D", "E", "F", "G", "H"]

        hdr_fill = PatternFill(start_color="1F497D", end_color="1F497D", fill_type="solid")
        hdr_font = Font(name="Calibri", size=11, bold=True, color="FFFFFF")
        sec_font = Font(name="Calibri", size=12, bold=True, color="1F497D")
        border = Border(left=Side(style="thin", color="D3D3D3"), right=Side(style="thin", color="D3D3D3"),
                        top=Side(style="thin", color="D3D3D3"), bottom=Side(style="thin", color="D3D3D3"))

        # 1. Layout Map
        ws_plate["A1"] = "1. Plate layout map (primers & dilutions):"
        ws_plate["A1"].font = sec_font
        for c in range(1, 13):
            cell = ws_plate.cell(row=3, column=c + 1, value=c)
            cell.fill = hdr_fill; cell.font = hdr_font; cell.alignment = Alignment(horizontal="center")
        for r_idx, r_let in enumerate(rows_letters, start=4):
            ws_plate.cell(row=r_idx, column=1, value=r_let).font = Font(bold=True)
            for c_idx in range(1, 13):
                cell = ws_plate.cell(row=r_idx, column=c_idx + 1, value=plate_tag_grid.loc[r_let, c_idx])
                cell.border = border; cell.alignment = Alignment(horizontal="center")

        # 2. Cq Map
        ws_plate["A15"] = "2. Bio-Rad instrument Cq values:"
        ws_plate["A15"].font = sec_font
        for c in range(1, 13):
            cell = ws_plate.cell(row=17, column=c + 1, value=c)
            cell.fill = hdr_fill; cell.font = hdr_font; cell.alignment = Alignment(horizontal="center")
        for r_idx, r_let in enumerate(rows_letters, start=18):
            ws_plate.cell(row=r_idx, column=1, value=r_let).font = Font(bold=True)
            for c_idx in range(1, 13):
                val = plate_cq_grid.loc[r_let, c_idx]
                cell = ws_plate.cell(row=r_idx, column=c_idx + 1)
                cell.border = border; cell.alignment = Alignment(horizontal="center")
                if pd.isna(val):
                    cell.value = "No Cq"
                    cell.font = Font(color="808080", italic=True)
                else:
                    cell.value = round(val, 2)
                    cell.number_format = "0.00"

        # Sheet 3: Raw Wells
        merged_df[["Well", "Row", "Col", "Primer", "Series", "Tag", "Cq", "Is_NTC"]].to_excel(writer, sheet_name="Raw_Wells", index=False)

        # Auto width
        for ws in [writer.sheets["Executive_Summary"], writer.sheets["Plate_96_Map"], writer.sheets["Raw_Wells"]]:
            for col in ws.columns:
                max_len = max(len(str(c.value or "")) for c in col)
                ws.column_dimensions[get_column_letter(col[0].column)].width = max(max_len + 3, 12)

    with open(excel_tmp.name, "rb") as ef:
        st.download_button(
            f"📊 Download complete Excel ({export_name}_report.xlsx)",
            data=ef.read(),
            file_name=f"{export_name}_report.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )

st.divider()

# ----------------- SECTION 2: 96-WELL HEATMAP WITH CONTRAST NUMBERS -----------------
st.markdown("### 🧪 96-well plate Cq heatmap")
st.caption("Each well displays its exact Cq value. Contrast colors auto-adapt (white text on dark blue, dark text on light blue).")

rows_order = ["H", "G", "F", "E", "D", "C", "B", "A"]
cols_order = list(range(1, 13))

cq_matrix = np.full((8, 12), np.nan)
hover_matrix = np.empty((8, 12), dtype=object)
annotations = []

for r_idx, r_let in enumerate(rows_order):
    for c_idx, c_num in enumerate(cols_order):
        well_id = f"{r_let}{c_num:02d}"
        match = merged_df[merged_df["Well"] == well_id]
        if not match.empty:
            val = match.iloc[0]["Cq"]
            tag = match.iloc[0]["Tag"]
            cq_matrix[r_idx, c_idx] = val
            
            # Smart contrast text annotation: white text on dark blue (val < 27), dark on light
            if pd.isna(val):
                txt = "—"
                font_col = "#888888"
                hover_matrix[r_idx, c_idx] = f"Well: <b>{well_id}</b><br>Sample: {tag}<br><i>No amplification</i>"
            else:
                txt = f"{val:.1f}"
                font_col = "#FFFFFF" if val < 27 else "#111111"
                hover_matrix[r_idx, c_idx] = f"Well: <b>{well_id}</b><br>Sample: {tag}<br>Cq: <b>{val:.2f}</b>"

            annotations.append(dict(
                x=str(c_num),
                y=r_let,
                text=txt,
                showarrow=False,
                font=dict(color=font_col, size=11, family="Arial")
            ))

fig_heat = go.Figure(data=go.Heatmap(
    z=cq_matrix,
    x=[str(c) for c in cols_order],
    y=rows_order,
    hovertext=hover_matrix,
    hoverinfo="text",
    colorscale="Blues_r",
    colorbar=dict(title="Cq"),
    xgap=4,
    ygap=4
))

fig_heat.update_layout(
    annotations=annotations,
    height=420,
    margin=dict(l=40, r=40, t=10, b=40)
)
fig_heat.update_xaxes(tickmode="linear", tick0=1, dtick=1, title="Column (1 - 12)")
fig_heat.update_yaxes(tickmode="linear", dtick=1, title="Row (A - H)")

st.plotly_chart(fig_heat, use_container_width=True)

st.divider()

# ----------------- SECTION 3: DETAILED PRIMER INSPECTOR -----------------
st.markdown("### 🔍 Detailed primer inspector & interactive optimizer")

cond_names = [f"{r['Primer']} ({r['Series']})" for _, r in conditions.iterrows()]
col_p_sel, _ = st.columns([2, 3])
with col_p_sel:
    selected_cond = st.selectbox("Choose primer to inspect:", cond_names)

sel_p, sel_s = selected_cond.split(" ")
sel_s = sel_s.replace("(", "").replace(")", "")
cond_data = merged_df[(merged_df["Primer"] == sel_p) & (merged_df["Series"] == sel_s)].copy()

valid_wells_sub = cond_data[~cond_data["Is_NTC"]].dropna(subset=["Cq"])
grp_pts = valid_wells_sub.groupby("Dilution_Factor")["Cq"].mean().reset_index()

s_init, int_init, r2_init, eff_init = calc_stats(-np.log10(grp_pts["Dilution_Factor"].values), grp_pts["Cq"].values) if len(grp_pts) >= 3 else (np.nan, np.nan, np.nan, np.nan)

# Smart recommendation
adv_msg = "Data fits optimal criteria."
rec_drop = None
if not np.isnan(eff_init) and not (90 <= eff_init <= 110 and r2_init >= 0.98):
    g_trim = grp_pts.iloc[:-1]
    if len(g_trim) >= 3:
        st_i, _, r2t_i, efft_i = calc_stats(-np.log10(g_trim["Dilution_Factor"].values), g_trim["Cq"].values)
        if 90 <= efft_i <= 110 and r2t_i >= 0.98:
            rec_drop = int(grp_pts.iloc[-1]["Dilution_Factor"])
            adv_msg = f"💡 **Smart advice:** Dilution 1/{rec_drop} reached the limit of quantification (plateau). Excluding it yields **Eff: {efft_i:.1f}%** and **R²: {r2t_i:.4f}** (optimal!)."

if rec_drop:
    st.warning(adv_msg)
else:
    st.success("✅ Standard curve satisfies standard guidelines (Eff 90-110%, R² >= 0.98).")

# Interactive Well Manager (OPEN BY DEFAULT)
with st.expander("🛠️ Interactive well manager (toggle replicates)", expanded=True):
    st.write("Uncheck any replicate or dilution point to observe real-time recalculation:")
    chk_cols = st.columns(5)
    active_well_list = []
    for i, (_, w_row) in enumerate(valid_wells_sub.iterrows()):
        col_c = chk_cols[i % 5]
        def_val = True
        if rec_drop and w_row["Dilution_Factor"] == rec_drop:
            def_val = False  # Auto uncheck recommended outlier

        is_checked = col_c.checkbox(
            f"{w_row['Well']} (1/{int(w_row['Dilution_Factor'])}): {w_row['Cq']:.2f}",
            value=def_val,
            key=f"w_chk_{w_row['Well']}"
        )
        if is_checked:
            active_well_list.append(w_row["Well"])

# Live Calculation
filtered_pts = valid_wells_sub[valid_wells_sub["Well"].isin(active_well_list)]
filtered_grp = filtered_pts.groupby("Dilution_Factor")["Cq"].agg(["mean", "std", "count"]).reset_index()

sl_live, int_live, r2_live, eff_live = np.nan, np.nan, np.nan, np.nan
if len(filtered_grp) >= 3:
    act_x = -np.log10(filtered_grp["Dilution_Factor"].values)
    act_y = filtered_grp["mean"].values
    sl_live, int_live, r2_live, eff_live = calc_stats(act_x, act_y)

# KPI Cards with neat delta
k1, k2, k3, k4 = st.columns(4)
eff_delta = f"{eff_live - eff_init:+.1f}% vs baseline" if not np.isnan(eff_init) and not np.isnan(eff_live) else None

k1.metric("PCR efficiency", f"{eff_live:.2f}%" if not np.isnan(eff_live) else "N/A", delta=eff_delta)
k2.metric("R-squared (R²)", f"{r2_live:.4f}" if not np.isnan(r2_live) else "N/A")
k3.metric("Slope", f"{sl_live:.3f}" if not np.isnan(sl_live) else "N/A")
k4.metric("Live status", "🟢 Optimal" if (90 <= eff_live <= 110 and r2_live >= 0.98) else "🔴 Suboptimal")

# 3-Panel Plotly Charts
c_p1, c_p2, c_p3 = st.columns(3)

# 1. Amplification
with c_p1:
    f_amp = go.Figure()
    if amp_df is not None:
        cycles_x = amp_df["Cycle"].values
        for _, w_r in cond_data.iterrows():
            w = w_r["Well"]
            if w in amp_df.columns:
                is_ntc = w_r["Is_NTC"]
                is_active = w in active_well_list
                col = "red" if is_ntc else ("#1f77b4" if is_active else "#dddddd")
                dash = "dash" if is_ntc else "solid"
                name = "NTC" if is_ntc else f"1/{int(w_r['Dilution_Factor'])} ({w})"
                f_amp.add_trace(go.Scatter(x=cycles_x, y=amp_df[w].values, mode="lines", name=name,
                                           line=dict(color=col, dash=dash, width=1.5 if is_active or is_ntc else 0.8)))
    f_amp.update_layout(title="Amplification curves", xaxis_title="Cycle", yaxis_title="RFU", height=380, showlegend=False)
    st.plotly_chart(f_amp, use_container_width=True)

# 2. Melt Peaks
with c_p2:
    f_melt = go.Figure()
    if melt_df is not None:
        temp_x = melt_df["Temperature"].values
        for _, w_r in cond_data.iterrows():
            w = w_r["Well"]
            if w in melt_df.columns:
                is_ntc = w_r["Is_NTC"]
                is_active = w in active_well_list
                col = "red" if is_ntc else ("#2ca02c" if is_active else "#dddddd")
                dash = "dash" if is_ntc else "solid"
                name = "NTC" if is_ntc else f"1/{int(w_r['Dilution_Factor'])} ({w})"
                f_melt.add_trace(go.Scatter(x=temp_x, y=melt_df[w].values, mode="lines", name=name,
                                            line=dict(color=col, dash=dash, width=1.5 if is_active or is_ntc else 0.8)))
    f_melt.update_layout(title="Melt peaks (-d(RFU)/dT)", xaxis_title="Temp (°C)", yaxis_title="-d(RFU)/dT", height=380, showlegend=False)
    st.plotly_chart(f_melt, use_container_width=True)

# 3. Standard Curve
with c_p3:
    f_std = go.Figure()
    if len(filtered_grp) >= 3:
        act_x = -np.log10(filtered_grp["Dilution_Factor"].values)
        act_y = filtered_grp["mean"].values
        act_sd = filtered_grp["std"].fillna(0).values

        f_std.add_trace(go.Scatter(
            x=act_x, y=act_y, mode="markers",
            error_y=dict(type="data", array=act_sd, visible=True),
            marker=dict(size=10, color="#1F4E79"),
            hovertext=[f"Dilution: 1/{int(d)}<br>Mean Cq: {m:.2f} ± {s:.2f}<br>Active wells: {int(c)}" 
                       for d, m, s, c in zip(filtered_grp["Dilution_Factor"], act_y, act_sd, filtered_grp["count"])],
            hoverinfo="text",
            name="Active points"
        ))

        x_line = np.linspace(min(act_x) - 0.2, max(act_x) + 0.2, 50)
        y_line = sl_live * x_line + int_live
        f_std.add_trace(go.Scatter(x=x_line, y=y_line, mode="lines", line=dict(color="#C00000", width=2), name="Fit"))

    f_std.update_layout(title=f"Standard curve (Eff={eff_live:.1f}%, R²={r2_live:.3f})", xaxis_title="Log10 rel conc", yaxis_title="Cq", height=380, showlegend=False)
    st.plotly_chart(f_std, use_container_width=True)