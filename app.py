import io
import os

import pandas as pd
import streamlit as st

# ---------------------- PAGE CONFIG ---------------------- #

st.set_page_config(page_title="BTP DEG Helper", layout="wide")

st.title("BTP DEG Helper: GEO2R + GPL → Up/Down Gene Lists")

st.markdown(
    """
Upload:

1. **GEO2R full results table** (Excel/CSV; columns: `ID`, `adj.P.Val`, `logFC`, …)  
2. **GPL platform table** (Excel/CSV; columns: `ID`, `ORF`, …)

The app will:

- Merge GEO2R results with the GPL annotation on **ID**
- Attach **Locus_tag** from the GPL `ORF` column
- Apply **adj.P.Val** and **logFC** thresholds
- Output **Test-up** and **Control-up** gene lists (SO_XXXX) and a full DEG table
"""
)

# ---------------------- SIDEBAR SETTINGS ---------------------- #

st.sidebar.header("DEG Settings")

fdr_cutoff = st.sidebar.number_input(
    "FDR (adj.P.Val) cutoff",
    min_value=0.0,
    max_value=1.0,
    value=0.05,
    step=0.01,
)

logfc_cutoff = st.sidebar.number_input(
    "Absolute logFC cutoff",
    min_value=0.0,
    max_value=10.0,
    value=1.0,
    step=0.25,
)

add_underscore = st.sidebar.checkbox(
    "Convert SOxxxx → SO_xxxx",
    value=True,
    help="If checked, ORF like 'SO1427' becomes 'SO_1427'.",
)

one_probe_per_gene = st.sidebar.checkbox(
    "Keep one probe per gene (max |logFC|)",
    value=True,
)

st.sidebar.markdown("---")
st.sidebar.markdown("Then upload your files on the main panel.")


# ---------------------- HELPERS ---------------------- #

def read_uploaded_file(uploaded_file: st.runtime.uploaded_file_manager.UploadedFile) -> pd.DataFrame:
    """
    Read an uploaded Excel/CSV/TXT file into a pandas DataFrame,
    using the first row as header and dropping any 'Unnamed:*' index columns.
    """
    if uploaded_file is None:
        return None

    filename = uploaded_file.name.lower()

    if filename.endswith(".xls") or filename.endswith(".xlsx"):
        df = pd.read_excel(uploaded_file, header=0)
    else:
        # autodetect separator for csv/tsv/txt
        df = pd.read_csv(uploaded_file, sep=None, engine="python", header=0)

    # drop auto-generated index columns
    drop_cols = [c for c in df.columns if str(c).startswith("Unnamed")]
    if drop_cols:
        df = df.drop(columns=drop_cols)

    return df


def df_to_tsv_bytes(df: pd.DataFrame) -> bytes:
    buf = io.StringIO()
    df.to_csv(buf, sep="\t", index=False)
    return buf.getvalue().encode()


def series_to_txt_bytes(series: pd.Series) -> bytes:
    buf = io.StringIO()
    series.to_csv(buf, index=False, header=False)
    return buf.getvalue().encode()


# ---------------------- FILE UPLOADS ---------------------- #

st.subheader("Step 1: Upload GEO2R full results table")
geo_file = st.file_uploader(
    "GEO2R results (Excel/CSV with columns: ID, adj.P.Val, logFC, ...)",
    type=["xls", "xlsx", "csv", "txt"],
    key="geo",
)

st.subheader("Step 2: Upload GPL platform table")
gpl_file = st.file_uploader(
    "GPL table (Excel/CSV with columns: ID, ORF, ...)",
    type=["xls", "xlsx", "csv", "txt"],
    key="gpl",
)

if geo_file is None or gpl_file is None:
    st.info("Upload both GEO2R and GPL files to start.")
    st.stop()

# ---------------------- READ FILES ---------------------- #

geo_df = read_uploaded_file(geo_file)
gpl_df = read_uploaded_file(gpl_file)

st.write("### Preview: GEO2R table (first 5 rows)")
st.dataframe(geo_df.head())

st.write("### Preview: GPL table (first 5 rows)")
st.dataframe(gpl_df.head())

# Drop any pre-existing Locus_tag in GEO2R (so we don't clash with new one)
if "Locus_tag" in geo_df.columns:
    geo_df = geo_df.drop(columns=["Locus_tag"])

required_geo_cols = {"ID", "adj.P.Val", "logFC"}
required_gpl_cols = {"ID", "ORF"}

missing_geo = required_geo_cols - set(geo_df.columns)
missing_gpl = required_gpl_cols - set(gpl_df.columns)

if missing_geo:
    st.error(f"GEO2R file is missing columns: {missing_geo}")
    st.stop()

if missing_gpl:
    st.error(f"GPL file is missing columns: {missing_gpl}")
    st.stop()

# Just in case: ensure each probe ID appears once in GEO2R
geo_df = geo_df.drop_duplicates(subset=["ID"])

# ---------------------- MERGE GEO2R + GPL ---------------------- #

st.subheader("Step 3: Merge GEO2R with GPL on ID")

merged = pd.merge(
    geo_df,
    gpl_df[["ID", "ORF"]],
    on="ID",
    how="left",
)

merged = merged.rename(columns={"ORF": "Locus_tag"})

before_drop = len(merged)
merged = merged.dropna(subset=["Locus_tag"])
after_drop = len(merged)

st.write(
    f"Dropped {before_drop - after_drop} rows without Locus_tag. "
    f"Remaining: {after_drop}"
)

# Fix SO IDs if requested
if add_underscore:
    def fix_so(tag: str) -> str:
        tag = str(tag)
        if tag.startswith("SO") and "_" not in tag and len(tag) > 2:
            return "SO_" + tag[2:]
        return tag

    merged["Locus_tag"] = merged["Locus_tag"].apply(fix_so)

# ---------------------- DEG FILTERING ---------------------- #

st.subheader("Step 4: Apply DEG thresholds")

deg_mask = (merged["adj.P.Val"] <= fdr_cutoff) & (merged["logFC"].abs() >= logfc_cutoff)
merged_deg = merged.loc[deg_mask].copy()

if merged_deg.empty:
    st.warning("No DEGs found with current thresholds. Try relaxing cutoffs.")
    st.stop()

# one probe per gene: keep the row with max |logFC|
merged_deg["abs_logFC"] = merged_deg["logFC"].abs()

if one_probe_per_gene:
    merged_deg = merged_deg.sort_values("abs_logFC", ascending=False)
    merged_deg = merged_deg.drop_duplicates(subset=["Locus_tag"])

# Split into test-up and control-up sets
test_up = merged_deg[merged_deg["logFC"] >= logfc_cutoff].copy()
control_up = merged_deg[merged_deg["logFC"] <= -logfc_cutoff].copy()

st.write(f"**Total DEGs after filtering & de-duplicating:** {len(merged_deg)}")
st.write(f"- Test-up (logFC ≥ {logfc_cutoff}): {len(test_up)} genes")
st.write(f"- Control-up (logFC ≤ -{logfc_cutoff}): {len(control_up)} genes")

st.write("### Preview: DEG table (first 10 rows)")
st.dataframe(merged_deg.head(10))

# ---------------------- DOWNLOADS ---------------------- #

base_name = os.path.splitext(os.path.basename(geo_file.name))[0]

deg_bytes = df_to_tsv_bytes(merged_deg)
test_up_bytes = series_to_txt_bytes(test_up["Locus_tag"].drop_duplicates())
control_up_bytes = series_to_txt_bytes(control_up["Locus_tag"].drop_duplicates())

st.subheader("Step 5: Download results")

st.download_button(
    label="Download full DEG table (TSV)",
    data=deg_bytes,
    file_name=f"{base_name}_DEG_table.tsv",
    mime="text/tab-separated-values",
)

st.download_button(
    label="Download TEST-UP locus tags (txt)",
    data=test_up_bytes,
    file_name=f"{base_name}_test_up_locus_tags.txt",
    mime="text/plain",
)

st.download_button(
    label="Download CONTROL-UP locus tags (txt)",
    data=control_up_bytes,
    file_name=f"{base_name}_control_up_locus_tags.txt",
    mime="text/plain",
)

st.markdown(
    """
Now you can:

- Use **TEST-UP** and **CONTROL-UP** locus tag lists directly in **Venny** and **ShinyGO/KEGG**  
- Use the **DEG table** in your BTP report for tables and pathway interpretation
"""
)

