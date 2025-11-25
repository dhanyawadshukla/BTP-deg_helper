import io
import os

import pandas as pd
import streamlit as st

# ---------------------- PAGE CONFIG ---------------------- #

st.set_page_config(page_title="BTP DEG Helper", layout="wide")

st.title("BTP DEG Helper: GEO2R → Up/Down Gene Lists")

st.markdown(
    """
This app takes **GEO2R results** (and optionally a **GPL platform file**) and gives you
ready-to-use **up- and down-regulated gene lists** (SO_XXXX) plus a DEG table.

It supports two modes:

1. **Single-file mode**: Your GEO2R file already has `ORF` or `Locus_tag`  
   → Just upload that one file.  
2. **Merge mode**: Your GEO2R file has only `ID` (no ORF)  
   → Upload GEO2R + GPL (with `ID` and `ORF`), and the app will merge them.
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

def read_uploaded_file(uploaded_file):
    """Read Excel/CSV/TXT into DataFrame, drop Unnamed index columns."""
    if uploaded_file is None:
        return None

    filename = uploaded_file.name.lower()

    if filename.endswith(".xls") or filename.endswith(".xlsx"):
        df = pd.read_excel(uploaded_file, header=0)
    else:
        df = pd.read_csv(uploaded_file, sep=None, engine="python", header=0)

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


def fix_so(tag: str) -> str:
    tag = str(tag)
    if tag.startswith("SO") and "_" not in tag and len(tag) > 2:
        return "SO_" + tag[2:]
    return tag


# ---------------------- FILE UPLOADS ---------------------- #

st.subheader("Step 1: Upload GEO2R results file (required)")
geo_file = st.file_uploader(
    "GEO2R results (Excel/CSV). Can already contain ORF/Locus_tag, or just ID.",
    type=["xls", "xlsx", "csv", "txt"],
    key="geo",
)

st.subheader("Step 2: Upload GPL platform file (optional)")
st.caption("O
