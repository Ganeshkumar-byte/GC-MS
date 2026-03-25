import streamlit as st
import pandas as pd
import numpy as np
import json
import zipfile
import requests
import io

# ---------------- LOAD DATABASE FROM GITHUB ZIP ---------------- #
@st.cache_data
def load_database():
    url = "https://raw.githubusercontent.com/Ganeshkumar-byte/GC-MS/main/MoNA-export-GC-MS_Spectra-json.zip"

    try:
        response = requests.get(url)
        response.raise_for_status()

        with zipfile.ZipFile(io.BytesIO(response.content)) as z:

            # ✅ Exact path inside your ZIP
            json_path = "MoNA-export-GC-MS_Spectra-json/MoNA-export-GC-MS_Spectra.json"

            if json_path not in z.namelist():
                st.error("❌ JSON file not found inside ZIP")
                st.write("Files inside ZIP:", z.namelist())
                return []

            with z.open(json_path) as f:
                data = json.load(f)

        return data

    except Exception as e:
        st.error(f"Error loading database: {e}")
        return []


# ---------------- MATCH FUNCTION ---------------- #
def calculate_match_factor(query_peaks, library_peaks):
    all_mz = set(query_peaks.keys()).union(set(library_peaks.keys()))

    numerator = 0
    sum_q_sq = 0
    sum_l_sq = 0

    for mz in all_mz:
        q_int = query_peaks.get(mz, 0)
        l_int = library_peaks.get(mz, 0)

        w_q = (mz ** 0.5) * (q_int ** 0.5)
        w_l = (mz ** 0.5) * (l_int ** 0.5)

        numerator += (w_q * w_l)
        sum_q_sq += (w_q ** 2)
        sum_l_sq += (w_l ** 2)

    if sum_q_sq == 0 or sum_l_sq == 0:
        return 0

    score = (numerator**2) / (sum_q_sq * sum_l_sq)
    return round(score * 1000)


# ---------------- PARSE FUNCTIONS ---------------- #
def parse_spectrum_string(spectrum_str):
    peaks_dict = {}

    if not isinstance(spectrum_str, str):
        return peaks_dict

    for peak_pair in spectrum_str.split(' '):
        if ':' in peak_pair:
            try:
                mz_str, intensity_str = peak_pair.split(':')
                mz = int(float(mz_str))
                intensity = float(intensity_str)
                peaks_dict[mz] = intensity
            except:
                continue

    return peaks_dict


def parse_user_input(input_text):
    peaks = {}
    pairs = input_text.split(',')

    for pair in pairs:
        if ':' in pair:
            try:
                mz, intensity = pair.split(':')
                peaks[int(float(mz.strip()))] = float(intensity.strip())
            except:
                continue

    return peaks


# ---------------- MATCH SEARCH ---------------- #
def find_top_matches(manual_data, database):

    results = []

    if isinstance(database, dict):
        compounds_to_process = database.values()
    elif isinstance(database, list):
        compounds_to_process = database
    else:
        st.error(f"Unexpected database type: {type(database)}")
        return pd.DataFrame()

    for compound_entry in compounds_to_process:

        if not isinstance(compound_entry, dict):
            continue

        try:
            # Extract compound details
            if 'compound' in compound_entry and len(compound_entry['compound']) > 0:
                actual = compound_entry['compound'][0]
            else:
                actual = {}

            # Name
            name = 'Unknown'
            if 'names' in actual and len(actual['names']) > 0:
                name = actual['names'][0].get('name', 'Unknown')

            # Formula
            formula = 'N/A'
            if 'metaData' in actual:
                for item in actual['metaData']:
                    if item.get('name') == 'molecular formula':
                        formula = item.get('value', 'N/A')
                        break

            # Spectrum
            spectrum_str = compound_entry.get('spectrum', '')
            library_peaks = parse_spectrum_string(spectrum_str)

            if not library_peaks:
                continue

            score = calculate_match_factor(manual_data, library_peaks)

            results.append({
                "Name": name,
                "Formula": formula,
                "Match Score": score
            })

        except:
            continue

    # ✅ Prevent crash if empty
    if len(results) == 0:
        st.warning("⚠️ No matches found. Check your input or database.")
        return pd.DataFrame(columns=["Name", "Formula", "Match Score"])

    df = pd.DataFrame(results)
    return df.sort_values(by="Match Score", ascending=False).head(5)


# ---------------- STREAMLIT UI ---------------- #
st.set_page_config(page_title="GC-MS Matcher", layout="centered")

st.title("🔬 GC-MS Spectral Matcher")

st.write("Enter peaks like: `43:100, 70:12, 61:11`")

# User input
user_input = st.text_area(
    "Enter m/z : intensity values",
    "43:100, 70:12, 61:11, 88:3, 41:8, 42:6"
)

# Load database
with st.spinner("Loading database..."):
    database = load_database()

# Check database
if not database:
    st.error("❌ Database failed to load. Check ZIP or URL.")
    st.stop()
else:
    st.success(f"✅ Database loaded ({len(database)} entries)")

# Run match
if st.button("Find Matches"):

    if not user_input:
        st.warning("Please enter peak values.")
    else:
        with st.spinner("Matching spectra..."):

            query_peaks = parse_user_input(user_input)

            if not query_peaks:
                st.error("Invalid input format!")
            else:
                matches = find_top_matches(query_peaks, database)

                st.success("Top Matches Found!")
                st.dataframe(matches)
