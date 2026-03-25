
import zipfile
import json
import requests
import io
import streamlit as st

@st.cache_data
def load_database():
    url = "https://raw.githubusercontent.com/Ganeshkumar-byte/GC-MS/main/MoNA-export-GC-MS_Spectra-json.zip"

    response = requests.get(url)
    response.raise_for_status()

    with zipfile.ZipFile(io.BytesIO(response.content)) as z:
        json_filename = [f for f in z.namelist() if f.endswith('.json')][0]

        with z.open(json_filename) as f:
            return json.load(f)

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
    else:
        compounds_to_process = database

    for compound_entry in compounds_to_process:

        if not isinstance(compound_entry, dict):
            continue

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

        score = calculate_match_factor(manual_data, library_peaks)

        results.append({
            "Name": name,
            "Formula": formula,
            "Match Score": score
        })

    df = pd.DataFrame(results)
    return df.sort_values(by="Match Score", ascending=False).head(5)


# ---------------- STREAMLIT UI ---------------- #
st.set_page_config(page_title="GC-MS Matcher", layout="centered")

st.title("🔬 GC-MS Spectral Matcher")

st.write("Enter peaks like: `43:100, 70:12, 61:11`")

# Input
user_input = st.text_area(
    "Enter m/z : intensity values",
    "43:100, 70:12, 61:11, 88:3, 41:8, 42:6"
)

# Load database once
database = load_database()

# Button
if st.button("Find Matches"):

    if not user_input:
        st.warning("Please enter peak values.")
    else:
        with st.spinner("Matching spectra..."):

            query_peaks = parse_user_input(user_input)
            matches = find_top_matches(query_peaks, database)

            st.success("Top Matches Found!")
            st.dataframe(matches)
