import streamlit as st
import pandas as pd

def main():
    st.title("Retention Time Analysis App")
    st.write("Upload your Excel database and enter a retention time to find the closest peaks.")

    # Upload Excel file
    uploaded_file = st.file_uploader("Choose an Excel file", type="xlsx")

    if uploaded_file is not None:
        df = pd.read_excel(uploaded_file)

        # Convert 'Retention Time' column to numeric, coercing errors to NaN
        # Assuming 'Retention Time' and 'Compound Name' are the correct column headers after loading
        if 'Retention Time' in df.columns and 'Compound Name' in df.columns:
            df['Retention Time'] = pd.to_numeric(df['Retention Time'], errors='coerce')
            df.dropna(subset=['Retention Time'], inplace=True)

            # Input retention time using Streamlit's number_input
            input_rt = st.number_input("Enter retention time:", value=0.0, format="%.4f")

            if input_rt > 0:
                # Calculate difference
                df['diff'] = (df['Retention Time'] - input_rt).abs()

                # Get 3 closest peaks
                closest_peaks = df.sort_values(by='diff').head(3)

                # Rename columns for display
                closest_peaks_display = closest_peaks.rename(columns={'Compound Name': 'NAME', 'Retention Time': 'RETENTION TIME'})

                st.subheader("Nearest 3 Peaks:")
                st.dataframe(closest_peaks_display[['NAME', 'RETENTION TIME']])

                # Check for exact match
                exact_match_found = (df['Retention Time'] == input_rt).any()

                if exact_match_found:
                    matched_entries = df[df['Retention Time'] == input_rt]
                    matched_names = matched_entries['Compound Name'].tolist()
                    st.info(f"Statement: The input retention time ({input_rt}) has an exact match for {', '.join(matched_names)}.")
                else:
                    st.warning(f"Statement: No exact retention time match found for {input_rt} in the database. The elution will vary if the specified GC instrument and column conditions are not met, or if no exact match exists.")

                st.markdown("\n--- Notes ---")
                st.write("This database was developed under the following conditions:")
                st.write("Conditions Instrument: GC-2010 Column: SH-200, 60 m, 0.32 mm ID, 1.00 μm (P/N: 227-36186-02) Injection: Split (split ratio: 50:1)Inj. Temp: 250 °C Carrier Gas: He, constant linear velocity mode, 25 cm/sec Oven Temp: 40 °C (0 min) to 310 °C at 4 °C/min Detector: FID, 330 °C")
            else:
                st.info("Please enter a retention time to analyze.")
        else:
            st.error("The uploaded Excel file must contain 'Retention Time' and 'Compound Name' columns.")
    else:
        st.info("Please upload your Excel file to begin.")

if __name__ == '__main__':
    main()
