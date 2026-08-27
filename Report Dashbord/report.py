import os
import pandas as pd
from Bio import PDB

# Set to True to load test data, or False to leave paths empty
test = True

# #import app inputs
if test:
    kraken_path = "test data/output_kraken2.txt"
    seqscreen_path = "test data/output_seqscreen.tsv"
    input_dir_ESMFold = "test data/esmfold_output"
else:
    kraken_path = ""
    seqscreen_path = ""
    input_dir_ESMFold = ""


# #from kraken2
if kraken_path:
    df_kraken = pd.read_csv(
        kraken_path,
        sep="\t",
        header=None,
        names=["Status", "Sequence_ID", "Taxonomy_ID", "Length", "Kmer_Mapping"],
    )
else:
    df_kraken = None  # or pd.DataFrame()


# #from seqscreen
if seqscreen_path:
    df_seqscreen = pd.read_csv(seqscreen_path, sep="\t")
else:
    df_seqscreen = None

# Kraken2 output report
total_sequences = len(df_kraken) if df_kraken is not None else 0
unclassified_df = (
    df_kraken[df_kraken["Status"] == "U"]
    if df_kraken is not None
    else pd.DataFrame()
)
classified_df = (
    df_kraken[df_kraken["Status"] != "U"]
    if df_kraken is not None
    else pd.DataFrame()
)

unclassified_count = len(unclassified_df)
classified_count = len(classified_df)
unclassified_pct = (
    (unclassified_count / total_sequences * 100) if total_sequences > 0 else 0
)
classified_pct = (
    (classified_count / total_sequences * 100) if total_sequences > 0 else 0
)

kraken_report_text = f"""Kraken2 Classification Report:
    - Total Sequences Evaluated: {total_sequences}
    - Classified Sequences: {classified_count} ({classified_pct:.2f}%)
    - Unclassified Sequences: {unclassified_count} ({unclassified_pct:.2f}%)"""

# SeqScreen output report
if df_seqscreen is not None:
    de_novo_condition = (
        df_seqscreen["top_hit_protein"].isna()
        | (df_seqscreen["top_hit_protein"] == "")
        | (df_seqscreen["top_hit_protein"].str.lower() == "unclassified")
        | (df_seqscreen["top_hit_protein"].str.lower() == "null")
    )
    de_novo_df = df_seqscreen[de_novo_condition].copy()
    annotated_df = df_seqscreen[~de_novo_condition].copy()
else:
    de_novo_df = pd.DataFrame()
    annotated_df = pd.DataFrame()

total_seqscreen = len(df_seqscreen) if df_seqscreen is not None else 0
de_novo_count = len(de_novo_df)
annotated_count = len(annotated_df)
de_novo_pct = (de_novo_count / total_seqscreen * 100) if total_seqscreen > 0 else 0
annotated_pct = (
    (annotated_count / total_seqscreen * 100) if total_seqscreen > 0 else 0
)

seqscreen_report_text = f"""SeqScreen Functional Characterization Report:
    - Total Input ORFs Evaluated: {total_seqscreen}
    - Functionally Annotated ORFs: {annotated_count} ({annotated_pct:.2f}%)
    - Unclassified Dark Matter ORFs: {de_novo_count} ({de_novo_pct:.2f}%)"""


## ESMFold output report
summary_records = []
plddt_threshold = 70.0

if os.path.exists(input_dir_ESMFold):
    for filename in os.listdir(input_dir_ESMFold):
        if filename.endswith(".pdb"):
            seq_id = filename.replace(".pdb", "")
            pdb_path = os.path.join(input_dir_ESMFold, filename)

            try:
                parser = PDB.PDBParser(QUIET=True)
                structure = parser.get_structure(seq_id, pdb_path)

                plddts = []
                residue_count = 0
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            if "CA" in residue:
                                residue_count += 1
                                val = residue["CA"].get_bfactor()
                                normalized_val = val * 100 if val <= 1.0 else val
                                plddts.append(normalized_val)

                mean_plddt = (sum(plddts) / len(plddts)) if plddts else 0.0
                status = (
                    "Passed" if mean_plddt > plddt_threshold else "Filtered Out"
                )

                summary_records.append(
                    {
                        "Sequence_ID": seq_id,
                        "Length_Residues": residue_count,
                        "Mean_pLDDT": round(mean_plddt, 2),
                        "Status": status,
                        "PDB_Filename": filename,
                    }
                )
            except Exception as e:
                print(f"Error parsing {filename}: {e}")

    esmfold_summary_df = pd.DataFrame(summary_records)
    total_evaluated = len(esmfold_summary_df)
    passed_df = esmfold_summary_df[
        esmfold_summary_df["Status"] == "Passed"
    ]
    passed_count = len(passed_df)
    passed_pct = (
        (passed_count / total_evaluated * 100) if total_evaluated > 0 else 0
    )
    filtered_out_count = total_evaluated - passed_count

    esmfold_report_text = f"""ESMFold Quality Filter Summary (pLDDT > {plddt_threshold}):
    - Total Structures Evaluated: {total_evaluated}
    - Passed Filter: {passed_count} ({passed_pct:.2f}%)
    - Filtered Out: {filtered_out_count}"""
else:
    esmfold_summary_df = pd.DataFrame()
    passed_df = pd.DataFrame()
    esmfold_report_text = "ESMFold Report: Directory not found."

# Final filtering & merging
final_filtered_ids = (
    passed_df["Sequence_ID"].astype(str).tolist()
    if not passed_df.empty
    else []
)

filtered_kraken = (
    unclassified_df[
        unclassified_df["Sequence_ID"].astype(str).isin(final_filtered_ids)
    ].copy()
    if not unclassified_df.empty
    else pd.DataFrame()
)
filtered_seqscreen = (
    de_novo_df[
        de_novo_df["seq_id"].astype(str).isin(final_filtered_ids)
    ].copy()
    if not de_novo_df.empty and "seq_id" in de_novo_df.columns
    else pd.DataFrame()
)
filtered_esmfold = (
    esmfold_summary_df[
        esmfold_summary_df["Sequence_ID"].astype(str).isin(final_filtered_ids)
    ].copy()
    if not esmfold_summary_df.empty
    else pd.DataFrame()
)

if not filtered_kraken.empty and not filtered_seqscreen.empty:
    combined_pipeline_df = pd.merge(
        filtered_kraken,
        filtered_seqscreen,
        left_on="Sequence_ID",
        right_on="seq_id",
        how="inner",
    )
    final_combined_df = pd.merge(
        combined_pipeline_df,
        filtered_esmfold,
        left_on="Sequence_ID",
        right_on="Sequence_ID",
        how="inner",
    )
else:
    final_combined_df = pd.DataFrame()

# Compile text outputs to import into Shiny
full_report_output = f"""{kraken_report_text}

{seqscreen_report_text}

{esmfold_report_text}

Pipeline Merge Complete. Candidates: {len(final_combined_df)}"""