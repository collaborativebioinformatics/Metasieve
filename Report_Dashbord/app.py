# --- app.py ---
import os
from pathlib import Path
import subprocess
from report import final_combined_df, full_report_output
from shiny import App, render, ui

# Retrieve path for ChimeraX if configured via main.py
CHIMERAX_EXE = os.environ.get("CHIMERAX_EXE")

PDB_DIR = Path("test data/esmfold output")

# Resolve absolute path for snapshot directory and static mounting
SNAPSHOT_DIR = Path("Result_Reports").resolve()
SNAPSHOT_DIR.mkdir(parents=True, exist_ok=True)

# Derive dropdown choices STRICTLY from final_combined_df
if not final_combined_df.empty:
    if "PDB_Filename" in final_combined_df.columns:
        pdb_files = final_combined_df["PDB_Filename"].dropna().unique().tolist()
    elif "Sequence_ID" in final_combined_df.columns:
        pdb_files = [f"{str(sid)}.pdb" for sid in final_combined_df["Sequence_ID"].dropna().unique()]
    else:
        pdb_files = []
else:
    pdb_files = []

if not pdb_files:
    pdb_files = ["No passing candidates found"]

# --- SHINY APP UI ---
app_ui = ui.page_fluid(
    ui.h2("Pipeline Result Summary", class_="my-3"),
    ui.hr(),
    # 1. Summary Report Section
    ui.layout_columns(
        ui.card(
            ui.card_header("Execution Summary Reports"),
            ui.output_text_verbatim("analysis_results"),
            full_screen=True,
        ),
        col_widths=12,
    ),
    ui.br(),
    # 2. Interactive Table Section
    ui.layout_columns(
        ui.card(
            ui.card_header("Filtered Candidate Seq"),
            ui.p(
                "Explore, sort, and search through your processed pipeline results:"
            ),
            ui.output_data_frame("combined_table"),
            full_screen=True,
        ),
        col_widths=12,
    ),
    ui.br(),
    # 3. ChimeraX Viewer Section (Displayed right under the table)
    ui.h3("Novel Protein Structure Prediction (ChimeraX )", class_="mt-4"),
    ui.layout_sidebar(
        ui.sidebar(
            ui.h4("Visualization Controls"),
            ui.input_select(
                "selected_pdb",
                "Choose Novel Protein:",
                choices=pdb_files,
                selected=pdb_files[0] if pdb_files and pdb_files[0] != "No passing candidates found" else None,
            ),
            ui.p(
                "Note: Selecting a protein triggers background ChimeraX high-res rendering."
            ),
        ),
        ui.card(
            ui.card_header("ChimeraX High-Accuracy Visual Output"),
            ui.output_ui("chimerax_render"),
        ),
    ),
    class_="p-4",
)


# --- SHINY APP SERVER ---
def server(input, output, session):

    @render.text
    def analysis_results():
        return full_report_output

    @render.data_frame
    def combined_table():
        return render.DataGrid(final_combined_df, filters=True)

    @output
    @render.ui
    def chimerax_render():
        chimerax_exec = os.environ.get("CHIMERAX_EXE") or CHIMERAX_EXE
        
        if not chimerax_exec or not Path(chimerax_exec).exists():
            return ui.p("ChimeraX is not installed or path is missing.", class_="text-danger")

        filename = input.selected_pdb()
        if not filename or filename in ["No passing candidates found", "No PDBs found", ""]:
            return ui.p("No structure selected.")

        pdb_path = PDB_DIR / filename
        if not pdb_path.exists():
            return ui.p(f"PDB file missing: {filename}", class_="text-danger")
            
        output_image_path = SNAPSHOT_DIR / f"{pdb_path.stem}_chimerax.png"

        pdb_file_str = str(pdb_path.resolve()).replace("\\", "/")
        output_img_str = str(output_image_path.resolve()).replace("\\", "/")

        # Always generate / overwrite regardless of whether the file exists
        chimerax_commands = (
            f"open '{pdb_file_str}';"
            "preset ribbon;"
            "color byattribute bfactor palette esmfold;"
            "set bg_color white;"
            "view;"
            "zoom 0.85;"
            "wait 20;"
            f"save '{output_img_str}' width 1920 height 1080;"
            "exit"
        )

        cmd_list = [chimerax_exec, "--cmd", chimerax_commands]

        try:
            print(f"Running ChimeraX render for: {filename} (Overwriting existing snapshot)")
            subprocess.run(
                cmd_list,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            print("ChimeraX render completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"ERROR: ChimeraX failed with exit code {e.returncode}")
            print(f"Stderr: {e.stderr}")
            return ui.p(f"Error rendering {filename} via ChimeraX.", class_="text-danger")
        except Exception as e:
            print(f"ERROR: Unexpected exception during subprocess run: {e}")
            return ui.p(f"Unexpected error: {e}", class_="text-danger")

        if output_image_path.exists() and output_image_path.stat().st_size > 0:
            file_mtime = output_image_path.stat().st_mtime
            img_url = f"/reports/{output_image_path.name}?t={file_mtime}"
            return ui.tags.img(
                src=img_url, 
                style="width: 100%; min-height: 450px; object-fit: contain; border-radius: 4px; background-color: white;"
            )
            
        return ui.p("Waiting for render output...")


# Mount Result_Reports absolute path directly to /reports
app = App(app_ui, server, static_assets={"/reports": SNAPSHOT_DIR})