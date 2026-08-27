# --- main.py ---
import os
from pathlib import Path
import subprocess
from install_chimerax import ensure_chimerax
from shiny import run_app

if __name__ == "__main__":
    print("Step 1: Running report.py backend pipeline...")
    try:
        subprocess.run(["python", "report.py"], check=True)
        print("report.py executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error while running report.py: {e}")
        exit(1)

    print("\nStep 2: Check ChimeraX environment...")
    chimerax_path = ensure_chimerax()

    if chimerax_path:
        print(f"ChimeraX ready at: {chimerax_path}")
        os.environ["CHIMERAX_EXE"] = chimerax_path
    else:
        print(
            "WARNING: ChimeraX could not be found or installed. Renders will fail"
            " until installed."
        )

    print("\nStep 3: Starting unified Shiny application...")
    # Pass the app module string filename directly to run_app
    # app.py itself handles the absolute static asset mounting now
    run_app("app.py", port=8000, launch_browser=True)