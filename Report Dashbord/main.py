# --- main.py ---
import os
from pathlib import Path
from install_chimerax import ensure_chimerax
from shiny import run_app

if __name__ == "__main__":
    print("Step 1: Check ChimeraX environment...")
    chimerax_path = ensure_chimerax()

    if chimerax_path:
        print(f"ChimeraX ready at: {chimerax_path}")
        os.environ["CHIMERAX_EXE"] = chimerax_path
    else:
        print(
            "WARNING: ChimeraX could not be found or installed. Renders will fail"
            " until installed."
        )

    print("\nStep 2: Starting unified Shiny application...")
    # Pass the app module string filename directly to run_app
    # app.py handles pipeline initialization on-demand inside the UI workflow
    run_app("app.py", port=8000, launch_browser=True)