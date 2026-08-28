import os
from pathlib import Path
from install_chimerax import ensure_chimerax
from shiny import run_app

if __name__ == "__main__":


    print("\nCheck ChimeraX environment...")
    chimerax_path = ensure_chimerax()

    if chimerax_path:
        print(f"ChimeraX ready at: {chimerax_path}")
        os.environ["CHIMERAX_EXE"] = chimerax_path
    else:
        print(
            "WARNING: ChimeraX could not be found or installed. Renders will fail"
            " until installed."
        )

    print("\nStarting unified Shiny application...")
    run_app("app.py", port=8000, launch_browser=True)
