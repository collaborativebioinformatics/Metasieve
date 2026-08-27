# --- install_chimerax.py ---
import os
import platform
from pathlib import Path
import shutil
import subprocess
import winreg

CWD = Path.cwd()
LOCAL_INSTALLER = CWD / "ChimeraX-1.12.exe"
LOCAL_CHIMERAX_DIR = CWD / "ChimeraX"
INSTALL_LOCK = CWD / ".chimerax_installed.lock"


def get_chimerax_from_registry():
  if platform.system() != "Windows":
    return None
  try:
    for root in [winreg.HKEY_LOCAL_MACHINE, winreg.HKEY_CURRENT_USER]:
      for subkey in [
          r"SOFTWARE\UCSF\ChimeraX",
          r"SOFTWARE\WOW6432Node\UCSF\ChimeraX",
      ]:
        try:
          with winreg.OpenKey(root, subkey) as key:
            install_path, _ = winreg.QueryValueEx(key, "")
            exe_path = Path(install_path) / "bin" / "ChimeraX.exe"
            if exe_path.exists():
              return str(exe_path.resolve())
        except FileNotFoundError:
          continue
  except Exception:
    pass
  return None


def ensure_chimerax():
  # 1. Check environment variable override
  env_path = os.environ.get("CHIMERAX_EXE")
  if env_path and Path(env_path).exists():
    return str(Path(env_path).resolve())

  # 2. Check local working directory installation
  local_exe = LOCAL_CHIMERAX_DIR / "bin" / "ChimeraX.exe"
  if local_exe.exists():
    return str(local_exe.resolve())

  # 3. Check Windows Registry
  reg_path = get_chimerax_from_registry()
  if reg_path:
    return reg_path

  # 4. Check standard paths
  system = platform.system()
  if system == "Windows":
    windows_paths = [
        Path(r"C:\Program Files\UCSF ChimeraX\bin\ChimeraX.exe"),
        Path(r"C:\Program Files (x86)\UCSF ChimeraX\bin\ChimeraX.exe"),
        Path.home()
        / r"AppData\Local\Programs\UCSF ChimeraX\bin\ChimeraX.exe",
    ]
    for p in windows_paths:
      if p.exists():
        return str(p.resolve())
  elif system == "Darwin":
    mac_path = Path("/Applications/ChimeraX.app/Contents/MacOS/ChimeraX")
    if mac_path.exists():
      return str(mac_path)
  elif system == "Linux":
    linux_paths = [
        Path("/usr/bin/chimerax"),
        Path("/usr/local/bin/chimerax"),
        Path.home() / ".local/bin/chimerax",
    ]
    for p in linux_paths:
      if p.exists():
        return str(p)

  # 5. Check PATH
  for name in ["chimerax", "ChimeraX"]:
    found_path = shutil.which(name)
    if found_path:
      return found_path

  # 6. Silent install if missing (Windows only)
  if (
      system == "Windows"
      and not local_exe.exists()
      and not INSTALL_LOCK.exists()
      and LOCAL_INSTALLER.exists()
  ):
    print(
        f"ChimeraX not detected. Performing automatic silent installation into"
        f" {LOCAL_CHIMERAX_DIR}..."
    )
    try:
      silent_args = [
          str(LOCAL_INSTALLER),
          "/VERYSILENT",
          "/SUPPRESSMSGBOXES",
          "/NORESTART",
          "/SP-",
          f"/DIR={LOCAL_CHIMERAX_DIR}",
      ]
      subprocess.run(silent_args, check=True)
      INSTALL_LOCK.touch()

      if local_exe.exists():
        return str(local_exe.resolve())
    except Exception as e:
      print(f"Automatic silent installation failed: {e}")

  return None