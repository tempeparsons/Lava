# check_install.py
import os, sys, importlib

print(f"Python: {sys.version.split()[0]}")
print(f"Executable: {sys.executable}")

in_conda = os.path.isdir(os.path.join(sys.prefix, "conda-meta"))
in_venv = sys.prefix != sys.base_prefix

if not in_conda and not in_venv:
    print("\nWARNING: you are running the system Python, not a virtual environment.")
    print("If you created an environment for this software, run the script with its")
    print(r"Python instead, e.g.  .\lava_env\Scripts\python.exe check_install.py")
    
print("Environment type:", "conda" if in_conda else "venv" if in_venv else "system Python (no environment)")

required = {
    "numpy": "numpy", "pandas": "pandas", "scipy": "scipy",
    "sklearn": "scikit-learn", "matplotlib": "matplotlib",
    "PySide6": "PySide6",
}
missing = []
for module, package in required.items():
    try:
        m = importlib.import_module(module)
        print(f"OK       {package} {getattr(m, '__version__', '')}")
    except ImportError as e:
        print(f"MISSING  {package}  ({e})")
        missing.append(package)

if missing:
    pkgs = " ".join(missing)
    print("\nTo fix, run this exact command:")
    print(f'  "{sys.executable}" -m pip install {pkgs}')
    if in_conda:
        print("(In a conda env you can also use: conda install -n <env_name> " + pkgs + ")")
    print("\nIf the pip command prints an error, copy it into your bug report.")
    sys.exit(1)

print("\nAll dependencies found.")