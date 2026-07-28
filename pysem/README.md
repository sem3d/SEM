# pysem

Tools for SEM3D snapshot parsing and mesh conversion.

## Requirements

- Python 3.11 or later
- An MPI implementation on your system (required to build `mpi4py`):
  - **Linux**: `sudo apt install libopenmpi-dev openmpi-bin` (Debian/Ubuntu) or the equivalent for your distro
  - **macOS**: `brew install open-mpi`
  - **Windows**: [Microsoft MPI](https://learn.microsoft.com/en-us/message-passing-interface/microsoft-mpi) (install both the runtime and SDK)

Below are three ways to install `pysem`: with [uv](https://docs.astral.sh/uv/) (recommended, fastest), with a plain `venv` + `pip`, or with `pip` directly into your current environment.

## Option A — Install with uv

Install `uv` once if you don't have it yet:

```bash
# Linux / macOS
curl -LsSf https://astral.sh/uv/install.sh | sh
```

```powershell
# Windows (PowerShell)
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

Then, from the repository root (same directory as `pyproject.toml`):

```bash
# Linux / macOS / Windows (any shell)
uv venv
uv pip install -e .
```

Or, for a plotting-enabled install:

```bash
uv pip install -e ".[plotting]"
```

`uv venv` creates a `.venv` directory using a Python 3.11+ interpreter (uv will download one automatically if none is available). Activate it with:

```bash
source .venv/bin/activate        # Linux / macOS
.venv\Scripts\activate           # Windows (cmd.exe)
.venv\Scripts\Activate.ps1       # Windows (PowerShell)
```

## Option B — Install with venv + pip

Create and activate a virtual environment with Python's built-in `venv` module, then install with `pip`.

**Linux / macOS**

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -e .
```

**Windows (PowerShell)**

```powershell
py -3.11 -m venv .venv
.venv\Scripts\Activate.ps1
pip install --upgrade pip
pip install -e .
```

**Windows (cmd.exe)**

```bat
py -3.11 -m venv .venv
.venv\Scripts\activate.bat
pip install --upgrade pip
pip install -e .
```

## Option C — Install with pip only

If you already have a Python 3.11+ environment active (system, conda, etc.) and don't want a dedicated virtual environment:

```bash
pip install -e .
```

This works the same way on Linux, macOS, and Windows.

## Optional extras

- `pip install -e ".[plotting]"` — pulls in `pyvista`, `matplotlib`, `ipywidgets` for interactive visualization.
- `pip install -e ".[rikpp]"` — pulls in `dask`, `pyproj`, `seaborn`, `utm`, required by the `pysem.rikpp` fault post-processing tools.

## Command-line tools

Installing the package registers the following commands:

| Command | Entry point |
| --- | --- |
| `parse-snapshots` | `pysem.parse_sem3d_snapshots:main` |
| `parse-traces` | `pysem.parse_sem3d_traces:main` |
| `generate-material` | `pysem.generate_h5_materials:main` |
| `create-stations` | `pysem.create_stations:main` |
| `create-source` | `pysem.sem_stf:main` |
| `create-xmf` | `pysem.create_xmf_from_mesh:main` |
| `create-material-input` | `pysem.create_material_input:main` |
