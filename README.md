# Pandora Target Database

This repository contains tools to index, search, and manage a database of targets observed and their corresponding `InfImg` or `VisSci` FITS files from the Pandora instrument directories.

## Database Schema (`pandora_observations.db`)

The SQLite database contains a `files` table structured as follows:

| Column Name | Data Type | Description |
| :--- | :--- | :--- |
| `id` | `INTEGER` | Primary key (auto-incremented). |
| `filepath` | `TEXT` | Absolute path to the FITS file (unique). |
| `filename` | `TEXT` | Name of the FITS file. |
| `target_id` | `TEXT` | Target name (sourced from FITS `TARG_ID` header, with filename fallback). |
| `file_type` | `TEXT` | Observation image type (`InfImg` or `VisSci`). |
| `instrument` | `TEXT` | Instrument name (`NIRDA` or `VISDA`). |
| `camera_id` | `TEXT` | Camera name (`H2rgCam` or `PcoCam`). |
| `exptime` | `REAL` | Exposure time. |
| `ra` | `REAL` | Target Right Ascension (`TARG_RA`). |
| `dec` | `REAL` | Target Declination (`TARG_DEC`). |
| `obs_time` | `TEXT` | UTC date/time of observation parsed from the filename (`YYYY-MM-DD HH:MM:SS`). |
| `month` | `TEXT` | Observation month subfolder component (e.g., `04`). |
| `day` | `TEXT` | Observation day subfolder component (e.g., `30`). |

---

## CLI Tool Usage

Always execute the script using the virtual environment's Python binary:

```bash
./.pandora_aux/bin/python3 pandora_get_targets.py <command> [args]
```

### 1. Re-index / Update Database

If new files are added to the observation directories, run the `index` command. It uses a cache to skip already indexed files, allowing updates to run in a fraction of a second.

```bash
./.pandora_aux/bin/python3 pandora_get_targets.py index
```

### 2. List All Observed Targets and File Counts

List all unique targets indexed in the database along with the count of `InfImg` and `VisSci` files for each:

```bash
./.pandora_aux/bin/python3 pandora_get_targets.py targets
```

You can search/filter targets containing a specific substring (e.g., to find all `WASP` planets):

```bash
./.pandora_aux/bin/python3 pandora_get_targets.py targets --search WASP
```

### 3. View / Export Files for a Specific Target

To display details and paths of all FITS files matching a target (e.g., `WASP-18b`):

```bash
./.pandora_aux/bin/python3 pandora_get_targets.py files --target WASP-18b
```

To filter by file type (`InfImg` or `VisSci`) and save the absolute file paths to a text file (useful for setting up automated processing pipelines for a specific target):

```bash
./.pandora_aux/bin/python3 pandora_get_targets.py files --target WASP-18b --type InfImg --output wasp18b_inf_files.txt
```

---

## Data Loading Tools (`pandora_tools.py`)

A set of helper functions is provided in pandora_tools.py to programmatically retrieve and load Pandora observation data.

### `read_InfImg(filepath, time_format="JD")`

Reads a single `InfImg` FITS file, extracts the science data cube, and computes chronological frame timestamps.

#### Parameters:
* **`filepath`** *(str)*: Path to the `InfImg` FITS file.
* **`time_format`** *(str)*: Format of the returned timestamps. 
  * `"JD"` (default): Julian Date.
  * `"MJD"`: Modified Julian Date.
  * `"seconds"`: Seconds since midnight J2000 calendar epoch (`2000-01-01T00:00:00Z`).

#### Returns:
* **`ramp_cube`** *(numpy.ndarray)*: 4D array with shape `(nint, ngroup, x, y)` representing the up-the-ramp readout.
* **`timestamps`** *(numpy.ndarray)*: 2D array with shape `(nint, ngroup)` containing chronological frame timestamps.

#### Usage Example:
```python
from pandora_tools import get_target_files, read_InfImg

# 1. Retrieve InfImg file paths for a target
files = get_target_files("WASP-18b", "InfImg")

# 2. Load the datacube and timestamps in MJD format
ramp_cube, mjd_times = read_InfImg(files[0], time_format="MJD")

print("Ramp Cube Shape (nint, ngroup, x, y) :", ramp_cube.shape)
print("Timestamps Shape (nint, ngroup)      :", mjd_times.shape)
```