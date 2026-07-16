# Pandora Observations Analysis Suite

This repository contains tools to index observations, load FITS datacubes, and perform data reduction (ramp fitting and persistence modeling) for the Pandora SmallSat mission.

---

## 1. Pandora Target Database

The database tool indexes Pandora FITS files from the specified observation directories into a SQLite database, allowing users to query observed targets and locate matching files.

### Database Schema (`pandora_observations.db`)

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

### CLI Database Management (`pandora_get_targets.py`)

Always execute the script using the virtual environment's Python binary:

```bash
./.pandora_aux/bin/python3 pandora_get_targets.py <command> [args]
```

#### Index / Update Database
Scans the observation directories and adds new FITS files to the index.
```bash
./.pandora_aux/bin/python3 pandora_get_targets.py index
```

#### List Observed Targets
Lists unique targets with file counts of `InfImg` and `VisSci`.
```bash
./.pandora_aux/bin/python3 pandora_get_targets.py targets
# Search/filter targets
./.pandora_aux/bin/python3 pandora_get_targets.py targets --search WASP
```

#### Locate / Export FITS Files
```bash
./.pandora_aux/bin/python3 pandora_get_targets.py files --target WASP-18b
# Export absolute paths of a specific type to a file:
./.pandora_aux/bin/python3 pandora_get_targets.py files --target WASP-18b --type InfImg --output wasp18b_inf_files.txt
```

### Programmatic Database Access

You can query file paths programmatically using `get_target_files` in `pandora_tools.py`:

```python
from pandora_tools import get_target_files

# Query InfImg FITS files for WASP-18b
files = get_target_files("WASP-18b", "InfImg")
print(f"Found {len(files)} files.")
```

---

## 2. Pandora Data Loading Tools

Helper functions are provided to load FITS files into multidimensional arrays and visualize them.

### `read_InfImg(filepath, time_format="JD")`

Reads a single `InfImg` FITS file, extracts the science data cube, and computes chronological frame timestamps.

#### Parameters:
* **`filepath`** *(str)*: Path to the `InfImg` FITS file.
* **`time_format`** *(str)*: Format of the returned timestamps.
  * `"JD"` (default): Julian Date.
  * `"MJD"`: Modified Julian Date relative to Jan 1, 2026 (`2026-01-01T00:00:00Z`).
  * `"seconds"`: Seconds since midnight J2000 calendar epoch (`2000-01-01T00:00:00Z`).

#### Returns:
* **`ramp_cube`** *(numpy.ndarray)*: 4D array with shape `(nint, ngroup, x, y)` representing the up-the-ramp readout.
* **`timestamps`** *(numpy.ndarray)*: 2D array with shape `(nint, ngroup)` containing chronological frame timestamps.

#### Timestamp Calculation Details:
The spacecraft clock coarse and fine start times are referenced to the J2000 calendar midnight epoch (`2000-01-01 00:00:00 UTC`).
- **Start time ($t_0$)**: $t_0 = \text{CORSTIME} + (\text{FINETIME} \times 10^{-9})$
- **Frame time ($t(i, g)$)**: $t(i, g) = t_0 + (i \times \text{ngroup} + g) \times (\text{FRMTIME} \times 10^{-3})$
- **JD**: $2451544.5 + t(i, g)/86400.0$
- **MJD**: $t(i, g)/86400.0 - 9497.0$ (days relative to Jan 1, 2026).

### `plot_ramp_cube(ramp_cube, integration_index=0, group_index=None, iraf_contrast=0.25, cmap="viridis", output_path=None)`

Plots a single 2D frame from the `ramp_cube` using IRAF-style `zscale` display scaling.

#### Parameters:
* **`ramp_cube`** *(numpy.ndarray)*: 4D array with shape `(nint, ngroup, x, y)`.
* **`integration_index`** *(int)*: Index of the integration to plot (default: `0`).
* **`group_index`** *(int)*: Index of the group to plot. If `None`, defaults to the last group (`ngroup - 1`).
* **`iraf_contrast`** *(float)*: Contrast parameter for `zscale` (default: `0.25`).
* **`cmap`** *(str)*: Colormap (default: `"viridis"`).
* **`output_path`** *(str)*: Path to save the plot.

### Loading and Displaying WASP-18b Data Example:

```python
from pandora_tools import get_target_files, read_InfImg, plot_ramp_cube

# 1. Get file path
files = get_target_files("WASP-18b", "InfImg")

# 2. Load the FITS datacube and timestamps
ramp_cube, jd_times = read_InfImg(files[0], time_format="JD")
nint, ngroup, _, _ = ramp_cube.shape

# 3. Plot the last group of the middle integration using IRAF display scale
plot_ramp_cube(
    ramp_cube=ramp_cube,
    integration_index=nint // 2,
    group_index=ngroup - 1,
    output_path="ramp_plot.png"
)
```

---

## 3. Pandora Data Reduction Tools

Functions to fit slopes for up-the-ramp reads and model/fit detector persistence.

### Up-the-Ramp Fitting: `fit_ramp(timestamps, ramp_cube, read_noise=10.0, gain=1.0)`

Fits an optimal weighted least squares line to each pixel's ramp across groups.

#### Parameters:
* **`timestamps`** *(numpy.ndarray)*: 1D array of shape `(ngroup)` or 2D array of shape `(nint, ngroup)`.
* **`ramp_cube`** *(numpy.ndarray)*: 4D array of shape `(nint, ngroup, x, y)`.
* **`read_noise`** *(float)*: Detector read noise in electrons (default: `10.0`).
* **`gain`** *(float)*: Detector gain in $e^-/\text{DN}$ (default: `1.0`).

#### Returns:
* **`slope_cube`** *(numpy.ndarray)*: 3D array of shape `(nint, x, y)` containing calculated slopes ($\text{DN}/\text{s}$).
* **`intercept_cube`** *(numpy.ndarray)*: 3D array of shape `(nint, x, y)` containing calculated intercepts ($\text{DN}$).

---

### Detector Persistence Modeling

Persistence is tracked using a discrete linear trap charge and release model across the continuous frame timeline.

#### Trapped Charge Model Equations
For each step $k$ in the flattened time series of duration $\Delta t_k = t_k - t_{k-1}$:
1. **Persistence Current Rate ($P_k$)**: $P_k = \frac{Q_{k-1} (1 - e^{-\Delta t_k / \tau})}{\Delta t_k}$ (DN/s)
2. **True Flux ($F_k$)**: $F_k = S_k - P_k$ (DN/s)
3. **Trapped Charge Update ($Q_k$)**: $Q_k = Q_{k-1} e^{-\Delta t_k / \tau} + \epsilon F_k \Delta t_k$ (DN)

#### `calculate_persistence(ramp_cube, timestamps, epsilon=0.18, tau=120.0, Q_init=0.0)`
Computes the persistence model variables ($P$, $F$, $Q$, and $S$) for each pixel.

#### `fit_persistence(ramp_cube, timestamps, epsilon=0.18, tau=120.0)`
Analytically solves for the 2D arrays of true flux $F(x, y)$ and initial trapped charge $Q_{\text{init}}(x, y)$ per pixel.

#### `plot_persistence_model(timestamps, S_cube, F_cube, P_cube, Q_cube, x_pixel, y_pixel, output_path=None)`
Plots model variables over time for a specific pixel.

---

### Data Reduction Example (Ramps and Persistence):

```python
from pandora_tools import get_target_files, read_InfImg, fit_ramp, fit_persistence, calculate_persistence, plot_persistence_model

# 1. Load FITS file
files = get_target_files("WASP-18b", "InfImg")
ramp_cube, times_sec = read_InfImg(files[0], time_format="seconds")

# 2. Fit Ramps (relative times to calculate bias intercept)
times_rel = times_sec - times_sec[:, 0:1]
slope_cube, intercept_cube = fit_ramp(times_rel, ramp_cube)

# 3. Fit Persistence to solve for 2D True Flux and Initial Trapped Charge
F_fit, Q_init_fit = fit_persistence(ramp_cube, times_sec, epsilon=0.18, tau=120.0)

# 4. Generate the full persistence model timeline using the fits
P_cube, F_cube, Q_cube, S_cube = calculate_persistence(
    ramp_cube=ramp_cube,
    timestamps=times_sec,
    epsilon=0.18,
    tau=120.0,
    Q_init=Q_init_fit
)

# 5. Plot the persistence model fit for pixel (40, 125)
plot_persistence_model(
    timestamps=times_sec,
    S_cube=S_cube,
    F_cube=F_cube,
    P_cube=P_cube,
    Q_cube=Q_cube,
    x_pixel=40,
    y_pixel=125,
    output_path="persistence_fit_plot.png"
)
```