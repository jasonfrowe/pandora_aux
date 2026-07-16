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

# Data Loading Tools (`pandora_tools.py`)

A set of helper functions is provided in pandora_tools.py to programmatically retrieve and load Pandora observation data.

## `read_InfImg(filepath, time_format="JD")`

Reads a single `InfImg` FITS file, extracts the science data cube, and computes chronological frame timestamps.

### Parameters:
* **`filepath`** *(str)*: Path to the `InfImg` FITS file.
* **`time_format`** *(str)*: Format of the returned timestamps. 
  * `"JD"` (default): Julian Date.
  * `"MJD"`: Modified Julian Date relative to Jan 1, 2026 (`2026-01-01T00:00:00Z`).
  * `"seconds"`: Seconds since midnight J2000 calendar epoch (`2000-01-01T00:00:00Z`).

### Returns:
* **`ramp_cube`** *(numpy.ndarray)*: 4D array with shape `(nint, ngroup, x, y)` representing the up-the-ramp readout.
* **`timestamps`** *(numpy.ndarray)*: 2D array with shape `(nint, ngroup)` containing chronological frame timestamps.

### Timestamp Calculation Details

The spacecraft clock coarse and fine start times are referenced to the J2000 calendar midnight epoch (`2000-01-01 00:00:00 UTC`).

1. **Calculate the start time ($t_0$)** in seconds since epoch:
   $$t_0 = \text{CORSTIME} + (\text{FINETIME} \times 10^{-9})$$

2. **Calculate the elapsed time for each frame** at indices $(i, g)$, where $i \in [0, \text{nint}-1]$ is the integration and $g \in [0, \text{ngroup}-1]$ is the group:
   $$\text{time}_{\text{seconds}}(i, g) = t_0 + (i \times \text{ngroup} + g) \times (\text{FRMTIME} \times 10^{-3})$$
   *(Note: `FRMTIME` in the header is in milliseconds, which is converted to seconds).*

3. **Convert to target format**:
   * **Julian Date (`JD`)**:
     $$\text{JD}(i, g) = 2451544.5 + \frac{\text{time}_{\text{seconds}}(i, g)}{86400.0}$$
     *(The Julian Date of midnight Jan 1, 2000 is  $2451544.5$).*
   * **Modified Julian Date (`MJD`) relative to Jan 1, 2026**:
     $$\text{MJD}(i, g) = \frac{\text{time}_{\text{seconds}}(i, g)}{86400.0} - 9497.0$$
     *(Where $9497.0$ is the number of days between midnight Jan 1, 2000 and midnight Jan 1, 2026).*



## `plot_ramp_cube(ramp_cube, integration_index=0, group_index=None, iraf_contrast=0.25, cmap="viridis", output_path=None)`

Plots a single 2D frame from the 4D `ramp_cube` using IRAF-style `zscale` for display scaling.

### Parameters:
* **`ramp_cube`** *(numpy.ndarray)*: 4D array with shape `(nint, ngroup, x, y)`.
* **`integration_index`** *(int)*: Index of the integration to plot (default: `0`).
* **`group_index`** *(int)*: Index of the group to plot. If `None`, defaults to the last group (`ngroup - 1`).
* **`iraf_contrast`** *(float)*: Contrast parameter for the `zscale` algorithm (default: `0.25`).
* **`cmap`** *(str)*: Colormap name (default: `"viridis"`).
* **`output_path`** *(str)*: Path to save the plot image. If `None`, displays the plot interactively.

---

### Usage Example:

```python
from pandora_tools import get_target_files, read_InfImg, plot_ramp_cube

# 1. Retrieve InfImg file paths for a target
files = get_target_files("WASP-18b", "InfImg")
first_file = files[0]

# 2. Load the datacube and timestamps
ramp_cube, jd_times = read_InfImg(first_file, time_format="JD")

# Get nint and ngroup sizes from the datacube shape directly
nint, ngroup, _, _ = ramp_cube.shape

# 3. Plot the last group of the middle integration
plot_ramp_cube(
    ramp_cube=ramp_cube,
    integration_index=nint // 2,
    group_index=ngroup - 1,
    output_path="ramp_plot.png"
)
```

## `fit_ramp(timestamps, ramp_cube, read_noise=10.0, gain=1.0)`

Fits an optimal weighted least squares line to each pixel's ramp across groups.

### Parameters:
* **`timestamps`** *(numpy.ndarray)*: 1D array of shape `(ngroup)` or 2D array of shape `(nint, ngroup)` representing the time coordinate of each read (in seconds).
* **`ramp_cube`** *(numpy.ndarray)*: 4D array of shape `(nint, ngroup, x, y)`.
* **`read_noise`** *(float)*: Detector read noise in electrons (default: `10.0`).
* **`gain`** *(float)*: Detector gain in $e^-/\text{DN}$ (default: `1.0`).

### Returns:
* **`slope_cube`** *(numpy.ndarray)*: 3D array of shape `(nint, x, y)` containing calculated slopes ($\text{DN}/\text{s}$).
* **`intercept_cube`** *(numpy.ndarray)*: 3D array of shape `(nint, x, y)` containing calculated intercepts ($\text{DN}$).

---

### Usage Example:

```python
from pandora_tools import get_target_files, read_InfImg, fit_ramp

# 1. Retrieve FITS file and load data
files = get_target_files("WASP-18b", "InfImg")
ramp_cube, times_sec = read_InfImg(files[0], time_format="seconds")

# 2. Convert absolute seconds to relative seconds from the start of each integration
# (This keeps the fit intercepts physically equal to the starting bias levels in DN)
times_rel = times_sec - times_sec[:, 0:1]

# 3. Compute optimal least squares fits in parallel
slope_cube, intercept_cube = fit_ramp(times_rel, ramp_cube, read_noise=10.0, gain=1.0)

print("Slope Cube shape:", slope_cube.shape)         # (89, 80, 250)
print("Intercept Cube shape:", intercept_cube.shape) # (89, 80, 250)
```

---

# Detector Persistence Model

A discrete linear trap charge and release model is implemented to track detector persistence on a frame-by-frame basis across the flattened time-series (covering both read ramps and resets).

## Trapped Charge Model Equations

For each step $k$ in the flattened time series of duration $\Delta t_k = t_k - t_{k-1}$:

1. **Persistence Current Rate ($P_k$)**:
   $$P_k = \frac{Q_{k-1} (1 - e^{-\Delta t_k / \tau})}{\Delta t_k}$$
   *(The rate in DN/s of trapped charge released during the step).*

2. **True Flux ($F_k$)**:
   $$F_k = S_k - P_k$$
   *(Where $S_k$ is the observed signal rate, estimated at each frame by taking differences, taking care at the reset frames).*

3. **Trapped Charge Update ($Q_k$)**:
   $$Q_k = Q_{k-1} e^{-\Delta t_k / \tau} + \epsilon F_k \Delta t_k$$
   *(The trapped charge at the end of the step).*

---

## `calculate_persistence(ramp_cube, timestamps, epsilon=0.18, tau=120.0, Q_init=0.0)`

Computes the trapped-charge persistence model for each pixel.

### Parameters:
* **`ramp_cube`** *(numpy.ndarray)*: 4D array of shape `(nint, ngroup, x, y)`.
* **`timestamps`** *(numpy.ndarray)*: 2D array of shape `(nint, ngroup)` containing absolute or relative timestamps (in seconds).
* **`epsilon`** *(float or numpy.ndarray)*: Trapping efficiency. Can be a scalar or a 2D array of shape `(x, y)` (default: `0.18`).
* **`tau`** *(float or numpy.ndarray)*: Trapping decay constant in seconds. Can be a scalar or a 2D array of shape `(x, y)` (default: `120.0`).
* **`Q_init`** *(float or numpy.ndarray)*: Initial trapped charge. Can be a scalar or a 2D array of shape `(x, y)` (default: `0.0`).

### Returns:
* **`P_cube`** *(numpy.ndarray)*: 3D array of shape `(nint * ngroup, x, y)` containing persistence rates ($\text{DN}/\text{s}$).
* **`F_cube`** *(numpy.ndarray)*: 3D array of shape `(nint * ngroup, x, y)` containing true fluxes ($\text{DN}/\text{s}$).
* **`Q_cube`** *(numpy.ndarray)*: 3D array of shape `(nint * ngroup, x, y)` containing trapped charges ($\text{DN}$).
* **`S_cube`** *(numpy.ndarray)*: 3D array of shape `(nint * ngroup, x, y)` containing observed rates ($\text{DN}/\text{s}$).

---

## `fit_persistence(ramp_cube, timestamps, epsilon=0.18, tau=120.0)`

Analytically fits the trapped-charge persistence model to each pixel's entire time series to solve for the 2D arrays of true incident flux ($F$) and initial trapped charge ($Q_{\text{init}}$) at $t = 0$.

#### Parameters:
* **`ramp_cube`** *(numpy.ndarray)*: 4D array of shape `(nint, ngroup, x, y)`.
* **`timestamps`** *(numpy.ndarray)*: 2D array of shape `(nint, ngroup)` containing absolute or relative timestamps (in seconds).
* **`epsilon`** *(float or numpy.ndarray)*: Trapping efficiency (default: `0.18`).
* **`tau`** *(float or numpy.ndarray)*: Trapping decay constant in seconds (default: `120.0`).

#### Returns:
* **`F_fit`** *(numpy.ndarray)*: 2D array of shape `(x, y)` containing the fitted true flux rate ($\text{DN}/\text{s}$).
* **`Q_init_fit`** *(numpy.ndarray)*: 2D array of shape `(x, y)` containing the fitted initial trapped charge ($\text{DN}$).

---

### Usage Example:

```python
from pandora_tools import get_target_files, read_InfImg, fit_persistence, calculate_persistence, plot_persistence_model

# 1. Load data
files = get_target_files("WASP-18b", "InfImg")
ramp_cube, times_sec = read_InfImg(files[0], time_format="seconds")

# 2. Fit the persistence model to solve for 2D arrays of F and Q_init
F_fit, Q_init_fit = fit_persistence(
    ramp_cube=ramp_cube,
    timestamps=times_sec,
    epsilon=0.18,
    tau=120.0
)

# 3. Calculate the full time-series model using the fitted Q_init
P_cube, F_cube, Q_cube, S_cube = calculate_persistence(
    ramp_cube=ramp_cube,
    timestamps=times_sec,
    epsilon=0.18,
    tau=120.0,
    Q_init=Q_init_fit
)

# 4. Plot the model fit for pixel at center (x=40, y=125) and save to image
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