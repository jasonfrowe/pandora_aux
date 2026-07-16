import os
import sqlite3
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.time import Time
from astropy.visualization import ZScaleInterval
from numba import njit, prange

def get_target_files(target_id, file_type):
    """
    Queries the pandora_observations database for FITS files of a specific 
    target and observation file type.
    
    Parameters:
      target_id (str): The name of the target (e.g., 'WASP-18b')
      file_type (str): The type of file ('InfImg' or 'VisSci')
      
    Returns:
      list: A list of absolute file paths sorted chronologically.
    """
    # Resolve the path to the database file in the same directory as this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    db_path = os.path.join(script_dir, "pandora_observations.db")
    
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Database not found at: {db_path}. Please index first.")
        
    # Connect to the SQLite database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Query for the file paths. 
    # Use LOWER() to ensure case-insensitivity on the target ID.
    query = """
        SELECT filepath 
        FROM files 
        WHERE LOWER(target_id) = LOWER(?) AND file_type = ?
        ORDER BY obs_time
    """
    
    cursor.execute(query, (target_id, file_type))
    
    # Extract file paths from the rows
    filepaths = [row[0] for row in cursor.fetchall()]
    
    # Clean up connection
    conn.close()
    
    return filepaths

def read_InfImg(filepath, time_format="JD"):
    """
    Reads a single InfImg FITS file and extracts the science data cube and timestamps.
    
    Parameters:
      filepath (str): Path to the InfImg FITS file.
      time_format (str): Format of the returned timestamps.
                         Options: 'JD' (Julian Date), 
                                  'MJD' (Modified Julian Date relative to Jan 1, 2026),
                                  or 'seconds' (seconds since 2000-01-01 00:00:00 UTC).
                         Default: 'JD'.
      
    Returns:
      tuple: (ramp_cube, timestamps)
        - ramp_cube (numpy.ndarray): 4D array with shape (nint, ngroup, x, y)
        - timestamps (numpy.ndarray): 2D array with shape (nint, ngroup) containing timestamps.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"FITS file not found at: {filepath}")
        
    with fits.open(filepath) as hdul:
        prim_hdr = hdul[0].header
        sci_hdr = hdul[1].header
        
        # Get dimensions from the primary and science headers
        nint = prim_hdr.get("INTEGRTS")
        ngroup = prim_hdr.get("GRPS")
        x = sci_hdr.get("NAXIS1")  # Dimension X (e.g., 80)
        y = sci_hdr.get("NAXIS2")  # Dimension Y (e.g., 250)
        
        # Read the science data
        # Astropy loads 3D FITS images as shape (NAXIS3, NAXIS2, NAXIS1) -> (nint * ngroup, y, x)
        data = hdul[1].data
        
        # Transpose the last two dimensions to change layout from (y, x) to (x, y)
        data_xy = np.transpose(data, (0, 2, 1))  # Shape becomes (nint * ngroup, x, y)
        
        # Reshape to 4D datacube with dimensions (nint, ngroup, x, y)
        ramp_cube = data_xy.reshape(nint, ngroup, x, y)
        
        # Calculate time in seconds since 2000-01-01 00:00:00 UTC
        # CORSTIME: coarse start time in seconds
        # FINETIME: fine start time in nanoseconds
        # FRMTIME: frame time in milliseconds
        t0_sec = prim_hdr.get("CORSTIME", 0.0) + prim_hdr.get("FINETIME", 0.0) * 1e-9
        frmtime_sec = prim_hdr.get("FRMTIME", 0.0) * 1e-3  # Convert milliseconds to seconds
        
        # Build 2D seconds-since-epoch array of shape (nint, ngroup)
        frame_indices = np.arange(nint * ngroup).reshape(nint, ngroup)
        times_sec = t0_sec + frame_indices * frmtime_sec
        
        # Convert to target format
        fmt = time_format.upper()
        if fmt == "JD":
            # JD epoch for 2000-01-01 00:00:00 UTC is 2451544.5
            timestamps = 2451544.5 + times_sec / 86400.0
        elif fmt == "MJD":
            # Return MJD relative to Jan 1, 2026 00:00:00 UTC (MJD offset: 61041.0)
            # J2000 midnight calendar epoch MJD is 51544.0.
            # Difference in days: 61041.0 - 51544.0 = 9497.0 days.
            timestamps = (times_sec / 86400.0) - 9497.0
        elif fmt == "SECONDS":
            timestamps = times_sec
        else:
            raise ValueError(f"Unsupported time_format: '{time_format}'. Choose from 'JD', 'MJD', or 'seconds'.")
        
    return ramp_cube, timestamps

def plot_ramp_cube(ramp_cube, integration_index=0, group_index=None, iraf_contrast=0.25, cmap="viridis", output_path=None):
    """
    Plots a single 2D frame from the 4D ramp_cube using IRAF-style zscale for optimal display scaling.
    
    Parameters:
      ramp_cube (numpy.ndarray): 4D array with shape (nint, ngroup, x, y).
      integration_index (int): Index of the integration to plot (default: 0).
      group_index (int): Index of the group to plot. If None, defaults to the last group (ngroup - 1).
      iraf_contrast (float): Contrast parameter for the zscale algorithm (default: 0.25).
      cmap (str): Matplotlib colormap (default: 'viridis').
      output_path (str): File path to save the plot. If None, display the plot interactively.
    """
    if ramp_cube.ndim != 4:
        raise ValueError(f"Expected a 4D ramp_cube; got array with shape {ramp_cube.shape}")
        
    nint, ngroup, x, y = ramp_cube.shape
    
    if group_index is None:
        group_index = ngroup - 1
        
    if not (0 <= integration_index < nint):
        raise IndexError(f"integration_index {integration_index} is out of bounds for nint={nint}")
    if not (0 <= group_index < ngroup):
        raise IndexError(f"group_index {group_index} is out of bounds for ngroup={ngroup}")
        
    # Extract the 2D frame (shape: x, y)
    frame = ramp_cube[integration_index, group_index]
    
    finite_mask = np.isfinite(frame)
    if not np.any(finite_mask):
        raise ValueError("Selected frame has no finite values for display scaling.")
        
    # Calculate IRAF-style zscale limits
    zscale = ZScaleInterval(contrast=iraf_contrast)
    vmin, vmax = zscale.get_limits(frame[finite_mask])
    
    # Calculate size matching the aspect ratio of the subarray
    ny_size, nx_size = frame.shape  # rows (x / spatial) and columns (y / dispersion)
    fig_width = 10.0
    fig_height = max(3.0, fig_width * (ny_size / max(nx_size, 1)))
    
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    
    # Matplotlib imshow maps row index to y-axis and column index to x-axis
    im = ax.imshow(
        frame,
        origin="lower",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        aspect="equal",
        interpolation="nearest"
    )
    
    ax.set_title(f"Ramp Frame (Int: {integration_index}, Group: {group_index})")
    ax.set_xlabel("Dispersion pixel")
    ax.set_ylabel("Spatial pixel")
    
    # Add colorbar
    fig.colorbar(im, ax=ax, label="Counts (DN)", pad=0.02, shrink=0.8)
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=150)
        plt.close(fig)
        print(f"Plot saved successfully to: {output_path}")
    else:
        plt.show()

@njit(parallel=True)
def _fit_ramp_numba(timestamps, ramp_cube, read_noise, gain):
    nint, ngroup, nx, ny = ramp_cube.shape
    slopes = np.zeros((nint, nx, ny), dtype=np.float64)
    intercepts = np.zeros((nint, nx, ny), dtype=np.float64)
    
    is_ts_2d = (timestamps.ndim == 2)
    
    for i in prange(nint):
        # Extract times for this integration
        if is_ts_2d:
            t = timestamps[i]
        else:
            t = timestamps
            
        for x in range(nx):
            for y in range(ny):
                sum_w = 0.0
                sum_wt = 0.0
                sum_ws = 0.0
                sum_wtt = 0.0
                sum_wts = 0.0
                
                for g in range(ngroup):
                    tg = t[g]
                    sg = ramp_cube[i, g, x, y]
                    
                    # Optimal weighting variance: read_noise^2 + Poisson noise
                    variance = (read_noise ** 2) + (max(sg, 0.0) / gain)
                    w = 1.0 / variance
                    
                    sum_w += w
                    sum_wt += w * tg
                    sum_ws += w * sg
                    sum_wtt += w * tg * tg
                    sum_wts += w * tg * sg
                    
                delta = sum_w * sum_wtt - sum_wt * sum_wt
                if delta != 0.0:
                    slope = (sum_w * sum_wts - sum_wt * sum_ws) / delta
                    intercept = (sum_ws - slope * sum_wt) / sum_w
                    slopes[i, x, y] = slope
                    intercepts[i, x, y] = intercept
                else:
                    slopes[i, x, y] = 0.0
                    intercepts[i, x, y] = 0.0
                    
    return slopes, intercepts

def fit_ramp(timestamps, ramp_cube, read_noise=10.0, gain=1.0):
    """
    Fits an optimal weighted least squares line to each pixel's ramp across groups.
    
    Parameters:
      timestamps (numpy.ndarray): 1D array of shape (ngroup) or 2D array of shape (nint, ngroup).
      ramp_cube (numpy.ndarray): 4D array of shape (nint, ngroup, x, y).
      read_noise (float): Read noise in electrons (default: 10.0).
      gain (float): Detector gain in e-/DN (default: 1.0).
      
    Returns:
      tuple: (slope_cube, intercept_cube)
        - slope_cube (numpy.ndarray): 3D array of shape (nint, x, y) containing the slopes (DN/s).
        - intercept_cube (numpy.ndarray): 3D array of shape (nint, x, y) containing the intercepts (DN).
    """
    if ramp_cube.ndim != 4:
        raise ValueError(f"Expected a 4D ramp_cube; got shape {ramp_cube.shape}")
        
    ramp_cube = np.asarray(ramp_cube, dtype=np.float64)
    timestamps = np.asarray(timestamps, dtype=np.float64)
    
    return _fit_ramp_numba(timestamps, ramp_cube, read_noise, gain)

@njit(parallel=True)
def _calc_persistence_numba(ramp_cube, times_flat, epsilon, tau, Q_init):
    nint, ngroup, nx, ny = ramp_cube.shape
    N = nint * ngroup
    
    # Outputs
    P_cube = np.zeros((N, nx, ny), dtype=np.float64)
    F_cube = np.zeros((N, nx, ny), dtype=np.float64)
    Q_cube = np.zeros((N, nx, ny), dtype=np.float64)
    S_cube = np.zeros((N, nx, ny), dtype=np.float64)  # observed signal rate
    
    # Time steps
    dt = np.zeros(N, dtype=np.float64)
    dt[0] = times_flat[1] - times_flat[0]
    for k in range(1, N):
        dt[k] = times_flat[k] - times_flat[k-1]
        
    for x in prange(nx):
        for y in range(ny):
            eps_val = epsilon[x, y]
            tau_val = tau[x, y]
            
            # Initial trapped charge
            Q_prev = Q_init[x, y]
            
            # 1. Calculate observed signal rate S for all steps
            for i in range(nint):
                k0 = i * ngroup
                C0 = ramp_cube[i, 0, x, y]
                C1 = ramp_cube[i, 1, x, y]
                dt1 = times_flat[k0 + 1] - times_flat[k0]
                S_cube[k0, x, y] = (C1 - C0) / dt1
                
                for g in range(1, ngroup):
                    k = k0 + g
                    C_curr = ramp_cube[i, g, x, y]
                    C_prev = ramp_cube[i, g-1, x, y]
                    dt_curr = times_flat[k] - times_flat[k-1]
                    S_cube[k, x, y] = (C_curr - C_prev) / dt_curr
            
            # 2. Run persistence filter across the flattened time series
            for k in range(N):
                dt_k = dt[k]
                S_k = S_cube[k, x, y]
                
                # Decay factor
                decay = np.exp(-dt_k / tau_val)
                
                # Persistence current rate (DN/s)
                P_k = Q_prev * (1.0 - decay) / dt_k
                P_cube[k, x, y] = P_k
                
                # True flux rate (DN/s)
                F_k = S_k - P_k
                F_cube[k, x, y] = F_k
                
                # Update trapped charge
                Q_k = Q_prev * decay + eps_val * F_k * dt_k
                Q_cube[k, x, y] = Q_k
                
                Q_prev = Q_k
                
    return P_cube, F_cube, Q_cube, S_cube

def calculate_persistence(ramp_cube, timestamps, epsilon=0.18, tau=120.0, Q_init=0.0):
    """
    Computes a trapped-charge persistence model for each pixel.
    
    Parameters:
      ramp_cube (numpy.ndarray): 4D array of shape (nint, ngroup, x, y).
      timestamps (numpy.ndarray): 2D array of shape (nint, ngroup) containing UTC timestamps in seconds.
      epsilon (float or numpy.ndarray): Trapping efficiency. Can be a scalar or a 2D array of shape (x, y).
      tau (float or numpy.ndarray): Trapping decay constant in seconds. Can be a scalar or a 2D array of shape (x, y).
      Q_init (float or numpy.ndarray): Initial trapped charge. Can be a scalar or a 2D array of shape (x, y).
      
    Returns:
      tuple: (P_cube, F_cube, Q_cube, S_cube)
        - P_cube (numpy.ndarray): 3D array of shape (nint * ngroup, x, y) containing persistence rates (DN/s).
        - F_cube (numpy.ndarray): 3D array of shape (nint * ngroup, x, y) containing true fluxes (DN/s).
        - Q_cube (numpy.ndarray): 3D array of shape (nint * ngroup, x, y) containing trapped charges (DN).
        - S_cube (numpy.ndarray): 3D array of shape (nint * ngroup, x, y) containing observed rates (DN/s).
    """
    nint, ngroup, nx, ny = ramp_cube.shape
    
    # Coerce epsilon, tau, and Q_init to (nx, ny) arrays
    if isinstance(epsilon, (int, float)):
        epsilon_arr = np.full((nx, ny), float(epsilon), dtype=np.float64)
    else:
        epsilon_arr = np.asarray(epsilon, dtype=np.float64)
        
    if isinstance(tau, (int, float)):
        tau_arr = np.full((nx, ny), float(tau), dtype=np.float64)
    else:
        tau_arr = np.asarray(tau, dtype=np.float64)
        
    if isinstance(Q_init, (int, float)):
        Q_init_arr = np.full((nx, ny), float(Q_init), dtype=np.float64)
    else:
        Q_init_arr = np.asarray(Q_init, dtype=np.float64)
        
    # Flatten timestamps
    times_flat = np.asarray(timestamps, dtype=np.float64).flatten()
    
    # Convert ramp_cube to float64 for fits
    ramp_cube_dbl = np.asarray(ramp_cube, dtype=np.float64)
    
    return _calc_persistence_numba(ramp_cube_dbl, times_flat, epsilon_arr, tau_arr, Q_init_arr)

def fit_persistence(ramp_cube, timestamps, epsilon=0.18, tau=120.0):
    """
    Fits the trapped-charge persistence model to the FITS ramp_cube data
    to analytically solve for the 2D arrays of true flux (F_fit) and 
    initial trapped charge (Q_init_fit) per pixel.
    
    Parameters:
      ramp_cube (numpy.ndarray): 4D array of shape (nint, ngroup, x, y).
      timestamps (numpy.ndarray): 2D array of shape (nint, ngroup) containing UTC timestamps in seconds.
      epsilon (float or numpy.ndarray): Trapping efficiency. Can be a scalar or a 2D array of shape (x, y).
      tau (float or numpy.ndarray): Trapping decay constant in seconds. Can be a scalar or a 2D array of shape (x, y).
      
    Returns:
      tuple: (F_fit, Q_init_fit)
        - F_fit (numpy.ndarray): 2D array of shape (x, y) containing the fitted true flux rate (DN/s).
        - Q_init_fit (numpy.ndarray): 2D array of shape (x, y) containing the fitted initial trapped charge (DN).
    """
    nint, ngroup, nx, ny = ramp_cube.shape
    N = nint * ngroup
    
    # Coerce epsilon and tau to (nx, ny) arrays
    if isinstance(epsilon, (int, float)):
        eps_arr = np.full((nx, ny), float(epsilon), dtype=np.float64)
    else:
        eps_arr = np.asarray(epsilon, dtype=np.float64)
        
    if isinstance(tau, (int, float)):
        tau_arr = np.full((nx, ny), float(tau), dtype=np.float64)
    else:
        tau_arr = np.asarray(tau, dtype=np.float64)
        
    # Flatten timestamps and calculate time steps dt
    times_flat = np.asarray(timestamps, dtype=np.float64).flatten()
    dt = np.zeros(N, dtype=np.float64)
    dt[0] = times_flat[1] - times_flat[0]
    dt[1:] = np.diff(times_flat)
    
    # Compute observed signal rate S_cube of shape (N, nx, ny)
    S_cube = np.zeros((N, nx, ny), dtype=np.float64)
    for i in range(nint):
        k0 = i * ngroup
        C0 = ramp_cube[i, 0].astype(np.float64)
        C1 = ramp_cube[i, 1].astype(np.float64)
        dt1 = times_flat[k0 + 1] - times_flat[k0]
        S_cube[k0] = (C1 - C0) / dt1
        
        for g in range(1, ngroup):
            k = k0 + g
            C_curr = ramp_cube[i, g].astype(np.float64)
            C_prev = ramp_cube[i, g-1].astype(np.float64)
            dt_curr = times_flat[k] - times_flat[k-1]
            S_cube[k] = (C_curr - C_prev) / dt_curr
            
    # Initialize linear least squares accumulation matrices
    M00 = np.zeros((nx, ny), dtype=np.float64)
    M01 = np.zeros((nx, ny), dtype=np.float64)
    M11 = np.zeros((nx, ny), dtype=np.float64)
    V0 = np.zeros((nx, ny), dtype=np.float64)
    V1 = np.zeros((nx, ny), dtype=np.float64)
    
    # a_k and b_k represent the linear coefficients of Q_k:
    # Q_k = a_k * Q_init + b_k * F
    a_prev = np.ones((nx, ny), dtype=np.float64)
    b_prev = np.zeros((nx, ny), dtype=np.float64)
    
    # Forward pass to accumulate coefficients
    for k in range(N):
        dt_k = dt[k]
        
        d_k = np.exp(-dt_k / tau_arr)
        h_k = (1.0 - d_k) / dt_k
        
        # Coefficients for S_model_k = A_k * Q_init + B_k * F
        A_k = a_prev * h_k
        B_k = 1.0 + b_prev * h_k
        
        # Update coefficients for the next step Q_k
        a_curr = a_prev * d_k
        b_curr = b_prev * d_k + eps_arr * dt_k
        
        # Accumulate sums
        S_k = S_cube[k]
        M00 += A_k * A_k
        M01 += A_k * B_k
        M11 += B_k * B_k
        V0 += A_k * S_k
        V1 += B_k * S_k
        
        a_prev = a_curr
        b_prev = b_curr
        
    # Solve the system of equations for each pixel:
    # [M00  M01] [Q_init] = [V0]
    # [M01  M11] [F     ]   [V1]
    det = M00 * M11 - M01 * M01
    
    # Avoid division by zero
    det_mask = (det > 1e-12)
    Q_init_fit = np.zeros((nx, ny), dtype=np.float64)
    F_fit = np.zeros((nx, ny), dtype=np.float64)
    
    # Solve where determinant is valid
    Q_init_fit[det_mask] = (M11[det_mask] * V0[det_mask] - M01[det_mask] * V1[det_mask]) / det[det_mask]
    F_fit[det_mask] = (M00[det_mask] * V1[det_mask] - M01[det_mask] * V0[det_mask]) / det[det_mask]
    
    # Enforce physical constraints: F >= 0, Q_init >= 0
    Q_init_fit = np.clip(Q_init_fit, 0.0, None)
    F_fit = np.clip(F_fit, 0.0, None)
    
    return F_fit, Q_init_fit

def plot_persistence_model(timestamps, S_cube, F_cube, P_cube, Q_cube, x_pixel, y_pixel, output_path=None):
    """
    Plots the observed signal rate, true flux, persistence current, and trapped charge 
    over time for a specific pixel.
    """
    times_flat = np.asarray(timestamps).flatten()
    # Relative time in seconds from start of observation
    time_rel = times_flat - times_flat[0]
    
    # Extract values for the selected pixel
    S = S_cube[:, x_pixel, y_pixel]
    F = F_cube[:, x_pixel, y_pixel]
    P = P_cube[:, x_pixel, y_pixel]
    Q = Q_cube[:, x_pixel, y_pixel]
    
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    
    # 1. Observed Signal and True Flux
    axes[0].plot(time_rel, S, label="Observed Signal Rate (S)", color="#1f77b4", alpha=0.8)
    axes[0].plot(time_rel, F, label="True Flux Rate (F)", color="#2ca02c", linestyle="--", alpha=0.8)
    axes[0].set_ylabel("Rate (DN/s)")
    axes[0].set_title(f"Persistence Model for Pixel (x={x_pixel}, y={y_pixel})")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # 2. Persistence Current
    axes[1].plot(time_rel, P, label="Persistence Current (P)", color="#ff7f0e", alpha=0.8)
    axes[1].set_ylabel("Persistence (DN/s)")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    # 3. Trapped Charge
    axes[2].plot(time_rel, Q, label="Trapped Charge (Q)", color="#d62728", alpha=0.8)
    axes[2].set_xlabel("Time (seconds since start of observation)")
    axes[2].set_ylabel("Trapped Charge (DN)")
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=150)
        plt.close(fig)
        print(f"Persistence plot saved successfully to: {output_path}")
    else:
        plt.show()

if __name__ == "__main__":
    # Demonstration of how to use get_target_files, read_InfImg, plot_ramp_cube, fit_ramp, calculate_persistence, and fit_persistence
    target = "WASP-18b"
    ftype = "InfImg"
    
    print(f"1. Querying database for {target} {ftype} files...")
    try:
        files = get_target_files(target, ftype)
        print(f"Found {len(files)} file(s).")
        
        if len(files) > 0:
            first_file = files[0]
            print(f"\n2. Reading first file: {os.path.basename(first_file)}...")
            with fits.open(first_file) as hdul:
                nint = hdul[0].header.get("INTEGRTS")
                ngroup = hdul[0].header.get("GRPS")
                
            ramp_cube, times_sec = read_InfImg(first_file, time_format="seconds")
            
            print(f"\n3. Datacube shape: {ramp_cube.shape}")
            
            # Save a sample plot
            output_png = "ramp_plot.png"
            mid_int = nint // 2
            print(f"\n4. Generating sample plot...")
            plot_ramp_cube(
                ramp_cube=ramp_cube,
                integration_index=mid_int,
                group_index=ngroup - 1,
                output_path=output_png
            )
            
            # Run Least Squares Ramp Fitting
            print(f"\n5. Fitting ramp data using optimal weighted least squares (relative times)...")
            times_rel = times_sec - times_sec[:, 0:1]
            slopes, intercepts = fit_ramp(times_rel, ramp_cube)
            
            print(f"\n6. Slope & Intercept Cubes properties:")
            print(f"   Slope cube shape      : {slopes.shape}")
            print(f"   Intercept cube shape  : {intercepts.shape}")
            print(f"   Slope range (DN/s)    : [{slopes.min():.3f}, {slopes.max():.3f}]")
            print(f"   Intercept range (DN)  : [{intercepts.min():.3f}, {intercepts.max():.3f}]")
            
            px, py = slopes.shape[1] // 2, slopes.shape[2] // 2
            print(f"\n7. Sample pixel (x={px}, y={py}) for integration {mid_int}:")
            print(f"   Slope                 : {slopes[mid_int, px, py]:.3f} DN/s")
            print(f"   Intercept (Bias)      : {intercepts[mid_int, px, py]:.3f} DN")
            
            # Run analytical persistence fitting to get 2D arrays F_fit and Q_init_fit
            print(f"\n8. Performing analytical linear least squares fit of the persistence model...")
            # We use epsilon=0.18 and tau=120.0
            F_fit, Q_init_fit = fit_persistence(ramp_cube, times_sec, epsilon=0.18, tau=120.0)
            
            print(f"\n9. Fitted Parameters Properties:")
            print(f"   Fitted True Flux F_fit shape      : {F_fit.shape}")
            print(f"   Fitted Initial Charge Q_init shape : {Q_init_fit.shape}")
            print(f"   F_fit range (DN/s)                 : [{F_fit.min():.3f}, {F_fit.max():.3f}]")
            print(f"   Q_init_fit range (DN)              : [{Q_init_fit.min():.3f}, {Q_init_fit.max():.3f}]")
            print(f"   Sample pixel (x={px}, y={py}) fitted:")
            print(f"     True Incident Flux (F_fit)       : {F_fit[px, py]:.3f} DN/s")
            print(f"     Initial Trapped Charge (Q_init)  : {Q_init_fit[px, py]:.3f} DN")
            
            # Generate the fitted forward model using our fitted parameters
            print(f"\n10. Calculating model using fitted parameters...")
            P_cube_fit, F_cube_fit, Q_cube_fit, S_cube_fit = calculate_persistence(
                ramp_cube=ramp_cube,
                timestamps=times_sec,
                epsilon=0.18,
                tau=120.0,
                Q_init=Q_init_fit
            )
            
            # Save the fitted persistence model plot
            fit_plot_png = "persistence_fit_plot.png"
            print(f"\n11. Saving fitted persistence model plot to {fit_plot_png}...")
            plot_persistence_model(
                timestamps=times_sec,
                S_cube=S_cube_fit,
                F_cube=F_cube_fit,
                P_cube=P_cube_fit,
                Q_cube=Q_cube_fit,
                x_pixel=px,
                y_pixel=py,
                output_path=fit_plot_png
            )
        else:
            print("No files found to read.")
    except Exception as e:
        print(f"An error occurred: {e}")
