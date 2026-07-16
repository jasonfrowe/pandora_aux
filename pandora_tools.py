import os
import sqlite3
import numpy as np
from astropy.io import fits
from astropy.time import Time

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
                         Options: 'JD' (Julian Date), 'MJD' (Modified Julian Date),
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
            # MJD epoch for 2000-01-01 00:00:00 UTC is 51544.0
            timestamps = 51544.0 + times_sec / 86400.0
        elif fmt == "SECONDS":
            timestamps = times_sec
        else:
            raise ValueError(f"Unsupported time_format: '{time_format}'. Choose from 'JD', 'MJD', or 'seconds'.")
        
    return ramp_cube, timestamps

if __name__ == "__main__":
    # Demonstration of how to use get_target_files and read_InfImg
    target = "WASP-18b"
    ftype = "InfImg"
    
    print(f"1. Querying database for {target} {ftype} files...")
    try:
        files = get_target_files(target, ftype)
        print(f"Found {len(files)} file(s).")
        
        if len(files) > 0:
            first_file = files[0]
            print(f"\n2. Reading first file (JD format): {os.path.basename(first_file)}...")
            ramp_cube, jd_times = read_InfImg(first_file, time_format="JD")
            
            print(f"\n3. Reading first file (MJD format):...")
            _, mjd_times = read_InfImg(first_file, time_format="MJD")
            
            # Print Properties
            print(f"\n4. Datacube and Timestamps Properties:")
            print(f"   Ramp Cube Shape (nint, ngroup, x, y) : {ramp_cube.shape}")
            print(f"   Timestamps Shape (nint, ngroup)      : {jd_times.shape}")
            print(f"   Data type                            : {ramp_cube.dtype}")
            print(f"   Data range                           : [{ramp_cube.min()}, {ramp_cube.max()}]")
            
            print(f"\n5. Timestamp Details:")
            print(f"   First Frame:")
            print(f"     JD   : {jd_times[0, 0]:.8f}")
            print(f"     MJD  : {mjd_times[0, 0]:.8f}")
            print(f"     UTC  : {Time(jd_times[0, 0], format='jd').iso}")
            print(f"   Last Frame:")
            print(f"     JD   : {jd_times[-1, -1]:.8f}")
            print(f"     MJD  : {mjd_times[-1, -1]:.8f}")
            print(f"     UTC  : {Time(jd_times[-1, -1], format='jd').iso}")
        else:
            print("No files found to read.")
    except Exception as e:
        print(f"An error occurred: {e}")
