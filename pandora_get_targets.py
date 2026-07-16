#!/usr/bin/env python3
"""
Pandora Target Database Tool

This script indexes Pandora FITS files from the specified observation directories 
into a SQLite database, allowing users to query observed targets and generate 
lists of corresponding InfImg or VisSci files for processing.

Usage:
  # 1. Index or update the database:
  python3 pandora_get_targets.py index

  # 2. List all observed targets and their file counts:
  python3 pandora_get_targets.py targets

  # 3. Filter targets containing a specific substring:
  python3 pandora_get_targets.py targets --search WASP

  # 4. List FITS files for a specific target:
  python3 pandora_get_targets.py files --target WASP-18b

  # 5. List only InfImg or VisSci FITS files for a target and save to a text file:
  python3 pandora_get_targets.py files --target WASP-18b --type InfImg --output wasp18b_inf_files.txt
"""

import os
import sys
import re
import glob
import sqlite3
import argparse
from astropy.io import fits

DEFAULT_DATA_DIR = "/opt/data2/rowe/pandora/2026"
DEFAULT_FOLDERS = ["01", "02", "03", "04", "05", "06", "07"]

def get_db_path():
    """Returns the absolute path to the SQLite database in the script's directory."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(script_dir, "pandora_observations.db")

def get_db_connection():
    """Establishes connection to the SQLite database."""
    db_path = get_db_path()
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    return conn

def init_db(conn):
    """Initializes the database schema and indexes."""
    cursor = conn.cursor()
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS files (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        filepath TEXT UNIQUE,
        filename TEXT,
        target_id TEXT,
        file_type TEXT,
        instrument TEXT,
        camera_id TEXT,
        exptime REAL,
        ra REAL,
        dec REAL,
        obs_time TEXT,
        month TEXT,
        day TEXT
    )
    """)
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_target_id ON files(target_id);")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_file_type ON files(file_type);")
    conn.commit()

def parse_filename(filename):
    """
    Parses metadata from the standard Pandora filename pattern:
    YYYY-MM-DD__HH-MM-SS_(InfImg|VisSci)_<Target>_<details>.fits
    """
    # Matches YYYY-MM-DD__HH-MM-SS_(InfImg|VisSci)_rest.fits
    m = re.match(r'^(\d{4})-(\d{2})-(\d{2})__(\d{2})-(\d{2})-(\d{2})_(InfImg|VisSci)_(.+)\.fits$', filename)
    if m:
        year, month, day, hour, minute, second, file_type, rest = m.groups()
        obs_time = f"{year}-{month}-{day} {hour}:{minute}:{second}"
        # Fallback target extraction: split on the first '_d' followed by a digit
        target_match = re.split(r'_d\d', rest)
        target_fallback = target_match[0] if target_match else None
        return {
            "obs_time": obs_time,
            "month": month,
            "day": day,
            "file_type": file_type,
            "target_fallback": target_fallback
        }
    return None

def index_files(data_dir, folders, force=False):
    """Scans directories and indexes FITS files into the database."""
    conn = get_db_connection()
    init_db(conn)
    
    # Get already indexed files to allow fast incremental updates
    cursor = conn.cursor()
    cursor.execute("SELECT filepath FROM files")
    existing_paths = {row["filepath"] for row in cursor.fetchall()}
    
    print(f"Database contains {len(existing_paths)} previously indexed files.")
    if force:
        print("Force flag is set. Re-indexing all files.")
        existing_paths.clear()

    # Discover all FITS files in requested folders
    fits_files = []
    print(f"Scanning FITS files in {data_dir} for folders: {', '.join(folders)}...")
    for folder in folders:
        folder_path = os.path.join(data_dir, folder)
        if not os.path.exists(folder_path):
            print(f"Warning: Folder {folder_path} does not exist. Skipping.")
            continue
        # Recursively find all *.fits files
        pattern = os.path.join(folder_path, "**", "*.fits")
        found = glob.glob(pattern, recursive=True)
        fits_files.extend(found)

    total_files = len(fits_files)
    print(f"Found a total of {total_files} FITS files on disk.")

    # Filter out files that are already indexed
    files_to_index = [f for f in fits_files if f not in existing_paths]
    to_index_count = len(files_to_index)
    
    if to_index_count == 0:
        print("All files are already indexed. Database is up to date.")
        conn.close()
        return

    print(f"Indexing {to_index_count} new FITS files...")
    
    new_records = []
    skipped_count = 0

    for idx, filepath in enumerate(files_to_index, 1):
        filename = os.path.basename(filepath)
        parsed = parse_filename(filename)
        
        if not parsed:
            # Not a matching pattern, skip or try header only
            parsed = {
                "obs_time": None,
                "month": None,
                "day": None,
                "file_type": "Unknown",
                "target_fallback": None
            }

        # Read FITS header for target name and other properties
        try:
            header = fits.getheader(filepath)
            target_id = header.get("TARG_ID")
            instrument = header.get("INSTRMNT")
            camera_id = header.get("CAMERAID")
            exptime = header.get("EXPTIME")
            ra = header.get("TARG_RA")
            dec = header.get("TARG_DEC")
        except Exception as e:
            # If FITS header can't be read, log warning and use fallback if possible
            print(f"\nWarning: Could not read FITS header for {filename}: {e}", file=sys.stderr)
            target_id = parsed.get("target_fallback")
            instrument = camera_id = exptime = ra = dec = None
            skipped_count += 1

        # Fallback to filename target if header didn't specify one
        if not target_id:
            target_id = parsed.get("target_fallback")

        # Fallback to filename for type if FITS instrument maps standard types
        file_type = parsed.get("file_type") or "Unknown"

        new_records.append((
            filepath,
            filename,
            target_id,
            file_type,
            instrument,
            camera_id,
            exptime,
            ra,
            dec,
            parsed.get("obs_time"),
            parsed.get("month"),
            parsed.get("day")
        ))

        # Show progress
        if idx % 100 == 0 or idx == to_index_count:
            print(f"Processed {idx}/{to_index_count} files...", end="\r", flush=True)

    print() # Newline after progress bar

    # Insert in bulk using a transaction
    if new_records:
        cursor.executemany("""
        INSERT OR REPLACE INTO files (
            filepath, filename, target_id, file_type, instrument, camera_id, 
            exptime, ra, dec, obs_time, month, day
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, new_records)
        conn.commit()

    print(f"Successfully indexed {len(new_records)} files (skipped {skipped_count} invalid files).")
    conn.close()

def print_table(headers, rows):
    """Helper function to print rows as a formatted text table."""
    if not rows:
        print("No data available.")
        return
        
    widths = [len(h) for h in headers]
    for row in rows:
        for idx, val in enumerate(row):
            widths[idx] = max(widths[idx], len(str(val if val is not None else "")))
    
    header_str = " | ".join(f"{str(val).ljust(widths[idx])}" for idx, val in enumerate(headers))
    print(header_str)
    print("-" * len(header_str))
    
    for row in rows:
        print(" | ".join(f"{str(val if val is not None else '').ljust(widths[idx])}" for idx, val in enumerate(row)))

def list_targets(search_query=None):
    """Lists observed targets with counts of corresponding FITS files."""
    conn = get_db_connection()
    init_db(conn)
    cursor = conn.cursor()
    
    query = """
        SELECT 
            target_id, 
            SUM(case when file_type='InfImg' then 1 else 0 end) as inf_count, 
            SUM(case when file_type='VisSci' then 1 else 0 end) as vis_count, 
            COUNT(*) as total_count 
        FROM files 
    """
    params = []
    if search_query:
        query += " WHERE target_id LIKE ? "
        params.append(f"%{search_query}%")
        
    query += " GROUP BY target_id ORDER BY target_id"
    
    cursor.execute(query, params)
    rows = cursor.fetchall()
    conn.close()
    
    if not rows:
        if search_query:
            print(f"No targets found matching '{search_query}'.")
        else:
            print("No targets found in the database. Run the 'index' command first.")
        return

    headers = ["Target ID", "InfImg Files", "VisSci Files", "Total Files"]
    table_rows = [(r["target_id"], r["inf_count"], r["vis_count"], r["total_count"]) for r in rows]
    
    print(f"\n--- Observed Targets ({len(rows)} unique targets found) ---")
    print_table(headers, table_rows)

def get_files_for_target(target, file_type=None, output_file=None):
    """Retrieves file paths for a target and optionally filters/saves them."""
    conn = get_db_connection()
    init_db(conn)
    cursor = conn.cursor()
    
    query = "SELECT filepath, filename, file_type, obs_time FROM files WHERE LOWER(target_id) = LOWER(?)"
    params = [target]
    
    if file_type:
        query += " AND LOWER(file_type) = LOWER(?)"
        params.append(file_type)
        
    query += " ORDER BY obs_time"
    
    cursor.execute(query, params)
    rows = cursor.fetchall()
    conn.close()
    
    if not rows:
        filter_str = f" of type {file_type}" if file_type else ""
        print(f"No files found for target '{target}'{filter_str}.")
        return

    # Extract paths
    filepaths = [r["filepath"] for r in rows]

    if output_file:
        try:
            with open(output_file, 'w') as f:
                for path in filepaths:
                    f.write(f"{path}\n")
            print(f"Successfully wrote {len(filepaths)} file paths to {output_file}")
        except Exception as e:
            print(f"Error writing to output file {output_file}: {e}", file=sys.stderr)
    else:
        print(f"\n--- Files for target '{target}' ({len(filepaths)} found) ---")
        headers = ["File Type", "Observation Time", "Filename"]
        table_rows = [(r["file_type"], r["obs_time"], r["filename"]) for r in rows]
        print_table(headers, table_rows)
        print("\nPaths (first 10 shown):")
        for path in filepaths[:10]:
            print(path)
        if len(filepaths) > 10:
            print(f"... and {len(filepaths) - 10} more paths.")

def main():
    parser = argparse.ArgumentParser(
        description="Manage the Pandora FITS observations database.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    subparsers = parser.add_subparsers(dest="command", help="Commands")

    # Index command
    index_parser = subparsers.add_parser("index", help="Scan and index FITS files")
    index_parser.add_argument("--data-dir", default=DEFAULT_DATA_DIR, help=f"Directory to scan (default: {DEFAULT_DATA_DIR})")
    index_parser.add_argument("--folders", default=",".join(DEFAULT_FOLDERS), help=f"Comma-separated list of subfolders (default: {','.join(DEFAULT_FOLDERS)})")
    index_parser.add_argument("--force", action="store_true", help="Force re-indexing of all files (skip cache)")

    # Targets command
    targets_parser = subparsers.add_parser("targets", help="List observed targets and their file counts")
    targets_parser.add_argument("--search", help="Filter target IDs by substring")

    # Files command
    files_parser = subparsers.add_parser("files", help="List or export files for a specific target")
    files_parser.add_argument("--target", required=True, help="Target ID (case-insensitive, e.g. WASP-18b)")
    files_parser.add_argument("--type", choices=["InfImg", "VisSci"], help="Filter by file type (InfImg or VisSci)")
    files_parser.add_argument("--output", help="Save the list of absolute file paths to a text file")

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(0)

    if args.command == "index":
        folders_list = [f.strip() for f in args.folders.split(",") if f.strip()]
        index_files(args.data_dir, folders_list, force=args.force)
    elif args.command == "targets":
        list_targets(search_query=args.search)
    elif args.command == "files":
        get_files_for_target(args.target, file_type=args.type, output_file=args.output)

if __name__ == "__main__":
    main()
