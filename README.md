### Database Schema in  pandora_observations.db 

  The database contains a  files  table structured as follows:

  •  id  (INTEGER PRIMARY KEY)
  •  filepath  (TEXT UNIQUE): Absolute path to the FITS file.
  •  filename  (TEXT): Name of the file.
  •  target_id  (TEXT): Target name (sourced from FITS  TARG_ID  with standard filename fallback).
  •  file_type  (TEXT):  InfImg  or  VisSci .
  •  instrument  (TEXT): Instrument name ( NIRDA  or  VISDA ).
  •  camera_id  (TEXT): Camera name ( H2rgCam  or  PcoCam ).
  •  exptime  (REAL): Exposure time.
  •  ra  /  dec  (REAL): Coordinates ( TARG_RA / TARG_DEC ).
  •  obs_time  (TEXT): UTC date/time of observation parsed from filename.
  •  month  /  day  (TEXT): Year subfolder components.
  ──────
  ### CLI Tool Usage

  Always run the script using the virtual environment's Python binary:
   ./.pandora_aux/bin/python3 pandora_get_targets.py <command> [args] 

  #### 1. Re-index/Update Database

  If new files are added to the observation directories, you can run the  index  command. It uses a cache set to skip files that are already
  indexed, making updates run in a fraction of a second.

    ./.pandora_aux/bin/python3 pandora_get_targets.py index

  #### 2. List All Observed Targets and File Counts

  Shows a formatted table of all targets, counting their  InfImg  and  VisSci  files.

    ./.pandora_aux/bin/python3 pandora_get_targets.py targets

  Tip: You can search/filter for targets containing a specific substring, e.g., to find all WASP planets:

    ./.pandora_aux/bin/python3 pandora_get_targets.py targets --search WASP

  #### 3. View/Export Files for a Specific Target

  To view details and paths of all FITS files matching a target (e.g.  WASP-18b ):

    ./.pandora_aux/bin/python3 pandora_get_targets.py files --target WASP-18b

  To filter by file type ( InfImg  or  VisSci ) and save the absolute file paths to a text file (useful for automated processing pipelines):

    ./.pandora_aux/bin/python3 pandora_get_targets.py files --target WASP-18b --type InfImg --output wasp18b_inf_files.txt