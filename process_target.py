# %%

# Enable autoreload when running in IPython/Jupyter without using % magics
# so this file remains valid Python for Pylance.
_get_ipython = globals().get("get_ipython")
if callable(_get_ipython):
    _ip = _get_ipython()
    if _ip is not None:
        _ip.run_line_magic("load_ext", "autoreload")
        _ip.run_line_magic("autoreload", "2")

import numpy as np
import pandora_tools as pandora

# %%
target = "WASP-18b"
ftype = "InfImg"
files = pandora.get_target_files(target, ftype)

print(f"Successfully retrieved {len(files)} {ftype} file(s) for {target}:\n")
for idx, filepath in enumerate(files, 1):
    print(f"{idx}: {filepath}")

# %%
ramp_cube, timestamps = pandora.read_InfImg(files[0])
print(f"\nRamp cube shape: {ramp_cube.shape}")
# %%
print(timestamps)
# %%
