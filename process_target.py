# %%
import matplotlib.pyplot as plt
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
ramp_cube, timestamps, start_timestamps = pandora.read_InfImg(files[0], time_format="MJD", return_start_times=True)
print(f"\nRamp cube shape: {ramp_cube.shape}")
# %%
print(timestamps)
# %%

pandora.plot_ramp_cube(
                ramp_cube=ramp_cube,
                integration_index=0,
                group_index=ramp_cube.shape[1] - 1,
            )
# %%
times_sec = timestamps * 24 * 3600
times_rel = times_sec - times_sec[:, 0:1]
slopes, intercepts = pandora.fit_ramp(times_rel, ramp_cube)
# %%
integration_index = 40
px = 40
py = 100

plt.plot(timestamps[integration_index, :], ramp_cube[integration_index, :, px, py], marker="o", linestyle="-")
model = slopes[integration_index, px, py] * times_rel[integration_index, :] + intercepts[integration_index, px, py]
plt.plot(timestamps[integration_index, :], model, color="red", label="Fitted line")
plt.xlabel("Time (MJD)")
plt.ylabel("Counts (DN)")
plt.show()
# %%
# Working on a persistence model for a single pixel (px, py) across all integrations and groups
times_sec_flat = times_sec.flatten()
ramp_cube_flat = ramp_cube.reshape(-1, ramp_cube.shape[2], ramp_cube.shape[3])

plt.plot(times_sec_flat - times_sec_flat[0], ramp_cube_flat[:, px, py], marker="o", linestyle="-")

S = np.zeros(ramp_cube_flat.shape[0])  # Observed signal for the pixel
Q = np.zeros(ramp_cube_flat.shape[0])  # Trapped charge for the pixel
F = np.zeros(ramp_cube_flat.shape[0])  # True charge rate for the pixel 
P = np.zeros(ramp_cube_flat.shape[0])  # Rate of persistence released  
M = np.zeros(ramp_cube_flat.shape[0])  # Model for the pixel

bias = ramp_cube_flat[0, px, py]  # Bias level for the pixel
Qo = 0.0  # Initial trapped charge for the pixel
tau = 120.0  # Decay time constant for the pixel
eps = 0.18  # Efficiency of the pixel

i = 0
dt = times_sec_flat[1] - times_sec_flat[0]  # Time difference between consecutive groups
Q[i] = Qo * np.exp(-dt / tau)  # Initial trapped charge for the pixel
S[i] = (ramp_cube_flat[1, px, py] - ramp_cube_flat[0, px, py]) / dt  # Initial flux for the pixel
P[i] = Qo * (1 - np.exp(-dt / tau))  # Initial persistence for the pixel

M[i] = S[i] + bias

plt.plot(times_sec_flat - times_sec_flat[0], M, color="red", label="Pixel model")

plt.xlim(0, times_sec_flat[12] - times_sec_flat[0])
plt.ylim(20000, 27000)

plt.xlabel("Time (s)")
plt.ylabel("Counts (DN)")
plt.show()


# %%
F_fit, Q_init_fit = pandora.fit_persistence(ramp_cube, times_sec, epsilon=0.18, tau=120.0)
P_cube, F_cube, Q_cube, S_cube = pandora.calculate_persistence(
    ramp_cube=ramp_cube,
    timestamps=times_sec,
    epsilon=0.18,
    tau=120.0,
    Q_init=Q_init_fit
)

pandora.plot_persistence_model(times_sec, S_cube, F_cube, P_cube, Q_cube, px, py)
# %%
