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
data = pandora.read_InfImg(files[0], time_format="MJD")
ramp_cube = data.ramp_cube
timestamps = data.timestamps
start_timestamps = data.start_timestamps
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

bias = intercepts[0, px, py]  # Guess at bias level for the pixel
Qo = 0.0  # Initial trapped charge for the pixel
tau = 120.0  # Decay time constant for the pixel
eps = 0.18  # Efficiency of the pixel

# Time from reset to Group 0 midpoint (s)
dt0 = (data.drops1 + (data.reads - 1) / 2.0) * data.frmtime_sec  

# 1. Observed signal rate (DN/sec)
S[0] = (ramp_cube_flat[0, px, py] - bias) / dt0  

# 2. Persistence rate (DN/sec)
decay = np.exp(-dt0 / tau)
P[0] = Qo * (1.0 - decay) / dt0  

# 3. True flux rate (DN/sec)
F[0] = S[0] - P[0]  

# 4. Trapped charge at the end of interval 0 (DN)
Q[0] = Qo * decay + eps * F[0] * dt0  

N = len(times_sec_flat)
for i in range(1, N):
    dt = times_sec_flat[i] - times_sec_flat[i-1]
    decay = np.exp(-dt / tau)

    if i % data.ngroup == 0:
        # Start of a new integration: calculate S[i] using the next group (i to i+1)
        dt_local = times_sec_flat[i+1] - times_sec_flat[i]
        S[i] = (ramp_cube_flat[i+1, px, py] - ramp_cube_flat[i, px, py]) / dt_local
    else:
        # Normal step within the same integration
        S[i] = (ramp_cube_flat[i, px, py] - ramp_cube_flat[i-1, px, py]) / dt

    P[i] = Q[i-1] * (1.0 - decay) / dt
    F[i] = S[i] - P[i]
    Q[i] = Q[i-1] * decay + eps * F[i] * dt



# # Extract the fitted true flux for this specific pixel
# F_pixel = F

# for k in range(N):
#     if k % data.ngroup == 0:
#         # Start of an integration: anchor the model to the first group readout
#         M[k] = ramp_cube_flat[k, px, py]
#     else:
#         # Within the same integration: accumulate signal
#         dt = times_sec_flat[k] - times_sec_flat[k-1]
        
#         # Modeled signal rate = True Flux + Modeled Persistence rate
#         S_model = F_pixel + P_cube[k, px, py]
        
#         # Integrate (counts at previous step + rate * time)
#         M[k] = M[k-1] + S_model * dt

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
