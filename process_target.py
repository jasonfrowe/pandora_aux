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

# --- Arrays & Parameters ---
N = len(times_sec_flat)
S = np.zeros(N)
Q = np.zeros(N)
F = np.zeros(N)
P = np.zeros(N)
M = np.zeros(N)
Mp = np.zeros(N)

bias = intercepts[0, px, py]   
tau = 120.0  
eps = 0.01  # Fraction of traps opening per physically generated photoelectron

# --- Initial State ---
dt0 = (data.drops1 + (data.reads - 1) / 2.0) * data.frmtime_sec  
star_rate = slopes[0, px, py]  

Q[0] = 0.0  
P[0] = 0.0
F[0] = star_rate
S[0] = star_rate 
M[0] = bias + S[0] * dt0       
Mp[0] = bias + F[0] * dt0      

time_since_reset = dt0 

# --- Up the Ramp Loop ---
for i in range(1, N):
    
    if i % data.ngroup == 0:
        # --- RESET BOUNDARY ---
        dt = (data.drops1 + (data.reads - 1) / 2.0) * data.frmtime_sec 
        is_reset = True
        time_since_reset = dt  # Restart the integration clock
    else:
        # --- NORMAL READ ---
        dt = times_sec_flat[i] - times_sec_flat[i-1]
        is_reset = False
        time_since_reset += dt # Advance the integration clock

    decay = np.exp(-dt / tau)
    
    # 1. Release charge from the trap pool (dominates EARLY in the ramp)
    released_charge = Q[i-1] * (1.0 - decay)
    P[i] = released_charge / dt 
    
    # 2. Trap new charge (dominates LATE in the ramp as signal accumulates)
    true_accumulated = star_rate * time_since_reset
    trapped_charge = eps * true_accumulated * dt
    
    # 3. Update the trapped pool
    Q[i] = Q[i-1] - released_charge + trapped_charge
    
    # 4. Calculate observed rate: (True) + (Released) - (Lost to traps)
    S[i] = star_rate + P[i] - (trapped_charge / dt)
    F[i] = star_rate
    
    # 5. Build the ramps
    if is_reset:
        M[i]  = bias + S[i] * dt
        Mp[i] = bias + F[i] * dt   # FIXED: Mp now resets correctly!
    else:
        M[i]  = M[i-1] + S[i] * dt
        Mp[i] = Mp[i-1] + F[i] * dt
        

Nplot = 100
# plt.plot(times_sec_flat - times_sec_flat[0], 0 * ramp_cube_flat[:, px, py], marker="o", linestyle="-", label="Observed")
# plt.plot(times_sec_flat - times_sec_flat[0], M - Mp, color="red", label="Pixel model", marker="o", linestyle="-")
# plt.plot(times_sec_flat - times_sec_flat[0], ramp_cube_flat[:, px, py], marker="o", linestyle="-", label="Observed")
plt.plot(times_sec_flat - times_sec_flat[0], Mp, color="blue", label="Pixel model", marker="o", linestyle="-")
plt.plot(times_sec_flat - times_sec_flat[0], M, color="red", label="Pixel model", marker="o", linestyle="-")
plt.xlim(times_sec_flat[0] - times_sec_flat[0] ,times_sec_flat[Nplot] - times_sec_flat[0])
plt.show()

# %%
F_fit, Q_init_fit = pandora.fit_persistence(ramp_cube, times_sec, epsilon=0.18, tau=120.0)
res = pandora.calculate_persistence(
    ramp_cube=ramp_cube,
    timestamps=times_sec,
    epsilon=0.18,
    tau=120.0,
    Q_init=Q_init_fit
)

pandora.plot_persistence_model(
    times_sec, 
    res.S_cube, 
    res.F_cube, 
    res.P_cube, 
    res.Q_cube, 
    px, 
    py, 
    ramp_cube=ramp_cube, 
    M_cube=res.M_cube, 
    Mp_cube=res.Mp_cube
)
# %%
