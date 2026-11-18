import numpy as np
import matplotlib.pyplot as plt
 
# ---- inputs ---------------------------------------------------------------
OMEGA = 733.04          # rad/s (7000 rpm)
R_IN, R_WALL = 0.030, 0.045   # m, free surface -> containment wall
P0 = 10e6               # Pa, core gas pressure at the free surface
P_BHAT = 14.589e6       # Pa, Bhat Table I inner-wall pressure (target)
N = 200                 # finite-difference steps
 
data = np.loadtxt("keese_fig4_2_fuel_T.csv", delimiter=",", skiprows=1)
r_data, T_data = data[:, 0], data[:, 1]
 
# 4th-order fit in normalized radius x = (r - 0.03)/0.015
FIT = [3634.08, -226.92, -2106.03, 3288.50, -3087.60]
def T_of_r(r):
    x = (r - R_IN) / (R_WALL - R_IN)
    return sum(c * x**k for k, c in enumerate(FIT))
 
def rho_U(T):
    """Liquid uranium density, Bhat Eq. 2 (Fischer 2000)."""
    return 17270 - 1.4485 * (T - 1408)
 
# ---- first-order forward-difference march of Eq. 1 ------------------------
r = np.linspace(R_IN, R_WALL, N + 1)
dr = r[1] - r[0]
T = np.interp(r, r_data, T_data)
rho = rho_U(T)
 
P = np.empty_like(r)
P[0] = P0
for i in range(N):
    P[i + 1] = P[i] + dr * OMEGA**2 * rho[i] * r[i]
 
err = (P[-1] - P_BHAT) / P_BHAT * 100
print(f"P_wall = {P[-1]/1e6:.3f} MPa  (Bhat: {P_BHAT/1e6:.3f} MPa, {err:+.2f}%)")

# ---- first-order forward-difference march of Eq. 1 ------------------------
r = np.linspace(R_IN, R_WALL, N + 1)
dr = r[1] - r[0]
T = np.interp(r, r_data, T_data)
rho = rho_U(T)
 
P = np.empty_like(r)
P[0] = P0
for i in range(N):
    P[i + 1] = P[i] + dr * OMEGA**2 * rho[i] * r[i]
 
err = (P[-1] - P_BHAT) / P_BHAT * 100
print(f"P_wall = {P[-1]/1e6:.3f} MPa  (Bhat: {P_BHAT/1e6:.3f} MPa, {err:+.2f}%)")
 
# ---- plots ----------------------------------------------------------------
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e6e5e1"
BLUE, ORANGE = "#2a78d6", "#eb6834"   # categorical slots 1, 2
 
plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlecolor": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
})
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
r_cm = r * 100
 
# (a) temperature: extracted data + fit
r_fine = np.linspace(R_IN, R_WALL, 400)
ax1.plot(r_fine * 100, T_of_r(r_fine), color=ORANGE, lw=2, label="4th-order fit", zorder=2)
ax1.plot(r_data[::3] * 100, T_data[::3], "o", ms=6, color=BLUE, mec="#fcfcfb", mew=1,
         label="Keese Fig. 4.2 (extracted)", zorder=3)
ax1.set(xlabel="Radius (cm)", ylabel="Temperature (K)", ylim=(1000, 4000),
        title="(a) Uranium temperature, input")
ax1.legend(frameon=False, loc="lower left", labelcolor=INK2)
ax1.annotate(f"Keese wall ≈ {T_data[-1]:.0f} K\n(Bhat BC: 1478.1 K)", (r_data[-1]*100, T_data[-1]),
             xytext=(3.75, 1700), color=INK2, ha="right",
             arrowprops=dict(arrowstyle="-", color=MUTED, lw=1))
 
# (b) pressure: FD result + Bhat target
ax2.plot(r_cm, P / 1e6, color=BLUE, lw=2, label=f"Forward Euler, N = {N}")
ax2.plot(R_WALL * 100, P_BHAT / 1e6, "D", ms=8, color=ORANGE, mec="#fcfcfb", mew=1.5,
         label="Bhat Table I (14.589 MPa)", zorder=3)
ax2.set(xlabel="Radius (cm)", ylabel="Pressure (MPa)", ylim=(10, 15),
        title="(b) Uranium pressure, Eq. 1")
ax2.legend(frameon=False, loc="upper left", labelcolor=INK2)
ax2.annotate(f"{P[-1]/1e6:.3f} MPa ({err:+.2f}%)", (4.5, P[-1] / 1e6),
             xytext=(3.35, 13.2), color=INK2,
             arrowprops=dict(arrowstyle="-", color=MUTED, lw=1))
 
fig.savefig("step0_uranium.png", dpi=180)
print("saved step0_uranium.png")