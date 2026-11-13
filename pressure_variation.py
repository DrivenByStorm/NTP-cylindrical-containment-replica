import numpy as np

omega = 733.04                      # rad/s
r = np.linspace(0.03, 0.045, 201)   # N = 200 steps
dr = r[1] - r[0]

T = T_of_r(r)                       # interpolate digitized Fig. 1
rho = 17270 - 1.4485*(T - 1408)

P = np.empty_like(r)
P[0] = 10e6
for i in range(len(r) - 1):
    P[i+1] = P[i] + dr * omega**2 * rho[i] * r[i]

print(P[-1]/1e6)  # target: ~14.589 MPa