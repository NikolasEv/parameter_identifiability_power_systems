"""Independent check of the FIM analysis using only NumPy and SciPy.

The Julia MBAM code (../src) holds x''_d fixed at its nominal value in the
model equations (see "Known issue" in ../README.md). This script:

  1. reproduces the saved Jacobian ../data/SG11_jacobian.txt by replicating
     that behavior, which confirms the model, parameterization and settings; and
  2. recomputes the FIM with x''_d varied consistently (x''_q = x''_d, since
     delta x_5 = 0), and reports its spectrum and the parameters spanning the
     six stiffest eigendirections.

Sensitivities are computed by central finite differences in log-parameters,
with sigma = 1 and all six states observed at t = 3:0.02:5.

Usage:  python fim_check.py      (takes about 15 seconds)
"""

from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

# Constants (src/pvalues.jl)
OMEGA_B = 2 * np.pi * 60
OMEGA_0 = 1.0
PM = 0.7
V = 1.09
VF0 = 4.14676
THETA = -0.0382548

# Initial condition (src/ic.jl)
Y0 = [0.5, 0.98, 2.13, 0.02, 1.98, 0.02]

# Independent parameters in the order used by src/SG11.jl
NAMES = ["H", "D", "dx1", "dx2", "dx3", "dx4", "x2d", "dTd", "dTq", "T2d0", "T2q0"]
P0 = np.array([2.53, 0.5, 0.12, 2.02, 1.932, 0.448, 0.48, 4.69, 1.29, 0.06, 0.21])
X2D_NOMINAL = 0.48

T_OBS = np.round(np.arange(3, 5 + 1e-9, 0.02), 10)


def physical(p):
    """Map independent parameters to physical ones (delta x_5 = 0)."""
    H, D, dx1, dx2, dx3, dx4, x2d, dTd, dTq, T2d0, T2q0 = p
    x2q = x2d
    x1d = x2q + dx4
    x1q = x1d + dx3
    xq = x1q + dx2
    xd = xq + dx1
    return dict(H=H, D=D, xd=xd, x1d=x1d, x2d=x2d, xq=xq, x1q=x1q, x2q=x2q,
                T1d0=dTd + T2d0, T2d0=T2d0, T1q0=dTq + T2q0, T2q0=T2q0)


def trajectory(logp, freeze_x2d):
    """Full 6th-order model (Sec. III of the paper), observed at T_OBS."""
    q = physical(np.exp(logp))
    x2d = X2D_NOMINAL if freeze_x2d else q["x2d"]

    def rhs(_, y):
        d, w, e1q, e1d, e2q, e2d = y
        vd, vq = V * np.sin(d - THETA), V * np.cos(d - THETA)
        i_d = (e2q - vq) / x2d
        i_q = (vd - e2d) / q["x2q"]
        Pg = vd * i_d + vq * i_q
        return [
            OMEGA_B * (w - OMEGA_0),
            (PM - Pg - q["D"] * (w - OMEGA_0)) / q["H"],
            (-e1q - (q["xd"] - q["x1d"]) * i_d + VF0) / q["T1d0"],
            (-e1d + (q["xq"] - q["x1q"]) * i_q) / q["T1q0"],
            (-e2q + e1q - (q["x1d"] - x2d) * i_d) / q["T2d0"],
            (-e2d + e1d + (q["x1q"] - q["x2q"]) * i_q) / q["T2q0"],
        ]

    sol = solve_ivp(rhs, (0, T_OBS[-1]), Y0, t_eval=T_OBS, method="Radau",
                    rtol=1e-12, atol=1e-12)
    return sol.y.ravel()  # state-major ordering, matching the saved Jacobian


def jacobian(freeze_x2d, h=1e-4):
    logp0 = np.log(P0)
    cols = []
    for i in range(len(P0)):
        e = np.zeros(len(P0))
        e[i] = h
        cols.append((trajectory(logp0 + e, freeze_x2d)
                     - trajectory(logp0 - e, freeze_x2d)) / (2 * h))
    return np.array(cols).T


def report(label, J):
    _, s, vt = np.linalg.svd(J, full_matrices=False)
    weight = (vt[:6] ** 2).sum(axis=0)
    print(f"\n{label}")
    print("  FIM eigenvalues:", np.array2string(s**2, precision=2, max_line_width=120))
    print("  Weight of each parameter in the six stiffest directions:")
    for name, wgt in sorted(zip(NAMES, weight), key=lambda x: -x[1]):
        print(f"    {name:5s} {wgt:.3f}")


def main():
    saved = np.loadtxt(Path(__file__).resolve().parent.parent / "data" / "SG11_jacobian.txt")

    J_frozen = jacobian(freeze_x2d=True)
    # The saved Jacobian is the residual Jacobian (data - model), hence the sign.
    rel = np.linalg.norm(-J_frozen - saved) / np.linalg.norm(saved)
    print(f"As coded (x''_d frozen): relative difference from saved Jacobian = {rel:.1e}")
    report("Saved Jacobian (as used in the paper)", saved)

    report("Corrected (x''_d varied consistently)", jacobian(freeze_x2d=False))


if __name__ == "__main__":
    main()
