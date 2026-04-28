"""Example 1: VLE flash — methanol + water system.

Demonstrates:
    - Loading components from the database
    - Setting binary interaction parameters (kij)
    - Running vle_flash_negative to trace a VLE curve
    - Verifying the component mass balance
    - Comparing computed vapor compositions against published experimental data

Reference:
    Ochi, K. et al. (1986) J. Chem. Eng. Data, measured at 322.91 K.
"""

import numpy as np

from pvttool import (
    FlashOptions,
    add_components,
    add_mixture,
    add_thermo,
    vle_flash,
    vle_flash_negative,
)

# ---------------------------------------------------------------------------
# System setup
# ---------------------------------------------------------------------------
comps = add_components(["Methanol", "Water"])
T = 322.91          # K — isothermal experiment
thermo = add_thermo()   # PR EOS, vdW mixing rule
opts = FlashOptions()

# kij = -0.07 (PR EOS binary interaction parameter for MeOH-H2O)
kij = np.array([[0.0, -0.07], [-0.07, 0.0]])

# Experimental data (Ochi et al. 1986)
p_exp   = [15932, 26131, 34474, 41369, 52142]    # Pa
y1_exp  = [0.2741, 0.6294, 0.7472, 0.8348, 0.9736]  # vapor methanol mole fraction

# ---------------------------------------------------------------------------
# Trace the VLE curve with equimolar feed using negative flash
# ---------------------------------------------------------------------------
print("=" * 60)
print("Methanol-Water VLE at 322.91 K  (PR EOS, kij = -0.07)")
print("=" * 60)
print(f"{'P [Pa]':>10}  {'V_frac':>8}  {'y_MeOH':>8}  {'x_MeOH':>8}  {'y_exp':>8}  {'|Δy|':>7}")
print("-" * 60)

for p, y_ref in zip(p_exp, y1_exp):
    mix = add_mixture(comps, T=T, p=p)
    mix.bip.eos_cons = kij

    y, x, V = vle_flash_negative(mix, thermo, opts)

    delta = abs(y[0] - y_ref) if 0.0 <= V <= 1.0 else float("nan")
    print(f"{p:>10d}  {V:>8.4f}  {y[0]:>8.4f}  {x[0]:>8.4f}  {y_ref:>8.4f}  {delta:>7.4f}")

# ---------------------------------------------------------------------------
# Standard vle_flash for a two-phase feed
# ---------------------------------------------------------------------------
print()
print("Standard vle_flash with z(MeOH) = 0.35 at P = 26131 Pa")
mix_tp = add_mixture(comps, T=T, p=26131)
mix_tp.mole_fraction = np.array([0.35, 0.65])
mix_tp.bip.eos_cons = kij

y, x, V = vle_flash(mix_tp, thermo, opts)

print(f"  Vapor fraction  : {V:.4f}")
print(f"  Vapor  y = {y}")
print(f"  Liquid x = {x}")

# Mass balance check
z = mix_tp.mole_fraction
z_calc = V * y + (1.0 - V) * x
print(f"  Mass balance error: {np.max(np.abs(z_calc - z)):.2e}  (should be < 1e-5)")
