"""Example 4: Residual thermodynamic properties.

Demonstrates:
    - Computing HR, SR, GR, VR, Cp_R, Cv_R from all three EOS
    - Comparing PR, SRK, and PR78 residual enthalpy for propane
    - Scanning temperature at fixed pressure to show the liquid-to-vapor
      transition in residual properties
    - Physical interpretation of sign conventions

System: pure propane (C3H8) scanned from 200 K (liquid) to 380 K (vapor)
        at P = 5 bar.  Critical point: Tc = 369.83 K, Pc = 42.48 bar.
"""

import numpy as np

from pvttool import (
    add_components,
    add_mixture,
    add_thermo,
    pr78eos,
    preos,
    srkeos,
)

comps = add_components(["Propane"])
P = 5e5  # 5 bar

# ---------------------------------------------------------------------------
# 1. Residual properties vs temperature (liquid and vapor branches)
# ---------------------------------------------------------------------------
print("=" * 70)
print("Residual properties of propane at P = 5 bar (PR EOS)")
print("=" * 70)
print(f"{'T [K]':>6}  {'Phase':>6}  {'HR [J/mol]':>12}  {'SR [J/molK]':>12}  "
      f"{'GR [J/mol]':>12}  {'VR [m3/mol]':>13}")
print("-" * 70)

thermo_liq = add_thermo()
thermo_liq.phase = 1
thermo_liq.fugacity_switch = 0

thermo_vap = add_thermo()
thermo_vap.phase = 2
thermo_vap.fugacity_switch = 0

for T in [200, 230, 260, 290, 320, 350, 380]:
    mix = add_mixture(comps, T=float(T), p=P)
    mix.mole_fraction = np.array([1.0])

    # Below Tc (369.83 K): show liquid Z-root; above: vapor Z-root
    thermo = thermo_liq if T < 369.83 else thermo_vap
    phase_label = "liquid" if T < 369.83 else "vapor"

    _, _, _, HR, props = preos(mix, thermo)
    print(f"{T:>6d}  {phase_label:>6}  {HR:>12.1f}  {props['SR']:>12.4f}  "
          f"{props['GR']:>12.1f}  {props['VR']:>13.2e}")

# ---------------------------------------------------------------------------
# 2. EOS comparison: residual enthalpy at 250 K, 10 bar (liquid propane)
# ---------------------------------------------------------------------------
print()
print("=" * 50)
print("EOS comparison: liquid propane at 250 K, 10 bar")
print("=" * 50)

mix_cmp = add_mixture(comps, T=250.0, p=10e5)
mix_cmp.mole_fraction = np.array([1.0])

for eos_fn, label in [(preos, "PR  "), (srkeos, "SRK "), (pr78eos, "PR78")]:
    t = add_thermo()
    t.eos = eos_fn
    t.phase = 1
    t.fugacity_switch = 0
    _, _, _, HR, props = eos_fn(mix_cmp, t)
    print(f"  {label}:  HR = {HR:>9.2f} J/mol   "
          f"SR = {props['SR']:>8.4f} J/(mol·K)   "
          f"VR = {props['VR']:>10.3e} m³/mol")

# ---------------------------------------------------------------------------
# 3. Cp_R and Cv_R for liquid propane
# ---------------------------------------------------------------------------
print()
print("=" * 50)
print("Heat-capacity departures at 250 K, 10 bar (liquid, PR EOS)")
print("=" * 50)
mix_cp = add_mixture(comps, T=250.0, p=10e5)
mix_cp.mole_fraction = np.array([1.0])
t_cp = add_thermo()
t_cp.phase = 1
t_cp.fugacity_switch = 0
_, _, _, _, props = preos(mix_cp, t_cp)
print(f"  Cp_R = {props['Cp_R']:>8.3f} J/(mol·K)")
print(f"  Cv_R = {props['Cv_R']:>8.3f} J/(mol·K)")
print(f"  GR   = {props['GR']:>8.1f} J/mol  (should equal HR - T*SR)")
GR_check = props["HR"] - 250.0 * props["SR"]
print(f"  HR - T·SR check: {GR_check:.1f} J/mol  (diff: {abs(GR_check - props['GR']):.2e})")
