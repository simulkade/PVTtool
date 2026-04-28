"""Example 5: Multi-component natural gas VLE.

Demonstrates:
    - A 6-component natural gas mixture (typical North Sea composition)
    - Running a stability test to confirm two-phase behaviour
    - VLE flash with custom FlashOptions
    - Bubble and dew point pressures at reservoir temperature
    - Comparing how K-values change from Wilson correlation to EOS values

Mixture (mole fractions):
    CH4   0.70
    C2H6  0.10
    C3H8  0.08
    n-C4  0.05
    n-C5  0.04
    n-C10 0.03
"""

import numpy as np

from pvttool import (
    FlashOptions,
    add_components,
    add_mixture,
    add_thermo,
    bubble_pressure,
    dew_pressure,
    stability_test,
    vle_flash,
)
from pvttool.flash._kvalue import kval_estimate

# ---------------------------------------------------------------------------
# Mixture definition
# ---------------------------------------------------------------------------
names = ["Methane", "Ethane", "Propane", "n-Butane", "n-Pentane", "n-Decane"]
z     = np.array([0.70, 0.10, 0.08, 0.05, 0.04, 0.03])

comps  = add_components(names)
thermo = add_thermo()
opts   = FlashOptions(accuracy=1e-7, iteration=200)

print("=" * 65)
print("Natural gas VLE — 6-component mixture")
print("=" * 65)
print("Composition (mol%): " + "  ".join(f"{n}: {zi*100:.0f}" for n, zi in zip(names, z)))

# ---------------------------------------------------------------------------
# 1. Reservoir conditions: 380 K, 15 MPa
# ---------------------------------------------------------------------------
T_res, P_res = 380.0, 15e6
mix = add_mixture(comps, T=T_res, p=P_res)
mix.mole_fraction = z.copy()

print(f"\n--- Reservoir conditions: T = {T_res} K, P = {P_res/1e6:.0f} MPa ---")

_, _, _, result = stability_test(mix, thermo, opts)
print(f"Stability: {result['overall']}")

if result["overall"] == "unstable":
    y, x, V = vle_flash(mix, thermo, opts)
    print(f"\nVLE Flash result:")
    print(f"  Vapor fraction V = {V:.4f}")
    print(f"\n  {'Component':<12} {'z':>6}  {'y (vapor)':>10}  {'x (liquid)':>10}  {'K = y/x':>10}")
    print("  " + "-" * 52)
    for i, name in enumerate(names):
        Ki = y[i] / x[i] if x[i] > 1e-10 else float("inf")
        print(f"  {name:<12} {z[i]:>6.3f}  {y[i]:>10.4f}  {x[i]:>10.4f}  {Ki:>10.4f}")
else:
    print("Single-phase at reservoir conditions.")

# ---------------------------------------------------------------------------
# 2. Separator conditions: 310 K, 3 MPa
# ---------------------------------------------------------------------------
T_sep, P_sep = 310.0, 3e6
mix_sep = add_mixture(comps, T=T_sep, p=P_sep)
mix_sep.mole_fraction = z.copy()

print(f"\n--- Separator conditions: T = {T_sep} K, P = {P_sep/1e6:.0f} MPa ---")

_, _, _, result_sep = stability_test(mix_sep, thermo, opts)
print(f"Stability: {result_sep['overall']}")

if result_sep["overall"] == "unstable":
    y_s, x_s, V_s = vle_flash(mix_sep, thermo, opts)
    print(f"  Vapor fraction V = {V_s:.4f}")

# ---------------------------------------------------------------------------
# 3. Phase envelope: bubble and dew pressure at 310 K
# ---------------------------------------------------------------------------
T_scan = 310.0
mix_scan = add_mixture(comps, T=T_scan, p=3e6)
mix_scan.mole_fraction = z.copy()

P_bub, y_bub, ok_bub = bubble_pressure(mix_scan, thermo, opts)
P_dew, x_dew, ok_dew = dew_pressure(mix_scan, thermo, opts)

print(f"\n--- Saturation pressures at T = {T_scan} K ---")
if ok_bub:
    print(f"  Bubble-point pressure : {P_bub/1e6:.4f} MPa")
    # For methane-rich mixtures the bubble point is near the critical region;
    # the incipient vapor composition converges toward the feed composition.
    print(f"  Incipient vapor       : " + "  ".join(f"{yi:.4f}" for yi in y_bub))
else:
    print("  Bubble-point: did not converge")
if ok_dew:
    print(f"  Dew-point pressure    : {P_dew/1e6:.4f} MPa")
    print(f"  Incipient liquid      : " + "  ".join(f"{xi:.4f}" for xi in x_dew))
else:
    print("  Dew-point: did not converge")

# ---------------------------------------------------------------------------
# 4. Wilson vs EOS K-values at separator conditions
# ---------------------------------------------------------------------------
print(f"\n--- K-value comparison at {T_sep} K, {P_sep/1e6:.0f} MPa ---")
mix_k = add_mixture(comps, T=T_sep, p=P_sep)
mix_k.mole_fraction = z.copy()

K_wilson = kval_estimate(mix_k)

if result_sep["overall"] == "unstable":
    K_eos = y_s / np.where(x_s > 1e-10, x_s, 1e-10)
    print(f"\n  {'Component':<12} {'K_Wilson':>10}  {'K_EOS':>10}")
    print("  " + "-" * 35)
    for i, name in enumerate(names):
        print(f"  {name:<12} {K_wilson[i]:>10.4f}  {K_eos[i]:>10.4f}")
else:
    print(f"\n  {'Component':<12} {'K_Wilson':>10}")
    for i, name in enumerate(names):
        print(f"  {name:<12} {K_wilson[i]:>10.4f}")
