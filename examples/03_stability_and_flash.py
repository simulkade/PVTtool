"""Example 3: Phase stability test followed by conditional VLE flash.

Demonstrates:
    - Running stability_test to determine whether a mixture will split
    - Using the result to decide whether a VLE flash is needed
    - Scanning T and P to find the two-phase boundary
    - Interpreting the result struct (flags, saturation values, message)

System: methane (60 mol%) + n-decane (40 mol%) — a highly asymmetric
hydrocarbon pair that is two-phase over a wide pressure range at 300 K.
"""

import numpy as np

from pvttool import (
    FlashOptions,
    add_components,
    add_mixture,
    add_thermo,
    stability_test,
    vle_flash,
)

# ---------------------------------------------------------------------------
# System setup
# ---------------------------------------------------------------------------
comps = add_components(["Methane", "n-Decane"])
z = np.array([0.60, 0.40])
thermo = add_thermo()
opts = FlashOptions()
flash_opts = FlashOptions(accuracy=1e-7, iteration=150)

# ---------------------------------------------------------------------------
# 1. Single-point stability test at 300 K, 5 MPa
# ---------------------------------------------------------------------------
T, P = 300.0, 5e6
mix = add_mixture(comps, T=T, p=P)
mix.mole_fraction = z.copy()

flags, SL, SV, result = stability_test(mix, thermo, opts)

print("=" * 60)
print("Stability test: 60% CH4 / 40% n-C10  at 300 K, 5 MPa")
print("=" * 60)
print(result["message"])
print(f"\nFlags : {flags}")
print(f"SV    : {SV:.4f}  (trial vapor saturation from Test 1)")
print(f"SL    : {SL:.4f}  (trial liquid saturation from Test 2)")

# ---------------------------------------------------------------------------
# 2. If unstable → run VLE flash
# ---------------------------------------------------------------------------
if result["overall"] == "unstable":
    print("\nMixture is unstable — running VLE flash …")
    y, x, V = vle_flash(mix, thermo, flash_opts)
    print(f"\n  Vapor fraction V = {V:.4f}")
    print(f"  Vapor phase  y  = {y}")
    print(f"  Liquid phase x  = {x}")

    # Verify equilibrium: Ki = yi/xi should equal phi_Li/phi_Vi
    K = y / np.where(x > 1e-10, x, 1e-10)
    print(f"\n  K-values (y/x) = {K}")

    z_check = V * y + (1.0 - V) * x
    print(f"  Mass balance error: {np.max(np.abs(z_check - z)):.2e}")

# ---------------------------------------------------------------------------
# 3. P-scan: find the two-phase region at T = 300 K
# ---------------------------------------------------------------------------
print()
print("=" * 60)
print("Pressure scan at 300 K (stability only)")
print("=" * 60)
print(f"{'P [MPa]':>10}  {'Overall':>14}  {'Flags':>10}")
print("-" * 40)

pressures = [0.5e6, 1e6, 2e6, 5e6, 10e6, 20e6, 30e6]
for p_val in pressures:
    mix_s = add_mixture(comps, T=300.0, p=p_val)
    mix_s.mole_fraction = z.copy()
    flags_s, _, _, res = stability_test(mix_s, thermo, opts)
    print(f"{p_val/1e6:>10.1f}  {res['overall']:>14}  {str(flags_s):>10}")
