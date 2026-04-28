"""Example 2: Phase envelope — methane + propane binary.

Demonstrates:
    - Tracing the bubble-point and dew-point pressure curves as a function
      of temperature for a binary hydrocarbon mixture
    - Using bubble_pressure and dew_pressure with a temperature loop
    - How the two curves meet at the cricondentherm (maximum temperature)
      and cricondenbar (maximum pressure)

The output is a table of (T, P_bub, P_dew) values that describe the
two-phase envelope.  If matplotlib is installed a P-T plot is also produced.
"""

import numpy as np

from pvttool import (
    FlashOptions,
    add_components,
    add_mixture,
    add_thermo,
    bubble_pressure,
    dew_pressure,
)

# ---------------------------------------------------------------------------
# System: 30 mol% methane + 70 mol% propane
# ---------------------------------------------------------------------------
comps = add_components(["Methane", "Propane"])
z = np.array([0.30, 0.70])
thermo = add_thermo()
opts = FlashOptions(accuracy=1e-6, iteration=150)

T_range = np.arange(200.0, 361.0, 10.0)   # K — from 200 K up to near cricondentherm

print("=" * 60)
print("Phase Envelope: 30% CH4 / 70% C3H8  (PR EOS, vdW mixing)")
print("=" * 60)
print(f"{'T [K]':>8}  {'P_bub [bar]':>12}  {'P_dew [bar]':>12}  {'Bub OK':>6}  {'Dew OK':>6}")
print("-" * 60)

T_bub_list, P_bub_list = [], []
T_dew_list, P_dew_list = [], []

for T in T_range:
    mix = add_mixture(comps, T=T, p=1e6)  # P is overridden by the solver
    mix.mole_fraction = z.copy()

    P_bub, y_bub, ok_bub = bubble_pressure(mix, thermo, opts)
    P_dew, x_dew, ok_dew = dew_pressure(mix, thermo, opts)

    p_bub_bar = P_bub / 1e5
    p_dew_bar = P_dew / 1e5
    print(f"{T:>8.1f}  {p_bub_bar:>12.3f}  {p_dew_bar:>12.3f}  {'yes' if ok_bub else 'no':>6}  {'yes' if ok_dew else 'no':>6}")

    if ok_bub:
        T_bub_list.append(T)
        P_bub_list.append(P_bub / 1e5)
    if ok_dew:
        T_dew_list.append(T)
        P_dew_list.append(P_dew / 1e5)

# ---------------------------------------------------------------------------
# Cricondenbar estimate (maximum bubble pressure)
# ---------------------------------------------------------------------------
if P_bub_list:
    idx = int(np.argmax(P_bub_list))
    print()
    print(f"Cricondenbar (approx): T ≈ {T_bub_list[idx]:.0f} K, "
          f"P ≈ {P_bub_list[idx]:.2f} bar")

# ---------------------------------------------------------------------------
# Optional matplotlib plot
# ---------------------------------------------------------------------------
try:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(T_bub_list, P_bub_list, "b-o", markersize=4, label="Bubble curve")
    ax.plot(T_dew_list, P_dew_list, "r--s", markersize=4, label="Dew curve")
    ax.set_xlabel("Temperature [K]")
    ax.set_ylabel("Pressure [bar]")
    ax.set_title("Phase Envelope: 30% CH₄ / 70% C₃H₈  (PR EOS)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("phase_envelope_ch4_c3h8.png", dpi=150)
    print("\nPlot saved to phase_envelope_ch4_c3h8.png")
    plt.show()
except ImportError:
    print("\n(Install matplotlib to generate the P-T plot)")
