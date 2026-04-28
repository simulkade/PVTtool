"""Universal Quasi-Chemical (UNIQUAC) activity coefficient model."""

from __future__ import annotations

import numpy as np


_R = 8.314
_z = 10  # coordination number


def uniquac(
    T: float,
    x: np.ndarray,
    components: list,
    bip,
) -> tuple[float, np.ndarray]:
    """UNIQUAC activity coefficient model.

    Parameters
    ----------
    T : float
        Temperature [K].
    x : np.ndarray, shape (n,)
        Mole fractions.
    components : list[Component]
        Pure-component objects (uses ``uniquacR`` and ``uniquacQ`` attributes).
    bip : BIP
        Binary interaction parameters.  Uses ``uniquac_cons``, ``uniquac_tdep``,
        ``uniquac_r`` and ``uniquac_q``.

    Returns
    -------
    g_ert : float
        Excess Gibbs energy divided by RT [dimensionless].
    gamma : np.ndarray, shape (n,)
        Activity coefficients.
    """
    N = len(x)
    r = bip.uniquac_r if np.any(bip.uniquac_r != 0) else np.array([c.uniquacR for c in components])
    q = bip.uniquac_q if np.any(bip.uniquac_q != 0) else np.array([c.uniquacQ for c in components])

    A = (bip.uniquac_cons + bip.uniquac_tdep * T)
    np.fill_diagonal(A, 0.0)
    taw = np.exp(-A / (_R * T))
    np.fill_diagonal(taw, 0.0)

    xr = np.dot(x, r)
    xq = np.dot(x, q)
    fay = x * r / xr
    teta = x * q / xq
    l = _z / 2.0 * (r - q) - (r - 1.0)

    if N == 2:
        x1, x2 = x[0], x[1]
        f1, f2 = fay[0], fay[1]
        t1, t2 = teta[0], teta[1]
        q1, q2_c = q[0], q[1]
        t21, t12 = taw[1, 0], taw[0, 1]

        gEcrt = (x1 * np.log(f1 / x1) + x2 * np.log(f2 / x2)
                 + _z / 2.0 * (q1 * x1 * np.log(t1 / f1) + q2_c * x2 * np.log(t2 / f2)))
        gErrt = (-q1 * x1 * np.log(t1 + t2 * t21)
                 - q2_c * x2 * np.log(t2 + t1 * t12))
        g_ert = float(gEcrt + gErrt)

        lngama = np.zeros(2)
        for i in range(2):
            j = 1 - i
            lngama[i] = (
                np.log(fay[i] / x[i])
                + _z / 2.0 * q[i] * np.log(teta[i] / fay[i])
                + fay[j] * (l[i] - r[i] / r[j] * l[j])
                - q[i] * np.log(teta[i] + teta[j] * taw[j, i])
                + teta[j] * q[i] * (taw[j, i] / (teta[i] + teta[j] * taw[j, i])
                                     - taw[i, j] / (teta[j] + teta[i] * taw[i, j]))
            )
        return g_ert, np.exp(lngama)

    # N > 2: vectorized
    g_ert = float(
        np.dot(x, np.log(fay / x))
        + _z / 2.0 * np.dot(x, q * np.log(teta / fay))
        - np.dot(x, q * np.log(teta @ taw.T))
    )
    xl = np.dot(x, l)
    teta_taw = teta @ taw.T        # shape (n,)
    lngama = (np.log(fay / x)
              + _z / 2.0 * q * np.log(teta / fay)
              + l - fay / x * xl
              - q * teta_taw
              + q
              - q * ((1.0 / teta_taw) * (teta @ taw.T)))
    return g_ert, np.exp(lngama)
