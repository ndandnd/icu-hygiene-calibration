"""Closed forms for the coefficient-extension study (no simulation, no frozen imports).

* Transition-defined (WHO-style) compliance. With per-boundary completion
  probability h and Markov persistence rho (stationary law (1-h, h); a miss is
  followed by a miss with probability rho + (1-rho)(1-h)), a zone transition
  (exit then entry) is missed only if both boundaries are missed:

      T(h, rho) = 1 - (1-h) (rho + (1-rho)(1-h)).

  For a target T = c, m = 1-h solves (1-rho) m^2 + rho m - (1-c) = 0 on (0, 1).
  The root is computed in 60-digit decimal arithmetic with the cancellation-free
  form m = 2(1-c) / (rho + sqrt(rho^2 + 4(1-rho)(1-c))) and then rounded once to
  the nearest double.

* Isolated k-opportunity bridge and the renewal comparator F_geo of supplement
  S4.2 ("The extended surface"): ``Pm``, ``Dm``, ``qk``, ``Fk``, ``estar``,
  ``sum_even``, ``Fgeo`` and ``Fgeo_vec`` are ported verbatim from
  ``manuscript_v2_2026-10-21/reviews_round1/analysis/round1_analysis.py``
  (section 4 and section 9), which cannot be imported because it runs every
  analysis at import time.

      S(h) = pi Gamma P Gamma {I - (P Gamma)^2}^{-1} 1,   F_geo = 1 - S(h1)/S(h0).
"""
import sys

sys.dont_write_bytecode = True

import decimal  # noqa: E402
from decimal import Decimal  # noqa: E402
import math  # noqa: E402

import numpy as np  # noqa: E402

DECIMAL_PRECISION = 60


# ------------------------------------------------------------ transition policy
def transition_compliance(h, rho):
    """WHO-style transition compliance: probability that exit or entry is performed."""
    return 1 - (1 - h) * (rho + (1 - rho) * (1 - h))


def transition_h_decimal(rho, target):
    """Exact per-boundary h (Decimal) with transition compliance ``target``.

    ``rho`` and ``target`` are read through ``str`` so that .5 and .90 are the
    exact decimals of the protocol, not their binary approximations.
    """
    with decimal.localcontext() as context:
        context.prec = DECIMAL_PRECISION
        r = Decimal(str(rho))
        miss = 1 - Decimal(str(target))
        if not (0 <= r < 1 and 0 < miss < 1):
            raise ValueError('rho must be in [0,1) and target in (0,1)')
        m = 2 * miss / (r + (r * r + 4 * (1 - r) * miss).sqrt())
        residual = (1 - r) * m * m + r * m - miss
        if abs(residual) > Decimal(10) ** (-(DECIMAL_PRECISION - 5)):
            raise ArithmeticError('transition root residual too large')
        return +(1 - m)


def transition_h(rho, target):
    """Nearest double to the exact root (``float`` of a Decimal is correctly rounded)."""
    return float(transition_h_decimal(rho, target))


def transition_root_record(rho, target):
    exact = transition_h_decimal(rho, target)
    value = float(exact)
    with decimal.localcontext() as context:
        context.prec = DECIMAL_PRECISION
        miss_exact = 1 - Decimal(str(target))
        miss_at_double = (1 - Decimal(value)) * (Decimal(str(rho)) + (1 - Decimal(str(rho))) * (1 - Decimal(value)))
        error = float(miss_at_double - miss_exact)
    return dict(rho=float(rho), target_transition_compliance=str(target),
                h_exact_decimal_40_digits=format(exact, '.40f'), h_double=value, h_double_repr=repr(value),
                h_rounded_6=round(value, 6),
                transition_compliance_at_double_float=transition_compliance(value, float(rho)),
                transition_miss_error_at_double_exact_arithmetic=error)


# ------------------------------------------------- round1_analysis.py, section 4
def Pm(h, rho):
    return np.array([[rho + (1 - rho) * (1 - h), (1 - rho) * h],
                     [(1 - rho) * (1 - h), rho + (1 - rho) * h]])


def Dm(e):
    return np.diag([1.0, 1.0 - e])


def qk(k, h, rho, e):
    pi = np.array([1 - h, h])
    M = Pm(h, rho) @ Dm(e)
    v = pi @ Dm(e)
    for _ in range(k - 1):
        v = v @ M
    return float(v.sum())


def Fk(k, h0, h1, rho, e):
    return 1 - qk(k, h1, rho, e) / qk(k, h0, rho, e)


def estar(h0, h1):
    s0 = math.sqrt(h0 * (1 - h0))
    s1 = math.sqrt(h1 * (1 - h1))
    return (s0 - s1) / (h1 * s0 - h0 * s1)


def sum_even(h, rho, e):
    """Renewal comparator S(h): one susceptible contact every two opportunities, unbounded stationary episode."""
    pi = np.array([1 - h, h])
    M = Pm(h, rho) @ Dm(e)
    S = pi @ Dm(e) @ Pm(h, rho) @ Dm(e) @ np.linalg.inv(np.eye(2) - M @ M)
    return float(S.sum())


def Fgeo(h0, h1, rho, e):
    """F_geo = 1 - S(h1)/S(h0) (supplement S4.2)."""
    return 1 - sum_even(h1, rho, e) / sum_even(h0, rho, e)


# ------------------------------------------------- round1_analysis.py, section 9
def Fgeo_vec(h0, h1, rho, e):
    rho = np.asarray(rho, float)
    e = np.asarray(e, float)

    def S(h):
        # closed form of pi D P D (I - (PD)^2)^{-1} 1 evaluated elementwise
        a11 = rho + (1 - rho) * (1 - h)
        a12 = (1 - rho) * h * (1 - e)
        a21 = (1 - rho) * (1 - h)
        a22 = (rho + (1 - rho) * h) * (1 - e)
        # M = P D = [[a11, a12],[a21, a22]]; M2 = M @ M
        m11 = a11 * a11 + a12 * a21
        m12 = a11 * a12 + a12 * a22
        m21 = a21 * a11 + a22 * a21
        m22 = a21 * a12 + a22 * a22
        # (I - M2)^{-1} 1
        A11, A12, A21, A22 = 1 - m11, -m12, -m21, 1 - m22
        det = A11 * A22 - A12 * A21
        x1 = (A22 - A12) / det
        x2 = (-A21 + A11) / det
        # v = pi D M
        v1 = (1 - h) * a11 + h * (1 - e) * a21
        v2 = (1 - h) * a12 + h * (1 - e) * a22
        return v1 * x1 + v2 * x2
    return 1 - S(h1) / S(h0)


def q2_response(h0, h1, rho, e):
    """Isolated two-opportunity response 1 - q2(h1)/q2(h0) (same as Fk(2, ...))."""
    return Fk(2, h0, h1, rho, e)


__all__ = ['Dm', 'Fgeo', 'Fgeo_vec', 'Fk', 'Pm', 'estar', 'q2_response', 'qk', 'sum_even',
           'transition_compliance', 'transition_h', 'transition_h_decimal', 'transition_root_record']
