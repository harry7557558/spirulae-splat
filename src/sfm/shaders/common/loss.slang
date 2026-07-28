#pragma once
// Robust loss functions. s is the squared residual norm; following Ceres,
// the total cost is sum(0.5 * rho(s)) and residual/Jacobian are scaled by
// sqrt(rho'(s)) (Triggs correction omitted -- fine for LM line quality).
//
// The loss is selected at shader-compile time with -DLOSS_*; the scale
// parameter (delta for Huber, c for Cauchy) comes in at runtime via push
// constant so it can be set per-problem without recompiling.

#include "real.slang"

interface ILoss {
    static Real cost(Real s, float param);
    static Real weight(Real s, float param); // rho'(s)
};

struct TrivialLoss : ILoss {
    static Real cost(Real s, float param) { return s; }
    static Real weight(Real s, float param) { return Real(1.0); }
};

struct HuberLoss : ILoss {
    // rho(s) = s                     if s <= d^2
    //          2 d sqrt(s) - d^2     otherwise
    static Real cost(Real s, float param) {
        Real d2 = Real(param) * Real(param);
        if (s <= d2) return s;
        return Real(2.0) * Real(param) * sqrt(s) - d2;
    }
    static Real weight(Real s, float param) {
        Real d2 = Real(param) * Real(param);
        if (s <= d2) return Real(1.0);
        return Real(param) / sqrt(s);
    }
};

struct CauchyLoss : ILoss {
    // rho(s) = c^2 log(1 + s/c^2);  rho'(s) = 1 / (1 + s/c^2)
    static Real weight(Real s, float param) {
        Real c2 = Real(param) * Real(param);
        return Real(1.0) / (Real(1.0) + s / c2);
    }
    static Real cost(Real s, float param) {
        Real c2 = Real(param) * Real(param);
        return c2 * rlog(Real(1.0) + s / c2);
    }
};

#if defined(LOSS_HUBER)
typealias Loss = HuberLoss;
#elif defined(LOSS_CAUCHY)
typealias Loss = CauchyLoss;
#else
typealias Loss = TrivialLoss;
#endif
