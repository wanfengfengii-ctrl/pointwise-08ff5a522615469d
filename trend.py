import math
from matrix import inverse, multiply

def fit(total):
    covariance = inverse(total["information"])
    coefficients = multiply(covariance, total["rhs"])
    intercept, slope = coefficients[0][0], coefficients[1][0]
    variance = covariance[1][1]
    if variance <= 0:
        raise ArithmeticError("nonpositive slope variance")
    se = math.sqrt(float(variance))
    return {
        "intercept": intercept, "slope": slope, "slope_variance": variance,
        "standard_error": se,
        "ci95": [float(slope) - 1.96*se, float(slope) + 1.96*se],
        "information_rank": total["information_rank"],
        "window_count": total["window_count"]
    }
