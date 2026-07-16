#ifndef _PML_HELPERS_HPP_
#define _PML_HELPERS_HPP_

#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

inline std::vector<double> compute_pml_offsets(int n, double total, const std::string& law, double parameter) {
    std::vector<double> off(n + 1, 0.0);
    if (n <= 1) {
        if (n == 1) off[1] = total;
        return off;
    }
    
    if (law == "power" || law == "powerlaw" || law == "pow") {
        double p = parameter; // Exponent
        for (int l = 1; l <= n; ++l) {
            off[l] = total * std::pow((double)l / n, p);
        }
    } else if (law == "linear" || law == "lin") {
        double R = parameter; // Ratio of last element size to first element size
        if (std::abs(R - 1.0) < 1e-9) {
            for (int l = 1; l <= n; ++l) off[l] = total * l / n;
        } else {
            double h1 = (2.0 * total) / (n * (R + 1.0));
            double d = (2.0 * total * (R - 1.0)) / (n * (n - 1) * (R + 1.0));
            for (int l = 1; l <= n; ++l) {
                off[l] = l * h1 + d * (l * (l - 1)) / 2.0;
            }
        }
    } else if (law == "geom" || law == "geometric" || law == "ratio") {
        double ratio = parameter; // Ratio between consecutive element sizes
        if (std::abs(ratio - 1.0) < 1e-9) {
            for (int l = 1; l <= n; ++l) off[l] = total * l / n;
        } else {
            double rn = std::pow(ratio, n);
            double h1 = total * (ratio - 1.0) / (rn - 1.0);
            double h = h1;
            for (int l = 1; l <= n; ++l) {
                off[l] = off[l - 1] + h;
                h *= ratio;
            }
        }
    } else { // Constant / Uniform (default)
        for (int l = 1; l <= n; ++l) {
            off[l] = total * l / n;
        }
    }
    return off;
}

#endif
