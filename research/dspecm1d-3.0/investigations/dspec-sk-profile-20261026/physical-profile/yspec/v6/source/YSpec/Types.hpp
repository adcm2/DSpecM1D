#pragma once

/**
 * @file Types.hpp
 * @brief The index and complex types used throughout.
 */

#include <complex>
#include <cstddef>

namespace YSpec {

/// The signed index type, as in GSHTrans.
using Int = std::ptrdiff_t;

/// The complex type of spectra and of the radial solutions.
using Complex = std::complex<double>;

}  // namespace YSpec
