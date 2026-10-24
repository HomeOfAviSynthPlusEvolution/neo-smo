#pragma once
#include "kernels/f3kdb.hpp"

namespace neo_smo::f3kdb {
// Plain C++ production backend, also used by Highway SCALAR and EMU128 targets.
Frame process_scalar(const Filter&, std::size_t frame, const std::vector<PlaneView>& input);
void process_row_scalar(const Row&, int mode, const double* thresholds, bool integer, bool blur, double angle,
                        double boost);
void weights_scalar(const double* responses[4], double* output, std::size_t count);
bool unpack_plane_scalar(const InputPlane&);
void pack_row_scalar(const OutputRow&);
void sobel_plane_scalar(const double* source, double* gx, double* gy, std::size_t width, std::size_t height);
void diffuse_stripe_scalar(const FsStripe&);
} // namespace neo_smo::f3kdb
