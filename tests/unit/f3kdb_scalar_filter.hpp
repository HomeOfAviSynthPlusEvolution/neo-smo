#pragma once
#include "kernels/f3kdb_scalar.hpp"

// Bypass dynamic dispatch so the differential baseline and scalar tests always
// execute the production C++ backend, even on a machine with hardware SIMD.
class ScalarFilter : public neo_smo::f3kdb::Filter {
public:
  using Filter::Filter;
  neo_smo::f3kdb::Frame process(std::size_t frame, const std::vector<neo_smo::f3kdb::PlaneView>& input) const {
    return neo_smo::f3kdb::process_scalar(*this, frame, input);
  }
};
