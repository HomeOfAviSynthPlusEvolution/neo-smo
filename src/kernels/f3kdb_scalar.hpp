#pragma once
#include "algorithms/f3kdb.hpp"

namespace neo_smo::f3kdb {
// Plain C++ production backend, also used by Highway SCALAR and EMU128 targets.
Frame process_scalar(const Filter&, std::size_t frame, const std::vector<PlaneView>& input);
} // namespace neo_smo::f3kdb
