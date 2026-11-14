#pragma once

#include <iomanip>
#include <limits>
#include <ostream>

namespace ocean_imu::sim {

// Preserve the input clock independently of the precision of sensor columns.
inline std::ostream& write_csv_time(std::ostream& stream, double time)
{
    const auto flags = stream.flags();
    const auto precision = stream.precision();
    stream << std::defaultfloat
           << std::setprecision(std::numeric_limits<double>::max_digits10) << time;
    stream.flags(flags);
    stream.precision(precision);
    return stream;
}

} // namespace ocean_imu::sim
