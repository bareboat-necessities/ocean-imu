#include "../common/FiniteMean.h"
#include <stdexcept>
#include <iostream>
int main() {
    FiniteMean m;
    if (!std::isnan(m.mean())) throw std::runtime_error("missing mean must be NaN");
    m.add(8);
    m.add(std::numeric_limits<double>::quiet_NaN());
    m.add(std::numeric_limits<double>::infinity());
    if (m.count != 1 || m.mean() != 8) throw std::runtime_error("invalid NIS diluted mean");
    m.add(4);
    if (m.count != 2 || m.mean() != 6) throw std::runtime_error("finite NIS mean mismatch");
    std::cout << "finite-only diagnostic means PASS\n";
}
