#pragma once
#include "StartupLearner.h"
#include <memory>
namespace startup_learning {
// Pimpl keeps the shipping implementation out of the candidate's translation
// unit. The baseline runs unchanged until its FIRST live sample, then stops.
class ShippingBaseline {
public:
    ShippingBaseline();
    ~ShippingBaseline();
    bool step(const Sample& s);
    Handover handover() const;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
}
