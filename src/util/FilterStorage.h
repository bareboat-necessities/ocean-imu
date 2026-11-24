#pragma once

#include <algorithm>
#include <memory>
#include <new>
#include <utility>

#if defined(ARDUINO_ARCH_ESP32) && defined(CONFIG_IDF_TARGET_ESP32S3) && defined(BOARD_HAS_PSRAM)
#include <esp_heap_caps.h>
#endif

namespace ocean_imu {

// The covariance matrix and its scratch space share the estimator allocation.
// Prefer internal SRAM on PSRAM-equipped S3 boards for the 200 Hz read/write path.
#if defined(ARDUINO_ARCH_ESP32) && defined(CONFIG_IDF_TARGET_ESP32S3) && defined(BOARD_HAS_PSRAM)
template<class Filter>
struct FilterDeleter {
    bool internal = false;
    void operator()(Filter* ptr) const {
        if (internal) {
            ptr->~Filter();
            heap_caps_free(ptr);
        } else {
            delete ptr;
        }
    }
};

template<class Filter>
using FilterPtr = std::unique_ptr<Filter, FilterDeleter<Filter>>;

template<class Filter, class... Args>
FilterPtr<Filter> makeFilter(Args&&... args) {
    std::unique_ptr<void, decltype(&heap_caps_free)> block(
        heap_caps_aligned_alloc(std::max(alignof(Filter), alignof(void*)), sizeof(Filter),
                                MALLOC_CAP_INTERNAL | MALLOC_CAP_8BIT), heap_caps_free);
    if (!block) return FilterPtr<Filter>(new Filter(std::forward<Args>(args)...));
    auto* ptr = ::new (block.get()) Filter(std::forward<Args>(args)...);
    block.release();
    return FilterPtr<Filter>(ptr, FilterDeleter<Filter>{true});
}
#else
template<class Filter>
using FilterPtr = std::unique_ptr<Filter>;

template<class Filter, class... Args>
FilterPtr<Filter> makeFilter(Args&&... args) {
    return std::make_unique<Filter>(std::forward<Args>(args)...);
}
#endif

}  // namespace ocean_imu
