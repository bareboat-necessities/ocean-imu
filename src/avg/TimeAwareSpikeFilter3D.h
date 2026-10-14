#ifndef TIME_AWARE_SPIKE_FILTER_3D_H
#define TIME_AWARE_SPIKE_FILTER_3D_H

/**
 * Copyright 2025, Mikhail Grushinskiy
 */

#include <math.h>

/**
 * 3D counterpart of TimeAwareSpikeFilter.
 *
 * Spike detection runs on the vector magnitude: the time derivative of
 * |v| for the newest sample is compared against the median |v| derivative
 * over the window. When the deviation exceeds the threshold the sample is
 * replaced by the window vector whose magnitude is the median (the mean of
 * the two middle vectors for an even window), so all three axes are replaced
 * together by an actually observed sample rather than per-axis medians.
 *
 * Non-finite inputs are passed through unchanged and are not stored, so a
 * single invalid sample cannot poison the window.
 */
class TimeAwareSpikeFilter3D {
  private:
    struct Sample {
      float x, y, z;
      float magnitude;
      float deltaTime; // Time since previous sample
    };

    Sample *window;      // Circular buffer for previous vectors and time deltas
    int windowSize;      // Size of the moving window
    int currentIndex;    // Current position in the circular buffer
    float threshold;     // Magnitude-derivative threshold for spike detection
    bool initialized;    // Flag to indicate if filter is initialized

    float *derivatives;
    float *temp;
    int *order;

  public:
    /**
       @brief Construct a new 3D Spike Filter object
       @param size Size of the moving window, clamped to at least 2
       @param thr Threshold on the magnitude derivative (signal units per second)
    */
    TimeAwareSpikeFilter3D(int size, float thr) : windowSize(size < 2 ? 2 : size), threshold(thr), initialized(false) {
      window = new Sample[windowSize];
      derivatives = new float[windowSize];
      temp = new float[windowSize];
      order = new int[windowSize];
      currentIndex = 0;
    }

    ~TimeAwareSpikeFilter3D() {
      delete[] order;
      delete[] temp;
      delete[] derivatives;
      delete[] window;
    }

    // The filter owns raw arrays, so copying it would leave two objects
    // pointing at the same storage and double-free on destruction.
    TimeAwareSpikeFilter3D(const TimeAwareSpikeFilter3D &) = delete;
    TimeAwareSpikeFilter3D &operator=(const TimeAwareSpikeFilter3D &) = delete;

    /**
       @brief Forget the history; the next sample re-initializes the window
    */
    void reset() {
      initialized = false;
      currentIndex = 0;
    }

    /**
       @brief Filtering when you already have the delta time
       @param x,y,z Raw input vector, overwritten with the filtered output
       @param deltaT Time since last sample
       @return true if the sample was flagged as a spike and replaced
    */
    bool filterWithDelta(float &x, float &y, float &z, float deltaT) {
      if (!isfinite(x) || !isfinite(y) || !isfinite(z) || !(deltaT > 0.0f)) {
        return false;
      }
      const float magnitude = sqrtf(x * x + y * y + z * z);

      if (!initialized) {
        // Initialize the window with the first vector
        for (int i = 0; i < windowSize; i++) {
          window[i] = Sample{x, y, z, magnitude, deltaT};
        }
        initialized = true;
        return false;
      }

      // Store the new vector and time delta, then advance the write cursor.
      // After the advance currentIndex points at the oldest sample, which is
      // the slot the next call overwrites.
      int newestIndex = currentIndex;
      window[newestIndex] = Sample{x, y, z, magnitude, deltaT};
      currentIndex = (currentIndex + 1) % windowSize;

      // Magnitude derivatives over consecutive sample pairs. A window of N
      // samples yields N - 1 derivatives.
      int derivativeCount = windowSize - 1;
      for (int i = 0; i < derivativeCount; i++) {
        int currentIdx = (currentIndex + i + 1) % windowSize;
        int prevIdx = (currentIndex + i) % windowSize;

        derivatives[i] = (window[currentIdx].magnitude - window[prevIdx].magnitude) / window[currentIdx].deltaTime;
      }

      // Find the median derivative (more robust than average)
      float medianDerivative = computeMedian(derivatives, derivativeCount);

      int prevIdx = (newestIndex - 1 + windowSize) % windowSize;
      float currentDerivative = (magnitude - window[prevIdx].magnitude) / deltaT;

      if (fabsf(currentDerivative - medianDerivative) > threshold) {
        // Spike detected - replace with the median-magnitude window vector
        sortIndicesByMagnitude();
        const Sample &lo = window[order[(windowSize - 1) / 2]];
        if (windowSize % 2 == 1) {
          x = lo.x;
          y = lo.y;
          z = lo.z;
        } else {
          const Sample &hi = window[order[windowSize / 2]];
          x = 0.5f * (lo.x + hi.x);
          y = 0.5f * (lo.y + hi.y);
          z = 0.5f * (lo.z + hi.z);
        }
        return true;
      }

      // No spike - keep the original vector
      return false;
    }

  private:
    float computeMedian(const float *arr, int n) {
      for (int i = 0; i < n; i++) {
        temp[i] = arr[i];
      }
      // Insertion sort; n is a small window size
      for (int i = 1; i < n; i++) {
        float key = temp[i];
        int j = i - 1;
        while (j >= 0 && temp[j] > key) {
          temp[j + 1] = temp[j];
          j--;
        }
        temp[j + 1] = key;
      }
      if (n % 2 == 1) {
        return temp[n / 2];
      } else {
        return (temp[n / 2 - 1] + temp[n / 2]) / 2.0f;
      }
    }

    void sortIndicesByMagnitude() {
      for (int i = 0; i < windowSize; i++) {
        order[i] = i;
      }
      for (int i = 1; i < windowSize; i++) {
        int key = order[i];
        int j = i - 1;
        while (j >= 0 && window[order[j]].magnitude > window[key].magnitude) {
          order[j + 1] = order[j];
          j--;
        }
        order[j + 1] = key;
      }
    }
};

#endif
