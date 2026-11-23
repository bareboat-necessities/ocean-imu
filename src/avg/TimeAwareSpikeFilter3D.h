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
 * Magnitude alone cannot see a sample that keeps |v| but points the wrong
 * way (e.g. a sign flip of the vertical axis). When a vertical axis is
 * selected, the same derivative-vs-median rule is also applied to that one
 * component; a vertical spike is replaced by the window vector whose vertical
 * component is the median. Only the chosen axis is checked, so ordinary
 * rotation that moves gravity between the horizontal axes is not flagged.
 *
 * Each sample may carry a companion gyroscope vector captured with the same
 * IMU read. It plays no part in detection, but when a sample is replaced its
 * gyro vector is replaced together with the accelerometer vector by the gyro
 * of the same chosen window sample(s): a corrupted accelerometer read is
 * evidence that the gyro read of that instant is suspect too.
 *
 * Non-finite inputs are passed through unchanged and are not stored, so a
 * single invalid sample cannot poison the window.
 */
class TimeAwareSpikeFilter3D {
  private:
    struct Sample {
      float x, y, z;
      float gx, gy, gz;    // Companion gyro, carried along with the vector
      float magnitude;
      float vertical;
      float magnitudeRate, verticalRate; // Finite differences from the preceding raw read
    };

    Sample *window;      // Circular buffer for previous vectors and time deltas
    int windowSize;      // Size of the moving window
    int currentIndex;    // Current position in the circular buffer
    int quietCount = 0;  // Consecutive slopes strictly below half their thresholds
    float threshold;     // Magnitude-derivative threshold for spike detection
    int verticalAxis;    // 0/1/2 = x/y/z component also checked; -1 = off
    float verticalThreshold; // Vertical-component derivative threshold
    bool initialized;    // Flag to indicate if filter is initialized

    float *temp;
    int *order;

  public:
    /**
       @brief Construct a new 3D Spike Filter object
       @param size Size of the moving window, clamped to at least 2
       @param thr Threshold on the magnitude derivative (signal units per second)
       @param vAxis Component (0 = x, 1 = y, 2 = z) additionally checked on its
              own; any other value disables the vertical check
       @param vThr Threshold on that component's derivative
    */
    TimeAwareSpikeFilter3D(int size, float thr, int vAxis = -1, float vThr = 0.0f)
        : windowSize(size < 2 ? 2 : size), threshold(thr),
          verticalAxis((vAxis >= 0 && vAxis <= 2) ? vAxis : -1), verticalThreshold(vThr),
          initialized(false) {
      window = new Sample[windowSize];
      temp = new float[windowSize];
      order = new int[windowSize];
      currentIndex = 0;
    }

    ~TimeAwareSpikeFilter3D() {
      delete[] order;
      delete[] temp;
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
      float gx = 0.0f, gy = 0.0f, gz = 0.0f;
      return filterWithDelta(x, y, z, gx, gy, gz, deltaT);
    }

    /**
       @brief Filtering of a vector together with its companion gyro sample
       @param x,y,z Raw input vector (detection), overwritten with the output
       @param gx,gy,gz Gyro sample of the same read, overwritten with the gyro
              of the window sample(s) chosen to replace a spike
       @param deltaT Time since last sample
       @return true if the sample was flagged as a spike and replaced
    */
    bool filterWithDelta(float &x, float &y, float &z,
                         float &gx, float &gy, float &gz, float deltaT) {
      if (!isfinite(x) || !isfinite(y) || !isfinite(z) ||
          !isfinite(gx) || !isfinite(gy) || !isfinite(gz) || !(deltaT > 0.0f)) {
        return false;
      }
      const float magnitude = sqrtf(x * x + y * y + z * z);
      const float vertical = verticalAxis == 0 ? x : (verticalAxis == 1 ? y : z);

      if (!initialized) {
        // The initial window repeats the first vector. Preserve the original
        // derivative even if a finite, out-of-range input overflowed its norm.
        const float initialRate = (magnitude - magnitude) / deltaT;
        for (int i = 0; i < windowSize; i++) {
          window[i] = Sample{x, y, z, gx, gy, gz, magnitude, vertical, initialRate, 0.0f};
        }
        quietCount = 0;
        initialized = true;
        return false;
      }

      // Keep the slope with its ending sample. Rejected samples also remain
      // in the raw history used by the next derivative and median.
      const Sample &previous = window[currentIndex ? currentIndex - 1 : windowSize - 1];
      const float magnitudeRate = (magnitude - previous.magnitude) / deltaT;
      const float verticalRate = verticalAxis >= 0 ? (vertical - previous.vertical) / deltaT : 0.0f;
      window[currentIndex] = Sample{x, y, z, gx, gy, gz, magnitude, vertical, magnitudeRate, verticalRate};
      if (++currentIndex == windowSize) currentIndex = 0;

      // If every slope has magnitude below half the threshold, their median
      // does too, so it cannot differ from the newest slope by the threshold.
      if (fabsf(magnitudeRate) < 0.5f * threshold &&
          (verticalAxis < 0 || fabsf(verticalRate) < 0.5f * verticalThreshold)) {
        if (quietCount < windowSize - 1) ++quietCount;
        if (quietCount == windowSize - 1) return false;
      } else {
        quietCount = 0;
      }

      if (isSpike(&Sample::magnitudeRate, magnitudeRate, threshold)) {
        replaceWithMedian(&Sample::magnitude, x, y, z, gx, gy, gz);
        return true;
      }
      if (verticalAxis >= 0 && isSpike(&Sample::verticalRate, verticalRate, verticalThreshold)) {
        replaceWithMedian(&Sample::vertical, x, y, z, gx, gy, gz);
        return true;
      }

      // No spike - keep the original vector
      return false;
    }

  private:
    // Derivative of the newest sample's key versus the median key derivative
    // over the window. A window of N samples yields N - 1 derivatives.
    bool isSpike(float Sample::*rate, float newestRate, float thr) {
      // The slope entering the oldest sample crosses the window boundary.
      // Take the remaining N-1 slopes in chronological order.
      int slot = currentIndex;
      for (int i = 0; i < windowSize - 1; ++i) {
        if (++slot == windowSize) slot = 0;
        temp[i] = window[slot].*rate;
      }
      return fabsf(newestRate - computeMedian(windowSize - 1)) > thr;
    }

    // Replace with the window vector whose key is the median (mean of the two
    // middle vectors for an even window). The companion gyro comes from the
    // same sample(s).
    void replaceWithMedian(float Sample::*key, float &x, float &y, float &z,
                           float &gx, float &gy, float &gz) {
      sortIndicesBy(key);
      const Sample &lo = window[order[(windowSize - 1) / 2]];
      if (windowSize % 2 == 1) {
        x = lo.x;
        y = lo.y;
        z = lo.z;
        gx = lo.gx;
        gy = lo.gy;
        gz = lo.gz;
      } else {
        const Sample &hi = window[order[windowSize / 2]];
        x = 0.5f * (lo.x + hi.x);
        y = 0.5f * (lo.y + hi.y);
        z = 0.5f * (lo.z + hi.z);
        gx = 0.5f * (lo.gx + hi.gx);
        gy = 0.5f * (lo.gy + hi.gy);
        gz = 0.5f * (lo.gz + hi.gz);
      }
    }

    float computeMedian(int n) {
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

    void sortIndicesBy(float Sample::*key) {
      for (int i = 0; i < windowSize; i++) {
        order[i] = i;
      }
      for (int i = 1; i < windowSize; i++) {
        int idx = order[i];
        int j = i - 1;
        while (j >= 0 && window[order[j]].*key > window[idx].*key) {
          order[j + 1] = order[j];
          j--;
        }
        order[j + 1] = idx;
      }
    }
};

#endif
