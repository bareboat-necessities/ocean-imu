#pragma once

/*
  Copyright (c) 2026  Mikhail Grushinskiy

  StartupVerticalBiasSeed - vertical accelerometer-bias seed for the OU
  filters' hand-over to Live, from a PII observer run on the startup proxy.

  Why this exists
  ---------------
  The OU filters go Live with a zero accelerometer-bias estimate, and bias
  learning is then held until the magnetic reference has been refined
  (SeaStateFusionFilter_OU_*::setAccBiasHold).  Until the hold releases, a few
  mg of vertical bias sits unexplained in the linear block, and the S
  pseudo-measurement turns it into a metres-sized one-sided displacement
  offset that lasts until the bias estimate converges.

  Before Live the filters already level the accelerometer with a private
  Mahony observer (VerticalAccelComplementary).  A PII observer fed that
  up-positive vertical acceleration rejects a constant input offset through
  its integral channel, so at equilibrium (p, v zero-mean)

      mean(a_up) = ks * S.

  The proxy reports a_up = -((R f)_z + g), and a body bias b_body enters
  f additively, so the world-down component of the bias is

      d = (R b_body)_z = -ks * S.

  That is the one bias component a vertical observer can see, and it is the
  one that drives the heave offset.

  Gains
  -----
  The PII core is a triple real pole at r, so ks * S follows the input mean
  through r^3 / (s + r)^3: it settles in about 7/r and passes wave
  acceleration at about (r/omega)^3.  The shipping PII schedules r near
  0.09 rad/s for heave, which settles in ~75 s -- longer than the ~30 s the
  filters spend before Live.  This observer therefore runs its own fixed r,
  below the 0.6-2.6 rad/s wave band but fast enough to settle before
  go-live, and reports a trailing one-pole average of -ks * S so that the
  wave-frequency ripple that does leak through is averaged out.  It is not
  adapted and has no bias-trend channel: its only job is this one mean.

  It starts after the proxy's own levelling transient, and the seed is only a
  mean: it sets the bias state at go-live and leaves the covariance, the
  hold, and every other state exactly as the deployed hand-over leaves them.

  Nothing here reads filter state; it is a pure function of the proxy's
  output, so the startup interconnection stays open.
*/

#include <cmath>

#include "pii_observer/VerticalPIIObserver.h"

class StartupVerticalBiasSeed {
public:
    using Observer = marine_obs::VerticalPIIObserver<float, false>;

    static constexpr float POLE_RATE_DEFAULT = 0.30f;     // r [rad/s]
    static constexpr float START_DELAY_SEC = 5.0f;        // proxy levelling transient
    static constexpr float AVERAGE_TAU_SEC = 8.0f;        // trailing average of -ks*S
    // The observer must have settled (~7/r) before the seed is trusted.
    static constexpr float MIN_RUN_SEC = 20.0f;
    // Larger than any healthy MEMS turn-on bias; a bigger value means the
    // observer has not settled, and the deployed zero seed is kept.
    static constexpr float MAX_ABS_MPS2 = 0.3f;

    StartupVerticalBiasSeed() { configure_(); }

    void setEnabled(bool on) noexcept { enabled_ = on; }
    bool enabled() const noexcept { return enabled_; }

    void setPoleRate(float r) {
        if (std::isfinite(r) && r > 0.0f) {
            pole_rate_ = r;
            configure_();
        }
    }

    void reset() {
        pii_.reset();
        elapsed_sec_ = 0.0f;
        run_sec_ = 0.0f;
        avg_ = 0.0f;
    }

    void update(float a_up_mps2, float dt) {
        if (!enabled_) return;
        if (!(dt > 0.0f) || !std::isfinite(dt) || !std::isfinite(a_up_mps2)) return;
        elapsed_sec_ += dt;
        if (elapsed_sec_ < START_DELAY_SEC) return;
        pii_.update(a_up_mps2, dt);
        run_sec_ += dt;
        const float d = -pii_.ks() * pii_.integral_displacement();
        avg_ += (dt / (AVERAGE_TAU_SEC + dt)) * (d - avg_);
    }

    // World-down accelerometer bias [m/s^2], NaN until trusted.
    float downBiasMps2() const {
        if (!enabled_ || run_sec_ < MIN_RUN_SEC) return NAN;
        if (!std::isfinite(avg_) || std::fabs(avg_) > MAX_ABS_MPS2) return NAN;
        return avg_;
    }

private:
    void configure_() {
        Observer::Config cfg{};
        cfg.r = pole_rate_;
        pii_.configure(cfg);
        pii_.set_adaptation_enabled(false);
        pii_.reset();
    }

    bool enabled_ = false;
    float pole_rate_ = POLE_RATE_DEFAULT;
    float elapsed_sec_ = 0.0f;
    float run_sec_ = 0.0f;
    float avg_ = 0.0f;
    Observer pii_{};
};
