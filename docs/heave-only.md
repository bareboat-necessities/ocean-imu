# Heave-only filters

`tests/heave_only/` implements three published heave-only filters as
comparison baselines (outside `src/`, which holds only shipping code). Each
estimates vertical displacement from the levelled vertical acceleration
alone, and all are compared with the OU filters on the same replay:

- `GodhavnHeaveFilter.h`: the standard adaptive heave filter of Godhavn
  (OCEANS'98). It is implemented as specified by Richter et al. (IFAC 2014),
  Sections 2.2–2.4, because the original paper is not openly available.
- `RichterHeaveFilter.h`: the zero-displacement heave filter of Richter et
  al. (IFAC 2014), Section 3.2.
- `KuchlerHeaveEKF.h`: the harmonic-mode EKF of Küchler et al. (IFAC 2011).

Two support headers:

- `AdaptiveHeaveFilter.h`: the engine the Godhavn and Richter filters share.
  It holds the cascaded sections, the identification of `w_p` and `A_p`, and
  cutoff retuning. Each filter supplies its own structure law.
- `AccelSpectrum.h`: the FFT identification front end.

The filters take the levelled vertical acceleration (z up, gravity removed)
and return heave. In `tests/heave_only` they share one attitude front
end with the PII observer, which runs alongside them. Equation numbers below
are those of the cited papers.

## Godhavn 1998 (as specified by Richter et al. 2014)

Heave is the vertical acceleration through (Eq. 3)

    H(s) = s^2 / (s^2 + 2 zeta w_c s + w_c^2)^2,   zeta = 1/sqrt(2)

This is a double integrator in series with two second-order Butterworth high
passes. It is realised as two identical sections, integrated with RK4. The
double zero rejects a constant bias. The cost is a phase lead,
`|1 - s^2 H(j w_p)| -> 2 sqrt(2) w_c / w_p` (Eq. 7).

The cutoff adapts as in Fig. 4. Every 30 s, the FFT of the buffered
acceleration (128 s at 4 Hz) gives the dominant heave frequency `w_p`. The
buffer also gives the mean wave height (Eq. 14):

    A_p = sqrt(2/N sum a_j^2 / w_p^4)

The cutoff then minimises the total-error bound (Eq. 11):

    J = 8 A_p^2 (w_c/w_p)^2 + sigma_n^2 / (2^(7/2) w_c^3)

The closed-form minimiser is (Eq. 12):

    w_c,opt = 2^(-3/2) (3 sigma_n^2 w_p^2 / A_p^2)^(1/5)

`sigma_n^2` is the datasheet white-noise density. The result is low-pass
filtered with a 60 s relaxation in `log w_c`.

As published, this law accounts only for white sensor noise. On the replay's
MEMS accelerometer it fails:
- It selects `w_c` = 0.010–0.03 rad/s in all but the smallest seas.
- The 5e-4 m/s²/√s bias random walk then passes through the filter's low
  frequency gain as a slow heave drift of metres.

Richter et al. validated the law on a test bench, where bias drift was
presumably negligible.

`godhavn_ext` adds two extensions that are not in either paper:
- **Bias-instability term.** `J` gains the variance
  `q_b / (2^(7/2) w_c^5)`, which a bias random walk of intensity `q_b` leaves
  after `H` (exact for this `H`). The extended `J` is minimised numerically.
- **Input-mean subtraction.** A 300 s running mean is removed ahead of `H`,
  so retuning does not move its DC state.

With both, the adapted cutoff runs from 0.17 rad/s (Hs 0.27 m) down to
0.032 rad/s (Hs 8.5 m).

## Richter et al. 2014: zero-displacement filter (Section 3.2)

Richter et al. identify the standard filter's phase lead as its main
real-time error, and propose moving one of its zeros off the origin (Eq. 23):

    Hzd(s) = s (s + a) / D(s)^2,   a = 2 sqrt(2) w_c (1 - w_c^2 / w_p^2)   (Eq. 25)

- **Error.** The error on a tone falls from `2 sqrt(2) w_c/w_p` to
  `4 (w_c/w_p)^2` (Eq. 26).
- **Cost.** It amplifies white noise about 9 times more (Eq. 27).
- **Cutoff.** The cutoff minimises Eq. 28, which gives
  `w_c^7 = (27 sqrt(2)/1024) sigma_n^2 w_p^4 / A_p^2` (Eq. 29).
- **Implementation.** `s/D^2` is the second section's position state, so the
  output is the standard heave plus `a` times that state.
- **Drift sensitivity.** With a single zero at the origin, Hzd passes far
  more bias drift than H. Its bias-instability variance is
  `q_b sqrt(2)/16 (w_c^-5 + 3 a^2 w_c^-7)`, exact for Hzd.

`richter_zd` is this filter with the same two extensions as `godhavn_ext`.

On their test bench, Richter et al. report this filter at 59.6 % of the
standard filter's RMS error (Table 1). Their pole-zero filter reaches 43.4 %.
Here, with a good vertical reference, it reaches 60 % (5.11 % against
8.45 % of Hs).

## Küchler et al. 2011

Heave is a sum of undamped modes with `x_j = [z_j, z_j', w_j]` and
`z_j'' = -w_j^2 z_j` (Eq. 5), where `w_j` is a random walk. A random-walk
offset state absorbs gravity and bias (Eq. 8). Its initial value is the
paper's `-g`, which is 0 after gravity removal. The measurement is the sum of
the mode accelerations plus the offset (Eq. 9). Each mode is discretised
exactly (Eq. 6).

The steps are:

- **Identification.** An FFT of the buffered acceleration gives the heave
  amplitude spectrum `|A_acc|/w^2` (Eq. 2). Its peaks set the modes. The EKF
  starts at the first identification.
- **Re-identification.** This runs every 30 s. A new peak extends the model,
  and only its states are initialised. A mode whose frequency approaches zero
  is removed. Converging modes are merged (Eq. 11).
- **Noise.** Process noise is larger for higher-frequency modes. `R` is the
  datasheet sensor noise.

Choices the paper leaves open:

- **Seeding.** New modes are seeded with amplitude and phase from a joint
  least-squares fit of the buffer, instead of reading the FFT bins.
- **Mode limits.** "Approaches zero" means below 0.04 Hz. Modes above 1 Hz
  are also removed. At most 4 modes are kept, and a peak twice as strong as
  the weakest mode replaces it.
- **Process noise.** Velocity process noise is
  `q_j = 4 zeta_q w_j^3 (A_j^2/2)`, with `zeta_q = 0.01` and a frequency
  random walk of `0.002 w_j/sqrt(s)`. These were picked by a sweep on this
  replay (mean 12.1–17.4 % of Hs over the range tried). Using 8 modes did not
  help.
- **Input.** The paper uses the raw body-z accelerometer and neglects roll
  and pitch. `HO_KUCHLER_BODY_Z=1` reproduces that: the mean error is
  13.7 %, against 12.6 % with levelled input. `HO_KUCHLER_OFFSET_FIT=1`
  seeds the offset from the buffer fit instead (12.4 %).

## Comparison on the OU replay

`tests/heave_only/heave_only-sim` runs every method through
`process_wave_file_for_tracker`, the noise template and seeds of the OU-II,
OU-III and TFG simulators. It uses:

- the 28 ft RAO dataset `v1.2.3`;
- an 8 mg accelerometer turn-on bias and a 5e-4 m/s²/√s bias random walk;
- a 25 Hz magnetometer.

It scores vertical displacement over the same trailing 900 s window. The
heave-only filters need a vertical reference (`--frontend`):

- **`mahony` (default).** The PII observer's Mahony AHRS as shipped: its
  adaptive gains, the heading-only magnetometer correction
  (`ahrs/MahonyHeadingMag.h`), and the firmware's sensor mapping, with
  gyro, accelerometer and magnetometer all mapped `[N, -E, -D]`.
- **`proxy`.** The same Mahony core, IMU-only, at the gains of the private
  Mahony observer that OU-II, OU-III and TFG share
  (`defaults::STARTUP_PROXY_TWO_KP = 0.2`, `STARTUP_PROXY_TWO_KI = 0.02`).
  No new tuning is introduced.
- **`truth`.** The measured accelerometer, still noisy and biased, levelled
  with the record's true attitude. This is an ideal VRU. PII is left out of
  this mode, because it carries its own attitude estimator.

Z RMS error, % of Hs, with the `proxy` front end:

| Record | OU-III | OU-II | TFG | PII | Richter zero-displacement | Godhavn ext | Küchler | Godhavn published law |
|---|---|---|---|---|---|---|---|---|
| JONSWAP 0.27 m | 4.25 | 6.02 | 4.38 | 3.71 | 4.58 | 6.83 | 4.79 | 22.05 |
| JONSWAP 1.5 m | 4.10 | 6.63 | 4.19 | 5.22 | 5.15 | 8.06 | 10.81 | 94.40 |
| JONSWAP 4 m | 4.02 | 6.35 | 4.03 | 5.81 | 5.03 | 8.58 | 13.61 | 223 |
| JONSWAP 8.5 m | 3.71 | 6.05 | 3.73 | 5.99 | 5.31 | 9.77 | 15.85 | 387 |
| PM-Stokes 0.27 m | 4.21 | 5.81 | 4.33 | 3.81 | 5.00 | 7.28 | 5.27 | 28.98 |
| PM-Stokes 1.5 m | 4.02 | 6.43 | 4.09 | 4.75 | 5.33 | 8.28 | 12.33 | 162 |
| PM-Stokes 4 m | 3.96 | 6.20 | 3.98 | 5.49 | 4.74 | 8.16 | 15.98 | 313 |
| PM-Stokes 8.5 m | 3.70 | 5.98 | 3.73 | 5.87 | 5.75 | 10.64 | 18.59 | 507 |
| **Mean** | **4.00** | **6.18** | **4.06** | **5.08** | **5.11** | **8.45** | **12.15** | **217** |
| Worst | 4.25 | 6.63 | 4.38 | 5.99 | 5.75 | 10.64 | 18.59 | 507 |

How each method changes with the front end. Tilt error is the angle between
the estimated and true gravity directions:

| Mean (worst) % of Hs | Tilt error, mean (worst) | PII | Richter zero-displacement | Godhavn ext | Küchler | Godhavn published law |
|---|---|---|---|---|---|---|
| `mahony`: PII Mahony as shipped | 1.53° (3.36°) | 5.64 (7.81) | 9.08 (18.46) | 11.88 (21.99) | 12.51 (19.88) | 203 (422) |
| `proxy`: shared IMU-only Mahony | 0.58° (1.18°) | 5.08 (5.99) | 5.11 (5.75) | 8.45 (10.64) | 12.15 (18.59) | 217 (507) |
| `truth`: true vertical | 0.00° (0.00°) | – | 5.00 (5.33) | 8.23 (9.65) | 12.61 (20.12) | 215 (492) |

### Mahony settings across the repository

Every Mahony instance here runs the same `ahrs/Mahony_AHRS.h` core. What
differs is the gains, whether the magnetometer is used, and what the
attitude feeds. The correction corner is about `2Kp/2` rad/s.

| User | 2Kp | 2Ki | Corner | Magnetometer | Role |
|---|---|---|---|---|---|
| PII observer (`AdaptiveVerticalPIIMahony`) | 1.2–1.7, adapted to sea state | 0.009–0.0125 | 0.6–0.85 rad/s (0.10–0.14 Hz) | yes | the attitude its heave channel is levelled with |
| OU-II, OU-III, TFG (`VerticalAccelComplementary`, `SeaStateFusionDefaults.h`) | 0.2 | 0.02 | 0.1 rad/s (0.016 Hz) | no | levels the wave-period/tuner channel and seeds the startup attitude; the MEKF supplies the attitude used for heave |
| TVG-NLO (`TimeVarGainNLO_Adapter`) | 0.35 | – | 0.18 rad/s | no | 2–6 s bootstrap only, then the NLO takes over |
| Heave-only filters, `proxy` | 0.2 | 0.02 | 0.1 rad/s | no | the attitude every heave-only filter is levelled with |

The wave band of these records is 0.11–0.42 Hz. The PII corner sits inside
it, so each wave's horizontal acceleration pulls the estimated vertical.
`SeaStateFusionDefaults.h` already states the rule the shared observer
follows: `2Kp` "must stay an order of magnitude below the wave band, or the
observer levels itself against the orbital specific force instead of
gravity".

A sweep of IMU-only fixed gains found heave error flat for `2Kp` from 0.03
to 0.2. At `2Kp = 0.01` the gyro carries attitude too long and tilt reaches
4.4°. A slower setting (`2Kp = 0.05`) gave 4.98 % for `richter_zd`,
against 5.11 % at the shared gains, which is not worth a separate tuning.

### Do these settings help the other filters?

Each filter was rerun on its own simulator (OU-II, OU-III, TFG, NLO at the
OU noise model; PII at its own) with only the Mahony settings changed. Z RMS
error is in % of Hs and roll/pitch RMS is in degrees, over the last 900 s,
averaged over the eight records.

| Filter | Mahony variant | Z mean (worst) | Roll / pitch | Effect |
|---|---|---|---|---|
| OU-III | private observer 0.2/0.02 (default) | 4.00 (4.25) | 0.25 / 0.13 | – |
| OU-III | 0.05/0.0003 | 3.98 (4.19) | 0.71 / 0.22 | attitude and accelerometer-bias gates fail |
| OU-III | 0.1/0.005 | 4.00 (4.25) | 0.23 / 0.13 | none |
| OU-II | 0.2/0.02 (default) | 6.18 (6.63) | 0.19 / 0.18 | – |
| OU-II | 0.05/0.0003 | 6.19 (6.66) | 0.33 / 0.21 | slightly worse |
| OU-II | 0.1/0.005 | 6.18 (6.63) | 0.18 / 0.18 | none |
| TFG | 0.2/0.02 (default) | 4.06 (4.38) | 0.29 / 0.06 | – |
| TFG | 0.05/0.0003 | 4.03 (4.27) | 1.13 / 0.19 | 3-D gates fail |
| TFG | 0.1/0.005 | 4.05 (4.35) | 0.56 / 0.09 | worse attitude |
| NLO | bootstrap 0.35 (default), 0.2, 0.05 | 6.96 (7.32) | 0.25 / 0.14 | none |
| PII | adaptive 1.2–1.7, textbook magnetometer (default) | 5.82 (8.46) | – | – |
| PII | adaptive 1.2–1.7, heading-only magnetometer | 5.52 (7.55) | – | better |
| PII | 0.2/0.02, textbook magnetometer | 5.04 (5.91) | – | better |
| PII | 0.2/0.02, no magnetometer | 5.01 (5.91) | – | better, but no heading |
| PII | 0.2/0.02, heading-only magnetometer | **5.01 (5.88)** | – | better, heading 2.7° RMS |

- **OU-II, OU-III, TFG.** The private observer only seeds the startup
  attitude and levels the wave-period/tuner input. Heave uses the MEKF
  attitude, so slower settings do not help, and the slowest one delays the
  startup attitude enough to hurt attitude and accelerometer-bias estimates.
- **NLO.** Its Mahony runs only during the 2–6 s bootstrap.
- **PII.** The one filter whose heave is levelled by its Mahony.
  - The gains are what matter: the shared 0.2/0.02 cut its mean error from
    5.82 % to about 5.0 % and its worst case by a third.
  - With the magnetometer in the right frame it costs at most 0.03 points.
    Restricting it to heading makes it free, while keeping a 2.7° heading.
  - The PII rows come from its own simulator with the sensor mapping
    corrected to the firmware's (next section). That simulator runs its own
    noise model (5 mg accelerometer bias, 20 Hz magnetometer), so its numbers
    differ slightly from the `heave_only` table.
  - The shipped PII defaults are not changed here.

### Why dropping the magnetometer first looked like an improvement

More information should not make the estimate worse. Two defects made it
look that way.

- **A sensor-frame bug in the simulator adapter.** The PII simulator
  (`tests/pii_observer/pii_observer-adaptive.cpp`) mapped gyro and
  accelerometer to `(-E, -N, -D)` but the magnetometer to `(N, -E, -D)`.
  That put the magnetometer in a body frame turned 90° about z, which then
  rotated inconsistently with the gyro once the boat rolled and pitched. A
  comment there recorded that this mapping was chosen to fix a yaw error, and
  the yaw report added twice the declination to compensate. This simulator
  first copied the same adapter.
  - The firmware (`atomS3R_ins_pii_observer.ino`) maps all three sensors
    `[N, -E, -D]` and reports heading as `-yaw`. With that mapping the
    configuration's tilt error, with the textbook magnetometer correction,
    falls from 2.94° to 1.25° (worst 6.73° to 2.75°), and its heading error
    is 2.6° RMS. `heave_only-sim` and the PII simulator now use the firmware
    mapping.
- **Textbook Mahony lets the magnetometer steer tilt.** Its magnetometer
  error enters all three axes with the same gain as the accelerometer error.
  The magnetometer carries little tilt information, but its residual
  calibration errors (2 µT hard iron, 1.5 % scale, 1 % cross-axis, 1°
  misalignment in this replay) tilt the estimate.
  - Projecting the magnetometer error onto the estimated vertical, so it
    corrects heading only, fixes this. At the shared gains, tilt error goes
    from 0.63° (textbook) to 0.56°, against 0.58° without the magnetometer.
    Richter's filter goes from 5.22 % to 5.09 % of Hs (5.11 % without).
  - The shipped PII observer applies this projection
    (`ahrs/MahonyHeadingMag.h`; `Mahony_AHRS.h` is unchanged). At the PII's
    own fast adaptive gains it is a trade-off. The magnetometer carries
    acceleration-free tilt information about the axis normal to the field,
    and the projection drops it, so the `mahony` front end's tilt error rises
    from 1.25° to 1.53° (worst 2.75° to 3.36°). Heave still improves for PII
    (5.96 % to 5.64 % of Hs), Richter (11.08 % to 9.08 %) and Godhavn ext
    (14.09 % to 11.88 %); Küchler's worst case rises from 18.87 % to 19.88 %.

### Where the error comes from

I separated the error sources offline by feeding the filters the record's own
vertical acceleration, then adding one corruption at a time.

- **Levelling.**
  - *What goes wrong.* The shipped Mahony settings are tuned for a fast AHRS.
    A ~1 s time constant lets every wave's horizontal acceleration tilt the
    vertical.
  - *The result.* Tilt error reaches 3.36°. Leaked horizontal acceleration
    leaves a low-frequency vertical error. The heave filters' gain there
    (about 235 m per m/s² at 0.01 Hz for Godhavn) turns it into metres of
    drift.
  - *The fix.* The repository's QMEKF AHRS does no better (3.9–4.3° at
    Hs 8.5 m). Using the OU filters' own private Mahony settings instead
    (IMU-only, corner below the wave band) brings tilt to 0.58°. That is as
    good as a true vertical for every heave filter here. The PII observer
    also improves, from 5.64 % to 5.08 % of Hs.
  - *Why OU-III differs.* OU-III estimates attitude and wave motion jointly.
- **Bias instability (the published Godhavn law).**
  - On exact acceleration the law gives under 2 % of Hs.
  - Adding the 5e-4 m/s²/√s bias random walk makes Eq. 12 choose
    0.010–0.03 rad/s, and the drift passes almost unattenuated into heave.
    A better vertical reference cannot help with that.
  - The bias term in `godhavn_ext` and `richter_zd` lands within about
    2 points of the best fixed cutoff for each record.
- **Phase lead (standard filter).** Even with a perfect vertical, the
  standard structure stays at about 8 % of Hs. Richter's zero-displacement
  filter removes most of that lead and reaches 5.0–5.1 %.
- **The model (Küchler).**
  - *The symptom.* On exact, noise-free acceleration the EKF fits the
    acceleration to under 1 %. Yet its heave error is 33–41 % (Hs 1.5 m) and
    58–69 % (Hs 8.5 m) of heave RMS, with 4 to 27 modes and any tuning tried.
    On a synthetic three-tone sea it reaches 3.3 %.
  - *Why.* This boat's heave is nearly as broadband as the sea. 78 % of heave
    energy is at 0.05–0.12 Hz, but 75 % of acceleration energy is at
    0.12–0.4 Hz. The `-w_j^2 z_j` measurement sees only each mode's spring
    acceleration. The process noise that lets modes follow a random sea adds
    acceleration the measurement never sees.
  - *Context.* Küchler et al. validated on a 175 m ship, whose heave spectrum
    is narrow. No front end changes this result.

### How the papers demonstrated their performance

- **Küchler et al. 2011.** They used two setups:
  - an MSS simulation of a 175 m, 24,609 t vessel in a JONSWAP sea with
    Hs 7 m;
  - a winch rig whose IMU hangs on the rope with its z axis vertical, with an
    encoder as reference.

  In both, the accelerometer is (close to) vertical, so tilt never enters.
  Results are shown as plots only, with no RMS figures.
- **Richter et al. 2014.** They used the Liebherr AHC test bench: a tripod
  driven by three hydraulic winches with MSS-simulated heave, roll and pitch
  (JONSWAP, Hs 6 m, `w_p` = 0.63 rad/s), with a rope encoder as reference.
  - *Levelling.* It uses an attitude estimate "similar to" Küchler et al.'s
    ACC 2011 paper and to Kim and Golnaraghi 2004.
  - *Noise level.* It is taken from the IMU datasheet.
  - *Results.* They are given only relative to the standard filter
    (Table 1, RMS): lead-lag 55.6 %, zero-displacement 59.6 %, pole-zero
    43.4 %.
- **Neither paper isolated the levelling error.** The bench's attitude
  motion and the IMU's bias stability were favourable, and no absolute error
  against Hs was reported. The ACC 2011 attitude paper (Küchler, Pregizer,
  Eberharter, Schneider, Sawodny, pp. 2411–2416) could not be retrieved here.

### Heave over the last 60 s

The charts below show the reference against every estimate over the last 60 s
of the replay. The legends give each method's RMS error over that 60 s window.

With the shared `proxy` Mahony:

![JONSWAP Hs 1.5 m](../reports/results/heave_only/heave_only_proxy_jonswap_medium.svg)

![JONSWAP Hs 8.5 m](../reports/results/heave_only/heave_only_proxy_jonswap_high.svg)

With the Mahony as shipped:

![PM-Stokes Hs 8.5 m](../reports/results/heave_only/heave_only_pmstokes_high.svg)

All records, for both front ends, are in
[`reports/results/heave_only/`](../reports/results/heave_only/).
`heave_only_proxy_*` use the `proxy` front end; the others use the
Mahony as shipped.

## Reproduce

```bash
make ensure-sim-data
for d in kalman_ou_iii kalman_ou_ii kalman_tfg heave_only; do
  (cd tests/$d && make build && W3D_WRITE_TIMESERIES=1 ./run_tests.sh)
done
cd plots/heave_only && ./draw_plots.sh   # all front ends; --png adds previews
```

`draw_plots.sh` also writes `<method>_<wave>_<group>_zkin.svg`, with heave
and heave rate per heave-only filter.

These variables change the defaults for sweeps:
- `HO_MAHONY_TWOKP` and `HO_MAHONY_TWOKI` set fixed front-end gains
  (`2Ki` defaults to 0);
- `HO_FIXED_WC`, `HO_BIAS_TERM` and `HO_SUBTRACT_MEAN` for the Godhavn and
  Richter filters;
- `HO_KUCHLER_ZETA_Q`, `HO_KUCHLER_OMEGA_RW`, `HO_KUCHLER_BODY_Z` and
  `HO_KUCHLER_OFFSET_FIT`.

The simulator gates each method on its own regression sentinel.

## References

- J.-M. Godhavn, "Adaptive tuning of heave filter in motion sensor",
  OCEANS'98 Conference Proceedings, vol. 1, pp. 174–178, 1998. No open
  copy; implemented as specified by Richter et al. 2014.
- M. Richter, K. Schneider, D. Walser, O. Sawodny, "Real-time heave motion
  estimation using adaptive filtering techniques", 19th IFAC World Congress,
  pp. 10119–10125, 2014.
  [PDF](https://skoge.folk.ntnu.no/prost/proceedings/ifac2014/media/files/0111.pdf)
- S. Küchler, J. K. Eberharter, K. Langer, K. Schneider, O. Sawodny, "Heave
  motion estimation of a vessel using acceleration measurements", 18th IFAC
  World Congress, pp. 14742–14747, 2011.
  [PDF](https://skoge.folk.ntnu.no/prost/proceedings/ifac11-proceedings/data/html/papers/1935.pdf)
