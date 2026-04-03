#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path


METRIC_METADATA = {
    "valid": ("Valid report", "1 if the spectral report passed internal validity checks, else 0."),
    "lowfreq_guard_hz": ("Low-frequency guard [Hz]", "Minimum trusted frequency guard used to avoid unstable low-frequency integration."),
    "observation_duration_s": ("Observation duration [s]", "Total equivalent time window represented by the spectrum."),
    "m_neg1": ("Spectral moment m-1", "Negative-first moment of the spectrum (emphasizes low frequencies)."),
    "m0": ("Spectral moment m0", "Zeroth spectral moment (total variance / energy proxy)."),
    "m1": ("Spectral moment m1", "First spectral moment (energy weighted by frequency)."),
    "m2": ("Spectral moment m2", "Second spectral moment (used in zero-crossing period metrics)."),
    "m3": ("Spectral moment m3", "Third spectral moment (high-frequency weighted moment)."),
    "m4": ("Spectral moment m4", "Fourth spectral moment (very high-frequency weighted moment)."),
    "mu2": ("Central moment μ2", "Second central moment of the spectral distribution."),
    "mu3": ("Central moment μ3", "Third central moment used for spectral skewness."),
    "mu4": ("Central moment μ4", "Fourth central moment used for spectral kurtosis."),
    "fp_hz": ("Peak frequency fp [Hz]", "Frequency at which spectral density is maximal."),
    "fp_rad": ("Peak angular frequency ωp [rad/s]", "Angular form of peak frequency."),
    "tp_s": ("Peak period Tp [s]", "Period corresponding to peak frequency."),
    "mean_frequency_hz": ("Mean frequency [Hz]", "Average frequency of the spectrum (energy-weighted)."),
    "mean_frequency_rad": ("Mean angular frequency [rad/s]", "Angular-frequency form of mean spectral frequency."),
    "t1_s": ("Mean period T1 [s]", "Mean wave period from low-order moments."),
    "tm02_s": ("Zero-crossing period Tm02 [s]", "Mean zero-upcrossing period from m0/m2."),
    "te_s": ("Energy period Te [s]", "Energy flux period often used for wave power."),
    "tm0m1_s": ("Period Tm-10 [s]", "Moment-derived period from m-1 and m0."),
    "tm1m1_s": ("Period Tm1-1 [s]", "Moment-derived period from m1 and m-1 equivalents in model."),
    "mean_group_period_s": ("Mean group period [s]", "Estimated mean duration of wave groups."),
    "rms_displacement_m": ("RMS displacement [m]", "Root-mean-square sea-surface elevation."),
    "hs_m": ("Significant wave height Hs [m]", "Estimated significant wave height (typically 4·sqrt(m0))."),
    "h1_10_crest_m": ("1/10 crest height [m]", "Representative average of highest 10% crest heights."),
    "h1_100_crest_m": ("1/100 crest height [m]", "Representative average of highest 1% crest heights."),
    "most_probable_max_crest_m": ("Most probable max crest [m]", "Most probable largest crest in the observation window."),
    "expected_max_crest_m": ("Expected max crest [m]", "Expected value of maximum crest over the window."),
    "upcrossing_rate_hz": ("Upcrossing rate [Hz]", "Estimated mean rate of zero upcrossings."),
    "downcrossing_rate_hz": ("Downcrossing rate [Hz]", "Estimated mean rate of zero downcrossings."),
    "wave_count.expected": ("Expected wave count", "Expected number of waves in the observation duration."),
    "wave_count.ci_lower": ("Wave count CI lower", "Lower bound of confidence interval for wave count."),
    "wave_count.ci_upper": ("Wave count CI upper", "Upper bound of confidence interval for wave count."),
    "wave_count.confidence": ("Wave count confidence", "Confidence level used for wave-count interval."),
    "rbw": ("Relative bandwidth", "Relative spectral bandwidth indicator."),
    "regularity_spec": ("Spectral regularity", "Regularity index based on spectrum shape."),
    "regularity_phase": ("Phase regularity", "Regularity index based on phase progression assumptions."),
    "rbw_phase_increment": ("RBW phase increment", "Phase-based increment used by RBW-related metric."),
    "narrowness_nu": ("Narrowness ν", "Spectral narrowness parameter."),
    "bandwidth_clh": ("Bandwidth CLH", "Cartwright-Longuet-Higgins style bandwidth metric."),
    "bandwidth_goda": ("Bandwidth Goda", "Goda spectral bandwidth metric."),
    "bandwidth_kuik": ("Bandwidth Kuik", "Kuik bandwidth metric."),
    "bandwidth_epsilon": ("Bandwidth ε", "Alternative epsilon-based spectral bandwidth."),
    "longuet_higgins_width": ("Longuet-Higgins width", "Longuet-Higgins spectral width parameter."),
    "spectral_bandwidth_hz": ("Spectral bandwidth [Hz]", "Bandwidth expressed in frequency units."),
    "spectral_narrowness_ratio": ("Spectral narrowness ratio", "Ratio indicating concentration around spectral peak."),
    "spectral_skewness": ("Spectral skewness", "Skewness of spectral distribution."),
    "spectral_kurtosis": ("Spectral kurtosis", "Kurtosis of spectral distribution."),
    "spectral_excess_kurtosis": ("Spectral excess kurtosis", "Kurtosis above Gaussian baseline (kurtosis-3)."),
    "peakedness_ochi_q": ("Ochi peakedness q", "Peakedness parameter in Ochi-style characterization."),
    "benassai_parameter": ("Benassai parameter", "Sea-state parameter used in selected nonlinearity characterization."),
    "peak_enhancement_gamma": ("Peak enhancement γ", "JONSWAP-like peak enhancement factor."),
    "peakedness_goda": ("Goda peakedness", "Peakedness metric following Goda formulation."),
    "spectrum_type": ("Spectrum type", "Sea-state class inferred from spectrum (e.g., Swell/Mixed/WindSea)."),
    "deep_water_wavelength_tp_m": ("Deep-water wavelength at Tp [m]", "Deep-water wavelength computed from Tp."),
    "wave_steepness": ("Wave steepness", "Representative steepness from height-to-length scaling."),
    "nonlinearity_parameter": ("Nonlinearity parameter", "Aggregate indicator of wave nonlinearity."),
    "ursell_number": ("Ursell number", "Nonlinearity-to-dispersion ratio used in wave theory."),
    "wave_age": ("Wave age", "Indicator of wave development/maturity."),
    "crest_exceed_prob_hs": ("Crest exceedance at Hs", "Probability crest exceeds an Hs-based threshold."),
    "crest_exceed_prob_tayfun_hs": ("Tayfun crest exceedance at Hs", "Nonlinear Tayfun-model exceedance probability at Hs threshold."),
    "pot_mean_excess_m": ("POT mean excess [m]", "Mean excess above selected peaks-over-threshold level."),
    "weibull_return_height_m": ("Weibull return height [m]", "Return height estimated from Weibull-tail fit."),
    "groupiness_factor": ("Groupiness factor", "Degree of wave grouping/intermittency."),
    "bfi": ("Benjamin-Feir index (BFI)", "Instability index linked to modulational instability likelihood."),
    "mean_group_length": ("Mean group length", "Average number or extent of waves per group (model-dependent)."),
    "group_height_factor": ("Group height factor", "Amplification factor of grouped-wave heights."),
    "expected_run_length_above_hs": ("Expected run length above Hs", "Expected consecutive waves exceeding Hs threshold."),
    "wave_energy_density_j_m2": ("Wave energy density [J/m²]", "Depth-integrated wave energy per sea-surface area."),
    "wave_power_w_m": ("Wave power [W/m]", "Wave energy transport per meter crest width."),
    "energy_flux_period_s": ("Energy flux period [s]", "Period representative for energy flux computations."),
    "breaking_probability": ("Breaking probability", "Estimated probability of wave breaking."),
    "depth_limited_breaking_height_m": ("Depth-limited breaking height [m]", "Depth-constrained breaking wave height estimate."),
    "bottom_orbital_velocity_mps": ("Bottom orbital velocity [m/s]", "Estimated near-bed wave orbital velocity."),
    "radiation_stress_xx_n_m": ("Radiation stress Sxx [N/m]", "Radiation stress tensor xx component."),
    "radiation_stress_yy_n_m": ("Radiation stress Syy [N/m]", "Radiation stress tensor yy component."),
    "radiation_stress_xy_n_m": ("Radiation stress Sxy [N/m]", "Radiation stress tensor xy/shear component."),
    "msdv": ("MSDV", "Motion sickness dose value indicator."),
    "seasickness_incidence_pct": ("Seasickness incidence [%]", "Estimated fraction of people likely seasick."),
    "motion_comfort_level_0_100": ("Motion comfort level [0-100]", "Composite comfort index where higher generally means more comfort."),
    "vertical_motion_intensity": ("Vertical motion intensity", "Relative intensity of vertical vessel/platform motion."),
    "motion_character": ("Motion character", "Qualitative or coded motion regime descriptor."),
    "time_to_onset_min": ("Time to onset [min]", "Estimated time to motion-sickness symptom onset."),
    "snr_db": ("SNR [dB]", "Signal-to-noise ratio of the derived spectrum/metrics."),
    "temporal_stability_hs": ("Temporal stability of Hs", "Stability index of Hs over time."),
    "temporal_stability_tp": ("Temporal stability of Tp", "Stability index of Tp over time."),
    "data_quality_0_1": ("Data quality [0-1]", "Overall data quality confidence score."),
    "wec_capture_width_ratio": ("WEC capture width ratio", "Wave-energy-converter capture-width performance ratio."),
    "sea_swell_partition": ("Sea/swell partition", "Partition indicator between wind sea and swell components."),
    "wave_age_class": ("Wave age class", "Categorical class derived from wave age."),
}


def tex_escape(text: str) -> str:
    return (
        text.replace("\\", r"\textbackslash{}")
        .replace("_", r"\_")
        .replace("%", r"\%")
        .replace("&", r"\&")
        .replace("#", r"\#")
    )


def pass_cell(v: str) -> str:
    return r"\textbf{PASS}" if str(v).strip() in ("1", "true", "True") else r"\textbf{FAIL}"


def describe_case(case_name: str) -> str:
    stem = Path(case_name).stem
    mode = ""
    if stem.startswith("spectrum_estimator_"):
        mode = "Goertzel"
        stem = stem[len("spectrum_estimator_"):]
    elif stem.startswith("spectrum_wavelets_"):
        mode = "Wavelets"
        stem = stem[len("spectrum_wavelets_"):]

    if "_" not in stem:
        return case_name

    wave_type, _, encoded = stem.partition("_")
    token_pattern = re.compile(r"([A-Z])([-+]?[0-9]*\.?[0-9]+)")
    values = {key: value for key, value in token_pattern.findall(encoded)}

    required = {"H", "L", "A", "P"}
    if not required.issubset(values):
        return case_name

    height = float(values["H"])
    length = float(values["L"])
    azimuth = float(values["A"])
    phase = float(values["P"])
    wave_type = wave_type.upper()
    details = f"{wave_type}: H={height:.2f} m, az={azimuth:.1f}°"
    return details


def build_table(rows, caption: str, label: str) -> str:
    lines = [
        r"\begin{table}[H]",
        r"\centering",
        rf"\caption{{{caption}}}",
        rf"\label{{{label}}}",
        r"\begin{tabular}{lrrrrrcl}",
        r"\toprule",
        r"Case & $H_{s,target}$ & $H_{s,est}$ & $H_s$ err [\%] & $T_{p,target}$ & $T_{p,est}$ \\",
        r"\midrule",
    ]

    for row in rows:
        line = (
            f"{tex_escape(describe_case(row['case_name']))} & "
            f"{float(row['hs_target_m']):.3f} & "
            f"{float(row['hs_est_m']):.3f} & "
            f"{float(row['hs_error_pct']):.2f} & "
            f"{float(row['tp_target_s']):.3f} & "
            f"{float(row['tp_est_s']):.3f} "
            r"\\"
        )
        lines.append(line)

    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}"])
    return "\n".join(lines) + "\n"


def build_full_metrics_block(full_metrics_rows, last_case_description: str | None = None) -> str:
    if not full_metrics_rows:
        return ""

    title = r"\subsection*{Full metrics}"
    if last_case_description:
        title = (
            r"\subsection*{Full metrics (last synthetic spectrum case: "
            + tex_escape(last_case_description)
            + ")}"
        )

    lines = [
        title,
        r"\begin{longtable}{p{0.23\linewidth}p{0.25\linewidth}p{0.35\linewidth}p{0.13\linewidth}}",
        r"\hline",
        r"Metric key & Metrics Name & Meaning & Value \\",
        r"\hline",
        r"\endfirsthead",
        r"\hline",
        r"Metric key & Metrics Name & Meaning & Value \\",
        r"\hline",
        r"\endhead",
    ]

    for row in full_metrics_rows:
        metric = row["metric"]
        value = format_metric_value(row["value"])
        display_name, meaning = METRIC_METADATA.get(
            metric,
            ("(unmapped metric)", "Metric emitted by report but not yet described in this script."),
        )
        lines.append(
            f"{tex_escape(metric)} & {tex_escape(display_name)} & {tex_escape(meaning)} & {tex_escape(value)} "
            r"\\"
        )

    lines.extend([r"\hline", r"\end{longtable}", ""])
    return "\n".join(lines)


def format_metric_value(raw_value: str) -> str:
    try:
        value = float(raw_value)
    except ValueError:
        return raw_value
    return f"{value:.6g}"


def parse_name_value_lines(path: Path):
    rows = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        metric, value = line.split("=", 1)
        rows.append({"metric": metric.strip(), "value": value.strip()})
    return rows


def build_fragment(rows, caption: str, label: str, full_metrics_rows=None) -> str:
    last_case_description = describe_case(rows[-1]["case_name"]) if rows else None
    return build_table(rows, caption, label) + build_full_metrics_block(
        full_metrics_rows or [],
        last_case_description,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", default="sea_metrics_from_spectrum_report.csv")
    parser.add_argument("--full-metrics-name-value", default="sea_metrics_from_spectrum_full_metrics.txt")
    parser.add_argument("--out-fragment", default="../../doc/spectrum/spectrum-sea_metrics.tex-part")
    parser.add_argument("--fragment-input", default="spectrum-sea_metrics.tex-part")
    parser.add_argument("--title", default="SeaMetricsFromSpectrum Test Report")
    parser.add_argument("--caption", default="SeaMetricsFromSpectrum validation from estimator spectra.")
    parser.add_argument("--label", default="tab:sea_metrics_from_spectrum")
    args = parser.parse_args()

    csv_path = Path(args.csv)
    fragment_path = Path(args.out_fragment)

    with csv_path.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))

    full_metrics_rows = []
    full_metrics_path = Path(args.full_metrics_name_value)
    if full_metrics_path.exists():
        full_metrics_rows = parse_name_value_lines(full_metrics_path)

    fragment_path.parent.mkdir(parents=True, exist_ok=True)
    fragment_path.write_text(build_fragment(rows, args.caption, args.label, full_metrics_rows), encoding="utf-8")

    print(f"Loaded CSV rows: {len(rows)} from {csv_path}")
    if full_metrics_rows:
        print(f"Loaded full metrics name=value rows: {len(full_metrics_rows)} from {full_metrics_path}")
    else:
        print(f"Full metrics name=value file not found (optional): {full_metrics_path}")
    print(f"Wrote LaTeX fragment: {fragment_path}")


if __name__ == "__main__":
    main()
