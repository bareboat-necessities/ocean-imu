// Host-only covariance observer for the OU-III exact-map certificate.
//
// This translation unit deliberately reuses CertificateAdapter verbatim so the
// exact-map trace, production configuration, startup logic, measurement order,
// and quality gates remain identical.  It only adds the full estimator
// covariance at each exact-map block boundary.  Production sources are not
// modified.

#define main ou3_certificate_sim_embedded_main
#include "ou3-certificate-sim.cpp"
#undef main

namespace {

// Standard-conforming explicit-instantiation access to two private observer
// members.  This is confined to the host certificate translation unit; it does
// not change CertificateAdapter or the deployed filter.
template<class Tag, typename Tag::type Member>
struct PrivateMemberAccess {
    friend typename Tag::type get_private_member(Tag) { return Member; }
};

struct FusionMemberTag {
    using type = CertificateAdapter::Fusion CertificateAdapter::*;
    friend type get_private_member(FusionMemberTag);
};
template struct PrivateMemberAccess<FusionMemberTag, &CertificateAdapter::fusion_>;

struct MapStartedMemberTag {
    using type = bool CertificateAdapter::*;
    friend type get_private_member(MapStartedMemberTag);
};
template struct PrivateMemberAccess<MapStartedMemberTag, &CertificateAdapter::map_block_started_>;

class InformationCertificateAdapter final : public IW3dFusionAdapter {
public:
    InformationCertificateAdapter(bool with_mag,
                                  const Vector3f& sigma_a_init,
                                  const Vector3f& sigma_g,
                                  const Vector3f& sigma_m)
        : inner_(with_mag, sigma_a_init, sigma_g, sigma_m)
    {
        const char* explicit_path = std::getenv("OU3_CERT_COV_TRACE");
        if (explicit_path && *explicit_path) {
            cov_path_ = explicit_path;
        } else {
            const char* map_path = std::getenv("OU3_CERT_MAP_TRACE");
            if (!map_path || !*map_path)
                throw std::runtime_error("OU3_CERT_MAP_TRACE is required for covariance trace derivation");
            cov_path_ = map_path;
            const std::string suffix = "_exact_maps.bin";
            const auto pos = cov_path_.rfind(suffix);
            if (pos != std::string::npos)
                cov_path_.replace(pos, suffix.size(), "_covariance.bin");
            else
                cov_path_ += ".covariance.bin";
        }

        cov_trace_.open(cov_path_, std::ios::binary);
        if (!cov_trace_) throw std::runtime_error("cannot open OU3 certificate covariance trace");
        write_cov_header();
    }

    ~InformationCertificateAdapter() override = default;

    void updateMag(const Vector3f& mag_body_ned) override
    {
        begin_cov_block_if_needed();
        const bool before = map_started();
        inner_.updateMag(mag_body_ned);
        const bool after = map_started();
        // updateMag never advances the 50-sample block counter, but retain this
        // defensive branch so a future source change cannot desynchronize the
        // two traces silently.
        if (before && !after) close_cov_block();
    }

    void update(float dt,
                const Vector3f& gyr_meas_ned,
                const Vector3f& acc_meas_ned,
                float temperature_c) override
    {
        begin_cov_block_if_needed();
        const bool before = map_started();
        inner_.update(dt, gyr_meas_ned, acc_meas_ned, temperature_c);
        time_s_ += dt;
        const bool after = map_started();
        if (before && !after) close_cov_block();
    }

    FilterSnapshot snapshot() const override
    {
        return inner_.snapshot();
    }

private:
    CertificateAdapter::Fusion& fusion()
    {
        const auto member = get_private_member(FusionMemberTag{});
        return inner_.*member;
    }

    const CertificateAdapter::Fusion& fusion() const
    {
        const auto member = get_private_member(FusionMemberTag{});
        return inner_.*member;
    }

    bool map_started() const
    {
        const auto member = get_private_member(MapStartedMemberTag{});
        return inner_.*member;
    }

    Matrix21f covariance() const
    {
        return fusion().raw().mekf().covariance_full();
    }

    void begin_cov_block_if_needed()
    {
        if (map_started() || cov_block_started_) return;
        cov_block_started_ = true;
        cov_start_time_ = time_s_;
        cov_start_ = covariance();
    }

    void close_cov_block()
    {
        if (!cov_block_started_)
            throw std::runtime_error("OU3 covariance trace lost map-block start alignment");
        const double t0 = static_cast<double>(cov_start_time_);
        const double t1 = static_cast<double>(time_s_);
        write_binary(cov_trace_, t0);
        write_binary(cov_trace_, t1);
        const Matrix21f end = covariance();
        for (int i = 0; i < kNX; ++i)
            for (int j = 0; j < kNX; ++j)
                write_binary(cov_trace_, cov_start_(i,j));
        for (int i = 0; i < kNX; ++i)
            for (int j = 0; j < kNX; ++j)
                write_binary(cov_trace_, end(i,j));
        cov_block_started_ = false;
    }

    void write_cov_header()
    {
        const char magic[8] = {'O','U','3','C','O','V','1','\0'};
        cov_trace_.write(magic, sizeof(magic));
        const std::uint32_t version = 1;
        const std::uint32_t nx = kNX;
        const std::uint32_t stride = static_cast<std::uint32_t>(
            env_positive_int("OU3_CERT_MAP_STRIDE", 50));
        write_binary(cov_trace_, version);
        write_binary(cov_trace_, nx);
        write_binary(cov_trace_, stride);
    }

    CertificateAdapter inner_;
    std::ofstream cov_trace_;
    std::string cov_path_;
    bool cov_block_started_ = false;
    float cov_start_time_ = 0.0f;
    float time_s_ = 0.0f;
    Matrix21f cov_start_ = Matrix21f::Zero();
};

void process_information_one(const std::string& filename,
                             bool with_mag,
                             const W3dRandomSeeds& seeds,
                             bool write_timeseries,
                             float validation_window_sec)
{
    auto result = process_wave_file_for_tracker<InformationCertificateAdapter>(
        filename, kDt, with_mag, true, kMagOdrHz,
        "_fusion_ou3_cert", "_fusion_ou3_cert_nomag", seeds, write_timeseries);
    if (!result) return;
    if (validation_window_sec > 0.0f)
        print_validation_metrics(*result, kDt, validation_window_sec, "OU_III_CERT");
    static constexpr W3dSummaryLabels labels{.target = "RS_target", .applied = "RS_applied"};
    print_summary_and_fail_if_needed(*result, kDt, Limits::value, labels);
}

} // namespace

int main(int argc, char** argv)
{
    bool with_mag = true;
    std::vector<std::string> files;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--nomag") with_mag = false;
        else if (arg == "--input" && i + 1 < argc) files.emplace_back(argv[++i]);
        else if (arg == "--help") {
            std::cout << "Usage: " << argv[0] << " [--nomag] [--input PATH]...\n";
            return 0;
        } else {
            std::cerr << "ERROR: unknown or incomplete argument: " << arg << "\n";
            return 2;
        }
    }
    if (files.empty()) files = collect_wave_data_files(".");

    W3dRandomSeeds seeds;
    try { seeds = w3d_random_seeds_from_env(); }
    catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }

    bool write_timeseries = true;
    if (const char* raw = std::getenv("W3D_WRITE_TIMESERIES"))
        write_timeseries = std::string(raw) != "0";
    float validation_window_sec = 900.0f;
    if (const char* raw = std::getenv("W3D_VALIDATION_WINDOW_SEC"))
        validation_window_sec = static_cast<float>(std::atof(raw));

    if (files.size() != 1u) {
        std::cerr << "ERROR: ou3-information-certificate-sim requires exactly one --input\n";
        return 2;
    }

    try { process_information_one(files.front(), with_mag, seeds, write_timeseries, validation_window_sec); }
    catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
    return w3d_any_quality_gate_failed() ? 1 : 0;
}
