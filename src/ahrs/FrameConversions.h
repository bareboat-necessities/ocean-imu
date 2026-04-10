#pragma once

/*
   Copyright 2025, Mikhail Grushinskiy
 */

#ifdef ARDUINO
#include <ArduinoEigenDense.h>
#else
#include <Eigen/Dense>
#endif

#include <cmath>
#include <numbers>

#ifdef FRAMECONV_TEST
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#endif

using Eigen::Vector3f;
using Eigen::Matrix3f;
using Eigen::Quaternionf;

// Coordinate conversions
//   Nautical Z-up (ENU: East, North, Up)
//   Aerospace NED (North, East, Down)
//
// Vector basis map:
//   v_ned = [ v_zu.y, v_zu.x, -v_zu.z ]

static inline Vector3f zu_to_ned(const Vector3f& v) {
    return Vector3f(v.y(), v.x(), -v.z());
}

static inline Vector3f ned_to_zu(const Vector3f& v) {
    return Vector3f(v.y(), v.x(), -v.z());
}

static inline Matrix3f basis_zu_to_ned() {
    Matrix3f S;
    S << 0.0f, 1.0f,  0.0f,
         1.0f, 0.0f,  0.0f,
         0.0f, 0.0f, -1.0f;
    return S;
}

// Quaternion helpers
//
// quat_from_euler() is just a ZYX matrix/quaternion builder:
//   R = Rz(yaw) * Ry(pitch) * Rx(roll)
//
// When used for aerospace Euler, interpret R as BODY->WORLD in NED.
// When used for nautical Euler in this file, interpret R as WORLD->BODY in ENU/Z-up.

static inline Quaternionf quat_from_euler(float roll_deg, float pitch_deg, float yaw_deg) {
    const float hr = roll_deg  * std::numbers::pi_v<float> / 180.0f * 0.5f;
    const float hp = pitch_deg * std::numbers::pi_v<float> / 180.0f * 0.5f;
    const float hy = yaw_deg   * std::numbers::pi_v<float> / 180.0f * 0.5f;

    const float cr = std::cos(hr);
    const float sr = std::sin(hr);
    const float cp = std::cos(hp);
    const float sp = std::sin(hp);
    const float cy = std::cos(hy);
    const float sy = std::sin(hy);

    Quaternionf q;
    q.w() = cy * cp * cr + sy * sp * sr;
    q.x() = cy * cp * sr - sy * sp * cr;
    q.y() = sy * cp * sr + cy * sp * cr;
    q.z() = sy * cp * cr - cy * sp * sr;
    return q.normalized();
}

static inline void matrix_to_euler_zyx_deg(const Matrix3f& R,
                                           float& roll_deg,
                                           float& pitch_deg,
                                           float& yaw_deg)
{
    const float pitch = std::atan2(
        -R(2, 0),
        std::sqrt(R(0, 0) * R(0, 0) + R(1, 0) * R(1, 0)));
    const float roll = std::atan2(R(2, 1), R(2, 2));
    const float yaw = std::atan2(R(1, 0), R(0, 0));

    constexpr float RAD2DEG = 180.0f / std::numbers::pi_v<float>;
    roll_deg = roll * RAD2DEG;
    pitch_deg = pitch * RAD2DEG;
    yaw_deg = yaw * RAD2DEG;
}

// Aerospace Euler convention in this file:
//   BODY->WORLD in NED.
static inline void quat_to_euler_aero(const Quaternionf& q,
                                      float& roll,
                                      float& pitch,
                                      float& yaw)
{
    const Matrix3f C_bw_ned = q.normalized().toRotationMatrix();
    matrix_to_euler_zyx_deg(C_bw_ned, roll, pitch, yaw);
}

// Convert BODY->WORLD NED rotation to WORLD->BODY ENU/Z-up rotation.
static inline Matrix3f rot_bw_ned_to_wb_zu(const Matrix3f& C_bw_ned) {
    const Matrix3f S = basis_zu_to_ned();
    return S * C_bw_ned.transpose() * S;
}

// Convert WORLD->BODY ENU/Z-up rotation to BODY->WORLD NED rotation.
static inline Matrix3f rot_wb_zu_to_bw_ned(const Matrix3f& C_wb_zu) {
    const Matrix3f S = basis_zu_to_ned();
    return S * C_wb_zu.transpose() * S;
}

// Nautical Euler convention in this file:
//   WORLD->BODY in ENU/Z-up.
//
// This preserves the same near-level sign behavior as the old helper:
//   roll_n  ~= -pitch_a
//   pitch_n ~= -roll_a
//   yaw_n   ~=  yaw_a
//
// but is correct for finite angles because it converts the full rotation first.
static inline void quat_to_euler_nautical(const Quaternionf& q_bw_ned,
                                          float& roll,
                                          float& pitch,
                                          float& yaw)
{
    const Matrix3f C_bw_ned = q_bw_ned.normalized().toRotationMatrix();
    const Matrix3f C_wb_zu = rot_bw_ned_to_wb_zu(C_bw_ned);
    matrix_to_euler_zyx_deg(C_wb_zu, roll, pitch, yaw);
}

// Finite-angle conversion: aerospace BODY->WORLD NED -> nautical WORLD->BODY ENU/Z-up.
static inline void aero_to_nautical(float& roll, float& pitch, float& yaw) {
    const Matrix3f C_bw_ned = quat_from_euler(roll, pitch, yaw).toRotationMatrix();
    const Matrix3f C_wb_zu = rot_bw_ned_to_wb_zu(C_bw_ned);
    matrix_to_euler_zyx_deg(C_wb_zu, roll, pitch, yaw);
}

// Finite-angle conversion: nautical WORLD->BODY ENU/Z-up -> aerospace BODY->WORLD NED.
static inline void nautical_to_aero(float& roll, float& pitch, float& yaw) {
    const Matrix3f C_wb_zu = quat_from_euler(roll, pitch, yaw).toRotationMatrix();
    const Matrix3f C_bw_ned = rot_wb_zu_to_bw_ned(C_wb_zu);
    matrix_to_euler_zyx_deg(C_bw_ned, roll, pitch, yaw);
}

// Magnetometer simulation using World Magnetic Model (WMM)
// Defaults: Statue of Liberty, USA, Sept 2025, elevation ~0 m

struct MagSim_WMM {
    static constexpr float default_declination_deg = -12.6f; // [deg] east is positive
    static constexpr float default_inclination_deg =  66.5f; // [deg] positive = down
    static constexpr float default_total_field_uT  = 52.0f;  // [uT]

    // World magnetic field in ENU (East, North, Up)
    // East = X, North = Y, Up = Z
    // Units: microteslas [uT]
    static Eigen::Vector3f mag_world_nautical(
        float declination_deg = default_declination_deg,
        float inclination_deg = default_inclination_deg,
        float total_field_uT  = default_total_field_uT)
    {
        const float dec_rad  = declination_deg * std::numbers::pi_v<float> / 180.0f;
        const float incl_rad = inclination_deg * std::numbers::pi_v<float> / 180.0f;

        const float h = std::cos(incl_rad);  // horizontal fraction
        const float v = -std::sin(incl_rad); // vertical in Z-up

        Eigen::Vector3f mag_world;
        mag_world.x() = h * std::sin(dec_rad); // East
        mag_world.y() = h * std::cos(dec_rad); // North
        mag_world.z() = v;                     // Up

        return mag_world * total_field_uT;
    }

    // World magnetic field vector in aerospace NED frame (North, East, Down)
    static Eigen::Vector3f mag_world_aero(
        float declination_deg = default_declination_deg,
        float inclination_deg = default_inclination_deg,
        float total_field_uT  = default_total_field_uT)
    {
        return zu_to_ned(mag_world_nautical(declination_deg, inclination_deg, total_field_uT));
    }

    // Simulate body-frame magnetometer [uT] from nautical Euler (deg)
    // Input: nautical Euler = WORLD->BODY in ENU/Z-up
    // Output: body-frame magnetometer in nautical frame (body ENU/Z-up)
    static Eigen::Vector3f simulate_mag_from_euler_nautical(
        float roll_deg, float pitch_deg, float yaw_deg,
        float declination_deg = default_declination_deg,
        float inclination_deg = default_inclination_deg,
        float total_field_uT  = default_total_field_uT)
    {
        const Matrix3f C_wb_zu = quat_from_euler(roll_deg, pitch_deg, yaw_deg).toRotationMatrix();
        const Eigen::Vector3f mag_world_zu =
            mag_world_nautical(declination_deg, inclination_deg, total_field_uT);
        return C_wb_zu * mag_world_zu;
    }

    // Simulate body-frame magnetometer [uT] from aerospace Euler (deg)
    // Input: aerospace Euler = BODY->WORLD in NED
    // Output: body-frame magnetometer in nautical frame (body ENU/Z-up)
    static Eigen::Vector3f simulate_mag_from_euler_aero(
        float roll_deg, float pitch_deg, float yaw_deg,
        float declination_deg = default_declination_deg,
        float inclination_deg = default_inclination_deg,
        float total_field_uT  = default_total_field_uT)
    {
        const Matrix3f C_bw_ned = quat_from_euler(roll_deg, pitch_deg, yaw_deg).toRotationMatrix();
        const Eigen::Vector3f mag_world_ned =
            mag_world_aero(declination_deg, inclination_deg, total_field_uT);
        const Eigen::Vector3f mag_body_ned = C_bw_ned.transpose() * mag_world_ned;
        return ned_to_zu(mag_body_ned);
    }
};

#ifdef FRAMECONV_TEST

inline void assert_close(float a, float b, float tol, const char* msg) {
    if (std::fabs(a - b) > tol) {
        std::cerr << "ASSERT FAIL: " << msg
                  << "  got=" << a << " expected=" << b
                  << " tol=" << tol << "\n";
        assert(false);
    }
}

inline void assert_mat_close(const Matrix3f& A, const Matrix3f& B, float tol, const char* msg) {
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            if (std::fabs(A(r, c) - B(r, c)) > tol) {
                std::cerr << "ASSERT FAIL: " << msg
                          << "  at (" << r << "," << c << ")"
                          << "  got=" << A(r, c)
                          << " expected=" << B(r, c)
                          << " tol=" << tol << "\n";
                assert(false);
            }
        }
    }
}

inline int test_frame_conversions() {
    const float tol_angle = 1e-3f;
    const float tol_vec   = 1e-6f;
    const float tol_mat   = 1e-5f;

    float r_n, p_n, y_n;
    float r_a, p_a, y_a;

    // Flat
    r_n = 0.0f; p_n = 0.0f; y_n = 0.0f;
    r_a = r_n; p_a = p_n; y_a = y_n;
    nautical_to_aero(r_a, p_a, y_a);
    aero_to_nautical(r_a, p_a, y_a);
    assert_close(r_a, r_n, tol_angle, "Flat roll");
    assert_close(p_a, p_n, tol_angle, "Flat pitch");
    assert_close(y_a, y_n, tol_angle, "Flat yaw");

    // Small-angle sign behavior preserved
    r_a = 1.0f; p_a = 0.0f; y_a = 0.0f;
    aero_to_nautical(r_a, p_a, y_a);
    assert_close(r_a, 0.0f, tol_angle, "Small-angle aero->nautical roll from roll");
    assert_close(p_a, -1.0f, tol_angle, "Small-angle aero->nautical pitch from roll");

    r_a = 0.0f; p_a = 1.0f; y_a = 0.0f;
    aero_to_nautical(r_a, p_a, y_a);
    assert_close(r_a, -1.0f, tol_angle, "Small-angle aero->nautical roll from pitch");
    assert_close(p_a, 0.0f, tol_angle, "Small-angle aero->nautical pitch from pitch");

    r_a = 0.0f; p_a = 0.0f; y_a = 1.0f;
    aero_to_nautical(r_a, p_a, y_a);
    assert_close(y_a, 1.0f, tol_angle, "Small-angle aero->nautical yaw");

    // General round-trip
    r_n = 30.0f; p_n = 20.0f; y_n = 45.0f;
    r_a = r_n; p_a = p_n; y_a = y_n;
    nautical_to_aero(r_a, p_a, y_a);
    aero_to_nautical(r_a, p_a, y_a);
    assert_close(r_a, r_n, tol_angle, "General roll");
    assert_close(p_a, p_n, tol_angle, "General pitch");
    assert_close(y_a, y_n, tol_angle, "General yaw");

    // Gravity vector round-trip
    Vector3f g_n(0.0f, 0.0f, -9.81f);
    Vector3f g_a = zu_to_ned(g_n);
    Vector3f g_back = ned_to_zu(g_a);
    assert_close(g_back.x(), g_n.x(), tol_vec, "Gravity x");
    assert_close(g_back.y(), g_n.y(), tol_vec, "Gravity y");
    assert_close(g_back.z(), g_n.z(), tol_vec, "Gravity z");

    // Magnetometer vector round-trip
    Vector3f m_n(1.0f, 0.0f, 0.0f);
    Vector3f m_a = zu_to_ned(m_n);
    Vector3f m_back = ned_to_zu(m_a);
    assert_close(m_back.x(), m_n.x(), tol_vec, "Mag x");
    assert_close(m_back.y(), m_n.y(), tol_vec, "Mag y");
    assert_close(m_back.z(), m_n.z(), tol_vec, "Mag z");

    // Quaternion / matrix consistency: nautical -> aerospace should match same physical rotation.
    {
        const float rn = 30.0f, pn = 20.0f, yn = 45.0f;
        float ra = rn, pa = pn, ya = yn;
        nautical_to_aero(ra, pa, ya);

        const Matrix3f C_wb_zu = quat_from_euler(rn, pn, yn).toRotationMatrix();
        const Matrix3f C_bw_from_n = rot_wb_zu_to_bw_ned(C_wb_zu);
        const Matrix3f C_bw_from_a = quat_from_euler(ra, pa, ya).toRotationMatrix();
        assert_mat_close(C_bw_from_n, C_bw_from_a, tol_mat, "Quaternion/matrix consistency");
    }

    // Randomized stress test
    for (int i = 0; i < 100; ++i) {
        float rn = (float(rand()) / RAND_MAX - 0.5f) * 180.0f;
        float pn = (float(rand()) / RAND_MAX - 0.5f) * 180.0f;
        float yn = (float(rand()) / RAND_MAX - 0.5f) * 360.0f;

        float ra = rn, pa = pn, ya = yn;
        nautical_to_aero(ra, pa, ya);
        aero_to_nautical(ra, pa, ya);

        assert_close(ra, rn, tol_angle, "Random roll");
        assert_close(pa, pn, tol_angle, "Random pitch");
        assert_close(ya, yn, tol_angle, "Random yaw");
    }

    // Quaternion -> Euler (nautical) direct test from aerospace quaternion input.
    {
        const float rn = 30.0f, pn = 20.0f, yn = 45.0f;
        float ra = rn, pa = pn, ya = yn;
        nautical_to_aero(ra, pa, ya);

        const Quaternionf q_bw_ned = quat_from_euler(ra, pa, ya);
        float r_e, p_e, y_e;
        quat_to_euler_nautical(q_bw_ned, r_e, p_e, y_e);
        assert_close(r_e, rn, tol_angle, "quat_to_euler_nautical roll");
        assert_close(p_e, pn, tol_angle, "quat_to_euler_nautical pitch");
        assert_close(y_e, yn, tol_angle, "quat_to_euler_nautical yaw");
    }

    // Magnetometer world-field tests
    Vector3f mag_enu = MagSim_WMM::mag_world_nautical();
    Vector3f mag_ned = MagSim_WMM::mag_world_aero();
    assert_close(
        mag_enu.z(),
        -std::sin(MagSim_WMM::default_inclination_deg * std::numbers::pi_v<float> / 180.0f)
            * MagSim_WMM::default_total_field_uT,
        1e-3f,
        "Mag ENU vertical");
    assert_close(
        mag_ned.z(),
        std::sin(MagSim_WMM::default_inclination_deg * std::numbers::pi_v<float> / 180.0f)
            * MagSim_WMM::default_total_field_uT,
        1e-3f,
        "Mag NED vertical");

    // Magnetometer body-frame test at zero Euler (should match world ENU)
    Vector3f mag_body0 = MagSim_WMM::simulate_mag_from_euler_nautical(0.0f, 0.0f, 0.0f);
    assert_close(mag_body0.x(), mag_enu.x(), 1e-3f, "Mag body0 east");
    assert_close(mag_body0.y(), mag_enu.y(), 1e-3f, "Mag body0 north");
    assert_close(mag_body0.z(), mag_enu.z(), 1e-3f, "Mag body0 up");

    // Yaw by 90 deg preserves horizontal magnitude.
    Vector3f mag_body_yaw90 = MagSim_WMM::simulate_mag_from_euler_nautical(0.0f, 0.0f, 90.0f);
    assert_close(mag_body0.head<2>().norm(),
                 mag_body_yaw90.head<2>().norm(),
                 1e-3f,
                 "Mag yaw90 preserves horizontal norm");

    // Yaw by 180 deg flips the horizontal vector.
    Vector3f mag_body_yaw180 = MagSim_WMM::simulate_mag_from_euler_nautical(0.0f, 0.0f, 180.0f);
    assert_close(mag_body_yaw180.x(), -mag_body0.x(), 1e-3f, "Mag yaw180 east");
    assert_close(mag_body_yaw180.y(), -mag_body0.y(), 1e-3f, "Mag yaw180 north");
    assert_close(mag_body_yaw180.z(),  mag_body0.z(), 1e-3f, "Mag yaw180 up");

    // Cross-check nautical and aerospace mag simulators for the same physical attitude.
    {
        const float rn = 30.0f, pn = 20.0f, yn = 45.0f;
        float ra = rn, pa = pn, ya = yn;
        nautical_to_aero(ra, pa, ya);

        const Vector3f m_from_n = MagSim_WMM::simulate_mag_from_euler_nautical(rn, pn, yn);
        const Vector3f m_from_a = MagSim_WMM::simulate_mag_from_euler_aero(ra, pa, ya);
        assert_close(m_from_n.x(), m_from_a.x(), 1e-3f, "Mag nautical/aero x");
        assert_close(m_from_n.y(), m_from_a.y(), 1e-3f, "Mag nautical/aero y");
        assert_close(m_from_n.z(), m_from_a.z(), 1e-3f, "Mag nautical/aero z");
    }

    std::cout << "All frame conversion tests passed\n";
    return 0;
}

static inline Quaternionf quat_wb_zu_from_csv(float w, float x, float y, float z) {
    Quaternionf q(w, x, y, z);
    if (q.norm() < 1e-6f) return Quaternionf::Identity();
    return q.normalized();
}

static inline void quat_wb_zu_to_euler_nautical(const Quaternionf& q_wb_zu,
                                                float& roll,
                                                float& pitch,
                                                float& yaw)
{
    const Matrix3f C_wb_zu = q_wb_zu.normalized().toRotationMatrix();
    matrix_to_euler_zyx_deg(C_wb_zu, roll, pitch, yaw);
}

static inline Vector3f mag_body_from_quat_wb_zu(
    const Quaternionf& q_wb_zu,
    float declination_deg = MagSim_WMM::default_declination_deg,
    float inclination_deg = MagSim_WMM::default_inclination_deg,
    float total_field_uT  = MagSim_WMM::default_total_field_uT)
{
    const Matrix3f C_wb_zu = q_wb_zu.normalized().toRotationMatrix();
    const Vector3f mag_world_zu =
        MagSim_WMM::mag_world_nautical(declination_deg, inclination_deg, total_field_uT);
    return C_wb_zu * mag_world_zu;
}

#endif // FRAMECONV_TEST
