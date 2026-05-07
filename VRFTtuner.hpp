#ifndef VRFT_TUNER_HPP
#define VRFT_TUNER_HPP

#include "arm_math.h"
#include <string.h>

class VRFTuner {
public:
    static const int SAMPLES = 2000; // 2秒分に増やして精度を上げる
    float32_t u_data[SAMPLES];
    float32_t y_data[SAMPLES];
    int count = 0;

    float32_t Td = 0.05f; 
    float32_t Ts = 0.001f;

    void record(float u, float y) {
        if (count < SAMPLES) {
            u_data[count] = u;
            y_data[count] = y;
            count++;
        }
    }

   bool solve(float& kp, float& ki, float& kd) {
    if (count < SAMPLES) return false;

    const int N = SAMPLES - 4;
    const float a1 = expf(-Ts / Td);
    const float b0 = 1.0f - a1;

    static float32_t phi_data[SAMPLES * 3];
    static float32_t u_target[SAMPLES];
    static float32_t pT_d[3 * SAMPLES];
    static float32_t ptp_d[9];
    static float32_t ptpI_d[9];
    static float32_t ptu_d[3];
    float32_t theta_res[3];

    memset(phi_data, 0, sizeof(phi_data));

    float fy_prev     = y_data[0];
    float integral_ev = 0.0f;

    for (int i = 2; i < N + 2; i++) {
        int idx = i - 2;

        // 🔴 修正1：正しい仮想誤差
        float fy = a1 * fy_prev + b0 * y_data[i];
        float ev = y_data[i] - fy;
        fy_prev  = fy;

        integral_ev += ev * Ts;

        // 🟡 修正3：5点SG微分
        float dev = (-y_data[i+2] + 8.0f*y_data[i+1]
                     - 8.0f*y_data[i-1] + y_data[i-2])
                    / (12.0f * Ts);

        phi_data[idx*3 + 0] = ev;
        phi_data[idx*3 + 1] = integral_ev;
        phi_data[idx*3 + 2] = dev;
        u_target[idx]        = u_data[i];
    }

    arm_matrix_instance_f32 Phi, U, PhiT, PhiTPhi, PhiTPhiInv, PhiTU, Theta;
    arm_mat_init_f32(&Phi,        N, 3, phi_data);
    arm_mat_init_f32(&U,          N, 1, u_target);
    arm_mat_init_f32(&PhiT,       3, N, pT_d);
    arm_mat_init_f32(&PhiTPhi,    3, 3, ptp_d);
    arm_mat_init_f32(&PhiTPhiInv, 3, 3, ptpI_d);
    arm_mat_init_f32(&PhiTU,      3, 1, ptu_d);
    arm_mat_init_f32(&Theta,      3, 1, theta_res);

    arm_mat_trans_f32(&Phi, &PhiT);
    arm_mat_mult_f32(&PhiT, &Phi, &PhiTPhi);
    if (arm_mat_inverse_f32(&PhiTPhi, &PhiTPhiInv) != ARM_MATH_SUCCESS)
        return false;
    arm_mat_mult_f32(&PhiT, &U, &PhiTU);
    arm_mat_mult_f32(&PhiTPhiInv, &PhiTU, &Theta);

    kp = theta_res[0];
    ki = theta_res[1];
    kd = theta_res[2];

    // Fit Rate で推定精度を数値化
    float ss_res = 0, ss_tot = 0, u_mean = 0;
    for (int i = 0; i < N; i++) u_mean += u_target[i];
    u_mean /= (float)N;
    for (int i = 0; i < N; i++) {
        float u_hat = phi_data[i*3+0]*kp
                    + phi_data[i*3+1]*ki
                    + phi_data[i*3+2]*kd;
        ss_res += (u_target[i] - u_hat) * (u_target[i] - u_hat);
        ss_tot += (u_target[i] - u_mean) * (u_target[i] - u_mean);
    }
    float fit = (ss_tot > 1e-6f)
              ? (1.0f - sqrtf(ss_res) / sqrtf(ss_tot)) * 100.0f
              : 0.0f;
    printf("Fit Rate: %d%%\n", (int)fit);

    if (!isfinite(kp)||!isfinite(ki)||!isfinite(kd)) return false;
    if (kp <= 0.0f || ki < 0.0f || kd < 0.0f)       return false;
    if (kp > 5000.0f || ki > 5000.0f || kd > 500.0f) return false;
    if (fit < 60.0f)                                  return false;

    return true;
}
};
#endif
