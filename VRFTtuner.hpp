#ifndef VRFT_TUNER_HPP
#define VRFT_TUNER_HPP

#include "arm_math.h"
#include <string.h>
#include <stdio.h>
#include <math.h>

class VRFTuner {
public:
    // -------------------------------------------------------
    // 定数
    // -------------------------------------------------------
    static const int SAMPLES = 2000;   // 100Hz × 20秒

    // -------------------------------------------------------
    // データバッファ
    // -------------------------------------------------------
    float32_t u_data[SAMPLES];
    float32_t y_data[SAMPLES];
    int       count = 0;

    // -------------------------------------------------------
    // パラメータ（main から設定する）
    // -------------------------------------------------------
    float Ts   = 0.01f;   // サンプリング周期 [s]
    float wc   = 40.0f;   // 参照モデルのカットオフ [rad/s]
                          // 目標帯域に合わせる（40rad/s目標なら10〜20が目安）

    // Td は削除（wc で直接指定する設計に変更）

    // -------------------------------------------------------
    // データ記録
    // -------------------------------------------------------
    void record(float u, float y) {
        if (count < SAMPLES) {
            u_data[count] = u;
            y_data[count] = y;
            count++;
        }
    }

    // -------------------------------------------------------
    // VRFT ソルバ
    // -------------------------------------------------------
    bool solve(float& kp, float& ki, float& kd) {

        // 有効サンプル数（前後2点を微分に使うので端を除く）
        const int N = count - 4;
        if (N < 50) {
            printf("ERR: Not enough samples (%d)\n", count);
            return false;
        }

        // --------------------------------------------------
        // 参照モデル M(z)：1次ローパス
        //   連続: M(s) = wc / (s + wc)
        //   離散: ym[k] = a1*ym[k-1] + b0*y[k]
        // --------------------------------------------------
        const float a1 = expf(-wc * Ts);
        const float b0 = 1.0f - a1;
        printf("VRFT: wc=%.1f Ts=%.4f a1=%.4f b0=%.4f\n", wc, Ts, a1, b0);

        // --------------------------------------------------
        // 静的バッファ（スタック節約）
        // --------------------------------------------------
        static float32_t phi_data[SAMPLES * 3];
        static float32_t u_target[SAMPLES];
        static float32_t pT_d   [3 * SAMPLES];
        static float32_t ptp_d  [9];
        static float32_t ptpI_d [9];
        static float32_t ptu_d  [3];
        float32_t         theta_res[3];

        memset(phi_data, 0, sizeof(float32_t) * N * 3);

        // --------------------------------------------------
        // 基底関数の構築
        //
        //  仮想誤差 ev[k] = M(z)*y[k] - y[k]
        //                 = ym[k] - y[k]
        //
        //  φ(k) = [ ev(k),  ∫ev·dt,  dev/dt ]
        //  u_target(k) = u(k)
        // --------------------------------------------------
        float ym_prev     = y_data[0];
        float integral_ev = 0.0f;

        for (int i = 2; i < N + 2; i++) {
            int idx = i - 2;

            // 参照モデル出力
            float ym = a1 * ym_prev + b0 * y_data[i];
            ym_prev  = ym;

            // 仮想誤差（目標応答との差）
            float ev = ym - y_data[i];

            // 積分
            integral_ev += ev * Ts;

            // 5点 Savitzky-Golay 微分（境界は2点確保済み）
            float dev = (-y_data[i+2] + 8.0f*y_data[i+1]
                         - 8.0f*y_data[i-1] + y_data[i-2])
                        / (12.0f * Ts);

            phi_data[idx*3 + 0] = ev;
            phi_data[idx*3 + 1] = integral_ev;
            phi_data[idx*3 + 2] = dev;
            u_target[idx]        = u_data[i];
        }

        // --------------------------------------------------
        // デバッグ：基底関数の統計を出力
        // --------------------------------------------------
        float ev_max=0, intev_max=0, dev_max=0;
        for (int i = 0; i < N; i++) {
            float e = fabsf(phi_data[i*3+0]);
            float n = fabsf(phi_data[i*3+1]);
            float d = fabsf(phi_data[i*3+2]);
            if (e > ev_max)    ev_max    = e;
            if (n > intev_max) intev_max = n;
            if (d > dev_max)   dev_max   = d;
        }
        printf("phi_max: ev=%.2f intev=%.2f dev=%.2f\n",
               ev_max, intev_max, dev_max);

        // --------------------------------------------------
        // 最小二乗法：θ = (Φ^T Φ)^{-1} Φ^T u
        // --------------------------------------------------
        arm_matrix_instance_f32 Phi, U, PhiT, PhiTPhi, PhiTPhiInv, PhiTU, Theta;
        arm_mat_init_f32(&Phi,        N, 3, phi_data);
        arm_mat_init_f32(&U,          N, 1, u_target);
        arm_mat_init_f32(&PhiT,       3, N, pT_d);
        arm_mat_init_f32(&PhiTPhi,    3, 3, ptp_d);
        arm_mat_init_f32(&PhiTPhiInv, 3, 3, ptpI_d);
        arm_mat_init_f32(&PhiTU,      3, 1, ptu_d);
        arm_mat_init_f32(&Theta,      3, 1, theta_res);

        arm_mat_trans_f32(&Phi, &PhiT);
        arm_mat_mult_f32(&PhiT, &Phi,  &PhiTPhi);
        if (arm_mat_inverse_f32(&PhiTPhi, &PhiTPhiInv) != ARM_MATH_SUCCESS) {
            printf("ERR: Matrix singular\n");
            return false;
        }
        arm_mat_mult_f32(&PhiT,      &U,      &PhiTU);
        arm_mat_mult_f32(&PhiTPhiInv, &PhiTU,  &Theta);

        kp = theta_res[0];
        ki = theta_res[1];
        kd = theta_res[2];

        printf("Raw gains: Kp=%d Ki=%d Kd=%d\n", (int)kp*100, (int)ki*100,(int) kd*100);

        // --------------------------------------------------
        // Fit Rate（NRMSE ベース）
        // --------------------------------------------------
        float ss_res = 0.0f, ss_tot = 0.0f, u_mean = 0.0f;
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

        // --------------------------------------------------
        // 妥当性チェック
        // --------------------------------------------------
        if (!isfinite(kp) || !isfinite(ki) || !isfinite(kd)) {
            printf("ERR: Non-finite gains\n");
            return false;
        }
        if (fit < 30.0f) {
            printf("ERR: Fit too low\n");
            return false;
        }
        if (kp <= 0.0f || ki < 0.0f || kd < 0.0f) {
            printf("ERR: Negative gains\n");
            return false;
        }
        if (kp > 50000.0f || ki > 50000.0f || kd > 5000.0f) {
            printf("ERR: Gains too large\n");
            return false;
        }

        return true;
    }
};

#endif
