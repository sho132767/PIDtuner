#include "Config.hpp"
#include "Encoder.hpp"
#include "Motor.hpp"
#include "pid.hpp"
#include "main.h"
#include "Timer.hpp"
#include "math.h"
#include "stm32f4xx_hal.h"
#include "VRFTturner.hpp"
#include <cstdint>
#include <stdio.h>
#include <stdlib.h>

// ============================================================
// サンプリング周期（100Hz）
// ============================================================
constexpr float dt = 0.01f;

// ============================================================
// 15ビット PRBS（最大長系列）
// ============================================================
uint32_t shift_reg = 0xACE1u;
float get_m_sequence() {
    uint32_t bit = ((shift_reg>>0) ^ (shift_reg>>2)
                  ^ (shift_reg>>3) ^ (shift_reg>>5)) & 1u;
    shift_reg = (shift_reg >> 1) | (bit << 14);  // 15bit
    return (bit ? 1.0f : -1.0f);
}

// ============================================================
int cpp_main()
// ============================================================
{
    main_timer::activate();

    Timer timer;
    timer.reset();
    timer.start();

    Encoder_SINGLE_Interrupt enc(GPIOB, GPIO_PIN_8,
                                 GPIOC, GPIO_PIN_11, 2048, NEGATIVE);

    Motor mot(&htim1, TIM_CHANNEL_4,
              &htim2, TIM_CHANNEL_4, POSITIVE, NEGATIVE);

    uint32_t resolution = mot.get_resolution();
    printf("resolution = %lu\n", (unsigned long)resolution);
    mot.set_limit((double)resolution, -(double)resolution);

    // --------------------------------------------------------
    // VRFTuner 設定
    //   wc = 目標帯域の 1/2 〜 1/3 程度から始める
    //   目標 40 rad/s → wc = 15 rad/s
    // --------------------------------------------------------
    VRFTuner tuner;
    tuner.Ts = dt;
    tuner.wc = 15.0f;   // ← fit が低ければ 5〜30 の範囲で調整

    printf("\n--- VRFT Auto Tuning System Start ---\n");
    printf("Ts=%.3f wc=%.1f SAMPLES=%d\n", tuner.Ts, tuner.wc, VRFTuner::SAMPLES);
    HAL_Delay(3000);

    // ========================================================
    // Phase 1: PRBS 励振 → データ収集
    // ========================================================
    printf("Phase 1: Recording...\n");

    const float amplitude = 4000.0f;
    const int   prbs_hold = 5;    // 5サンプル（50ms）ごとに切替
    int   clock_div  = 0;
    float prbs_val   = 1.0f;
    int   last_print = -1;

    float t = timer.read();

    while (tuner.count < VRFTuner::SAMPLES) {

        // --- エンコーダ更新 ---
        enc.update();
        float omega = (float)enc.get_omega();

        // --- PRBS 切替 ---
        if (++clock_div >= prbs_hold) {
            clock_div = 0;
            prbs_val  = get_m_sequence();
        }

        // --- モータ出力 & 記録 ---
        float u = prbs_val * amplitude;
        mot.set_value((int32_t)u);
        tuner.record(u, omega);

        // --- 進捗表示 ---
        if (tuner.count % 200 == 0 && tuner.count != last_print) {
            last_print = tuner.count;
            printf("Data: %d/%d  omega=%.2f\n",
                   tuner.count, VRFTuner::SAMPLES, omega);
        }

        // --- 正確な周期待ち ---
        while (timer.read() - t < dt);
        t = timer.read();
    }

    mot.set_value(0);
    printf("Phase 1: Done.\n");
    HAL_Delay(500);

    // ========================================================
    // Phase 2: VRFTゲイン計算
    // ========================================================
    printf("Phase 2: Calculating...\n");

    float kp_new, ki_new, kd_new;

    if (tuner.solve(kp_new, ki_new, kd_new)) {
        printf(">>> Tuning Success! <<<\n");
        printf("Kp=%.4f Ki=%.4f Kd=%.4f\n", kp_new, ki_new, kd_new);
    } else {
        printf(">>> Tuning Failed. Default gains. <<<\n");
        kp_new = 15.0f;
        ki_new =  1.0f;
        kd_new =  0.0f;
    }

    Vel_PID pid(kp_new, ki_new, kd_new);
    pid.set_limits(-(double)resolution, (double)resolution);

    printf("Phase 3: Start (Target: 40.0 rad/s)\n");
    HAL_Delay(2000);

    // ========================================================
    // Phase 3: 速度制御
    // ========================================================
    timer.reset();
    timer.start();
    t = timer.read();

    const double target_omega = 40.0;

    while (true) {
        enc.update();
        double omega = enc.get_omega();

        pid.Input(target_omega, omega);
        mot.set_value(static_cast<int32_t>(pid.Output()));

        while (timer.read() - t < dt);
        t = timer.read();
    }

    return 0;
}
