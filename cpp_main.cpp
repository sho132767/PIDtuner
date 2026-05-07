#include "Config.hpp"
#include "Encoder.hpp"
#include "Motor.hpp"
#include "pid.hpp"
#include "main.h"
#include "Timer.hpp"
#include "math.h"
#include "stm32f4xx_hal.h"
#include "VRFTtuner.hpp"
#include <cstdint>
#include <stdio.h>
#include <stdlib.h>


constexpr double dt = 0.005;

uint32_t shift_reg = 0xACE1u;
float get_m_sequence() {
    uint32_t bit = ((shift_reg>>0)^(shift_reg>>2)
                   ^(shift_reg>>3)^(shift_reg>>5)) & 1;
    shift_reg = (shift_reg >> 1) | (bit << 15);
    return (bit ? 1.0f : -1.0f);
}

int cpp_main()
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
    printf("resolution = %lu\n", (unsigned long)resolution); // 要確認
    mot.set_limit(resolution * 1.0, resolution * -1.0);

    VRFTuner tuner;
    tuner.Td = 0.05f;
    tuner.Ts = (float)dt;  // 0.001f

    printf("\n--- VRFT Auto Tuning System Start ---\n");
    HAL_Delay(3000);

    // ---------------------------------------------------------
    // Phase 1: データ収集
    // ---------------------------------------------------------
    printf("Phase 1: Recording...\n");
   

    float amplitude = 4000.0f;
    int clock_div   = 0;
    float prbs_val  = 1.0f;
    int last_printed = -1;  // 🔴 修正2：同じcountで連続表示しない

    double t = timer.read();

    while (tuner.count < VRFTuner::SAMPLES) {
        double current_t = timer.read();
       
        if (current_t - t < dt) {
            t = current_t;

            enc.update();
            float omega = (float)enc.get_omega();

            if (++clock_div >= 10) {
                clock_div = 0;
                prbs_val = get_m_sequence();
            }
            float u = prbs_val * amplitude;
            mot.set_value((int32_t)u);
            tuner.record(u, omega);

           
            if (tuner.count % 100 == 0 && tuner.count != last_printed) {
                last_printed = tuner.count;
                printf("Data: %d/%d\n", tuner.count, VRFTuner::SAMPLES);
            }
            
            

        }
    }
    mot.set_value(0);
    printf("Phase 1: Done.\n");

    // ---------------------------------------------------------
    // Phase 2: PIDゲイン計算
    // ---------------------------------------------------------
    printf("Phase 2: Calculating...\n");
    float kp_new, ki_new, kd_new;

    if (tuner.solve(kp_new, ki_new, kd_new)) {
        printf(">>> Tuning Success! <<<\n");
        printf("Kp:%d Ki:%d Kd:%d (x0.001)\n",
               (int)(kp_new*1000),
               (int)(ki_new*1000),
               (int)(kd_new*1000));
    } else {
        printf(">>> Tuning Failed. Default gains. <<<\n");
        kp_new = 15.0f; ki_new = 1.0f; kd_new = 0.0f;
    }

    Vel_PID pid(kp_new, ki_new, kd_new);
    pid.set_limits(resolution * -1.0, resolution * 1.0);

    printf("Phase 3: Start (Target: 40.0 rad/s)\n");
    HAL_Delay(2000);

    // ---------------------------------------------------------
    // Phase 3: 速度制御
    // ---------------------------------------------------------
    timer.reset();
    timer.start();
    t = timer.read();  

    while (true) {
        
           

            enc.update();
            double omega = enc.get_omega();
            double target_omega = 40.0;

            pid.Input(target_omega, omega);
            mot.set_value(static_cast<int32_t>(pid.Output()));

            printf("%d, %d\n",
                   (int)(omega * 100),
                   (int)(target_omega * 100));
                while(timer.read()-t<dt);
                 t = timer.read();
                 
     }
           
    return 0;
}
