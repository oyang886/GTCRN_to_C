#include <stdio.h>
#include <stdlib.h>
#include "kiss_fft.h"
#include "test_data.h"

#define BLOCK 512
#define HALF BLOCK/2
#define DEBUG_MODE 1

// 时域带噪音频
float td_noisy[HALF];
// 上一帧的带噪音频
float last_x1[HALF];
// 组合起来加窗分帧
float xx1[BLOCK];
float last_output[HALF];
cfloat X1[BLOCK];
float input_real[HALF + 1];
float input_imag[HALF + 1];

typedef struct
{
    kiss_fft_state* kfft;
} CommonState;

CommonState common;

static void forward_transform(cfloat* out, const float* in)
{
    int i;
    cfloat x[NFFT];
    for (i = 0; i < NFFT; i++) {
        x[i].r = in[i];
        x[i].i = 0.0F;
    }
    opus_fft(common.kfft, x, out, 0);
}

static void inverse_transform(float* out, const cfloat* in)
{
    int i;
    cfloat x[NFFT];
    cfloat y[NFFT];
    for (i = 0; i < (NFFT >> 1) + 1; i++)
        x[i] = in[i];
    for (i = (NFFT >> 1) + 1; i < NFFT; i++) {
        x[i].r = x[NFFT - i].r;
        x[i].i = -x[NFFT - i].i;
    }
    opus_fft(common.kfft, x, y, 0);
    /* output in reverse order for IFFT. */
    out[0] = NFFT * y[0].r;
    for (i = 1; i < NFFT; i++)
        out[i] = NFFT * y[NFFT - i].r;
}

// 启动先初始化fft
void init_fft()
{
    int i;
    common.kfft = opus_fft_alloc_twiddles(NULL, NULL, NULL, 0);
}

void init_parm()
{
    int i = 0;
    for (i = 0; i < BLOCK; i++)
    {
        X1[i].r = 0.0f;
        X1[i].i = 0.0f;
    }
    memset(last_x1, 0, sizeof(last_x1));
    memset(last_output, 0, sizeof(last_output));
}

// 加窗分帧
void apply_window(float *x, int LEN)
{
    int i = 0;
    for (i = 0; i < LEN; i++)
    {
        x[i] *= hann_window[i];
    }
}

void test_function()
{
    int i;
    memcpy(&xx1[0], &last_x1[0], HALF * sizeof(float));
    memcpy(&xx1[HALF], &td_noisy[0], HALF * sizeof(float));
    // 保存此帧到下一帧用
    memcpy(last_x1, td_noisy, sizeof(td_noisy));
    apply_window(xx1, BLOCK);
    // fft
    forward_transform(X1, xx1);
    for (i = 0; i < HALF + 1; i++)
    {
        input_real[i] = X1[i].r;
        input_imag[i] = X1[i].i;
    }

#if DEBUG_MODE
    // 对比测试fft准确性
    for (i = 0; i < HALF + 1; i++)
    {
        printf("real%d: %f\n", i + 1, input_real[i]);
    }
    printf("\n");
    for (i = 0; i < HALF + 1; i++)
    {
        printf("imag%d: %f\n", i + 1, input_imag[i]);
    }
#endif
    
}

int main()
{
    init_fft();
    init_parm();
    test_function();
}
