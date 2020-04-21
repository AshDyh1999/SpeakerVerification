#include "stdio.h"
#include "math.h"
#include "kfft.c"

#define PI 3.1415926535

main()
{
    int i, j;
    double pr[64], pi[64], fr[64], fi[64], t[64];
    for (i = 0; i <= 63; i++)  //生成输入信号
    {
        t[i] = i * 0.001;
        pr[i] = 1.2 + 2.7 * cos(2 * PI * 33 * t[i]) + 5 * cos(2 * PI * 200 * t[i] + PI / 2); pi[i] = 0.0;
    }

    kfft(pr, pi, 64, 6, fr, fi);  //调用FFT函数
    for (i = 0; i < 64; i++)
    {
        printf("%d\t%lf\n", i, pr[i]); //输出结果
    }
}
