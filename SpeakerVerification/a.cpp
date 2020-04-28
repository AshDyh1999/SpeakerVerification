#include <stdlib.h>
#include <string.h>
#include<stdio.h>
#include<math.h>
#pragma warning(disable:4996)
#define PI 3.1415926535
#define NFFT 512
#define CEPLIFTER 22

//fft相关函数
typedef struct {       //定义一个结构体表示复数的类型
    float real;
    float imag;
}complex;
complex x[NFFT], * W;   //定义输入序列和旋转因子
int size = 512;   //定义数据长度

void change()
{
    complex temp;
    unsigned short i = 0, j = 0, k = 0;
    double t;
    for (i = 0; i < size; i++)
    {
        k = i;
        j = 0;
        t = (log(size) / log(2));
        while ((t--) > 0)
        {
            j = j << 1;
            j |= (k & 1);
            k = k >> 1;
        }
        if (j > i)
        {
            temp = x[i];
            x[i] = x[j];
            x[j] = temp;
        }
    }
    // output();
}

void transform()
{
    int i;
    W = (complex*)malloc(sizeof(complex) * size);
    for (i = 0; i < size; i++)
    {
        W[i].real = cos(2 * PI / size * i);
        W[i].imag = -1 * sin(2 * PI / size * i);
    }
}

void add(complex a, complex b, complex* c)
{
    c->real = a.real + b.real;
    c->imag = a.imag + b.imag;
}

void sub(complex a, complex b, complex* c)
{
    c->real = a.real - b.real;
    c->imag = a.imag - b.imag;
}

void mul(complex a, complex b, complex* c)
{
    c->real = a.real * b.real - a.imag * b.imag;
    c->imag = a.real * b.imag + a.imag * b.real;
}

void fft()
{
    int i = 0, j = 0, k = 0, m = 0;
    complex q, y, z;
    change();
    for (i = 0; i < log(size) / log(2); i++)
    {
        m = 1 << i;
        for (j = 0; j < size; j += 2 * m)
        {
            for (k = 0; k < m; k++)
            {
                mul(x[k + j + m], W[size * k / 2 / m], &q);
                add(x[j + k], q, &y);
                sub(x[j + k], q, &z);
                x[j + k] = y;
                x[j + k + m] = z;
            }
        }
    }
}

float* dct(float* feat) {
    float* f = (float*)malloc(sizeof(float) * 13);
    int i, k, n;
    float A, s;
    for (k = 0; k < 13; k++)
    {
        s = 0;
        if (k == 0)
            A = sqrt(1.0 / 26); //计算k=0时的系数
        else
            A = sqrt(2.0 / 26); //计算k!=0时的系数
        for (n = 0; n < 26; n++)
        {
            float tmp = (*(feat + + n)) * cos((PI * (2 * n + 1) * k) / (2 * 26));
            s = s + tmp;	//累加求和
        }
        (*(f + + k)) = A * s;	//X[k]等于累和结果s乘以系数A
    }
    return f;
}

float* lifter(float* feat, int ceplifter) {
    float lift[13] = { 1.00000000,2.56546330,4.09905815,5.56956530,6.94704914,8.20346832,9.31324577,10.2537889,11.0059519,11.5544224,11.8880358,12.0000000,11.8880358 };
    int i, j;
    for (j = 0; j < 13; j++)
    {
        *(feat + i * 13 + j) = (*(feat + i * 13 + j)) * lift[j];

    }
    return feat;
}


float* mfcc_frame(short* signal_frame, int winlenth, float* winfunc, float* last_point) {
	int i = 0;
	float* new_signal = NULL;
	new_signal = (float*)malloc(sizeof(float) * size);
	*new_signal = *signal_frame;
	for (i = 1; i < winlenth; i++) {
		*(new_signal + 1) = (float)(*(signal_frame + i)) - (float)(*(signal_frame + i - 1)) * 0.97;
	}

    //快速傅里叶变化=======================================================
	int j;
	float* complex_spec = NULL;
	float abs_fft;
	float* energy = NULL;
	float sum = 0.0;

    if (winlenth < NFFT) {
        for (j = 0; j < winlenth; j++) {
            //temp[j] = *(frame + i * winlenth + j);
            x[j].real = *(new_signal + i * winlenth + j);
            x[j].imag = 0.0;
        }
        for (j = 0; j < NFFT - winlenth; j++) {
            //temp[winlenth + j] = 0.0;
            x[winlenth + j].real = 0.0;
            x[winlenth + j].imag = 0.0;
        }
    }
    else {
        for (j = 0; j < NFFT; j++) {
            //temp[j] = *(frame + i * winlenth + j);
            x[j].real = *(new_signal + i * winlenth + j);
            x[j].imag = 0.0;
        }
    }
    transform();//变换序列顺序
    fft();//蝶形运算

    float* complex_spec = (float*)malloc(sizeof(float) * ((NFFT / 2 + 1)));
    for (j = 0; j < NFFT / 2 + 1; j++)
    {
        //*(complex_spec + j + i * (NFFT / 2 + 1)) = *(temp + j);
        abs_fft = sqrt(x[j].real * x[j].real + x[j].imag * x[j].imag);
        //printf("<%d>%6f|\t", j, abs_fft*abs_fft/NFFT);
        *(complex_spec + j) = abs_fft * abs_fft / NFFT;
        sum += abs_fft * abs_fft / NFFT;
    }

    //读取mel滤波==========================================================
    float* FEAT;

    //log=============================================================
    float* feat;
    feat = (float*)malloc(sizeof(float)* 26);

    for (j = 1; j <= 26; j++)
    {
        *(feat + + (j - 1)) = (float)log((double)*(FEAT+j));
    }

    //dct=============================================================
    float* feat_mfcc;
    feat_mfcc = dct(feat);

    //用倒谱提升器，得到倒谱的矩阵。这有增加高频DCT系数的作用
    feat_mfcc = lifter(feat_mfcc, CEPLIFTER);
    //用帧能量的对数替换第一个倒谱系数
    feat_mfcc = appendEnergy(feat_mfcc, energy, frame_num);
    printf("mfcc[0] = %f\n", feat_mfcc[0]);



	
}