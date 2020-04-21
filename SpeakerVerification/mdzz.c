#include<stdio.h>
#include<math.h>
#include "MyMatrix.h"
#include "mfcc.h"
#pragma warning(disable:4996)
#define PI 3.1415926535
#define NFFT 512
#define CEPLIFTER 22

//-----------------------------------------------------------fft--------------
typedef struct {       //定义一个结构体表示复数的类型
	float real;
	float imag;
}complex;
complex x[NFFT], *W;   //定义输入序列和旋转因子
int size = 512;   //定义数据长度
void output()
{
	int i;
	for (i = 0; i < size; i++)
	{
		printf("%.4f", x[i].real);
		if (x[i].imag >= 0.0001)
		{
			printf("+%.4fj\n", x[i].imag);
		}
		else if (fabs(x[i].imag) < 0.0001)
		{
			printf("\n");
		}
		else
		{
			printf("%.4fj\n", x[i].imag);
		}
	}
}
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
//--------------------------------------------------------------------------

float* preemphasis(short *signal, float preemph, int size) {
	float* new_signal = NULL;
	new_signal = (float*)malloc(sizeof(float) * size);
	int i = 0;
	*new_signal = *signal;
	for (i = 1; i < size; i++) {
		(float)(*(new_signal + i)) = (float)(*(signal + i)) - (float)(*(signal + i - 1)) * preemph;
	}
	printf("1: %.2f\t2: %.2f\n", new_signal[0], new_signal[1]);
	return new_signal;
}

float* framesig(float* signal, float frame_len, float frame_step, float* winfunc, int slen) {
	//slen为signal长度
	int framelen = (int)(frame_len + 0.5);//帧长 16000*0.025
	int framestep = (int)(frame_step + 0.5);//帧移动 16000*0.01
	int numframes, padlen; //帧数，总长度

	if (slen <= framelen)
		numframes = 1;
	else
		numframes = 1 + (int)(ceil((1.0 * slen - framelen) / framestep));//帧数，95

	padlen = (int)((numframes - 1) * framestep + framelen); //总长度15440

	float* padsignal = (float*)malloc(sizeof(float) * (padlen));
	float* psignal = (float*)malloc(sizeof(float) * (numframes * framelen));
	for (int i = 0; i < slen; i++) {
		*(padsignal + i) = *(signal + i);
	}

	for (int i = slen; i < padlen; i++) { //一开始是15360后面补0
		*(padsignal + i) = 0;
	}

	for (int i = 0; i < numframes; i++) { //对每帧进行加窗
		for (int j = 0; j < framelen; j++) {
			*(psignal + i * framelen + j) = *(padsignal + i * framestep + j);
		}
	}

	return psignal;
}

float* framesig2(float* signal, float frame_len, float frame_step, float* winfunc, int slen, int* frame_num) {
	//slen为signal长度
	int framelen = (int)(frame_len + 0.5);//帧长 16000*0.025
	int framestep = (int)(frame_step + 0.5);//帧移动 16000*0.01
	int numframes, padlen; //帧数，总长度

	if (slen <= framelen)
		numframes = 1;
	else
		numframes = 1 + (int)(ceil((1.0 * slen - framelen) / framestep));//帧数，95
	*frame_num = numframes;
	padlen = (int)((numframes - 1) * framestep + framelen); //总长度15440

	float* padsignal = (float*)malloc(sizeof(float) * (padlen));
	float* psignal = (float*)malloc(sizeof(float) * (numframes * framelen));
	for (int i = 0; i < slen; i++) {
		*(padsignal + i) = *(signal + i);
	}

	for (int i = slen; i < padlen; i++) { //一开始是15360后面补0
		*(padsignal + i) = 0;
	}

	for (int i = 0; i < numframes; i++) { //对每帧进行加窗
		for (int j = 0; j < framelen; j++) {
			*(psignal + i * framelen + j) = *(padsignal + i * framestep + j);
		}
	}

	return psignal;
}

float* dct(float* feat, int numframes) {
	float* f = (float*)malloc(sizeof(float) * numframes * 13);
	int i,k,n;
	float A, s;
	for ( i = 0; i < numframes; i++)
	{
		for (k = 0; k < 13; k++)
		{
			s = 0;
			if (k == 0)
				A = sqrt(1.0 / 26); //计算k=0时的系数
			else
				A = sqrt(2.0 / 26); //计算k!=0时的系数
			for (n = 0; n < 26; n++)
			{
				float tmp = (*(feat+i*26+n)) * cos((PI*(2 * n + 1)*k) / (2 * 26));
				s = s + tmp;	//累加求和
			}
			(*(f + i * 13 + k)) = A * s;	//X[k]等于累和结果s乘以系数A
		}
	}
	return f;
}

float* lifter(float* feat, int frame_num, int ceplifter) {
	int nframes = frame_num;
	float lift[13]={1.00000000,2.56546330,4.09905815,5.56956530,6.94704914,8.20346832,9.31324577,10.2537889,11.0059519,11.5544224,11.8880358,12.0000000,11.8880358};
	int i, j;
	for (i = 0; i < nframes; i++) {
		for (j = 0; j < 13; j++)
		{
			*(feat + i * 13 + j) = (*(feat + i * 13 + j)) * lift[j];

		}
	}
	return feat;
}

float* appendEnergy(float* feat, float* energy, int frame_num) {
	int i;
	for (i = 0; i < frame_num; i++)
	{
		*(feat + i * 13) = log(*(energy + i));
	}
	return feat;
}


void main() {
	printf("program start!\n");
	FILE* fp = NULL;
	int count, fsize = 0;
	short* sig;

	//打开文件
	fp = fopen("testwav", "rb");

	//计算所需要的空间
	fseek(fp, 0, SEEK_END);
	fsize = ftell(fp) / 2;
	printf("%d\n", fsize);
	rewind(fp);

	//动态分配内存
	sig = (short*)malloc(sizeof(short) * fsize);

	//循环读取到内存（堆）
	while (!feof(fp)) {
		count = fread(sig, sizeof(short), fsize, fp);
		int n = feof(fp);
	}
	fclose(fp);
	//printf("%d\n", sig[0]);
	//===================================================================
	//参数
	int samplerate = 16000;
	float winlen = 0.025;
	float winstep = 0.01;
	float* winfunc = NULL;

	//预加重==============================================================
	float* pre_sig;
	pre_sig = preemphasis(sig, 0.97, fsize);

	//加窗分帧============================================================
	float* frame;
	int frame_num;
	frame = framesig2(pre_sig, winlen * samplerate, winstep * samplerate, winfunc, fsize, &frame_num);

	//快速傅里叶变化=======================================================
	int i, j, winlenth = 0;
	float* complex_spec = NULL;
	float abs_fft;
	float* energy = NULL;
	float sum = 0.0;

	//float temp[NFFT];
	winlenth = winlen * samplerate;
	complex_spec = (float*)malloc(sizeof(float) * (frame_num*(NFFT / 2 + 1)));
	energy = (float*)malloc(sizeof(float) * frame_num);
	for (i = 0; i < frame_num; i++)
	{
		if (winlenth < NFFT) {
			for (j = 0; j < winlenth; j++) {
				//temp[j] = *(frame + i * winlenth + j);
				x[j].real = *(frame + i * winlenth + j);
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
				x[j].real = *(frame + i * winlenth + j);
				x[j].imag = 0.0;
			}
		}
		//rfft(temp, NFFT);
		transform();//变换序列顺序
		fft();//蝶形运算
		//printf("输出FFT后的结果\n");
		//output();
		//printf("\n第%d帧\n", i);
		for (j = 0; j < NFFT / 2 + 1; j++)
		{
			//*(complex_spec + j + i * (NFFT / 2 + 1)) = *(temp + j);
			abs_fft = sqrt(x[j].real * x[j].real + x[j].imag * x[j].imag);
			//printf("<%d>%6f|\t", j, abs_fft*abs_fft/NFFT);
			*(complex_spec + j + i * (NFFT / 2 + 1)) = abs_fft * abs_fft / NFFT;
			sum += abs_fft * abs_fft / NFFT;
		}
		//printf("%d\t%f\n", i, sum);
		*(energy + i) = sum;
		sum = 0.0;
	}

	//读取mel滤波
	Matrix A, A_t, B, FEAT;
	A = Create_Matrix(26, 257);
	SetData_Matrix(A, mel);
	printf("A(1, 1) = %f, A(1, 2) = %f, A(1, 3) = %f\n", PickInMat(A, 1, 1), PickInMat(A, 1, 2), PickInMat(A, 1, 3));
	A_t = Trans_Matrix(A);
	printf("row = %d, col = %d\n", A_t->row, A_t->column);
	B = Create_Matrix(frame_num, 257);
	SetData_Matrix(B, complex_spec);
	printf("B(1, 1) = %f, B(1, 2) = %f, B(1, 3) = %f\n", PickInMat(B, 1, 1), PickInMat(B, 1, 2), PickInMat(B, 1, 3));
	FEAT = Mult_Matrix(B, A_t);
	printf("feat(1, 1) = %f, feat(1, 2) = %f, feat(95, 26) = %f\n", PickInMat(FEAT, 1, 1), PickInMat(FEAT, 1, 2), PickInMat(FEAT, 95, 26));
	printf("FEAT->row = %d\t FEAT->col = %d\n", FEAT->row, FEAT->column);
	Free_Matrix(A);
	Free_Matrix(A_t);
	Free_Matrix(B);

	//取对数
	float* feat;
	float* feat2;
	feat = (float*)malloc(sizeof(float) * frame_num * 26);
	feat2= (float*)malloc(sizeof(float) * frame_num * 13);
	for (i = 1; i <= frame_num; i++)
	{
		for (j = 1; j <= 26; j++)
		{
			*(feat + (i - 1) * 26 + (j - 1)) = (float)log((double)PickInMat(FEAT, i, j));
		}
	}
	//dct
	feat2 = dct(feat, frame_num);


	//用倒谱提升器，得到倒谱的矩阵。这有增加高频DCT系数的作用
	feat2 = lifter(feat2, frame_num, CEPLIFTER);
	//用帧能量的对数替换第一个倒谱系数
	feat2 = appendEnergy(feat2, energy, frame_num);
	printf("program end!\n");
	//mfcc = mfcc(signal, rate);
	return;
}


float* mfcc(short* signal, int samplerate, float winlen, float winstep, int numcep) {

	return signal;
}
