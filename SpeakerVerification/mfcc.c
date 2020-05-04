#include<stdio.h>
#include<math.h>
#include "MyMatrix.h"
#include "mfcc.h"
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

//
float* preemphasis(short *signal, float preemph, int size) {
	float* new_signal = NULL;
	new_signal = (float*)malloc(sizeof(float) * size);
	int i = 0;
	*new_signal = *signal;
	for (i = 1; i < size; i++) {
		(float)(*(new_signal + i)) = (float)(*(signal + i)) - (float)(*(signal + i - 1)) * preemph;
	}
	printf("1: %.2f\t2: %.2f\n", new_signal[0], new_signal[1]);
    //free(new_signal);
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
    free(padlen);
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
    int i, k, n;
    float A, s;
    for (i = 0; i < numframes; i++)
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
                float tmp = (*(feat + i * 26 + n)) * cos((PI * (2 * n + 1) * k) / (2 * 26));
                s = s + tmp;	//累加求和
            }
            (*(f + i * 13 + k)) = A * s;	//X[k]等于累和结果s乘以系数A
        }
    }
    return f;
}

float* dct_frame(float* feat) {
    float* f = (float*)malloc(sizeof(float)* 13);
    int k, n;
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
            //printf("%f,", cos((PI * (2 * n + 1) * k) / (2 * 26)));
            float tmp = (*(feat + n)) * dct_arg[k * 26 + n];
            s = s + tmp;	//累加求和
        }
        (*(f + k)) = A * s;	//X[k]等于累和结果s乘以系数A
        //printf("\n");
    }
    return f;
}

float* lifter(float* feat, int frame_num, int ceplifter) {
    int nframes = frame_num;
    float lift[13] = { 1.00000000,2.56546330,4.09905815,5.56956530,6.94704914,8.20346832,9.31324577,10.2537889,11.0059519,11.5544224,11.8880358,12.0000000,11.8880358 };
    int i, j;
    for (i = 0; i < nframes; i++) {
        for (j = 0; j < 13; j++)
        {
            *(feat + i * 13 + j) = (*(feat + i * 13 + j)) * lift[j];

        }
    }
    return feat;
}

float* lifter_frame(float* feat, int ceplifter) {
    float lift[13] = { 1.00000000,2.56546330,4.09905815,5.56956530,6.94704914,8.20346832,9.31324577,10.2537889,11.0059519,11.5544224,11.8880358,12.0000000,11.8880358 };
    int j;
    for (j = 0; j < 13; j++)
    {
        *(feat + j) = (*(feat + j)) * lift[j];

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

float* appendEnergy_frame(float* feat, float energy) {
    *(feat) = log(energy);
    return feat;
}

float unit_gaussian(float* data, float* tmp_mu, float* tmp_cov) {
    Matrix x, mu, sigma;
    x = Create_Matrix(1, 13);
    mu = Create_Matrix(1, 13);
    sigma = Create_Matrix(13, 13);
    SetData_Matrix(x, data);
    SetData_Matrix(mu, tmp_mu);
    SetData_Matrix(sigma, tmp_cov);
    Matrix inv_cov = LUInverse_Matrix(sigma);
    Matrix x_mu = AddorSub_Matrix(x, mu, 1);
    float D = mu->column;
    float z = (float)(1 / ((sqrt(pow(2 * PI, D))) * (sqrt(Det_Matrix(sigma)))));
    Matrix tmp = Mult_Matrix(Mult_Matrix(x_mu, inv_cov), Trans_Matrix(x_mu));
    float exponent = (float)(exp(-0.5 * (**(tmp->data))));

    return z * exponent;
}

float calculate_likelihood(int N, int K, float* data, float* mu_k, float* cov_k, float* pi_k) {
    float log_likelihood = 0;
    float temp = 0;
    float tmp_mu[13], tmp_cov[169];
    for (int i = 0; i < N; i++)
    {
        temp = 0;
        for ( int k= 0;  k< K; k++)
        {
            temp +=pi_k[k] * unit_gaussian(data + i * 13, mu_k + k * 13, cov_k + k * 13 * 13);
        }
        log_likelihood += log(temp);
    }
    //log_likelihood += log(temp);
    /*for (int i = 0; i < K; i++) {
        for (int j = 0; j < 13; j++) {
            tmp_mu[j] = mu_k[i * 13 + j];
        }
        for (int j = 0; j < 169; j++) {
            tmp_cov[j] = cov_k[i * 169 + j];
        }


        temp += pi_k[i] * (unit_gaussian(data, tmp_mu, tmp_cov));
    }*/
    return log_likelihood;

}

short* openwav(char* filename, int* fsize) {
    FILE* fp = NULL;
    int count = 0;
    short* sig;

    //打开文件
    fp = fopen("testwav", "rb");

    //计算所需要的空间
    fseek(fp, 0, SEEK_END);
    *fsize = ftell(fp) / 2;
    printf("%d\n", fsize);
    rewind(fp);

    //动态分配内存
    sig = (short*)malloc(sizeof(short) * *(fsize));

    //循环读取到内存（堆）
    while (!feof(fp)) {
        count = fread(sig, sizeof(short), *(fsize), fp);
        int n = feof(fp);
    }
    fclose(fp);
    return sig;
}

void* load_model(char* filepath, float* mu, float* cov, float* weights) {
    FILE* fp = NULL;
    int i;
    //打开文件
    fp = fopen(filepath, "rb");
    //提取数据
    fread(mu, sizeof(float), 208, fp);
    fread(cov, sizeof(float), 2704, fp);
    fread(weights, sizeof(float), 16, fp);
    //关闭文件
    fclose(fp);
}

void* save_model(char* filepath, float* mu, float* cov, float* weights) {
    FILE* fp = NULL;
    int i;
    //创建文件
    fp = fopen(filepath, "wb");
    //存储数据
    fwrite(mu, sizeof(float), 208, fp);
    fwrite(cov, sizeof(float), 2704, fp);
    fwrite(weights, sizeof(float), 16, fp);
    //关闭文件
    fclose(fp);
}

float* mfcc(short* signal, int fsize, int samplerate, float winlen, float winstep, int numcep, float* winfunc, float preemph) {
    //预加重==============================================================
    float* pre_sig;
    pre_sig = preemphasis(signal, preemph, fsize);

    //加窗分帧============================================================
    float* frame;
    int frame_num;
    frame = framesig2(pre_sig, winlen * samplerate, winstep * samplerate, winfunc, fsize, &frame_num);

    //快速傅里叶变化=======================================================
    int i, j, winlenth;
    float* complex_spec = NULL;
    float abs_fft;
    float* energy = NULL;
    float sum = 0.0;

    winlenth = winlen * samplerate;
    complex_spec = (float*)malloc(sizeof(float) * (frame_num * (NFFT / 2 + 1)));
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

    //读取mel滤波==========================================================
    Matrix A, A_t, B, FEAT;
    A = Create_Matrix(26, 257);
    SetData_Matrix(A, mel_odd);
    printf("A(1, 1) = %f, A(1, 2) = %f, A(1, 3) = %f\n", PickInMat(A, 1, 1), PickInMat(A, 1, 2), PickInMat(A, 1, 3));
    A_t = Trans_Matrix(A);
    printf("row = %d, col = %d\n", A_t->row, A_t->column);
    B = Create_Matrix(frame_num, 257);
    SetData_Matrix(B, complex_spec);
    printf("B(1, 1) = %f, B(1, 2) = %f, B(1, 3) = %f\n", PickInMat(B, 1, 1), PickInMat(B, 1, 2), PickInMat(B, 1, 3));
    FEAT = Mult_Matrix(B, A_t);
    printf("feat(1, 1) = %f, feat(1, 2) = %f, feat(95, 26) = %f\n", PickInMat(FEAT, 1, 1), PickInMat(FEAT, 1, 2), PickInMat(FEAT, frame_num, 26));
    printf("FEAT->row = %d\t FEAT->col = %d\n", FEAT->row, FEAT->column);

    Free_Matrix(A);
    Free_Matrix(A_t);
    Free_Matrix(B);

    //log=============================================================
    float* feat;
    feat = (float*)malloc(sizeof(float) * frame_num * 26);
    for (i = 1; i <= frame_num; i++)
    {
        for (j = 1; j <= 26; j++)
        {
            *(feat + (i - 1) * 26 + (j - 1)) = (float)log((double)PickInMat(FEAT, i, j));
        }
    }
    //dct=============================================================
    float* feat_mfcc;
    feat_mfcc = dct(feat, frame_num);


    //用倒谱提升器，得到倒谱的矩阵。这有增加高频DCT系数的作用
    feat_mfcc = lifter(feat_mfcc, frame_num, CEPLIFTER);
    //用帧能量的对数替换第一个倒谱系数
    feat_mfcc = appendEnergy(feat_mfcc, energy, frame_num);
    printf("mfcc[0] = %f\n", feat_mfcc[0]);


    free(pre_sig);
    free(frame);
    free(feat);
    free(complex_spec);
    free(energy);

    return feat_mfcc;
}

float* mfcc_frame(short* signal_frame, float* new_signal, int winlenth, float* winfunc, short last_point) {
    //预加重
    int i = 0;
    //new_signal = (float*)malloc(sizeof(float) * winlenth);
    *new_signal = *(signal_frame)-last_point * 0.97;

    for (i = 1; i < winlenth; i++) {
        (float)(*(new_signal + i)) = (float)(*(signal_frame + i)) - (float)(*(signal_frame + i - 1)) * 0.97;
    }
    //快速傅里叶变化=======================================================
    int j;
    float abs_fft;
    float energy = 0.0;

    if (winlenth < NFFT) {
        for (j = 0; j < winlenth; j++) {
            //temp[j] = *(frame + i * winlenth + j);
            x[j].real = *(new_signal +  j);
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
            x[j].real = *(new_signal + j);
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
        energy += abs_fft * abs_fft / NFFT;
    }

    //读取mel滤波===mel_odd,mel_even==========================================================
    float* FEAT = (float*)malloc(sizeof(float) * 26);
    float temp = 0.0;
    short pre_flag = 1, flag = 0;
    j = 0;
    for ( i = 0; i < 257; i++)
    {
        //确定当前flag状态,flag = 1 表示间隔
        if (mel_even[i] == 0.) {
            flag = 1;
        }
        else {
            flag = 0;
        }
        //对应flag操作
        if (flag == 0) {
            temp += mel_even[i] * *(complex_spec + i);
        }
        else if (flag == 1 && pre_flag == 0) {
            *(FEAT + j) = log(temp);
            //printf("%f", FEAT);
            j += 2;
            temp = 0.0;
        }
        pre_flag = flag;
    }
    pre_flag = 1;
    flag = 0;
    j = 1;
    for (i = 0; i < 257; i++)
    {
        //确定当前flag状态,flag = 1 表示间隔
        if (mel_odd[i] == 0.) {
            flag = 1;
        }
        else {
            flag = 0;
        }
        //对应flag操作
        if (flag == 0) {
            temp += mel_odd[i] * *(complex_spec + i);
        }
        else if (flag == 1 && pre_flag == 0) {
            *(FEAT + j) = log(temp);
            //printf("%f", FEAT);
            j += 2;
            temp = 0.0;
        }
        pre_flag = flag;
    }

    //dct=============================================================
    float* feat_mfcc;
    feat_mfcc = dct_frame(FEAT);
    //用倒谱提升器，得到倒谱的矩阵。这有增加高频DCT系数的作用
    feat_mfcc = lifter_frame(feat_mfcc, CEPLIFTER);
    //用帧能量的对数替换第一个倒谱系数
    *(feat_mfcc) = log(energy);
    //printf("mfcc[0] = %f\n", feat_mfcc[0]);
    return feat_mfcc;
}

void map_adaptation_frame(float* mfcc, float* mu, float* pi , float* cov, short threshold, short max_iterations, short K, short rf) {
    float old_likelihood = 9999;
    float new_likelihood = 0;
    short iterations = 0;
    
    float* num = (float*)malloc(sizeof(float) * K);
    float* new_mu = (float*)malloc(sizeof(float) * K * 13);
    
    int i, j;
    float sum = 0.0;
    float adaption = 0.;
    float temp = 0.;

    while (fabs((double)old_likelihood - (double)new_likelihood) > threshold && iterations < max_iterations)
    {
        iterations += 1;
        old_likelihood = new_likelihood;
        printf("epoch:%d\n", iterations);
        //float* num = (float*)malloc(sizeof(float) * K);
        for ( i = 0; i < K; i++)
        {
            *(num + i) = *(pi + i) * unit_gaussian(mfcc, mu + i*13, cov + i*13*13);
            sum += *(num + i);
        }
        //float* new_mu = (float*)malloc(sizeof(float) * K * 13);

        for ( i = 0; i < K; i++)
        {
            *(num + i) /= sum;
            for (j = 0; j < 13; j++)
            {
                *(new_mu + i * 13 + j) += *(num + i) + *(mfcc + j);
            }
        }
        //float* adaptation = (float*)malloc(sizeof(float) * K);
        for (i = 0; i < K; i++)
        {
            //*(adaptation + i) = *(num + i) / (*(num + i) + rf);
            adaption = *(num + i) / (*(num + i) + rf);
            for ( j = 0; j < 13; j++)
            {
                //temp = (1 - *(adaptation + i)) * (*(mu + i * 13 * j));
                temp = (1 - adaption) * (*(mu + i * 13 * j));
                //*(mu + i * 13 * j) = temp + (*(adaptation + i) * (*(new_mu + i * 13 + j)));
                *(mu + i * 13 * j) = temp + adaption * (*(new_mu + i * 13 + j));
            }
        }
        //new_likelihood = calculate_likelihood(K, mfcc, mu, cov, pi);
    }

    free(new_mu);
    free(num);
    //save_model("map", mu, cov, pi);
}

void map_adaptation(char* mfcc_file, int frames, float* mu, float* pi, float* cov, short threshold, short max_iterations, short K, short rf) {
    //读取mfcc
    float* mfcc = (float*)malloc(sizeof(float) * 13 * frames);
    FILE* fp = NULL;
    fp = fopen(mfcc_file, "rb");
    fread(mfcc, sizeof(float), 13 * frames, fp);

    //
    float old_likelihood = 9999;
    float new_likelihood = 0;
    short iterations = 0;

    float* num = (float*)malloc(sizeof(float) * K);
    float* new_mu = (float*)malloc(sizeof(float) * K * 13);
    float* z_n_k = (float*)malloc(sizeof(float) * K * frames);
    float* n_k = (float*)malloc(sizeof(float) * frames);

    int i, j;
    float* sum = (float*)malloc(sizeof(float) * frames);
    float* temp = (float*)malloc(sizeof(float) * 13);
    float* adaptation_coefficient = (float*)malloc(sizeof(float) * K);
    float adaption = 0.;
    float temp2 = 0.;

    while (fabs((double)old_likelihood - (double)new_likelihood) > threshold && iterations < max_iterations)
    {
        iterations += 1;
        old_likelihood = new_likelihood;
        printf("epoch:%d\n", iterations);

        for (j = 0; j < frames; j++)
        {
            *(sum + j) = 0.0;
            for (i = 0; i < K; i++)
            {
                *(num + i) = *(pi + i) * unit_gaussian((mfcc + j * 13), mu + i * 13, cov + i * 13 * 13);
                *(sum + j) += *(num + i);
            }
            for (i = 0; i < K; i++)
            {
                *(num + i) /= *(sum + j);
                *(z_n_k + j * 16 + i) = *(num + i);
            }
        }

        memset(n_k, 0, sizeof(n_k)*K);
        for (j = 0; j < K; j++)
        {
            for (i = 0; i < frames; i++)
            {
                *(n_k + j) += *(z_n_k + j + i * 16);
            }
            *(n_k + j) += 1e-10;
        }



        for (i = 0; i < K; i++)
        {
            for (j = 0; j < 13; j++)
            {
                *(temp + j) = 0.0;
            }

            for (j = 0; j < frames; j++)
            {
                for (int k = 0; k < 13; k++)
                {
                    *(temp + k) += *(z_n_k + j * 16 + i) * (*(mfcc + j * 13 + k));
                }
            }

            for (j = 0; j < 13; j++)
            {
                *(new_mu + i * 13 + j) = *(temp + j) / (*(n_k + i));
            }
        }


        for (i = 0; i < K; i++)
        {
            *(adaptation_coefficient + i) = (*(n_k + i)) / ((*(n_k + i)) + rf);
        }




        for (i = 0; i < K; i++)
        {
            for (j = 0; j < 13; j++)
            {

                temp2 = (1 - (*(adaptation_coefficient + i))) * (*(mu + i * 13 + j));

                *(mu + i * 13 + j) = temp2 + (*(adaptation_coefficient + i)) * (*(new_mu + i * 13 + j));
            }
        }
        new_likelihood = calculate_likelihood(frames, K, mfcc, mu, cov, pi);
        printf("calculate likelihood:%f", new_likelihood);
    }

    free(new_mu);
    free(num);
    //save_model("map", mu, cov, pi);
}

float test_model(short start, short end, char* ubm_filepath, char* map_filepath, char* mfcc_filepath) {
    FILE* fp = NULL;
    float score;
    short n;
    n = end - start;
    //加载ubm数据
    float* ubm_mu = (float*)malloc(sizeof(float) * 208);
    float* ubm_cov = (float*)malloc(sizeof(float) * 2704);
    float* ubm_we = (float*)malloc(sizeof(float) * 16);
    load_model(ubm_filepath, ubm_mu, ubm_cov, ubm_we);
    //加载map数据
    float* map_mu = (float*)malloc(sizeof(float) * 208);
    float* map_cov = (float*)malloc(sizeof(float) * 2704);
    float* map_we = (float*)malloc(sizeof(float) * 16);
    load_model(map_filepath, map_mu, map_cov, map_we);
    //加载mfcc数据
    fp = fopen(mfcc_filepath, "wb");
    float* data = (float*)malloc(sizeof(float) * n * 13);
    fread(data, sizeof(float), n * 13, fp);
    //计算
    score = calculate_likelihood(n, 16, data, map_mu, map_cov, map_we) - calculate_likelihood(n, 16, data, ubm_mu, ubm_cov, ubm_we);
    //关闭文件并返回
    fclose(fp);

    free(data);
    free(ubm_cov);
    free(ubm_mu);
    free(ubm_we);
    free(map_cov);
    free(map_mu);
    free(map_we);
    return score;
}

//void main() {
//	printf("program start!\n");
//    //读音频=============================================================
//    char* filename = "testwav";
//    short* sig;
//    int fsize;
//    sig = openwav(filename, &fsize);
//	//printf("%d\n", sig[0]);
//
//	//提取mfcc===========================================================
//	int samplerate = 16000;
//	float winlen = 0.025;
//    float winstep = 0.01;
//    float preemph = 0.97;
//	float* winfunc = NULL;
//    int numcep = 13;
//
//    float* feat;
//    feat = mfcc(sig, fsize, samplerate, winlen, winstep, numcep, winfunc, preemph);
//	
//    printf("program end!\n");
//	return;
//}
//

//void main() {
//    printf("program start!\n");
//    //读音频=============================================================
//    char* filename = "testwav";
//    short* sig;
//    int fsize;
//    sig = openwav(filename, &fsize);
//    //printf("%d\n", sig[0]);
//    
//    //单帧处理 fsize=15360 
//    int frames = 0;
//    short flag = 0;
//    fsize -= 240;
//    frames = fsize / 160;
//    flag = fsize % 160;
//
//    if (flag) {
//        frames += 1;
//    }
//    int start = 0, i, j;
//    int winlength = 400;
//    short last_point = 0;
//    short* sig_frame;
//    float* new_signal;
//    float* feat = NULL;
//    sig_frame = (short*)malloc(sizeof(short) * 400);
//    new_signal = (float*)malloc(sizeof(float) * 400);
//
//    FILE* fp = NULL;
//    fp = fopen("mfcc", "wb");
//    float* mu = (float*)malloc(sizeof(float) * 208);
//    float* cov = (float*)malloc(sizeof(float) * 2704);
//    float* pi = (float*)malloc(sizeof(float) * 16);
//    load_model("ubm", mu, cov, pi);
//    for (i = 0; i < frames; i++)
//    {
//        if (i > 0) {
//            last_point = *(sig + start - 1);
//        }
//        for (j = 0; j < winlength; j++)
//        {
//            *(sig_frame + j) = *(sig + start + j);
//        }
//        if (i == frames - 1) {
//            for (j = 0; j < flag; j++)
//            {
//                *(sig_frame + winlength - flag + j) = (short)0;
//            }
//        }
//        start += 160;
//        feat = mfcc_frame(sig_frame, new_signal, winlength, NULL, last_point);
//        fwrite(feat, sizeof(float), 13, fp);
//        //map_adaptation_frame(feat, mu, pi, cov, 1, 10, 16, 1);
//    }
//    save_model("map", mu, cov, pi);
//
//    printf("program end!\n");
//    return;
//}

void main() {
    printf("program start!\n");
    //读音频=============================================================
    char* filename = "testwav";
    short* sig;
    int fsize;
    sig = openwav(filename, &fsize);
    //printf("%d\n", sig[0]);

    //单帧处理 fsize=15360 
    int frames = 0;
    short flag = 0;
    fsize -= 240;
    frames = fsize / 160;
    flag = fsize % 160;

    if (flag) {
        frames += 1;
    }
    int start = 0, i, j;
    int winlength = 400;
    short last_point = 0;
    short* sig_frame;
    float* new_signal;
    float* feat = NULL;
    sig_frame = (short*)malloc(sizeof(short) * 400);
    new_signal = (float*)malloc(sizeof(float) * 400);

    FILE* fp = NULL;
    fp = fopen("mfcc", "wb");
    float* mu = (float*)malloc(sizeof(float) * 16 * 13);
    float* cov = (float*)malloc(sizeof(float) * 16 * 13 * 13);
    float* pi = (float*)malloc(sizeof(float) * 16);
    load_model("ubm", mu, cov, pi);
    for (i = 0; i < frames; i++)
    {
        if (i > 0) {
            last_point = *(sig + start - 1);
        }
        for (j = 0; j < winlength; j++)
        {
            *(sig_frame + j) = *(sig + start + j);
        }
        if (i == frames - 1) {
            for (j = 0; j < flag; j++)
            {
                *(sig_frame + winlength - flag + j) = (short)0;
            }
        }
        start += 160;
        feat = mfcc_frame(sig_frame, new_signal, winlength, NULL, last_point);
        fwrite(feat, sizeof(float), 13, fp);
    }
    fclose(fp);
    map_adaptation("mfcc", frames, mu, pi, cov, 1, 100, 16, 16);
    //save_model("map", mu, cov, pi);

    printf("program end!\n");
    return;
}

//[1,2,3,4,5,6,7,8] --> [1, 2-0.97*1, 3-0.97*2, 4-3*0.97, 5-4*0.97, 6 - 5 * 0.97, 7 - 6 * 0.97, 8 - 7 * 0.97]
//[1, 2, 3, 4] --> 1 - 0.97 * 0, 2 - 1 * 0.97, 3 - 3 * 0.97, 4 - 3 * 0.97
//[2, 3, 4, 5] --> 2 - 0.97 * 1, 3 - 2 * 0.97, 4 - 3 * 0.97, 5 - 4 * 0.97
//[3, 4, 5, 6] --> 3 - 0.97 * 2, 4 - 3 * 0.97, 5 - 4 * 0.97, 6 - 5 * 0.97
//[4, 5, 6, 7]
//
//[2-0.97*1,3]


