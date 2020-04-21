/*
 * MyMatrix.c
 *
 *  Created on: Jul 13, 2019
 *      Author: xuuyann
 */

#include "MyMatrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void Init_Matrix(Matrix mat)
{
	int i, j;
	for (i = 0; i < mat->row; i++) {
		for (j = 0; j < mat->column; j++) {
			mat->data[i][j] = 0;
		}
	}
}

// 释放申请的矩阵空间
void Free_Matrix(Matrix mat)
{
	for (int i = 0; i < mat->row; i++)
		free(mat->data[i]); // 释放行指针
	free(mat->data); // 释放头指针
	free(mat); // 释放结构体指针
}

Matrix Create_Matrix(int row, int col)
{
	Matrix mat;
	mat = (Matrix)malloc(sizeof(struct MNode)); //　分配结构体指针
	if (row <= 0 || col <= 0) {
		printf("ERROR, in creat_Matrix the row or col <= 0\n");
		exit(1);
	}
	if (row > 0 && col > 0) {
		mat->row = row;
		mat->column = col;
		mat->data = (float**)malloc(row * sizeof(float*));// 分配头指针
		if (mat->data == NULL) {
			printf("ERROR, in creat_Matrix the mat->data == NULL\n");
			exit(1);
		}
		int i;
		for (i = 0; i < row; i++) {
			*(mat->data + i) = (float*)malloc(col * sizeof(float)); //　分配每行的指针
			if (mat->data[i] == NULL) {
				printf("ERROR, in create_Matrix the mat->data[i] == NULL\n");
				exit(1);
			}
		}
		Init_Matrix(mat);
	}
	return mat;
}

void Show_Matrix(Matrix mat, char* s)
{
	int i, j;
	printf("%s\n", s);
	for (i = 0; i < mat->row; i++) {
		for (j = 0; j < mat->column; j++)
			printf("%.6f\t", mat->data[i][j]);
		printf("\n");
	}
	printf("\n");
}

void SetData_Matrix(Matrix mat, float data[])
{
	int i, j;
	for (i = 0; i < mat->row; i++) {
		for (j = 0; j < mat->column; j++) {
			mat->data[i][j] = data[i * mat->column + j];
		}
	}
}

//flag = 0, add; flag = 1, sub
Matrix AddorSub_Matrix(Matrix mat_1, Matrix mat_2, int flag)
{
	Matrix rst_mat;
	if (mat_1->column != mat_2->column) {
		printf("ERROR in AddorSub, column !=\n");
		exit(1);
	}
	if (mat_1->row != mat_2->row) {
		printf("ERROR in AddorSub, row !=\n");
		exit(1);
	}
	int i, j;
	rst_mat = Create_Matrix(mat_1->row, mat_1->column);
	for (i = 0; i < mat_1->row; i++) {
		for (j = 0; j < mat_1->column; j++)
			rst_mat->data[i][j] = mat_1->data[i][j] + pow(-1, flag) * mat_2->data[i][j];
	}
	return rst_mat;
}

//转置
Matrix Trans_Matrix(Matrix mat)
{
	Matrix mat_;
	int i, j;
	mat_ = Create_Matrix(mat->column, mat->row);
	for (i = 0; i < mat->row; i++) {
		for (j = 0; j < mat->column; j++)
			mat_->data[j][i] = mat->data[i][j];
	}
	return mat_;
}

Matrix Mult_Matrix(Matrix mat_1, Matrix mat_2)
{
	Matrix rst_mat;
	int i, j, m;
	if (mat_1->column != mat_2->row) {
		printf("ERROR in Mult_Matrix, column != row\n");
		exit(1);
	}
	else {
		rst_mat = Create_Matrix(mat_1->row, mat_2->column);
		for (i = 0; i < mat_1->row; i++) {
			for (j = 0; j < mat_2->column; j++) {
				for (m = 0; m < mat_1->column; m++)
					rst_mat->data[i][j] += mat_1->data[i][m] * mat_2->data[m][j];
			}
		}
	}
	return rst_mat;
}

Matrix eye(int n)
{
	Matrix E;
	int i, j;
	if (n <= 0) {
		printf("ERROR in eye\n");
		exit(1);
	}
	E = Create_Matrix(n, n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i == j)
				E->data[i][j] = 1;
			else
				E->data[i][j] = 0;
		}
	}
	return E;
}

float PickInMat(Matrix mat, int r, int c)
{
	float rst;
	rst = mat->data[r - 1][c - 1];
	return rst;
}

//　顺序高斯消去法
Matrix Gauss_shunxu(Matrix A, Matrix b)
{
	int n;
	n = b->row;
	Matrix x;
	x = Create_Matrix(b->row, b->column);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			if (A->data[i][j] == 0) {
				printf("can't use the Gauss_shunxu\n");
				exit(1);
			}
	}
	// 消元
	double* L;
	L = (double*)malloc(sizeof(double) * n);
	for (int k = 0; k < n - 1; k++) {
		for (int i = k + 1; i < n; i++) {
			L[i] = A->data[i][k] / A->data[k][k];
			for (int j = k + 1; j < n; j++)
				A->data[i][j] = A->data[i][j] - L[i] * A->data[k][j];
			b->data[i][0] = b->data[i][0] - L[i] * b->data[k][0];
		}
	}
	for (int i = 0; i < n; i++) {
		if (A->data[i][i] == 0) {
			printf("can't use the Gauss_shunxu\n");
			exit(1);
		}
	}
	// 回代
	x->data[n - 1][0] = b->data[n - 1][0] / A->data[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--) {
		double sum_a = 0;
		for (int j = i + 1; j < n; j++)
			sum_a += A->data[i][j] * x->data[j][0];
		x->data[i][0] = (b->data[i][0] - sum_a) / A->data[i][i];
	}
	free(L);
	return x;
}

//　列主元高斯消去法
Matrix Gauss_lie(Matrix A, Matrix b)
{
	int n;
	n = b->row;
	Matrix x;
	x = Create_Matrix(b->row, b->column);
	double* L;
	L = (double*)malloc(sizeof(double) * n);
	for (int k = 0; k < n - 1; k++) {
		int cnt = k;
		double a_max = fabs(A->data[k][k]);
		for (int i = k; i < n; i++) {
			if (fabs(A->data[i][k]) > a_max) {
				a_max = fabs(A->data[i][k]);
				cnt = i; // 确定下标i
			}
		}
		if (A->data[cnt][k] == 0) {
			printf("Gauss_lie: no unique solution\n");
			exit(1);
		}
		// 换行
		if (cnt != k) {
			double t = 0, s = 0;
			for (int j = k; j < n; j++) {
				t = A->data[k][j];
				A->data[k][j] = A->data[cnt][j];
				A->data[cnt][j] = t;
				s = b->data[cnt][0];
				b->data[cnt][0] = b->data[k][0];
				b->data[k][0] = s;
			}
		}
		// step 5
		for (int i = k + 1; i < n; i++) {
			L[i] = A->data[i][k] / A->data[k][k];
			for (int j = k + 1; j < n; j++)
				A->data[i][j] = A->data[i][j] - L[i] * A->data[k][j];
			b->data[i][0] = b->data[i][0] - L[i] * b->data[k][0];
		}
	}
	for (int i = 0; i < n; i++) {
		if (A->data[i][i] == 0.0) {
			printf("Gauss_lie: no unique solution\n");
			exit(1);
		}
	}
	// 回代
	x->data[n - 1][0] = b->data[n - 1][0] / A->data[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--) {
		double sum_a = 0;
		for (int j = i + 1; j < n; j++)
			sum_a += A->data[i][j] * x->data[j][0];
		x->data[i][0] = (b->data[i][0] - sum_a) / A->data[i][i];
	}
	free(L);
	return x;
}

//// 3×3向量叉乘
//Matrix Cross(Matrix a, Matrix b)
//{
//	Matrix C;
//	C = Create_Matrix(a->row, a->column);
//	double ax, ay, az;
//	double bx, by, bz;
//	double c[a->row];
//	ax = PickInMat(a, 1, a->column);
//	ay = PickInMat(a, 2, a->column);
//	az = PickInMat(a, 3, a->column);
//	bx = PickInMat(b, 1, b->column);
//	by = PickInMat(b, 2, b->column);
//	bz = PickInMat(b, 3, b->column);
//	c[0] = ay * bz - az * by;
//	c[1] = az * bx - ax * bz;
//	c[2] = ax * by - ay * bx;
//	SetData_Matrix(C, c);
//	return C;
//}

// 矩阵求逆，利用矩阵LU分解求逆（直接三角分解）
Matrix LUInverse_Matrix(Matrix mat)
{
	Matrix inv;
	if (mat->column != mat->row) {
		printf("ERROR in inverse, the row != the column\n");
		exit(1);
	}
	if (Det_Matrix(mat) == 0) {
		printf("The Matrix is not invertible\n");
		exit(1);
	}
	// Ｌ矩阵和Ｕ矩阵
	Matrix L, U;
	int n = mat->column;
	inv = Create_Matrix(n, n);
	L = Create_Matrix(n, n);
	U = Create_Matrix(n, n);
	// 计算Ｕ的第一行元素
	for (int j = 0; j < n; j++)
		U->data[0][j] = mat->data[0][j];
	//　计算Ｌ的第一列元素
	for (int i = 0; i < n; i++) {
		L->data[i][0] = mat->data[i][0] / U->data[0][0];
		L->data[i][i] = 1.0;
	}
	double sum_u = 0, sum_l = 0;
	//　分别计算Ｕ和Ｌ的第２到ｎ行、列元素
	for (int k = 1; k < n; k++) {
		// 求U，U矩阵从第２行迭代到第n行，且U矩阵先于L矩阵一个节拍
		for (int j = k; j < n; j++) {
			sum_u = 0;
			for (int m = 0; m <= k - 1; m++)
				sum_u += L->data[k][m] * U->data[m][j];
			U->data[k][j] = mat->data[k][j] - sum_u;
		}
		//　求Ｌ，L矩阵从第２列迭代到第n列
		for (int i = k + 1; i < n; i++) {
			sum_l = 0;
			for (int m = 0; m <= k - 1; m++)
				sum_l += L->data[i][m] * U->data[m][k];
			L->data[i][k] = (mat->data[i][k] - sum_l) / U->data[k][k];
		}
	}
	Show_Matrix(U, "U");
	Show_Matrix(L, "L");
	// 分别求下三角矩阵Ｌ和上三角矩阵Ｕ的逆矩阵
	Matrix L_, U_;
	L_ = Create_Matrix(n, n);
	U_ = Create_Matrix(n, n);
	// 求矩阵Ｕ的逆
	double sum_u_;
	for (int i = 0; i < n; i++) {
		U_->data[i][i] = 1.0 / U->data[i][i];// 对角线元素的值直接取倒数
		for (int j = i - 1; j >= 0; j--) {
			sum_u_ = 0;
			for (int k = j + 1; k <= i; k++)
				sum_u_ += U->data[j][k] * U_->data[k][i];
			U_->data[j][i] = -sum_u_ / U->data[j][j]; // 迭代计算，按列倒序依次得到每一个值
		}
	}
	// 求Ｌ的逆
	for (int i = 0; i < n; i++) {
		L_->data[i][i] = 1; // 对角线元素的值直接取倒数，这里为１
		for (int k = i + 1; k < n; k++) {
			for (int j = i; j <= k - 1; j++)
				L_->data[k][i] -= L->data[k][j] * L_->data[j][i]; // 迭代计算，按列顺序依次得到每一个值
		}
	}
	//Show_Matrix(U_, "U_");
	//Show_Matrix(L_, "L_");
	// 已经得到Ｌ和Ｕ的逆矩阵
	inv = Mult_Matrix(U_, L_);

	Free_Matrix(L_);
	Free_Matrix(U_);
	Free_Matrix(L);
	Free_Matrix(U);
	return inv;
}

// 矩阵第n行与第m行互换
void Swap_row(Matrix mat, int n, int m)
{
	double temp;
	for (int i = 0; i < mat->column; i++) {
		temp = mat->data[n - 1][i];
		mat->data[n - 1][i] = mat->data[m - 1][i];
		mat->data[m - 1][i] = temp;
	}
}

// 矩阵求逆，利用初等行变换求逆
Matrix EleTransInv_Matrix(Matrix mat)
{
	if (mat->column != mat->row) {
		printf("ERROR in inverse, the row != the column\n");
		exit(1);
	}
	if (Det_Matrix(mat) == 0) {
		printf("The Matrix is not invertible\n");
		exit(1);
	}
	int n = mat->row;
	// 为防止改变原矩阵，此处进行复制处理
	Matrix mat_ = Copy_Matrix(mat);
	//　创建单位矩阵
	Matrix inv_eye = eye(n);
	double e, a_max;
	int i, j, k, t, cnt;
	for (k = 0; k < n - 1; k++) {
		a_max = fabs(mat_->data[k][k]);
		cnt = k;
		// 选主元
		for (i = k; i < n; i++) {
			if (fabs(mat_->data[i][k]) > a_max) {
				a_max = fabs(mat_->data[i][k]);
				cnt = i;
			}
		}
		//　换行，原矩阵换行的同时，单位矩阵也换行
		double temp, temp_e;
		if (cnt != k) {
			for (j = 0; j < n; j++) {
				temp = mat_->data[k][j];
				mat_->data[k][j] = mat_->data[cnt][j];
				mat_->data[cnt][j] = temp;
				temp_e = inv_eye->data[k][j];
				inv_eye->data[k][j] = inv_eye->data[cnt][j];
				inv_eye->data[cnt][j] = temp_e;
			}
		}
		// 消元
		for (i = k + 1; i < n; i++) {
			e = mat_->data[i][k] / mat_->data[k][k];
			for (j = 0; j < n; j++) {
				mat_->data[i][j] = mat_->data[i][j] - e * mat_->data[k][j];
				inv_eye->data[i][j] = inv_eye->data[i][j] - e * inv_eye->data[k][j];
			}
		}
	}
	// 将矩阵每行的行首元素化为１
	for (i = 0; i < n; i++) {
		e = 1.0 / mat_->data[i][i];
		for (j = 0; j < n; j++) {
			mat_->data[i][j] = mat_->data[i][j] * e;
			inv_eye->data[i][j] = inv_eye->data[i][j] * e;
		}
	}
	// 从最后一排开始消元，把增广矩阵左边的矩阵化为单位矩阵
	for (i = n - 1; i > 0; i--) {
		for (j = i - 1; j >= 0; j--) {
			e = mat_->data[j][i] / mat_->data[i][i];
			for (t = 0; t < n; t++) {
				mat_->data[j][t] = mat_->data[j][t] - e * mat_->data[i][t];
				inv_eye->data[j][t] = inv_eye->data[j][t] - e * inv_eye->data[i][t];
			}
		}
	}
	Free_Matrix(mat_);
	return inv_eye;
}

/*
 * 通过行变换（列主元消去法）将矩阵变换成上三角矩阵
 * 注意行变换会改变行列式的符号
 */
double Det_Matrix(Matrix mat)
{
	double det = 1.0;
	if (mat->row != mat->column) {
		printf("ERROR in Det_Matrix, the row != the column\n");
		exit(1);
	}
	if (mat->row == 1 && mat->column == 1)
		return mat->data[0][0];
	// 为防止改变原矩阵，此处进行复制处理
	Matrix mat_ = Copy_Matrix(mat);
	int n = mat_->row, s = 0, cnt;
	double L, a_max;
	for (int k = 0; k < n - 1; k++) {
		cnt = k;
		a_max = fabs(mat_->data[k][k]);
		for (int i = k; i < n; i++) {
			if (fabs(mat_->data[i][k]) > a_max) {
				a_max = fabs(mat_->data[i][k]);
				cnt = i; // 确定下标i
			}
		}
		//　换行
		double temp;
		if (cnt != k) {
			s++; // 换行次数
			for (int j = 0; j < n; j++) {
				temp = mat_->data[cnt][j];
				mat_->data[cnt][j] = mat_->data[k][j];
				mat_->data[k][j] = temp;
			}
		}
		// 消元计算
		for (int i = k + 1; i < n; i++) {
			L = mat_->data[i][k] / mat_->data[k][k];
			for (int j = k + 1; j < n; j++)
				mat_->data[i][j] = mat_->data[i][j] - L * mat_->data[k][j];
		}
	}
	if (s % 2 == 0)
		s = 1;
	else
		s = -1;
	for (int i = 0; i < n; i++)
		det *= mat_->data[i][i];
	det *= s;
	Free_Matrix(mat_); //释放掉复制矩阵的内存
	return det;
}


//　对一个矩阵进行复制
Matrix Copy_Matrix(Matrix mat)
{
	Matrix copy_mat = Create_Matrix(mat->row, mat->column);
	for (int i = 0; i < mat->row; i++) {
		for (int j = 0; j < mat->column; j++)
			copy_mat->data[i][j] = mat->data[i][j];
	}
	return copy_mat;
}
