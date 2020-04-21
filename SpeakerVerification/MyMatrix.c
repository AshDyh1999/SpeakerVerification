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

// �ͷ�����ľ���ռ�
void Free_Matrix(Matrix mat)
{
	for (int i = 0; i < mat->row; i++)
		free(mat->data[i]); // �ͷ���ָ��
	free(mat->data); // �ͷ�ͷָ��
	free(mat); // �ͷŽṹ��ָ��
}

Matrix Create_Matrix(int row, int col)
{
	Matrix mat;
	mat = (Matrix)malloc(sizeof(struct MNode)); //������ṹ��ָ��
	if (row <= 0 || col <= 0) {
		printf("ERROR, in creat_Matrix the row or col <= 0\n");
		exit(1);
	}
	if (row > 0 && col > 0) {
		mat->row = row;
		mat->column = col;
		mat->data = (float**)malloc(row * sizeof(float*));// ����ͷָ��
		if (mat->data == NULL) {
			printf("ERROR, in creat_Matrix the mat->data == NULL\n");
			exit(1);
		}
		int i;
		for (i = 0; i < row; i++) {
			*(mat->data + i) = (float*)malloc(col * sizeof(float)); //������ÿ�е�ָ��
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

//ת��
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

//��˳���˹��ȥ��
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
	// ��Ԫ
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
	// �ش�
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

//������Ԫ��˹��ȥ��
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
				cnt = i; // ȷ���±�i
			}
		}
		if (A->data[cnt][k] == 0) {
			printf("Gauss_lie: no unique solution\n");
			exit(1);
		}
		// ����
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
	// �ش�
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

//// 3��3�������
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

// �������棬���þ���LU�ֽ����棨ֱ�����Ƿֽ⣩
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
	// �̾���ͣվ���
	Matrix L, U;
	int n = mat->column;
	inv = Create_Matrix(n, n);
	L = Create_Matrix(n, n);
	U = Create_Matrix(n, n);
	// ����յĵ�һ��Ԫ��
	for (int j = 0; j < n; j++)
		U->data[0][j] = mat->data[0][j];
	//������̵ĵ�һ��Ԫ��
	for (int i = 0; i < n; i++) {
		L->data[i][0] = mat->data[i][0] / U->data[0][0];
		L->data[i][i] = 1.0;
	}
	double sum_u = 0, sum_l = 0;
	//���ֱ����պ̵ͣĵڣ������С���Ԫ��
	for (int k = 1; k < n; k++) {
		// ��U��U����ӵڣ��е�������n�У���U��������L����һ������
		for (int j = k; j < n; j++) {
			sum_u = 0;
			for (int m = 0; m <= k - 1; m++)
				sum_u += L->data[k][m] * U->data[m][j];
			U->data[k][j] = mat->data[k][j] - sum_u;
		}
		//����̣�L����ӵڣ��е�������n��
		for (int i = k + 1; i < n; i++) {
			sum_l = 0;
			for (int m = 0; m <= k - 1; m++)
				sum_l += L->data[i][m] * U->data[m][k];
			L->data[i][k] = (mat->data[i][k] - sum_l) / U->data[k][k];
		}
	}
	Show_Matrix(U, "U");
	Show_Matrix(L, "L");
	// �ֱ��������Ǿ���̺������Ǿ���յ������
	Matrix L_, U_;
	L_ = Create_Matrix(n, n);
	U_ = Create_Matrix(n, n);
	// �����յ���
	double sum_u_;
	for (int i = 0; i < n; i++) {
		U_->data[i][i] = 1.0 / U->data[i][i];// �Խ���Ԫ�ص�ֱֵ��ȡ����
		for (int j = i - 1; j >= 0; j--) {
			sum_u_ = 0;
			for (int k = j + 1; k <= i; k++)
				sum_u_ += U->data[j][k] * U_->data[k][i];
			U_->data[j][i] = -sum_u_ / U->data[j][j]; // �������㣬���е������εõ�ÿһ��ֵ
		}
	}
	// ��̵���
	for (int i = 0; i < n; i++) {
		L_->data[i][i] = 1; // �Խ���Ԫ�ص�ֱֵ��ȡ����������Ϊ��
		for (int k = i + 1; k < n; k++) {
			for (int j = i; j <= k - 1; j++)
				L_->data[k][i] -= L->data[k][j] * L_->data[j][i]; // �������㣬����˳�����εõ�ÿһ��ֵ
		}
	}
	//Show_Matrix(U_, "U_");
	//Show_Matrix(L_, "L_");
	// �Ѿ��õ��̺ͣյ������
	inv = Mult_Matrix(U_, L_);

	Free_Matrix(L_);
	Free_Matrix(U_);
	Free_Matrix(L);
	Free_Matrix(U);
	return inv;
}

// �����n�����m�л���
void Swap_row(Matrix mat, int n, int m)
{
	double temp;
	for (int i = 0; i < mat->column; i++) {
		temp = mat->data[n - 1][i];
		mat->data[n - 1][i] = mat->data[m - 1][i];
		mat->data[m - 1][i] = temp;
	}
}

// �������棬���ó����б任����
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
	// Ϊ��ֹ�ı�ԭ���󣬴˴����и��ƴ���
	Matrix mat_ = Copy_Matrix(mat);
	//��������λ����
	Matrix inv_eye = eye(n);
	double e, a_max;
	int i, j, k, t, cnt;
	for (k = 0; k < n - 1; k++) {
		a_max = fabs(mat_->data[k][k]);
		cnt = k;
		// ѡ��Ԫ
		for (i = k; i < n; i++) {
			if (fabs(mat_->data[i][k]) > a_max) {
				a_max = fabs(mat_->data[i][k]);
				cnt = i;
			}
		}
		//�����У�ԭ�����е�ͬʱ����λ����Ҳ����
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
		// ��Ԫ
		for (i = k + 1; i < n; i++) {
			e = mat_->data[i][k] / mat_->data[k][k];
			for (j = 0; j < n; j++) {
				mat_->data[i][j] = mat_->data[i][j] - e * mat_->data[k][j];
				inv_eye->data[i][j] = inv_eye->data[i][j] - e * inv_eye->data[k][j];
			}
		}
	}
	// ������ÿ�е�����Ԫ�ػ�Ϊ��
	for (i = 0; i < n; i++) {
		e = 1.0 / mat_->data[i][i];
		for (j = 0; j < n; j++) {
			mat_->data[i][j] = mat_->data[i][j] * e;
			inv_eye->data[i][j] = inv_eye->data[i][j] * e;
		}
	}
	// �����һ�ſ�ʼ��Ԫ�������������ߵľ���Ϊ��λ����
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
 * ͨ���б任������Ԫ��ȥ����������任�������Ǿ���
 * ע���б任��ı�����ʽ�ķ���
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
	// Ϊ��ֹ�ı�ԭ���󣬴˴����и��ƴ���
	Matrix mat_ = Copy_Matrix(mat);
	int n = mat_->row, s = 0, cnt;
	double L, a_max;
	for (int k = 0; k < n - 1; k++) {
		cnt = k;
		a_max = fabs(mat_->data[k][k]);
		for (int i = k; i < n; i++) {
			if (fabs(mat_->data[i][k]) > a_max) {
				a_max = fabs(mat_->data[i][k]);
				cnt = i; // ȷ���±�i
			}
		}
		//������
		double temp;
		if (cnt != k) {
			s++; // ���д���
			for (int j = 0; j < n; j++) {
				temp = mat_->data[cnt][j];
				mat_->data[cnt][j] = mat_->data[k][j];
				mat_->data[k][j] = temp;
			}
		}
		// ��Ԫ����
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
	Free_Matrix(mat_); //�ͷŵ����ƾ�����ڴ�
	return det;
}


//����һ��������и���
Matrix Copy_Matrix(Matrix mat)
{
	Matrix copy_mat = Create_Matrix(mat->row, mat->column);
	for (int i = 0; i < mat->row; i++) {
		for (int j = 0; j < mat->column; j++)
			copy_mat->data[i][j] = mat->data[i][j];
	}
	return copy_mat;
}
