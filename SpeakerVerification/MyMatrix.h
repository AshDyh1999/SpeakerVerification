/*
 * MyMatrix.h
 *
 *  Created on: Jul 13, 2019
 *      Author: xuuyann
 */

#ifndef HEADER_MYMATRIX_H_
#define HEADER_MYMATRIX_H_
 /*
  * �����¼��
  *��* �����ڴ����õ�������Ϊÿ�δ������󶼻�����һ�����ڴ棬��û��free��
  *��* LU�ֽ����治̫��
  */
typedef struct MNode* PtrToMNode;
struct MNode
{
    int row;
    int column;
    float** data;
};
typedef PtrToMNode Matrix;
// �ͷ�����ľ���ռ�
void free_Matrix(Matrix mat);
//����������
Matrix Create_Matrix(int row, int column);
// ��ʼ������(������Ԫ�س�ʼ��Ϊ��)
void Init_Matrix(Matrix mat);
// ������ÿ��Ԫ�ظ�ֵ
void SetData_Matrix(Matrix mat, double data[]);
//����ӡ����
void Show_Matrix(Matrix mat, char* s);
//������Ӽ���
Matrix AddorSub_Matrix(Matrix mat_1, Matrix mat_2, int flag);
//������ת��
Matrix Trans_Matrix(Matrix mat);
// ����˷�
Matrix Mult_Matrix(Matrix mat_1, Matrix mat_2);
//��������λ����
Matrix eye(int n);
//��ȡ������ĳ��ĳ�е�Ԫ��
float PickInMat(Matrix mat, int r, int c);
//��˳���˹��ȥ��
Matrix Gauss_shunxu(Matrix A, Matrix b);
//������Ԫ��˹��ȥ��
Matrix Gauss_lie(Matrix A, Matrix b);
// �������
Matrix Cross(Matrix a, Matrix b);

/*
 * �������棬���þ���LU�ֽ����棨ֱ�����Ƿֽ⣩
 * �÷�����Դ��˳���˹��Ԫ����û�н���ѡ��Ԫ�����ܲ�����ֵ���ȶ�������
 */
Matrix LUInverse_Matrix(Matrix mat);

// �������棬���ó����б任����
Matrix EleTransInv_Matrix(Matrix mat);
// �����n�����m�л���
void Swap_row(Matrix mat, int n, int m);

/*
 * ���������ʽ��ֵ��ͨ���б任������Ԫ��ȥ����������任�������Ǿ���
 * ע���б任��ı�����ʽ�ķ���
 */
double Det_Matrix(Matrix mat);
//���������(��ûд)
Matrix Adj_Matrix(Matrix mat);
//����һ��������и���
Matrix Copy_Matrix(Matrix mat);

/*
 * ������ת����Ĺ�ʽ��������α任�������
 * ������ѧר��(Ŀǰ��ûд)
 */
Matrix Rot_inv(Matrix T);


#endif /* HEADER_MYMATRIX_H_ */

#pragma once
