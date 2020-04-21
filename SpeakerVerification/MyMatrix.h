/*
 * MyMatrix.h
 *
 *  Created on: Jul 13, 2019
 *      Author: xuuyann
 */

#ifndef HEADER_MYMATRIX_H_
#define HEADER_MYMATRIX_H_
 /*
  * 问题记录：
  *　* 存在内存滥用的现象，因为每次创建矩阵都会申请一块新内存，并没有free掉
  *　* LU分解求逆不太熟
  */
typedef struct MNode* PtrToMNode;
struct MNode
{
    int row;
    int column;
    float** data;
};
typedef PtrToMNode Matrix;
// 释放申请的矩阵空间
void free_Matrix(Matrix mat);
//　创建矩阵
Matrix Create_Matrix(int row, int column);
// 初始化矩阵(将所有元素初始化为０)
void Init_Matrix(Matrix mat);
// 给矩阵每个元素赋值
void SetData_Matrix(Matrix mat, double data[]);
//　打印矩阵
void Show_Matrix(Matrix mat, char* s);
//　矩阵加减法
Matrix AddorSub_Matrix(Matrix mat_1, Matrix mat_2, int flag);
//　矩阵转置
Matrix Trans_Matrix(Matrix mat);
// 矩阵乘法
Matrix Mult_Matrix(Matrix mat_1, Matrix mat_2);
//　创建单位矩阵
Matrix eye(int n);
//　取出矩阵某行某列的元素
float PickInMat(Matrix mat, int r, int c);
//　顺序高斯消去法
Matrix Gauss_shunxu(Matrix A, Matrix b);
//　列主元高斯消去法
Matrix Gauss_lie(Matrix A, Matrix b);
// 向量叉乘
Matrix Cross(Matrix a, Matrix b);

/*
 * 矩阵求逆，利用矩阵LU分解求逆（直接三角分解）
 * 该方法来源于顺序高斯消元法，没有进行选主元，可能产生数值不稳定的现象
 */
Matrix LUInverse_Matrix(Matrix mat);

// 矩阵求逆，利用初等行变换求逆
Matrix EleTransInv_Matrix(Matrix mat);
// 矩阵第n行与第m行互换
void Swap_row(Matrix mat, int n, int m);

/*
 * 求矩阵行列式的值，通过行变换（列主元消去法）将矩阵变换成上三角矩阵
 * 注意行变换会改变行列式的符号
 */
double Det_Matrix(Matrix mat);
//　伴随矩阵(还没写)
Matrix Adj_Matrix(Matrix mat);
//　对一个矩阵进行复制
Matrix Copy_Matrix(Matrix mat);

/*
 * 利用旋转矩阵的公式来计算齐次变换矩阵的逆
 * 机器人学专用(目前还没写)
 */
Matrix Rot_inv(Matrix T);


#endif /* HEADER_MYMATRIX_H_ */

#pragma once
