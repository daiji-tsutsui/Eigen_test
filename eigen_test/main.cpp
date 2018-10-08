//
//  main.cpp
//  eigen_test
//
//  Created by 筒井 大二 on 2018/07/29.
//  Copyright © 2018年 筒井 大二. All rights reserved.
//

#include <iostream>
#include <time.h>
#include <Eigen/Core>
#define PRINT_MAT(X) cout << #X << ":\n" << X << endl << endl
#define PRINT_MAT2(X,DESC) cout << DESC << ":\n" << X << endl << endl
#define PRINT_FNC    cout << "[" << __func__ << "]" << endl
using namespace std;
using namespace Eigen;
void Matrix_Init_Test1();
void Matrix_Init_Test2(int n);
void Matrix_Init_Test3();
void Vector_Init_Test();
void Matrix_RowCol_Major_Test();
void Matrix_speed_check();
void Matrix_Copy_Swap_Test();
MatrixXd Matrix_Subroutine_Survival_Test1();
MatrixXd* Matrix_Subroutine_Survival_Test2();

int main(int argc, const char * argv[]) {
    
//    MatrixXf A = MatrixXf::Zero(2,2);
//    A(0,0) = 2;
//    A(1,1) = 5;
//    MatrixXf B(2,2); B.setZero();
//    B(0,1) = 2;
//    B(1,0) = 5;
    
//    cout << A << endl;
	
	
	MatrixXd A1 = Matrix_Subroutine_Survival_Test1();
	MatrixXd* A2 = Matrix_Subroutine_Survival_Test2();
	PRINT_MAT(A1);
	PRINT_MAT(*A2);
	//どちらでも正しく返すことができる -> 解放ができていないのでは？
	
    return 0;
}

void Matrix_Init_Test1(){
    /* 行列クラスの初期化 */
    Matrix3i A;  // 3x3 の int 型行列
    
    // (1) カンマ演算による初期化
    A << 1,2,3,
        4,5,6,
        7,8,9;
    
    PRINT_MAT(A);
}

void Matrix_Init_Test2(int n){
    // (2) 特殊型のベクトルで割り当て
    MatrixXd A1 = MatrixXd::Zero(n,n);          // すべて0
    MatrixXd A2 = MatrixXd::Ones(n,n);          // すべて1
    MatrixXd A3 = MatrixXd::Constant(n,n,2);    // 定数ベクトル
    MatrixXd A4 = MatrixXd::Random(n,n);        // ランダム
    MatrixXd A5 = MatrixXd::Identity(5,5);      // 単位行列（正方行列でなくてもOK）
    
//    Matrix3d B = Matrix3d::Identity();          // 固定サイズのときはサイズ指定必要なし。
    
    PRINT_MAT(A1);
    PRINT_MAT(A2);
    PRINT_MAT(A3);
    PRINT_MAT(A4);
    PRINT_MAT(A5);
}

void Matrix_Init_Test3(){
    MatrixXd A(3,3);
    for(int i=0;i<3;++i){
        for(int j=0;j<3;++j) A(i,j)=i+j;
    }
    PRINT_MAT(A);
}

void Vector_Init_Test(){
    /* ベクトルクラスの初期化
     * 行列クラスと同様に配列サイズの固定バージョン Vector3d などがある */
    VectorXd u(4);
    
    // (1) カンマ演算による初期化
    u << 1,2,3,4;
    
    PRINT_MAT(u);
    
    // (2) 特殊型のベクトルで割り当て
    int n = 4;
    VectorXd v1 = VectorXd::Zero(n);          // すべて0
    VectorXd v2 = VectorXd::Ones(n);          // すべて1
    VectorXd v3 = VectorXd::Constant(n,2);    // 定数ベクトル
    VectorXd v4 = VectorXd::Random(n);        // ランダム
    VectorXd v5 = VectorXd::LinSpaced(n,1,2); // 線形補間
    VectorXd v6 = VectorXd::Unit(n,2);        // 単位ベクトル
    
    PRINT_MAT(v1);
    PRINT_MAT(v2);
    PRINT_MAT(v3);
    PRINT_MAT(v4);
    PRINT_MAT(v5);
    PRINT_MAT(v6);
}

void Matrix_RowCol_Major_Test(){
    PRINT_FNC;
    
    Matrix<double,2,3,RowMajor> Arow;
    Matrix<double,2,3,ColMajor> Acol;
    Arow << 1,2,3,
            4,5,6;
    Acol << 1,2,3,
            4,5,6;
    
    // 行列の表示
    PRINT_MAT(Arow);
    PRINT_MAT(Acol);
    
    // ストレージの格納順に表示
    cout << "Storage ordering of Arow: " << endl;
    for(int i=0;i<6;++i){
        cout << Arow(i) << " ";
    }
    cout << endl;
    
    cout << "Storage ordering of Acol: " << endl;
    for(int i=0;i<6;++i){
        cout << Acol(i) << " ";
    }
    cout << endl;
    /* 注意！
     * デフォルトは ColMajor になっている．
     * カンマオペレータによる初期化は「ストレージの方向」に依存しない！
     * （常に横方向に伸びる，見た目が直感的になるから？）
     * ストレージの格納方法が異なることは上の Print storage の部分から
     * 確認できる．
     */
    cout << endl;
}

void Matrix_speed_check(){
	int k;
	clock_t start = clock();
	MatrixXd A = MatrixXd::Random(500,500);        // ランダム
	MatrixXd B(500,500);
	for(k = 0; k < 100; k++){
		B = A * A;
	}
	clock_t end = clock();
	cout << "time = " << (double)(end - start) / CLOCKS_PER_SEC << "[sec]\n" << "\n";
}

void Matrix_Copy_Swap_Test(){
	/* コピーおよびスワップ（交換）*/
	PRINT_FNC;
	
	MatrixXd A_orig = MatrixXd::Random(3,3);
	MatrixXd A_dump = A_orig;  /* ディープコピー（参照でない！） */
	
	A_orig(0,0) = 1111;
	
	/* 一方を変更しても他方には影響しない */
	PRINT_MAT(A_orig);
	PRINT_MAT(A_dump);
	
	/* スワッピング！ */
	cout << "  (Swapping A_orig <--> A_dump)" << endl << endl;
	A_orig.swap(A_dump);
	PRINT_MAT(A_orig);
	PRINT_MAT(A_dump);
}

MatrixXd Matrix_Subroutine_Survival_Test1(){
	PRINT_FNC;
	
	MatrixXd A1 = MatrixXd::Random(3,3);
	PRINT_MAT(A1);
	
	return A1;
}

MatrixXd* Matrix_Subroutine_Survival_Test2(){
	PRINT_FNC;
	
	MatrixXd* A2; A2 = new MatrixXd(3,3);
	MatrixXd B = MatrixXd::Random(3,3); *A2 = B;
	PRINT_MAT(*A2);
	
	return A2;
}

