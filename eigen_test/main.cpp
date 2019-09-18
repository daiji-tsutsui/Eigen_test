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
#include <Eigen/Eigenvalues>
//#include <Eigen/Geometry>
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
void Matrix_speed_check(int n);
void Matrix_Copy_Swap_Test();
MatrixXd Matrix_Subroutine_Survival_Test1();
MatrixXd* Matrix_Subroutine_Survival_Test2();
void Matrix_memory_check(int n);
void Array_Calculation_Test();
void Matrix_To_Vector_Test();
void Matrix_Rowwise_Test();
void Matrix_Block_Test();
void Matrix_Block_Subst_Test();
void Matrix_Row_Col_Subst_Test();
void Matrix_Row_Col_Calc_Test();
void Matrix_log_product_Test();
void Matrix_colwise_initialize_Test();
void Matrix_Row_Vector_Test();
void Matrix_Append_Test();
void Matrix_Comp_OP_Test1();
void Matrix_Comp_OP_Test2();

int main(int argc, const char * argv[]) {
    
//    MatrixXf A = MatrixXf::Zero(2,2);
//    A(0,0) = 2;
//    A(1,1) = 5;
//    MatrixXf B(2,2); B.setZero();
//    B(0,1) = 2;
//    B(1,0) = 5;
//    cout << A << endl;
    
    
//    MatrixXd A1 = Matrix_Subroutine_Survival_Test1();
//    MatrixXd* A2 = Matrix_Subroutine_Survival_Test2();
//    PRINT_MAT(A1);
//    PRINT_MAT(*A2);
//    //どちらでも正しく返すことができる -> 解放ができていないのでは？
    
//    Matrix_memory_check(100000);
//    //なぜだかわからないが，メモリの解放はできているらしい．
    
	
    Matrix_Block_Test();
    
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

void Matrix_speed_check(int n){
    int k;
    clock_t start = clock();
    MatrixXd A = MatrixXd::Random(500,500);        // ランダム
    MatrixXd B(500,500);
    for(k = 0; k < n; k++){
        B = A * A;
    }
    clock_t end = clock();
    cout << "time = " << (double)(end - start) / CLOCKS_PER_SEC << "[sec]\n" << "\n";
    getchar();
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

void Matrix_memory_check(int n){
    int k;
    clock_t start = clock();
    MatrixXd A;
    for(k = 0; k < n; k++){
        A = Matrix_Subroutine_Survival_Test1();    //代入処理のときに元のAのメモリは解放される？
    }
    clock_t end = clock();
    cout << "time = " << (double)(end - start) / CLOCKS_PER_SEC << "[sec]\n" << "\n";
    getchar();
}

void Array_Calculation_Test(){
    MatrixXd A = MatrixXd::Random(3,3);
    PRINT_MAT(A);
    PRINT_MAT(A * A);
    
    MatrixXd B;
    B = A.array() * A.array();
    PRINT_MAT(B);
    
    PRINT_MAT(A.array().pow(2.0));
    
    PRINT_MAT(A);
    PRINT_MAT(A.array().min(0.5).max(-0.5));
    PRINT_MAT((A.array() > 0.0));
    
    ArrayXXd X = A;
    PRINT_MAT(X);
    PRINT_MAT(X.min(0.5).max(-0.5));
    PRINT_MAT((X > 0.0));
    
    B = A.transpose() * A;
    PRINT_MAT2(B,"B = A^T * A");
//    PRINT_MAT(B.sqrt());    //Undefined template
    SelfAdjointEigenSolver<MatrixXd> es(B);
    PRINT_MAT(es.operatorSqrt());
    PRINT_MAT(es.operatorSqrt() * es.operatorSqrt());
}

void Matrix_To_Vector_Test(){
    MatrixXd A;
    VectorXd v = VectorXd::Random(3);
	
    PRINT_MAT(v);
    A = v;
    PRINT_MAT(A);
    
    A = MatrixXd::Random(3,2);
    PRINT_MAT(A);
//    v = A;    //error
    Map<VectorXd> w(A.data(), A.size());
    PRINT_MAT(w);
    
    for(int i=0; i<A.cols(); ++i){
        A.col(i) = w.block(A.rows()*i,0,A.rows(),1);
    }
    PRINT_MAT(A);
    
//    MatrixXd B(w);
//    PRINT_MAT(B);
}

void Matrix_Rowwise_Test(){
    MatrixXd A = MatrixXd::Random(3,4);
    VectorXd v(4);
    v << 1, 2, 3, 4;
    
    PRINT_MAT(A);
    A.rowwise() += v.transpose();    // rowwiseで足したければ，横ベクトルでなくてはならない．
    PRINT_MAT(A);
}

void Matrix_Block_Test(){
    MatrixXd m(4,4);
    m <<  1, 2, 3, 4,
        5, 6, 7, 8,
        9,10,11,12,
        13,14,15,16;
    PRINT_MAT(m);
    cout << "Block in the middle" << endl;
    cout << m.block<2,2>(1,1) << endl << endl;
    for (int i = 1; i <= 3; ++i) {
        cout << "Block of size " << i << "x" << i << endl;
        cout << m.block(1,0,i,i) << endl << endl;
    }
}

void Matrix_Block_Subst_Test(){
//    ArrayXXd a = ArrayXXd::Constant(4,4,0.1);    //Arrayでも同様．
    MatrixXd a = MatrixXd::Constant(4,4,0.1);
    PRINT_MAT(a);
    
    MatrixXd m(2,2);
    m <<  1, 2,
        3, 4;
    PRINT_MAT(m);
    
    a.block<2,2>(1,1) = m;
    PRINT_MAT(a);
    
    a.block(0,0,2,3) = a.block(2,1,2,3);
    PRINT_MAT(a);
}

void Matrix_Row_Col_Subst_Test(){
    MatrixXd m(3,3);
    m << 1,2,3,
        4,5,6,
        7,8,9;
    PRINT_MAT2(m,"Here is the matrix m:");
    PRINT_MAT(m.row(1));
    m.col(2) += 3 * m.col(0);
    PRINT_MAT2(m,"m.col(2) += 3 * m.col(0)");
}

void Matrix_Row_Col_Calc_Test(){
	VectorXd v(3); v << 1,2,3;
	VectorXd exp_v = (v.array() - v.maxCoeff()).exp();
	VectorXd sm_v = exp_v / exp_v.sum();
	PRINT_MAT(v);
	PRINT_MAT(exp_v);
	PRINT_MAT(sm_v);
	
	MatrixXd A(3,3);
	A << 1, 2, 3,
		4, 5, 6,
		7, 8, 9;
	MatrixXd exp_A = (A.colwise() - A.rowwise().maxCoeff()).array().exp();
	MatrixXd sm_A = exp_A.array().colwise() * (exp_A.rowwise().sum()).array().inverse();
	PRINT_MAT(A);
	PRINT_MAT(exp_A);
	PRINT_MAT(sm_A);
}

void Matrix_log_product_Test(){
	MatrixXd A = MatrixXd::Random(4,4);
	A = A.array() - A.minCoeff() + 0.1;
	VectorXd v = VectorXd::Random(4);
	v = v.array() - v.minCoeff() + 0.1;
	VectorXd u;
	
	PRINT_MAT(A);
	PRINT_MAT(v);
	PRINT_MAT(A*v);
	
	MatrixXd log_A = A.array().log();
	VectorXd log_v = v.array().log();
	PRINT_MAT(((log_A.rowwise() + log_v.transpose()).array().exp()).rowwise().sum());
}

void Matrix_colwise_initialize_Test(){
	Matrix<double,3,1> v[3];
	for(int i = 0; i < 3; i++){
		v[i] << 1,2,3;
	}
	Matrix<double,3,3> a;
	a << v[0], v[1], v[2];
	PRINT_MAT(a);
}

void Matrix_Row_Vector_Test(){
    MatrixXd A = MatrixXd::Random(3,10);
    PRINT_MAT(A);
    PRINT_MAT(A(5));
    PRINT_MAT(A.block(0,0,1,10));
    PRINT_MAT(A.row(0));
}

void Matrix_Append_Test(){
    MatrixXd A(0,10);
    PRINT_MAT(A);
    
    for(int i=0; i<3; ++i){
        VectorXd v = VectorXd::Random(10);
        A.conservativeResize(A.rows()+1, A.cols());
        A.row(A.rows()-1) = v;
        PRINT_MAT(A);
    }
}

// ``C++で線形代数を''より
void Matrix_Comp_OP_Test1(){
//    // 行列の各成分に対して条件判定して bool 型の行列を返す
//    PRINT_FNC;
//
//    MatrixXd A = MatrixXd::Identity(3,3);
//
//    Matrix&lt;bool, dynamic, dynamic&gt; B = A.array()&lt;1.0; // 比較
//
//    PRINT_MAT(A);
//    PRINT_MAT2(B, "A<1.0");
}

void Matrix_Comp_OP_Test2(){
    /* 行列の各要素が条件をみたすかどうかチェック
     * (all) すべての要素に対して満たすか？
     * (any) どれかひとつの要素に対して満たすか？
     * (count) 条件を満たす要素数
     */
    PRINT_FNC;
    
    MatrixXd A = MatrixXd::Identity(3,3);
    
    PRINT_MAT(A);
    cout << "Does A(i,j)<1.0 hold for all elements? : " << (A.array()<1.0).all() << endl;
    cout << "Does A(i,j)<2.0 hold for all elements? : " << (A.array()<2.0).all() << endl;
    cout << "Does A(i,j)<1.0 hold for any elements? : " << (A.array()<1.0).any() << endl;
    cout << "How many elements do A(i,j)<1.0 hold?  : " << (A.array()<1.0).count() << endl;
    
    cout << endl;
}
