//
// Created by zhest on 02.10.2022.
//

#ifndef SLAE_CLASS_H
#define SLAE_CLASS_H
#include <vector>
class Matrix{
    double** matrix= nullptr;
    int dimension=0;

public:
    void fill();
    void del();
    void setDimension();
    int getDimension();
    void print();
    void initialization();
    void newMatrix();
    void rowSwap(int a,int b);
    friend Matrix operator*(Matrix &A,Matrix&B);
    friend Matrix operator+(Matrix &A,Matrix&B);
    std::vector<double> X(std::vector<double> &y);
    std::vector<double> Y(Matrix &P, std::vector<double> &b);
    friend Matrix operator*(double &a,Matrix&B);
    void E();
    friend double Norm1(Matrix &A);
    friend double Norm3(Matrix &A);
    Matrix& operator=(const Matrix&B);
    Matrix transposition();
    double getElement(int a,int b);
    void newMatrix(int n);
    typedef double(*fptr)(Matrix &A);
    bool diagonal_predominance();
    fptr simpleChoice(Matrix &A,std::vector<double>&b,std::vector<double>&c);
    Matrix::fptr simpleChoiceNorm();
    void positively_define(std::vector<double> &b);
    friend std::vector<double> Seidel(Matrix &A,std::vector<double> &b,double eps, int& k);
    friend std::vector<double> operator*(Matrix&A,std::vector<double>& vector);
    friend std::vector<double> simpleIterations( std::vector<double> &vector, Matrix &matrix,double eps, int& k);
    friend std::vector<double> LU(Matrix& A,std::vector<double> &b);
    friend Matrix doubleW(std::vector<double>w1,std::vector<double>w2);
    void insertion(Matrix&Q);
    void newMatrix1(std::vector<double> &b,std::vector<double>& solution);
    void newMatrix2(std::vector<double> &b,std::vector<double>& solution);
    void newMatrix3(std::vector<double> &b,std::vector<double>& solution);
    void newMatrix4(std::vector<double> &b,std::vector<double>& solution);
    void newMatrix5(std::vector<double> &b,std::vector<double>& solution);
    void newMatrix6(int n, double eps,std::vector<double>&b,std::vector<double>&solution);
    friend std::vector<double> QR(Matrix& A, std::vector<double> &b);
};

#endif //SLAE_CLASS_H
