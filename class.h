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
    friend Matrix operator*(Matrix &A,Matrix&B);
    friend Matrix operator+(Matrix &A,Matrix&B);
    friend Matrix operator*(double &a,Matrix&B);
    void E();
    friend double Norm1(Matrix &A);
    friend double Norm3(Matrix &A);
    Matrix& operator=(const Matrix&B);
    Matrix transposition();
    double getElement(int a,int b);
    void newMatrix(int n);
    typedef double(*fptr)(Matrix &A);
    fptr simpleChoice(Matrix &A,std::vector<double>&b,std::vector<double>&c);
    Matrix::fptr simpleChoiceNorm();
    friend std::vector<double> operator*(Matrix&A,std::vector<double>& vector);
    friend std::vector<double> simpleIterations( std::vector<double> &vector, Matrix &matrix,double eps);
    friend std::vector<double> LU(Matrix& A,std::vector<double> &b);
};

#endif //SLAE_CLASS_H
