//
// Created by zhest on 02.10.2022.
//
#include "class.h"
#include "iostream"
#include "cmath"
#include "functions.h"

void Matrix::setDimension() {
    std::cout << "Type dimension od matrix" << std::endl;
    std::cin >> dimension;
}
void Matrix::initialization() {
    matrix = new double *[dimension];
    for (int i = 0; i < dimension; i++) {
        matrix[i] = new double[dimension];
        for (int j = 0; j < dimension; j++) {
            matrix[i][j] = 0;
        }
    }
}
Matrix::fptr Matrix::simpleChoiceNorm() {
    if (Norm1(*this) < 1) {
        return (&Norm1);
    }
    if (Norm3(*this) < 1) {
        return (&Norm3);
    }
    return nullptr;//Если никакая норма не подходит
}
Matrix::fptr Matrix::simpleChoice(Matrix &A,std::vector<double>&b,std::vector<double> &c){
    Matrix E;
    E.dimension=A.dimension;
    E.initialization();
    E.E();
    double mu=0;
    mu=-(1/(Norm3(A)));
    Matrix Res;
    Res.dimension=A.dimension;
    Res.initialization();
    Res=((mu)*(A));
    (*this)=E+Res;
    mu=-mu;
    c =mu*b;
    if (this->simpleChoiceNorm()!=nullptr){
        Res.del();
        E.del();
        return this->simpleChoiceNorm();
    }
    Matrix AT;
    AT.dimension=A.dimension;
    AT.initialization();
    AT=A.transposition();
    b=AT*b;
    A=AT*A;
    mu=-(1/(Norm3(A)));
    Res=((mu)*(A));
    (*this)=E+Res;
    mu=-mu;
    c=mu*b;
    A.print();
    std::cout<<'\n';
    for (int i = 0; i < b.size(); ++i) {
        std::cout<<b[i]<<" ";
    }
    std::cout<<'\n';
    if (this->simpleChoiceNorm()!= nullptr){
        Res.del();
        E.del();
        return this->simpleChoiceNorm();
    }
    return nullptr;
}
void Matrix::fill() {
    std::cout << "Type matrix elements(move by rows)" << std::endl;
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            std::cout << "a" << i << j << ' ';
            std::cin >> matrix[i][j];
        }
    }
}
void Matrix::del() {
    for (int i = 0; i < dimension; i++) {
        delete[]matrix[i];
    }
}

void Matrix::newMatrix() {
    this->setDimension();
    this->initialization();
    this->fill();
}
void Matrix::newMatrix(int n) {
    this->dimension=n;
    this->initialization();
}
void Matrix::print() {
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

Matrix operator+(Matrix &A, Matrix &B) {
    Matrix C;
    C.dimension = A.dimension;
    C.initialization();
    for (int i = 0; i < C.dimension; i++) {
        for (int j = 0; j < C.dimension; j++) {
            C.matrix[i][j] = A.matrix[i][j] + B.matrix[i][j];
        }
    }
    return C;
}
double Matrix::getElement(int a,int b){
    return this->matrix[a][b];
}
Matrix Matrix::transposition() {
    Matrix A;
    A.dimension = this->dimension;
    A.initialization();
    for (int i = 0; i < this->dimension; i++) {
        for (int j = 0; j < this->dimension; j++) {
            A.matrix[i][j] = this->matrix[j][i];
        }
    }
    return A;
}
Matrix operator*(Matrix &A, Matrix &B) {
    Matrix C;
    C.dimension = A.dimension;
    C.initialization();
    double S = 0;
    for (int row = 0; row < C.dimension; row++) {
        for (int intoRow = 0; intoRow < C.dimension; intoRow++) {
            for (int column = 0; column < C.dimension; column++) {
                S += A.matrix[row][column] * B.matrix[column][intoRow];
            }
            C.matrix[row][intoRow] = S;
            S = 0;
        }
    }
    return C;

}
Matrix& Matrix::operator=(const Matrix&B){
    this->dimension=B.dimension;
    for (int i = 0; i < B.dimension; ++i) {
        for (int j = 0; j < B.dimension; ++j) {
            this->matrix[i][j]=B.matrix[i][j];
        }
    }
    return *this;
}
Matrix operator*(double &a, Matrix &A) {
    Matrix C;
    C.dimension = A.dimension;
    C.initialization();
    for (int i = 0; i < C.dimension; i++) {
        for (int j = 0; j < C.dimension; j++) {
            C.matrix[i][j] = a * (A.matrix[i][j]);
        }
    }
    return C;
}
std::vector<double> operator*(Matrix&A,std::vector<double>& vector) {
    double s=0;
    std::vector<double> a;
    a.resize(A.dimension);
    for (int i = 0; i < A.dimension; ++i) {
        for (int j = 0; j < A.dimension; ++j) {
            s=s+((A.matrix[i][j])*vector[j]);
        }
        a[i]=s;
        s=0;
    }
    return a;
}
void Matrix::E(){
    int j =0;
    for (int i = 0; i < this->dimension; ++i,++j) {
        matrix[i][j]=1;
    }
}
int Matrix::getDimension(){
    return this->dimension;
}
std::vector<double> simpleIterations( std::vector<double> &vector, Matrix &A,double eps) {
    Matrix B;
    B.dimension=A.dimension;
    B.initialization();
    bool Flag=true;
    std::vector<double> xk;
    std::vector<double> xK;
    std::vector<double> x0;
    std::vector<double> c;
    Matrix::fptr finalNorm=B.simpleChoice(A,vector,c);
    xk=c;
    if (finalNorm!=nullptr) {
        std::vector<double> xKk;
        double Num= (finalNorm(B) / (1 - finalNorm(B)));
        while (Flag) {
            x0 = (B * xk);
            xK = x0 + c;
            xKk = xK - xk;
            xk=xK;
            if ((Num * Norm(xKk)) < eps) {
                Flag = false;
            }
        }
    }else{
        std::vector<double> Result;
        while (Flag) {
            x0 = (B * xk);
            xK = x0 + c;
            Result=A*xK;
            Result=(Result-vector);
            xk=xK;
            if (Norm(Result)<eps) {
                Flag = false;
            }
        }
    }
    B.del();
    return xK;
}
std::vector<double> LU(Matrix& A,std::vector<double> &b){
    Matrix P;
    P.newMatrix(A.dimension);
    P.E();
    Matrix LU;
    LU=A;
    


}
