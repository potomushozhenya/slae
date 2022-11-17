//
// Created by zhest on 02.10.2022.
//
#include "class.h"
#include "iostream"
#include "cmath"
#include "float.h"
#include "functions.h"
#include <iomanip>

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
            std::cout << "a" << i << j <<" ";
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
std::vector<double> simpleIterations( std::vector<double> &vector, Matrix &A,double eps, int& k) {
    k=0;
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
            ++k;
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
            ++k;
            if (Norm(Result)<eps) {
                Flag = false;
            }
        }
    }
    B.del();
    return xK;
}
void Matrix::rowSwap(int a,int b){
    double* tmpRow=this->matrix[a];
    this->matrix[a]=this->matrix[b];
    this->matrix[b]=tmpRow;
}
std::vector<double> Matrix::Y(Matrix &P, std::vector<double> &b){
    std::vector<double> pb=P*b;
    std::vector<double> y;
    y.resize(pb.size());
    double s;
    for (int i = 0; i < this->dimension; ++i) {
        s=0;
        for (int j = 0; j < i; ++j) {
            s+=y[j]*(this->matrix[i][j]);
        }
        y[i]=pb[i]-s;
    }
    return y;
}
std::vector<double> Matrix::X(std::vector<double> &y){
    std::vector<double> x;
    x.resize(y.size());
    double s;
    for (int i = (this->dimension)-1; i > -1; --i) {
        s=0;
        int j;
        for (j = (this->dimension)-1; j > i; --j) {
            s+=x[j]*(this->matrix[i][j]);
        }
        x[i]=(y[i]-s)/(this->matrix[i][j]);
    }
    return x;
}
std::vector<double> LU(Matrix& A,std::vector<double> &b){
    Matrix P;
    P.newMatrix(A.dimension);
    P.E();
    Matrix LU;
    LU.newMatrix(A.dimension);
    LU=A;
    //LU.rowChange(0,1);
    int maxElemIndex=0;
    double maxElem=-1;
    for (int i = 0; i < LU.dimension; ++i) {
        maxElemIndex=0;
        maxElem=-1;
        for (int j = i; j < LU.dimension; ++j) {
            if (std::abs(LU.matrix[j][i])>maxElem){
                maxElem=std::abs(LU.matrix[j][i]);
                maxElemIndex=j;
            }
        }
        if (maxElemIndex!=i){
            LU.rowSwap(i,maxElemIndex);
            P.rowSwap(i,maxElemIndex);
        }
        for (int j = i+1; j < LU.dimension; ++j) {
            LU.matrix[j][i]=(LU.matrix[j][i])/(LU.matrix[i][i]);
            for (int k = i+1; k < LU.dimension; ++k) {
                LU.matrix[j][k]=LU.matrix[j][k]-(LU.matrix[j][i]*LU.matrix[i][k]);
            }
        }
    }
    std::vector<double> y= LU.Y(P,b);
    std::vector<double> x=LU.X(y);
    P.del();
    LU.del();
    return x;
}
bool Matrix::diagonal_predominance(){
    int i;
    double s;
    for (i = 0; i < this->dimension; ++i) {
        s=0;
        for (int j = 0; j < this->dimension; ++j) {
            if(i!=j) {
                s += std::abs(this->matrix[i][j]);
            }
        }
        if(std::abs(this->matrix[i][i])<=s){
            break;
        }
    }
    if(i==(this->dimension)){
        return true;
    } else{
        return false;
    }
}
void Matrix::positively_define(std::vector<double> &b){
    Matrix AT;
    AT.newMatrix(this->dimension);
    AT=this->transposition();
    (*this)=AT*(*this);
    b=AT*b;
    AT.del();
}
std::vector<double> Seidel(Matrix &A,std::vector<double> &b,double eps,int& k){
    if(!A.diagonal_predominance()){
        A.positively_define(b);
    }
    k=0;
    bool Flag=true;
    std::vector<double> xK;
    std::vector<double> d;
    d.resize(A.dimension);
    for (int i = 0; i < d.size(); ++i) {
        d[i]=(b[i])/(A.matrix[i][i]);
    }
    xK=d;
    while(Flag){
        double s;
        ++k;
        for (int i = 0; i < A.dimension; ++i) {
            s=0;
            for (int j = 0; j < A.dimension; ++j) {
                if(i!=j){
                    s+=((-(A.matrix[i][j]))/(A.matrix[i][i]))*xK[j];
                }
            }
            s+=d[i];
            xK[i]=s;
        }
        std::vector<double> tmpVec=A*xK;
        tmpVec=tmpVec-b;
        if (Norm(tmpVec)<eps){
            Flag=false;
        }
    }
    return xK;
}
Matrix doubleW(std::vector<double>w1,std::vector<double>w2){
    Matrix matrix;
    matrix.newMatrix(w1.size());
    for (int i = 0; i < w1.size(); ++i) {
        for (int j = 0; j < w1.size(); ++j) {
            matrix.matrix[i][j]=-2*(w1[i]*w1[j]);
        }
    }
    return matrix;
}
void Matrix::insertion(Matrix& Q){
    int diff=(this->dimension-Q.dimension);
    for (int i = (this->dimension-1); i > diff-1; --i) {
        for (int j = (this->dimension-1); j >diff-1; --j) {
            this->matrix[i][j]=Q.matrix[i-diff][j-diff];
        }
    }
};
std::vector<double> QR(Matrix &A, std::vector<double> &b){
    Matrix R;
    R.newMatrix(A.dimension);
    R=A;
    Matrix QSummary;
    Matrix Q;
    QSummary.newMatrix(A.dimension);
    Q.newMatrix(A.dimension);
    Q.E();
    QSummary=Q;
    std::vector<double> w;
    std::vector<double> y;
    std::vector<double> z;
    y.resize(A.dimension);
    z.resize(A.dimension);
    z[0]=1;
    w.resize(A.dimension);
    for (int i = 0; i < (A.dimension-1); ++i) {
        //find y
        y.resize(A.dimension-i);
        z.resize(A.dimension-i);
        for (int j = (A.dimension-1); j > i-1; --j) {
            y[j-i]=R.matrix[j][i];
            std::cout<<R.matrix[j][i]<<' '<<i<<' '<<j<<std::endl;
        }
        std::cout<<"Vector y"<<std::endl;
        for (int j = 0; j < y.size(); ++j) {
            std::cout<<y[j]<<" ";
        }
        std::cout<<'\n';
        double a=Norm(y);
        std::vector<double>tmpVec= (-a)*z;
        tmpVec=y+tmpVec;
        //y-(a*z)
        double ro=Norm(tmpVec);
        tmpVec=(1/ro)*y;
        w=(-a/ro)*z;
        w=w+tmpVec;
        std::cout<<"Vector w"<<std::endl;
        for (int j = 0; j < w.size(); ++j) {
            std::cout<<w[j]<<" ";
        }
        std::cout<<'\n';
        Matrix E;
        E.newMatrix(A.dimension-i);
        E.E();
        Matrix W;
        Matrix QAdd;
        QAdd.newMatrix(A.dimension);
        QAdd.E();
        W.newMatrix(A.dimension-i);
        W= doubleW(w,w);
        Q=E+W;
        Q.print();
        QAdd.insertion(Q);
        std::cout<<"Matrix QAdd"<<std::endl;
        QAdd.print();
        std::cout<<'\n';
        R=QAdd*R;
        std::cout<<"Matrix R"<<std::endl;
        R.print();
        std::cout<<'\n';
        QSummary=QSummary*QAdd;
        W.del();
        QAdd.del();
        E.del();
    }
    QSummary.print();
    std::cout<<'\n';
    R.print();
    QSummary=QSummary.transposition();
    std::vector<double> ResAdd=QSummary*b;
    std::vector<double> result;
    result.resize(A.dimension);
    for (int i = (A.dimension-1); i > -1 ; --i) {
        double s=0;
        for (int j = (A.dimension-1); j > i ; --j) {
            s+=result[j]*R.matrix[i][j];
        }
        result[i]=(ResAdd[i]-s)/R.matrix[i][i];
    }
    Q.del();
    QSummary.del();
    R.del();
    return result;
}
void Matrix::newMatrix1(std::vector<double> &b,std::vector<double>& solution){
    (*this).newMatrix(3);
    b.resize(this->dimension);
    solution.resize(this->dimension);
    this->matrix[0][0]=0;
    this->matrix[0][1]=2;
    this->matrix[0][2]=3;
    this->matrix[1][0]=1;
    this->matrix[1][1]=2;
    this->matrix[1][2]=4;
    this->matrix[2][0]=4;
    this->matrix[2][1]=5;
    this->matrix[2][2]=6;
    b[0]=13;
    b[1]=17;
    b[2]=32;
    solution[0]=1;
    solution[1]=2;
    solution[2]=3;
}
void Matrix::newMatrix2(std::vector<double> &b,std::vector<double>& solution){
    (*this).newMatrix(3);
    b.resize(this->dimension);
    solution.resize(this->dimension);
    this->matrix[0][0]=7;
    this->matrix[0][1]=1;
    this->matrix[0][2]=1;
    this->matrix[1][0]=1;
    this->matrix[1][1]=9;
    this->matrix[1][2]=1;
    this->matrix[2][0]=1;
    this->matrix[2][1]=1;
    this->matrix[2][2]=11;
    b[0]=9;
    b[1]=11;
    b[2]=13;
    solution[0]=1;
    solution[1]=1;
    solution[2]=1;
}
void Matrix::newMatrix3(std::vector<double> &b,std::vector<double>& solution){
    (*this).newMatrix(3);
    b.resize(this->dimension);
    solution.resize(this->dimension);
    this->matrix[0][0]=-7;
    this->matrix[0][1]=1;
    this->matrix[0][2]=1;
    this->matrix[1][0]=1;
    this->matrix[1][1]=-9;
    this->matrix[1][2]=1;
    this->matrix[2][0]=1;
    this->matrix[2][1]=1;
    this->matrix[2][2]=-11;
    b[0]=-9;
    b[1]=-11;
    b[2]=-13;
    solution[0]=(1.7228915662650602);
    solution[1]=(1.5783132530120481);
    solution[2]=(1.4819277108433734);
}
void Matrix::newMatrix4(std::vector<double> &b,std::vector<double>& solution){
    (*this).newMatrix(3);
    b.resize(this->dimension);
    solution.resize(this->dimension);
    this->matrix[0][0]=-7;
    this->matrix[0][1]=8;
    this->matrix[0][2]=9;
    this->matrix[1][0]=10;
    this->matrix[1][1]=-9;
    this->matrix[1][2]=6;
    this->matrix[2][0]=9;
    this->matrix[2][1]=10;
    this->matrix[2][2]=-11;
    b[0]=9;
    b[1]=11;
    b[2]=13;
    solution[0]=(1.4940029985007496);
    solution[1]=(1.1799100449775112);
    solution[2]=(1.1131934032983508);
}
void Matrix::newMatrix5(std::vector<double> &b,std::vector<double>& solution){
    (*this).newMatrix(3);
    b.resize(this->dimension);
    solution.resize(this->dimension);
    this->matrix[0][0]=7;
    this->matrix[0][1]=6;
    this->matrix[0][2]=6;
    this->matrix[1][0]=6;
    this->matrix[1][1]=9;
    this->matrix[1][2]=6;
    this->matrix[2][0]=6;
    this->matrix[2][1]=6;
    this->matrix[2][2]=11;
    b[0]=9;
    b[1]=11;
    b[2]=13;
    solution[0]=(0.0196078431372549);
    solution[1]=(0.6732026143790849);
    solution[2]=(0.8039215686274509);
}
void Matrix::newMatrix6(int n, double epsilon,std::vector<double>&b,std::vector<double>&solution) {
    (*this).dimension=n;
    (*this).initialization();
    b.resize(this->dimension);
    solution.resize(this->dimension);
    double eps=epsilon*5;
    for (int i = 0; i < this->dimension; ++i) {
        for (int j = 0; j < this->dimension; ++j) {
            if(i<j) {
                this->matrix[i][j] = (-eps) - 1;
                continue;
            }
            if (i>j){
                this->matrix[i][j]=eps;
                continue;
            } else{
                this->matrix[i][j]=eps+1;
            }
        }
    }
    for (int i = 0; i < b.size(); ++i) {
        if(i!=b.size()-1){
            b[i]=-1;
        } else{
            b[i]=1;
        }
    }
    for (int i = 0; i < solution.size(); ++i) {
        if(i!=b.size()-1){
            solution[i]=0;
        } else{
            solution[i]=0.995024875621890547;
        }
    }
}