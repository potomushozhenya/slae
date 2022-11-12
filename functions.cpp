#include "iostream"
#include "class.h"
#include "cmath"
double Norm1(Matrix &A) {
    double max = 0;
    double s = 0;
    for (int i = 0; i < A.dimension; i++) {
        for (int j = 0; j < A.dimension; j++) {
            s += std::abs(A.matrix[i][j]);
        }
        if (s > max) {
            max = s;
            s = 0;
        } else {
            s = 0;
        }
    }
    return max;
}
double Norm3(Matrix &A) {
    double max = 0;
    double s = 0;
    for (int i = 0; i < A.dimension; i++) {
        for (int j = 0; j < A.dimension; j++) {
            s += std::abs(A.matrix[j][i]);
        }
        if (s > max) {
            max = s;
            s = 0;
        } else {
            s = 0;
        }
    }
    return max;
}
std::vector<double> operator*(double a, std::vector<double>& vector){
    std::vector<double> c;
    c.resize(vector.size());
    for (int i = 0; i < vector.size(); ++i) {
        c[i]=vector[i]*a;
    }
    return c;
}
std::vector<double> operator+(std::vector<double>&a,std::vector<double>&b){
    std::vector<double> c;
    c.resize(a.size());
    for (int i = 0; i < a.size(); ++i) {
        c[i]=a[i]+b[i];
    }
    return c;
}
std::vector<double> operator-(std::vector<double>&a,std::vector<double>&b){
    std::vector<double> c;
    c.resize(a.size());
    for (int i = 0; i < a.size(); ++i) {
        c[i]=a[i]-b[i];
    }
    return c;
}
double Norm(std::vector<double>&a){
    double s=0;
    for (int i = 0; i < a.size(); ++i) {
        s+=(a[i]*a[i]);
    }
    return std::sqrt(s);
}
