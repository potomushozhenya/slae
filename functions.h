//
// Created by zhest on 19.10.2022.
//
#include "class.h"
#ifndef SLAE_FUNCTIONS_H
#define SLAE_FUNCTIONS_H
double Norm1(Matrix &A);
double Norm3(Matrix &A);
std::vector<double> operator*(double a, std::vector<double>& vector);
std::vector<double> operator+(std::vector<double>&a,std::vector<double>&b);
std::vector<double> operator-(std::vector<double>&a,std::vector<double>&b);
double Norm(std::vector<double>&a);

#endif //SLAE_FUNCTIONS_H
