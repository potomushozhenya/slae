#include <iostream>
#include "class.h"
#include <iomanip>
//#include "functions.cpp"

int main() {
    int var;
    std::cout << "Choose one of 4 options" << std::endl;
    std::cout << "1- Simple Iterations" << std::endl;
    std::cout << "2- Seidel's Method" << std::endl;
    std::cout << "3- LU" << std::endl;
    std::cout << "4- QR" << std::endl;
    std::cin >> var;
    Matrix P;
    std::vector<double> b;
    P.newMatrix();
    b.resize(P.getDimension());
    std::cout<<"Type vector elements "<<std::endl;
    for (double & i : b) {
        std::cin>>i;
    }
    int k;
    /*int n=10;
    double epsilon=0.000001;
    std::vector<double> solution;
    P.newMatrix6(n,epsilon,b,solution);*/
    switch (var) {
        case (1): {
            double eps;
            std::cout<<"Type epsilon\t";
            std::cin>>eps;
            std::vector<double> otvet= simpleIterations(b,P,eps,k);
            for (int i = 0; i < otvet.size(); ++i) {
                std::cout << std::setprecision(16)<<"x"<<i+1<<"="<<otvet[i] << " ";
            }
            /*std::cout<<std::endl<<"k:"<<k<<std::endl;
            for (int i = 0; i < solution.size(); ++i) {
                std::cout << std::setprecision(16)<<"x"<<i+1<<"="<<std::abs(otvet[i]-solution[i]) << " ";
            }
            std::cout<<std::endl;
            for (int i = 0; i < solution.size(); ++i) {
                std::cout << std::setprecision(16)<<"x"<<i+1<<"="<<solution[i]<< " ";
            }*/
            break;
        }
        case (2): {
            double eps;
            std::cout<<"Type epsilon\t";
            std::cin>>eps;
            std::vector<double> otvet= Seidel(P,b,eps,k);
            for (int i = 0; i < otvet.size(); ++i) {
                std::cout << std::setprecision(16)<<"x"<<i+1<<"="<<otvet[i] << " ";
            }
            std::cout<<std::endl<<"k:"<<k<<std::endl;
            /*for (int i = 0; i < solution.size(); ++i) {
                std::cout << std::setprecision(16)<<"x"<<i+1<<"="<<std::abs(otvet[i]-solution[i]) << " ";
            }
            std::cout<<std::endl;
            for (int i = 0; i < solution.size(); ++i) {
                std::cout << std::setprecision(16)<<"x"<<i+1<<"="<<solution[i]<< " ";
            }*/
            break;
        }
        case (3): {
            std::vector<double> otvet= LU(P,b);
            for (int i = 0; i < otvet.size(); ++i) {
                std::cout << std::setprecision(16)<<"x"<<i+1<<"="<<otvet[i] << " ";
            }
            /*std::cout<<std::endl;
            for (int i = 0; i < solution.size(); ++i) {
                std::cout << std::setprecision(16)<<"x"<<i+1<<"="<<std::abs(otvet[i]-solution[i]) << " ";
            }
            std::cout<<std::endl;
            for (int i = 0; i < solution.size(); ++i) {
                std::cout << std::setprecision(16)<<"x"<<i+1<<"="<<solution[i]<< " ";
            }*/
            break;
        }
        case (4): {
            std::vector<double> otvet= QR(P,b);
            for (int i = 0; i < otvet.size(); ++i) {
                std::cout << std::setprecision(16)<<"x"<<i+1<<"="<<otvet[i] << " ";
            }
            /*std::cout<<std::endl;
            for (int i = 0; i < solution.size(); ++i) {
                std::cout << std::setprecision(16)<<"x"<<i+1<<"="<<std::abs(otvet[i]-solution[i]) << " ";
            }
            std::cout<<std::endl;
            for (int i = 0; i < solution.size(); ++i) {
                std::cout << std::setprecision(16)<<"x"<<i+1<<"="<<solution[i]<< " ";
            }*/
            break;
        }
    }
    P.del();
    return 0;
}
