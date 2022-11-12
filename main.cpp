#include <iostream>
#include "class.h"
//#include "functions.cpp"

int main() {
    int var;
    std::cout << "Choose one of 4 options" << std::endl;
    std::cout << "1- Simple Iterations" << std::endl;
    std::cout << "2- Seidel's Method" << std::endl;
    std::cout << "3- LU" << std::endl;
    std::cout << "4- QR" << std::endl;
    std::cin >> var;
    Matrix A;
    A.newMatrix();
    std::vector<double> b;
    b.resize(A.getDimension());
    std::cout<<"Type vector elements ";
    for (double & i : b) {
        std::cin>>i;
    }
    switch (var) {
        case (1): {
            double eps;
            std::cout<<"Type epsilon\t";
            std::cin>>eps;
            std::vector<double> otvet= simpleIterations(b,A,eps);
            for (int i = 0; i < otvet.size(); ++i) {
                std::cout << otvet[i] << " ";
            }
        }
        case (2): {

        }
        case (3): {

        }
        case (4): {

        }
    }
    A.del();
    return 0;
}
