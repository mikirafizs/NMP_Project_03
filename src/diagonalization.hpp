#ifndef DIAGONALIZATION_HPP
#define DIAGONALIZATION_HPP
#include <lapacke.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <cstddef>



//delta
inline double node_dist(double x_max, size_t n){
    return (2*x_max)/(n-1);
}

inline std::vector<double> h_diagonal(double delta, double x_max, size_t n){
    std::vector<double> h(n, 0.0);
    for(size_t i = 0; i < n; i++){
        double x_i = (-1)*x_max+delta*i;
        h[i] = 1/(delta*delta)+(x_i*x_i)/2;
    }
    return h;
}

inline std::vector<double> h_subdiagonal(double delta, size_t n){
    std::vector<double> h(n, 0.0);
    for(size_t i = 0; i < n-1; i++){
        h[i] = (-1)/(2*delta*delta);
    }
    return h;
}

inline void diagonalization(std::vector<double> &d, std::vector<double> &e, std::vector<double> &z, size_t n){
    if(z.size() != n*n){
        std::cout<<"Size does not match.";
        return;
    }
    const lapack_int info = LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', n, &d[0], &e[0], &z[0], n);
    if(info != 0){
        std::cout << "Error in using LAPACKE";
        return;
    }
}

inline void eigen_norm(std::vector <double> &z, double delta, std::size_t n){
    for(std::size_t k = 0; k < n; k++){
        double sum = 0;
        for(std::size_t i = 0; i < n; i++){
            sum += z[i + k*n]*z[i + k*n]*delta;
        }
        double c = std::sqrt(sum);
        for(std::size_t i = 0; i < n; i++){
            z[i + k*n] /= c;
        }
    }
}

inline double hermite_psi(double x, std::size_t k){
    double a = 0, b = 1, c = 1;
    for(size_t i = 1; i <= k; i++){
        c = 2*x*b-2*(i-1)*a;
        a = b;
        b = c;
    }
    return c*std::exp((-1*x*x)/2);
}


#endif