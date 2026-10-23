#ifndef DIAGONALIZATION_HPP
#define DIAGONALIZATION_HPP
#include <lapacke.h>
#include <cmath>
#include <vector>
#include <iostream>


//delta
inline double node_dist(double x_max, size_t n){
    return (2*x_max)/(n-1);
}

inline std::vector<double> h_diagonal(double node_dist, double x_max, size_t n){
    std::vector<double> h(n, 0.0);
    for(size_t i = 0; i < n; i++){
        double x_i = (-1)*x_max+node_dist*i;
        h[i] = 1/(node_dist*node_dist)+(x_i*x_i)/2;

        return h;
    }
}

inline std::vector<double> h_subdiagonal(double node_dist, size_t n){
    std::vector<double> h(n, 0.0);
    for(size_t i = 0; i < n-1; i++){
        h[i] = (-1)/(2*node_dist*node_dist);
    }

    return h;
}

inline void diagonalization(std::vector<double> &d, std::vector<double> &e, std::vector<double> &z, size_t n){
    const lapack_int info = LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', n, &d[0], &e[0], &z[0], n);
    if(info != 0){
        std::cout << "Succesfully used LAPACKE_dstev routine.";
    }
}

inline void eigen_norm(){
    
}



#endif