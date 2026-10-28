#include "diagonalization.hpp"
#include <fstream>
#include <iomanip>

const size_t n = 100;
const double x_max = 10;

int main(){
    double delta = node_dist(x_max, n);
    const std::vector<double> d_original = h_diagonal(delta, x_max, n);
    const std::vector<double> e_original = h_subdiagonal(delta, n);
    std::vector<double> z(n*n, 0.0);

    auto d = d_original, e = e_original;

    diagonalization(d, e, z, n);
    std::ofstream eigenvalues;
    eigenvalues.open("../data/eigenvalues.csv");
    for(std::size_t i = 0; i < n; i++){
        eigenvalues<<d[i]<<"\n";
    }
    eigenvalues.close();
    std::ofstream eigenvectors;
    eigenvectors.open("../data/eigenvectors.csv");
    eigen_norm(z, delta, n);
    for(std::size_t k = 0; k < n; k++){
        for(std::size_t i = 0; i < n; i++){
            eigenvectors<<z[i+k*n]<<"\n";
        }
    }
    eigenvectors.close();


    return 0;
}