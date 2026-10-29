#include "diagonalization.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

const size_t n = 100;
const double x_max = 10;

int main(){
    double delta = node_dist(x_max, n);
    std::vector<double> d = h_diagonal(delta, x_max, n);
    std::vector<double> e = h_subdiagonal(delta, n);
    std::vector<double> z(n*n, 0.0);


    diagonalization(d, e, z, n);
    eigen_norm(z, delta, n);


    //control print
    std::cout << "E0, E1, E2: " << d[0] << " " << d[1] << " " << d[2] << "\n";
    const std::vector<std::size_t> to_check = {0, 1, 50};
    for(std::size_t k : to_check){
        double sum = 0;
        for(std::size_t i = 0; i < n; i++){
            sum += z[i + k*n]*z[i + k*n]*delta;
        }
        std::cout << "column norm " << k << ": " << sum << "\n";
    }

    std::ofstream spectrum;
    spectrum.open("../data/spectrum.csv");
    if(!spectrum.is_open()){
        std::cerr << "Cannot open spectrum.csv \n";
        return 1;
    }
    spectrum << std::setprecision(12);
    for(std::size_t k = 0; k < n; k++){
        spectrum << k << " " << d[k] << " " << k+0.5 <<"\n";
    }
    spectrum.close();

    // According to a task we are interested in states {1,2,3,15,20,30}
    const std::vector<std::size_t> states = {0, 1, 2, 14, 19, 29};
    std::vector<std::vector<double>> exact(states.size(), std::vector<double>(n));

    for(std::size_t s = 0; s < states.size(); s++){
        const std::size_t k = states[s];
        for(std::size_t i = 0; i < n; i++){
            const double x_i = -x_max + delta*i;
            exact[s][i] = hermite_psi(x_i, k);
        }
        normalize_vector(exact[s].data(), n, delta);

        // check of the sign, dot product with exact and numeric vector (if sign is negative there is a sign change in numeric value)
        double dot = 0;
        for(std::size_t i = 0; i < n; i++){
            dot += z[i + k*n]*exact[s][i];
        }
        if(dot < 0){
            for(std::size_t i = 0; i < n; i++){
                z[i + k*n] = -z[i + k*n];
            }
        }
    }

    std::ofstream eigenvectors;
    eigenvectors.open("../data/eigenvectors.csv");
    if(!eigenvectors.is_open()){
        std::cerr << "Cannot open eigenvectors.csv \n";
        return 1;
    }
    eigenvectors << std::setprecision(12);
    for(std::size_t i = 0; i < n; i++){
        double x_i = -x_max+delta*i;
        eigenvectors << x_i;
        for (std::size_t s = 0; s < states.size(); s++){
            eigenvectors << " " << z[i + states[s]*n] << " " << exact[s][i];
        }
        eigenvectors << "\n";
    }
    eigenvectors.close();


    return 0;
}