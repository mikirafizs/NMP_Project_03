#include "diagonalization.hpp"
#include <fstream>
#include <iomanip>

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

    const std::vector<size_t> stages = {0,1,2,14,19,29}; // According to a task we are interested in stages {1,2,3,15,20,30}
    std::ofstream eigenvectors;
    eigenvectors.open("../data/eigenvectors.csv");
    if(!eigenvectors.is_open()){
        std::cerr << "Cannot open eigenvectors.csv \n";
        return 1;
    }
    for(std::size_t i = 0; i < n; i++){
        double x_i = -x_max+delta*i;
        eigenvectors << x_i;
        for (std::size_t k : stages){
            eigenvectors << " " << z[i+k*n];
        }
        eigenvectors << "\n";
    }
    eigenvectors.close();


    return 0;
}