#include<iostream>
#include <cmath>
#include<vector>
#include<fstream>
/**
 * metoda jakobiego - iteracyjnie przybizamy dokładnosc rowziazania, w  kazdej iteracji porpawia dokładnosc poprzedniego rozwiania
*/

void writeToFile(const std::vector<double>& t, const std::vector<double>& data, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        int i =0;
        for (double val : data) {
            file << t[i] << " " << val << "\n";
            i++;
        }
        file.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}


int main(){
    // Parameters
    constexpr double h = 0.02f;
    constexpr int N = 1000;
    double V0 = 0.0f;
    double x0 = 1.0f;
    double omega = 1.0f;
    double beta, F0, Omega;

    std::vector<double> t(N, 0); // siatka
    for(int i = 0; i<N; i++) t[i] = h*i;

    // CASE 1
    
    beta  = 0.0;
    F0    = 0.0;
    Omega = 0.8;


    // CASE 2
    beta  = 0.4;
    F0    = 0.0;
    Omega = 0.8;


    // CASE 3
    beta  = 0.4;
    F0    = 0.1;
    Omega = 0.8;

    double a1 = 1;
    double a2 = omega*omega * h*h -2 -beta*h;
    double a3 = 1 + beta*h;


    // Wyrazy Wolne 
    std::vector<double> b(N, 0);
    for( int i = 0; i< N; i++)
        b[i] = F0* std::sin( omega*h*i ) * ( h * h );

    std::vector<double> d0(N, 1.0);
    std::vector<double> d1(N, 1.0);
    std::vector<double> d2(N, 1.0);

    d0[0] =  0.;
    d0[1] =  0.;

    d1[0] =  0.;
    d1[1] = -1.;

    d2[0] =  0.;
    d2[1] =  0.;


    for( int i =2; i < N; i++)
    {
        d0[i] = a3;
        d1[i] = a2;
        d2[i] = a1;
    }
    

    // WEKTOR ROZWIAZAN
    std::vector<double> x_n(N, 0.0); 
    x_n[0] = 1.0;
    x_n[1] = 1.0;

    for( int i =2; i < N; i++)
        x_n[i] = ( 1/d0[i] ) * ( b[i] - d1[i]*x_n[i-1] - d2[i]*x_n[i-2] ); 
    

    writeToFile(t, x_n, "wyniki3.txt");
    
  
    return 0;
}
