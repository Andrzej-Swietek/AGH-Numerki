#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>


using namespace std;

#define PI M_PI 
#define TWO_PI (2.0 * M_PI)


double signal_fn(double i, double omega) { return sin(omega*i) + sin(2*omega*i) + sin(3*omega*i); }


void generate_periodic_signal(std::vector<double>& signal, int N) {
     double omega = TWO_PI / (double)N;
    for (int i = 0; i < N; ++i) signal[i] = signal_fn(i, omega);  
}


void add_noise(std::vector<double>& signal, int N) {
    for (int i = 0; i < N; ++i) {
        double delta = 2 * (static_cast<double>(rand()) / RAND_MAX - 0.5);
        signal[i] += delta;
    }
}


void prepare_data_array(const std::vector<double>& signal, std::vector<double>& data) {
    int N = signal.size();
    for (int i = 0; i < N; ++i) {
        data[2 * i] = signal[i];
        data[2 * i + 1] = 0.0; // Imaginary part is zero
    }
}


void filter_signal(std::vector<double>& data, int N, double threshold) {
    for (int i = 0; i < N; ++i) {
        double real_part = data[2 * i];
        double imag_part = data[2 * i + 1];
        double magnitude = std::sqrt(real_part * real_part + imag_part * imag_part);
        if (magnitude < threshold) {
            data[2 * i] = 0.0;
            data[2 * i + 1] = 0.0;
        }
    }
}

void save_signal(std::vector<double> signal, std::string file_name){
    ofstream file("data/fft_" + file_name + ".csv");
    file << "i signal\n";
    for (double i =0; i < signal.size(); i++ )
        file << i << " " << signal[i] << endl;
    
    file.close();
}

void save_result(std::vector<double> data, int N, std::string file_name){
    ofstream file("data/fft_result_" + file_name + ".csv");
    file << "i Re Im\n";
    for (int i = 0; i < N; ++i) {
        std::cout << data[2 * i] << " "; // real part
        std::cout << data[2 * i +1] << "\n"; // imaginary part
        file << i << " " << data[2 * i] << " " << data[2*i+1] << "\n";
    }
    file << endl;
    file.close();
}

int main() {
    srand(static_cast<unsigned int>(time(0)));

    // Define N = 2^k for k = 8, 10, 12
    const int K_VALUES[] = {8, 10, 12};
    const int NUM_K_VALUES = sizeof(K_VALUES) / sizeof(K_VALUES[0]);

    for (int k_idx = 0; k_idx < NUM_K_VALUES; ++k_idx) {
        int k = K_VALUES[k_idx];
        int N = 1 << k; // N = 2^k

        // CLEAN SIGNAL
        std::vector<double> signal(N);
        generate_periodic_signal(signal, N);
        save_signal(signal, "clean_signal_"+to_string(k));

        // NOISED SIGNAL
        add_noise(signal, N);
        save_signal(signal, "noiced_signal_"+to_string(k));

         // Prepare data array for FFT REAL + IMAGINARY
        std::vector<double> data(2 * N); // 2x bigger size (2N)
        prepare_data_array(signal, data);


      // Perform FFT
      gsl_fft_complex_radix2_forward(data.data(), 1, N);

      save_result(data, N, "transformed_"+to_string(k));

      // Apply the threshold filter
        double max_magnitude = 0.0;
        for (int i = 0; i < N; ++i) {
            double real_part = data[2 * i];
            double imag_part = data[2 * i + 1];
            double magnitude = std::sqrt(real_part * real_part + imag_part * imag_part);

            // magnitude -> file

            if (magnitude > max_magnitude) {
                max_magnitude = magnitude;
            }
        }
        double threshold = max_magnitude / 2.0;
        filter_signal(data, N, threshold);

        // INVERSE FFT
        gsl_fft_complex_radix2_backward(data.data(), 1, N);

        // NORMALIZATION OF THE SIGNAL
        for (int i = 0; i < N; ++i) {
            data[2 * i] /= (double)N;
            data[2 * i + 1] /= (double)N;
        }

        save_result(data, N, "wynik"+to_string(k));

    }
    return 0;
}