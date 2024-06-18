#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>
#include<fstream>
#include<sstream>

#include <cmath>
#include <random>
#include <iomanip>

// Parametry
const double k = 1.38e-23;                  // Stała Boltzmanna
const double T = 100.0;                     // Temperatura w Kelvinach
const double m = 40.0 * 1.66e-27;           // Masa cząsteczki (Argon) w kg
const double sigma = std::sqrt(k * T / m);  // Sigma dla rozkładu normalnego

double generateNormal(std::mt19937 &gen, std::uniform_real_distribution<> &dist) {
    double x1 = dist(gen);
    double x2 = dist(gen);
    return sigma * std::sqrt(-2.0 * std::log(x1)) * std::cos(2.0 * M_PI * x2);
}


// Funkcja do generowania prędkości o rozkładzie Maxwella
std::vector<double> generateMaxwell(int NL, std::mt19937 &gen, std::uniform_real_distribution<> &dist) {
    std::vector<double> velocities(NL);
    for (int i = 0; i < NL; ++i) {
        double Vx = generateNormal(gen, dist);
        double Vy = generateNormal(gen, dist);
        double Vz = generateNormal(gen, dist);
        velocities[i] = std::sqrt(Vx * Vx + Vy * Vy + Vz * Vz);
    }
    return velocities;
}


double calculateMean(const std::vector<double> &velocities) {
    double sum = 0.0;
    for (double v : velocities) {
        sum += v;
    }
    return sum / velocities.size();
}



double calculateStandardDeviation(const std::vector<double> &velocities, double mean) {
    double sum = 0.0;
    for (double v : velocities) {
        sum += (v - mean) * (v - mean);
    }
    return std::sqrt(sum / (velocities.size() - 1));
}



void createHistogram(const std::vector<double> &velocities, int bins, double maxVelocity) {
    std::vector<int> histogram(bins, 0);
    double binSize = maxVelocity / bins;
    for (double v : velocities) {
        int bin = std::min(static_cast<int>(v / binSize), bins - 1);
        histogram[bin]++;
    }

    std::cout << "Histogram:\n";
    for (int i = 0; i < bins; ++i) {
        double binStart = i * binSize;
        double binEnd = binStart + binSize;
        // std::cout << std::setw(5) << (binStart + binEnd)/2<< ";" << < std::setw(5) << histogram[i] << '\n';
        std::cout << "[" << std::setw(5) << binStart << ", " << std::setw(5) << binEnd << "): " << histogram[i] << '\n';
    }
}

int main() {
    int NL = 10000; 
    int bins = 30;  // Liczba przedziałów histogramu
    double maxVelocity = 5.0 * sigma;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0.0, 1.0);

    std::vector<double> velocities = generateMaxwell(NL, gen, dist);

    createHistogram(velocities, bins, maxVelocity);

    double mean = calculateMean(velocities);
    double stddev = calculateStandardDeviation(velocities, mean);

    std::cout << "Średnia prędkość: " << mean << " m/s\n";
    std::cout << "Odchylenie standardowe: " << stddev << " m/s\n";

    return 0;
}