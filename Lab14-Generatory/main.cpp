#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>
#include<fstream>
#include<sstream>

using namespace std;

struct Vector3D {
    double x, y, z;
};


const unsigned long long M3 = (1ULL << 32) - 5;

double gen_U1() {
    static long int x=10;
    int a = 17;
    int c = 0;
    long int m = pow(2,13)-1;
    x = (a*x + c) % m;

    return (double)x / (m + 1.0);
}

double gen_U2() {
    static long int x=10;
    int a = 85;
    int c = 0;
    long int m = pow(2,13)-1;
    x = (a*x + c) % m;

    return (double)x / (double)(m + 1.0);
}


double gen_U3() {
    static long long x0 = 10, x1 = 10, x2 = 10;
    long long next = (1176 * x0 + 1476 * x1 + 1776 * x2) % M3;
    x2 = x1;
    x1 = x0;
    x0 = next;
    return static_cast<double>(next) / static_cast<double>(M3 + 1.0);
}


void saveToFile(const std::vector<double>& numbers, const std::string& filename) {
    std::ofstream outFile(filename);
    int i =0;
    for (double num : numbers) {
        outFile << i++ << "\t" << num << "\n";
    }
    outFile.close();
}

void saveVectorsToFile(const std::vector<Vector3D>& vectors, const std::string& filename) {
    std::ofstream outFile(filename);
    for (const auto& vec : vectors) {
        outFile << vec.x << " " << vec.y << " " << vec.z << "\n";
    }
    outFile.close();
}



std::vector<Vector3D> generateVectors3D(int N) {
    std::vector<Vector3D> vectors;

    for (int i = 0; i < N; ++i) {
        double u1 = gen_U3();
        double u2 = gen_U3();
        double u3 = gen_U3();
        double u4 = gen_U3();

        double x = sqrt(-2.0 * log(1.0 - u1)) * cos(2 * M_PI * u2);
        double y = sqrt(-2.0 * log(1.0 - u1)) * sin(2 * M_PI * u2);
        double z = sqrt(-2.0 * log(1.0 - u3)) * cos(2 * M_PI * u4);

        vectors.push_back({ x, y, z });
    }

    return vectors;
}


void normalizeToSphereVolume(std::vector<Vector3D>& vectors, int dimensions = 3) {

    for (auto& vec : vectors) {
        double u = gen_U3();
        double s = pow(u, 1.0 / dimensions);
        double length = sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
        vec.x = (vec.x / length) * s;
        vec.y = (vec.y / length) * s;
        vec.z = (vec.z / length) * s;
    }
}


void analyzeDensity(const std::vector<Vector3D>& vectors, int N) {
    const int K = 10;
    std::vector<int> counts(K, 0);
    double delta = 1.0 / K;

    for (const auto& vec : vectors) {
        double r = sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
        int j = static_cast<int>(r / delta);
        if (j < K) {
            counts[j]++;
        }
    }

    std::ofstream outFile("data/density_analysis.dat");
    for (int j = 0; j < K; ++j) {
        double Rj = delta * (j + 1);
        double Rj1 = delta * j;
        double Vj = (4.0 / 3.0) * M_PI * pow(Rj, 3);
        double Vj1 = (4.0 / 3.0) * M_PI * pow(Rj1, 3);
        double V = Vj - Vj1;
        double g = counts[j] / V;
        outFile << j + 1 << " " << counts[j] << " " << g << "\n";
    }
    outFile.close();
}

void save_data_for_generator_chart(vector<double> random, int N, std::string filename){
    std::ofstream outFile1(filename + "r_1.csv");
    outFile1 << "x_1 x_2\n";
    for (int i = 0; i < N-1; i++) {
        outFile1 << random[i] << " " << random[i+1] << "\n";
    }
    outFile1.close();

    std::ofstream outFile2(filename+ "r_2.csv");
    outFile2 << "x_1 x_2\n";
    for (int i = 0; i < N-2; i++) {
        outFile2 << random[i] << " " << random[i+2] << "\n";
    }
    outFile2.close();

    std::ofstream outFile3(filename + "r_3.csv");
    outFile3 << "x_1 x_2\n";
    for (int i = 0; i < N-3; i++) {
        outFile3 << random[i] << " " << random[i+3] << "\n";
    }
    outFile3.close();
}

int main() {
    constexpr int N = 2000;
    std::vector<double> random_U1;
    std::vector<double> random_U2;
    std::vector<double> random_U3;

    for(int i = 0; i < N; i++){
        random_U1.push_back(gen_U1());
        random_U2.push_back(gen_U2());
        random_U3.push_back(gen_U3());
    }

    saveToFile(random_U1, "data/generator_1.csv");
    saveToFile(random_U2, "data/generator_2.csv");
    saveToFile(random_U3, "data/generator_3.csv");
    

    // KULA
    auto vectors = generateVectors3D(N);

    normalizeToSphereVolume(vectors);

    saveVectorsToFile(vectors, "data/vectors_volume.dat");

    analyzeDensity(vectors, N);


    save_data_for_generator_chart(random_U1, N, "data/wykres_gen_1");
    save_data_for_generator_chart(random_U2, N, "data/wykres_gen_2");
    save_data_for_generator_chart(random_U3, N, "data/wykres_gen_3");

    return 0;
}