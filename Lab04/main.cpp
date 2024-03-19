#include<iostream>
#include<cmath>
#include<vector>
#include <fstream>

using namespace std;

void initMatrix(double A[][7], int N)
{
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            A[i][j] = 1 / sqrt(2 + abs(i - j));
        }
    }
}

void transpose(double** A, int n)
{
    for(int i = 0; i<n; i++)
        for(int j= 0; j < n; j++)
            if(i != j) A[i][j] = A[j][i];
  
}

double norm(const vector<double>&v)
{
    double sum = 0.0;
    for(double val : v) {
        sum += val * val;
    }
    return sqrt(sum);
}

void normalize(vector<double>& v)
{
    double n = norm(v);
    for( double& val : v)
    {
        val /= n;
    }
}

void printMatrix(double A[][7], int N)
{
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
}


void powerMethod(double A[][7], double W[][7], int Kval, int N, int IT_MAX, ofstream& file)
{
    for(int k =0; k < Kval; k++)
    {
        vector<double> x_k(N, 1.0); // Inicjalizacja wektora startowego
        double lambda = 0;

        for (int i = 0; i < IT_MAX; i++) {
            // Obliczanie iloczynu macierz-wektor
            vector<double> result(N, 0.0);
            for (int x = 0; x < N; ++x)
                for (int y = 0; y < N; ++y)
                    result[x] += W[x][y] * x_k[y];

            // Normalizacja wyniku
            normalize(result);

            // Obliczenie lambda
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j)
                    lambda += W[i][j] * x_k[i] * x_k[j];

            // Zapisanie wartości lambda do pliku
            file << lambda << " ";

            x_k = result;
        }

        for(int i =0; i < N; ++i)
            for(int j =0; j < N; ++j)
                W[i][j] = W[i][j] - lambda * x_k[i] * x_k[j];
    }
}


int main(){
    constexpr int N = 7;
    constexpr int Kval = N;
    constexpr int IT_MAX = 12;
    
    ofstream eigenvalues_file("eigenvalues.txt");
    ofstream matrix_D_file("matrix_D.txt");

    double A[N][N];
    double W[N][N];

    initMatrix(A,N);
    initMatrix(W,N);

    cout << "Matrix A:" << endl;
    printMatrix(A, N);
    cout << "\n";
    

    powerMethod(A, W, Kval, N, IT_MAX, eigenvalues_file);

    eigenvalues_file.close();
    matrix_D_file.close();


    return 0;
}

