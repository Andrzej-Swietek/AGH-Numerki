#include<iostream>
#include<cmath>
#include<vector>
#include <fstream>

using namespace std;

void initMatrix(double A[][7], int N)
{
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            A[i][j] = 1.0 / sqrt(2 + abs(i - j));
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


vector<double> MatrixDotVector(double W[][7], vector<double>& x_k, int N) {
    vector<double> result(N, 0.0);
    for (int x = 0; x < N; ++x)
        for (int y = 0; y < N; ++y)
            result[x] += W[x][y] * x_k[y];
    return result;
}


double VectorDotProduct(vector<double>& x, vector<double>& y, int N) {
    double prod = 0;
    for(int i = 0; i < N; i++) prod += x[i] * y[i];
    return prod;
}


void powerMethod( double W[][7], int Kval, int N, int IT_MAX, ofstream& file)
{
    for(int k =0; k < Kval; k++)
    {
        vector<double> x_k(N, 1.0); // Inicjalizacja wektora startowego
        double lambda = 0;

        for (int i = 0; i < IT_MAX; i++) {
            // Obliczanie iloczynu macierz-wektor
            vector<double> result(N, 0.0);
            result = MatrixDotVector(W, x_k, N);    // x_k+1

            // Obliczenie lambda
            double new_product = VectorDotProduct(result, x_k, N);
            double old_product = VectorDotProduct(x_k, x_k, N);
            lambda = new_product / old_product;     // nowy iloczyn skalany / starny
            old_product = new_product;              // zamiana mianownika

            // Zapisanie wartości lambda do pliku
            file << lambda << " ";

            // Normalizacja wyniku
            normalize(result);

            x_k = result;
        }
        file << "\n";

        for(int i =0; i < N; i++)
            for(int j =0; j < N; j++)
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
    

    powerMethod(W, Kval, N, IT_MAX, eigenvalues_file);

    eigenvalues_file.close();
    matrix_D_file.close();


    return 0;
}

