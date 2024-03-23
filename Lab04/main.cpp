#include<iostream>
#include<cmath>
#include<vector>
#include <fstream>

using namespace std;
using Matrix = vector<vector<double>>;

void initMatrix(double A[][7], int N)
{
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            A[i][j] = 1.0 / sqrt(2 + abs(i - j));
        }
    }
}

Matrix transpose(const Matrix& A, int n)
{
    Matrix transposed = A;
    for(int i = 0; i<n; i++)
        for(int j= 0; j < n; j++)
            if(i != j) transposed[i][j] = A[j][i];
    return transposed;
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


struct PowerMethodResult {
    vector<double> eigen_values;
    vector<vector<double>> eigen_vector_matrix;
};


PowerMethodResult powerMethod( double W[][7], int Kval, int N, int IT_MAX, ofstream& file)
{
    PowerMethodResult data;

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
        data.eigen_values.push_back(lambda);
        data.eigen_vector_matrix.push_back(x_k);
        file << "\n";

        for(int i =0; i < N; i++)
            for(int j =0; j < N; j++)
                W[i][j] = W[i][j] - lambda * x_k[i] * x_k[j];
    }

    return data;
}


Matrix multiplyMatrices(const Matrix& matrix1, const Matrix& matrix2) {
    int rows1 = matrix1.size();
    int cols1 = matrix1[0].size();
    int rows2 = matrix2.size();
    int cols2 = matrix2[0].size();

    if (cols1 != rows2) {
        cerr << "Error: Matrix dimensions mismatch for multiplication" << endl;
        return vector<vector<double>>();
    }

    vector<vector<double>> result(rows1, vector<double>(cols2, 0.0));

    for (int i = 0; i < rows1; ++i) {
        for (int j = 0; j < cols2; ++j) {
            for (int k = 0; k < cols1; ++k) {
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }

    return result;
}

vector<vector<double>> arrayToVector(const double A[][7], int N) {
    vector<vector<double>> result(N, vector<double>(N));

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            result[i][j] = A[i][j];

    return result;
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


    PowerMethodResult eigen_data = powerMethod(W, Kval, N, IT_MAX, eigenvalues_file);

    cout << "[ Eigen Values ]: \n";
    for ( double val : eigen_data.eigen_values )
        cout << val << " ";
    cout << "\n\n";


    // X^T * A

    Matrix transposed_X = transpose( eigen_data.eigen_vector_matrix , N);

    Matrix D = multiplyMatrices(
            multiplyMatrices( eigen_data.eigen_vector_matrix, arrayToVector(A, N)), // Ponieważ push_back vectory wczesniej to dostajemy macierz juz ztransponowana
            transposed_X
    );

    cout << "[ Diagonal Matrix ]: \n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << D[i][j] << " ";
            matrix_D_file << D[i][j] << " ";
        }
        cout << "\n";
        matrix_D_file << "\n";
    }


    eigenvalues_file.close();
    matrix_D_file.close();

    return 0;
}

