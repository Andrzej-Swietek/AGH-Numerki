/**
 * gcc *.c -lgsl -lgslcblas -lm
 * g++ *.c -lgsl -lgslcblas -lm
*/
#include <iostream>
#include <iomanip>
// #include "/usr/include/gsl/gsl_math.h"
// #include "/usr/include/gsl/gsl_linalg.h"

using namespace std;

void displayMatrix(double (*matrix)[4], int columns, int rows, std::string name = "");

int main() {
    int i,j,k;

    // cout << "Podaj Liczbe rownan" << endl;
    // cin >> n;

    constexpr int n = 3;

    /**
     * Deklaracja tablicy n + n+1
    */
    double mat[n][n+1];

    /**
     * Macierz wynikowa 
    */
    double wynik[n];


    //------------------------------------------------------------------------------------------------
    // col = 0                  col = 1                col = 2                  col = 3 - wyrazy wolne
    mat[0][0] = -1;         mat[0][1] = 2;          mat[0][2] = 1;        /*   |  */       mat[0][3] = -1;
    mat[1][0] = 1;         mat[1][1] = -3;          mat[1][2] = -2;       /*   |  */       mat[1][3] = -1;
    mat[2][0] = 3;         mat[2][1] = -1;          mat[2][2] = -1;       /*   |  */       mat[2][3] = 4;


    // ZEROWANIE- METODA ELIMINACJI GAUSA
    
    // ETAP I eliminacja zmiennych i macierz trójkątna
    for ( int i = 0; i < n-1; i++ )
       for ( int j = i+1; j < n; j++ ) 
       {
            double factor = mat[j][i] / mat[i][i];

            for ( k = 0; k < n+1; k++ )
                mat[j][k] -= factor * mat[i][k]; 
       }

    displayMatrix(mat,4, 3);


    // ETAP II wyznaczenie w najniższym wierszu i podstawianie wyzej
    for ( int i = n-1; i >= 0; i-- )
    {
        wynik[i] = mat[i][n];
        for ( int j = i+1; j < n; j++ ) 
            if( i != j ) {
                wynik[i] = wynik[i] - mat[i][j] * wynik[j];
            }

        wynik[i] /=  mat[i][i];
    }


    // ETAP III - Prezentacja wyniku
    for ( int i = 0; i < n; i++ )
            cout << wynik[i] << endl;




    return 0;
}


void displayMatrix(double (*matrix)[4], int columns, int rows, std::string name) {
    int lineLength = columns * 12;
    cout << string(lineLength, '-') << endl;
    if (!name.empty())
        cout << setw(10) << name << ":" << "\n" << string(lineLength, '-') << endl << endl;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            cout << setw(10) << matrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << string(lineLength, '-') << endl;
}