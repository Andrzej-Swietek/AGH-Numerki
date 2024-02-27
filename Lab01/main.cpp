/**
 * gcc *.c -lgsl -lgslcblas -lm
*/
#include <iostream>
// #include "/usr/include/gsl/gsl_math.h"
// #include "/usr/include/gsl/gsl_linalg.h"

using namespace std;

int main() {
    int i,j,k,n;

    // cout << "Podaj Liczbe rownan" << endl;
    // cin >> n;

    n = 3;

    /**
     * Deklaracja tablicy n + n+1
    */
    float mat[n][n+1];

    /**
     * Macierz wynikowa 
    */
    float wynik[n];


    // for ( int row = 0; row < n; row++ )
    //     for ( int col = 0; col < n+1; col++ )
    //     {
    //         float tmp = 0.0;
    //         cin >> tmp;
    //         mat[row][col] = tmp;
    //     }

    //------------------------------------------------------------------------------------------------
    // col = 0                  col = 1                col = 2                  col = 3 - wyrazy wolne
    mat[0][0] = -1;         mat[0][1] = 2;          mat[0][2] = 1;       /*   |  */       mat[0][3] = -1;
    mat[1][0] = 1;         mat[1][1] = -3;          mat[1][2] = -2;       /*   |  */       mat[1][3] = -1;
    mat[2][0] = 3;         mat[2][1] = -1;          mat[2][2] = -1;       /*   |  */       mat[2][3] = 4;


    // ZEROWANIE- METODA ELIMINACJI GAUSA
    // for( int col; col < n-1; col++ )
    // {
    //     for ( int row = 1; row < n; row++)
    //     {
    //         // 0 = mat[col][row] - factor * mat[col][0] => factor = mat[col][row]  / mat[col][0]
    //         float factor = mat[col][row] / mat[col][col];
    //     }
    // }

    // ETAP I eliminacja zmiennych i macierz trójkątna
    for ( int i = 0; i < n-1; i++ )
       for ( int j = i+1; j < n; j++ ) 
       {
            float factor = mat[j][i] / mat[i][i];

            for ( k = 0; k < n+1; k++ )
                mat[j][k] -= factor * mat[i][k]; 
       }


    // ETAP II wyznaczenie w najniższym wierszu i podstawianie wyzej
    for ( int i = n-1; i > 0; i-- )
    {
        wynik[i] = mat[i][n];
        for ( int j = i+1; j < n; j++ ) 
            if( i != j ) 
                wynik[i] = wynik[i] - mat[i][j] * wynik[j];
    }


    // ETAP III - Prezentacja wyniku

    for ( int i = 0; i < n; i++ )
            cout << wynik[i] << endl;




    return 0;
}