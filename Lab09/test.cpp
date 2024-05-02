#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;
#include "/usr/include/gsl/gsl_math.h"
#include "/usr/include/gsl/gsl_linalg.h"
#include "/usr/include/gsl/gsl_eigen.h"

double f(double x){
    return cos(x);
}

int N = 4;
int M = 4;

int main(void){

    int n = N+M;

    gsl_vector *ck = gsl_vector_calloc(n+1);
    double silnia = 1;
    double minus = 1.0;
    for(int i=0;i<n;i++){
        double temp = 0.0;
        if(i%2==0){
            temp = minus;
            if(minus>0){minus = -1.0;}
            else{minus = 1.0;}
        }
        else{
            temp = 0.0;
        }
        // cout<<temp<<" ";
        if(i!=0){
            silnia = silnia*i;
        }
        
        // cout<<silnia<<" ";
        gsl_vector_set(ck, i, temp/silnia);

    }
    cout<<endl;

    // for(int i=0;i<n;i++){
    //     cout<<gsl_vector_get(ck, i)<<" ";
    // }

    cout<<endl;

    gsl_matrix *A = gsl_matrix_calloc(M, M);
    gsl_vector *y = gsl_vector_calloc(M);
    int signum;

    for(int i=0;i<M;i++){
        for(int j=0;j<M;j++){
            gsl_matrix_set(A, i, j, gsl_vector_get(ck, N-M+i+j+1));
        }
        double tempo = gsl_vector_get(ck, N+1+i);
        gsl_vector_set(y, i, -tempo);
    }

    gsl_permutation *p = gsl_permutation_alloc(M);
    gsl_vector *x = gsl_vector_calloc(M);
    gsl_linalg_LU_decomp(A, p, &signum);
    gsl_linalg_LU_solve(A, p, y, x);
    // for(int i=0;i<n;i++){
    //     cout<<gsl_vector_get(ck, i)<<" ";
    // }
    // cout<<endl;

    for(int i=0;i<M;i++){
        cout<<gsl_vector_get(x, i)<<" ";
    }
    cout<<endl;

    double a[N];
    cout<<"a: ";
    for(int i=0;i<N;i++){
        double suma = 0.0;
        for(int j=0;j<=i;j++){
            suma+= gsl_vector_get(ck, i-j)*gsl_vector_get(x, j);
        }
        a[i] = suma;
        cout<<a[i]<<" ";
    }
    cout<<endl;
    
    double b[M];
    b[0] = 1;
    cout<<"b: ";
    for(int i=M-1;i>=1;i--){
        b[i] = gsl_vector_get(x, i);
        // cout<<b[i]<<" ";
    }
    for(int i=0;i<M;i++){

        cout<<b[i]<<" ";
    }
    cout<<endl;

    ofstream file;
    file.open("sx.txt");


    double xmax = 5;
    double xmin = -5;

    double h = (xmax - xmin)/(100);

    double sumaa = 0;
    double sumab = 0;
    for(int i=0;i<100;i++){
        sumaa=0;
        sumab=0;
        double x = xmin + h*(i-1);
        for(int j=0;j<N;j++){
            sumaa+= a[j]*pow(x, j);
            sumab+= b[j]*pow(x,j);
        }
        double f = sumaa/sumab;
        cout<<f<<" ";
        file<<x<<" "<<f<<endl;
        // cout<<x<<endl;
    }
    cout<<endl;

    return 0;
}