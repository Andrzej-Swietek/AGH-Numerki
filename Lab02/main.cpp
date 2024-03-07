#include<iostream>
#include <vector>


constexpr int poisson_ro(double x, double Xa, double Xb) {
	if ( x == 0.0 || (x >= -Xb && x <= -Xa) || (x >= Xa && x <= Xb) ) return 0;
	if ( x >= -Xa && x < 0.0 ) return 1;
	if ( x > 0.0 && x <= Xa ) return -1;
}



int main() {
	constexpr int N = 500;
	double matrix_A[N][N];
	double matrix_B[N][N];
	double result[N];

	constexpr int N = 500;
	constexpr double Xa = 0.5;
	constexpr double Xb = 2.0;
	constexpr double h = 2 * Xb / (N-1);

	double vector_d[N]; // -2/h^2
	double vector_a[N]; // 1/h^2
	double vector_c[N]; // 1/h^2
	double vector_x_i[N];
	double vector_ro[N];

	// Wypełnienie wektora Di
	for(int i = 0; i < N; i++)
		vector_d[i] = -2/ (h*h);

	// Wypełnienie wektora Ai i Ci
	for(int i = 0; i < N; i++)
	{
		vector_a[i] = -1/ (h*h);
		vector_c[i] = -1/ (h*h);
	}

	for (int i=0; i< N; i++)
	{
		// vector_x_i[i] = -Xb + h * (i-1);
		vector_x_i[i] = -Xb + h * (i-1 -1);
	}

	for ( int i = 0; i < N; i++ )
		vector_ro[i] = -poisson_ro((double)i, Xa, Xb); 


	// Warunki Brzegowe 
	vector_d[0]  = 1;	vector_d[N-1] = 1;
	vector_c[0]  = 0;	vector_c[0]   = 0;
	vector_ro[0] = 0;	vector_ro[0]  = 0;


	// 4. Zadanie LU
	double l[N];
	double u[N];
	u[0] = vector_d[0];
	for(int i = 1; i< N; i++)
	{
		l[i] = vector_a[i] / u[i-1];
		u[i] = vector_d[i] - l[i] * vector_c[i-1];
	}


	// Ly = b - b to ro
	/**
	 * y1 = b1
	 * yi = bi - l_i*y_i-1
	*/
	double vector_y[N]; 
	vector_y[0] = vector_ro[0];
	for(int i = 1; i < N; i++)
		vector_y[i] = vector_ro[i] - l[i] * vector_y[i-1];

	double vector_v[N];
	vector_v[N-1] = vector_y[N-1]/u[N-1];
	for (int i = N-1+1; )


	return 0;
}
