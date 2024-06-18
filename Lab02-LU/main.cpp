#include<iostream>
#include <vector>
#include<fstream>

double poisson_ro(double x, double Xa, double Xb) {
    if ( x == 0.0 || (x >= -Xb && x < -Xa) || (x > Xa && x <= Xb)) return 0.;
    if ( x >= -Xa && x < 0.0 ) return 1.;
    if ( x > 0.0 && x <= Xa ) return -1.;
    else return 0.0;
}


double V_exact_solution(double x, double Xa, double Xb) {
    if (-Xb <= x && x <= -Xa) {
        return x/16.0 + 1.0/8.0;
    }
    if (-Xa < x && x < 0) {
        return -x*x / 2.0 - (7.0/16.0)*x;
    }
    if (0 <= x && x <= Xa) {
        return x*x / 2.0 - (7.0/16.0)*x;
    } else {
        return x/16.0 - 1.0/8.0;
    }
}



int main() {
	constexpr int N = 500;

	constexpr double Xa = 0.5;
	constexpr double Xb = 2.0;
	constexpr double h = 2.0 * Xb / (double)(N-1);

	double vector_d[N]; // -2/h^2
	double vector_a[N]; // 1/h^2
	double vector_c[N]; // 1/h^2
	double vector_x_i[N];
	double vector_ro[N];

	// Wypełnienie wektora Di
	for(int i = 0; i < N; i++)
		vector_d[i] = -2.0/ (h*h);

	// Wypełnienie wektora Ai i Ci
	for(int i = 0; i < N; i++)
	{
		vector_a[i] = 1.0/ (h*h);
		vector_c[i] = 1.0/ (h*h);
	}

    // Wypełnienie wektora x_i
    for (int i = 0; i < N; i++)
        vector_x_i[i] = -Xb + h * (double)i;                 //         x_i = −Xb + h ∗ (i − 1), i = 1, 2, . . . N


    // Wypełnienie wyrazów wolnych

    for (int i = 0; i < N; i++)
        vector_ro[i] = -poisson_ro(vector_x_i[i], Xa, Xb);



	// Warunki Brzegowe 
	vector_d[0]  = 1;	vector_d[N-1]   = 1;
	vector_c[0]  = 0;	vector_c[N-1]   = 0;
	vector_ro[0] = 0;	vector_ro[N-1]  = 0;


	// 4. Zadanie LU
	double l[N];
	double u[N];
	u[0] = vector_d[0];                 // u_1 = d_1
	for(int i = 1; i< N; i++)           // n = 2,3 ..N
	{
		l[i] = vector_a[i] / u[i-1];
		u[i] = vector_d[i] - l[i] * vector_c[i-1];
	}


//    for(int i = 1; i< N; i++)
//    {
//        std::cout << l[i] << " " << u[i] << std::endl;
//    }




    // Ly = b - b to ro
	/**
	 * y1 = b1
	 * yi = bi - l_i*y_i-1
	*/
	double vector_y[N]; for (int i =0 ; i < N; i++) vector_y[i] = 0;
	vector_y[0] = vector_ro[0];
	for(int i = 1; i < N; i++)
		vector_y[i] = vector_ro[i] - l[i] * vector_y[i-1];   // b_i − l_i *y_{i−1}


    double vector_v[N];
	vector_v[N-1] = vector_y[N-1]/u[N-1];   // v_n = y_n/u_n
	for (int i = (N-1)-1; i >= 0; i--)
        vector_v[i] = (vector_y[i] - vector_c[i]*vector_v[i+1]) / u[i];   // vi = (yi − c_i*v_{i+1})/u_i , i = n − 1, n − 2, . . . , 1



    // ROZWIAZANIE

    std::ofstream outputFile("data/output.txt");

    if (!outputFile.is_open()) {
        std::cerr << "Error opening the file!" << std::endl;
        return 1;
    }

    for ( int i = 0; i < N; i++ )
    {
//        std::cout << "v[ " << i << " ] = " << vector_v[i] << "\n";
        outputFile << i  << " "  << vector_v[i] << "\n";
    }

    outputFile.close();


    std::ofstream solutionFile("data/exact_solution.txt");
    if (!solutionFile.is_open()) { std::cerr << "Error opening the file!" << std::endl;return 1; }
    for ( int i = 0; i < N; i++ )
    {
        solutionFile << i  << " "  << V_exact_solution(vector_x_i[i], Xa, Xb) << "\n";
    }
    solutionFile.close();

	return 0;
}
