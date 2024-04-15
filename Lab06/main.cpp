#include <stdio.h>
#include<iostream>
#include <math.h>
// #include "/usr/include/gsl/gsl_math.h"
// #include "/usr/include/gsl/gsl_linalg.h"
// #include "/usr/include/gsl/gsl_eigen.h"


/**
 * Wyznaczanie pierwiastków równań nie liniowych 
 * metoda iteracyjna: Metoda siecznych
 * f(x) = x^4 - 7.899x^3 + 23.281114x^2 - 30.33468152x+ 14.73866033
 * 
 * Kometarz:
 * Przy odpowiednio małym kroku metody mogą znalezc pierwiastki podwojne
*/



double f(double x) {
    double a = 1.0;
    double b = -7.899;
    double c = 23.281114;
    double d = -30.33468152;
    double e = 14.73866033;

    return a*(x*x*x*x) + b * (x*x*x) + c*(x*x) + d*x + e; 
}

double dfdx( double x ) {
    const double Delta_x = 0.001;

    return ( f(x + Delta_x) - f(x - Delta_x) ) 
                    / (2*Delta_x);
}

double u( double x ) {
    return f(x) / dfdx(x);
}

void Find_Solution(double x1, double x2, double x3, FILE* file)
{

    int iter_count = 0;
    double epsilon = abs(x3 - x1);
    while( abs(x3 - x1) > pow(10, -6)) {
        double tmp1 = x2;
        double tmp2 = x3;

        epsilon = abs(x3 - x1);
        x3 = x1 - f(x1)*(x1-x2)/(f(x1) - f(x2)); 

        x1 = tmp1;
        x2 = tmp2;
        fprintf(file, "%d\t%lf\t%lf\n", iter_count, x3, epsilon);
        iter_count++;
    }
    std::cout<<x3<<std::endl;
}


void Find_Solution_Modified(double x1, double x2, double x3, FILE* file)
{

    int iter_count = 0;
    double epsilon = abs(x3 - x1);
    while( abs(x3 - x1) > pow(10, -6)) {
        double tmp1 = x2;
        double tmp2 = x3;

        epsilon = abs(x3 - x1);
        x3 = x1 - u(x1)*(x1-x2)/(u(x1) - u(x2)); 

        x1 = tmp1;
        x2 = tmp2;
        fprintf(file, "%d\t%lf\t%lf\n", iter_count, x3, epsilon);
        iter_count++;
    }
    std::cout<<x3<<std::endl;
}





int main() {

    // f(x) --> Pod wykres
    // FILE *fp_function_values = fopen("data/fx.txt", "w");
    // fprintf(fp_function_values, "x\tf_x\n");
    // for (double x = 1.5; x < 2.4; x+= 0.01)
    //     fprintf(fp_function_values, "%lf\t%lf\n", x, f(x));

    // fclose(fp_function_values);


/// Pierwszy Pierwsiastek

    FILE *fp_solutions_1 = fopen("data/pierwiastek_1.csv", "w");
    double x1 = 1.5; double x2 = 1.75; double x3 = 1.65;
    fprintf(fp_solutions_1, "iteracja\tx\tepsilon\n");

    Find_Solution(x1, x2, x3, fp_solutions_1);

    fclose(fp_solutions_1);

/// Drugi Pierwiastek

    FILE *fp_solutions_2 = fopen("data/pierwiastek_2.csv", "w");
    fprintf(fp_solutions_2, "iteracja\tx\tepsilon\n");

 
    x1 = 1.85;
    x2 = 2.0;
    x3 = 1.9;
    Find_Solution(x1, x2, x3, fp_solutions_2);

    fclose(fp_solutions_2);

/// Trzeci Pierwiastek
   FILE *fp_solutions_3 = fopen("data/pierwiastek_3.csv", "w");
    fprintf(fp_solutions_3, "iteracja\tx\tepsilon\n");


    x1 = 2.0;
    x2 = 2.3; 
    x3 = 2.4;
    Find_Solution(x1, x2, x3, fp_solutions_3);

    fclose(fp_solutions_3);


    /// Metoda Zmodyfikowana


std::cout << "===================== Modified Method =====================\n";

/// Pierwszy Pierwsiastek

    FILE *fp_solutions_1_m = fopen("data/pierwiastek_1_m.csv", "w");
    x1 = 1.6; 
    x2 = 1.8; 
    x3 = 1.7;
    fprintf(fp_solutions_1, "iteracja\tx\tepsilon\n");

    Find_Solution_Modified(x1, x2, x3, fp_solutions_1_m);

    fclose(fp_solutions_1_m);

/// Drugi Pierwiastek

    FILE *fp_solutions_2_m = fopen("data/pierwiastek_2_m.csv", "w");
    fprintf(fp_solutions_2_m, "iteracja\tx\tepsilon\n");

 
    x1 = 1.85;
    x2 = 2.0;
    x3 = 1.9;
    Find_Solution_Modified(x1, x2, x3, fp_solutions_2_m);

    fclose(fp_solutions_2_m);

/// Trzeci Pierwiastek
   FILE *fp_solutions_3_m = fopen("data/pierwiastek_3_m.csv", "w");
    fprintf(fp_solutions_3_m, "iteracja\tx\tepsilon\n");


    x1 = 2.0;
    x2 = 2.3; 
    x3 = 2.4;
    Find_Solution_Modified(x1, x2, x3, fp_solutions_3_m);

    fclose(fp_solutions_3_m);


    return 0;
}
