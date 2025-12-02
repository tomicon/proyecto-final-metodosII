#include <stdio.h>
#include <math.h>
#include <complex.h>

#define WIDTH 800
#define HEIGHT 800
#define X_MIN -1.5
#define X_MAX 1.5
#define Y_MIN -1.5
#define Y_MAX 1.5
#define MAX_ITER 50
#define TOL 1e-6
#define TOL_SQUARED (TOL * TOL)

// Funciones del fractal z^3 - 1
double complex f(double complex z) {
    return z * z * z - 1.0;
}
double complex df(double complex z) {
    return 3.0 * z * z;
}

int newton_rawson(double complex z, const double complex roots[3], int *iterations); //devuelve el índice de la raíz convergida

int main() {
    // Definir las 3 raíces
    const double complex roots[3] = {
        1.0 + 0.0 * I,
        -0.5 + (sqrt(3.0) / 2.0) * I,
        -0.5 - (sqrt(3.0) / 2.0) * I
    };

    // Abrir archivo CSV en modo texto ("w")
    FILE *fp = fopen("fractal_data.csv", "w");
    if (!fp) {
        fprintf(stderr, "Error al crear el archivo CSV.\n");
        return 1;
    }

    // Escribir cabecera (Header)
    fprintf(fp, "root_idx,iterations\n");

    // Recorremos la imagen
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            
            // Mapeo al plano complejo
            double re = X_MIN + (double)x / (WIDTH - 1) * (X_MAX - X_MIN);
            double im = Y_MIN + (double)y / (HEIGHT - 1) * (Y_MAX - Y_MIN);
            double complex z = re + im * I;

            int iterations = 0;
            int root_idx;

            // Método de Newton
            root_idx = newton_rawson(z, roots, &iterations);
            // Escribir datos en el archivo CSV
            fprintf(fp, "%d,%d\n", root_idx, iterations);
        }
    }

    fclose(fp);
    printf("Proceso completado. Archivo 'fractal_data.csv' creado.\n");
    return 0;
}

int newton_rawson(double complex z, const double complex roots[3], int *iterations) {
    for (int k = 0; k < MAX_ITER; k++) {
        double complex deriv = df(z);
        if (cabs(deriv) < 1e-14) return -1;

        z = z - f(z) / deriv;
        *iterations += 1;

        int converged = 0;
        for (int r = 0; r < 3; r++) {
            double complex diff = z - roots[r];
            if (creal(diff)*creal(diff) + cimag(diff)*cimag(diff) < TOL_SQUARED) {
                return r;
            }
        }
    }
}

