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

    printf("Generando CSV (%d pixeles)...\n", WIDTH * HEIGHT);

    // Recorremos la imagen
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            
            // Mapeo al plano complejo
            double re = X_MIN + (double)x / (WIDTH - 1) * (X_MAX - X_MIN);
            double im = Y_MIN + (double)y / (HEIGHT - 1) * (Y_MAX - Y_MIN);
            double complex z = re + im * I;

            int iterations = 0;
            int root_idx = -1;

            // Método de Newton
            for (int k = 0; k < MAX_ITER; k++) {
                double complex deriv = df(z);
                if (cabs(deriv) < 1e-14) break;

                z = z - f(z) / deriv;
                iterations++;

                int converged = 0;
                for (int r = 0; r < 3; r++) {
                    double complex diff = z - roots[r];
                    // Distancia al cuadrado para optimizar
                    if (creal(diff)*creal(diff) + cimag(diff)*cimag(diff) < TOL_SQUARED) {
                        root_idx = r;
                        converged = 1;
                        break;
                    }
                }
                if (converged) break;
            }

            // Escribir en el CSV: índice de la raíz y número de iteraciones
            // Formato: 0,15 (ejemplo)
            fprintf(fp, "%d,%d\n", root_idx, iterations);
        }
    }

    fclose(fp);
    printf("Proceso completado. Archivo 'fractal_data.csv' creado.\n");
    return 0;
}