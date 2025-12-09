#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "modelo.h"
#include "calculo_raices.h"

#define WIDTH 800
#define HEIGHT 800
#define LIMIT 1.5
#define MAX_ITER 35
#define TOL 1e-6
#define TOL_SQUARED (TOL * TOL)

int newton_raphson(double complex z, const double complex roots[], int roots_count, int *iterations); //devuelve el índice de la raíz convergida

int main() {
    double complex z_min = -LIMIT - LIMIT * I;
    double complex z_max = LIMIT + LIMIT * I;
    RootStore store = {.count = 0};

    encontrar_todas_las_raices(z_min, z_max, &store); 

    // Abrir archivo CSV en modo texto ("w")
    FILE *fp = fopen("fractal_data.csv", "w");
    if (!fp) {
        fprintf(stderr, "Error al crear el archivo CSV.\n");
        return 1;
    }

    // Escribir cabecera (Header)
    fprintf(fp, "root_idx,iterations\n");

    double largo_total = 2.0 * LIMIT;

    // Recorremos la imagen
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            
            // Mapeo al plano complejo
            double re = -LIMIT + (double)x / (WIDTH - 1) * largo_total;
            double im = -LIMIT + (double)y / (HEIGHT - 1) * largo_total;
            double complex z = re + im * I;

            int iterations = 0;
            int root_idx;

            // Método de Newton
            root_idx = newton_raphson(z, store.roots, store.count, &iterations);
            // Escribir datos en el archivo CSV
            fprintf(fp, "%d,%d\n", root_idx, iterations);
        }
    }

    fclose(fp);
    printf("Proceso completado. Archivo 'fractal_data.csv' creado.\n");
    printf("Raíces encontradas: %d\n", store.count);
    return 0;
}

int newton_raphson(double complex z, const double complex roots[], int roots_count, int *iterations) {
    for (int k = 0; k < MAX_ITER; k++) {
        double complex deriv = df(z);
        if (cabs(deriv) < 1e-14) return -1;

        z = z - f(z) / deriv;
        *iterations += 1;

        int converged = 0;
        for (int r = 0; r < roots_count; r++) {
            double complex diff = z - roots[r];
            if (creal(diff)*creal(diff) + cimag(diff)*cimag(diff) < TOL_SQUARED) {
                return r;
            }
        }
    }
    return -1; // Retorna -1 si no convergió después de MAX_ITER
}

