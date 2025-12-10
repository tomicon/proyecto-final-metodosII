#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "modelo.h"
#include "calculo_raices.h"

#define LADO_CUADRICULA 800
#define LIMIT 1.5
#define MAX_ITER 100
#define TOL 1e-6
#define TOL_CUADRADA (TOL * TOL)

int newton_raphson(double complex z, const double complex roots[], int roots_count, int *iteraciones); //devuelve el índice de la raíz convergida

int main() {
    double complex z_min = -LIMIT - LIMIT * I;
    double complex z_max = LIMIT + LIMIT * I;
    RootStore store = {.count = 0};

    encontrar_todas_las_raices(z_min, z_max, &store); 

    FILE *fp = fopen("fractal_data.csv", "w");
    if (!fp) {
        fprintf(stderr, "Error al crear el archivo CSV.\n");
        return 1;
    }

    fprintf(fp, "root_idx,iteraciones\n");

    double largo_total = 2.0 * LIMIT;

    // Recorremos la imagen
    for (int y = 0; y < LADO_CUADRICULA; y++) {
        for (int x = 0; x < LADO_CUADRICULA; x++) {
            
            // Mapeo al plano complejo
            double re = -LIMIT + (double)x / (LADO_CUADRICULA - 1) * largo_total;
            double im = -LIMIT + (double)y / (LADO_CUADRICULA - 1) * largo_total;
            double complex z = re + im * I;

            int iteraciones = 0;
            int root_idx;

            // Método de Newton
            root_idx = newton_raphson(z, store.roots, store.count, &iteraciones);
            // Escribir datos en el archivo CSV
            fprintf(fp, "%d,%d\n", root_idx, iteraciones);
        }
    }

    fclose(fp);
    printf("Se creo el archivo 'fractal_data.csv'.\n");
    printf("Raices encontradas: %d\n", store.count);
    return 0;
}

int newton_raphson(double complex z, const double complex roots[], int roots_count, int *iteraciones) {
    for (int k = 0; k < MAX_ITER; k++) {
        double complex deriv = df(z);
        if (cabs(deriv) < 1e-14) return -1;  // Si se aproxima a 0 diverge

        z = z - f(z) / deriv;
        *iteraciones += 1;

        for (int r = 0; r < roots_count; r++) {
            double complex diff = z - roots[r];
            if (creal(diff)*creal(diff) + cimag(diff)*cimag(diff) < TOL_CUADRADA) {
                return r;
            }
        }
    }
    return -1; // Retorna -1 si no convergió después de MAX_ITER
}

