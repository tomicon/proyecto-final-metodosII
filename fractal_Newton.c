#include <stdio.h>
#include <stdlib.h> // Necesario para atoi (texto a entero) y atof (texto a double)
#include <math.h>
#include <complex.h>

#define MAX_ITER 50
#define TOL 1e-6
#define TOL_SQUARED (TOL * TOL)

// Funciones del fractal
double complex f(double complex z) {
    return z * z * z - 1.0;
}
double complex df(double complex z) {
    return 3.0 * z * z;
}

// Main ahora acepta argumentos: argc (cantidad) y argv (valores en texto)
int main(int argc, char *argv[]) {
    // 1. Verificación de seguridad
    // Esperamos 7 argumentos: el nombre del programa + 6 números
    if (argc != 7) {
        fprintf(stderr, "Uso: %s WIDTH HEIGHT X_MIN X_MAX Y_MIN Y_MAX\n", argv[0]);
        return 1;
    }

    // 2. Convertir texto (argv) a números (int/double)
    int width = atoi(argv[1]);
    int height = atoi(argv[2]);
    double x_min = atof(argv[3]);
    double x_max = atof(argv[4]);
    double y_min = atof(argv[5]);
    double y_max = atof(argv[6]);

    const double complex roots[3] = {
        1.0 + 0.0 * I,
        -0.5 + (sqrt(3.0) / 2.0) * I,
        -0.5 - (sqrt(3.0) / 2.0) * I
    };

    FILE *fp = fopen("fractal_data.csv", "w");
    if (!fp) {
        fprintf(stderr, "Error al crear CSV.\n");
        return 1;
    }
    
    // Escribimos la cabecera
    fprintf(fp, "root_idx,iterations\n");

    // 3. Usamos las variables width y height en los bucles
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            
            // 4. Usamos los límites dinámicos (x_min, etc.)
            double re = x_min + (double)x / (width - 1) * (x_max - x_min);
            double im = y_min + (double)y / (height - 1) * (y_max - y_min);
            double complex z = re + im * I;

            int iterations = 0;
            int root_idx = -1;

            for (int k = 0; k < MAX_ITER; k++) {
                double complex deriv = df(z);
                if (cabs(deriv) < 1e-14) break;

                z = z - f(z) / deriv;
                iterations++;

                int converged = 0;
                for (int r = 0; r < 3; r++) {
                    double complex diff = z - roots[r];
                    if (creal(diff)*creal(diff) + cimag(diff)*cimag(diff) < TOL_SQUARED) {
                        root_idx = r;
                        converged = 1;
                        break;
                    }
                }
                if (converged) break;
            }

            fprintf(fp, "%d,%d\n", root_idx, iterations);
        }
    }

    fclose(fp);
    return 0;
}