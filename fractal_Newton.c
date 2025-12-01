#include <stdio.h>
#include <stdlib.h> // Necesario para atoi (texto a entero) y atof (texto a double)
#include <math.h>
#include <complex.h>

#define MAX_ITER 50 // maximo de iteraciones
#define TOL 1e-6 // tolerancia para convergencia de raiz
#define TOL_SQUARED (TOL * TOL) // cuadtado de la tolerancia, para comparar, para evitar calcular raiz

// Funciones del fractal
// funcion
double complex f(double complex z) {
    return z * z * z - 1.0;
}
// derivada de la funcion
double complex df(double complex z) {
    return 3.0 * z * z;
}

// Main ahora acepta argumentos: argc (cantidad) y argv (valores en texto)
int main(int argc, char *argv[]) {
    // verificacion de seguridad
    // Esperamos 7 argumentos: el nombre del programa + 6 números
    if (argc != 7) {
        fprintf(stderr, "Uso: %s WIDTH HEIGHT X_MIN X_MAX Y_MIN Y_MAX\n", argv[0]);
        return 1;
    }

    // convierte los arg de texto (argv) a números (int/double)
    int width = atoi(argv[1]);
    int height = atoi(argv[2]);
    double x_min = atof(argv[3]);
    double x_max = atof(argv[4]);
    double y_min = atof(argv[5]);
    double y_max = atof(argv[6]);

    // arreglo con las raices de la funcion
    const double complex roots[3] = {
        1.0 + 0.0 * I,
        -0.5 + (sqrt(3.0) / 2.0) * I,
        -0.5 - (sqrt(3.0) / 2.0) * I
    };

    // abre el archivo en modo escritura
    FILE *fp = fopen("fractal_data.csv", "w");
    if (!fp) {
        fprintf(stderr, "Error al crear CSV.\n");
        return 1;
    }
    
    // escribimos la cabecera
    fprintf(fp, "root_idx,iterations\n");

    // usamos las variables width y height en los bucles
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            
            // usamos los límites dinámicos (x_min, etc.)
            double re = x_min + (double)x / (width - 1) * (x_max - x_min); // parte real del numero
            double im = y_min + (double)y / (height - 1) * (y_max - y_min); // parte imaginaria del numero
            double complex z = re + im * I; // define el complejo

            int iterations = 0; // contador de iteraciones
            int root_idx = -1; // indice de raiz -1, no convergencia

            // metodo de newton
            for (int k = 0; k < MAX_ITER; k++) {
                double complex deriv = df(z);
                if (cabs(deriv) < 1e-14) break; // si la derivada es casi cero, detiene la iteracion

                z = z - f(z) / deriv; // formula de newton
                iterations++;

                int converged = 0; // controlar si el punto converge

                // recorre las tres raices definidas
                for (int r = 0; r < 3; r++) {
                    double complex diff = z - roots[r]; // diferencia entre punto actual y raiz r
                    
                    // comprueba la convergencia
                    if (creal(diff)*creal(diff) + cimag(diff)*cimag(diff) < TOL_SQUARED) {
                        root_idx = r;
                        converged = 1;
                        break;
                    }
                }
                if (converged) break;
            }

            fprintf(fp, "%d,%d\n", root_idx, iterations); // resultados por pixel
        }
    }

    fclose(fp);
    return 0;
}