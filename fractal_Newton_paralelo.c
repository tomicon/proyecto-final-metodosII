#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <mpi.h>  // Biblioteca MPI
#include "modelo.h"
#include "calculo_raices.h"

#define WIDTH 800
#define HEIGHT 800
#define LIMIT 1.5
#define MAX_ITER 35
#define TOL 1e-6
#define TOL_SQUARED (TOL * TOL)

int newton_rawson(double complex z, const double complex roots[3], int *iterations);

int main(int argc, char *argv[]) {
    int rank, size;
    double start_time, end_time;
    
    // Inicializar MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // ID del proceso
    MPI_Comm_size(MPI_COMM_WORLD, &size);  // Número total de procesos
    
    start_time = MPI_Wtime();  // Tiempo inicial
    
    double complex z_min = -LIMIT - LIMIT * I;
    double complex z_max = LIMIT + LIMIT * I;
    RootStore store = {.count = 0};

    // Solo el proceso maestro busca las raíces
    if (rank == 0) {
        printf("Buscando raíces...\n");
        encontrar_todas_las_raices(z_min, z_max, &store);
        printf("Raices encontradas: %d\n", store.count);
    }
    
    //Como el proceso maestro solo calculo el numero de raices, este lo comparte con los demas procesos
    MPI_Bcast(&store.count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    //MPI NO puede enviar un numero complejo directamente, por lo que separamos las partes real e imaginaria
    double raices_parte_real[MAX_ROOTS];
    double raices_parte_imaginaria[MAX_ROOTS];
    
    if (rank == 0) {
        for (int i = 0; i < store.count; i++) {
            raices_parte_real[i] = creal(store.roots[i]);
            raices_parte_imaginaria[i] = cimag(store.roots[i]);
        }
    }
    
    MPI_Bcast(raices_parte_real, store.count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(raices_parte_imaginaria, store.count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    //Reconstruir las raíces complejas en procesos no-maestros
    if (rank != 0) {
        for (int i = 0; i < store.count; i++) {
            store.roots[i] = raices_parte_real[i] + raices_parte_imaginaria[i] * I;
        }
    }

    double largo_total = 2.0 * LIMIT;
    int total_pixels = WIDTH * HEIGHT;
    
    // *** DISTRIBUCIÓN DEL TRABAJO ***
    // Cada proceso calculará una porción de filas
    int filas_por_procesos = HEIGHT / size;
    int filas_sobrantes = HEIGHT % size;
    
    int indice_inicio, indice_fin;
    if (rank < filas_sobrantes) {
        indice_inicio = rank * (filas_por_procesos + 1);
        indice_fin = indice_inicio + filas_por_procesos + 1;
    } else {
        indice_inicio = rank * filas_por_procesos + filas_sobrantes;
        indice_fin = indice_inicio + filas_por_procesos;
    }
    
    int local_pixels = (indice_fin - indice_inicio) * WIDTH;
    
    
    // Arrays locales para cada proceso
    int *convergencia_raices_locales = (int *)malloc(local_pixels * sizeof(int));
    int *cant_iteraciones_locales = (int *)malloc(local_pixels * sizeof(int));
    
    if (!convergencia_raices_locales || !cant_iteraciones_locales) {
        fprintf(stderr, "Proceso %d: Error al asignar memoria.\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    // *** CALCULO PARALELO (cada proceso trabaja en sus filas) ***
    int local_idx = 0;
    for (int y = indice_inicio; y < indice_fin; y++) {
        for (int x = 0; x < WIDTH; x++) {
            // Mapeo al plano complejo
            double re = -LIMIT + (double)x / (WIDTH - 1) * largo_total;
            double im = -LIMIT + (double)y / (HEIGHT - 1) * largo_total;
            double complex z = re + im * I;

            int iterations = 0;
            int root_idx;

            // Método de Newton
            root_idx = newton_rawson(z, store.roots, &iterations);
            
            convergencia_raices_locales[local_idx] = root_idx;
            cant_iteraciones_locales[local_idx] = iterations;
            local_idx++;
        }
    }
    
    // *** RECOLECCION DE RESULTADOS EN EL PROCESO MAESTRO ***
    int *todas_convergencia_raices = NULL;
    int *final_cant_iteraciones = NULL;
    int *recvcounts = NULL;
    int *displs = NULL;
    
    if (rank == 0) {
        todas_convergencia_raices = (int *)malloc(total_pixels * sizeof(int));
        final_cant_iteraciones = (int *)malloc(total_pixels * sizeof(int));
        recvcounts = (int *)malloc(size * sizeof(int));
        displs = (int *)malloc(size * sizeof(int));
        
        // Calcular cuántos datos recibirá de cada proceso
        int offset = 0;
        for (int i = 0; i < size; i++) {
            int proc_start_row, proc_end_row;
            if (i < filas_sobrantes) {
                proc_start_row = i * (filas_por_procesos + 1);
                proc_end_row = proc_start_row + filas_por_procesos + 1;
            } else {
                proc_start_row = i * filas_por_procesos + filas_sobrantes;
                proc_end_row = proc_start_row + filas_por_procesos;
            }
            recvcounts[i] = (proc_end_row - proc_start_row) * WIDTH;
            displs[i] = offset;
            offset += recvcounts[i];
        }
    }
    
    // Gather de los resultados
    MPI_Gatherv(convergencia_raices_locales, local_pixels, MPI_INT,
                todas_convergencia_raices, recvcounts, displs, MPI_INT,
                0, MPI_COMM_WORLD);
    
    MPI_Gatherv(cant_iteraciones_locales, local_pixels, MPI_INT,
                final_cant_iteraciones, recvcounts, displs, MPI_INT,
                0, MPI_COMM_WORLD);
    
    // ============================================================
    // ESCRITURA DEL ARCHIVO (solo proceso maestro)
    // ============================================================
    if (rank == 0) {
        printf("Proceso maestro escribiendo resultados al archivo...\n");
        
        FILE *fp = fopen("fractal_data.csv", "w");
        if (!fp) {
            fprintf(stderr, "Error al crear el archivo CSV.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        fprintf(fp, "root_idx,iterations\n");
        
        for (int i = 0; i < total_pixels; i++) {
            fprintf(fp, "%d,%d\n", todas_convergencia_raices[i], final_cant_iteraciones[i]);
        }
        
        fclose(fp);
        
        free(todas_convergencia_raices);
        free(final_cant_iteraciones);
        free(recvcounts);
        free(displs);
        
        end_time = MPI_Wtime();
        
        printf("\n=== PROCESO COMPLETADO ===\n");
        printf("Archivo 'fractal_data.csv' creado.\n");
        printf("Tiempo de ejecucion: %.2f segundos\n", end_time - start_time);
        printf("Pixeles procesados: %d\n", total_pixels);
        printf("Procesos MPI: %d\n", size);
        printf("Velocidad: %.0f pixeles/segundo\n", total_pixels / (end_time - start_time));
    }
    
    // Liberar memoria local
    free(convergencia_raices_locales);
    free(cant_iteraciones_locales);
    
    MPI_Finalize();
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
    return -1;
}
