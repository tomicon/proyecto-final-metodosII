#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <mpi.h>
#include "modelo.h"
#include "calculo_raices.h"
#include <string.h>

#define LADO_CUADRICULA 800
#define LIMIT 1.5
#define MAX_ITER 100
#define TOL 1e-6
#define TOL_CUADRADA (TOL * TOL) // Nos evita calcular una raíz cuadrada


int newton_raphson(double complex z, const double complex roots[], int roots_count, int *iteraciones);

int main(int argc, char *argv[]) {
    int rank, size;
    double tiempo_de_inicio, tiempo_final, end_time_computacion;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    tiempo_de_inicio = MPI_Wtime();
    
    //Definimos los limites del area del plano complejo
    double complex z_min = -LIMIT - LIMIT * I;
    double complex z_max = LIMIT + LIMIT * I;
    RootStore store = {.count = 0};

    //Solo el proceso 0 realiza la busqueda de raices, luego las comparte a los demas procesos
    if (rank == 0) {
        printf("Buscando raíces...\n");
        encontrar_todas_las_raices(z_min, z_max, &store);
        printf("Raices encontradas: %d\n", store.count);
    }
    
    //Se comparten las raices a los demás procesos
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
    
    // Se reparten las raíces a los demás procesos
    MPI_Bcast(raices_parte_real, store.count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(raices_parte_imaginaria, store.count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    //Reconstruir las raíces complejas en procesos que no son el 0
    if (rank != 0) {
        for (int i = 0; i < store.count; i++) {
            store.roots[i] = raices_parte_real[i] + raices_parte_imaginaria[i] * I;
        }
    }

    double largo_total = 2.0 * LIMIT;
    int total_pixeles = LADO_CUADRICULA * LADO_CUADRICULA;
    
    // Balance de cargas
    // Cada proceso calculara una porcion de filas
    int filas_por_procesos = LADO_CUADRICULA / size;
    int filas_sobrantes = LADO_CUADRICULA % size;
    
    // Se distribuyen las filas sobrantes entre los primeros filas_sobrantes procesos
    int indice_inicio, indice_fin;
    if (rank < filas_sobrantes) {
        indice_inicio = rank * (filas_por_procesos + 1);
        indice_fin = indice_inicio + filas_por_procesos + 1;
    } else {
        indice_inicio = rank * filas_por_procesos + filas_sobrantes;
        indice_fin = indice_inicio + filas_por_procesos;
    }
    
    int cant_pixeles_locales = (indice_fin - indice_inicio) * LADO_CUADRICULA;
    
    //Cada proceso reserva memoria segun la cantidad de pixeles que procesa
    int *convergencia_raices_locales = (int *)malloc(cant_pixeles_locales * sizeof(int));
    int *cant_iteraciones_locales = (int *)malloc(cant_pixeles_locales * sizeof(int));
    
    if (!convergencia_raices_locales || !cant_iteraciones_locales) {
        fprintf(stderr, "Proceso %d: Error al asignar memoria.\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    int indice_local = 0;
    for (int y = indice_inicio; y < indice_fin; y++) {
        for (int x = 0; x < LADO_CUADRICULA; x++) {
            // Mapeo al plano complejo
            double re = -LIMIT + (double)x / (LADO_CUADRICULA - 1) * largo_total;
            double im = -LIMIT + (double)y / (LADO_CUADRICULA - 1) * largo_total;
            double complex z = re + im * I;

            int iteraciones = 0;
            int root_idx;

            // Metodo de Newton
            root_idx = newton_raphson(z, store.roots, store.count, (int *)&iteraciones);
            
            convergencia_raices_locales[indice_local] = root_idx;
            cant_iteraciones_locales[indice_local] = iteraciones;
            indice_local++;
        }
    }
    
    end_time_computacion = MPI_Wtime();  // Tiempo despues del calculo

    // escritura paralela del archivo (cada proceso escribe su parte)    
    char *buffer = (char *)malloc(cant_pixeles_locales * 30 * sizeof(char)); // 30 caractere maximos por línea 
    if (!buffer) {
        fprintf(stderr, "Proceso %d: Error al asignar buffer\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int buffer_pos = 0; //almacena la posicion actual en el buffer. contiene el tamaño total de datos escritos en el buffer
    for (int i = 0; i < cant_pixeles_locales; i++) {

        buffer_pos += sprintf(buffer + buffer_pos, "%d,%d\n", 
                             convergencia_raices_locales[i], 
                             cant_iteraciones_locales[i]);
    }

    // Calcular offsets exactos mediante comunicación colectiva
    int mi_tamanio = buffer_pos;
    int *tamanios_todos = NULL;
    //El proceso 0 reserva memoria para almacenar los tamaños de escritura de todos los procesos
    if (rank == 0) {
        tamanios_todos = (int *)malloc(size * sizeof(int));
    }
    
    //Si no es proceso 0 manda su tamaño de escritura al proceso 0
    MPI_Gather(&mi_tamanio, 1, MPI_INT, tamanios_todos, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    char header[] = "root_idx,iteraciones\n";
    MPI_Offset longitud_header = strlen(header);
            
    // Calcular offsets acumulados (usando MPI_Offset para archivos grandes), donde empieza a escribir cada proceso
    MPI_Offset *offsets = NULL;
    if (rank == 0) {
        offsets = (MPI_Offset *)malloc(size * sizeof(MPI_Offset));
        offsets[0] = longitud_header; //Cabecera del archivo ocupa 20 bytes
        for (int i = 1; i < size; i++) {
            offsets[i] = offsets[i-1] + (MPI_Offset)tamanios_todos[i-1]; //aqui se almacena la posicion de escritura de cada proceso
        }
    } else {
        offsets = (MPI_Offset *)malloc(1 * sizeof(MPI_Offset));
    }
    
    // Distribuir el offset a cada proceso
    MPI_Offset mi_offset;
    MPI_Scatter(offsets, 1, MPI_LONG_LONG, &mi_offset, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD); //se indica a cada proceso donde debe escribir en el archivo
    
    //Abrir archivo con MPI-IO (todos los procesos)
    MPI_File fh;
    MPI_Status status;
    
    int result = MPI_File_open(MPI_COMM_WORLD, "fractal_data.csv",
                               MPI_MODE_WRONLY | MPI_MODE_CREATE,
                               MPI_INFO_NULL, &fh);
    
    if (result != MPI_SUCCESS) {
        fprintf(stderr, "Proceso %d: Error al abrir archivo con MPI_File_open\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    //Escribir cabecera (solo rank 0) y datos (todos)
    if (rank == 0) {
        char header[] = "root_idx,iteraciones\n";
        MPI_Offset longitud_header = strlen(header);
        MPI_File_write_at(fh, 0, header, longitud_header, MPI_CHAR, &status);
    }
    
    //Los procesos escriben sus datos en la posicion correspondiente
    MPI_File_write_at(fh, mi_offset, buffer, buffer_pos, MPI_CHAR, &status);
    
    MPI_File_close(&fh);

    free(buffer);
    if (rank == 0) {
        free(tamanios_todos);
    }
    free(offsets);
    
    if (rank == 0) {
        tiempo_final = MPI_Wtime();
        double tiempo_total = tiempo_final - tiempo_de_inicio;
        double tiempo_computacion = end_time_computacion - tiempo_de_inicio;
        double tiempo_io = tiempo_final - end_time_computacion;
        
        printf("\n=== PROCESO COMPLETADO (ESCRITURA PARALELA MPI-IO) ===\n");
        printf("Archivo 'fractal_data.csv' creado.\n");
        printf("Tiempo TOTAL: %.2f segundos\n", tiempo_total);
        printf("  - Tiempo de CALCULO: %.2f segundos (%.1f%%)\n", 
               tiempo_computacion, 100.0 * tiempo_computacion / tiempo_total);
        printf("  - Tiempo de E/S: %.2f segundos (%.1f%%)\n", 
               tiempo_io, 100.0 * tiempo_io / tiempo_total);
        printf("Pixeles procesados: %d\n", total_pixeles);
        printf("Procesos MPI: %d\n", size);
        printf("Velocidad: %.0f pixeles/segundo\n", total_pixeles / tiempo_total);
    }
    
    free(convergencia_raices_locales);
    free(cant_iteraciones_locales);
    
    MPI_Finalize();
    return 0;
}

int newton_raphson(double complex z, const double complex roots[], int roots_count, int *iteraciones) {
    for (int k = 0; k < MAX_ITER; k++) {
        double complex deriv = df(z);
        if (cabs(deriv) < 1e-14) return -1;

        z = z - f(z) / deriv;
        *iteraciones += 1;

        for (int r = 0; r < roots_count; r++) {
            double complex diff = z - roots[r];
            if (creal(diff)*creal(diff) + cimag(diff)*cimag(diff) < TOL_CUADRADA) {
                return r;
            }
        }
    }
    return -1;
}
