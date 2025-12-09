#ifndef ROOTS_H
#define ROOTS_H

#include <complex.h>
#include <stdbool.h>

#define MAX_ROOTS 20       // Límite de raíces a encontrar
#define NEWTON_TOL 1e-10      // Precisión deseada para la raíz
#define MISMA_ROOT_TOL 1e-5    // Distancia para considerar dos raíces como "la misma"
#define N_INTEGRACION 100    // Pasos de integración por lado del rectángulo
#define MAX_RECURSION_PROFUNDIDAD 20  // Profundidad máxima de recursión

// Estructura para almacenar las raíces encontradas
typedef struct {
    double complex roots[MAX_ROOTS];
    int count;
} RootStore;

void encontrar_todas_las_raices(double complex z_min, double complex z_max, RootStore *store);

#endif