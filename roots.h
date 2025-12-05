#ifndef ROOTS_H
#define ROOTS_H

#include <complex.h>
#include <stdbool.h>

// Configuración de precisión
#define MAX_ROOTS 20       // Límite de raíces a encontrar (ajustable)
#define NEWTON_TOL 1e-9       // Precisión deseada para la raíz
#define SAME_ROOT_TOL 1e-4    // Distancia para considerar dos raíces como "la misma"
#define INTEGRAL_STEPS 50     // Pasos de integración por lado del rectángulo (calidad vs velocidad)

// Estructura para almacenar las raíces encontradas
typedef struct {
    double complex roots[MAX_ROOTS];
    int count;
} RootStore;

void encontrar_todas_las_raices(double complex z_min, double complex z_max, RootStore *store);

#endif