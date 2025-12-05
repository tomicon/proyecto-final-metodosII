
#ifndef MODELO_H
#define MODELO_H

#include <complex.h>
#include <tgmath.h> 
// tgmath.h es CRUCIAL: permite que 'pow', 'sin', etc. funcionen 
// automáticamente con números complejos (mapea pow -> cpow)

// Función: z**3-1
double complex f(double complex z) {
    return pow(z, 3) - 1;
}

// Derivada: 3*z**2
double complex df(double complex z) {
    return 3*pow(z, 2);
}


#endif // MODELO_H
