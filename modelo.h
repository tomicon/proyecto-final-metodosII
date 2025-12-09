
#ifndef MODELO_H
#define MODELO_H

#include <complex.h>
#include <tgmath.h> 
// tgmath.h es CRUCIAL: permite que 'pow', 'sin', etc. funcionen 
// automáticamente con números complejos (mapea pow -> cpow)

// Función: sin(z**2+1)
double complex f(double complex z) {
    return sin(pow(z, 2) + 1);
}

// Derivada: 2*z*cos(z**2 + 1)
double complex df(double complex z) {
    return 2*z*cos(pow(z, 2) + 1);
}


#endif // MODELO_H
