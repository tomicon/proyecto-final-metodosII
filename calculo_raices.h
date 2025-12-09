#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "roots.h"
#include "modelo.h"

//utilizada para encontrar una raíz en una región dada
double complex newton_refinado(double complex z_inicial) {
    double complex z = z_inicial;
    for (int i = 0; i < 50; i++) { 
        double complex f_val = f(z);
        double complex df_val = df(z);
        
        if (cabs(df_val) < 1e-14) break; // Si la derivada se aproxima a 0 el método diverge
        
        double complex diff = f_val / df_val;
        z = z - diff;
        
        if (cabs(diff) < NEWTON_TOL) return z;
    }
    return z; //Si no converge, retorna el último valor calculado ya que sabemos que había una raíz en esa región
}

// Agrega una raíz al almacén si es única y válida
void agregar_root(double complex root, RootStore *store) {
    if (store->count >= MAX_ROOTS) return;

    for (int i = 0; i < store->count; i++) {
        if (cabs(root - store->roots[i]) < MISMA_ROOT_TOL) {
            return; 
        }
    }
    
    if (isfinite(creal(root)) && isfinite(cimag(root))) {
        store->roots[store->count] = root;
        store->count++;
    }
}


double complex integral_curvilinea(double complex p1, double complex p2) {
    double complex sum = 0;
    double complex paso = (p2 - p1) / (double)N_INTEGRACION;
    
    double complex z = p1;
    for (int i = 0; i < N_INTEGRACION; i++) {
        double complex z_siguiente = z + paso;
        
        double complex f_val1 = f(z);
        double complex f_val2 = f(z_siguiente);
        
        // Evitar división por cero
        if (cabs(f_val1) > 1e-15 && cabs(f_val2) > 1e-15) {
            double complex val1 = df(z) / f_val1;   // f'(z) / f(z) es la función a integrar según el principio del argumento
            double complex val2 = df(z_siguiente) / f_val2;
            sum += (val1 + val2) * 0.5 * paso;
        }
        z = z_siguiente;
    }
    return sum;
}

// Se busca recursivamente las regiones que contienen raíces, dividiendolas hasta aislar cada raíz
void rastreo_recursivo(double complex z_min, double complex z_max, RootStore *store, int profundidad) {
    if (profundidad > MAX_RECURSION_PROFUNDIDAD) {
        double complex center = (z_min + z_max) / 2.0;
        agregar_root(newton_refinado(center), store);
        return;
    }
    
    // Si la región es muy pequeña, refinar directamente
    double tamanio_region = cabs(z_max - z_min);
    if (tamanio_region < 1e-4) {
        double complex centro = (z_min + z_max) / 2.0;
        agregar_root(newton_refinado(centro), store);
        return;
    }

    double complex c1 = z_min; 
    double complex c2 = creal(z_max) + cimag(z_min)*I; 
    double complex c3 = z_max; 
    double complex c4 = creal(z_min) + cimag(z_max)*I;

    //Se divide la integral en los 4 lados del rectángulo para encontrar el número de raíces dentro
    double complex integral = 0;
    integral += integral_curvilinea(c1, c2); 
    integral += integral_curvilinea(c2, c3); 
    integral += integral_curvilinea(c3, c4); 
    integral += integral_curvilinea(c4, c1); 

    double complex N_complejo = integral / (2.0 * M_PI * I);
    int N = (int)round(creal(N_complejo));   //nos aseguramos de redondear a un número entero
    
    
    double parte_imaginaria = fabs(cimag(N_complejo));  //Esto nos permite saber si el calculo de la integral tiene sentido
    if (parte_imaginaria > 0.5) {
        // Integral con mucho ruido, subdividir más
        N = -1;
    }
    /*Cabe aclarar que en esta parte del código tal vez se esté descartando una región con raíces si es que también había polos
    Pero como este proyecto está pensado para raíces que no tengan polos, no debería haber problema*/
    if (N <= 0) { 
        return; 
    } 
    else if (N == 1) {   // Si hay exactamente una raíz, refinar con Newton y agregar la raiz encontrada
        double complex centro = (z_min + z_max) / 2.0;
        double complex root = newton_refinado(centro);
        agregar_root(root, store);
    } 
    else {
        double medio_r = (creal(z_min) + creal(z_max)) / 2.0;
        double medio_i = (cimag(z_min) + cimag(z_max)) / 2.0;
        double complex medio = medio_r + medio_i * I;

        // Dividir en 4 sub-regiones y rastrear recursivamente
        rastreo_recursivo(z_min, medio, store, profundidad + 1); 
        rastreo_recursivo(creal(medio) + cimag(z_min)*I, creal(z_max) + cimag(medio)*I, store, profundidad + 1); 
        rastreo_recursivo(creal(z_min) + cimag(medio)*I, creal(medio) + cimag(z_max)*I, store, profundidad + 1); 
        rastreo_recursivo(medio, z_max, store, profundidad + 1); 
    }
}

// función de seguridad por si al rastreo inteligente se le escapó alguna raíz
// divide la cuadrícula en una malla fina y prueba cada punto
void encontrar_todas_las_raices(double complex z_min, double complex z_max, RootStore *store) {
    store->count = 0;
    rastreo_recursivo(z_min, z_max, store, 0);
    
    int count_inicial = store->count;
    int tamanio_cuadricula = 20;
    double re_min = creal(z_min);
    double re_max = creal(z_max);
    double im_min = cimag(z_min);
    double im_max = cimag(z_max);
    
    for (int i = 0; i < tamanio_cuadricula; i++) {
        for (int j = 0; j < tamanio_cuadricula; j++) {
            double re = re_min + (re_max - re_min) * i / (tamanio_cuadricula - 1);
            double im = im_min + (im_max - im_min) * j / (tamanio_cuadricula - 1);
            double complex prueba = re + im * I;
            double complex root = newton_refinado(prueba);
            
            if (cabs(f(root)) < 1e-6) {
                agregar_root(root, store);
            }
            
            if (store->count >= MAX_ROOTS) {
                return;
            }
        }
    }
}