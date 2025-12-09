#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "roots.h"
#include "modelo.h"

double complex newton_refine(double complex z_guess) {
    double complex z = z_guess;
    for (int i = 0; i < 50; i++) { 
        // CAMBIO AQUI: f y df
        double complex f_val = f(z);
        double complex df_val = df(z);
        
        if (cabs(df_val) < 1e-14) break; 
        
        double complex diff = f_val / df_val;
        z = z - diff;
        
        if (cabs(diff) < NEWTON_TOL) return z;
    }
    return z;
}

// 2. Agregar raíz única (Sin cambios, solo lógica de arrays)
void add_unique_root(double complex root, RootStore *store) {
    if (store->count >= MAX_ROOTS) return;

    for (int i = 0; i < store->count; i++) {
        if (cabs(root - store->roots[i]) < SAME_ROOT_TOL) {
            return; 
        }
    }
    
    if (isfinite(creal(root)) && isfinite(cimag(root))) {
        store->roots[store->count] = root;
        store->count++;
    }
}

// 3. Integración Numérica: df(z)/f(z)
double complex line_integral(double complex p1, double complex p2) {
    double complex sum = 0;
    double complex step = (p2 - p1) / (double)INTEGRAL_STEPS;
    
    double complex z = p1;
    for (int i = 0; i < INTEGRAL_STEPS; i++) {
        double complex z_next = z + step;
        
        // CAMBIO AQUI: f y df con protección contra singularidades
        double complex f_val1 = f(z);
        double complex f_val2 = f(z_next);
        
        // Evitar división por cero
        if (cabs(f_val1) > 1e-15 && cabs(f_val2) > 1e-15) {
            double complex val1 = df(z) / f_val1;
            double complex val2 = df(z_next) / f_val2;
            sum += (val1 + val2) * 0.5 * step;
        }
        z = z_next;
    }
    return sum;
}

// 4. Recursión (Sin cambios estructurales, llama a funciones internas)
void scan_recursive(double complex z_min, double complex z_max, RootStore *store, int depth) {
    if (depth > MAX_RECURSION_DEPTH) {
        double complex center = (z_min + z_max) / 2.0;
        add_unique_root(newton_refine(center), store);
        return;
    }
    
    // Si la región es muy pequeña, refinar directamente
    double region_size = cabs(z_max - z_min);
    if (region_size < 1e-4) {
        double complex center = (z_min + z_max) / 2.0;
        add_unique_root(newton_refine(center), store);
        return;
    }

    double complex c1 = z_min; 
    double complex c2 = creal(z_max) + cimag(z_min)*I; 
    double complex c3 = z_max; 
    double complex c4 = creal(z_min) + cimag(z_max)*I; 

    double complex integral = 0;
    integral += line_integral(c1, c2); 
    integral += line_integral(c2, c3); 
    integral += line_integral(c3, c4); 
    integral += line_integral(c4, c1); 

    double complex N_complex = integral / (2.0 * M_PI * I);
    int N = (int)round(creal(N_complex));
    
    // Validar que el número de raíces sea razonable
    double imag_part = fabs(cimag(N_complex));
    if (imag_part > 0.5) {
        // Integral con mucho ruido, subdividir más
        N = -1;
    }

    if (N <= 0) {
        return; 
    } 
    else if (N == 1) {
        double complex center = (z_min + z_max) / 2.0;
        double complex refined = newton_refine(center);
        // Agregar la raíz si está dentro de la región o cerca de ella
        add_unique_root(refined, store);
    } 
    else {
        double mid_r = (creal(z_min) + creal(z_max)) / 2.0;
        double mid_i = (cimag(z_min) + cimag(z_max)) / 2.0;
        double complex mid = mid_r + mid_i * I;

        scan_recursive(z_min, mid, store, depth + 1); 
        scan_recursive(creal(mid) + cimag(z_min)*I, creal(z_max) + cimag(mid)*I, store, depth + 1); 
        scan_recursive(creal(z_min) + cimag(mid)*I, creal(mid) + cimag(z_max)*I, store, depth + 1); 
        scan_recursive(mid, z_max, store, depth + 1); 
    }
}

// Función Pública con método de respaldo
void encontrar_todas_las_raices(double complex z_min, double complex z_max, RootStore *store) {
    store->count = 0;
    scan_recursive(z_min, z_max, store, 0);
    
    // Método de muestreo en cuadrícula (complementario)
    // Se ejecuta siempre para asegurar que no se pierdan raíces
    int initial_count = store->count;
    int grid_size = 20;  // Aumentado para mejor cobertura
    double re_min = creal(z_min);
    double re_max = creal(z_max);
    double im_min = cimag(z_min);
    double im_max = cimag(z_max);
    
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
            double re = re_min + (re_max - re_min) * i / (grid_size - 1);
            double im = im_min + (im_max - im_min) * j / (grid_size - 1);
            double complex guess = re + im * I;
            double complex root = newton_refine(guess);
            
            // Verificar que sea realmente una raíz
            if (cabs(f(root)) < 1e-6) {
                add_unique_root(root, store);
            }
            
            // Detener si ya encontramos el máximo
            if (store->count >= MAX_ROOTS) {
                return;
            }
        }
    }
}