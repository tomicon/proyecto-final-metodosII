# Proyecto Final - Métodos Numéricos II
## Generación de Fractales de Newton

Este proyecto implementa la generacion de fractales de Newton utilizando el metodo de Newton-Raphson para encontrar las raices de funciones complejas. Incluye implementaciones en version serial y paralela (usando MPI) para comparar el rendimiento.

---

## Contenido de los Archivos

### Archivos de Cabecera (.h)

#### `modelo.h`
Contiene la definición de la funcion compleja a analizar y su derivada. Actualmente implementa:
- **Función**: `f(z) = sin(z² + 1)`
- **Derivada**: `f'(z) = 2z·cos(z² + 1)`

Utiliza la biblioteca `<tgmath.h>` para permitir operaciones con numeros complejos de forma automatica (mapea `pow` → `cpow`, `sin` → `csin`, etc.).

#### `roots.h`
Define las estructuras y constantes para el manejo de raices:
- **Constantes**:
  - `MAX_ROOTS`: Límite maximo de raices a encontrar (20)
  - `NEWTON_TOL`: Tolerancia para convergencia de Newton-Raphson (1e-10)
  - `MISMA_ROOT_TOL`: Distancia minima para considerar dos raices como distintas (1e-5)
  - `N_INTEGRACION`: Numero de pasos para integracion numerica (100)
  - `MAX_RECURSION_PROFUNDIDAD`: Profundidad maxima del algoritmo recursivo (20)
- **Estructura `RootStore`**: Almacena las raices encontradas y su contador

#### `calculo_raices.h`
Implementa los algoritmos de busqueda de raices:
- **`newton_refinado()`**: Aplica el metodo de Newton-Raphson para refinar una raiz en una region dada
- **`agregar_root()`**: Agrega una raiz al almacen verificando que sea unica y valida
- **`integral_curvilinea()`**: Calcula la integral de contorno para aplicar el teorema del argumento
- **`rastreo_recursivo()`**: Subdivide recursivamente las regiones del plano complejo para aislar raices usando el principio del argumento
- **`encontrar_todas_las_raices()`**: Funcion principal que combina rastreo inteligente con busqueda exhaustiva en malla fina

### Archivos de Codigo C (.c)

#### `fractal_Newton.c`
Implementacion **serial** del generador de fractales:
- Define una cuadricula de 800x800 píxeles sobre el plano complejo (region: -1.5 a +1.5 en ambos ejes)
- Para cada pixel, aplica el metodo de Newton-Raphson iterativamente
- Registra a qué raiz converge cada punto y el número de iteraciones necesarias
- Genera un archivo CSV (`fractal_data.csv`) con dos columnas:
  - `root_idx`: Índice de la raiz a la que converge
  - `iteraciones`: Numero de iteraciones hasta convergencia

#### `fractal_Newton_paralelo.c`
Implementación **paralela con MPI** del generador de fractales:
- Distribuye el trabajo entre multiples procesos MPI
- El proceso 0 busca las raíces y las distribuye a todos los procesos
- Cada proceso calcula una porcion de filas de la cuadrícula (balance de carga)
- Utiliza **MPI-IO** para escritura paralela del archivo CSV
- Ventajas:
  - Reduce significativamente el tiempo de computo
  - Escala eficientemente con el numero de procesos
  - Ideal para cuadriculas de alta resolucion
- Muestra estadisticas de rendimiento: tiempo total, tiempo de calculo, tiempo de E/S y velocidad en pixeles/segundo

### Notebooks de Jupyter (.ipynb)

#### `generador.ipynb`
Script para generar automaticamente el archivo `modelo.h`:
- Solicita al usuario una función matematica en notación Python (ej: `sin(z**2 + 1)`)
- Utiliza **SymPy** para:
  - Calcular simbolicamente la derivada
  - Convertir las expresiones a codigo C valido
- Genera el archivo `modelo.h` con las funciones `f(z)` y `df(z)`
- Facilita experimentar con diferentes funciones sin editar codigo manualmente

#### `graficar_fractal.ipynb`
Script de visualizacion de los fractales generados:
- Lee el archivo `fractal_data.csv` generado por los programas en C
- Crea una visualizacion colorida del fractal donde:
  - Cada raiz tiene un color distinto (paleta `tab10` de matplotlib)
  - El brillo varia según el numero de iteraciones (convergencia más rapida = más brillante)
- Genera una imagen de alta resolución con ejes etiquetados
- Muestra el numero de raíces detectadas y la resolución de la imagen

---

## Guía de Ejecución Paso a Paso

### Prerrequisitos

#### Para Windows:
1. **Compilador C**: MinGW-w64 o Visual Studio con soporte C
2. **MPI**: Microsoft MPI (MS-MPI)
3. **Python 3.x** con las siguientes bibliotecas:
   ```powershell
   pip install sympy pandas numpy matplotlib
   ```

### Paso 1: Generar el Archivo de Función (Opcional)

Si deseas cambiar la funcion a analizar ejecutar el archivo:

```powershell
generador.ipynb
```

Ejecuta la celda y proporciona tu funcion cuando se solicite (ejemplo: `sin(z**2 + 1)`). Esto generará/actualizará el archivo `modelo.h`.

### Paso 2: Compilar la Versión Serial

```powershell
gcc -o fractal_Newton.exe fractal_Newton.c -lm
```

**Nota**: Si usas MinGW, asegurarse de tener las bibliotecas matematicas complejas.

### Paso 3: Compilar la Versión Paralela

```powershell
mpicc -o fractal_Newton_paralelo.exe fractal_Newton_paralelo.c -lm
```

### Paso 4: Ejecutar el Programa

#### Versión Serial:
```powershell
.\fractal_Newton.exe
```

Esto generará el archivo `fractal_data.csv` con los resultados.

#### Versión Paralela (con 4 procesos):
```powershell
mpirun -np 4 .\fractal_Newton_paralelo.exe
```

Puedes ajustar el numero despues de `-np` segun la cantidad de nucleos de tu CPU.

### Paso 5: Visualizar el Fractal

```powershell
graficar_fractal.ipynb
```

Ejecuta las celdas para generar la visualización del fractal. La imagen mostrara:
- Diferentes colores para cada raíz
- Variaciones de brillo según velocidad de convergencia
- Ejes del plano complejo etiquetados

---


## Detalles Técnicos

### Método de Busqueda de Raices
El proyecto utiliza el **Teorema del Argumento** (tambien conocido como Principio del Argumento):
- Divide el plano complejo en regiones rectangulares
- Calcula integrales de contorno para determinar el número de raices en cada region
- Subdivide recursivamente las regiones que contienen multiples raices
- Refina las raices aisladas usando Newton-Raphson

### Generación del Fractal
- Mapea cada pixel de la imagen a un punto del plano complejo
- Aplica iterativamente el metodo de Newton: `z_{n+1} = z_n - f(z_n)/f'(z_n)`
- Determina a que raíz converge cada punto inicial
- Colorea segun la raíz de convergencia y numero de iteraciones

---

## Autores
Atim, Guadalupe

Castillo, Santiago Emanuel

Consoli, Tomas


Proyecto desarrollado para la materia de "**Métodos Numéricos II**" - Facultad de Ciencias Exacatas y Tecnoliga - UNT