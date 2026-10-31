## A.1 El mecanismo de la falla
**Universidad católica del Maule** \
**Profesor:** Juan Cabezas Contardo

#### 1.¿Que es un **pivote** en la eliminación Gaussiana y qué papel juega en el calculo de los multiplicadores?.

El pivote es el número utilizado para hacer ceros todos los demás elementos que estan por debajo de él en la misma columna. Este además es el denominador en la formula de los multiplicadores:

$$
m_{ij} = \frac{a_{ij}}{a_{jj}} \tag{1}
$$

#### 2. ¿Qué pasa cuando el pivote es exactamente cero? ¿Es lo mismo que cuando es muy pequeño pero distinto de cero?.

A partir de la función ya definida (1), podemos notar que si el $a_{ii}$ vale exactamente 0, el multiplicador se indefine por división por cero. Además, si el pivote fuera muy pequeño pero no cero, podemos notar como el multiplicador se hace extremadamente grande lo cual dificulta el cálculo a mano y en entornos donde se usen cifras significativas puede generar indeterminación o fallas de precisión. Burden y Faires ilustran este caso con un pivote de $0.003$ y aritmética de 4 cifras, donde la solución sale incorrecta pese a que el sistema no es singular [1, sec. 6.2]. Por lo tanto, no son lo mismo: con pivote cero el algoritmo falla de forma explícita (error aritmético), mientras que con pivote pequeño el algoritmo "funciona" pero entrega un resultado erróneo sin avisar (pérdida de precisión). Esta segunda situación es más peligrosa precisamente porque no genera ningún error visible.

#### 3. Define cancelación catastrófica y muestra con el sistema del problema motivador en qué línea del algoritmo ocurre.

**Cancelación catastrófica:** Es un fenómeno en este caso de la computación en donde una operación aritmetica de dos números genera una pérdida drástica de digitos significativos en el resultado.

En el algoritmo `src/pivoteo.py` ocurriría en la función r5(x), la cual reduce las cifras significativas de los cálculos, esto con el fin de forzar el error como si se hiciera un cálculo a mano o un cálculo limitado.



## A.2 Las estrategias

#### 4. Describe el pivoteo parcial (partial pivoting): qué busca, dónde busca y qué intercambia.

El pivoteo parcial es otro tipo de pivoteo que busca eliminar el problema que teníamos en (1) cuando el pivote vale 0, este lo que hace es cambiar la fila de la diagonal que vale 0, con la fila que tenga el valor mas grande en su misma columna, así:

$$
A =\begin{pmatrix} 
0 & 1 \\ 1 & 1 
\end{pmatrix}, \qquad
b = 
\begin{pmatrix} 1 \\ 2 
\end{pmatrix}
$$
Pasa a:
$$
A = \begin{pmatrix} 
1 & 1 \\ 0 & 1 
\end{pmatrix}, \qquad
b = 
\begin{pmatrix} 2 \\ 1 
\end{pmatrix}
$$
Como podemos ver, se intercambiaron las filas F1 ⇔ F2, lo cual también afecta a las filas de $b$.

#### 5. Describe el pivoteo completo (full pivoting) y el pivoteo parcial escalado (scaled partial pivoting). Para cada uno: ¿qué problema resuelve que el parcial no resuelve, y qué cuesta de más?

**Pivoteo completo (full pivoting):** En el paso $k$ no se limita a la columna $k$, sino que busca el elemento de mayor valor absoluto en **toda la submatriz restante** (filas $i \geq k$ y columnas $j \geq k$), y lo lleva a la posición $a_{kk}$ intercambiando tanto filas como columnas:

$$
|a_{pq}| = \max_{i,j \geq k} |a_{ij}| \tag{2}
$$

Como intercambiar columnas equivale a intercambiar el orden de las incógnitas, hay que llevar un registro de esa permutación para reordenar $x$ al final.

- **Qué resuelve que el parcial no:** el parcial solo compara dentro de una columna, así que puede quedar con un pivote relativamente pequeño si el mejor candidato está en otra columna. El completo garantiza el pivote más grande posible en cada paso, lo que da la mejor estabilidad numérica y permite detectar el rango de la matriz.
- **Qué cuesta de más:** buscar el máximo en la submatriz de tamaño $(n-k+1)^2$ en cada paso suma del orden de $n^3/3$ comparaciones en total, es decir, del mismo orden que las operaciones aritméticas de la eliminación. Además hay que llevar la permutación de columnas.

**Pivoteo parcial escalado (scaled partial pivoting):** Antes de eliminar, calcula para cada fila un factor de escala, que es su mayor valor absoluto:

$$
s_i = \max_{j} |a_{ij}| \tag{3}
$$

Luego, en el paso $k$, elige la fila que maximiza el cociente $|a_{ik}| / s_i$ y la intercambia con la fila $k$. Es decir, compara el tamaño de cada candidato **relativo a su propia fila**.

- **Qué resuelve que el parcial no:** el parcial compara valores absolutos, y por eso depende de las unidades o la escala de cada ecuación. Si una fila tiene números enormes, el parcial puede elegir como pivote un valor que es grande en términos absolutos pero pequeño comparado con el resto de su fila. Ejemplo clásico (Burden y Faires):

$$
A = \begin{pmatrix} 30 & 591400 \\ 5.291 & -6.130 \end{pmatrix}, \qquad
b = \begin{pmatrix} 591700 \\ 46.78 \end{pmatrix}
$$

El parcial elige $30$ porque $30 > 5.291$. Pero los factores de escala son $s_1 = 591400$ y $s_2 = 6.130$, así que los cocientes son:

$$
\frac{30}{591400} \approx 5.1 \times 10^{-5}, \qquad \frac{5.291}{6.130} \approx 0.863
$$

El escalado elige la fila 2, que es la decisión correcta, y con pocas cifras significativas obtiene la solución exacta $x_1 = 10$, $x_2 = 1$.

- **Qué cuesta de más:** calcular los $s_i$ una sola vez cuesta del orden de $n^2$ operaciones, y en cada paso hay que hacer una división por fila candidata. Es mucho más barato que el pivoteo completo, aunque algo más que el parcial.


#### 6. ¿Por qué el pivoteo parcial es el que usan en la práctica las bibliotecas reales (LAPACK, NumPy, MATLAB) y no el completo, si el completo es más robusto?

## A.3 Teoría y límites

#### 7. Investiga qué es el factor de crecimiento (growth factor) y qué relación tiene con la estabilidad del método. Menciona cuál es su cota teórica con pivoteo parcial y por qué en la práctica casi nunca se alcanza.

#### 8. ¿Existen matrices en las que se puede aplicar eliminación gaussiana sin pivoteo con total seguridad? Nombra al menos dos familias y justifica por qué.

$$
\begin{cases}
0.0001\,x_1 + x_2 = 1 \\
x_1 + x_2 = 2
\end{cases}
$$

$$
A_1 = 
\begin{pmatrix} 
0.0001 & 1 \\ 1 & 1 
\end{pmatrix}, \qquad
b_1 = 
\begin{pmatrix} 1 \\ 2 
\end{pmatrix}
$$