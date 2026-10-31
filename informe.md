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

En el sistema del problema motivador, el multiplicador del primer paso es:

$$
m_{21} = \frac{1}{0.0001} = 10^4 \tag{2}
$$

La cancelación ocurre en la actualización del elemento $a_{22}$, que en el código corresponde a la línea de `eliminar_columna`:

```python
A[j][k] = r(A[j][k] - r(m * A[i][k]))
```

Con los valores del sistema:

$$
a_{22} \leftarrow 1 - 10^4 \cdot 1 = -9999 \tag{3}
$$

Con 3 cifras significativas, $-9999$ se redondea a $-1.00 \times 10^4$ y el $1$ original desaparece. Lo mismo pasa con el lado derecho: $b_2 \leftarrow 2 - 10^4 \cdot 1 = -9998 \approx -1.00 \times 10^4$. Entonces $x_2 = 1.00$, y al sustituir hacia atrás $x_1 = (1 - 1.00)/0.0001 = 0$, cuando lo correcto es $1.0001$. La función `r` solo simula la precisión limitada; la información se destruye en esa resta, porque $10^4$ es tan grande que absorbe al $1$.


## A.2 Las estrategias

#### 4. Describe el pivoteo parcial (partial pivoting): qué busca, dónde busca y qué intercambia.

El pivoteo parcial busca, en cada paso $i$, el elemento de mayor valor absoluto en la columna $i$ entre las filas $i, i+1, \dots, n$ (de la diagonal hacia abajo), y lo intercambia con la fila $i$, junto con el elemento correspondiente de $b$. Se aplica **siempre**, no solo cuando el pivote vale 0. Resuelve dos problemas: evita la división por cero y, sobre todo, evita pivotes pequeños. Como el pivote es el mayor de la columna, todos los multiplicadores cumplen $|m_{ij}| \leq 1$, lo que impide que los errores de redondeo se amplifiquen.

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
|a_{pq}| = \max_{i,j \geq k} |a_{ij}| \tag{4}
$$

Como intercambiar columnas equivale a intercambiar el orden de las incógnitas, hay que llevar un registro de esa permutación para reordenar $x$ al final.

**Qué resuelve que el parcial no:** el parcial solo compara dentro de una columna, así que puede quedar con un pivote relativamente pequeño si el mejor candidato está en otra columna. El completo garantiza el pivote más grande posible en cada paso, lo que da la mejor estabilidad numérica y permite detectar el rango de la matriz.

**Qué cuesta de más:** buscar el máximo en la submatriz de tamaño $(n-k+1)^2$ en cada paso suma del orden de $n^3/3$ comparaciones en total, es decir, del mismo orden que las operaciones aritméticas de la eliminación. Además hay que llevar la permutación de columnas.

**Pivoteo parcial escalado (scaled partial pivoting):** Antes de eliminar, calcula para cada fila un factor de escala, que es su mayor valor absoluto:

$$
s_i = \max_{j} |a_{ij}| \tag{5}
$$

Luego, en el paso $k$, elige la fila que maximiza el cociente $|a_{ik}| / s_i$ y la intercambia con la fila $k$. Es decir, compara el tamaño de cada candidato **relativo a su propia fila**.

**Qué resuelve que el parcial no:** el parcial compara valores absolutos, y por eso depende de las unidades o la escala de cada ecuación. Si una fila tiene números enormes, el parcial puede elegir como pivote un valor que es grande en términos absolutos pero pequeño comparado con el resto de su fila. Ejemplo clásico (Burden y Faires):

$$
A = \begin{pmatrix} 30 & 591400 \\ 5.291 & -6.130 \end{pmatrix}, \qquad
b = \begin{pmatrix} 591700 \\ 46.78 \end{pmatrix}
$$

El parcial elige $30$ porque $30 > 5.291$. Pero los factores de escala son $s_1 = 591400$ y $s_2 = 6.130$, así que los cocientes son:

$$
\frac{30}{591400} \approx 5.1 \times 10^{-5}, \qquad \frac{5.291}{6.130} \approx 0.863
$$

El escalado elige la fila 2, que es la decisión correcta, y con pocas cifras significativas obtiene la solución exacta $x_1 = 10$, $x_2 = 1$.

**Qué cuesta de más:** calcular los $s_i$ una sola vez cuesta del orden de $n^2$ operaciones, y en cada paso hay que hacer una división por fila candidata. Es mucho más barato que el pivoteo completo, aunque algo más que el parcial.


#### 6. ¿Por qué el pivoteo parcial es el que usan en la práctica las bibliotecas reales (LAPACK, NumPy, MATLAB) y no el completo, si el completo es más robusto?

Las bibliotecas reales usan pivoteo parcial porque ofrece casi la misma estabilidad que el completo a un costo mucho menor. La diferencia está en la búsqueda del pivote: el parcial revisa una columna por paso, mientras que el completo revisa toda la submatriz restante.

$$
\underbrace{\sum_{i=1}^{n} (n-i+1)}_{\text{parcial}} \approx \frac{n^2}{2}
\qquad \text{frente a} \qquad
\underbrace{\sum_{i=1}^{n} (n-i+1)^2}_{\text{completo}} \approx \frac{n^3}{3} \tag{6}
$$

## A.3 Teoría y límites

#### 7. Investiga qué es el factor de crecimiento (growth factor) y qué relación tiene con la estabilidad del método. Menciona cuál es su cota teórica con pivoteo parcial y por qué en la práctica casi nunca se alcanza.

El **factor de crecimiento** mide cuánto pueden crecer los elementos de la matriz durante la eliminación, comparado con el tamaño original:

$$
\rho = \frac{\max_{i,j,k} |a_{ij}^{(k)}|}{\max_{i,j} |a_{ij}|} \tag{7}
$$

donde $a_{ij}^{(k)}$ es el elemento después del paso $k$.

**Relación con la estabilidad:** el error de retroceso de la eliminación gaussiana es proporcional a $\rho$. Si $\rho$ es pequeño, el método es estable; si es grande, los elementos intermedios crecen y arrastran errores de redondeo que contaminan la solución [3, cap. 3], [2, lec. 22]. Sin pivoteo, $\rho$ puede ser arbitrariamente grande: en el Sistema 1, el multiplicador $10^4$ hace crecer $a_{22}$ desde $1$ hasta $-9999$.

**Cota con pivoteo parcial:** como $|m_{ij}| \leq 1$, cada paso puede, como máximo, duplicar el tamaño de los elementos, por lo que:

$$
\rho \leq 2^{n-1} \tag{8}
$$

**Por qué casi nunca se alcanza:** para llegar a $2^{n-1}$ se necesita que en cada paso todos los signos y magnitudes se alineen de forma perfecta, una estructura que las matrices que aparecen en la práctica casi nunca tienen. En matrices aleatorias y en las de aplicaciones reales, $\rho$ suele ser pequeño (del orden de $n^{2/3}$ o menos en experimentos) [2, lec. 22], [4, cap. 9]. Por eso el pivoteo parcial se considera estable en la práctica, aunque su garantía teórica sea mala.



#### 8. ¿Existen matrices en las que se puede aplicar eliminación gaussiana sin pivoteo con total seguridad? Nombra al menos dos familias y justifica por qué.

Existen dos familias:

- **Matrices Estrictamente Diagonalmente Dominantes por Filas** \
Una matriz la cual el valor absoluto del elemento en la diagonal es estrictamente mayor que la suma de los valores absolutos de todos los demás elementos de esa misma fila:

$$
|{a_{ii}}| > \sum_{j\neq i}|a_{ij}|,\ \forall i =1,...n
$$
- **Matrices Definidas Positivas** \
Una matriz simétrica (o hermítica en el caso complejo) $A$ es definida positiva si para cualquier vector no nulo $x$, se cumple que $x^T A x > 0$.

## Matriz 4x4

- Empezamos por una matriz "sana", con enteros entre $-1$ y $5$ elegidos a mano. Para calcular $B_3$ elegimos primero un $x_{exacto}$ simple y calculamos el producto matricial $B_3 = A_3\,x$. Por ejemplo, la primera fila da $1\cdot1 + 2\cdot2 + 1\cdot3 + 3\cdot4 = 20$. $B_3$ no se redondea, para que la solución exacta sea exactamente $x$.

$$
A_3 = \begin{pmatrix}
1 & 2 & 1 & 3 \\
1 & 1 & -1 & 2 \\
2 & -1 & 4 & 1 \\
3 & 1 & 2 & 5
\end{pmatrix}, \qquad
x = \begin{pmatrix} 1 \\ 2 \\ 3 \\ 4 \end{pmatrix}, \qquad
B_3 = \begin{pmatrix} 20 \\ 8 \\ 16 \\ 31 \end{pmatrix}
$$

- Luego hacemos pequeños dos elementos de la diagonal ($a_{11}$ y $a_{22}$). Esto funciona porque un pivote pequeño genera multiplicadores enormes $(1000, 2000, 3000)$, y al restar $a_{jk} - m\,a_{ik}$ el término grande absorbe al pequeño (el mismo mecanismo de la pregunta 3). El resto de la matriz no cambia, y $B_3$ se recalcula con el mismo $x$.

$$
A_3 = \begin{pmatrix}
0.001 & 2 & 1 & 3 \\
1 & 0.01 & -1 & 2 \\
2 & -1 & 4 & 1 \\
3 & 1 & 2 & 5
\end{pmatrix}, \qquad
B_3 = \begin{pmatrix} 19.001 \\ 6.02 \\ 16 \\ 31 \end{pmatrix}
$$

- Los valores $0.001$ y $0.01$ y las cifras de $r$ se eligieron probando: con 4 cifras significativas el sin pivoteo da un resultado equivocado sin lanzar excepción (con 3 cifras lanza `ZeroDivisionError`).

| Matriz | Sin pivoteo | Parcial | Escalado |
|---|---|---|---|
| Sana | $0$ | $1.4\times10^{-3}$ | $6.6\times10^{-4}$ |
| Con pivotes pequeños | $1.37$ | $8.7\times10^{-4}$ | $8.7\times10^{-4}$ |

En la matriz modificada, parcial y escalado empatan porque las filas tienen magnitudes parecidas y eligen los mismos pivotes.

## Analisis

### ¿En qué casos el pivoteo parcial y el escalado dieron el mismo resultado? ¿Por qué?

En los sistemas 1, 2 y 3 dieron la misma solución, el mismo error y el mismo número de intercambios (1, 2 y 2). Esto ocurre porque en esos sistemas las filas tienen magnitudes parecidas: los factores de escala son (1, 1), (2, 3, 2) y (3, 2, 4, 5). Al dividir cada candidato entre un factor similar, el orden de los candidatos no cambia, así que ambos métodos eligen el mismo pivote.

<!-- Con la aritmética de 4 cifras, en el sistema 1 los tres métodos coinciden (error 7,07·10⁻⁵), porque el resultado exacto 1,0001 se redondea a 1. Sin pivoteo, el sistema 2 falla (división por cero, ya que a₁₁ = 0) y en el sistema 3 el error relativo es 1,37, frente a 8,7·10⁻⁴ con pivoteo. -->

### ¿Hay algún sistema donde el pivoteo escalado ganó claramente?

En los sistemas dados, no. Se construyó uno donde sí ocurre, mezclando filas de magnitudes muy distintas:

    30,00·x + 591400·y = 591700
    5,291·x − 6,130·y = 46,78        (solución exacta: x = 10, y = 1)

Con aritmética de 4 cifras:

| Método | Solución | Error relativo | Intercambios |
|---|---|---|---|
| Sin pivoteo | x = −10, y = 1,001 | 1,99 | 0 |
| Pivoteo parcial | x = −10, y = 1,001 | 1,99 | 0 |
| Pivoteo escalado | x = 10, y = 1 | 2·10⁻¹³ | 1 |

El pivoteo parcial compara solo 30,00 y 5,291, y elige el 30,00. Pero ese valor es diminuto en comparación con el 591400 de su fila: los cocientes del escalado son 30/591400 ≈ 5·10⁻⁵ contra 5,291/6,130 ≈ 0,86, y por eso intercambia las filas y obtiene la solución correcta.

### Tiempo de ejecución en una matriz aleatoria de 200×200

| Método | Tiempo (s) | Diferencia vs. sin pivoteo |
|---|---|---|
| Sin pivoteo | 9.003 s | referencia |
| Pivoteo parcial | 8.259 s | −8.3 % |
| Pivoteo escalado | 8.542 s | −5.1 % |

Los tiempos se midieron en una sola ejecución por método. Las diferencias son negativas, es decir, los métodos con pivoteo salieron más rápidos que el que no busca pivote, lo cual no puede ser un efecto real, porque hacen el mismo trabajo más la búsqueda. Concluimos que la diferencia medida está dentro del ruido de la medición (orden de ejecución, estado de la CPU, otros procesos), y que el costo del pivoteo es demasiado pequeño para distinguirlo con esta técnica. Entre los dos métodos con pivoteo, el escalado tardó un 3.4 % más que el parcial.

Este resultado coincide con la estimación teórica: la eliminación crece como $n^3/3$ (unas 2.7 millones de operaciones para $n = 200$), mientras que la búsqueda del pivote parcial crece como $n^2/2$ (unas 20 mil comparaciones), es decir, cerca de un 0.75 % del trabajo de la eliminación. Además, en nuestra implementación cada operación aritmética pasa por la función `r`, que formatea el número como texto y es mucho más cara que una resta normal, lo que hace aún más pequeña la fracción que representa la búsqueda. El escalado cuesta algo más que el parcial porque calcula los factores de escala y hace una división por cada candidato.

**¿Justifica el costo el beneficio?** Sí. Un costo menor al 1 % del tiempo (no medible en nuestro experimento) evita errores de orden de magnitud: sin pivoteo, el Sistema 2 falla con `ZeroDivisionError` y el Sistema 3 tiene un error relativo de 137 %. El pivoteo parcial suele bastar; el escalado se justifica cuando las filas del sistema tienen escalas muy distintas, como en el ejemplo de Burden y Faires.