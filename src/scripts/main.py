import csv
import math
import numpy as np


def resolver_referencia(A, B):
    return np.linalg.solve(np.array(A, dtype=float), np.array(B, dtype=float)).tolist()


def error_norma2(x_exacto, x_calc):
    num = math.sqrt(sum((c - e) ** 2 for c, e in zip(x_calc, x_exacto)))
    den = math.sqrt(sum(e ** 2 for e in x_exacto))
    return num / den


def ejecutar(sistema, metodo, funcion, A, B, x_exacto):
    A = [fila[:] for fila in A]
    B = B[:]
    try:
        x_calc, hist = funcion(A, B)
        return {"sistema": sistema,
                "metodo": metodo,
                "solucion": " ".join(f"{v:.6g}" for v in x_calc),
                "error_relativo": error_norma2(x_exacto, x_calc),
                "intercambios": len(hist) - 1}
    except ZeroDivisionError as error:
        return {"sistema": sistema,
                "metodo": metodo,
                "solucion": f"ZeroDivisionError: {error}",
                "error_relativo": "",
                "intercambios": 0}


def guardar_csv(resultados, ruta="./labU3/pivoteo-eliminacion-gaussiana/src/data/resultados.csv"):
    with open(ruta, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(resultados[0].keys()))
        w.writeheader()
        w.writerows(resultados)


def r(x, n=4):
    return float(f"{x:.{n-1}e}")


def error_relativo(x_exacto: list[float], x_calc: list[float]):
    return [f"{r(abs(e - c) / abs(e) * 100)}%" for e, c in zip(x_exacto, x_calc)]


def eliminar_columna(A, B, i):
    n = len(A)
    for j in range(i + 1, n):
        m = r(A[j][i] / A[i][i])
        for k in range(i, n):
            A[j][k] = r(A[j][k] - r(m * A[i][k]))
        B[j] = r(B[j] - r(m * B[i]))


def eliminacion_adelante(A, B):
    for i in range(len(A)):
        eliminar_columna(A, B, i)


def gauss_sin_pivoteo(A, B):
    hist = []
    hist.append(([fila[:] for fila in A], B[:]))
    eliminacion_adelante(A, B)
    return sustitucion_regresiva(A, B), hist


def gauss_pivoteo_parcial(A, B):
    hist = []
    hist.append(([fila[:] for fila in A], B[:]))
    for i in range(len(A)):
        idx = i
        for j in range(i + 1, len(A)):
            if abs(A[j][i]) > abs(A[idx][i]):
                idx = j
        if idx != i:
            A[i], A[idx] = A[idx], A[i]
            B[i], B[idx] = B[idx], B[i]
            hist.append(([fila[:] for fila in A], B[:]))
        eliminar_columna(A, B, i)
    return sustitucion_regresiva(A, B), hist


def gauss_pivoteo_escalado(A, B):
    hist = []
    hist.append(([fila[:] for fila in A], B[:]))
    s = [max(abs(v) for v in fila) for fila in A]   # una sola vez
    for i in range(len(A)):
        idx = i
        for j in range(i + 1, len(A)):
            if abs(A[j][i]) / s[j] > abs(A[idx][i]) / s[idx]:
                idx = j
        if idx != i:
            A[i], A[idx] = A[idx], A[i]
            B[i], B[idx] = B[idx], B[i]
            s[i], s[idx] = s[idx], s[i]
            hist.append(([fila[:] for fila in A], B[:]))
        eliminar_columna(A, B, i)
    return sustitucion_regresiva(A, B), hist


def sustitucion_regresiva(A, B):
    n = len(A)
    x = [0] * n
    for i in range(n - 1, -1, -1):
        suma = 0
        for k in range(i + 1, n):
            suma = suma + r(A[i][k] * x[k])
        x[i] = r((B[i] - suma) / A[i][i])
    return x


def imprimir_matriz(A, B):
    for i in range(len(A)):
        for j in range(len(A)):
            print(A[i][j], end="\t")
        print(B[i])


if __name__ == "__main__":
    A1, B1 = [[0.0001, 1], [1, 1]], [1, 2]
    A2, B2 = [[0, 2, 1], [1, -1, 3], [2, 1, -1]], [5, 6, 1]
    A3, B3 = [[0.001, 2, 1, 3], [1, 0.01, -1, 2], [2, -1, 4, 1], [3, 1, 2, 5]], [19.001, 6.02, 16, 31]

    sistemas = {1: (A1, B1), 2: (A2, B2), 3: (A3, B3)}  
    metodos = {"Sin pivoteo": gauss_sin_pivoteo,
               "Pivoteo parcial": gauss_pivoteo_parcial,
               "Pivoteo escalado": gauss_pivoteo_escalado}

    resultados = []
    for s, (A, B) in sistemas.items():
        x_exacto = resolver_referencia(A, B)
        for nombre, funcion in metodos.items():
            resultados.append(ejecutar(s, nombre, funcion, A, B, x_exacto))

    guardar_csv(resultados)
    print("Listo: resultados.csv")