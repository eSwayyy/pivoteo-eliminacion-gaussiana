def valor_absoluto(A):
    for i in range(len(A)):
        abs(A[i])
    return A


def r(x, n=4):
    return float(f"{x:.{n-1}e}")


def error_relativo(x_exacto: list[float], x_calc: list[float]):
    return [f"{r(abs(e - c) / abs(e) * 100)}%" for e, c in zip(x_exacto, x_calc)]


def gauss_sin_pivoteo(A, B):
    hist = []
    hist.append(([fila[:] for fila in A], B[:]))
    eliminacion_adelante(A, B)
    return sustitucion_regresiva(A, B), hist



def gauss_pivoteo_parcial(A, B):
    hist = []
    hist.append(([fila[:] for fila in A], B[:]))
    for i in range(len(A)):
        if A[i][i] == 0:
            pivote = abs(A[i][i])
            idx = i
            for j in range(i+1, len(A)):
                if abs(A[i][j]) > abs(A[idx][j]):
                    idx = i
            if idx != j:
                A[j], A[idx] = A[idx], A[j]
                B[j], B[idx] = B[idx], B[j]
                hist.append(([f[:] for f in A], B[:]))
        eliminar_columna(A, B, j)
    return sustitucion_regresiva(A, B), hist
    


def gauss_pivoteo_escalado(A, B):
    hist = []
    hist.append(([fila[:] for fila in A], B[:]))
    s = [max(abs(v) for v in fila) for fila in A]
    for i in range(len(A)):
        pivote = r(abs(A[i][i]) / max(valor_absoluto(A[i])))
        idx = i
        for j in range(i+1, len(A)):
            if abs(A[j][i]) / s[i] > abs(A[idx][j] / s[idx]):
                idx = i
        if idx != j:
            A[j], A[idx] = A[idx], A[j]
            B[j], B[idx] = B[idx], B[j]
            s[j], s[idx] = s[idx], s[j]
            hist.append(([fila[:] for fila in A], B[:]))
    eliminacion_adelante(A, B)
    return sustitucion_regresiva(A, B), hist
    

def eliminar_columna(A, B, k):
    n = len(A)
    for j in range(k + 1, n):
        m = r(A[j][k] / A[k][k])
        for c in range(k, n):
            A[j][c] = r(A[j][c] - r(m * A[k][c]))
        B[j] = r(B[j] - r(m * B[k]))


def eliminacion_adelante(A, B):
    for k in range(len(A)):
        eliminar_columna(A, B, k)


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
    A1 = [
        [0.0001, 1],
        [1, 1]
    ]
    B1 = [1, 2]

    A2 = [
        [0.003000, 59.14],
        [5.291,  -6.130]
    ]
    B2 = [59.17, 46.78]
    x_calc, historial = gauss_pivoteo_escalado(A2, B2)
    for i in historial:
        print(i)
    print(x_calc)
    print(error_relativo([10, 1], x_calc))
    print("\n")