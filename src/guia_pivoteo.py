def valor_absoluto(A):
    for i in range(len(A)):
        abs(A[i])
    return A


def r5(x):
    return float(f"{x:.5g}")


def error_relativo(x_exacto: list[float], x_calc: list[float]):
    a = []*len(x_calc)
    for i, j in x_exacto, x_calc:
        a.append(str(r5(abs(i- j) / i * 100)) + "%")
    return a


def gauss_sin_pivoteo(A, B):
    for i in range(len(A)):
        if A[i][i] == 0:
            pivote = r5(abs(A[i][i]))
            idx = i
            for j in range(i+1, len(A)):
                if A[j][i] != 0:
                    print("cambio")
                    idx = j
                    break 
            A[i], A[idx] = A[idx], A[i]
            B[i], B[idx] = B[idx], B[i]


def gauss_pivoteo_parcial(A, B):
    for i in range(len(A)):
        if A[i][i] == 0:
            pivote = abs(A[i][i])
            idx = i
            for j in range(i+1, len(A)):
                if abs(A[j][i]) > pivote:
                    print("cambio")
                    pivote = abs(A[j][i]) 
                    idx = j
            A[i], A[idx] = A[idx], A[i]
            B[i], B[idx] = B[idx], B[i]


def gauss_pivoteo_escalado(A, B):
    for i in range(len(A)):
        pivote = r5(abs(A[i][i]) / max(valor_absoluto(A[i])))
        idx = i
        for j in range(i+1, len(A)):
            if abs(A[j][i]) / max(valor_absoluto(A[j])) > pivote:
                print("cambio")
                pivote = r5(abs(A[j][i]) / max(valor_absoluto(A[j])))
                idx = j
        A[i], A[idx] = A[idx], A[i]
        B[i], B[idx] = B[idx], B[i] 
    

def sustitucion_regresiva(A, B):
    n = len(A)
    for i in range(n):
        for j in range(i + 1, n):
                m = r5(A[j][i] / A[i][i])
                for k in range(i, n):
                    A[j][k] = r5(A[j][k] - r5(m * A[i][k]))
                B[j] = r5(B[j] - r5(m * B[i]))


def sustitucion_hacia_atras(A, B):
    n = len(A)
    x = [0] * n
    for i in range(n - 1, -1, -1):
        suma = 0
        for k in range(i + 1, n):
            suma = suma + r5(A[i][k] * x[k])
        x[i] = r5((B[i] - suma) / A[i][i])
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
    gauss_sin_pivoteo(A1, B1)
    sustitucion_regresiva(A1, B1)
    imprimir_matriz(A1, B1)
    x_calc = sustitucion_hacia_atras(A1, B1)
    print(x_calc)
    print(error_relativo([1.00010001, 0.99989999], x_calc))
    print("\n")