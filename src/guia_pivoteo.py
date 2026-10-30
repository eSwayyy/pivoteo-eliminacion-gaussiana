# Eliminacion Gaussiana

def valor_absoluto(A):
    for i in range(len(A)):
        abs(A[i])
    return A


def gauss_sin_pivoteo(A, B):
    for i in range(len(A)):
        if A[i][i] == 0:
            pivote = abs(A[i][i])
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
        pivote = abs(A[i][i]) / max(valor_absoluto(A[i]))
        idx = i
        for j in range(i+1, len(A)):
            if abs(A[j][i]) / max(valor_absoluto(A[j])) > pivote:
                print("cambio")
                pivote = abs(A[j][i]) / max(valor_absoluto(A[j]))
                idx = j
        A[i], A[idx] = A[idx], A[i]
        B[i], B[idx] = B[idx], B[i] 
    

def sustitucion_regresiva(A, B):
    n = len(A)
    for i in range(n):
        for j in range(i + 1, n):
                m = A[j][i] / A[i][i]
                for k in range(i, n):
                    A[j][k] = A[j][k] - m * A[i][k]
                B[j] = B[j] - m * B[i]


def sustitucion_hacia_atras(A, B):
    n = len(A)
    x = [0] * n
    for i in range(n - 1, -1, -1):
        suma = 0
        for k in range(i + 1, n):
            suma = suma + A[i][k] * x[k]
        x[i] = (B[i] - suma) / A[i][i]
    return x


def imprimir_matriz(A, B):
    for i in range(len(A)):
        for j in range(len(A)):
            print(A[i][j], end="\t")
        print(B[i])


if __name__ == "__main__":
    A1 = [
        [0.00000001, 1],
        [1, 1]
    ]
    B1 = [1, 2]
    # imprimir_matriz(A1, B1)
    gauss_pivoteo_escalado(A1, B1)
    sustitucion_regresiva(A1, B1)
    imprimir_matriz(A1, B1)
    print(sustitucion_hacia_atras(A1, B1))
    print("\n")

    A2 = [
        [1, 1, 1],
        [2, 2, 5],
        [4, 6, 8]
    ]
    B2 = [3, 3, 10]
    # pivoteo_trivial(A2, B2)
    # print()
    # imprimir_matriz(A2, B2)