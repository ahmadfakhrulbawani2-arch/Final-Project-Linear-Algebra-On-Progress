# file dari obeSpl.py

import numpy as np

def calcDet(A):
    D = np.linalg.det(A)
    return D

def printMatrix(M):
    for row in M:
        formatted = ["{:.2f}".format(x) for x in row]
        print("   [" + ", ".join(formatted) + "]")
    print()

def gaussStep(matrix3):
    n = len(matrix3)
    A = [row[:] for row in matrix3]   # salin matriks
    step = 1

    print("Proses Eliminasi Gauss-Jordan:\n")
    
    for i in range(n):
        pivot = A[i][i]

        # Jika pivot = 0, cari baris lain untuk swap
        if pivot == 0:
            for k in range(i+1, n):
                if A[k][i] != 0:
                    print(f"Iterasi-{step}: Tukar R{i+1} <-> R{k+1}")
                    A[i], A[k] = A[k], A[i]
                    printMatrix(A)
                    step += 1
                    pivot = A[i][i]
                    break

        # Jika tetap 0 maka lanjut (tidak bisa normalisasi)
        if pivot == 0:
            continue

        # Normalisasi pivot: jadikan 1
        print(f"Iterasi-{step}: Normalisasi R{i+1} (bagi {pivot:.2f})")
        for c in range(len(A[0])):
            A[i][c] /= pivot
        printMatrix(A)
        step += 1

        # Eliminasi ke semua baris lain (Gauss-Jordan)
        for r in range(n):
            if r == i:
                continue
            if A[r][i] == 0:
                continue

            factor = A[r][i]
            print(f"Iterasi-{step}: R{r+1} = R{r+1} - ({factor:.2f}) * R{i+1}")

            for c in range(len(A[0])):
                A[r][c] -= factor * A[i][c]

            printMatrix(A)
            step += 1

    print("Matriks diagonal akhir:")
    printMatrix(A)

    return A

def solveFromDiagonal(matrix3):
    row = len(matrix3)
    col = len(matrix3[0])
    vars = col - 1

    solutions = []

    for i in range(vars):
        # Jika diagonal = 0
        if abs(matrix3[i][i]) < 1e-9:
            if abs(matrix3[i][col-1]) < 1e-9:
                return None, "infinite"
            else:
                return None, "nosolution"

        solutions.append(matrix3[i][col-1])

    return solutions, "ok"

def showSol(matrix3):
    row = len(matrix3)
    col = len(matrix3[0])

    if row < col - 1:
        print("Jumlah persamaan tidak cukup untuk solusi unik.\n")
        return
    
    print("Matriksnya adalah:\n")
    printMatrix(matrix3)
    print("Solusi:\n")

    finalM = gaussStep(matrix3)
    sol, status = solveFromDiagonal(finalM)

    print("Maka solusinya adalah:\n")

    if status == "nosolution":
        print("Tidak ada solusi (persamaan tidak konsisten)")
    elif status == "infinite":
        print("Solusi tak hingga (sistem kurang determinan)")
    else:
        for i, val in enumerate(sol):
            print(f"x{i+1} = {val:.4f}")


# Pemanggilan fungsi

A1 = [
    [2, 1, -1, 8],
    [-3, -1, 2, -11],
    [-2, 1, 2, -3],
]

A = [
    [2, 3, 1, 4],
    [1, 0, -1, 2],
    [3, 1, 2, 0],
    [2, 4, 1, 3]
]

B = [
    [2,  3,  1],
    [4, -1,  2],
    [5,  0,  3]
]

C = [
    [2, 1],
    [5, 3]
]

# Matriks tidak persegi / tidak lengkap
D = [
    [1, 3, 2],
    [4, 6, 5]
]

print("\n=== Penyelesaian SPL dengan Gauss-Jordan OBE ===\n")
showSol(A1)
