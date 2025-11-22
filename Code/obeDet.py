import numpy as np

swap = 0

def calcDet(A):
    D = np.linalg.det(A)
    return D

def printMatrix(M):
    for row in M:
        print("   ", row)
    print()

def gaussStep(matrix2):
    global swap   

    n = len(matrix2)
    A = [row[:] for row in matrix2]   # copy matrix
    step = 1

    print("Proses eliminasi OBE:\n")

    for i in range(n - 1):
        pivot = A[i][i]

        # pivot nol maka tukar baris
        if pivot == 0:
            for k in range(i + 1, n):
                if A[k][i] != 0:
                    print(f"Iterasi {step}: Tukar R{i+1} <-> R{k+1}")
                    A[i], A[k] = A[k], A[i]
                    printMatrix(A)
                    step += 1
                    pivot = A[i][i]
                    swap += 1       
                    break

        # eliminasi
        for r in range(i + 1, n):
            if A[r][i] == 0:
                continue

            factor = A[r][i] / pivot
            print(f"Iterasi {step}: R{r+1} = R{r+1} - ({factor:.2f}) * R{i+1}")

            for c in range(n):
                A[r][c] = A[r][c] - factor * A[i][c]

            printMatrix(A)
            step += 1

    print("Akhir proses OBE (bentuk segitiga atas):")
    printMatrix(A)

    return A


def showStep2(matrix2):
    row = len(matrix2)
    col = len(matrix2[0])

    print("Matriks A:")
    printMatrix(matrix2)

    if row != col:
        print(f"So sorry, there's no valid determinant of A {row} x {col}")
        print("Please input determinant with order of n x n")
        print()
    else:
        print(f"Operasi OBE untuk matriks A {row}x{col}:\n")
        matrix2 = gaussStep(matrix2)
        det = 1
        diag = []

        for i in range(row):
            diag.append(f"({matrix2[i][i]})")
            det *= matrix2[i][i]

        swapSign = pow(-1, swap)
        print("Determinant is:", " * ".join(diag))
        print(f"Determinant result: {swapSign * det:.2f}, because we swap {swap} row(s)")
        print()

# test sample
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
    [2, 3],
    [4, 5]
]

# undefine determinant handling

D = [
    [1, 3, 2],
    [4, 6, 5]
]

print("\n=== Determinant with Elementary Row Operation(s) ===\n")
showStep2(D)
