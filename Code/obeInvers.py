# from obeInvers.py

import numpy as np

# functions
def cheksDets(M):
    D = np.linalg.det(M)
    return D

def printMatrix(M):
    for row in M:
        formatted = ["{:.2f}".format(x) for x in row]
        print("   [" + ", ".join(formatted) + "]")
    print()

def gaussStep(matrix4, B):
    n = len(matrix4)
    A = [row[:] for row in matrix4]   # copy matrix
    step = 1

    print("Process of Gauss-Jordan Elimination:\n")
    
    for i in range(n):
        pivot = A[i][i]

        # Jika pivot = 0 maka cari baris lain
        if pivot == 0:
            for k in range(i+1, n):
                if A[k][i] != 0:
                    print(f"Iteration-{step}: Swapping R{i+1} <-> R{k+1}")
                    A[i], A[k] = A[k], A[i]
                    B[i], B[k] = B[k], B[i]
                    print("Matrix A becomes:")
                    printMatrix(A)
                    print("Matrix inverse of A becomes:")
                    printMatrix(B)
                    step += 1
                    pivot = A[i][i]
                    break

        # Jika tetap 0 maka tidak bisa ditentukan
        if pivot == 0:
            continue

        # Normalisasi pivot: yaitu jadikan 1 utama
        print(f"Iteration-{step}: Normalizing R{i+1} (divide by {pivot:.2f})")
        for c in range(len(A[0])):
            A[i][c] /= pivot
            B[i][c] /= pivot
        print("Matrix A becomes:")
        printMatrix(A)
        print("Matrix inverse of A becomes:")
        printMatrix(B)
        step += 1

        # Eliminasi ke semua baris lain (Gauss-Jordan)
        for r in range(n):
            if r == i:
                continue
            if A[r][i] == 0:
                continue

            factor = A[r][i]
            print(f"Iteration-{step}: R{r+1} = R{r+1} - ({factor:.2f}) * R{i+1}")

            for c in range(len(A[0])):
                A[r][c] -= factor * A[i][c]
                B[r][c] -= factor * B[i][c]
            print("Matrix A becomes:")
            printMatrix(A)
            print("Matrix inverse of A becomes:")
            printMatrix(B)
            step += 1

    return A, B

def showStep(matrix4):
    print("Matrix A is:")
    printMatrix(matrix4)
    row = len(matrix4)  
    col = len(matrix4)

    if row != col:                                  # because determinant doesn't exist
        print("No inverse exist for matrix A")
        return

    checkDet = cheksDets(matrix4)
    if checkDet == 0:                               # because determinant should not be zero
        print("No inverse exist for matrix A")
        return 
    else:
        B = []
        for i in range (row):
            B.append([])                            # simply A = [[1 if i == j else 0 for j in range(col)] for i in range(row)]
            for j in range (col):
                if i == j:
                    B[i].append(1)
                else:
                    B[i].append(0)
    
    matrix4, B = gaussStep(matrix4, B)

    print("So matrix A:")
    printMatrix(matrix4)
    print("Have an inverse matrix:")
    printMatrix(B)

# function calls

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

# undefine determinant handling

D = [
    [1, 3, 2],
    [4, 6, 5]
]

print("=== Matrix Inverse of A using Gauss-Jordan ERO ===")
showStep(C)