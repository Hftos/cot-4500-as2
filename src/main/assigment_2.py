import numpy as np

def neville_interpolation(x,val, w):
    n = len(x)
    neville = [[0 for _ in range(n)] for _ in range(n)]

    for i in range(n):
        neville[i][0]= val[i]
    
    for i in range(1,n):
        for j in range(1, i+1):
            term1 = (w - x[i-j]) * neville[i][j-1]
            term2 = (w - x[i]) * neville[i-1][j-1]
            neville[i][j] = (term1 - term2) / (x[i] - x[i-j])
    
    return neville

def newton_forward_method(x, fx):

    n = len(x)
    h = x[1] - x[0]
    
    diffs = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        diffs[i][0] = fx[i]
    
    for i in range(1, n):
        for j in range(1, i + 1):
            diffs[i][j] = diffs[i][j-1] - diffs[i-1][j-1]
    
    return diffs

def get_coefficients(x, diffs, degree):

    h = x[1] - x[0] 
    return diffs[degree][degree] / (factorial(degree) * h**degree)

def factorial(n):
    if n == 0 or n == 1:
        return 1
    return n*(factorial(n-1))

def evaluate_polynomial(x, fx, eval_point, degree):
    diffs = newton_forward_method(x, fx)
    h = x[1] - x[0]
    
    p = fx[0]
    
    for k in range(1, degree + 1):
        term = diffs[k][k] / (factorial(k) * h**k)
        for i in range(k):
            term *= (eval_point - x[i])
        p += term
    
    return p

def hermite_divided_differences(x, fx, fpx):
    n = len(x)
    m = 2 * n 
    table = [[0.0 for _ in range(5)] for _ in range(m)]
    
    row = 0
    for i in range(n):
        table[row][0]     = x[i]
        table[row][1]     = fx[i]
        table[row + 1][0] = x[i]
        table[row + 1][1] = fx[i]
        row += 2

    for i in range(n):
        table[2*i][2] = fpx[i]
        if i < n - 1:
            table[2*i + 1][2] = (fx[i+1] - fx[i]) / (x[i+1] - x[i])
        else:
            table[2*i + 1][2] = 0.0
    
    
    for col in range(3, 5):
        for i in range(m - (col - 1)):
            denom = table[i + (col - 1)][0] - table[i][0]
            if abs(denom) < 1e-15:
                table[i][col] = 0.0
            else:
                table[i][col] = (table[i+1][col-1] - table[i][col-1]) / denom

    return table

def build_natural_spline_system(x, y):
    h0 = x[1] - x[0]
    h1 = x[2] - x[1]
    h2 = x[3] - x[2]

    
    A = np.array([
        [1,           0,          0, 0],
        [h0, 2*(h0 + h1),        h1, 0],
        [0,           h1, 2*(h1 + h2), h2],
        [0,           0,          0, 1]
    ], dtype=float)

    diff_1_0 = (y[1] - y[0]) / h0 
    diff_2_1 = (y[2] - y[1]) / h1  
    diff_3_2 = (y[3] - y[2]) / h2 

    B = np.array([
        0,
        3 * (diff_2_1 - diff_1_0),
        3 * (diff_3_2 - diff_2_1), 
        0
    ], dtype=float)

    return A, B


def main():
    # Question 1 input
    x1 = [3.6, 3.8, 3.9]
    val1 = [1.675, 1.436, 1.318]
    w1 = 3.7
    
    result1 = neville_interpolation(x1, val1, w1)
    interpolated_value1 = result1[2][2]
    print(f"{interpolated_value1}\n")

    # Question 2 input
    x2 =[7.2, 7.4, 7.5, 7.6]
    fx2 = [23.5492, 25.3913, 26.8224, 27.4589]

    diffs = newton_forward_method(x2, fx2)

    for degree in range(1, 4):
        coeff = get_coefficients(x2, diffs, degree)
        print(f"{coeff}")

    # Question 3
    eval_point = 7.3
    degree = 3
    approx_value = evaluate_polynomial(x2, fx2, eval_point, degree)
    print(f"\n{approx_value}\n")

    # Question 4 input
    x4 = [3.6, 3.8, 3.9]
    fx4 = [1.675, 1.436, 1.318]
    fpx4 = [-1.195, -1.188, -1.182]
    
    table = hermite_divided_differences(x4, fx4, fpx4)

    for row in table:
        print("[ ", end="")
        for j in range(5):
            print(f"{row[j]:12.6e}", end=" ")
        print("]")
    
    print("\n")
    # Question 5 input
    x5 = [2.0, 5.0, 8.0, 10.0]
    y5 = [3.0, 5.0, 7.0,  9.0]

    A, B = build_natural_spline_system(x5, y5)
    M = np.linalg.solve(A, B)

    for row in A:
        print(row)
    print(B)
    print(M)

if __name__ == "__main__":
    main()