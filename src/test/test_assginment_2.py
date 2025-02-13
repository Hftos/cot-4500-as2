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
    degree = 3  # Use the highest degree polynomial (degree 3)
    approx_value = evaluate_polynomial(x2, fx2, eval_point, degree)
    print(f"\n{approx_value}\n")

if __name__ == "__main__":
    main()