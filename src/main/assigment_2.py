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

def main():
    # Question 1 input
    x1 = [3.6, 3.8, 3.9]
    val1 = [1.675, 1.436, 1.318]
    w1 = 3.7
    
    result1 = neville_interpolation(x1, val1, w1)
    interpolated_value1 = result1[2][2]
    print(f"{interpolated_value1}\n")

if __name__ == "__main__":
    main()