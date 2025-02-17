import pytest
import numpy as np

from src.main.assignment_2 import (
    neville_interpolation,
    newton_forward_method,
    get_coefficients,
    evaluate_polynomial,
    hermite_divided_differences,
    build_natural_spline_system
)

def test_neville_interpolation():

    x = [3.6, 3.8, 3.9]
    val = [1.675, 1.436, 1.318]
    w = 3.7
    neville_table = neville_interpolation(x, val, w)
    approx_value = neville_table[-1][-1]
    expected_value = 1.5550
    tolerance = 1e-2
    assert abs(approx_value - expected_value) < tolerance

def test_newton_forward_method():

    x = [7.2, 7.4, 7.5, 7.6]
    fx = [23.5492, 25.3913, 26.8224, 27.4589]
    diffs = newton_forward_method(x, fx)
    expected = 25.3913 - 23.5492
    assert abs(diffs[1][1] - expected) < 1e-7

def test_get_coefficients():
    
    x = [7.2, 7.4, 7.5, 7.6]
    fx = [23.5492, 25.3913, 26.8224, 27.4589]
    diffs = newton_forward_method(x, fx)
    coeff_deg1 = get_coefficients(x, diffs, 1)
    assert coeff_deg1 != 0

def test_evaluate_polynomial():

    x = [7.2, 7.4, 7.5, 7.6]
    fx = [23.5492, 25.3913, 26.8224, 27.4589]
    val = evaluate_polynomial(x, fx, 7.3, 3)
    assert 24.0 < val < 26.0

def test_hermite_divided_differences():
    
    x = [3.6, 3.8, 3.9]
    fx = [1.675, 1.436, 1.318]
    fpx = [-1.195, -1.188, -1.182]
    table = hermite_divided_differences(x, fx, fpx)
    assert abs(table[0][0] - 3.6) < 1e-12
    assert abs(table[5][1] - 1.318) < 1e-12

def test_build_natural_spline_system():
    
    x = [2.0, 5.0, 8.0, 10.0]
    y = [3.0, 5.0, 7.0, 9.0]
    A, B = build_natural_spline_system(x, y)

    assert np.allclose(A[1], [3.0, 12.0, 3.0, 0.0])
    assert np.allclose(A[2], [0.0, 3.0, 10.0, 2.0])
    assert np.allclose(B, [0.0, 0.0, 1.0, 0.0])

    M = np.linalg.solve(A, B)

    assert abs(M[0]) < 1e-12
    assert abs(M[3]) < 1e-12
    assert abs(M[1] + 0.027) < 0.001
    assert abs(M[2] - 0.108) < 0.001
