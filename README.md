# cot-4500-as2
#### Description
This repository contains the solutions to Assignment 2 for COT-4500.
It implements the following interpolation techniques:
- Neville's Method
- Newton's Forward Method
- Polynomial Approximation
- Hermite Divided Differences
- Cubic Spline Interpolation

#### Requirements
This program uses Python 3.x and requires:
- numpy>=1.18.0
- pytest>=7.0.0
Install the dependencies with the following:
```sh
pip install -r requirements.txt
```
#### How to Run
1. Clone the repository
```sh
git clone git@github.com:Hftos/cot-4500-as2.git
```
2. Navigate to cot-4500-as2/src/main
3. Run
```sh
python assigment_2.py
```

#### How to test
The repository can be tested using pytest.
1. Install pytest
pytest is in the requirements, but if that failed then use the following
```sh
pip install pytest
```
2. From the root directory, run the following
```sh
pytest
```