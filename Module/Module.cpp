#ifndef MATH_OPERATIONS_HPP
#define MATH_OPERATIONS_HPP

#include <cmath>
#include <iostream>
#include <vector>
namespace MathOperations {
    double sin(double x) {
        return std::sin(x);
    }

    double cos(double x) {
        return std::cos(x);
    }

    double tan(double x) {
        return std::tan(x);
    }

    double cot(double x) {
        return 1.0 / std::tan(x);
    }

    double power(double base, double exponent) {
        return std::pow(base, exponent);
    }
}

#endif // MATH_OPERATIONS_HPP
#ifndef MATRIX_HPP
#define MATRIX_HPP

class Matrix {
private:
    std::vector<std::vector<double>> data;
    size_t rows;
    size_t cols;

public:
    Matrix() : rows(0), cols(0) {}

    Matrix(size_t rows, size_t cols) : rows(rows), cols(cols) {
        data.resize(rows, std::vector<double>(cols, 0.0));
    }

    double& operator()(size_t row, size_t col) {
        return data[row][col];
    }

    const double& operator()(size_t row, size_t col) const {
        return data[row][col];
    }

    size_t getRows() const {
        return rows;
    }

    size_t getCols() const {
        return cols;
    }

    void print() const {
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                std::cout << data[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Matrix dimensions must match for addition");
        }

        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result(i, j) = data[i][j] + other(i, j);
            }
        }
        return result;
    }

    Matrix operator-(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Matrix dimensions must match for subtraction");
        }

        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result(i, j) = data[i][j] - other(i, j);
            }
        }
        return result;
    }

    Matrix operator*(double scalar) const {
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result(i, j) = data[i][j] * scalar;
            }
        }
        return result;
    }

    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows) {
            throw std::invalid_argument("Number of columns in first matrix must match number of rows in second matrix for multiplication");
        }

        Matrix result(rows, other.cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < other.cols; ++j) {
                for (size_t k = 0; k < cols; ++k) {
                    result(i, j) += data[i][k] * other(k, j);
                }
            }
        }
        return result;
    }

    Matrix transpose() const {
        Matrix result(cols, rows);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result(j, i) = data[i][j];
            }
        }
        return result;
    }
};

#endif // MATRIX_HPP
int main() {
    double x = 3.14159; // Пі
    std::cout << "sin(" << x << ") = " << MathOperations::sin(x) << std::endl;
    std::cout << "cos(" << x << ") = " << MathOperations::cos(x) << std::endl;
    std::cout << "tan(" << x << ") = " << MathOperations::tan(x) << std::endl;
    std::cout << "cot(" << x << ") = " << MathOperations::cot(x) << std::endl;
    std::cout << "2^3 = " << MathOperations::power(2, 3) << std::endl;

    Matrix A(2, 2);
    A(0, 0) = 1; A(0, 1) = 2;
    A(1, 0) = 3; A(1, 1) = 4;

    Matrix B(2, 2);
    B(0, 0) = 5; B(0, 1) = 6;
    B(1, 0) = 7; B(1, 1) = 8;

    std::cout << "Matrix A:\n";
    A.print();

    std::cout << "\nMatrix B:\n";
    B.print();

    Matrix C = A + B;
    std::cout << "\nMatrix A + B:\n";
    C.print();

    Matrix D = A * 2;
    std::cout << "\nMatrix A * 2:\n";
    D.print();

    Matrix E = A * B;
    std::cout << "\nMatrix A * B:\n";
    E.print();

    Matrix F = A.transpose();
    std::cout << "\nTranspose of Matrix A:\n";
    F.print();

    return 0;
}