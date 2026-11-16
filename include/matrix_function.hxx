// matrix_function.hxx

#pragma once

#include <Eigen/Dense>

#include <vector>
#include <stdexcept>
#include <cstddef>

class MatrixFunction {
public:
    using Matrix = Eigen::MatrixXcd;

    // Constructor
    MatrixFunction(const std::vector<double>& x, std::size_t nrows, std::size_t ncols)
        : x_(x), vals_(x.size(), Matrix::Zero(nrows, ncols)), nrows_(nrows), ncols_(ncols)
    {
        if (nrows == 0 || ncols == 0) throw std::invalid_argument("dimensions must be positive");
    }
    MatrixFunction(const std::vector<double>& x, std::size_t n) : MatrixFunction(x, n, n) {}

    // Accessors
    std::size_t size() const { return x_.size(); }
    std::size_t nrows() const { return nrows_; }
    std::size_t ncols() const { return ncols_; }
    const std::vector<double>& x() const { return x_; }

    const Matrix& operator[](std::size_t i) const { return vals_.at(i); }
    Matrix& operator[](std::size_t i) { return vals_.at(i); }

    // Right matrix multiplication
    MatrixFunction matmul_right(const Matrix& C) const {
        if (C.rows() != static_cast<Eigen::Index>(ncols_) || C.cols() == 0) throw std::invalid_argument("Matrix dimensions do not match");

        MatrixFunction result(x_, nrows_, C.cols());
        for (std::size_t i = 0; i < size(); ++i) result[i] = vals_[i] * C;

        return result;
    }

    // Left matrix multiplication
    MatrixFunction matmul_left(const Matrix& C) const {
        if (C.cols() != static_cast<Eigen::Index>(nrows_) || C.rows() == 0) throw std::invalid_argument("Matrix dimensions do not match");

        MatrixFunction result(x_, C.rows(), ncols_);
        for (std::size_t i = 0; i < size(); ++i) result[i] = C * vals_[i];

        return result;
    }

private:
    std::vector<double> x_;
    std::vector<Matrix> vals_;
    std::size_t nrows_;
    std::size_t ncols_;

};