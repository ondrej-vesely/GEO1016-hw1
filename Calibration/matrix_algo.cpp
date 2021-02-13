/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */


#include "matrix_algo.h"
#include <iostream>
#include <3rd_party/Eigen/Dense>


namespace easy3d {

    bool inverse(const Matrix<double> &A, Matrix<double> &invA) {
        const int m = A.rows();
        const int n = A.cols();
        if (m != n) {
            std::cout << "could not compute inverse (not a square matrix)" << std::endl;
            return false;
        }

        Eigen::MatrixXd C(m, n);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j)
                C(i, j) = A(i, j);
        }

        invA.resize(m, n);
        const Eigen::MatrixXd invC = C.inverse();
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j)
                invA(i, j) = invC(i, j);
        }
        return true;
    }


    void svd_decompose(const Matrix<double> &A, Matrix<double> &U, Matrix<double> &S, Matrix<double> &V)
    {
        const int m = A.rows();
        const int n = A.cols();
        Eigen::MatrixXd C(m, n);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j)
                C(i, j) = A(i, j);
        }

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(C, Eigen::ComputeFullU | Eigen::ComputeFullV);
        const Eigen::MatrixXd& eU = svd.matrixU();
        const Eigen::MatrixXd& eV = svd.matrixV();
        const auto& eS = svd.singularValues();

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m; ++j)
                U(i, j) = eU(i, j);
        }

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j)
                V(i, j) = eV(i, j);
        }

        S.load_zero();
        for (int i = 0; i < std::min(m, n); ++i)
            S(i, i) = eS(i);
    }
}