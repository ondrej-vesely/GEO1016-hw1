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

#ifndef EASY3D_MATRIX_ALGOTHMS_H
#define EASY3D_MATRIX_ALGOTHMS_H


#include "matrix.h"

namespace easy3d {

    /**
     * Compute the inverse of a square matrix. This is a wrapper around Eigen's inverse function.
     * @param A The input matrix.
     * @param invA The inverse of A.
     * @return false if failed (failure is reported only if the input matrix is not square).
     */
    bool inverse(const Matrix<double> &A, Matrix<double> &invA);


    /**
     * Compute the Singular Value Decomposition (SVD) of an M by N matrix. This is a wrapper around Eigen's JacobiSVD.
     *
     * For an m-by-n matrix A, the singular value decomposition is an m-by-m orthogonal matrix U, an m-by-n diagonal
     * matrix S, and an n-by-n orthogonal matrix V so that A = U*S*V^T.
     *
     * The singular values, s[k] = S[k][k], are sorted in decreasing order, i.e., sigma[i] >= sigma[i+1] for any i.
     *
     * The singular value decomposition always exists, so the decomposition will never fail.
     *
     * @param A The input matrix to be decomposed, which can have an arbitrary size.
     * @param U The left side M by M orthogonal matrix.
     * @param S The middle M by N diagonal matrix, with zero elements outside of its main diagonal。
     * @param V The right side N by N orthogonal matrix V.
     */
    void svd_decompose(const Matrix<double> &A, Matrix<double> &U, Matrix<double> &S, Matrix<double> &V);

}

#endif // EASY3D_MATRIX_ALGOTHMS_H
