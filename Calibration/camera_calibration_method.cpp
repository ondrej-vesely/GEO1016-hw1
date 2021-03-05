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

#include "camera_calibration.h"
#include "matrix_algo.h"

using namespace easy3d;


/**
 * TODO: Finish this function for calibrating a camera from the corresponding 3D-2D point pairs.
 *       You may define a few functions for some sub-tasks.
 *
 * @param points_3d   An array of 3D points.
 * @param points_2d   An array of 2D points.
 * @return True on success, otherwise false. On success, the camera parameters are returned by
 *           - fx and fy: the focal length (in our slides, we use 'alpha' and 'beta'),
 *           - cx and cy: the principal point (in our slides, we use 'u0' and 'v0'),
 *           - skew:      the skew factor ('-alpha * cot_theta')
 *           - R:         the 3x3 rotation matrix encoding camera orientation.
 *           - t:         a 3D vector encoding camera location.
 */
bool CameraCalibration::calibration(
        const std::vector<vec3>& points_3d,
        const std::vector<vec2>& points_2d,
        float& fx, float& fy,
        float& cx, float& cy,
        float& skew,
        mat3& R,
        vec3& t)
{
    std::cout << std::endl;
    std::cout << "TODO: I am going to implement the calibration() function in the following file:" << std::endl
              << "\t" << __FILE__ << std::endl;
    std::cout << "TODO: After implementing the calibration() function, I will disable all unrelated output ...\n\n";

    // check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points must match)
    if (points_2d.size() >= 6 && points_2d.size() == points_3d.size()) {
        std::cout << "\t"<< "Input is valid" << std::endl;
    }
    else {
        std::cout << "\t" << "Input is not valid" << std::endl;
        return false;
    }

    // Construct the P matrix (so P * m = 0).
    const int m = 2 * points_3d.size();
    const int n = 12;
    Matrix<double> P(m, n, 0.0);

    for (int i = 0; i < 2*points_2d.size(); i++) {
        const int h = i/2;              // half of i
        const int shift = (i % 2) * 4;  // if row number is odd, shift first 4 by 4 to left

        P(i, 0 + shift) = points_3d[h][0];
        P(i, 1 + shift) = points_3d[h][1];
        P(i, 2 + shift) = points_3d[h][2];
        P(i, 3 + shift) = 1.0;
        P(i, 8) = points_3d[h][0] * -1.0 * points_2d[h][i % 2];
        P(i, 9) = points_3d[h][1] * -1.0 * points_2d[h][i % 2];
        P(i, 10) = points_3d[h][2] * -1.0 * points_2d[h][i % 2];
        P(i, 11) = -1.0 * points_2d[h][i % 2];
    }
    std::cout << "P-Matrix: \n" << P << std::endl;

    
    // TODO: solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.
    // Compute the SVD decomposition of P-matrix
    Matrix<double> U(m, m, 0.0);
    Matrix<double> S(m, n, 0.0); 
    Matrix<double> V(n, n, 0.0); 
    svd_decompose(P, U, S, V);

    // Check the results
    std::cout << "U: \n" << U << std::endl;
    std::cout << "S: \n" << S << std::endl;
    std::cout << "V: \n" << V << std::endl;

    // Check 1: U is orthogonal, so U * U^T must be identity
    std::cout << "U*U^T: \n" << U * transpose(U) << std::endl;

    // Check 2: V is orthogonal, so V * V^T must be identity
    std::cout << "V*V^T: \n" << V * transpose(V) << std::endl;

    // Check 3: S must be a diagonal matrix
    std::cout << "S: \n" << S << std::endl;

    // Check 4: according to the definition, P = U * S * V^T
    std::cout << "P - U * S * V^T: \n" << P - U * S * transpose(V) << std::endl;



    // Initialise M
    Matrix<double> M(3, 4, 0.0);

    // Fill M with last column of V
    int k = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            M(i, j) = V(k, 11);
            k++;
        }
    }

    std::cout << "M: \n" << M << std::endl;

    //   Optional: you can check if your M is correct by applying M on the 3D points. If correct, the projected point
    //             should be very close to your input images points.


    // Extract intrinsic parameters from M.
    // Define a1, a2, a3
    vec3 a1 = vec3(M(0, 0), M(0, 1), M(0, 2));
    vec3 a2 = vec3(M(1, 0), M(1, 1), M(1, 2));
    vec3 a3 = vec3(M(2, 0), M(2, 1), M(2, 2));

    double rho = 1 / a3.length();
    cx = pow(rho, 2) * dot(a1, a3);         // cx and cy have been declared in line 50. Is this the right way to "fill in" the values?
    cy = pow(rho, 2) * dot(a2, a3);
    
    // theta
    double theta = acos(-(dot(cross(a1, a3), cross(a2, a3)) / (norm(cross(a1, a3)) * norm(cross(a2, a3)))));
    // fx
    double alpha = pow(rho, 2) * cross(a1, a3).length() * sin(theta);
    fx = (float) alpha;
    // fy
    double beta = pow(rho, 2) * cross(a2, a3).length() * sin(theta);
    fy = (float) beta;
    // skewness
    skew = (float)(- fx * cos(theta));
    std::cout << "Intrinsic parameters: " " cx: " << cx << " cy: " << cy << " theta: " << theta << " fx: " << fx << " fy: " << fy << " skewness: " << skew << std::endl;
    // TODO: extract extrinsic parameters from M.
    // b1, b2, b3
    double b1 = M(0, 3);
    double b2 = M(1, 3);
    double b3 = M(2, 3);
    // r1, r2, r3
    auto r1 = (cross(a2, a3)) / (cross(a2, a3).length());
    auto r3 = rho * a3;
    auto r2 = cross(r3, r1);

    // TODO: uncomment the line below to return true when testing your algorithm and in you final submission.
    return false;


    
    // TODO: The following code is just an example showing you SVD decomposition, matrix inversion, and some related.
    // TODO: Delete the code below (or change "#if 1" in the first line to "#if 0") in you final submission.
#if 1
    std::cout << "[Liangliang:] Camera calibration requires computing the SVD and inverse of matrices.\n"
                 "\tIn this assignment, I provide you with a Matrix data structure for storing matrices of arbitrary\n"
                 "\tsizes (see matrix.h). I also wrote the example code to show you how to:\n"
                 "\t\t- use the dynamic 1D array data structure 'std::vector' from the standard C++ library;\n"
                 "\t\t  The points (both 3D and 2D) are stored in such arrays;\n"
                 "\t\t- use the template matrix class (which can have an arbitrary size);\n"
                 "\t\t- compute the SVD of a matrix;\n"
                 "\t\t- compute the inverse of a matrix;\n"
                 "\t\t- compute the transpose of a matrix.\n"
                 "\tThe following are just the output of these examples. You should delete ALL unrelated code and\n"
                 "\tavoid unnecessary output in you final submission.\n\n";

    // This is a 1D array of 'double' values. Alternatively, you can use 'double mat[25]' but you cannot change it
    // length. With 'std::vector', you can do append/delete/insert elements, and much more. The 'std::vector' can store
    // not only 'double', but also any other types of objects. In case you may want to learn more about 'std::vector'
    // check here: https://en.cppreference.com/w/cpp/container/vector
    std::vector<double> array = {1, 3, 3, 4, 7, 6, 2, 8, 2, 8, 3, 2, 4, 9, 1, 7, 3, 23, 2, 3, 5, 2, 1, 5, 8, 9, 22};
    array.push_back(5); // append 5 to the array (so the size will increase by 1).
    array.insert(array.end(), 10, 3);  // append ten 3 (so the size will grow by 10).

    // To access its values
    for (int i=0; i<array.size(); ++i)
        std::cout << array[i] << " ";  // use 'array[i]' to access its i-th element.
    std::cout << std::endl;

    // Define an m-by-n double valued matrix.
    // Here I use the above array to initialize it. You can also use A(i, j) to initialize/modify/access its elements.
    //const int m = 2*points_3d.size(), n = 12;
    Matrix<double> A(m, n, 0.0);    // 'array.data()' returns a pointer to the array.
    std::cout << "M: \n" << A << std::endl;

    //Matrix<double> U(m, m, 0.0);   // initialized with 0s
    //Matrix<double> S(m, n, 0.0);   // initialized with 0s
    //Matrix<double> V(n, n, 0.0);   // initialized with 0s

    // Compute the SVD decomposition of A
    //svd_decompose(A, U, S, V);

    // Now let's check if the SVD result is correct

    // Check 1: U is orthogonal, so U * U^T must be identity
    std::cout << "U*U^T: \n" << U * transpose(U) << std::endl;

    // Check 2: V is orthogonal, so V * V^T must be identity
    std::cout << "V*V^T: \n" << V * transpose(V) << std::endl;

    // Check 3: S must be a diagonal matrix
    std::cout << "S: \n" << S << std::endl;

    // Check 4: according to the definition, A = U * S * V^T
    std::cout << "M - U * S * V^T: \n" << A - U * S * transpose(V) << std::endl;

    // Define a 5 by 5 square matrix and compute its inverse.
    Matrix<double> B(5, 5, array.data());    // Here I use part of the above array to initialize B
    // Compute its inverse
    Matrix<double> invB(5, 5);
    inverse(B, invB);
    // Let's check if the inverse is correct
    std::cout << "B * invB: \n" << B * invB << std::endl;

    return false;
    // TODO: delete the above code in you final submission (which are just examples).
#endif
}

















