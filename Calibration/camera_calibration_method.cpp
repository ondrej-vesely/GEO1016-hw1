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

    // Check if input is valid
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

    for (int i = 0; i < m; i++) {
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
    
    // Compute the SVD decomposition of P-matrix
    Matrix<double> U(m, m, 0.0);
    Matrix<double> S(m, n, 0.0); 
    Matrix<double> V(n, n, 0.0); 
    svd_decompose(P, U, S, V);

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

    /*
    // Check the results
    std::cout << "P: \n" << P << std::endl;
    std::cout << "U: \n" << U << std::endl;
    std::cout << "S: \n" << S << std::endl;
    std::cout << "V: \n" << V << std::endl;
    std::cout << "M: \n" << M << std::endl;

    // Check 1: U is orthogonal, so U * U^T must be identity
    std::cout << "U*U^T: \n" << U * transpose(U) << std::endl;

    // Check 2: V is orthogonal, so V * V^T must be identity
    std::cout << "V*V^T: \n" << V * transpose(V) << std::endl;

    // Check 3: S must be a diagonal matrix
    std::cout << "S: \n" << S << std::endl;

    // Check 4: according to the definition, P = U * S * V^T
    std::cout << "P - U * S * V^T: \n" << P - U * S * transpose(V) << std::endl;
    */

    //   Optional: you can check if your M is correct by applying M on the 3D points. If correct, the projected point
    //             should be very close to your input images points.
    std::cout << "\n" << "Reconstructed 2D coordinates:" << "\n";

    for (int i = 0; i < points_2d.size(); i++) {
        vec3 p3 = points_3d[i];
        vec2 p2 = points_2d[i];

        double coords[]{ p3[0], p3[1], p3[2], 1 };
        Matrix<double> pt3_(4, 1, coords);
        Matrix<double> pt2_ = M * pt3_;
        vec2 pt2d(pt2_(0, 0) / pt2_(2, 0), pt2_(1, 0) / pt2_(2, 0));
        vec2 error = p2 - pt2d;
        std::cout << pt2d << "\t <- error: " << error.length() << "\n";
    }


    // Extract intrinsic parameters from M.
    // Define a1, a2, a3
    vec3 a1 = vec3(M(0, 0), M(0, 1), M(0, 2));
    vec3 a2 = vec3(M(1, 0), M(1, 1), M(1, 2));
    vec3 a3 = vec3(M(2, 0), M(2, 1), M(2, 2));

    // rho
    double rho = 1.0 / a3.length();
    // theta
    double theta = acos(-(dot(cross(a1, a3), cross(a2, a3)) / (norm(cross(a1, a3)) * norm(cross(a2, a3)))));
    
    // cx - u0
    cx = pow(rho, 2) * dot(a1, a3);
    // cy = v0
    cy = pow(rho, 2) * dot(a2, a3);
    // fx - alpha
    double alpha = pow(rho, 2) * cross(a1, a3).length() * sin(theta);
    fx = (float) alpha;
    // fy - beta
    double beta = pow(rho, 2) * cross(a2, a3).length() * sin(theta);
    fy = (float) beta;
    // skew
    skew = (float)(- fx * cos(theta));
    
    std::cout << "\n" << "Intrinsic parameters:"
        << "\n" << "cx / u0: " << "\t" << cx
        << "\n" << "cy / v0: " << "\t" << cy
        << "\n" << "fx / alpha: " << "\t" << fx
        << "\n" << "fy / beta: " << "\t" << fy
        << "\n" << "skew: " << "\t\t" << skew
        << "\n" << "theta: " << "\t\t" << theta
        << "\n" << "ro: " << "\t\t" << rho
        << "\n\n";

    
    // Extract extrinsic parameters from M. 
    // r1, r2, r3 rows
    vec3 r1 = (cross(a2, a3)) / (cross(a2, a3).length());
    vec3 r3 = rho * a3;
    vec3 r2 = cross(r3, r1);

    // Set rows of output R matrix
    R = mat3(
        r1[0], r1[1], r1[2],
        r2[0], r2[1], r2[2],
        r3[0], r3[1], r3[2]
    );

    // Create K matrix
    Matrix<double> K(3, 3, std::vector<double> {
        fx,   -fx * cos(theta)/sin(theta),   cx,
        0,     fy / sin(theta),              cy,
        0,     0,                            1
    }.data());

    // Calculate K inverse
    Matrix<double> K_inv(3, 3, 0.0);
    inverse(K, K_inv);

    // Extract b vector from M matrix
    Matrix<double> b(3, 1, std::vector<double> {
        M(0, 3), M(1, 3), M(2, 3)
    }.data());

    // Calculate T matrix
    auto T = rho * K_inv * b;

    // Set t output matrix
    for (int i = 0; i < 3; i++) { 
        t[i] = (float) T(i, 0); 
    }

    std::cout << "Extrinsic parameters: "
        << "\n" << "R: " << "\n" << R
        << ""   << "T: " << "\n" << t
        << "\n\n";


    return true;
}

















