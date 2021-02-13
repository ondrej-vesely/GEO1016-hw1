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

#ifndef GEO1016_CAMERA_CALIBRATION_H
#define GEO1016_CAMERA_CALIBRATION_H

#include <easy3d/viewer/viewer.h>


class CameraCalibration : public easy3d::Viewer
{
public:
    CameraCalibration(const std::string& title, const std::string& model_file);

protected:
    bool open() override;
    std::string usage() const override ;
    bool key_press_event(int key, int modifiers) override;
    bool mouse_press_event(int x, int y, int button, int modifiers) override;
    bool mouse_drag_event(int x, int y, int dx, int dy, int button, int modifiers) override;

    void create_cameras_drawable();

    bool calibration(
            const std::vector<easy3d::vec3>& points_3d,
            const std::vector<easy3d::vec2>& points_2d,
            float& fx, float& fy,
            float& cx, float& cy,
            float& skew,
            easy3d::mat3& R,
            easy3d::vec3& t);

private:
    std::vector<easy3d::vec3>  points_3d_;
    std::vector<easy3d::vec2>  points_2d_;
};


#endif // GEO1016_CAMERA_CALIBRATION_H
