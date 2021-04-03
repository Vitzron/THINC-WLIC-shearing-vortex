/*
 * File:    THINC_WLIC.h
 * Desc:    header file of class THINC_WLIC
 */

#ifndef THINC_WLIC_
#define THINC_WLIC_

#include <Eigen/Dense>

class THINC_WLIC
{
    public:
        THINC_WLIC  () = default;
        ~THINC_WLIC () = default;
        void update (Eigen::ArrayXXd &phi, const Eigen::ArrayXXd &velo_x, const Eigen::ArrayXXd &velo_y,
                     const double &dt, const double &dx, const double &dy);
    private:
        void flux   (double &f, const Eigen::RowVector4d &phi, const double &u,
                     const Eigen::RowVector2d &nx, const Eigen::RowVector2d &nn,
                     const double &dt, const double &dx) const;
        void normal (Eigen::ArrayXXd &nx, Eigen::ArrayXXd &ny, const Eigen::ArrayXXd &phi,
                     const double &dx, const double &dy) const;
        void bc     (Eigen::ArrayXXd &phi) const;
        size_t step = 0;
};

#endif