/*
 * File: Grid_Utilities.h
 * Desc: header file of class Grid_Utiliites
 *       auxiliary functions of class Grid
 *       staggered arrangement on the grid:
 *       **********fy**********
 *       *         *         *
 *       *         *         *
 *       fx        ct        fx
 *       *         *         *
 *       *         *         *
 *       **********fy**********
 *       ct: cell center
 *       fx: cell faces perpendicular to x-direction
 *       fy: cell faces perpendicular to y-direction
 */

#ifndef GRID_UTILITIES_
#define GRID_UTILITIES_

#include "Grid.h"

class Grid_Utilities
{
    public:
        Grid_Utilities              () = default;
        Grid_Utilities              (const Grid &g);
        ~Grid_Utilities             () = default;
        Eigen::ArrayXXd get_X_fx    () const { return X.bottomRows(dims(1)); };
        Eigen::ArrayXXd get_Y_fx    () const { return Y.bottomRows(dims(1)) + dx(1) / 2.0; };
        Eigen::ArrayXXd get_X_fy    () const { return X.leftCols(dims(0)) + dx(0) / 2.0; };
        Eigen::ArrayXXd get_Y_fy    () const { return Y.leftCols(dims(0)); };
        Eigen::ArrayXXd get_X_ct    () const { return X.bottomLeftCorner(dims(1), dims(0)) + dx(0) / 2.0; };
        Eigen::ArrayXXd get_Y_ct    () const { return Y.bottomLeftCorner(dims(1), dims(0)) + dx(1) / 2.0; };
    private:
        Eigen::Array2i  dims;
        Eigen::Array2d  dx;
        Eigen::Array2d  lo;
        Eigen::Array2d  hi;
        Eigen::ArrayXXd X;
        Eigen::ArrayXXd Y;
};

#endif