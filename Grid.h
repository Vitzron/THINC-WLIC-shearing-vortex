/*
 * File: Grid.h
 * Desc: header file of class Grid
 *       staggered arrangement on the grid:
 *       **********S**********
 *       *         *         *
 *       *         *         *
 *       W         P         E
 *       *         *         *
 *       *         *         *
 *       **********N**********
 *       the scalar variables (pressure, density, volume of fraction etc.) are stored in the cell centers (P)
 *       the velocity variables are located at the cell faces (N, W, S, E)
 */

#ifndef GRID_
#define GRID_

#include <iostream>
#include <Eigen/Dense>

class Grid
{
    friend std::ostream &operator<< (std::ostream &os, const Grid &g);
    public:
        Grid                        () = default;
        Grid                        (const Eigen::Array2i &x_dims, const Eigen::Array2d &x_lo,
                                     const Eigen::Array2d &x_hi, const size_t num_gc = 1);
        ~Grid                       () = default;
        size_t get_ng               () const { return num_ghostcells; };
        Eigen::Array2d get_lo       () const { return lo; };
        Eigen::Array2d get_hi       () const { return hi; };
        Eigen::Array2i get_dims     () const { return dims; };
        Eigen::Array2d get_dx       () const { return dx; };
    private:
        size_t                      num_ghostcells;
        Eigen::Array2d              lo;
        Eigen::Array2d              hi;
        Eigen::Array2i              dims;
        Eigen::Array2d              dx;
};

#endif