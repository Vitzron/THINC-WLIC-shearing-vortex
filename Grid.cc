#include "Grid.h"

Grid::Grid(const Eigen::Array2i &x_dims, const Eigen::Array2d &x_lo,
        const Eigen::Array2d &x_hi, const size_t num_gc)
{
    num_ghostcells = num_gc;
    dims = x_dims + 2 * num_gc;

    dx(0) = (x_hi(0) - x_lo(0)) / x_dims(0);
    dx(1) = (x_hi(1) - x_lo(1)) / x_dims(1);
    lo = x_lo - dx * num_ghostcells;
    hi = x_hi + dx * num_ghostcells;
}

std::ostream &operator<<(std::ostream &os, const Grid &g)
{
    os << "Geometric limits: [" << g.lo(0) << ", " << g.hi(0) << "] x ["
       << g.lo(1) << ", " << g.hi(1) << "]" << std::endl;
    os << "Grid resolution: " << g.dims(1) << " x " << g.dims(0) << std::endl;
    os << "Grid space: dx = " << g.dx(0) << ", " << g.dx(1) << std::endl;
    os << "Number of ghostcells: " << g.num_ghostcells << std::endl;
	return os;
}