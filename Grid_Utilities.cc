#include "Grid_Utilities.h"
using namespace Eigen;

Grid_Utilities::Grid_Utilities(const Grid &g) : dims(g.get_dims()), dx(g.get_dx()), lo(g.get_lo()), hi(g.get_hi())
{
    VectorXd x = VectorXd::LinSpaced(dims(0) + 1, lo(0), hi(0));
    VectorXd y = VectorXd::LinSpaced(dims(1) + 1, hi(1), lo(1));
    X = (VectorXd::LinSpaced(dims(1) + 1, 1, 1) * x.transpose()).array();
    Y = (y * VectorXd::LinSpaced(dims(0) + 1, 1, 1).transpose()).array();
}