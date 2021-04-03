#include <iostream>
#include <cmath>
#include "THINC_WLIC.h"
using namespace std;
using namespace Eigen;

const double beta = 3.5;
const double eps = 1.0e-10;

void THINC_WLIC::update(Eigen::ArrayXXd &phi, const Eigen::ArrayXXd &velo_x, const Eigen::ArrayXXd &velo_y,
        const double &dt, const double &dx, const double &dy)
{
    if (phi.rows() != velo_x.rows() || phi.cols() + 1 != velo_x.cols() ||
        phi.rows() + 1 != velo_y.rows() || phi.cols() != velo_y.cols())
    {
        cerr << "INCOMPATIBLE shapes in THINC_SW!" << endl;
        exit(EXIT_FAILURE);
    }

    Index rows = phi.rows() - 2;
    Index cols = phi.cols() - 2;

    ArrayXXd nx, ny, nn;
    ArrayXXd phi_flux;

    ArrayXXd phi_0 = phi.block(1, 1, rows, cols);
    if (step % 2 == 0)
    {
        // x-direction
        normal(nx, ny, phi, dx, dy);
        nn = nx + ny;
        phi_flux = ArrayXXd::Zero(rows, cols + 1);
        for (Index j = 0; j != cols - 1; ++j)
        {
            for (Index i = 0; i != rows; ++i)
            {
                flux(phi_flux(i, j + 1), phi.block(i + 1, j, 1, 4).matrix(), velo_x(i + 1, j + 2), nx.block(i, j, 1, 2).matrix(), nn.block(i, j, 1, 2), dt, dx);
            }
        }
        phi.block(1, 1, rows, cols) = phi.block(1, 1, rows, cols)
            - (phi_flux.rightCols(cols) - phi_flux.leftCols(cols)) / dx
            + phi_0 * (velo_x.block(1, 2, rows, cols) - velo_x.block(1, 1, rows, cols)) * dt / dx;
        bc(phi);

        // y-direction
        normal(nx, ny, phi, dx, dy);
        nn = nx + ny;
        phi_flux = ArrayXXd::Zero(rows + 1, cols);
        for (Index j = 0; j != cols; ++j)
        {
            for (Index i = 0; i != rows - 1; ++i)
            {
                flux(phi_flux(i + 1, j), phi.block(i, j + 1, 4, 1).transpose().reverse().matrix(), velo_y(i + 2, j + 1),
                    ny.block(i, j, 2, 1).transpose().reverse().matrix(), nn.block(i, j, 2, 1).transpose().reverse().matrix(), dt, dy);
            }
        }
        phi.block(1, 1, rows, cols) = phi.block(1, 1, rows, cols) 
            - (phi_flux.topRows(rows) - phi_flux.bottomRows(rows)) / dx
            + phi_0 * (velo_y.block(1, 1, rows, cols) - velo_y.block(2, 1, rows, cols)) * dt / dy;
        bc(phi);
    }
    else
    {
        // y-direction
        normal(nx, ny, phi, dx, dy);
        nn = nx + ny;
        phi_flux = ArrayXXd::Zero(rows + 1, cols);
        for (Index j = 0; j != cols; ++j)
        {
            for (Index i = 0; i != rows - 1; ++i)
            {
                flux(phi_flux(i + 1, j), phi.block(i, j + 1, 4, 1).transpose().reverse().matrix(), velo_y(i + 2, j + 1),
                    ny.block(i, j, 2, 1).transpose().reverse().matrix(), nn.block(i, j, 2, 1).transpose().reverse().matrix(), dt, dy);
            }
        }
        phi.block(1, 1, rows, cols) = phi.block(1, 1, rows, cols) 
            - (phi_flux.topRows(rows) - phi_flux.bottomRows(rows)) / dx
            + phi_0 * (velo_y.block(1, 1, rows, cols) - velo_y.block(2, 1, rows, cols)) * dt / dy;
        bc(phi);

        // x-direction
        normal(nx, ny, phi, dx, dy);
        nn = nx + ny;
        phi_flux = ArrayXXd::Zero(rows, cols + 1);
        for (Index j = 0; j != cols - 1; ++j)
        {
            for (Index i = 0; i != rows; ++i)
            {
                flux(phi_flux(i, j + 1), phi.block(i + 1, j, 1, 4).matrix(), velo_x(i + 1, j + 2), nx.block(i, j, 1, 2).matrix(), nn.block(i, j, 1, 2), dt, dx);
            }
        }
        phi.block(1, 1, rows, cols) = phi.block(1, 1, rows, cols)
            - (phi_flux.rightCols(cols) - phi_flux.leftCols(cols)) / dx
            + phi_0 * (velo_x.block(1, 2, rows, cols) - velo_x.block(1, 1, rows, cols)) * dt / dx;
        bc(phi);
    }
    ++step;
}

void THINC_WLIC::bc(Eigen::ArrayXXd &phi) const
{
    Index rows = phi.rows();
    Index cols = phi.cols();

    phi.row(0) = phi.row(1);
    phi.row(rows - 1) = phi.row(rows - 2);
    phi.col(0) = phi.col(1);
    phi.col(cols - 1) = phi.col(cols - 2);
}

void THINC_WLIC::flux(double &f, const Eigen::RowVector4d &phi, const double &u,
        const Eigen::RowVector2d &nx, const Eigen::RowVector2d &nn,
        const double &dt, const double &dx) const
{
    if (u == 0.0)
        f = 0.0;
    else
    {
        Index iup = u > 0.0 ? 1 : 2;
        if (phi(iup) > 1.0 - eps || phi(iup) < eps || nn(iup - 1) == 0.0)
            f = phi(iup) * u * dt;
        else
        {
            double gamma  = phi(iup - 1) <= phi(iup + 1) ? 1.0 : -1.0;
            double lambda = u < 0.0 ? 0.0 : 1.0;
            double x_m = 0.5 / beta * log((exp(beta / gamma * (1.0 + gamma - 2.0 * phi(iup))) - 1.0) 
                / (1.0 - exp(beta / gamma * (1.0 - gamma - 2.0 * phi(iup)))));
            double fxx = 0.5 * (u * dt - gamma * dx / beta * log((cosh(beta * (lambda - x_m - u * dt / dx))) 
                / (cosh(beta * (lambda - x_m)))));
            double weight = nx(iup - 1) / nn(iup - 1);
            f = fxx * weight + phi(iup) * u * dt * (1.0 - weight);
        }
    }
}

void THINC_WLIC::normal(Eigen::ArrayXXd &nx, Eigen::ArrayXXd &ny, const Eigen::ArrayXXd &phi,
        const double &dx, const double &dy) const
{
    Index rows = phi.rows() - 2;
    Index cols = phi.cols() - 2;

    ArrayXXd mx_v = (phi.block(0, 1, rows + 1, cols + 1) + phi.block(1, 1, rows + 1, cols + 1)
        - phi.block(0, 0, rows + 1, cols + 1) - phi.block(1, 0, rows + 1, cols + 1)) / dx / 2.0;
    ArrayXXd my_v = (phi.block(0, 0, rows + 1, cols + 1) + phi.block(0, 1, rows + 1, cols + 1)
        - phi.block(1, 0, rows + 1, cols + 1) - phi.block(1, 1, rows + 1, cols + 1)) / dy / 2.0;
    nx = (mx_v.block(0, 0, rows, cols) + mx_v.block(0, 1, rows, cols)
          + mx_v.block(1, 0, rows, cols) + mx_v.block(1, 1, rows, cols)).abs() / 4.0;
    ny = (my_v.block(0, 0, rows, cols) + my_v.block(0, 1, rows, cols)
          + my_v.block(1, 0, rows, cols) + my_v.block(1, 1, rows, cols)).abs() / 4.0;
    
}