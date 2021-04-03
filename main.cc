#include "Grid.h"
#include "Grid_Utilities.h"
#include "THINC_WLIC.h"
#include "Eigen_Utilities.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;
using namespace Eigen;

Eigen::ArrayXXd Heaviside(const Eigen::ArrayXXd &phi, const double &eps);
double Heaviside(const double &phi, const double &eps);

int main(int argc, char* argv[])
{
    size_t rows = atoi(argv[1]);
    size_t cols = atoi(argv[1]);
    Array2i dims(cols, rows);
    Array2d lo(0., 0.);
    Array2d hi(1.0, 1.0);
    Grid grid(dims, lo, hi);
    double dx = grid.get_dx()(0);
    double dy = grid.get_dx()(1);
    cout << "**************************grid**************************"
         << endl << grid 
         << "**************************grid**************************"
         << endl;
    
    Grid_Utilities gu(grid);
    ArrayXXd X = gu.get_X_ct();
    ArrayXXd Y = gu.get_Y_ct();
    ArrayXXd X_fx = gu.get_X_fx();
    ArrayXXd Y_fx = gu.get_Y_fx();
    ArrayXXd X_fy = gu.get_X_fy();
    ArrayXXd Y_fy = gu.get_Y_fy();

    ofstream output;
    const char *X_file = "X.dat";
    const char *Y_file = "Y.dat";

    output.open(X_file, std::ios::out);
    if (!output)
    {
        cerr << "CANNOT open file " << X_file << "!" << endl;
        exit(EXIT_FAILURE);
    }
    output << X;
    output.flush();
    output.close();
    cout << "SUCCESSFULLY write X into file " << X_file << "!" << endl;
    output.open(Y_file, std::ios::out);
    if (!output)
    {
        cerr << "CANNOT open file " << Y_file << "!" << endl;
        exit(EXIT_FAILURE);
    }
    output << Y;
    output.flush();
    output.close();
    cout << "SUCCESSFULLY write Y into file " << Y_file << "!" << endl;

    ArrayXXd alpha = ArrayXXd::Zero(rows + 2, cols + 2);
    ArrayXXd u = ArrayXXd::Zero(rows + 2, cols + 3);
    ArrayXXd v = ArrayXXd::Zero(rows + 3, cols + 2);
    ArrayXXd p = ArrayXXd::Zero(rows + 2, cols + 2);
    ArrayXXd rho = ArrayXXd::Zero(rows + 2, cols + 2);

    const double x_lo_w = 5.0;
    const double x_hi_w = 15.0;
    const double y_lo_w = 2.0;
    const double y_hi_w = 10.0;
    const double c_x = 0.50;
    const double c_y = 0.75;
    const double radius = 0.15;
    Eigen_Utilities eu = Eigen_Utilities();
    alpha = 1.0 - eu.positive(X - x_lo_w) * eu.positive(x_hi_w - X) * eu.positive(Y - y_lo_w) * eu.positive(y_hi_w - Y);
    alpha = Heaviside(((X - c_x).square() + (Y - c_y).square()).sqrt() - radius, dx * 1.5);

    THINC_WLIC thinc = THINC_WLIC();

    const double CFL = 0.1;
    const double dt = CFL * min(dx, dy);
    cout << "dx: " << dx << "; dy: " << dy << "; dt: " << dt << endl;

    double T = 0.0;
    const double T_total = 8.0;
    ArrayXd t_records(5);
    t_records << 0.0, 2.0, 4.0, 6.0, 8.0;
    size_t i = 0;
    size_t index_record = 0;

    vector<string> phi_files;
    phi_files.push_back("phi_0.dat");
    phi_files.push_back("phi_2.dat");
    phi_files.push_back("phi_4.dat");
    phi_files.push_back("phi_6.dat");
    phi_files.push_back("phi_8.dat");


    while(T <= T_total)
    {
        // cout << "step: " << i << "; dt: " << dt << "; T: " << T;
        // cout << "; min of alpha: " << alpha.minCoeff() << "; max of alpha: " << alpha.maxCoeff() << "; initial total mass: " << (rows * cols - alpha.block(1, 1, rows, cols).sum()) * dx * dy << endl;

        if (T <= t_records(index_record) && T + dt > t_records(index_record))
        {
            output.open(phi_files[index_record], std::ios::out);
            if (!output)
            {
                cerr << "CANNOT open file " << phi_files[index_record] << "!" << endl;
                exit(EXIT_FAILURE);
            }
            output << alpha;
            output.flush();
            output.close();
            ++index_record;
            cout << "**********************************************************" << endl;
            cout << "**************Write data down at time " << T << "****************" << endl;
            cout << "**********************************************************" << endl;
        }
        
        u = (2.0 * EIGEN_PI * Y_fx).sin() * (EIGEN_PI * X_fx).sin().square() * cos(EIGEN_PI * T / T_total);
        v = -(2.0 * EIGEN_PI * X_fy).sin() * (EIGEN_PI * Y_fy).sin().square() * cos(EIGEN_PI * T / T_total);
        thinc.update(alpha, u, v, dt, dx, dy);

        T += dt;
        ++i;
    }

    return 0;
}


Eigen::ArrayXXd Heaviside(const Eigen::ArrayXXd &phi, const double &eps)
{
    ArrayXXd H = phi;
    for (Index j = 0; j != phi.cols(); ++j)
    {
        for (Index i = 0; i != phi.rows(); ++i)
        {
            H(i, j) = Heaviside(phi(i, j), eps);
        }
    }
    return H;
}

double Heaviside(const double &phi, const double &eps)
{
    if (phi < -eps)
    {
        return 0.0;
    }
    else if (phi > eps)
    {
        return 1.0;
    }
    else
    {
        return 0.5 + phi / 2.0 / eps + 1.0 / 2.0 / EIGEN_PI * sin(EIGEN_PI * phi / eps);
    }
}