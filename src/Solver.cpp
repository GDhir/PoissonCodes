#include "Solver.hpp"
#include "utils.hpp"
#include <iostream>

using std::cout;

void Solver::solveChannelFlow(Grid &grid, std::vector<double> &dataValue, const ConstantData<double> &cvals, const double tol)
{

    auto neta = grid.getNeta();
    auto nzeta = grid.getNzeta();
    double k = 0;
    double zeta = 0;
    double eta = 0;
    double c1 = 0, c2 = 0, c3 = 0, c4 = 0;

    if (dataValue.size() != static_cast<size_t>(neta * nzeta))
        throw("Incorrect Initialization");

    double error = pow(10, 6);

    vector<double> prevValue;

    int iter = 0;

    while ( error > tol )
    {
        iter += 1;
        prevValue = dataValue;

        // Implement Neumann BCs
        for (int etaID = 1; etaID < neta; etaID++)
        {

            dataValue[etaID * nzeta + 0] = (-dataValue[(etaID)*nzeta + 2] +
                                            4 * dataValue[(etaID)*nzeta + 1]) /
                                           3.0;
        }

        int j = neta - 1;
        int i = nzeta - 2;

        zeta = grid.calcXFromIndex(i, nzeta);
        k = (2 * cvals.a * zeta) / (2 * cvals.a + cvals.b);

        dataValue[j * nzeta + i] = (4 * k * dataValue[(j)*nzeta + i + 1] - k * dataValue[(j)*nzeta + i + 2] + 4 * dataValue[(j - 1) * nzeta + i] - dataValue[(j - 2) * nzeta + i]) / (3.0 * k + 3);

        double spacing = grid.getH();

        for (int j = 1; j < neta - 1; j++)
        {
            for (int i = 1; i < nzeta - 1; i++)
            {

                eta = grid.calcXFromIndex(j, neta);
                zeta = grid.calcXFromIndex(i, nzeta);

                double cdenom = 2 * cvals.a * eta + cvals.b;

                c1 = 8 * pow(cvals.a, 2) * zeta / (pow(cvals.h * cdenom, 2));

                c2 = 4 * (pow(cvals.a * zeta, 2) + pow(cvals.h, 2)) / (pow(cdenom * cvals.h, 2));

                c3 = -4 * cvals.a * zeta / (cvals.h * cvals.h * cdenom);

                c4 = 1 / (cvals.h * cvals.h);

                dataValue[j * nzeta + i] = (c1 * (dataValue[j * nzeta + i + 1] - dataValue[j * nzeta + i - 1]) / (2 * spacing) +
                                            c2 * (dataValue[j * nzeta + i + 1] + dataValue[j * nzeta + i - 1]) / (spacing * spacing) +
                                            c3 * (dataValue[(j + 1) * nzeta + i + 1] - dataValue[(j - 1) * nzeta + i + 1] - dataValue[(j + 1) * nzeta + i - 1] + dataValue[(j - 1) * nzeta + i - 1]) / (4 * spacing * spacing) +
                                            c4 * (dataValue[(j + 1) * nzeta + i] + dataValue[(j - 1) * nzeta + i]) / (spacing * spacing) + 1) /
                                           (2 * c2 / (spacing * spacing) + 2 * c4 / (spacing * spacing));
            }
        }

        for (int i = 1; i < nzeta - 2; i++)
        {

            zeta = grid.calcXFromIndex(i, nzeta);
            k = (2 * cvals.a * zeta) / (2 * cvals.a + cvals.b);

            int j = neta - 1;

            dataValue[j * nzeta + i] = (4 * k * dataValue[(j)*nzeta + i + 1] - k * dataValue[(j)*nzeta + i + 2] + 4 * dataValue[(j - 1) * nzeta + i] - dataValue[(j - 2) * nzeta + i]) / (3.0 * k + 3);
        }

        if (std::isnan(dataValue[neta * nzeta / 2]))
        {
            throw("NaNs in output ");
        }

        error = 0;

        // if( iter%50 == 0 ) {

        //     for( int i = 0; i < dataValue.size(); i++ ) {

        //         error += pow( ( dataValue[i] - prevValue[i] ), 2 );

        //     }
        //     cout << error << "\t" << iter << "\n";
        // }
        // else
        //     error = 23;

        if (iter % 50 == 0)
        {

            for (int etaID = 1; etaID < neta; etaID++)
            {

                error += std::pow(dataValue[etaID * nzeta + 0] - (-dataValue[(etaID)*nzeta + 2] +
                                                                  4 * dataValue[(etaID)*nzeta + 1]) /
                                                                     3.0,
                                  2);
            }

            int j = neta - 1;
            int i = nzeta - 2;

            zeta = grid.calcXFromIndex(i, nzeta);
            k = (2 * cvals.a * zeta) / (2 * cvals.a + cvals.b);

            error += pow(dataValue[j * nzeta + i] - (4 * k * dataValue[(j)*nzeta + i + 1] - k * dataValue[(j)*nzeta + i + 2] + 4 * dataValue[(j - 1) * nzeta + i] - dataValue[(j - 2) * nzeta + i]) / (3.0 * k + 3), 2);

            double spacing = grid.getH();

            for (int j = 1; j < neta - 1; j++)
            {
                for (int i = 1; i < nzeta - 1; i++)
                {

                    eta = grid.calcXFromIndex(j, neta);
                    zeta = grid.calcXFromIndex(i, nzeta);

                    double cdenom = 2 * cvals.a * eta + cvals.b;

                    c1 = 8 * pow(cvals.a, 2) * zeta / (pow(cvals.h * cdenom, 2));

                    c2 = 4 * (pow(cvals.a * zeta, 2) + pow(cvals.h, 2)) / (pow(cdenom * cvals.h, 2));

                    c3 = -4 * cvals.a * zeta / (cvals.h * cvals.h * cdenom);

                    c4 = 1 / (cvals.h * cvals.h);

                    error += std::pow(dataValue[j * nzeta + i] - (c1 * (dataValue[j * nzeta + i + 1] - dataValue[j * nzeta + i - 1]) / (2 * spacing) +
                                                                  c2 * (dataValue[j * nzeta + i + 1] + dataValue[j * nzeta + i - 1]) / (spacing * spacing) +
                                                                  c3 * (dataValue[(j + 1) * nzeta + i + 1] - dataValue[(j - 1) * nzeta + i + 1] - dataValue[(j + 1) * nzeta + i - 1] + dataValue[(j - 1) * nzeta + i - 1]) / (4 * spacing * spacing) +
                                                                  c4 * (dataValue[(j + 1) * nzeta + i] + dataValue[(j - 1) * nzeta + i]) / (spacing * spacing) + 1) /
                                                                     (2 * c2 / (spacing * spacing) + 2 * c4 / (spacing * spacing)),
                                      2);
                }
            }

            for (int i = 1; i < nzeta - 2; i++)
            {

                zeta = grid.calcXFromIndex(i, nzeta);
                k = (2 * cvals.a * zeta) / (2 * cvals.a + cvals.b);

                int j = neta - 1;

                error += std::pow(dataValue[j * nzeta + i] - (4 * k * dataValue[(j)*nzeta + i + 1] - k * dataValue[(j)*nzeta + i + 2] + 4 * dataValue[(j - 1) * nzeta + i] - dataValue[(j - 2) * nzeta + i]) / (3.0 * k + 3), 2);
            }

            error = sqrt(error);
            cout << error << "\n";
        }
        else
        {
            error = 23;
        }
    }

    cout << "end"
         << "\n";
}

vector<pair<double, int>> Solver::solveJacobiPoisson(Grid &grid, std::vector<double> &dataValue, vector<double> &fval, const double w, const double tol, int maxiter, cmptype comparator)
{

    double h = grid.getH();
    int N = grid.getNx();

    double res = pow(10, 5);
    vector<pair<double, int>> resConvergence;

    std::vector<double> prevValue;

    int iter = 0;

    while ( comparator( res, tol, iter, maxiter ) )
    {

        prevValue = dataValue;

        for (int j = 1; j < N - 1; j++)
        {

            for (int i = 1; i < N - 1; i++)
            {

                dataValue[j * N + i] = w * (h * h * fval[j * N + i] + prevValue[j * N + i + 1] + prevValue[j * N + i - 1] + prevValue[(j + 1) * N + i] + prevValue[(j - 1) * N + i]) / 4.0 + (1 - w) * dataValue[j * N + i];
            }
        }

        if (iter % 10 == 0)
        {

            res = calculateResidual(dataValue, fval, "Jacobi", h, N);

            resConvergence.push_back(make_pair(res, iter));
            //cout << res << "\n";
        }

        iter += 1;
    }

    return resConvergence;
}

vector<pair<double, int>> Solver::solveGSPoisson(Grid &grid, std::vector<double> &dataValue, vector<double> &fval, const double w, const double tol, int maxiter, cmptype comparator )
{

    double h = grid.getH();
    int N = grid.getNx();

    double res = pow(10, 5);
    vector<pair<double, int>> resConvergence;

    int iter = 0;

    while ( comparator( res, tol, iter, maxiter ) )
    {

        for (int j = 1; j < N - 1; j++)
        {

            for (int i = 1; i < N - 1; i++)
            {

                dataValue[j * N + i] = w * (h * h * fval[j * N + i] + dataValue[j * N + i + 1] + dataValue[j * N + i - 1] + dataValue[(j + 1) * N + i] + dataValue[(j - 1) * N + i]) / 4.0 + (1 - w) * dataValue[j * N + i];
            }
        }

        if (iter % 10 == 0)
        {

            res = calculateResidual(dataValue, fval, "GS", h, N);

            resConvergence.push_back(make_pair(res, iter));
            //cout << res << "\n";
        }

        iter += 1;
    }

    return resConvergence;
}

void Solver::multigridVCycle(Grid refgrid, std::vector<double> &dataValue, vector<double> &fval, int currlevel,
                                      const double w, const double tol, string method, int& totalIter, int niter )
{

    int N = refgrid.getNx();
    double h = refgrid.getH();
    auto comparator = []( double res, double tol, int iter, int maxiter ){ if( iter < maxiter ) return true; else return false; };
        
    if (currlevel == 1)
    {

        solveGSPoisson(refgrid, dataValue, fval, w, tol, niter, comparator);
        totalIter += niter;
        return;
    }

    int coarseN = ceil(N / 2.0);
    Grid coarsegrid{coarseN, coarseN};

    auto currRes = calculateResidualArray(dataValue, fval, method, h, N);
    auto restrictedRes = refgrid.restrictVal(currRes, coarsegrid);

    vector<double> errorVal(coarseN * coarseN, 0);
    multigridVCycle(coarsegrid, errorVal, restrictedRes, currlevel - 1, w, tol, method, totalIter, niter);

    auto interpolatedError = refgrid.interpolateVal(errorVal, coarsegrid);
    addError(dataValue, interpolatedError, N);

    solveGSPoisson(refgrid, dataValue, fval, w, tol, niter, comparator);
    totalIter += niter;

}

vector<pair<double, int>> Solver::computeMultigrid( Grid refgrid, std::vector<double> &dataValue, vector<double> &fval, const int maxlevels,
                                      const double w, const double tol, string method ) {


    int niter = 50;
    int N = refgrid.getNx();
    double h = refgrid.getH();
    int totalIter{0};
    vector<pair<double, int>> resConvergence;

    auto comparator = []( double res, double tol, int iter, int maxiter ){ if( iter < maxiter ) return true; else return false; };

    double res =  calculateResidual(dataValue, fval, method, h, N);
    resConvergence.push_back( make_pair( res, totalIter ) );

    int currlevel = maxlevels;

    while( res > tol ) {

        auto resHistory = solveGSPoisson(refgrid, dataValue, fval, w, tol, niter, comparator); 
        totalIter += niter;
        multigridVCycle(refgrid, dataValue, fval, currlevel, w, tol, method, totalIter, niter);

        res = calculateResidual(dataValue, fval, method, h, N);
        resConvergence.push_back( make_pair( res, totalIter ) );
        //cout << res << "\n";

    }

    return resConvergence;

}

Solver::Solver() {

    funcMap.insert( std::unordered_map<std::string, mytype>::value_type( "GS", &Solver::solveGSPoisson ) );
    funcMap.insert( std::unordered_map<std::string, mytype>::value_type( "Jacobi", &Solver::solveJacobiPoisson ) );

}