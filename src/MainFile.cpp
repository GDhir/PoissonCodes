#include "Grid.hpp"
#include "Solver.hpp"
#include "ConstantData.hpp"
#include "interpolation.h"
#include "utils.hpp"
using namespace alglib;
using namespace std;
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

void computeChannelFlow() {

    double l = 3, h = 1, b = 1.0;
    double tol = pow( 10, -4 );
    vector<int> gridsize{ 41 };
    vector<Grid> grids;

    string prefixstr = "/home/gaurav/plots/plotdir_b=1.0/";   
    writeOutput<int>( gridsize, prefixstr + "gridsize.txt" );

    for( auto& sz: gridsize ) {
        grids.push_back( Grid( sz, sz ) );
    }

    vector<double> QVals;
    vector<double> l2Error;
    vector<double> linfError;

    ConstantData<double> cvals( l, h, b );

    int N = grids[0].getNeta();
    vector<double> refsoln( N*N, 0);
    Solver sjacobi;

    sjacobi.solveChannelFlow( grids[0], refsoln, cvals, tol );
    auto sum = simpsons( refsoln, grids[0], cvals );
    QVals.push_back( sum );
    l2Error.push_back(0);
    linfError.push_back(0);

    using funcType = std::pair<double, double>(Grid::*)(int, int, ConstantData<double>& );
    funcType fn = &Grid::getXY;

    //cout << sum << "\n";
    writeSoln( refsoln, grids[0], prefixstr + "solnsN=" + to_string( N ) + ".txt", fn, cvals );

    for( int i = 1; i < static_cast<int>( grids.size() ); i++ ) {

        auto grid = grids[i];

        int N = grid.getNeta();

        vector<double> dataValue( N*N, 0 );
        
        sjacobi.solveChannelFlow( grid, dataValue, cvals, tol );

        auto sum = simpsons( dataValue, grid, cvals );
        QVals.push_back( sum );

        l2Error.push_back( calculateL2Error( dataValue, refsoln, grids[i], grids[0] ) );
        linfError.push_back( calculateLinfError( dataValue, refsoln, grids[i], grids[0] ) );

        writeSoln( dataValue, grids[i], prefixstr + "solnsN=" + to_string( N ) + ".txt", fn, cvals );

        //std::cout << sum << "\t" << l2Error[i] << "\t" << linfError[i] << "\n";

    }

    writeOutput<double>( QVals, prefixstr + "flowrate.txt" );
    writeOutput<double>( l2Error, prefixstr + "l2Error.txt" );
    writeOutput<double>( linfError, prefixstr + "linfError.txt" );
    writeOutput<int>( gridsize, prefixstr + "gridsize.txt" );

}

void computePareto() {

    double l = 3;
    vector<double> hvals, bvals;
    double tol = pow( 10, -6 );
    Grid grid( 41, 41 );

    for( double k = 0.1; k <= 1; k+= 0.1  ) {
        hvals.push_back( k );
    }

    for( double k = 0.0; k <= 1; k+= 0.1  ) {
        bvals.push_back( k );
    }

    string prefixstr = "/home/gaurav/plots/pareto/";

    vector<vector<double>> QVals;

    using funcType = std::pair<double, double>(Grid::*)(int, int, ConstantData<double>& );
    funcType fn = &Grid::getXY;
    Solver sjacobi;
    ofstream qvalpareto(prefixstr + "qvalpareto.txt");

    for( auto& h: hvals) {
        for( auto& b: bvals ) {

            int N = grid.getNeta();
            ConstantData<double> cvals( l, h, b );

            vector<double> dataValue( N*N, 0 );
            
            sjacobi.solveChannelFlow( grid, dataValue, cvals, tol );

            auto sum = simpsons( dataValue, grid, cvals );
            qvalpareto << h << "\t" << b << "\t" << sum << "\n";

            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << b;
            std::string bstr = ss.str();

            std::stringstream hs;
            hs << std::fixed << std::setprecision(2) << h;
            std::string hstr = hs.str();

            writeSoln( dataValue, grid, prefixstr + "solnsN=" + to_string( N ) + "b=" + bstr + "h=" + hstr + ".txt", fn, cvals );

        }

    }

}

void computePoissonIterative( string method ) {

    vector<double> wvals;

    // for( double val = 0.1; val < 2; val += 0.01 ) {

    //     wvals.push_back(val);

    // }
    wvals.push_back( 0.5 );
    wvals.push_back( 0.8 );
    wvals.push_back(1);

    auto comparator = []( double res, double tol, int iter, int maxiter ){ if( res > tol ) return true; else return false; };

    for( auto&w: wvals ) {

        int N{49}, maxiter{1000000};
        Grid grid( N, N );

        string prefixstr = "/home/gaurav/plots/PoissonMG/" + method + "/";

        Solver slv;

        vector<double> dataValue( N*N, 0 );

        setType sourcePos( { make_pair( 2, 1 ), make_pair(2, 3), make_pair( 3, 3 ), make_pair( 4, 1 ) } );
        auto fval = getFval( sourcePos, N );

        double tol = pow( 10, -6 );

        auto slvmethod = slv.funcMap[method];

        auto residual = (slv.*slvmethod)( grid, dataValue, fval, w, tol, maxiter, comparator );

        std::stringstream ss;
        ss << std::fixed << std::setprecision(2) << w;
        std::string wstr = ss.str();

        writeOutputPair( residual, prefixstr + "residualN=" + to_string( N ) + "w=" + wstr + ".txt" );

        using funcType = std::pair<double, double>(Grid::*)(int, int );
        funcType fn = &Grid::getXYRegular;

        writeSoln( dataValue, grid, prefixstr + "solnsN=" + to_string( N ) + "w=" + wstr + ".txt", fn );
    
    }

}

void computePoissonMG( int Nval, int nlevels ) {

    vector<double> wvals;

    wvals.push_back( 0.5 );
    wvals.push_back( 0.8 );
    wvals.push_back( 1.0 );

    int N{Nval}, maxlevels{nlevels};
    Grid grid( N, N );
    string method = "GS";
    vector<pair<double, int>> resConvergence;

    for( auto&w: wvals ) {
        cout << w;

        string prefixstr = "/home/gaurav/plots/PoissonMG/MG/";

        Solver slv;

        vector<double> dataValue( N*N, 0 );

        setType sourcePos( { make_pair( 2, 1 ), make_pair(2, 3), make_pair( 3, 3 ), make_pair( 4, 1 ) } );
        auto fval = getFval( sourcePos, N );

        double tol = pow( 10, -6 );

        resConvergence = slv.computeMultigrid( grid, dataValue, fval, maxlevels, w, tol, method );

        std::stringstream ss;
        ss << std::fixed << std::setprecision(2) << w;
        std::string wstr = ss.str();

        writeOutputPair( resConvergence, prefixstr + "residualN=" + to_string( N ) + "w=" + wstr + "maxlevels=" + to_string(maxlevels) + ".txt" );

        using funcType = std::pair<double, double>(Grid::*)(int, int );
        funcType fn = &Grid::getXYRegular;

        writeSoln( dataValue, grid, prefixstr + "solnsN=" + to_string( N ) + "w=" + wstr + "maxlevels=" + to_string(maxlevels) + ".txt", fn );
    
    }

}

void getAllSources( int i, int j, vector<setType>& sourceArr, setType& currArr ) {

    if( currArr.size() == 4 ) {
        sourceArr.push_back( currArr );
        return;
    }
    else {

        currArr.insert( make_pair( i, j ) );

        if( i + 1 <= 4 ) {

            getAllSources( i + 1, j, sourceArr, currArr );

        }
        else if( j + 1 <= 4 ) {

            getAllSources( 1, j + 1, sourceArr, currArr );

        }
        else {

            if( currArr.size() == 4 ) {
                sourceArr.push_back( currArr );
            }

        }

        currArr.erase( make_pair( i, j ) );

        if( i + 1 <= 4 ) {

            getAllSources( i + 1, j, sourceArr, currArr );

        }
        else if( j + 1 <= 4 ) {

            getAllSources( 1, j + 1, sourceArr, currArr );

        }

    }

}

double getNormalError( vector<double>& dataValue, vector<vector<double>>& actnormalDerivative, Grid& grid ) {

    double error{0}, derivative{0};

    int N = grid.getNx();
    double h = grid.getH();

    for( int i = 1; i < N - 1; i++ ) {

        derivative = ( -3*dataValue[ 0*N + i ] + 4*dataValue[ 1*N + i ] - dataValue[ 2*N + i ] )/( 2*h );

        error += abs( derivative - actnormalDerivative[0][i]);

        derivative = ( 3*dataValue[ (N - 1)*N + i ] - 4*dataValue[ (N - 2)*N + i ] + dataValue[ (N - 3)*N + i ] )/( 2*h );

        error += abs( derivative - actnormalDerivative[1][i] );

    }

    for( int j = 1; j < N - 1; j++ ) {

        derivative = ( -3*dataValue[ j*N + 0 ] + 4*dataValue[ j*N + 1 ] - dataValue[ j*N + 2 ] )/( 2*h );

        error += abs( derivative - actnormalDerivative[2][j] );

        derivative = ( 3*dataValue[ j*N + N - 1 ] - 4*dataValue[ j*N + N - 2 ] + dataValue[ j*N + N - 3 ] )/( 2*h );

        error += abs( derivative - actnormalDerivative[3][j] );

    }

    return error;
}

double getAlternateNormalError( vector<double>& dataValue, vector<vector<double>>& actnormalDerivative, Grid& grid ) {

    double error{0}, derivative{0};

    int N = grid.getNx();
    double h = grid.getH();

    for( int i = 1; i < N - 1; i++ ) {

        derivative = ( 4*dataValue[ 0*N + i ] - 2*dataValue[ 1*N + i ] - dataValue[ 0*N + i + 1 ] - dataValue[ 0*N + i - 1 ] )/( 2*h );

        error += abs( derivative - actnormalDerivative[0][i]);

        derivative = ( 4*dataValue[ (N - 1)*N + i ] - 2*dataValue[ (N - 2)*N + i ] - dataValue[ (N - 1)*N + i + 1 ] - dataValue[ (N - 1)*N + i - 1 ] )/( 2*h );

        error += abs( derivative - actnormalDerivative[1][i]);

    }

    for( int j = 1; j < N - 1; j++ ) {

        derivative = ( 4*dataValue[ j*N + 0 ] - 2*dataValue[ j*N + 1 ] - dataValue[ (j + 1)*N + 0 ] - dataValue[ (j - 1)*N + 0 ] )/( 2*h );

        error += abs( derivative - actnormalDerivative[2][j] );

        derivative = ( 4*dataValue[ j*N + N - 1 ] - 2*dataValue[ j*N + N - 2 ] - dataValue[ (j + 1)*N + N - 1 ] - dataValue[ (j - 1)*N + N - 1 ] )/( 2*h );

        error += abs( derivative - actnormalDerivative[3][j] );

    }

    return error;
}

void findSources() {

    double w = 1;
    int N = 25, maxlevels = 2, maxiter = 1000000;
    Grid grid( N, N );
    string method = "GS";

    vector<setType> sourceArr;
    setType currArr;

    getAllSources( 1, 1, sourceArr, currArr );

    Solver slv;

    setType bestSourcePos;
    double bestError{1e5};

    vector<vector<double>> normalDerivative;

    normalDerivative.push_back( vector<double>{ 0.0000, 0.0076, 0.0154, 0.0233, 0.0314, 0.0398, 0.0484, 0.0565, 0.0635, 0.0686, 0.0709, 0.0703, 0.0670, 0.0619, 0.0556, 0.0492, 0.0429, 0.0369, 0.0313, 0.0259, 0.0206, 0.0154, 0.0103, 0.0051, 0.0000 } );
    normalDerivative.push_back( vector<double>{ 0.0000, 0.0099, 0.0198, 0.0297, 0.0396, 0.0494, 0.0588, 0.0676, 0.0748, 0.0800, 0.0823, 0.0816, 0.0784, 0.0734, 0.0674, 0.0611, 0.0548, 0.0486, 0.0424, 0.0360, 0.0292, 0.0222, 0.0150, 0.0075, 0.0000 } );
    normalDerivative.push_back( vector<double>{ 0.0000, 0.0075, 0.0150, 0.0222, 0.0292, 0.0360, 0.0424, 0.0486, 0.0548, 0.0611, 0.0674, 0.0734, 0.0784, 0.0816, 0.0823, 0.0800, 0.0748, 0.0676, 0.0588, 0.0494, 0.0396, 0.0297, 0.0198, 0.0099, 0.0000 } );
    normalDerivative.push_back( vector<double>{ 0.0000, 0.0051, 0.0103, 0.0154, 0.0206, 0.0259, 0.0313, 0.0369, 0.0429, 0.0492, 0.0556, 0.0619, 0.0670, 0.0703, 0.0709, 0.0686, 0.0635, 0.0565, 0.0484, 0.0398, 0.0314, 0.0233, 0.0154, 0.0076, 0.0000 } );
    auto comparator = []( double res, double tol, int iter, int maxiter ){ if( res > tol ) return true; else return false; };

    for( auto& sourcePos: sourceArr ) {

        vector<double> dataValue( N*N, 0 );

        auto fval = getFval( sourcePos, N );

        double tol = pow( 10, -6 );

        slv.computeMultigrid( grid, dataValue, fval, maxlevels, w, tol, method );
        //slv.solveGSPoisson( grid, dataValue, fval, w, tol, maxiter, comparator );

        auto error = getNormalError( dataValue, normalDerivative, grid );

        cout << error << "\n";

        if( error < bestError ) {
            bestError = error;
            bestSourcePos = sourcePos;
        }

    }

    ofstream dfval( "sourcePos.txt" );
    for( auto& s: bestSourcePos ) {

        dfval << s.first + (s.second - 1)*4 << "\t";

    }

    dfval << "\n" << bestError << "\n";


}

int main() {

    cout << "a" << "\n";
    //computePoissonMG( 97, 2 );
    //computePoissonMG( 97, 3 );
    //computePoissonMG( 97, 4 );
    //computePoissonMG(97);
    //computePoissonIterative( "GS" );
    //computePoissonIterative( "GS" );
    computePareto();
    vector<setType> sourceArr;
    setType currArr;

    //getAllSources( 1, 1, sourceArr, currArr );

    //findSources();

    // ofstream dfval( "sources.txt" );

    // for( auto& elem: sourceArr ) {

    //     for( auto& s: elem ) {

    //         dfval << s.first + (s.second - 1)*4 << "\t";

    //     }

    //     dfval << "\n";

    // }

    // dfval << sourceArr.size() << "\n";

    return 0;
}