#ifndef utilsdef
#define utilsdef

#include <vector>
#include<Grid.hpp>
#include<ConstantData.hpp>
#include "interpolation.h"
#include <string>
#include <fstream>
#include <unordered_set>

using namespace alglib;
using namespace std;

struct pair_hash {
    inline std::size_t operator()(const std::pair<int,int> & v) const {
        return v.first*31+v.second;
    }
};

using setType = unordered_set< pair<int, int>, pair_hash >;

void buildInterpolant( vector<double>& zeta, vector<double>& eta, vector<double>& fcn, spline2dinterpolant& interpolator );

double simpsons( vector<double>& dataValue, Grid& grid, ConstantData<double>& cvals );

double calculateL2Error( vector<double>& fcn, vector<double>& refFcn, Grid& grid, Grid& refgrid );

template<typename T>
void writeOutput( vector<T>& vals, string filename ) {

    std::ofstream outFile( filename );
    // the important part
    for (const auto &e : vals) outFile << e << "\n";

    outFile.close();

}

void writeOutputPair( vector< pair<double, int> >& vals, string filename );

template<class Fn, typename... Ts>
void writeSoln( vector<double>& soln, Grid& grid, string filename, Fn fn, Ts&... cvals ) {

    ofstream outfile( filename );

    auto N = grid.getNeta();

    for( int j = 0; j < N; j++ ) {

        for( int i = 0; i < N; i++ ) {

            auto xyval = (grid.*fn)( i, j, cvals... );

            outfile << xyval.first << "\t" << xyval.second << "\t" << soln[ j*N + i ] << "\n";

        }

    }

    outfile.close();

}

double calculateLinfError( vector<double>& fcn, vector<double>& refFcn, Grid& grid, Grid& refgrid );

double getF( const setType& sourcePos, int i, int j, int N );

double calculateResidual( std::vector<double> &dataValue, vector<double>& fval, string method, double h, int N );

vector<double> getFval( const setType& sourcePos, int N );

vector<double> calculateResidualArray( std::vector<double> &dataValue, vector<double>& fval, string method, double h, int N );

void addError( vector<double>& dataValue, vector<double>& interpolatedError, int N );

#endif