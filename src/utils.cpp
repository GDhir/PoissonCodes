#include "utils.hpp"

double simpsons( vector<double>& dataValue, Grid& grid, ConstantData<double>& cvals ) {

    vector<double> f( dataValue );

    int neta = grid.getNeta();
    int nzeta = grid.getNzeta();

    double spacing = grid.getH();

    for( int j = 0; j < neta; j++ ) {
        for( int i = 0; i < nzeta; i++ ) {

            double eta = static_cast<double>(j)/static_cast<double>( neta - 1 );

            f[ j*nzeta + i ] = f[ j*nzeta + i ]*( 2*cvals.a*eta + cvals.b )*cvals.h/2;

        }

    }

    vector<double> svals( neta, 0 );

    for( int j = 0; j < neta; j++ ) {

        svals[j] += f[ j*nzeta + 0 ] + f[ j*nzeta + nzeta - 1 ];

        for( int i = 1; i < nzeta - 1; i++ ) {

            if( i%2 != 0 ) {

                svals[j] += 4*f[ j*nzeta + i ];

            }
            else {
                
                svals[j] += 2*f[ j*nzeta + i ];

            }

        }
    }

    double sum{0};

    sum += svals[ 0 ] + svals[ neta - 1 ];

    for( int j = 1; j < neta - 1; j++ ) {

        if( j%2 != 0 ) {

            sum += 4*svals[j];

        }
        else {

            sum += 2*svals[j];

        }

    }

    return sum*spacing*spacing/9;

}

void buildInterpolant( vector<double>& zeta, vector<double>& eta, vector<double>& fcn, spline2dinterpolant& interpolator ) {

    int N = static_cast<int>( zeta.size() );

    alglib::real_1d_array zetaArr;
    zetaArr.attach_to_ptr( zeta.size(), zeta.data() );

    alglib::real_1d_array etaArr;
    etaArr.attach_to_ptr( eta.size(), eta.data() );

    alglib::real_1d_array fArr;
    fArr.attach_to_ptr( fcn.size(), fcn.data() );

    spline2dbuildbicubicv(zetaArr, N, etaArr, N, fArr, 1, interpolator);

}

double calculateL2Error( vector<double>& fcn, vector<double>& refFcn, Grid& grid, Grid& refgrid ) {

    spline2dinterpolant interpolator;

    int N = grid.getNeta();
    int Nref = refgrid.getNeta();

    vector<double> eta( N, 0 );

    for( int j = 0; j < N; j++ ) {

        eta[j] = grid.calcXFromIndex( j, N );

    }

    vector<double> zeta( eta );

    buildInterpolant( zeta, eta, fcn, interpolator );

    double error = 0;

    for( int j = 0; j < Nref; j++ ) {

        auto eta = refgrid.calcXFromIndex( j, Nref );

        for( int i = 0; i < Nref; i++ ) {

            auto zeta = refgrid.calcXFromIndex( i, Nref );
            error += pow( refFcn[ j*Nref + i ] - spline2dcalc( interpolator, zeta, eta ), 2 );

        }

    }

    error = sqrt( error );

    return error;

}

double calculateLinfError( vector<double>& fcn, vector<double>& refFcn, Grid& grid, Grid& refgrid ) {

    spline2dinterpolant interpolator;

    int N = grid.getNeta();
    int Nref = refgrid.getNeta();

    vector<double> eta( N, 0 );

    for( int j = 0; j < N; j++ ) {

        eta[j] = grid.calcXFromIndex( j, N );

    }

    vector<double> zeta( eta );

    buildInterpolant( zeta, eta, fcn, interpolator );

    double error = 0;

    for( int j = 0; j < Nref; j++ ) {

        auto etaval = refgrid.calcXFromIndex( j, Nref );

        for( int i = 0; i < Nref; i++ ) {

            auto zetaval = refgrid.calcXFromIndex( i, Nref );

            error = std::max( error, std::abs( refFcn[ j*Nref + i ] - spline2dcalc( interpolator, zetaval, etaval ) ) );

        }

    }

    return error;

}

double getF( const setType& sourcePos, int i, int j, int N ){

    int nPerBlk = ( N - 1 )/6;
    int blkIdx = i/nPerBlk;
    int blkIdy = j/nPerBlk;

    auto blkIndices = make_pair( blkIdx, blkIdy );

    if( sourcePos.find( blkIndices ) != sourcePos.end() ) {
        return 1;
    }
    else {

        if( i%nPerBlk == 0 && blkIdx > 1 && blkIdy > 0 ) {

            if( sourcePos.find( make_pair( blkIdx - 1, blkIdy ) ) != sourcePos.end() ) {
                return 1;
            }

        }

        if( j%nPerBlk == 0 && blkIdy > 1 && blkIdx > 0 ) {

            if( sourcePos.find( make_pair( blkIdx, blkIdy - 1 ) ) != sourcePos.end() ) {
                return 1;
            }

        }

        if( i%nPerBlk == 0 && blkIdx > 1 && j%nPerBlk == 0 && blkIdy > 1 ) {

            if( sourcePos.find( make_pair( blkIdx - 1, blkIdy - 1 ) ) != sourcePos.end() ) {
                return 1;
            }

        }

        return 0;

    }

}

double calculateResidual( std::vector<double> &dataValue, vector<double>& fval, string method, double h, int N ) {

    double res{0};

    if( method == "GS" ) {

            for (int j = 1; j < N - 1; j++)
            {

                for (int i = 1; i < N - 1; i++)
                {

                    res += pow(fval[ j*N + i ] + ( -4*dataValue[j * N + i] + dataValue[j * N + i + 1] + dataValue[j * N + i - 1] +
                                                       dataValue[(j + 1) * N + i] + dataValue[(j - 1) * N + i] )/(h*h),
                               2);
                }
            }

            res = sqrt(res*h);

    }
    else {

            for (int j = 1; j < N - 1; j++)
            {

                for (int i = 1; i < N - 1; i++)
                {

                    res += pow(fval[ j*N + i ] + ( -4*dataValue[j * N + i] + dataValue[j * N + i + 1] + dataValue[j * N + i - 1] +
                                                       dataValue[(j + 1) * N + i] + dataValue[(j - 1) * N + i] )/(h*h),
                               2);
                }
            }

            res = sqrt(res * h);

    }

    return res;

}

vector<double> calculateResidualArray( std::vector<double> &dataValue, vector<double>& fval, string method, double h, int N ) {

    vector<double> res( N*N, 0 );

    if( method == "GS" ) {

            for (int j = 1; j < N - 1; j++)
            {

                for (int i = 1; i < N - 1; i++)
                {

                    res[ j*N + i ] = fval[ j*N + i ] + ( -4*dataValue[j * N + i] + dataValue[j * N + i + 1] + dataValue[j * N + i - 1] +
                                                       dataValue[(j + 1) * N + i] + dataValue[(j - 1) * N + i] )/(h*h);
                }
            }

    }
    else {

            for (int j = 1; j < N - 1; j++)
            {

                for (int i = 1; i < N - 1; i++)
                {

                    res[ j*N + i ] = (dataValue[j * N + i] - (h * h *fval[j*N + i] + dataValue[j * N + i + 1] + dataValue[j * N + i - 1] +
                                                       dataValue[(j + 1) * N + i] + dataValue[(j - 1) * N + i]) /
                                                          4.0 );
                }
            }

    }

    return res;

}

vector<double> getFval( const setType& sourcePos, int N ) {

    vector<double> fval( N*N, 0 );

    for (int j = 1; j < N - 1; j++)
    {

        for (int i = 1; i < N - 1; i++)
        {

            fval[ j*N + i ] = getF(sourcePos, i, j, N);
        }
    }

    return fval;

}

void addError( vector<double>& dataValue, vector<double>& interpolatedError, int N ){

    for( int j = 1; j < N - 1; j++ ) {
        for( int i = 0; i < N - 1; i++ ) {

            dataValue[ j*N + i ] += interpolatedError[ j*N + i ];

        }
    }

}

void writeOutputPair( vector< pair<double, int> >& vals, string filename ) {

    std::ofstream outFile( filename );
    // the important part
    for (const auto &e : vals) outFile << e.first << "\t" << e.second << "\n";

    outFile.close();

}