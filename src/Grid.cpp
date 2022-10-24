#include "Grid.hpp"

int Grid::getNeta()
{
    return this->neta;
}

int Grid::getNzeta()
{
    return this->nzeta;
}

int Grid::getNx()
{
    return this->nzeta;
}

int Grid::getNy()
{
    return this->neta;
}

double Grid::calcXFromIndex(int idx, int npoints)
{

    return static_cast<double>(idx) / static_cast<double>((npoints - 1));
}

Grid::Grid(int _neta, int _nzeta) : neta(_neta), nzeta(_nzeta) {}

double Grid::getH()
{

    if (neta == nzeta)
    {
        return 1 / (static_cast<double>(neta - 1));
    }
    else
        throw("neta != nzeta");
}

std::pair<double, double> Grid::getXY(int i, int j, ConstantData<double> &cvals)
{

    double eta = calcXFromIndex(j, neta);
    double zeta = calcXFromIndex(i, nzeta);

    std::pair<double, double> xyval;

    xyval.first = zeta * (2 * cvals.a * eta + cvals.b) / 2;
    xyval.second = eta * cvals.h;

    return xyval;
}

std::pair<double, double> Grid::getXYRegular(int i, int j)
{

    double y = calcXFromIndex(j, neta);
    double x = calcXFromIndex(i, nzeta);

    return std::make_pair(x, y);
}

vector<double> Grid::restrictVal(vector<double> &fineVals, Grid &coarsegrid)
{

    int coarseN = coarsegrid.getNx();
    int fineN = this->getNx();

    vector<double> coarseVals(coarseN * coarseN, 0);

    for (int j = 1; j < coarseN - 1; j++)
    {

        for (int i = 1; i < coarseN - 1; i++)
        {

            coarseVals[j * coarseN + i] = fineVals[2 * j * fineN + 2 * i] / 4 + fineVals[2 * j * fineN + 2 * i + 1] / 8 + fineVals[2 * j * fineN + 2 * i - 1] / 8 +
                                          fineVals[(2 * j + 1) * fineN + 2 * i] / 8 + fineVals[(2 * j + 1) * fineN + 2 * i + 1] / 16 + fineVals[(2 * j + 1) * fineN + 2 * i - 1] / 16 +
                                          fineVals[(2 * j - 1) * fineN + 2 * i] / 8 + fineVals[(2 * j - 1) * fineN + 2 * i + 1] / 16 + fineVals[(2 * j - 1) * fineN + 2 * i - 1] / 16;
        }
    }

    return coarseVals;
}

vector<double> Grid::interpolateVal(vector<double> &coarseVals, Grid &coarsegrid)
{

    int coarseN = coarsegrid.getNx();
    int fineN = this->getNx();

    vector<double> fineVals(fineN*fineN, 0);

    for (int j = 1; j < fineN - 1; j++)
    {
        for (int i = 1; i < fineN - 1; i++)
        {

            if (i % 2 == 0 && j % 2 == 0)
            {
                fineVals[j * fineN + i] = coarseVals[(j / 2) * coarseN + (i / 2)];
            }
            else if (i % 2 != 0 && j % 2 != 0)
            {
                fineVals[j * fineN + i] = (coarseVals[(j / 2) * coarseN + (i / 2)] + coarseVals[(j / 2) * coarseN + (i / 2) + 1] +
                                           coarseVals[(j / 2 + 1) * coarseN + (i / 2)] + coarseVals[(j / 2 + 1) * coarseN + (i / 2 + 1)]) /
                                          4;
            }
            else if (i % 2 == 0)
            {
                fineVals[j * fineN + i] = (coarseVals[(j / 2) * coarseN + (i / 2)] + coarseVals[(j / 2 + 1) * coarseN + (i / 2)]) / 2;
            }
            else if (j % 2 == 0)
            {
                fineVals[j * fineN + i] = (coarseVals[(j / 2) * coarseN + (i / 2)] + coarseVals[(j / 2) * coarseN + (i / 2) + 1]) / 2;
            }
        }
    }

    return fineVals;
}