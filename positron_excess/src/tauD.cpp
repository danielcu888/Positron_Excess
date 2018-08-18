/*
 * tauD.cpp
 *
 *  Created on: Feb 7, 2012
 *      Author: danielcumberbatch
 */

#include "tauD.h"
#include "constants.h"
#include <vector>
#include <cmath>
#include <stdexcept>
#include "integration.h"

using std::vector; using std::pow;
using std::log; using std::out_of_range;

/**
 * Returns the value of I[r,z,deltav], for specified Galacto-centric radius, specified azimuthal
 * distance from the Galactic disc and specified deltav, as defined by Eq.(25) in Baltz & Edsjo.
 *
 * @param r - specified Galacto-centric radius, in kpc
 *
 * @param z - specified azimuthal distance, in kpc.
 *
 * @param deltav - specified value of v - v'
 *
 * @return the value of I[r,z,deltav].
 *
 */
double I(double r, double z, double deltav, double beta_min, double beta_max, int hm)
{
    return r_integral(r, deltav, beta_min, beta_max, hm) * z_integral(z, deltav)/ 4.0 / Con::K_0 / Con::TAU_E / deltav;
}

/**
 * Returns the value of tau_D[r,z,deltav], for specified Galacto-centric radius, specified azimuthal
 * distance from the Galactic disc and specified deltav, as defined by Eq.(25) in Baltz & Edsjo.
 *
 * @param r - specified Galacto-centric radius, in kpc
 *
 * @param z - specified azimuthal distance, in kpc
 *
 * @param beta_min - specified minimum power of dummy Galacto-centric radius, given by r_min' = r*pow(10.0, beta_min)
 *
 * @param beta_max - specified maximum power of dummy Galacto-centric radius, given by r_max' = r*pow(10.0, beta_max)
 *
 * @param deltav - specified value of v - v'
 *
 * @param hm - halo model flag
 *
 * @return the value of tau_D[r,z,deltav], in s.
 */
double tauD(double r, double z, double deltav, double beta_min, double beta_max, int hm)
{
    return Con::TAU_E * I(r, z, deltav, beta_min, beta_max, hm);
}

/**
 * Returns a 2-column table containing the values of deltav (column 1) vs. tauD(r, z, deltav) (column 2) for specified
 * range of the dummy Galactocentric radius r (for dummy r-integration), and azimuthal distance from the Galactic disc, z,
 * and specified range of log10(deltav).
 *
 * @param r - specified Galacto-centric radius, in kpc
 *
 * @param z - specified azimuthal distance, in kpc
 *
 * @param deltav_ind_min - specified minimum log10(deltav) of the deltav stored in table, i.e., ret[0][0] = deltav_min = 10^deltav_ind_min
 *
 * @param deltav_ind_max - specified maximum log10(deltav) of the deltav stored in table, i.e., ret[ret.size()-1][0] = deltav_max = 10^deltav_ind_max
 *
 * @param beta_min - specified minimum power of dummy Galacto-centric radius, given by r_min' = r*pow(10.0, beta_min)
 *
 * @param beta_max - specified maximum power of dummy Galacto-centric radius, given by r_max' = r*pow(10.0, beta_max)
 *
 * @param hm - halo model flag
 *
 * @return a 2-column table containing values of deltav (column 1) vs. tauD(r, z, deltav) (column 2).
 *
 */
vector<vector<double> > tauDTab(double r, double z, double deltav_ind_min, double deltav_ind_max, double beta_min, double beta_max, int hm)
{
    vector<vector<double> > ret;
    for(double deltav_ind = deltav_ind_min; deltav_ind <= deltav_ind_max; deltav_ind += Con::DDELTAV_IND)
        {
            double deltav = pow(10.0, deltav_ind);
            vector<double> tmp(2);
            tmp[0] = deltav;
            tmp[1] = tauD(r, z, deltav, beta_min, beta_max, hm);
            ret.push_back(tmp);
        }
    return ret;
}

/**
 * Returns the value of tauD(r, z, deltav) for a specified deltav by linearly interpolating the values in a
 * specified table of tauD vs. deltav.
 *
 * @param deltav - specified value of v - v'
 *
 * @param tauDtable - specified table of deltav vs. tauD
 *
 * @return tauD(deltav)
 */
double tauDinterpol(double deltav, const vector<vector<double> >& tauDtable)
{
    if(deltav < tauDtable[0][0])
        throw out_of_range("Error: deltav too low for interpolation");
    else if(deltav > tauDtable[tauDtable.size() - 1][0])
        throw out_of_range("Error: deltav too large for interpolation");

    vector<vector<double> >::size_type i = 0;
    while(i < tauDtable.size() - 1 && tauDtable[i + 1][0] < deltav)
        ++i;
    if(i == tauDtable.size() - 1)
        throw out_of_range("Error: deltav too large for interpolation");

    double x0 = log(tauDtable[i][0]), x1 = log(tauDtable[i+1][0]);
    double y0 = tauDtable[i][1], y1 = tauDtable[i+1][1];
    double grad = (y1 - y0)/(x1 - x0);
    return y0 + grad * (log(deltav) - x0);
}

