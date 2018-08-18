/*
 * dnde.cpp
 *
 *  Created on: Feb 7, 2012
 *      Author: danielcumberbatch
 */

#include "dnde.h"
#include "constants.h"
#include "utility_funcs_model_A.h"
#include <vector>
#include <cmath>
#include "tauD.h"
#include "spectrum.h"
#include "halo.h"

using std::vector; using std::log; using std::exp;

/**
 * Returns the number density of positrons per energy interval at specified energy, for a specified epsilon, WIMP positron spectrum per annihilation, dphide,
 * halo model, and table, tauD, of deltav vs. tauD(r, z, deltav).
 *
 * @param epsilon - scaled positron energy, E/GeV.
 *
 * @param ind - vector<double> containing coefficients of the positron spectrum, dphi/dE, per WIMP annihilation, in GeV^-1.
 *
 * @param hm - halo model flag.
 *
 * @param tauDtable - 2-column table containing values of deltav (column 1) and tauD(r,z,deltav) (column 2) for specified r, z.
 *
 * @returns the number density, in cm^-3, of positrons per energy interval at specified location and energy.
 *
 */
double dnde(double epsilon, const vector<double>& ind, int hm, const vector<vector<double> >& tauDtable)
{
    double v_min = e2v(Con::M_X); //minimum v'
    double v_max = e2v(epsilon); //maximum v'
    double lnv_min = log(v_min); //minimum ln(v')
    double lnv_max = log(v_max); //maximum ln(v')

    double integrand = 0.0; //initialize integrand

    for(double lnv = lnv_min; lnv < lnv_max; lnv += Con::DLNV)
    {
        double v = exp(lnv); //v'
        double deltav = v_max - v; //definition of deltav

        integrand += v * w(v2e(v)) * //weight factor
                     dphide(ind, v2e(v)) * //spectrum
                     tauDinterpol(deltav, tauDtable); //tauD, interpolated to v'
    }

    return pow(n0(hm), 2.0) * //normalising WIMP number density
    		Con::SIGMA_V * //annihilation cross-section
            pow(epsilon, -2.0) * integrand * Con::DLNV;
}



