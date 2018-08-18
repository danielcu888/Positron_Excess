/*
 * flux.cpp
 *
 *  Created on: Feb 15, 2012
 *      Author: danielcumberbatch
 */

#include "flux.h"
#include "constants.h"
#include "dnde.h"
#include <vector>

using std::vector;

/**
 * Returns the flux density of positrons per energy interval at specified energy, for a specified epsilon, WIMP positron spectrum per annihilation, dphide,
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
 * @returns the flux density, in cm^-2 s^-1 sr^-1, of positrons per energy interval at specified location and energy.
 *
 */
double flux(double epsilon, const vector<double>& ind, int hm, const vector<vector<double> >& tauDtable)
{
	return dnde(epsilon, ind, hm, tauDtable)* Con::C / 4.0 / Con::PI;
}



