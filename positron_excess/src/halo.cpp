/*
 * halo.cpp
 *
 *  Created on: Feb 7, 2012
 *      Author: danielcumberbatch
 */

#include "constants.h"
#include "halo.h"
#include <stdexcept>
#include <cmath>

using std::invalid_argument; using std::pow;

/**
 * Returns the scale number density of WIMPs, in cm^-3 for the
 * specified halo model.
 *
 * @param hm - halo model flag (0 = ISO, 1 = NFW)
 *
 * @throws invalid_argument
 *
 * @return the scale number density of WIMPs, in cm^-3.
 *
 */
double n0(int hm)
{
    if(hm == 0)
        return Con::RHO_0_ISO/Con::M_X;
    else if(hm == 1)
        return Con::RHO_0_NFW/Con::M_X;
    else
        throw invalid_argument("illegal halo model.");
}

/**
 * Returns the scaled relic density in an NFW virialized halo at specified
 * Galacto-centric radius and azimuthal distance from the Galactic disc.
 *
 * @param r - Galacto-centric radius, in kpc.
 *
 * @param z - Azimuthal distance from Galactic disc, in kpc.
 *
 * @throws invalid_argument
 *
 * @return the scaled relic density in an NFW virialized halo.
 *
 */
double g_NFW(double r, double z)
{
    if(r <= 0 || Con::A_NFW <= 0) throw invalid_argument("g_NFW - illegal arguments");

    double rsph = sqrt(pow(r, 2.0) + pow(z, 2.0));
    return pow(Con::A_NFW, 3.0)/rsph/pow((rsph + Con::A_NFW), 2.0);
}

/**
 * Returns the scaled relic density in an isothermal virialized halo at specified
 * Galacto-centric radius and azimuthal distance from the Galactic disc.
 *
 * @param r - Galacto-centric radius, in kpc.
 *
 * @param z - Azimuthal distance from Galactic disc, in kpc.
 *
 * @throws invalid_argument
 *
 * @return the scaled relic density in an isothermal virialized halo.
 *
 */
double g_ISO(double r, double z)
{
    if(r <= 0 || Con::A_ISO <= 0) throw invalid_argument("g_ISO - illegal arguments");

    return pow(Con::A_ISO, 2.0)/(pow(r, 2.0) + pow(z, 2.0) + pow(Con::A_ISO, 2.0));
}

/**
 * Returns the azimuthal-averaged value of the scaled relic density at specified
 * Galacto-centric radius, for a specified halo model.
 *
 * @param r - Galacto-centric radius, in kpc.
 *
 * @param hm - halo model flag (0 = ISO, 1 = NFW)
 *
 * @return the azimuthal-averaged value of the scaled relic density at specified
 * Galacto-centric radius, for a specified halo model.
 *
 */
double f(double r, int hm)
{
    double integrand = 0.0;
    int num_steps = 100;
    double dz = 2*Con::L/num_steps;

    if(hm == 0)
    {
        for(double z = -Con::L; z < Con::L; z += dz)
            integrand += pow(g_NFW(r, z), 2.0);
    }
    else if(hm == 1)
    {
        for(double z = -Con::L; z < Con::L; z += dz)
            integrand += pow(g_ISO(r, z), 2.0);
    }

    return integrand * dz/2.0/Con::L;
}
