/*
 * integration.cpp
 *
 *  Created on: Feb 7, 2012
 *      Author: danielcumberbatch
 */

#include "integration.h"
#include "constants.h"
#include "halo.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_erf.h>
#include <cmath>

using std::pow; using std::exp;

/**
 * Returns the summation over image charges as specified in
 * Eq.(25) of Baltz & Edsjo, for a specified azimuthal distance
 * from the Galactic disc, and a specified deltav.
 *
 * @param z - specified azimuthal distance from the Galactic disc, in kpc.
 *
 * @param deltav - specified value of v - v'
 *
 * @return the z-summation over image charges.
 *
 */
double z_integral(double z, double deltav)
{
    double integrand = 0.0;
    for(int n = -Con::N_MAX; n <= Con::N_MAX; ++n)
    {
        double arg1 = ((pow(-1.0, n)*Con::L) + (2.0*Con::L*n + z))/sqrt(4.0*Con::K_0*Con::TAU_E*deltav),
                arg2 = ((pow(-1.0, n)*Con::L) + (2.0*Con::L*n - z))/sqrt(4.0*Con::K_0*Con::TAU_E*deltav);
        integrand += (gsl_sf_erf(arg1) + gsl_sf_erf(arg2));
    }
    return integrand;
}

/**
 * Returns the result of the integral over dummy Galacto-centric radii, as specified in
 * Eq.(25) of Baltz & Edsjo, for a specified Galacto-centric radius, specified range of
 * dummy Galacto-centric radii, specified value of deltav, and specified halo model.
 *
 * @param r - specified Galacto-centric radius, in kpc.
 *
 * @param deltav - specified value of v - v'
 *
 * @param beta_min - specified minimum power of dummy Galacto-centric radius, given by r_min' = r*pow(10.0, beta_min)
 *
 * @param beta_max - specified maximum power of dummy Galacto-centric radius, given by r_max' = r*pow(10.0, beta_max)
 *
 * @return the result of the integral over dummy Galacto-centric radii, in kpc^2
 *
 */
double r_integral(double r, double deltav, double beta_min, double beta_max, int hm)
{
    const double c = 2.0*r*r/4.0/Con::K_0/Con::TAU_E;
    double dbeta = 0.01;
    double integrand = 0.0;
    for(double beta = beta_min; beta < beta_max; beta += dbeta)
        {
            double tmp = pow(10.0, 2.0*beta)
                         * f((r*pow(10.0, beta)), hm)
                         * gsl_sf_bessel_i1_scaled(pow(10.0, beta) * c / deltav)
                         * exp(-pow((1-pow(10.0, beta)), 2.0) * c / 2.0 / deltav);
            integrand += tmp;
        }
    return pow(r, 2.0) * log(10.0) * dbeta * integrand;
}


