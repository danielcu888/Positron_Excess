/*
 * utility_funcs_model_A.cpp
 *
 *  Created on: Feb 7, 2012
 *      Author: danielcumberbatch
 */

#include "utility_funcs_model_A.h"
#include "constants.h"
#include <cmath>

using std::pow;

/**
 * Returns e = E/1GeV corresponding to specified v for Baltz & Edsjo model A
 *
 * @param v - specified v
 *
 * @return v(e)
 */
double v2e(double v)
{
    double mantissa = (1 - Con::ALPHA) * v;
    double exponent = -1.0/(1.0 - Con::ALPHA);
    return pow(mantissa, exponent);
}

/**
 * Returns v corresponding to specified epsilon = E/1GeV for Baltz & Edsjo model A
 *
 * @param epsilon - specified epsilon
 *
 * @return v(epsilon)
 */
double e2v(double epsilon)
{
    return pow(epsilon, Con::ALPHA - 1.0)/(1.0 - Con::ALPHA);
}

/**
 * Returns weight factor for specified value of epsilon = E/1GeV for model A of Baltz and Edsjo
 *
 * @param epsilon - specified epsilon
 *
 * @return weight factor for specified e, for maodel A of Baltz and Edsjo
 */
double w(double epsilon)
{
    double mantissa = (1.0 - Con::ALPHA) * e2v(epsilon);
    double exponent = -1.0 * (2.0 - Con::ALPHA)/(1.0 - Con::ALPHA);
    return pow(mantissa, exponent);
}
