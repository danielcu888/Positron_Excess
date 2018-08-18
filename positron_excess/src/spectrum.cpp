/*
 * spectrum.cpp
 *
 *  Created on: Feb 7, 2012
 *      Author: danielcumberbatch
 */

#include "spectrum.h"
#include <vector>
#include <cmath>

using std::vector; using std::pow;

/**
 * Returns the number of positrons at specified scaled energy, epsilon = E/1GeV,
 * per WIMP annihilation for given parameterised spectrum.
 *
 * @param ind - a vector<double>, v, where v[i] is such that the parameterised
 * spectrum is given by the sum \Sigma_(i)^(v.size()-1) {v[i]*pow(epsilon, (double)i)}
 *
 * @param epsilon - the specified scaled energy
 *
 * @return the number of positrons at specified scaled energy.
 *
 */
double dphide(const vector<double>& ind, double epsilon)
{
    double ret = 0.0;
    for(vector<double>::size_type i = 0; i != ind.size(); ++i)
        ret += ind[i] * pow(epsilon, static_cast<double>(i));
    return ret;
}


