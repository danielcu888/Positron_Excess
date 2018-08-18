/*
 * tauD.h
 *
 *  Created on: Feb 7, 2012
 *      Author: danielcumberbatch
 */

#ifndef TAUD_H_
#define TAUD_H_

#include <vector>

double I(double, double, double, double, double, int);
double tauD(double, double, double, double, double, int);
std::vector<std::vector<double> > tauDTab(double, double, double, double, double, double, int);
double tauDinterpol(double, const std::vector<std::vector<double> >&);

#endif /* TAUD_H_ */
