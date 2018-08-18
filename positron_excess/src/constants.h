/*
 * constants.h
 *
 *  Created on: Jan 11, 2012
 *      Author: danielcumberbatch
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <cmath>

namespace Con
{
	//Integration constants
	const int N_MAX = 10; //max number of image charges
	const double DELTAV_IND_MIN = -5.0;
	const double DELTAV_IND_MAX = 1.5;
	const double DDELTAV_IND = 0.1;
	const double BETA_MIN = -10.0;
	const double BETA_MAX = 3.0;
    const double DLNV = 0.1; //increment in ln(v') for v' integral in double dnde()

	//Conversion Factors
    const double KPC = 3.08568025 * std::pow(10.0, 21.0);

    //Diffusion parameters
    const double K_0 = 3.0 * std::pow(10.0, 27.0) / std::pow(Con::KPC, 2.0); //spatial diffusion constant in kpc^2 s^-1
    const double TAU_E = 1.0 * std::pow(10.0, 16.0); //energy loss time scale, in s
    const double L = 4.0; //diffusion zone width, in kpc
    const double ALPHA = 0.6; //spatial diffusion energy dependence index

    //DM Halo parameters
    const double R_LOC = 8.5; //local galacto-centric radius, in kpc
    const double Z_LOC = 0.0; //local azimuthal distance, in kpc
    const double RHO_LOC = 0.3; //DM denisty at R_LOC, in GeV cm^-3
    const int HM = 1; //Halo Model Flag

    const double A_ISO = 5.0; //scale radius for isothermal DM halo, in kpc
    const double RHO_0_ISO = RHO_LOC * (std::pow(R_LOC/A_ISO, 2.0) + 1.0); //scale density for isothermal DM halo model, in GeV cm^-3
    const double A_NFW = 25.0; //scale radius for NFW DM halo, in kpc
    const double RHO_0_NFW = RHO_LOC * R_LOC * std::pow((R_LOC + A_NFW), 2.0)/std::pow(A_NFW, 3.0); //scale density for isothermal DM halo model, in GeV cm^-3

    //WIMP parameters
    const double M_X = 50; //WIMP mass, in GeV
    const double SIGMA_V = 2.7 * std::pow(10.0, -26.0); //thermally averaged product of WIMP annihilation cross-section and relative velocity, in cm^3 s^-1
    const double SPEC_IND_0 = 0.38183;
    const double SPEC_IND_1 = -1.5962;
    const double SPEC_IND_2 = -0.96026;
    const double SPEC_IND_3 = 0.27938;

    //Scientific Constants
    const double C = 3.0 * pow(10.0, 10.0); //speed of light in a vacuum, in cm s^-1
    const double PI = 3.141592654;
}

#endif /* CONSTANTS_H_ */
