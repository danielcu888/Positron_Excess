//============================================================================
// Name        : main.cpp
// Author      : D. T. Cumberbatch
// Version     : 1.0
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "constants.h"
#include "flux.h"
#include "tauD.h"
#include <iostream>
#include <vector>
#include <iomanip>

using std::vector; using std::scientific;
using std::cout; using std::endl;

int main()
{
    vector<vector<double> > tab = tauDTab(Con::R_LOC, Con::Z_LOC, Con::DELTAV_IND_MIN, Con::DELTAV_IND_MAX, Con::BETA_MIN, Con::BETA_MAX, Con::HM); //construct tauD table for desired location r, z

    //construct spectrum
    double indices[] = {Con::SPEC_IND_0, Con::SPEC_IND_1, Con::SPEC_IND_2, Con::SPEC_IND_3};
    vector<double> spec(indices, indices + sizeof(indices)/sizeof(*indices));

    /*for(double deltav_ind = -2; deltav_ind <= 1.5; deltav_ind += 0.1)
    {
        double deltav = pow(10.0, deltav_ind);
        cout << scientific << deltav << "\t"
             //<< tauD(r, z, deltav, beta_min, beta_max, 0) << "\t"
             << tauD(r, z, deltav, beta_min, beta_max, 1)
             << endl;
    }*/
    /*
    for(vector<vector<double> >::size_type i = 0; i != tab.size(); ++i)
        {
            for(vector<double>::size_type j = 0; j != tab[i].size(); ++j)
                cout << scientific << tab[i][j] << "\t";
            cout << endl;
        }
    */

    //display dnde vs e
    for(int i = 1; i != 100; ++i)
    {
        double epsilon = (double)i * Con::M_X / 100.0;
        cout << scientific
             << "E = " << epsilon << " GeV, \t"
             << "dnde(r,z,e) = " << flux(epsilon, spec, Con::HM, tab) << " GeV^-1 cm^-2 s^-1 sr^-1"
             << endl;
    }

    return 0;
}
