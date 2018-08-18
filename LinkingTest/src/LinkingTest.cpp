//============================================================================
// Name        : LinkingTest.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <gsl/gsl_sf_bessel.h>
#include <cmath>

using namespace std;

int main()
{
	for(double i = 0.0; i < 100.0; i += 1)
	{
		cout << "i = " << i << scientific << "\t " << gsl_sf_bessel_i0_scaled(i) << endl;
	}

	return 0;
}
