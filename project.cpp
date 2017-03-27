#include <cstdlib>
#include <cmath>
#include <fstream>
#include <ctime>
#include <boost/random.hpp>
#include <stdlib.h>
#include <iostream>

using namespace std;
using boost::uniform_01;
using boost::mt19937;

const double pi = 4.0*atan(1.0);
const double h = 6.626E-34;
const double c = 2.998E8;
const double T = 5778;
const double kb = 1.381E-23;

/* random number generator to use later.

double random01(mt19937 generator){
	static uniform_01<mt19937> dist(generator);
	return dist();
}
*/

int main(int argc, char *argv[]) {

	double R,I;
	//double xmin=atof(argv[1]), xmax=atof(argv[2]), dx=10E-9;
	double xmin=400E-9, xmax=700E-9, dx=10E-9; //visible wavelength range
	double f=1E25; //scaling factor

	ofstream file;
	file.open("data.txt");

	for(double x=xmin; x<xmax; x+=dx){

		R=(2.0*h*c*c*pow(x,-5.0))*pow(exp(h*c/(x*kb*T))-1,-1.0);
		I=pow(x,-4.0); //only wavelength dependence for now
		file << x << "\t" << f*R << "\t" << I*R << "\n";
	}

	file.close();

	return 0;


}

/*

TO DO:

- Add dependency on radius of air molecules into intensity function. 
  Randomly scatter off N2 80% of the time and O2 20%. 
  We will have to estimate N2 and O2 atoms to be spheres. 

- Randomly generate wavelengths proportional to spectral distribution 
  (Planck's Law) in the visible wavelength range.

- Write something that will plot this data.

*/