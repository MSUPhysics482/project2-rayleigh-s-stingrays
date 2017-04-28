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

const double h = 6.626E-34;           // planck's constant
const double c = 2.998E8;             // speed of light
const double T = 5778;                // absolute temperature of the Sun
const double T0 = 288.15;             // standard temperature of atmosphere at sea level
const double kb = 1.381E-23;          // boltzmann constant
const double b = 2.898E-3;            // wien's displacement constant
const double R = 8.3145;              // universal gas constant
const double p0 = 101325;             // atmospheric pressure
const double n0 = p0*pow(kb*T0,-1.0); // loschmidt's number
const double g = 9.80665;             // acceleration due to gravity

// molar masses:
const double M_N2 = 0.028013; 
const double M_O2 = 0.031999; 
const double M_Ar = 0.039948; 

// proportionality constants for scattering cross-section: 
const double a_N2 = 4.42299E-56;
const double a_O2 = 3.95023E-56;
const double a_Ar = 3.79322E-56;

// generates random number between 0 and 1
double random01(mt19937 generator){
	static uniform_01<mt19937> dist(generator);
	return dist();
}

// spectral radiance of sun as a function of wavelength
// proportionality factors are omitted
double radiance(double wavelength){
	return pow(wavelength,-5.0)*pow(exp(h*c/(wavelength*kb*T))-1,-1.0);
}

// rayleigh scattering cross-section 
double cross_section(double a, double wavelength){
	return a*pow(wavelength,-4.0);
}

// number of gas particles in the atmosphere above a unit area
// i.e. integral of n(h)dh from 0 to inf, where n(h)=n0*exp(-Mgh/(RT0))
// is the particle density as a function of h=height
double part_density(double M){
	return n0*R*T0/(M*g);
}

// probability of light with a given wavelength is scattered
// while traveling through a given gas
double probability(double wavelength, double a, double M){
	return part_density(M)*cross_section(a, wavelength);
}


int main(int argc, char *argv[]) {

    // spectral radiance, maximum spectral radiance
	double y,ymax=0;

	//visible wavelength range
	double xmin=400E-9, xmax=700E-9, xdiff=xmax-xmin, xmid; 

	// fraction of scattered light
	double f,f_N2,f_O2,f_Ar;

	// scaling factor for intensity function
	double scale=4.2E-26;

	// set up for histogram
	int ntot=1E6, nbins=100, ibin;
	double before[nbins]={0.0}, after[nbins]={0.0}, binsize=xdiff/nbins;
	double scattered[nbins]={0.0}, unscattered[nbins]={0.0};

	// give names to random numbers 
	double rx,ry;

	// set up random number generator
	mt19937 generator(time(0));

	// use wien's displacement law to find maximum spectral radiance
	ymax=radiance(b/T);

	// open a data file
	ofstream file;
	file.open("data.dat"); 

	// generate ntot random wavelengths proportional to Planck's distribution
	for(int n=0; n<ntot; n++){

		BEGINNING:

		// throw a random wavelength in visible range
		rx=xmin+xdiff*random01(generator);

		// throw a random spectral radiance value from 0 to ymax
		ry=ymax*random01(generator);

		// calculate the theoretical spectral radiance at rx
		y=radiance(rx);

		// reject rx if ry is larger than theoretical value y
		if(ry>y) goto BEGINNING;

		// keep rx if ry is less than theoretical value y
		// add 1 to bin corresponding to rx
		else{
			ibin=lrint(floor((rx-xmin)/binsize));
			before[ibin]+=1.0;
		}

	}

	// calculate fraction of light that is scattered by atmosphere
	// apply intensity change to scattered portion of light only
	for(ibin=0; ibin<nbins; ibin++){

		// calculate wavelength that corresponds to the center of each bin
		xmid=xmin+(ibin+0.5)*binsize;

		// calculate fraction of light that is scattered by each gas that makes up 
		// the atmosphere and scale by relative percentages of composition
		f_N2=0.78*probability(xmid,a_N2,M_N2);
		f_O2=0.21*probability(xmid,a_O2,M_O2);
		f_Ar=0.01*probability(xmid,a_Ar,M_Ar);

		// total fraction of light scattered
		f=f_N2+f_O2+f_Ar;

		// unscattered light
		unscattered[ibin]=(1.0-f)*before[ibin];

		// scattered light
		scattered[ibin]=f*before[ibin];

		// sun's wavelength spectrum after rayleigh scattering
		// apply intensity change to scattered portion of light only
		after[ibin]=(1.0-f)*before[ibin]+f*scale*pow(xmid,-4.0)*before[ibin];

		// print everything to file
		file << xmid*10E8 << "\t" << before[ibin] << "\t" << after[ibin] 
		<< "\t" << unscattered[ibin] << "\t" << scattered[ibin] << "\n";
	}

	file.close();	

	return 0;
}

