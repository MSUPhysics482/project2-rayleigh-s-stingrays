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
const double h = 6.626E-34; // planck constant
const double c = 2.998E8; // speed of light
const double T = 5778; // absolute temperature of the Sun
const double T0 = 288.15; // standard temperature of atmosphere at sea level
const double kb = 1.381E-23; // boltzmann constant
const double b = 2.898E-3; // wien's displacement constant
const double R = 8.3145; // universal gas constant
const double p0 = 101325; // atmospheric pressure
const double n0 = p0*pow(kb*T0,-1.0); // particle density of ideal gas at p0 and T0
const double g = 9.80665; // acceleration due to gravity

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

double cross_section(double a, double wavelength){
	return a*pow(wavelength,-4.0);
}

double part_density(double M){
	return n0*R*T0/(M*g);
}

// fraction of light with given wavelength that is scattered
// while traveling through a given gas
double probability(double wavelength, double a, double M){
	return part_density(M)*cross_section(a, wavelength);
}


int main(int argc, char *argv[]) {

    // spectral radiance, maximum spectral radiance, intensity;
	double y,ymax=0,I;
	double f,f_N2,f_O2,f_Ar;

	//visible wavelength range
	double xmin=400E-9, xmax=700E-9, xdiff=xmax-xmin, xmid; 

	// give names to random numbers 
	double rx,ry;

	// set up for histogram
	int ntot=1E6, nbins=100, ibin;
	double before[nbins]={0.0}, after[nbins]={0.0}, binsize=xdiff/nbins;

	// set up random number generator
	mt19937 generator(time(0));

	// use wien's displacement law to find maximum spectral radiance
	ymax=radiance(b/T);

	// open a data file
	ofstream file;
	file.open("data.txt"); 

	// generate ntot random wavelengths proportional to Planck's distribution
	for(int n=0; n<ntot; n++){

		BEGINNING:

		// throw a random wavelength in visible range
		rx=xmin+xdiff*(random01(generator));

		// throw a random spectral radiance value from 0 to ymax
		ry=ymax*(random01(generator));

		// calculate the theoretical spectral radiance at rx
		y=radiance(rx);

		if(ry>y) goto BEGINNING;
		else{
			ibin=lrint(floor((rx-xmin)/binsize));
			before[ibin]+=1.0;
		}

	}

	for(ibin=0; ibin<nbins; ibin++){
		xmid=xmin+(ibin+0.5)*binsize;
		f_N2=0.78*probability(xmid,a_N2,M_N2);
		f_O2=0.21*probability(xmid,a_O2,M_O2);
		f_Ar=0.01*probability(xmid,a_Ar,M_Ar);
		f=f_N2+f_O2+f_Ar;
		cout << f << endl;
		after[ibin]=(1.0-f)*before[ibin];
		after[ibin]+=f*before[ibin]*pow(xmid,-4.0)/2.4E25;
		file << xmid*10E8 << "\t" << before[ibin] << "\t" << after[ibin] << "\n";
	}

	cout << n0 << endl;

	file.close();	


	return 0;

}

