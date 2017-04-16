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
const double b = 2.898E-3;


double random01(mt19937 generator){
	static uniform_01<mt19937> dist(generator);
	return dist();
}

double SR(double wavelength){
	 return (2.0*h*c*c*pow(wavelength,-5.0))*pow(exp(h*c/(wavelength*kb*T))-1,-1.0);
}


int main(int argc, char *argv[]) {

    // spectral radiance, maximum spectral radiance, intensity;
	double y,yf,ymax=0,I; 
	//visible wavelength range
	double xmin=400E-9, xmax=700E-9, xdiff=xmax-xmin, xmid; 

	// use wein's law to find maximum spectral radiance
	ymax=SR(b/T);

	// set up for histogram
	int ntot=10E4, nbins=50, ibin;
	double before[nbins]={0.0}, after[nbins]={0.0}, binsize=xdiff/nbins;

	// set up random number generator
	mt19937 generator(time(0));

	// open a data file
	ofstream file;
	file.open("data.txt"); 

	// generate ntot random wavelengths proportional to 
	// Planck's distribution
	for(int n=0; n<ntot; n++){

		// give names to random numbers 
		double rx,ry;

		// throw a random wavelength in visible range
		rx=xmin+xdiff*(random01(generator));

		// calculate the theoretical spectral radiance at rx
		y=SR(rx);

		ry=ymax*(random01(generator));

		if(ry<y){
			//file << rx << "\t" << ry << endl;
			ibin=lrint(floor((rx-xmin)/binsize));
			before[ibin]+=1.0;

		}

	}


	for(ibin=0; ibin<nbins; ibin++){
		xmid=xmin+(ibin+0.5)*binsize;
		y=(2.0*h*c*c*pow(xmid,-5.0))*pow(exp(h*c/(xmid*kb*T))-1,-1.0);
		after[ibin]=before[ibin]*pow(xmid,-4.0);
		file << xmid*10E8 << "\t" << 1.5*pow(10,25)*before[ibin] << "\t" << after[ibin] << "\n";
	}



	file.close();	


	return 0;

}



/*

TO DO:

- Make sure everything is scaled correctly.

- Add dependency on radius of air molecules into intensity function. 
  Randomly scatter off N2 80% of the time and O2 20%. 
  We will have to estimate N2 and O2 atoms to be spheres. 

- Write something that will plot this data.

*/