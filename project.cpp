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

/*
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
*/



int main(int argc, char *argv[]) {

    // spectral radiance, maximum spectral radiance, intensity;
	double y,ymax=0,I; 
	//visible wavelength range
	double xmin=400E-9, xmax=700E-9, xdiff=xmax-xmin, xmid, dx=1E-9; 

	// use wein's law to find maximum spectral radiance
	ymax=(2.0*h*c*c*pow(b/T,-5.0))*pow(exp(h*c/(b*kb))-1,-1.0);

	// set up for histogram
	int ntot=10E5, nbins=30, ibin;
	double before[nbins]={0.0}, after[nbins]={0.0}, binsize=xdiff/nbins;

	// open a data file
	ofstream file;
	file.open("data.txt"); 

	// generate ntot random wavelengths proportional to 
	// Planck's distribution
	for(int n=0; n<ntot; n++){

		// give names to random numbers 
		double rx,ry;

		// set up random number generator
		mt19937 generator(time(0));

		// throw a random wavelength in visible range
		rx=xmin+xdiff*(random01(generator));

		// calculate the theoretical spectral radiance at rx
		y=(2.0*h*c*c*pow(rx,-5.0))*pow(exp(h*c/(rx*kb*T))-1,-1.0);

		// "keep=true" means we found a probable wavelength
		bool keep=false;
		while(keep==false){

			// throw a random spectral radiance ranging from 0 to ymax
			ry=ymax*(random01(generator));

			// if ry<y, then keep the random wavelength rx
			if(ry<y) {
				keep=true;
				//file << rx << "\t" << ry << endl;
				ibin=lrint(floor((rx-xmin)/binsize));
				before[ibin]+=1.0;
			}

			// if ry>y, we go back to the beginning of this loop
			else continue;
		}
	}

	/*

	up to this point, the program is doing what we expect
	we pick a random rx in the visible wavelength range 400-700 nm
	then we pick a random ry in the range 0-ymax
	so we have a random coordinate in the plane (rx,ry)
	we compare (rx,ry) with (rx,B(rx)), where B is Planck's distribution function
	we end up with a bunch of points (rx,ry) that only lie under the curve B

	the value of ry doesn't really matter after this point
	all that matters, is the number of rx's in a given range [x,x+dx]
	if we made a histogram of all the rx coordinates we found, 
	the histogram should be proportional to Planck's distribution
	but it's not, so our next goal is to figure this out


	*/



	for(ibin=0; ibin<nbins; ibin++){
		xmid=xmin+(ibin+0.5)*binsize;
		y=(2.0*h*c*c*pow(xmid,-5.0))*pow(exp(h*c/(xmid*kb*T))-1,-1.0);
		after[ibin]=before[ibin]*pow(xmid,-4.0);
		file << xmid << "\t" << before[ibin] << "\t" << after[ibin] << "\n";
	}


	file.close();	


	return 0;

}



/*

TO DO:

- ERROR CHECK: random wavelengths are not proportional to Planck's distribution
  need to go through each step and make sure it's working like we expect.

- Add dependency on radius of air molecules into intensity function. 
  Randomly scatter off N2 80% of the time and O2 20%. 
  We will have to estimate N2 and O2 atoms to be spheres. 

- Write something that will plot this data.

*/