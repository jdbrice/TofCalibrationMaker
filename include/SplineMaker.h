#ifndef SPLINE_MAKER_H
#define SPLINE_MAKER_H

#include "TH1D.h"
#include "TMath.h"
#include "TGraph.h"
#include "Math/Interpolator.h"

using namespace ROOT::Math;
using namespace std;


class SplineMaker
{
public:
	const static int AlignLeft = 0;
	const static int AlignCenter = 1;
	const static int AlignRight = 2;

	// pass through constructor
	SplineMaker( const vector< double > &x, const vector< double > &y, Interpolation::Type type = Interpolation::kAKIMA );
	
	// from histogram
	SplineMaker( TH1D* hist, int place = SplineMaker::AlignLeft, Interpolation::Type type = Interpolation::kCSPLINE, int firstBin = 1, int lastBin = -1 );

	TGraph* graph( double xmin, double xmax, double step );
	//void draw( TH1D* hist, double xmin, double xmax, double step );

	double eval( double x );

	Interpolator* getSpline() { return spline; }

	~SplineMaker();



	

private:
	Interpolator* spline;
	double domainMin, domainMax;
};




#endif