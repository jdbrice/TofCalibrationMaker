#include "SplineMaker.h"

#include <iostream>
using namespace std;
SplineMaker::SplineMaker( const vector< double > &x, const vector< double > &y, Interpolation::Type type ){
	
	/*for ( int i = 0; i < x.size(); i++ ){
		cout << "x[" << i << " ] = " << x[i] << " => " << y[ i ] << endl;
	}*/

	spline = NULL;
	spline = new Interpolator( x, y, type);

	domainMin = x[ 0 ];
	domainMax = x[ x.size() - 1 ];
	//cout << "Domain : ( " << domainMin << ", " << domainMax << " ) " << endl;

	

}

SplineMaker::SplineMaker( TH1D* hist, int place, Interpolation::Type type , int firstBin , int lastBin ){
	spline = NULL;
	if ( !hist ) 
		return;

	int nBins = hist -> GetNbinsX();
	int fBin = firstBin;
	int lBin = lastBin;

	
	if ( lBin <= -1 || lBin >= nBins )
		lBin = nBins;
	if ( fBin < 0 || fBin >= nBins )
		fBin = 1;

	vector<double> x((lBin - fBin) + 1 + 2);
	vector<double> y((lBin - fBin) + 1 + 2);


	int j = 1;
	x[ 0 ] = hist->GetBinLowEdge( fBin );
	y[ 0 ] = hist->GetBinContent( fBin );
	for ( int i = fBin; i <= lBin; i++){

		double bEdge = hist->GetBinLowEdge( i );
		double bWidth = hist->GetBinWidth( i );

		double _y = hist->GetBinContent( i );
		double _x = bEdge;
		if ( SplineMaker::AlignLeft == place )
			_x = bEdge + .0000001; 		// makes sure there are no troubles with doubles on the edge
		else if ( SplineMaker::AlignCenter == place )
			_x = (bEdge + (bWidth/2.0));
		else if ( SplineMaker::AlignRight == place )
			_x = ( bEdge + bWidth - .0000001 );// makes sure there are no troubles with doubles on the edge

		x[ j ] = _x;
		y[ j ] = _y;

		j++;
	}
	x[ j ] = hist->GetBinLowEdge( lBin ) + hist->GetBinWidth( lBin );
	y[ j ] = hist->GetBinContent( lBin );

	domainMin = x[ 0 ];
	domainMax = x[ x.size() - 1 ];

	spline = new Interpolator( x, y, type);

}

SplineMaker::~SplineMaker(){

	if ( spline )
		delete spline;
	spline = NULL;
}


TGraph * SplineMaker::graph( double xmin, double xmax, double step ){
	
	if ( xmin < domainMin )
		xmin = domainMin;
	if ( xmax > domainMax )
		xmax = domainMax;

   	const Int_t n = ( (xmax - xmin ) / step) + 1 ;
   	Int_t i = 0;
   	Float_t xcoord[n], ycoord[n];

   	for ( double xi = xmin; xi < xmax; xi += step) { 
   		//cout << "\tx "  << xi << endl;
      	xcoord[i] = xi;
      	ycoord[i] = spline->Eval(xi);
      	i++; 
   	}
   	
   	TGraph *gr = new TGraph( n, xcoord, ycoord );

	return gr;
}


double SplineMaker::eval( double x ){

	double ex = x;
	if ( ex < domainMin )
		ex = domainMin;
	else if ( ex > domainMax )
		ex = domainMax;

	if ( !spline )
		return 0;

	return spline->Eval( ex );
}