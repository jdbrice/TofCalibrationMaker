#ifndef TOF_CORRECTION_H
#define TOF_CORRECTION_H

#include <memory>
#include "HistoBins.h"
using namespace jdb;
using namespace std;

#include "SplineMaker.h"


class TofCorrection
{
protected:

	bool useSplineForTot;
	bool useSplineForZ;
	double t0;
	HistoBins * totBins;
	vector<double> totCorrs;
	SplineMaker* totSpline;

	HistoBins * zBins;
	vector<double> zCorrs;
	SplineMaker* zSpline;

public:
	TofCorrection( vector<double> totBinEdges, vector<double> zBinEdges, bool _splineTot = false, bool _splineZ = false ){

		totBins = ( new HistoBins( totBinEdges ) );
		zBins = ( new HistoBins( zBinEdges ) );
		t0 = 0;

		for ( int i = 0; i < totBins->nBins(); i++  ){
			totCorrs.push_back( 0.0 );
		}

		for ( int i = 0; i < zBins->nBins(); i++  ){
			zCorrs.push_back( 0.0 );
		}

		useSplineForZ = _splineZ;
		useSplineForTot = _splineTot;
		totSpline = NULL;
		zSpline = NULL;
		updateTotSpline();
		updateZSpline();

	}

	void makeTotBins( vector<double> nTotBins ){

		if ( totBins )
			delete totBins;
		totCorrs.clear();

		totBins = ( new HistoBins( nTotBins ) );

		for ( int i = 0; i < totBins->nBins(); i++  ){
			totCorrs.push_back( 0.0 );
		}
	}
	void makeZBins( vector<double> nZBins ){

		if ( zBins )
			delete zBins;
		zCorrs.clear();

		zBins = ( new HistoBins( nZBins ) );

		for ( int i = 0; i < zBins->nBins(); i++  ){
			zCorrs.push_back( 0.0 );
		}
	}

	~TofCorrection(){
		if ( totBins )
			delete totBins;
		if ( zBins )
			delete zBins;
		if ( totSpline )
			delete totSpline;
		if ( zSpline )
			delete zSpline;
	}

	double tof( double raw, double tot, double z ){
		return raw - getTot( tot ) - getZ( z ) - t0;
	}

	double tofForT0( double raw, double tot, double z ){
		return raw - getTot( tot ) - getZ( z );
	}

	double tofForTot( double raw, double z ){
		return raw - t0 - getZ( z );
	}
	double tofForZ( double raw, double tot ){
		return raw - getTot( tot ) - t0;
	}

	void setT0( double _t0 ){ t0 = _t0; }
	void setTot( int bin, double val ){ if ( bin >= 0 && bin <= totCorrs.size() ) totCorrs[ bin ] = val; }
	void setZ( int bin, double val ){ if ( bin >= 0 && bin <= zCorrs.size() ) zCorrs[ bin ] = val; }

	double getTot( double tot, bool forceBin = false ) const { 
		bool firstOrLast = (	totBins->findBin( tot ) == 0 || 
							totBins->findBin( tot ) >= totBins->nBins()-1 );
		if ( !useSplineForTot || NULL == totSpline || forceBin ){
			if ( !forceBin ){
				int bin = totBins->findBin( tot );
				if ( bin >= 0 && bin < totCorrs.size() ) 
					return totCorrs[ bin ]; 
			} else {
				int bin = (int)tot;
				if ( bin >= 0 && bin < totCorrs.size() ) 
					return totCorrs[ bin ]; 
			}
		} else {
			if ( !firstOrLast )
				return totSpline->eval( tot );
			else { // splines arent garanteed to exist at the boundaries
				int bin = totBins->findBin( tot );
				if ( bin >= 0 && bin < totCorrs.size() ) 
					return totCorrs[ bin ]; 
			}
		}
		return 0;
	}
	double getZ( double z, bool forceBin = false ) const { 
		if ( !useSplineForZ || NULL == zSpline || forceBin){
			if (!forceBin){
				int bin = zBins->findBin( z );
				if ( bin >= 0 && bin < zCorrs.size() ) 
					return zCorrs[ bin ]; 		
			} else {
				int bin = (int) z;
				if ( bin >= 0 && bin < zCorrs.size() ) 
					return zCorrs[ bin ]; 		
			}
		} else {

			return zSpline->eval( z );
		}
		return 0;
	}
	double getT0( ) const { return t0; }

	HistoBins * getTotBins() const { return totBins; }
	HistoBins * getZBins() const { return zBins; }
	SplineMaker * getTotSpline() const { return totSpline; }
	SplineMaker * getZSpline() const { return zSpline; }

	void updateTotSpline() {
		if ( totSpline )
			delete totSpline;

		vector<double> x = totBins->getBins();
		vector<double> y = totCorrs;

		x.push_back( x[ x.size() - 1 ]+.1 );
		y.push_back( y[ y.size() - 1 ] );
		
		totSpline = ( new SplineMaker( x, y, Interpolation::Type::kLINEAR) );
	}

	void updateZSpline() {
		if ( zSpline )
			delete zSpline;

		
		vector<double> x = zBins->getBins();
		vector<double> y = zCorrs;

		x.push_back( x[ x.size() - 1 ] );
		y.push_back( y[ y.size() - 1 ] );

		zSpline = ( new SplineMaker( x, y ) );
		
	}

	
};

#endif