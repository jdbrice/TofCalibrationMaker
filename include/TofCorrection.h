#ifndef TOF_CORRECTION_H
#define TOF_CORRECTION_H

#include <memory>
#include "HistoBins.h"
using namespace jdb;
using namespace std;


class TofCorrection
{
protected:

	double t0;
	HistoBins * totBins;
	vector<double> totCorrs;

	HistoBins * zBins;
	vector<double> zCorrs;

public:
	TofCorrection( vector<double> totBinEdges, vector<double> zBinEdges ){

		totBins = ( new HistoBins( totBinEdges ) );
		zBins = ( new HistoBins( zBinEdges ) );
		t0 = 0;

		for ( int i = 0; i < totBins->nBins(); i++  ){
			totCorrs.push_back( 0.0 );
		}

		for ( int i = 0; i < zBins->nBins(); i++  ){
			zCorrs.push_back( 0.0 );
		}

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
	}

	double tof( double raw, double tot, double z ){
		int tBin = totBins->findBin( tot );
		int zBin = zBins->findBin( z );
		return raw - getTot( tBin ) - t0;
	}

	double tofForT0( double raw, double tot, double z ){
		int tBin = totBins->findBin( tot );
		int zBin = zBins->findBin( z );
		return raw - getTot( tBin );
	}

	double tofForTot( double raw, double z ){
		int zBin = zBins->findBin( z );
		return raw - t0;
	}

	void setT0( double _t0 ){ t0 = _t0; }
	void setTot( int bin, double val ){ if ( bin >= 0 && bin <= totCorrs.size() ) totCorrs[ bin ] = val; }
	void setZ( int bin, double val ){ if ( bin >= 0 && bin <= zCorrs.size() ) zCorrs[ bin ] = val; }

	double getTot( int bin ) const { if ( bin >= 0 && bin < totCorrs.size() ) return totCorrs[ bin ]; }


	
};

#endif