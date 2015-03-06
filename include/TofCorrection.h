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
	HistoBins * zBins;

public:
	TofCorrection( vector<double> totBinEdges, vector<double> zBinEdges ){

		totBins = ( new HistoBins( totBinEdges ) );
		zBins = ( new HistoBins( zBinEdges ) );
		t0 = 0;

	}
	~TofCorrection(){
		if ( totBins )
			delete totBins;
		if ( zBins )
			delete zBins;
	}

	double tof( double rawTof ){
		return rawTof - t0;
	}

	void setT0( double _t0 ){ t0 = _t0; }



	
};

#endif