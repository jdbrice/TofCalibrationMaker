#ifndef TOF_CALIBRATION_MAKER_H
#define TOF_CALIBRATION_MAKER_H

#include "TreeAnalyzer.h"
#include "XmlRange.h"

#include "TofCorrection.h"

class TofCalibrationMaker : public TreeAnalyzer
{
public:
	virtual char* classname() const { return "TofCalibrationMaker"; }
	TofCalibrationMaker() {}
	~TofCalibrationMaker() {}

	virtual void init(  XmlConfig &_config, string _nodePath="", int _jobIndex = -1 );

	void alignT0();
	void fillTot();
	void correctTot();
	
	void fillZLocal();
	void correctZLocal();
	void inverseBeta();

	static const int nTrays;
	static const int nModules;
	static const int nCells;
protected:

	virtual void overrideConfig() {
		DEBUG( classname(), "" );

		int _jobIndex = config.getInt( "jobIndex", -1 );

		if ( _jobIndex >= 0 ){
			map<string, string> opts;
			INFO( classname(), "Overriding config" );
			INFO( classname(), "jobIndex = " << _jobIndex );
			opts[ nodePath + ".Trays:min" ] = ts( _jobIndex + 1 );
			opts[ nodePath + ".Trays:max" ] = ts( _jobIndex + 1 );

			config.applyOverrides( opts );

		}
	}

	virtual void make() {
		
		INFO( classname(), "Making" );
		if ( config.exists( nodePath+".histograms" ) )
			book->makeAll( nodePath+".histograms" );

	    int nIterations = config.getInt( nodePath + ":nIterations", 3 );
	    for ( int i = 0; i < nIterations; i++ ){
	    	INFO( classname(), "Iteration " << i );
	        inverseBeta();   
	        
	        alignT0();
	        
	        fillTot();
	        correctTot();
	        // reportTot();
	        
	        fillZLocal();
	        correctZLocal();
	        // reportZLocal();

	        iteration++;
	    }

	    alignT0();
	    inverseBeta();

	    string eName = "t_" + ts( (int)trayRange->min ) + "_" + ts((int)trayRange->max) +
	                    "m_" + ts( (int)moduleRange->min ) + "_" + ts((int)moduleRange->max) +
	                    "c_" + ts( (int)cellRange->min ) + "_" + ts((int)cellRange->max);
	    exportT0Params( eName + "__t0.dat" );
	    exportTotParams( eName + "__tot.dat" );
	    exportZParams( eName + "__z.dat" );


	}


	int iteration = 0;	
	static const double cLight; //= 29.9792458 cm / ns
	static const double mPi; // pi mass = 0.13957 in GeV / c^2

	string splitMode;
	unique_ptr<XmlRange> trayRange;
	unique_ptr<XmlRange> moduleRange;
	unique_ptr<XmlRange> cellRange;
	vector< unique_ptr<TofCorrection> > corrections;

	int nElements; // the number of separate calibrated elements 
	

	bool keepEvent();
	bool keepTrack( int iHit );

	int absIndex( int tray, int module = 1, int cell = 1 );
	int relIndex( int tray, int module = 1, int cell = 1 );
	vector<int> fromRelIndex( int id ); // [0] = tray, [1] = module, [2] = cell

	void importZParams( string pFile );
	void importTotParams( string totFile );

	void exportT0Params( string pFile );
	void exportTotParams( string pFile );
	void exportZParams( string pFile );

	double expectedTof( double length, double p ){
		return TMath::Sqrt( length*length / (cLight*cLight) * ( 1 + mPi*mPi / (p*p) ) );
	}


	void makeBins( string var, int nBins, double min, double max );

	void reportTot();
	void reportZLocal();

	string nameFor( int rID ){

		vector<int> place = fromRelIndex( rID );
		if ( "tray" == splitMode )
			return "Tray "+ts(place[0]);
		else if ( "board" == splitMode )
			return "Tray "+ts(place[0]) + " Module " + ts( place[1] );
		else if ( "module" == splitMode )
			return "T"+ts(place[0]) + "M" + ts( place[1] );
		else if ( "cell" == splitMode )
			return "T"+ts(place[0]) + "M" + ts( place[1] )+ "C" + ts( place[2] );
		return "unknown";
	}

};


#endif