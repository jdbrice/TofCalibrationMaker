#ifndef TOF_CALIBRATION_MAKER
#define TOF_CALIBRATION_MAKER


/**
 * JDB LIB
 */
#include "XmlConfig.h"
#include "Logger.h"
#include "DataSource.h"
#include "HistoBook.h"
#include "LoggerConfig.h"
#include "Reporter.h"
#include "ConfigRange.h"
using namespace jdb;

/**
 * STL
 */
#include <memory>
using namespace std;

/**
 * ROOT
 */
#include "TMath.h"

/**
 * Local
 */
#include "TofCorrection.h"



class TofCalibrationMaker
{
protected:

	/**
	 * Config
	 */
	XmlConfig * cfg;
	string nodePath;

	unique_ptr<Logger> logger;

	unique_ptr<DataSource> ds;
	int nEventsToProcess;

	unique_ptr<HistoBook> book;
	unique_ptr<Reporter> reporter;


	string splitMode;
	unique_ptr<ConfigRange> trayRange;
	unique_ptr<ConfigRange> moduleRange;
	unique_ptr<ConfigRange> cellRange;
	vector< unique_ptr<TofCorrection> > corrections;

	int nElements; // the number of separate calibrated elements 

	int iteration;

	const double cLight; //= 29.9792458 cm / ns
	const double mPi; // pi mass = 0.13957 in GeV / c^2

public:
	TofCalibrationMaker( XmlConfig * config, string np, string fileList ="", string jobPrefix = "" );
	~TofCalibrationMaker();

	void make();

	void alignT0();

	void FillTot();
	//void slewing();

	void FillZLocal();

	void binTot();

	static const int nTrays;
	static const int nModules;
	static const int nCells;

protected:

	int absIndex( int tray, int module = 1, int cell = 1 );
	int relIndex( int tray, int module = 1, int cell = 1 );

	void importZParams( string pFile );
	void importTotParams( string totFile );

	double expectedTof( double length, double p ){
		return TMath::Sqrt( length*length / (cLight*cLight) * ( 1 + mPi*mPi / (p*p) ) );
	}

	bool keepEvent();
	bool keepTrack( int iHit );

	
	
};

#endif