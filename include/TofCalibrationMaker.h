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


public:
	TofCalibrationMaker( XmlConfig * config, string np, string fileList ="", string jobPrefix = "" );
	~TofCalibrationMaker();

	void make();

	void alignT0();

	void FillZLocal();

	void binTot();

protected:

	int absIndex( int tray, int module = 1, int cell = 1 );
	int relIndex( int tray, int module = 1, int cell = 1 );

	
	
};

#endif