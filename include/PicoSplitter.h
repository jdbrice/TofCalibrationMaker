#ifndef PICO_SPLITTER_H
#define PICO_SPLITTER_H


/**
 * JDB LIB
 */
#include "XmlConfig.h"
#include "Logger.h"
#include "DataSource.h"
#include "HistoBook.h"
#include "LoggerConfig.h"
#include "Reporter.h"
#include "Utils.h"

#include "TofCellData.h"
using namespace jdb;

/**
 * STL
 */
#include <memory>
using namespace std;

/**
 * ROOT
 */
#include "TTree.h"




class PicoSplitter
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

	// for split modes
	vector< TTree* > forest;
	vector< TofCellData > data;
	vector< TFile* > files;
	vector< int > nHits;
	string splitMode;
	int nTrees;

	// process only some trays at a time to make job manageable
	int firstTray, lastTray;

	unique_ptr<HistoBook> book;



public:
	PicoSplitter( XmlConfig * config, string np, string fileList ="", string jobPrefix = "" );
	~PicoSplitter();

	int absIndex( int tray, int module = 1, int cell = 1 );

	int relIndex( int tray, int module = 1, int cell = 1 );

	void setupTrees();

	void make();

	void addTrack( int id, int iHit );

	bool keepEvent();

	
};

#endif