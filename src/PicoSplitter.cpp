#include "PicoSplitter.h"



PicoSplitter::PicoSplitter( XmlConfig * config, string np, string fileList, string jobPrefix ){


	cfg = config;
	nodePath = np;

	/**
	 * Logger
	 */
	logger = unique_ptr<Logger>(LoggerConfig::makeLogger( cfg, np + "Logger" ));
	Logger::setGlobalLogLevel( logger->getLogLevel() );
	logger->info(__FUNCTION__) << logger->getLogLevel() << endl;
	logger->setClassSpace( "TofCalibrationMaker" );
	logger->info(__FUNCTION__) << "Got config with nodePath = " << np << endl;
	

    /**
     * Sets up the input, should switch seemlessly between chain only 
     * and a DataSource 
     */
    if ( cfg->exists( np+"DataSource" ) ){
    	ds = unique_ptr<DataSource>(new DataSource( cfg, np + "DataSource", fileList ) );
    } else {
    	logger->error(__FUNCTION__) << "No DataSource given " << endl;
    }


    splitMode = cfg->getString( "SplitBy", "tray" );


    firstTray = cfg->getInt( "Trays:min", 1 );
    lastTray = cfg->getInt( "Trays:max", 120 );
    int absFirstIndex = absIndex( firstTray );
    logger->info(__FUNCTION__) << "Tree First Index " << absFirstIndex << endl;

    int nTraysToProcess = lastTray - firstTray + 1;
    nTrees = nTraysToProcess;
	
	if ( "module" == splitMode )
		nTrees = nTraysToProcess * 32;
	else if ( "cell" == splitMode )
		nTrees = nTraysToProcess * 32 * 6;
	else if ( "board" == splitMode )
		nTrees = nTraysToProcess * 8;

    // make the trees
    for ( int iTree = 0; iTree < nTrees; iTree++ ){
    	TofCellData tcd;

    	files.push_back( (new TFile(( cfg->getString( "output:path", "" ) +  "tree_" + ts( absFirstIndex + 1 + iTree ) + ".root").c_str(), "RECREATE") ) );
    	forest.push_back( ( new TTree( "tof", "Split Tof Trees" ) ) );
    	data.push_back( tcd );
    	nHits.push_back( 0 );

    }

    setupTrees();

    book = unique_ptr<HistoBook>( new HistoBook( "info.root", cfg) );
    book->makeAll( "histograms" );

}

PicoSplitter::~PicoSplitter(){

	logger->info(__FUNCTION__) << endl;

	for ( int iTree = 0; iTree < forest.size(); iTree++ ){
		
		logger->debug(__FUNCTION__) << "Saving " << ("tree_" + ts(iTree) + ".root") << endl;
		if ( files[ iTree ] ){
			files[ iTree ]->Write();
			files[ iTree ]->Close();
		}else {
			cout << " BAD FILE" << endl;
		}

	}

}

int PicoSplitter::absIndex( int tray, int module, int cell ){

	if ( "tray" == splitMode )
		return (tray - 1) ;
	if ( "module" == splitMode )
		return ( (tray - 1) * 32 ) + (module-1);
	if ( "board" == splitMode )
		return ( (tray - 1) * 8 ) + ((module-1) / 4);
	if ( "cell" == splitMode )
		return ( (tray - 1) * (32 * 6 ) ) + ( (module-1) * 6 ) + (cell-1);
	return -1;
}

int PicoSplitter::relIndex( int tray, int module, int cell ){

	if ( tray < firstTray || tray > lastTray )
		return -2;

	if ( "tray" == splitMode )
		return (tray - firstTray) ;
	if ( "module" == splitMode )
		return ( (tray - firstTray) * 32 ) + (module-1);
	if ( "board" == splitMode )
		return ( (tray - firstTray) * 8 ) + ((module-1) / 4);
	if ( "cell" == splitMode )
		return ( (tray - firstTray) * (32 * 6 ) ) + ( (module-1) * 6 ) + (cell-1);
	return -1;
}

void PicoSplitter::setupTrees(){

	for ( int iTree = 0; iTree < nTrees; iTree++ ){

		forest[ iTree ]->SetAutoSave(1000);
					
		forest[ iTree ]->Branch("vertexX",				&data[iTree].vertexX,"vertexX/F");
		forest[ iTree ]->Branch("vertexY",				&data[iTree].vertexY,"vertexY/F");
		forest[ iTree ]->Branch("vertexZ",				&data[iTree].vertexZ,"vertexZ/F");
		
		forest[ iTree ]->Branch("numberOfVpdEast",		&data[iTree].numberOfVpdEast,"numberOfVpdEast/I");
		forest[ iTree ]->Branch("numberOfVpdWest",		&data[iTree].numberOfVpdWest,"numberOfVpdWest/I");
		

		forest[ iTree ]->Branch("tStart",				&data[iTree].tStart,"tStart/D");
		forest[ iTree ]->Branch("vpdVz",				&data[iTree].vpdVz,"vpdVz/F");

		forest[ iTree ]->Branch("nTofHits",				&data[iTree].nTofHits,"nTofHits/I");
		forest[ iTree ]->Branch("tray",					&data[iTree].tray,"tray[nTofHits]/I");
		forest[ iTree ]->Branch("module",				&data[iTree].module,"module[nTofHits]/I");
		forest[ iTree ]->Branch("cell",					&data[iTree].cell,"cell[nTofHits]/I");
		forest[ iTree ]->Branch("leTime",				&data[iTree].leTime,"leTime[nTofHits]/D");
		forest[ iTree ]->Branch("tot",					&data[iTree].tot,"tot[nTofHits]/D");
		forest[ iTree ]->Branch("matchFlag",			&data[iTree].matchFlag,"matchFlag[nTofHits]/I");

		forest[ iTree ]->Branch("zLocal",				&data[iTree].zLocal,"zLocal[nTofHits]/F");
		
		forest[ iTree ]->Branch("charge",				&data[iTree].charge,"charge[nTofHits]/I");
		forest[ iTree ]->Branch("pt",					&data[iTree].pt,"pt[nTofHits]/F");
		forest[ iTree ]->Branch("eta",					&data[iTree].eta,"eta[nTofHits]/F");
		
		forest[ iTree ]->Branch("length",				&data[iTree].length,"length[nTofHits]/F");

		forest[ iTree ]->Branch("nHitsFit",				&data[iTree].nHitsFit,"nHitsFit[nTofHits]/I");
		forest[ iTree ]->Branch("nHitsDedx",			&data[iTree].nHitsDedx,"nHitsDedx[nTofHits]/I"); 
		
		forest[ iTree ]->Branch("dedx",					&data[iTree].dedx,"dedx[nTofHits]/F"); 
		forest[ iTree ]->Branch("nSigE",				&data[iTree].nSigE,"nSigE[nTofHits]/F");
		forest[ iTree ]->Branch("nSigPi",				&data[iTree].nSigPi,"nSigPi[nTofHits]/F");
		forest[ iTree ]->Branch("nSigK",				&data[iTree].nSigK,"nSigK[nTofHits]/F");
		forest[ iTree ]->Branch("nSigP",				&data[iTree].nSigP,"nSigP[nTofHits]/F");
		
	}


}


void PicoSplitter::make(){

	TaskTimer t;
	t.start();

	Int_t nEvents = (Int_t)ds->getEntries();
	
	nEventsToProcess = cfg->getInt( nodePath+"DataSource:maxEvents", nEvents );

	if ( nEventsToProcess > nEvents )
		nEventsToProcess = nEvents;
	
	logger->info(__FUNCTION__) << "Loaded: " << nEventsToProcess << " events " << endl;
	
	TaskProgress tp( "Event Loop", nEventsToProcess );

	// loop over all events
	for(Int_t i=0; i<nEventsToProcess; i++) {
    	ds->getEntry(i);

    	tp.showProgress( i );

    	book->fill( "events", "All" );
    	if ( !keepEvent() )
    		continue;
    	book->fill( "events", "Good" );

    	for ( int iTree = 0; iTree < nTrees; iTree++ ){
    		nHits[ iTree ] = 0;
    	}

    	int nTofHits = ds->getInt( "nTofHits" );
    	for ( int iHit = 0; iHit < nTofHits; iHit++ ){

    		int tray = 		ds->getInt( "tray", 	iHit );
    		int module = 	ds->getInt( "module", 	iHit );
    		int cell = 		ds->getInt( "cell", 	iHit );
    		//logger->info( __FUNCTION__ ) << " ( " << tray << ", " << module << ", " << cell << " ) " << endl;	
    		int id = relIndex( tray, module, cell );
    		
    		book->fill( "trayHits", id );

    		if ( id >= 0 && nHits[ id ] < kMaxHits )
    			addTrack( id, iHit );
    		else if ( id >= 0 && nHits[ id ] >= kMaxHits )
    			logger->debug(__FUNCTION__) << "Losing hit on ( " << tray << ", " << module << ", " << cell << " ) " << endl;

    	} // loop on tofHits


    	for ( int iTree = 0; iTree < nTrees; iTree++ ){
    		if ( nHits[ iTree ] > 0 ){
    			forest[ iTree ]->Fill();
    		}
    	}


	} // end loop on events
	logger->info(__FUNCTION__) << "Completed in " << t.elapsed() << endl;

}




void PicoSplitter::addTrack( int id, int iHit ){

	data[ id ].vertexX = ds->get( "vertexX" );
	data[ id ].vertexY = ds->get( "vertexY" );
	data[ id ].vertexZ = ds->get( "vertexZ" );

	data[ id ].numberOfVpdWest 	= ds->getInt( "numberOfVpdWest" );
	data[ id ].numberOfVpdEast 	= ds->getInt( "numberOfVpdEast" );

	data[ id ].tStart 	= ds->get( "tStart" );
	data[ id ].vpdVz 	= ds->get( "vpdVz" );
	
	data[ id ].nTofHits = nHits[ id ] + 1;

	data[ id ].tray[ nHits[ id ] ] 		= 		ds->getInt( "tray", iHit );
	data[ id ].module[ nHits[ id ] ] 	= 		ds->getInt( "module", iHit );
	data[ id ].cell[ nHits[ id ] ] 		= 		ds->getInt( "cell", iHit );
	data[ id ].leTime[ nHits[ id ] ] 	= 		ds->get( "leTime", iHit );
	data[ id ].tot[ nHits[ id ] ] 		=	 	ds->get( "tot", iHit );
	data[ id ].matchFlag[ nHits[ id ] ] = 		ds->get( "matchFlag", iHit );
	data[ id ].zLocal[ nHits[ id ] ] 	= 		ds->get( "zLocal", iHit );

	data[ id ].charge[ nHits[ id ] ] 	= 		ds->get( "charge", iHit );
	data[ id ].pt[ nHits[ id ] ] 		= 		ds->get( "pt", iHit );
	data[ id ].eta[ nHits[ id ] ] 		= 		ds->get( "eta", iHit );

	data[ id ].length[ nHits[ id ] ] 	= 		ds->get( "length", iHit );

	data[ id ].nHitsFit[ nHits[ id ] ] 	= 		ds->get( "nHitsFit", iHit );
	data[ id ].nHitsDedx[ nHits[ id ] ] = 		ds->get( "nHitsDedx", iHit );
	data[ id ].dedx[ nHits[ id ] ] 		= 		ds->get( "dedx", iHit );

	data[ id ].nSigE[ nHits[ id ] ] 	= 		ds->get( "nSigE", iHit );
	data[ id ].nSigPi[ nHits[ id ] ] 	= 		ds->get( "nSigPi", iHit );
	data[ id ].nSigP[ nHits[ id ] ] 	= 		ds->get( "nSigP", iHit );
	data[ id ].nSigK[ nHits[ id ] ] 	= 		ds->get( "nSigK", iHit );
		
	nHits[ id ]++;
}



bool PicoSplitter::keepEvent( ){
	
	if ( 0 >= ds->getInt( "numberOfVpdEast" ) || 0 >= ds->getInt( "numberOfVpdWest" ) ) 
		return false;
    if ( ds->get( "vR" ) > 1.0 ) 
    	return false;
    if ( TMath::Abs( ds->get( "vpdVz" )  - ds->get("vertexZ") ) > 6.0 ) 
    	return false;

    return true;
}









