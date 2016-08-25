
#ifndef TOF_PICO_SPLITTER_H
#define TOF_PICO_SPLITTER_H


#include "TreeAnalyzer.h"
#include "TofCellData.h"

class TofPicoSplitter : public TreeAnalyzer
{
public:
	virtual char* classname() const { return "TofPicoSplitter"; }
	TofPicoSplitter() {}
	~TofPicoSplitter() {}

	virtual void init( XmlConfig &_config, string _nodePath="", int _jobIndex = -1 ){
		
		TreeAnalyzer::init( _config, _nodePath, _jobIndex );

	    splitMode = config.getString( nodePath +".SplitBy", "tray" );


	    firstTray = config.getInt( nodePath + ".Trays:min", 1 );
	    lastTray = config.getInt( nodePath + ".Trays:max", 120 );
	    int absFirstIndex = absIndex( firstTray );
	    INFO( classname(), "Tree First Index " << absFirstIndex );

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

	    	files.push_back( (new TFile(( config.getString( nodePath + ".output:path", "" ) +  "tree_" + ts( absFirstIndex + 1 + iTree ) + ".root").c_str(), "RECREATE") ) );
	    	forest.push_back( ( new TTree( "tof", "Split Tof Trees" ) ) );
	    	data.push_back( tcd );
	    	nHits.push_back( 0 );

	    }


	    INFO( classname(), "Setting Up Trees" );
	    setupTrees();

	}



protected:

	// for split modes
	vector< TTree* > forest;
	vector< TofCellData > data;
	vector< TFile* > files;
	vector< int > nHits;
	string splitMode;
	int nTrees;

	// process only some trays at a time to make job manageable
	int firstTray, lastTray;


	int absIndex( int tray, int module = 1, int cell = 1 ){
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
	int relIndex( int tray, int module = 1, int cell = 1 ){

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
	
	void setupTrees(){

		for ( int iTree = 0; iTree < nTrees; iTree++ ){
			INFO( classname(), "Setting up " << iTree );
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
	} // setupTrees


	void addTrack( int id, int iHit ){

		// INFO( classname(), "Add Track " << id << ", " << iHit );

		data[ id ].vertexX = 		ds->get( "vertexX" );
		data[ id ].vertexY = 		ds->get( "vertexY" );
		data[ id ].vertexZ = 		ds->get( "vertexZ" );

		data[ id ].numberOfVpdWest 	= ds->get( "numberOfVpdWest" );
		data[ id ].numberOfVpdEast 	= ds->get( "numberOfVpdEast" );

		data[ id ].tStart 	= ds->get( "tStart" );
		data[ id ].vpdVz 	= ds->get( "vpdVz" );
		
		data[ id ].nTofHits = nHits[ id ] + 1;

		data[ id ].tray[ nHits[ id ] ] 		= 		ds->get( "tray", iHit );
		data[ id ].module[ nHits[ id ] ] 	= 		ds->get( "module", iHit );
		data[ id ].cell[ nHits[ id ] ] 		= 		ds->get( "cell", iHit );
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
	} // addTrack

	bool keepEvent( ){
		
		
		if ( 0 >= ds->get( "numberOfVpdEast" ) || 0 >= ds->get( "numberOfVpdWest" ) ) 
			return false;
		
	    if ( ds->get( "vR" ) > 1.0 ) 
	    	return false;
	    
	    if ( TMath::Abs( ds->get( "vpdVz" )  - ds->get("vertexZ") ) > 6.0 ) 
	    	return false;
	    

	    return true;
	} // keepEvent



	virtual void analyzeEvent(){

		for ( int iTree = 0; iTree < nTrees; iTree++ ){
			nHits[ iTree ] = 0;
		}

		int nTofHits = ds->get( "nTofHits" );
		for ( int iHit = 0; iHit < nTofHits; iHit++ ){

			int tray = 		ds->get( "tray", 	iHit );
			int module = 	ds->get( "module", 	iHit );
			int cell = 		ds->get( "cell", 	iHit );
			// INFO( classname(), " ( " << tray << ", " << module << ", " << cell << " ) " );	
			int id = relIndex( tray, module, cell );
			// INFO( classname(), "id = " << id );
			// book->fill( "trayHits", id );

			if ( id >= 0 && nHits[ id ] < kMaxHits )
				addTrack( id, iHit );
			else if ( id >= 0 && nHits[ id ] >= kMaxHits ){
				DEBUG( classname(), "Losing hit on ( " << tray << ", " << module << ", " << cell << " ) " );
			}

		} // loop on tofHits


		for ( int iTree = 0; iTree < nTrees; iTree++ ){
			if ( nHits[ iTree ] > 0 ){
				// INFO( classname(), "Filling " << iTree );
				forest[ iTree ]->Fill();
			}
		}


	}

	virtual void postMake() {
		for ( auto f : files ){
			f->Close();
		}
	}

	// virtual void make(){
	// 	DEBUG( classname(), "" );
	// 	TaskTimer t;
	// 	t.start();

	// 	Int_t nEvents = (Int_t)ds->getEntries();
		
	// 	nEventsToProcess = config.getInt( nodePath + ".DataSource:maxEvents", nEvents );

	// 	if ( nEventsToProcess > nEvents )
	// 		nEventsToProcess = nEvents;
		
	// 	DEBUG( classname(), "Loaded: " << nEventsToProcess << " events " );
		
	// 	TaskProgress tp( "Event Loop", nEventsToProcess );

	// 	// loop over all events
	// 	for(Int_t i=0; i<nEventsToProcess; i++) {
	//     	ds->getEntry(i);

	//     	tp.showProgress( i );

	//     	book->fill( "events", "All" );
	//     	if ( !keepEvent() )
	//     		continue;
	//     	book->fill( "events", "Good" );

	//     	for ( int iTree = 0; iTree < nTrees; iTree++ ){
	//     		nHits[ iTree ] = 0;
	//     	}

	//     	int nTofHits = ds->get( "nTofHits" );
	//     	for ( int iHit = 0; iHit < nTofHits; iHit++ ){

	//     		int tray = 		ds->get( "tray", 	iHit );
	//     		int module = 	ds->get( "module", 	iHit );
	//     		int cell = 		ds->get( "cell", 	iHit );
	//     		//logger->info( __FUNCTION__ ) << " ( " << tray << ", " << module << ", " << cell << " ) " << endl;	
	//     		int id = relIndex( tray, module, cell );
	    		
	//     		book->fill( "trayHits", id );

	//     		if ( id >= 0 && nHits[ id ] < kMaxHits )
	//     			addTrack( id, iHit );
	//     		else if ( id >= 0 && nHits[ id ] >= kMaxHits ){
	//     			DEBUG( classname(), "Losing hit on ( " << tray << ", " << module << ", " << cell << " ) " );
	//     		}

	//     	} // loop on tofHits


	//     	for ( int iTree = 0; iTree < nTrees; iTree++ ){
	//     		if ( nHits[ iTree ] > 0 ){
	//     			forest[ iTree ]->Fill();
	//     		}
	//     	}


	// 	} // end loop on events
	// 	INFO( classname(), "Completed in " << t.elapsed() );

	// }
	
};


#endif
