#include "TofCalibrationMaker.h"


#include "TMath.h"
#include <math.h>

const int TofCalibrationMaker::nTrays = 120;
const int TofCalibrationMaker::nModules = 32;
const int TofCalibrationMaker::nCells = 6;


/**
 * Constructor
 * sets up the data source, config, logger
 */
TofCalibrationMaker::TofCalibrationMaker( XmlConfig * config, string np, string fileList, string jobPrefix ){

	cfg = config;
	nodePath = np;

	/**
	 * Logger
	 */
	logger = unique_ptr<Logger>(LoggerConfig::makeLogger( cfg, np + "Logger" ));
	Logger::setGlobalLogLevel( logger->getLogLevel() );
	logger->setClassSpace( "TofCalibrationMaker" );
	logger->info(__FUNCTION__) << "Got config with nodePath = " << np << endl;
	

	// create the book
    logger->info(__FUNCTION__) << " Creating book " << config->getString( np + "output.data", "TreeAnalyzer" ) << endl;
    book = unique_ptr<HistoBook>(new HistoBook( jobPrefix + config->getString( np + "output.data", "TreeAnalyzer" ), config, "", "" ) );
    	    
    if ( "" == jobPrefix && cfg->exists( np+"Reporter.output:url" ) ) {
	    reporter = unique_ptr<Reporter>(new Reporter( cfg, np+"Reporter.", jobPrefix ) );
	    logger->info(__FUNCTION__) << "Creating report " << config->getString( np+"Reporter.output:url" ) << endl;
    }

    /**
     * Sets up the input, should switch seemlessly between chain only 
     * and a DataSource 
     */
    if ( cfg->exists( np+"DataSource" ) ){
    	ds = unique_ptr<DataSource>(new DataSource( cfg, np + "DataSource", fileList ) );
    } else {
    	logger->error(__FUNCTION__) << "No DataSource given " << endl;
    }


    splitMode = cfg->getString( "TofCalibration:splitMode", "tray" );
    trayRange = unique_ptr<ConfigRange>( new ConfigRange( cfg, "TofCalibration.Trays", 1, 120 ) );
    moduleRange = unique_ptr<ConfigRange>( new ConfigRange( cfg, "TofCalibration.Modules", 1, 32 ) );
    cellRange = unique_ptr<ConfigRange>( new ConfigRange( cfg, "TofCalibration.Cells", 1, 6 ) );


    int nTraysToProcess = trayRange->max - trayRange->min + 1;
    int nModsToProcess = moduleRange->max - moduleRange->min + 1;
    int nCellsToProcess = cellRange->max - cellRange->min + 1;
    nElements = nTraysToProcess;
    if ( "module" == splitMode )
        nElements = nTraysToProcess * nModsToProcess;
    else if ( "cell" == splitMode )
        nElements = nTraysToProcess * nModsToProcess * nCellsToProcess;
    else if ( "board" == splitMode )
        nElements = nTraysToProcess * 8;

    /*************************
        working here when stopped
    */
    logger->info(__FUNCTION__) << "Making " << nElements << " zerod corrections" << endl;
    HistoBins totBins( cfg, "TofCalibration.Bins.tot" );
    HistoBins zBins( cfg, "TofCalibration.Bins.z" );
    for ( int i = 0; i < nElements; i++ ){
        corrections.push_back( unique_ptr<TofCorrection>( new TofCorrection( totBins.getBins(), zBins.getBins() ) ) );
    }
    iteration = 0;
    logger->info(__FUNCTION__) << "Making zerpo corrections" << endl;


    /**
     * Param import
     */
    if ( cfg->exists( "TofCalibration.Import.TotParams:url" ) ){
        importTotParams( cfg->getString( "TofCalibration.Import.TotParams:url" ) );
    }


}


TofCalibrationMaker::~TofCalibrationMaker(){


}


void TofCalibrationMaker::make(){

    book->makeAll( "histograms" );

    int nIterations = cfg->getInt( "TofCalibration:nIterations", 3 );
    for ( int i = 0; i < nIterations; i++ ){
        logger->info(__FUNCTION__) << "T0 Step " << i << endl;
        alignT0();
        
        FillTot();



        iteration++;
    }


}

void TofCalibrationMaker::FillTot(  ){

    if ( !ds ){
        logger->error(__FUNCTION__) << "Invalid DataSource " << endl;
        return;
    }


    /**
     * Make histos
     */
    
    book->cd( "totStep_"+ts(iteration) );
    for ( int i = 0; i < nElements; i++ ){
        book->clone( "/", "dtVsTot", "totStep_"+ts(iteration), "tot_" + ts(i) );
    }

    const Double_t c_light = 29.9792458;
    TF1 * piBetaCut = new TF1( "piBetaCut", "TMath::Sqrt([0]*[0]+x*x)/x",0,5);
    piBetaCut->SetParameter(0, 0.4 );

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

        /**
         * Select good events
         */
        if ( 0 >= ds->getInt( "numberOfVpdEast" ) || 0 >= ds->getInt( "numberOfVpdWest" ) ) continue;
        //if ( ds->get( "vR" ) > 1.0 ) continue;
        if ( TMath::Abs( ds->get( "vpdVz" )  - ds->get("vertexZ") ) > 6.0 ) continue;

        int nTofHits = ds->getInt( "nTofHits" );
        for ( int iHit = 0; iHit < nTofHits; iHit++ ){

            // Start
            double p = ds->get( "pt", iHit ) * TMath::CosH( ds->get( "eta", iHit ) );
            double dEdx = ds->get("dedx", iHit );

            if ( 0.3 > p || 1.0 < p ) continue;
            if ( p < 0.6 && dEdx > 2.8 ) continue;
            if ( 25 < ds->getInt( "nHitsFit", iHit ) ) continue;

            //cout << " ( " << ds->getInt( "tray", iHit ) << ", " << ds->getInt( "module", iHit ) << ", " << ds->getInt( "cell", iHit ) << " ) " << endl;
            int tray = ds->getInt( "tray", iHit );
            int module = ds->getInt( "module", iHit );
            int cell = ds->getInt( "cell", iHit );
            int id = relIndex( tray, module, cell );

            if ( iteration == 0 ){
                book->cd();
                book->fill( "traysHit", tray );
                book->fill( "modulesHit", module );
                book->fill( "cellsHit", cell );
                book->cd( "totStep_"+ts(iteration) );
            }

            if ( id < 0 )
                continue;


            double tot = ds->get( "tot", iHit );
            double tLength = ds->get( "length", iHit );
            double M = 0.13957;
            double bGamma = p / 0.13957; // pi mass in GeV / c^2
            double velocity = TMath::Sqrt(1.0/(1.0/bGamma/bGamma+1.0))*c_light;
            double piTof = tLength / velocity;
            double tofExpected = TMath::Sqrt( tLength*tLength / (c_light*c_light) * ( 1 + M*M / (p*p) ) ); // in nanoseconds  
            double rawTof = ds->get( "leTime", iHit ) - ds->get( "tStart" );
            rawTof = corrections[id]->tofForTot( rawTof, 0 ); // aplly all corrections except tot
            double corrTof = corrections[id]->tof( rawTof, tot, 0 );

            double dt = rawTof - tofExpected;
            double iBeta = (corrTof / tLength)*c_light;

            // cut on the inverse beta curve
            //if ( 0 < iteration && (iBeta > piBetaCut->Eval( p ) || iBeta < 1 ) ) continue;   

            book->fill( "tot_" + ts(id), ds->get( "tot", iHit ), dt );

        } // loop on tofHits
        
    } // end loop on events
    logger->info(__FUNCTION__) << "Completed in " << t.elapsed() << endl;



    //reporter->newPage();
    //book->style( "iBeta" )->set( "draw", "colz" )->draw();
    //piBetaCut->Draw("same");
    //reporter->savePage();


}

void TofCalibrationMaker::FillZLocal(  ){

	if ( !ds ){
		logger->error(__FUNCTION__) << "Invalid DataSource " << endl;
		return;
	}


	/**
	 * Make histos
	 */
	
	book->makeAll( "histograms" );
    for ( int iTray = 0; iTray < 120; iTray++ ){
        //for ( int iMod = 0; iMod < 32; iMod++ ){
            //for ( int iCell = 0; iCell < 6; iCell++ ){
                int id = (iTray+1);
                book->clone( "dtVsZ", "dtZ_" + ts(id) );
                book->clone( "dtVsTot", "dtTot_" + ts(id) );
            //}
        //}
    }

	const Double_t c_light = 29.9792458;
	TF1 * piBetaCut = new TF1( "piBetaCut", "TMath::Sqrt([0]*[0]+x*x)/x",0,5);
  	piBetaCut->SetParameter(0, 0.4 );

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



    	/**
    	 * Select good events
    	 */
    	if ( 0 >= ds->getInt( "numberOfVpdEast" ) || 0 >= ds->getInt( "numberOfVpdWest" ) ) continue;
    	if ( ds->get( "vR" ) > 1.0 ) continue;
    	if ( TMath::Abs( ds->get( "vpdVz" )  - ds->get("vertexZ") ) > 6.0 ) continue;

    	int nTofHits = ds->getInt( "nTofHits" );
    	for ( int iHit = 0; iHit < nTofHits; iHit++ ){

    		double p = ds->get( "pt", iHit ) * TMath::CosH( ds->get( "eta", iHit ) );
    		double dEdx = ds->get("dedx", iHit );

    		if ( 0.3 > p || 1.0 < p ) continue;
    		if ( p < 0.6 && dEdx > 2.8 ) continue;
    		if ( 25 < ds->getInt( "nHitsFit", iHit ) ) continue;

    		//cout << " ( " << ds->getInt( "tray", iHit ) << ", " << ds->getInt( "module", iHit ) << ", " << ds->getInt( "cell", iHit ) << " ) " << endl;
    		int tray = ds->getInt( "tray", iHit );
    		int module = ds->getInt( "module", iHit );
    		int cell = ds->getInt( "cell", iHit );

            int id = (tray );
    		book->fill( "traysHit", tray );
    		book->fill( "modulesHit", module );
    		book->fill( "cellsHit", cell );


    		double tLength = ds->get( "length", iHit );
    		double bGamma = p / 0.13957; // pi mass in GeV / c^2
    		double vel = bGamma / TMath::Sqrt( 1.0 + bGamma*bGamma );
            vel *= c_light;

    		double piTof = tLength / vel;
    		double tof =  ds->get( "tStart" ) - ds->get( "leTime", iHit );
    		

    		double dt = tof - piTof;
    		double iBeta = tof / tLength*c_light;

    		// cut on the inverse beta curve
    		if ( iBeta > piBetaCut->Eval( p ) ) continue;
    		book->fill( "iBeta", p, iBeta );

    		
            book->fill( "dtVsZ", ds->get( "zLocal", iHit ), dt );
            book->fill( "dtVsTot", ds->get( "tot", iHit ), dt );
    		book->fill( "dt", dt );

            
            book->fill( "dtZ_" + ts(id), ds->get( "zLocal", iHit ), dt );
            book->fill( "dtTot_" + ts(id), ds->get( "tot", iHit ), dt );
            book->fill( "dt", dt );



    		


    	} // loop on tofHits


    	
    	
    	
	} // end loop on events
	logger->info(__FUNCTION__) << "Completed in " << t.elapsed() << endl;



	reporter->newPage();

	book->style( "iBeta" )->set( "draw", "colz" )->draw();
	piBetaCut->Draw("same");
	reporter->savePage();


}


void TofCalibrationMaker::alignT0(){

    if ( !ds ){
        logger->error(__FUNCTION__) << "Invalid DataSource " << endl;
        return;
    }


    /**
     * Make histos
     */
    book->cd( "step_"+ts(iteration) );
    book->clone( "/", "t0Tray", "step_"+ts(iteration), "t0Tray" );
    book->clone( "/", "iBeta", "step_"+ts(iteration), "iBeta" );

    book->make2D( "elementT0", "t0", nElements, -0.5, nElements-0.5, 200, -100, 100 );
    if ( 0 < iteration )
        book->make2D( "elementT0Corr", "corrected t0", nElements, -0.5, nElements-0.5, 200, -10, 10 );

    const Double_t c_light = 29.9792458;
    TF1 * piBetaCut = new TF1( "piBetaCut", "TMath::Sqrt([0]*[0]+x*x)/x",0,5);
    piBetaCut->SetParameter(0, 0.4 );

    TaskTimer t;
    t.start();

    double globalT0 = 0;
    double nPi = 0;

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


        /**
         * Select good events
         */
        if ( 0 >= ds->getInt( "numberOfVpdEast" ) || 0 >= ds->getInt( "numberOfVpdWest" ) ) continue;
        if ( ds->get( "vR" ) > 1.0 ) continue;
        if ( TMath::Abs( ds->get( "vpdVz" )  - ds->get("vertexZ") ) > 6.0 ) continue;

        int nTofHits = ds->getInt( "nTofHits" );
        for ( int iHit = 0; iHit < nTofHits; iHit++ ){

            double p = ds->get( "pt", iHit ) * TMath::CosH( ds->get( "eta", iHit ) );
            double dEdx = ds->get("dedx", iHit );

            if ( 0.3 > p || 1.0 < p ) continue;
            if ( p < 0.6 && dEdx > 2.8 ) continue;
            if ( 25 < ds->getInt( "nHitsFit", iHit ) ) continue;

            //cout << " ( " << ds->getInt( "tray", iHit ) << ", " << ds->getInt( "module", iHit ) << ", " << ds->getInt( "cell", iHit ) << " ) " << endl;
            int tray = ds->getInt( "tray", iHit );
            int module = ds->getInt( "module", iHit );
            int cell = ds->getInt( "cell", iHit );
            int id = relIndex( tray, module, cell );

            if ( id < 0 )
                continue;
            //cout << " ID "<< id << endl;

            //book->fill( "traysHit", tray );
            //book->fill( "modulesHit", module );
            //book->fill( "cellsHit", cell );


            double tLength = ds->get( "length", iHit );
            double M = 0.13957;
            double bGamma = p / 0.13957; // pi mass in GeV / c^2
            double velocity = TMath::Sqrt(1.0/(1.0/bGamma/bGamma+1.0))*c_light;
            double piTof = tLength / velocity;
            double tofExpected = TMath::Sqrt( tLength*tLength / (c_light*c_light) * ( 1 + M*M / (p*p) ) ); // in nanoseconds  
            double rawTof = ds->get( "leTime", iHit ) - ds->get( "tStart" );
            //double rawTof = ds->get( "tofCorr", iHit );

            double corrTof = corrections[id]->tof( rawTof, ds->get("tot", iHit), 0 );

            double dt = corrTof - tofExpected;
            double dtRaw = rawTof - tofExpected;
            double iBeta = (corrTof / tLength)*c_light;

            // cut on the inverse beta curve
            if ( 0 < iteration && (iBeta > piBetaCut->Eval( p ) || iBeta < 1 ) ) continue;
            
            book->fill( "iBeta", p, iBeta );    
            book->fill( "elementT0", id, dtRaw );
            if ( 0 < iteration )
                book->fill( "elementT0Corr", id, dt );


        } // loop on tofHits
    } // end loop on events
    logger->info(__FUNCTION__) << "Completed in " << t.elapsed() << endl;

    TH2 * hT0 = book->get2D( "elementT0" );
    for ( int i = 0; i < nElements; i++ ){
        TH1D * hTmp = (TH1D*)hT0->ProjectionY( "_tmp", i+1, i+1 );
        double newT0 = hTmp->GetMean(1);

        hTmp->SetDirectory( 0 );

        corrections[ i ]->setT0( newT0 );
    }



    reporter->newPage();

    reporter->savePage();

}




int TofCalibrationMaker::absIndex( int tray, int module, int cell ){

    if ( trayRange && moduleRange && cellRange ){
        if ( tray < trayRange->min || tray > trayRange->max )
            return -2;
        if ( module < moduleRange->min || module > moduleRange->max )
            return -2;
        if ( cell < cellRange->min || cell > cellRange->max )
            return -2;
    } else {
        logger->error(__FUNCTION__) << endl;
    }

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

int TofCalibrationMaker::relIndex( int tray, int module, int cell ){

    if ( trayRange && moduleRange && cellRange ){
        if ( tray < trayRange->min || tray > trayRange->max )
            return -2;
        if ( module < moduleRange->min || module > moduleRange->max )
            return -2;
        if ( cell < cellRange->min || cell > cellRange->max )
            return -2;
    } else {
        logger->error(__FUNCTION__) << endl;
    }
    logger->debug( __FUNCTION__ ) << "( " << tray << ", " << module << ", " << cell << " ) " << endl;
    logger->debug( __FUNCTION__ ) << trayRange->min << " -> " << trayRange->max << endl;

    if ( "tray" == splitMode )
        return (tray - trayRange->min);
    if ( "module" == splitMode )
        return ( (tray - trayRange->min) * 32 ) + (module - moduleRange->min);
    if ( "board" == splitMode )
        return ( (tray - trayRange->min) * 8 ) + ((module - moduleRange->min) / 4);
    if ( "cell" == splitMode )
        return ( (tray - trayRange->min) * (32 * 6 ) ) + ( (module - moduleRange->min) * 6 ) + (cell - cellRange->min);
    return -1;
}


void TofCalibrationMaker::importTotParams( string totFile ){

    logger->info(__FUNCTION__) << endl;
    ifstream params( totFile.c_str() );

    if ( !params.good() ){
        logger->error( __FUNCTION__ ) << "Invalid tot parameter file " << totFile << endl;
        return;
    }

    int trayId, moduleId, cellId, boardId;
    int nbin;
    int iCalibType;
    
    params >> iCalibType;
    logger->info(__FUNCTION__) << "CalibType : " << iCalibType << endl;

    for( int i = 0; i < nTrays; i++ ) {
        for( int j = 0; j < nModules; j++ ) {
            for( int l = 0; l < nCells; l++ ){
                
                params >> trayId >> moduleId >> cellId;
                params >> nbin;

                int ri = relIndex( trayId, moduleId, cellId );
                if ( ri >= 0 ){
                    logger->info( __FUNCTION__ ) << "Index : " << ri << endl;
                    logger->info(__FUNCTION__) << "( " << trayId << ", " << moduleId << ", " << cellId << " ) " << endl;
                    logger->info(__FUNCTION__ ) << "#Bins = " << nbin << endl;    
                }
                
                vector<double> binEdges;
                for(int k = 0; k <= nbin; k++ ) {
                    double bEdge = 0;
                    params >> bEdge;
                    binEdges.push_back( bEdge );
                }

                if ( ri >= 0 ){
                    corrections[ ri ]->makeTotBins( binEdges );
                }
                
                vector<double> cors;
                
                for(int k = 0; k <= nbin; k++) {
                    double bCont = 0;
                    params >> bCont;
                    if ( ri >= 0 )
                        corrections[ ri ]->setTot( k, bCont );
                }
            }//cell
        }//module
    }//tray

    logger->info( __FUNCTION__ ) << "Complete" << endl;



}




/*
void calib::binTot( bool variableBinning ) {

    logger->info(__FUNCTION__) << "Starting " << endl;

    if ( variableBinning )
        logger->info(__FUNCTION__ ) << "Variable Binning ToT Range :  " << minTOT << " -> " << maxTOT << endl;
    else
        logger->info(__FUNCTION__ ) << "Fixed Binning ToT Range :  " << minTOT << " -> " << maxTOT << endl;

    logger->info(__FUNCTION__ ) << "Using " << numTOTBins << " bins for TOT" << endl;

    startTimer();


    Int_t nevents = (int)_chain->GetEntries();
    vector<double> tots[ constants::nChannels];

    logger->info(__FUNCTION__ ) << "] Processing " <<  nevents << " events" << endl;

    for(Int_t i=0; i<nevents; i++) {
        _chain->GetEntry(i);

        progressBar( i, nevents, 75 );
        Int_t numEast = pico->numberOfVpdEast;
        Int_t numWest = pico->numberOfVpdWest;
     
        if( numWest > constants::minHits){
            
            for(Int_t j = 0; j < constants::endWest; j++) {
                Double_t tot = getX( j );
              
                if(tot > minTOT && tot < maxTOT ) 
                    tots[j].push_back(tot);
            }

        }

        if( numEast > constants::minHits ){
    
            for(Int_t j = constants::startEast; j < constants::endEast; j++) {
                Double_t tot = getX( j );
      
                if( tot > minTOT && tot < maxTOT) 
                    tots[j].push_back(tot);
            }

        }

    } // lopp events    

    // get a threshold for a dead detector
    int threshold = 0;
    for(Int_t i=0; i<constants::nChannels; i++) {
        Int_t size = tots[i].size();
        threshold += size;
    }
    threshold /= (double)constants::nChannels; // the average of all detectors
    threshold *= .25;

    // loop through the channels and determine binning
    for(Int_t i=0; i<constants::nChannels; i++) {
      
        Int_t size = tots[i].size();
        cout << "[calib.binTOT] Channel[ " << i << " ] : " << size << " hits" << endl;
        
        if( size < threshold ) { // check for dead channels
            
            Double_t step = ( maxTOT - minTOT ) / numTOTBins;

            for(Int_t j=0; j <= numTOTBins; j++) {

                totBins[ i ][ j ] = ( step * j ) + minTOT; 
            }
            cout  << "[calib.binTOT] VPD Channel [ " << i << " ] is dead! " << "( " << size << " hits)" <<endl;
            
            // set this detector to dead
            deadDetector[ i ] = true;

        } else { // channel not dead

            deadDetector[ i ] = false;

            if ( variableBinning ){
                
                Int_t step = size / (numTOTBins + 1 ); 
        
                // sort into ascending order
                std::sort( tots[i].begin(), tots[i].end());
                
                totBins[ i ][0] = minTOT;
                totBins[ i ][ numTOTBins ] = maxTOT;
                
                for( Int_t j = 1; j < numTOTBins ; j++) {

                    double d1 = tots[i].at( step * j );
                    totBins[ i ][ j ] = d1;
                    
                }   // loop over tot bins
            }   // end variable binning 
            else { // fixed binning

                for ( int s = 0; s <= numTOTBins; s++ ){
                    double edge = ((maxTOT - minTOT) / (double) numTOTBins) * s;
                    edge += minTOT;
                    totBins[ i ][ s ] = edge;
                }
            
            }

      } // end channle not dead

    } // end loop channles
    
    for(Int_t i = 0; i< constants::nChannels; i++) {
        tots[i].clear();
    }   
    
    logger->info(__FUNCTION__ ) << "] completed in " << elapsed() << " seconds " << endl;
}
*/

/*
tof->Draw("leTime - tStart>>h1(10, -50, 50)", "numberOfVpdEast>0 && numberOfVpdWest>0 && abs(vertexZ-vpdVz) < 6.0 && TMath::Sqrt(vertexX*vertexX + vertexY*vertexY )< 1.0 && nHitsFit > 25& pt*TMath::CosH( eta ) < 1.0 & pt*TMath::CosH( eta ) > 0.3 && dedx < 2.8")
 */











