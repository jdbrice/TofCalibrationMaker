#include "TofCalibrationMaker.h"

#include "ChainLoader.h"


#include <math.h>
#include "TProfile.h"

const int TofCalibrationMaker::nTrays = 120;
const int TofCalibrationMaker::nModules = 32;
const int TofCalibrationMaker::nCells = 6;


/**
 * Constructor
 * sets up the data source, config, logger
 */
TofCalibrationMaker::TofCalibrationMaker( XmlConfig * config, string np, int jobId, string fileList, string jobPrefix )
    : cLight(29.9792458 /* [cm/ns] */), mPi( 0.13957 /* [GeV/c^2] */ )
{

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
    	    
    if ( cfg->exists( np+"Reporter.output:url" ) ) {
	    reporter = unique_ptr<Reporter>(new Reporter( cfg, np+"Reporter.", jobPrefix ) );
	    logger->info(__FUNCTION__) << "Creating report " << config->getString( np+"Reporter.output:url" ) << endl;
    }

    /**
     * Sets up the input, should switch seemlessly between chain only 
     * and a DataSource 
     */
    /*if ( cfg->exists( np+"DataSource" ) && jobId <= 0 ){
    	ds = unique_ptr<DataSource>(new DataSource( cfg, np + "DataSource", fileList ) );
    } if ( cfg->exists( np+"DataSource" ) && jobId >= 1 ){
        ds = unique_ptr<DataSource>(new DataSource( cfg, np + "DataSource", cfg->getString( np + "DataSource:url") + ts(jobId)+".root" ) );
    } else {
    	logger->error(__FUNCTION__) << "No DataSource given " << endl;
    }*/

    string fName = cfg->getString( np + "DataSource:url") + "tree_"+ ts(jobId)+".root";
    chain = new TChain( "tof" );
    ChainLoader::load( chain, fName );

    tuple = new TofTuple( chain );



    splitMode = cfg->getString( "TofCalibration:splitMode", "tray" );
    trayRange = unique_ptr<ConfigRange>( new ConfigRange( cfg, "TofCalibration.Trays", 1, 120 ) );
    moduleRange = unique_ptr<ConfigRange>( new ConfigRange( cfg, "TofCalibration.Modules", 1, 32 ) );
    cellRange = unique_ptr<ConfigRange>( new ConfigRange( cfg, "TofCalibration.Cells", 1, 6 ) );
    logger->info(__FUNCTION__) << "SplitMode " << splitMode << endl;
    logger->info(__FUNCTION__) << "Tray Range " << trayRange->toString() << endl;
    logger->info(__FUNCTION__) << "Module Range " << moduleRange->toString() << endl;
    logger->info(__FUNCTION__) << "Cell Range " << cellRange->toString() << endl;

    if ( jobId >= 1  ){
        trayRange->min = (int)(jobId / 8) + 1;
        trayRange->max = (int)(jobId / 8) + 1;
        int iB = jobId - ( trayRange->min - 1 ) * 8;
        
        moduleRange->min = ( iB - 1 ) * 4 + 1;
        moduleRange->max = ( iB ) * 4;

        cellRange->min = 1;
        cellRange->max = 6;
    }

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

    logger->info(__FUNCTION__) << "Calibrating on " << nElements << " elements" << endl;

    
    HistoBins totBins( cfg, "TofCalibration.Bins.tot" );
    HistoBins zBins( cfg, "TofCalibration.Bins.zLocal" );
    for ( int i = 0; i < nElements; i++ ){
        corrections.push_back( unique_ptr<TofCorrection>( 
            new TofCorrection( totBins.getBins(), 
                                zBins.getBins(), 
                                cfg->getBool( "TofCalibration.Spline:tot", false ),
                                cfg->getBool( "TofCalibration.Spline:zLocal", false ) ) 
            ) 
        );
    }
    iteration = 0;
    logger->info(__FUNCTION__) << "Setup null corrections for zLocal and tot" << endl;

    if ( cfg->getBool( "TofCalibration.Bins.tot:variable", false ) )
        makeBins( "tot", cfg->getDouble( "TofCalibration.Bins.tot:nBins" ), cfg->getDouble( "TofCalibration.Bins.tot:min" ), cfg->getDouble( "TofCalibration.Bins.tot:max" ) );
    if ( cfg->getBool( "TofCalibration.Bins.zLocal:variable", false ) )
        makeBins( "zLocal", cfg->getDouble( "TofCalibration.Bins.zLocal:nBins" ), cfg->getDouble( "TofCalibration.Bins.zLocal:min" ), cfg->getDouble( "TofCalibration.Bins.zLocal:max" ) );
    


    /**
     * Param import
     */
    if ( cfg->exists( "TofCalibration.Import.TotParams:url" ) ){
        importTotParams( cfg->getString( "TofCalibration.Import.TotParams:url" ) );
    }


}


TofCalibrationMaker::~TofCalibrationMaker(){


}


bool TofCalibrationMaker::keepEvent( ){

    if ( 0 >= tuple->numberOfVpdEast || 0 >= tuple->numberOfVpdWest ) 
        return false;
    double vX = tuple->vertexX;
    double vY = tuple->vertexY;
    double vR = TMath::Sqrt( vX * vX + vY * vY  );
    if ( vR > 1.0 ) 
        return false;
    if ( TMath::Abs( tuple->vpdVz  - tuple->vertexZ ) > 6.0 ) 
        return false;
    return true;
}
bool TofCalibrationMaker::keepTrack( int iHit ){
    double p = tuple->pt[ iHit ] * TMath::CosH( tuple->eta[ iHit ] );
    double nSigPi = tuple->nSigPi[ iHit ];

    if ( 0.3 > p || 0.6 < p ) 
        return false;
    if ( nSigPi > 2.0 ) 
        return false;
    if ( 25 > tuple->nHitsFit[ iHit ] ) 
        return false;
    return true;
}

void TofCalibrationMaker::make(){

    book->makeAll( "histograms" );

    int nIterations = cfg->getInt( "TofCalibration:nIterations", 3 );
    for ( int i = 0; i < nIterations; i++ ){
        inverseBeta();   
        
        alignT0();
        
        fillTot();
        correctTot();
        reportTot();
        
        fillZLocal();
        correctZLocal();
        reportZLocal();

        iteration++;
    }

    alignT0();
    inverseBeta();

    string eName =  "t_" + ts( (int)trayRange->min ) + "_" + ts((int)trayRange->max) +
                    "m_" + ts( (int)moduleRange->min ) + "_" + ts((int)moduleRange->max) +
                    "c_" + ts( (int)cellRange->min ) + "_" + ts((int)cellRange->max);
    exportT0Params( eName + "__t0.dat" );
    exportTotParams( eName + "__tot.dat" );
    exportZParams( eName + "__z.dat" );

}

void TofCalibrationMaker::fillTot(  ){

    if ( !chain ){
        logger->error(__FUNCTION__) << "Invalid DataSource " << endl;
        return;
    }


    /**
     * Make histos
     */
    
    book->cd( "totStep_"+ts(iteration) );
    HistoBins dtBins( cfg, "b.dt" );
    for ( int i = 0; i < nElements; i++ ){
        HistoBins * hb = corrections[ i ]->getTotBins();
        book->make2D( "tot_" + ts( i ), "dt vs. Tot", hb->nBins(), hb->getBins().data(), dtBins.nBins(), dtBins.getBins().data() );
        book->make2D( "corrTot_" + ts( i ), "dt vs. Tot", hb->nBins(), hb->getBins().data(), dtBins.nBins(), dtBins.getBins().data() );
    }

    const Double_t c_light = 29.9792458;
    TF1 * piBetaCut = new TF1( "piBetaCut", "TMath::Sqrt([0]*[0]+x*x)/x",0,5);
    piBetaCut->SetParameter(0, 0.4 );

    TaskTimer t;
    t.start();

    Int_t nEvents = (Int_t)chain->GetEntries();
    
    nEventsToProcess = cfg->getInt( nodePath+"DataSource:maxEvents", nEvents );

    if ( nEventsToProcess > nEvents )
        nEventsToProcess = nEvents;
    
    logger->info(__FUNCTION__) << "Loaded: " << nEventsToProcess << " events " << endl;
    
    TaskProgress tp( "Plotting Tot", nEventsToProcess );
    logger->info( __FUNCTION__ ) << "Iteration : " << iteration << endl;
    // loop over all events
    for(Int_t i=0; i<nEventsToProcess; i++) {
        chain->GetEntry(i);

        tp.showProgress( i );

        /**
         * Select good events (should already be done in splitter )
         */
        if ( !keepEvent() )
            continue;

        int nTofHits = tuple->nTofHits;
        for ( int iHit = 0; iHit < nTofHits; iHit++ ){

            if ( !keepTrack( iHit ) )
                continue;

            int tray = tuple->tray[ iHit ];
            int module = tuple->module[ iHit ];
            int cell = tuple->cell[ iHit ];
            int id = relIndex( tray, module, cell );

            if ( id < 0 )
                continue;

            const double tLength = tuple->length[ iHit ];
            const double p = tuple->pt[ iHit ] * TMath::CosH( tuple->eta[ iHit ] );
            const double tofExp = expectedTof( tLength, p );
            const double tot = tuple->tot[ iHit ];;
            
            double rawTof = tuple->leTime[ iHit ] - tuple->tStart;
            double corrTof = corrections[id]->tof( rawTof, tot, 0 ); // all corrections
            double tof = corrections[id]->tofForTot( rawTof, 0 ); // aplly all corrections except tot
            

            double dt = tof - tofExp;
            
            double iBeta = (corrTof / tLength ) * cLight;
            // cut on the inverse beta curve
            if ( 0 < iteration && (iBeta > piBetaCut->Eval( p ) ) ) continue;   

            book->fill( "tot_" + ts(id), tuple->tot[ iHit ], dt );
            book->fill( "corrTot_" + ts(id), tuple->tot[ iHit ], (corrTof - tofExp) );

        } // loop on tofHits
        
    } // end loop on events
    logger->info(__FUNCTION__) << "Completed in " << t.elapsed() << endl;


}
void TofCalibrationMaker::fillZLocal(  ){

    if ( !chain ){
        logger->error(__FUNCTION__) << "Invalid DataSource " << endl;
        return;
    }


    /**
     * Make histos
     */
    
    book->cd( "zLocalStep_"+ts(iteration) );
    HistoBins dtBins( cfg, "b.dtZ" );
    for ( int i = 0; i < nElements; i++ ){
        HistoBins * hb = corrections[ i ]->getZBins();
        book->make2D( "zLocal_" + ts( i ), "dt vs. zLocal", hb->nBins(), hb->getBins().data(), dtBins.nBins(), dtBins.getBins().data() );
        book->make2D( "corrZLocal_" + ts( i ), "dt vs. zLocal", hb->nBins(), hb->getBins().data(), dtBins.nBins(), dtBins.getBins().data() );
    }

    
    TF1 * piBetaCut = new TF1( "piBetaCut", "TMath::Sqrt([0]*[0]+x*x)/x",0,5);
    piBetaCut->SetParameter(0, 0.4 );

    TaskTimer t;
    t.start();

    Int_t nEvents = (Int_t)chain->GetEntries();
    
    nEventsToProcess = cfg->getInt( nodePath+"DataSource:maxEvents", nEvents );

    if ( nEventsToProcess > nEvents )
        nEventsToProcess = nEvents;
    
    logger->info(__FUNCTION__) << "Loaded: " << nEventsToProcess << " events " << endl;
    
    TaskProgress tp( "Event Loop", nEventsToProcess );
    logger->info( __FUNCTION__ ) << "Iteration : " << iteration << endl;
    // loop over all events
    for(Int_t i=0; i<nEventsToProcess; i++) {
        chain->GetEntry(i);

        tp.showProgress( i );

        /**
         * Select good events (should already be done in splitter )
         */
        if ( !keepEvent() )
            continue;

        int nTofHits = tuple->nTofHits;
        for ( int iHit = 0; iHit < nTofHits; iHit++ ){

            if ( !keepTrack( iHit ) )
                continue;

            int tray = tuple->tray[ iHit ];
            int module = tuple->module[ iHit ];
            int cell = tuple->cell[ iHit ];
            int id = relIndex( tray, module, cell );

            if ( id < 0 )
                continue;

            const double zLocal = tuple->zLocal[ iHit ];
            const double tLength = tuple->length[ iHit ];
            const double p = tuple->pt[ iHit ] * TMath::CosH( tuple->eta[ iHit ] );
            const double tofExp = expectedTof( tLength, p );
            const double tot = tuple->tot[ iHit ];;
            
            double rawTof = tuple->leTime[ iHit ] - tuple->tStart;
            double corrTof  = corrections[id]->tof( rawTof, tot, zLocal ); // all corrections
            double tof      = corrections[id]->tofForZ( rawTof, tot ); // aplly all corrections except tot
            

            double dt = tof - tofExp;
            
            double iBeta = (corrTof / tLength ) * cLight;
            // cut on the inverse beta curve
            if ( 0 < iteration && (iBeta > piBetaCut->Eval( p ) ) ) continue;     

            book->fill( "zLocal_" + ts(id), zLocal, dt );
            book->fill( "corrZLocal_" + ts(id), zLocal, (corrTof - tofExp) );

        } // loop on tofHits
        
    } // end loop on events
    logger->info(__FUNCTION__) << "Completed in " << t.elapsed() << endl;


}

void TofCalibrationMaker::correctZLocal(  ){

    /**
     * Make histos
     */
    book->cd( "zLocalStep_"+ts(iteration) );
    TaskTimer t;
    t.start();
    TaskProgress tp( "Determing zLocal Corrections", nEventsToProcess );

    for ( int i = 0; i < nElements; i++ ){

        double x1 = cfg->getDouble( "TofCalibration.Bins.zLocal:min" );
        double x2 = cfg->getDouble( "TofCalibration.Bins.zLocal:max" );

        tp.showProgress( i );

        TH2 * tot2D = book->get2D( "zLocal_"+ts(i) );
        TProfile * profile = tot2D->ProfileX();
        TF1 * pol1 = new TF1( "pol1", "pol1", x1, x2 );

        profile->SetDirectory( 0 );
        book->add( "zLocalProf_"+ts(i), (TH1D*)(profile->Clone( ("zLocalProf_"+ts(i)).c_str() )) );

        profile->Fit( pol1, "QR" );

        int l = corrections[ i ]->getZBins()->size();
        double zThere = x1;
        for ( int j = 0; j < l; j ++ ){
            if ( !cfg->getBool( "TofCalibration.Spline:zLocalPol1", false ) )
                corrections[ i ]->setZ( j, profile->GetBinContent( j+1 ) );
            else {
                zThere = x1 + j * (( x2 - x1 ) / (double)l);
                cout << " zThere " << zThere << endl;
                corrections[ i ]->setZ( j, pol1->Eval( zThere ) );
            }
        }
        corrections[ i ]->updateZSpline();
    } 
    logger->info(__FUNCTION__) << "Completed in " << t.elapsed() << endl;


}

void TofCalibrationMaker::correctTot(  ){

    /**
     * Make histos
     */
    book->cd( "totStep_"+ts(iteration) );
    TaskTimer t;
    t.start();
    TaskProgress tp( "Determing tot Corrections", nEventsToProcess );

    for ( int i = 0; i < nElements; i++ ){

        tp.showProgress( i );

        TH2 * tot2D = book->get2D( "tot_"+ts(i) );
        TProfile * profile = tot2D->ProfileX();
        profile->SetDirectory( 0 );
        book->add( "totProf_"+ts(i), (TH1D*)(profile->Clone( ("totProf_"+ts(i)).c_str() )) );

        int l = corrections[ i ]->getTotBins()->size();
        for ( int j = 0; j < l; j ++ ){
            corrections[ i ]->setTot( j, profile->GetBinContent( j+1 ) );
        }
        corrections[ i ]->updateTotSpline();

    } 
    logger->info(__FUNCTION__) << "Completed in " << t.elapsed() << endl;
}



void TofCalibrationMaker::alignT0(){

    if ( !chain ){
        logger->error(__FUNCTION__) << "Invalid DataSource " << endl;
        return;
    }

    /**
     * Make histos
     */
    book->cd( "t0Step_"+ts(iteration) );

    if ( 0 == iteration )
       book->make2D( "elementT0", "t0", nElements, -0.5, nElements-0.5, 200, -100, 100 ); // bin With = 1ns
   else if ( 1 == iteration )
        book->make2D( "elementT0", "t0", nElements, -0.5, nElements-0.5, 200, -10, 10 ); // bin Width = 0.1ns
    else 
        book->make2D( "elementT0", "t0", nElements, -0.5, nElements-0.5, 2000, -5, 5 ); // bin Width = 0.005ns
    



    TF1 * piBetaCut = new TF1( "piBetaCut", "TMath::Sqrt([0]*[0]+x*x)/x",0,5);
    piBetaCut->SetParameter(0, 0.4 );

    TaskTimer t;
    t.start();

    Int_t nEvents = (Int_t)chain->GetEntries();
    
    nEventsToProcess = cfg->getInt( nodePath+"DataSource:maxEvents", nEvents );

    if ( nEventsToProcess > nEvents )
        nEventsToProcess = nEvents;
    
    logger->info(__FUNCTION__) << "Loaded: " << nEventsToProcess << " events " << endl;
    
    TaskProgress tp( "T0 Alignment", nEventsToProcess );
    logger->info( __FUNCTION__ ) << "Iteration : " << iteration << endl;
    // loop over all events
    for(Int_t i=0; i<nEventsToProcess; i++) {
        chain->GetEntry(i);

        tp.showProgress( i );


        /**
         * Select good events
         */
        if ( !keepEvent() )
            continue;

        int nTofHits = tuple->nTofHits;
        for ( int iHit = 0; iHit < nTofHits; iHit++ ){

            if ( !keepTrack( iHit ) )
                continue;

            //cout << " ( " << ds->getInt( "tray", iHit ) << ", " << ds->getInt( "module", iHit ) << ", " << ds->getInt( "cell", iHit ) << " ) " << endl;
            int tray = tuple->tray[ iHit ];
            int module = tuple->module[ iHit ];
            int cell = tuple->cell[ iHit ];
            int id = relIndex( tray, module, cell );

            if ( id < 0 )
                continue;

            const double zLocal = tuple->zLocal[ iHit ];
            const double tLength = tuple->length[ iHit ];
            const double p = tuple->pt[ iHit ] * TMath::CosH( tuple->eta[ iHit ] );
            const double tofExp = expectedTof( tLength, p );
            const double tot = tuple->tot[ iHit ];;
            
            double rawTof = tuple->leTime[ iHit ] - tuple->tStart;
            const double corrTof = corrections[id]->tof( rawTof, tot, 0 ); // all corrections
            const double iBeta = ( corrTof / tLength ) * cLight;
    
            double dt = corrTof - tofExp;

            // cut on the inverse beta curve
            if ( 0 < iteration && (iBeta > piBetaCut->Eval( p ) ) ) continue;
            
            //book->fill( "iBeta", p, iBeta );    
            book->fill( "elementT0", id, dt );

        } // loop on tofHits
    } // end loop on events
    logger->info(__FUNCTION__) << "Completed in " << t.elapsed() << endl;

    
    TH2 * hT0 = book->get2D( "elementT0" );
    for ( int i = 0; i < nElements; i++ ){
        TH1D * hTmp = (TH1D*)hT0->ProjectionY( "_tmp", i+1, i+1 );
        double newT0 = hTmp->GetMean(1);
        if ( iteration >= 2 ){
            TF1 * gg = new TF1( "gaus", "gaus", -.25, .25 );
            hTmp->Fit( gg, "QRN" );
            newT0 = gg->GetParameter( 1 );    
        }
        
        hTmp->SetDirectory( 0 );

        double cT0 = corrections[ i ]->getT0();
        corrections[ i ]->setT0( cT0 + newT0 );
    }


    reporter->newPage();

    book->style( "elementT0" )->
    set( "title", "T0 : " + nameFor( 0 ) + " -> " + nameFor( nElements-1 ) )->
    set( "y", "#Delta TOF_{m} - TOF_{exp}")->
    set( "draw", "colz" )->set( "logz", 1 )->draw();

    reporter->savePage();

}


void TofCalibrationMaker::inverseBeta(  ){

    if ( !chain ){
        logger->error(__FUNCTION__) << "Invalid DataSource " << endl;
        return;
    }

    /**
     * Make histos
     */
    book->cd( "iBetaStep_"+ts(iteration) );
    HistoBins pBins( cfg, "b.p" );
    if ( 0 == iteration ){
            HistoBins ibBins( cfg, "b.iBetaFirst" ); 
            book->make2D( "inverseBeta" , "1/beta", pBins.nBins(), pBins.getBins().data(), ibBins.nBins(), ibBins.getBins().data() );
    } else {
            HistoBins ibBins( cfg, "b.iBeta" ); 
            book->make2D( "inverseBeta" , "1/beta", pBins.nBins(), pBins.getBins().data(), ibBins.nBins(), ibBins.getBins().data() );
    }
    
    TaskTimer t;
    t.start();

    Int_t nEvents = (Int_t)chain->GetEntries();
    nEventsToProcess = cfg->getInt( nodePath+"DataSource:maxEvents", nEvents );

    if ( nEventsToProcess > nEvents )
        nEventsToProcess = nEvents;
    
    logger->info(__FUNCTION__) << "Loaded: " << nEventsToProcess << " events " << endl;
    
    TaskProgress tp( "Plotting 1/beta", nEventsToProcess );
    logger->info( __FUNCTION__ ) << "Iteration : " << iteration << endl;
    // loop over all events
    for(Int_t i=0; i<nEventsToProcess; i++) {
        chain->GetEntry(i);

        tp.showProgress( i );

        /**
         * Select good events (should already be done in splitter )
         */
        if ( !keepEvent() )
            continue;

        int nTofHits = tuple->nTofHits;
        for ( int iHit = 0; iHit < nTofHits; iHit++ ){

            int tray = tuple->tray[ iHit ];
            int module = tuple->module[ iHit ];
            int cell = tuple->cell[ iHit ];
            int id = relIndex( tray, module, cell );

            if ( id < 0 )
                continue;

            const double zLocal = tuple->zLocal[ iHit ];
            const double tLength = tuple->length[ iHit ];
            const double p = tuple->pt[ iHit ] * TMath::CosH( tuple->eta[ iHit ] );
            const double tofExp = expectedTof( tLength, p );
            const double tot = tuple->tot[ iHit ];;
            
            double rawTof = tuple->leTime[ iHit ] - tuple->tStart;
            double corrTof = corrections[id]->tof( rawTof, tot, zLocal ); // all corrections
            
            double iBeta = (corrTof / tLength )*cLight;

            book->fill( "inverseBeta", p, iBeta );
            
        } // loop on tofHits
    } // end loop on events
    logger->info(__FUNCTION__) << "Completed in " << t.elapsed() << endl;

    reporter->newPage( 1, 1 );
    book->style( "inverseBeta" )->
    set( "title", "#beta^{-1} : " + nameFor( 0 ) + " -> " + nameFor( nElements-1 ) )->
    set( "x", "p [GeV]" )->set( "y", "#beta^{-1}")->
    set( "draw", "colz" )->set( "logz", 1 )->draw();
    reporter->savePage();
    logger->info(__FUNCTION__) << "Report Written" << endl;
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


vector<int> TofCalibrationMaker::fromRelIndex( int id ){

    vector<int> res(3);

    int tm = trayRange->min;
    int mm = moduleRange->min;
    int cm = cellRange->min;

    if ( "tray" == splitMode ){
        res[ 0 ] = id + tm;
        res[ 1 ] = 0;
        res[ 2 ] = 0;
    } else if ( "board" == splitMode ){
        res[ 0 ] = id / 8 + tm;
        res[ 1 ] = (id - ( res[0] - tm) * 8 ) * 4 + mm;
        res[ 2 ] = 0;
    } else if ( "module" == splitMode ){
        res[ 0 ] = id / 32 + tm;
        res[ 1 ] = ( id - (res [ 0 ] - tm) * 32 + mm );
        res[ 3 ] = 0;
    } else if ( "cell" == splitMode ){
        res[ 0 ] = ( id / (32 * 6) ) + tm;
        res[ 1 ] = (id - ( res[ 0 ] - tm)*32*6 ) / 6 + mm;
        res[ 2 ] = id - ( res[ 0 ] - tm )*32*6 - ( res[1] - mm )*6 + cm;
    }
    return res;
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

    for ( int i = 0; i < nElements; i++ ){
        corrections[ i ]->updateTotSpline();    
    }
    
    logger->info( __FUNCTION__ ) << "Complete" << endl;



}


void TofCalibrationMaker::makeBins( string var, int nBins, double min, double max ){


    vector<double> * vals = new vector<double>[ nElements ];

    TaskTimer t;
    t.start();

    Int_t nEvents = (Int_t)chain->GetEntries();
    
    nEventsToProcess = cfg->getInt( nodePath+"DataSource:maxEvents", nEvents );

    if ( nEventsToProcess > nEvents )
        nEventsToProcess = nEvents;
    
    logger->info(__FUNCTION__) << "Loaded: " << nEventsToProcess << " events " << endl;
    
    TaskProgress tp( "Binning " + var , nEventsToProcess );

    // loop over all events
    for(Int_t i=0; i<nEventsToProcess; i++) {
        chain->GetEntry(i);

        tp.showProgress( i );


        /**
         * Select good events
         */
        if ( !keepEvent() )
            continue;

        int nTofHits = tuple->nTofHits;
        for ( int iHit = 0; iHit < nTofHits; iHit++ ){

            if ( !keepTrack( iHit ) )
                continue;

            //cout << " ( " << ds->getInt( "tray", iHit ) << ", " << ds->getInt( "module", iHit ) << ", " << ds->getInt( "cell", iHit ) << " ) " << endl;
            int tray = tuple->tray[ iHit ];
            int module = tuple->module[ iHit ];
            int cell = tuple->cell[ iHit ];
            int id = relIndex( tray, module, cell );

            if ( id < 0 )
                continue;


            double val = 0;
            if ( "tot" == var ) 
                val = tuple->tot[ iHit ];
            if ( "zLocal" == var ) 
                val = tuple->zLocal[ iHit ];

            if ( val < min || val > max )
                continue;

            vals[ id ].push_back( val );

        } // loop on tracks
    } // loop on events


    for ( int i = 0; i < nElements; i++  ){

        Int_t size = vals[i].size();
        Int_t step = size / (nBins + 1 ); 
        
        // sort into ascending order
        std::sort( vals[i].begin(), vals[i].end());
        vector<double> binEdges;
        binEdges.push_back( min );
        

        for( Int_t j = 1; j < nBins ; j++) {

            double d1 = vals[ i ].at( step * j );
            binEdges.push_back( d1 );
            
        }   // loop over tot bins

        binEdges.push_back( max );

        // use the bin edges here
        if ( "tot" == var )
            corrections[ i ]->makeTotBins( binEdges );

    }
    

    delete[] vals;
    logger->info(__FUNCTION__) << "Completed in " << t.elapsed() << endl;

    return;

}



void TofCalibrationMaker::reportTot(){

    logger->info(__FUNCTION__) << endl;
    gStyle->SetOptStat( 1111 );
    book->cd( "totStep_"+ts(iteration) );
    reporter->newPage( 4, 4);
    for ( int i = 0; i < nElements; i++ ){
        if ( !corrections[ i ]->getTotBins() )
            continue;
        //logger->info(__FUNCTION__) << corrections[ i ]->getTotBins()->nBins() << endl;
        double x1 = corrections[ i ]->getTotBins()->getBins()[ 1 ];
        double x2 = corrections[ i ]->getTotBins()->getBins()[ corrections[ i ]->getTotBins()->nBins() - 1 ];
        logger->info(__FUNCTION__) << "Graphing Spline From ( " << x1 << ", " << x2 << " ) " << endl;
        book->style( "tot_"+ts(i) )->
            set( "draw", "colz" )->set( "title", nameFor( i ) )->
            set( "x", "ToT [ns]")->set( "y", "#Delta TOF_{meas} - TOF_{exp}" )->
            draw();
        logger->info(__FUNCTION__) << "Updating Spline" << endl;
        corrections[ i ]->updateTotSpline();
        TGraph * g = corrections[ i ]->getTotSpline()->graph( x1, x2, .5 );
        g->SetLineColor( kRed );
        g->Draw("same");

        reporter->next();
    }
    reporter->savePage();

    /**
     * Draw the corrected ones
     */
    if ( 1 <= iteration ){
        reporter->newPage( 4, 4);
        for ( int i = 0; i < nElements; i++ ){

            book->style( "corrTot_"+ts(i) )->
                set( "draw", "colz" )->set( "title",  nameFor( i ) + " : Corrected ToT" )->
                set( "x", "ToT [ns]")->set( "y", "#Delta TOF_{meas} - TOF_{exp}" )->
                draw();

            reporter->next();
        }
        reporter->savePage();
    }


}

void TofCalibrationMaker::reportZLocal(){

    book->cd( "zLocalStep_"+ts(iteration) );
    reporter->newPage( 3, 3);
    for ( int i = 0; i < nElements; i++ ){

        double x1 = corrections[ i ]->getZBins()->getBins()[ 0 ];
        double x2 = corrections[ i ]->getZBins()->getBins()[ corrections[ i ]->getZBins()->nBins() ];
        logger->info(__FUNCTION__) << "Graphing Spline From ( " << x1 << ", " << x2 << " ) " << endl;
        
        
        book->style( "zLocal_"+ts(i) )->
            set( "draw", "colz" )->set( "title", nameFor( i ) )->
            set( "x", "zLocal [cm]")->set( "y", "#Delta TOF_{meas} - TOF_{exp}" )->
            draw();
        corrections[ i ]->updateTotSpline();
        TGraph * g = corrections[ i ]->getZSpline()->graph( x1, x2, .5 );
        g->SetLineColor( kRed );
        g->Draw("same");
        
        reporter->next();
    }
    reporter->savePage();

    /**
     * Draw the corrected ones
     */
    if ( 1 <= iteration ){
        reporter->newPage( 4, 4);
        for ( int i = 0; i < nElements; i++ ){

            book->style( "corrZLocal_"+ts(i) )->
                set( "draw", "colz" )->set( "title",  nameFor( i ) + " : Corrected zLocal" )->
                set( "x", "zLocal [cm]")->set( "y", "#Delta TOF_{meas} - TOF_{exp}" )->
                draw();

            reporter->next();
        }
        reporter->savePage();
    }


}


void TofCalibrationMaker::exportT0Params( string pFile ) {

    ofstream params( pFile.c_str() );

    for ( int iTray = 1; iTray <= nTrays; iTray++ ){
        for ( int iMod = 1; iMod <= nModules; iMod++ ){
            for ( int iCell = 1; iCell <= nCells; iCell++ ){
                int id = relIndex( iTray, iMod, iCell );
                if ( id < 0 )
                    continue;

                params << iTray << "\t" << iMod << "\t" << iCell << endl;
                params << corrections[ id ]->getT0() << endl;

            }
        }
    }

    params.close();

}

void TofCalibrationMaker::exportTotParams( string pFile ) {

    ofstream params( pFile.c_str() );

    if ( trayRange->min == 1 ){
        if ( "cell" == splitMode )
           params << (nTrays * nModules * nCells) << endl;
        else if ( "module" == splitMode )
           params << (nTrays * nModules ) << endl;
        else if ( "board" == splitMode )
           params << (nTrays * 8) << endl;
    }
    for ( int iTray = 1; iTray <= nTrays; iTray++ ){
        for ( int iMod = 1; iMod <= nModules; iMod++ ){
            for ( int iCell = 1; iCell <= nCells; iCell++ ){
                int id = relIndex( iTray, iMod, iCell );
                if ( id < 0 )
                    continue;

                params << iTray << "\t" << iMod << "\t" << iCell << endl;
                int nBins = corrections[ id ]->getTotBins()->nBins();
                params << nBins-1 << endl;
                for ( int i = 0; i < nBins; i++ )
                    params << corrections[ id ]->getTotBins()->getBins()[ i ] << " ";
                params << endl;
                for ( int i = 0; i < nBins; i++ )
                    params << corrections[ id ]->getTot( i, true ) << " ";
                params << endl;

            }
        }
    }

    params.close();

}

void TofCalibrationMaker::exportZParams( string pFile ) {

    ofstream params( pFile.c_str() );

    if ( trayRange->min == 1 ){
        if ( "cell" == splitMode )
           params << (nTrays * nModules * nCells) << endl;
        else if ( "module" == splitMode )
           params << (nTrays * nModules ) << endl;
        else if ( "board" == splitMode )
           params << (nTrays * 8) << endl;
    }
    for ( int iTray = 1; iTray <= nTrays; iTray++ ){
        for ( int iMod = 1; iMod <= nModules; iMod++ ){
            for ( int iCell = 1; iCell <= nCells; iCell++ ){
                int id = relIndex( iTray, iMod, iCell );
                if ( id < 0 )
                    continue;

                params << iTray << "\t" << iMod << "\t" << iCell << endl;
                int nBins = corrections[ id ]->getZBins()->nBins();
                params << nBins-1 << endl;
                for ( int i = 0; i < nBins; i++ )
                    params << corrections[ id ]->getZBins()->getBins()[ i ] << " ";
                params << endl;
                for ( int i = 0; i < nBins; i++ )
                    params << corrections[ id ]->getZ( i, true ) << " ";
                params << endl;

            }
        }
    }

    params.close();

}










/*
tof->Draw("leTime - tStart>>h1(10, -50, 50)", "numberOfVpdEast>0 && numberOfVpdWest>0 && abs(vertexZ-vpdVz) < 6.0 && TMath::Sqrt(vertexX*vertexX + vertexY*vertexY )< 1.0 && nHitsFit > 25& pt*TMath::CosH( eta ) < 1.0 & pt*TMath::CosH( eta ) > 0.3 && dedx < 2.8")
 */











