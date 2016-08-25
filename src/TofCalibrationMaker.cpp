#include "TofCalibrationMaker.h"
#include "TProfile.h"


const int TofCalibrationMaker::nTrays   = 120;
const int TofCalibrationMaker::nModules = 32;
const int TofCalibrationMaker::nCells   = 6;

const double TofCalibrationMaker::cLight = 29.9792458; //cm / ns
const double TofCalibrationMaker::mPi = 0.13957; //in GeV / c^2


void TofCalibrationMaker::init( XmlConfig &_config, string _nodePath, int _jobIndex ){
	TreeAnalyzer::init( _config, _nodePath, _jobIndex );


	splitMode   = config.getString( nodePath + ":splitMode", "tray" );
	trayRange   = unique_ptr<XmlRange>( new XmlRange( &config, nodePath + ".Trays", 1, 120 ) );
	moduleRange = unique_ptr<XmlRange>( new XmlRange( &config, nodePath + ".Modules", 1, 32 ) );
	cellRange   = unique_ptr<XmlRange>( new XmlRange( &config, nodePath + ".Cells", 1, 6 ) );
    
    INFO( classname(), "SplitMode "    << splitMode );
    INFO( classname(), "Tray Range "   << trayRange->toString() );
    INFO( classname(), "Module Range " << moduleRange->toString() );
    INFO( classname(), "Cell Range "   << cellRange->toString() );

    int nTraysToProcess = trayRange->max - trayRange->min + 1;
    int nModsToProcess  = moduleRange->max - moduleRange->min + 1;
    int nCellsToProcess = cellRange->max - cellRange->min + 1;

    nElements = nTraysToProcess;
    if ( "module" == splitMode )
        nElements = nTraysToProcess * nModsToProcess;
    else if ( "cell" == splitMode )
        nElements = nTraysToProcess * nModsToProcess * nCellsToProcess;
    else if ( "board" == splitMode )
        nElements = nTraysToProcess * 8;

    INFO( classname(), "Calibrating on " << nElements << " elements" );

    
    HistoBins totBins( config, nodePath + ".Bins.tot" );
    HistoBins zBins( config, nodePath + ".Bins.zLocal" );
    for ( int i = 0; i < nElements; i++ ){
    	INFO( classname(), "Creating corrections for" << i );
        corrections.push_back( unique_ptr<TofCorrection>( 
            
            new TofCorrection( 	totBins.getBins(), 
                                zBins.getBins(), 
                                config.getBool( nodePath + ".Spline:tot", false ),
                                config.getBool( nodePath + ".Spline:zLocal", false ) ) 
            ) 
        );
    }
    iteration = 0;
    INFO( classname(), "Setup null corrections for zLocal and tot" );

    if ( config.getBool( nodePath + ".Bins.tot:variable", false ) )
        makeBins( "tot", config.getDouble( nodePath + ".Bins.tot:nBins" ), config.getDouble( nodePath + ".Bins.tot:min" ), config.getDouble( nodePath + ".Bins.tot:max" ) );
    if ( config.getBool( nodePath + ".Bins.zLocal:variable", false ) )
        makeBins( "zLocal", config.getDouble( nodePath + ".Bins.zLocal:nBins" ), config.getDouble( nodePath + ".Bins.zLocal:min" ), config.getDouble( nodePath + ".Bins.zLocal:max" ) );
    
    INFO( classname(), corrections[ 0 ]->getZBins()->toString() );
    INFO( classname(), vts( corrections[ 0 ]->getZBins()->bins ) );
    INFO( classname(), "Done Setting up" );
    INFO( classname(), "Done Setting up" );
    INFO( classname(), "Done Setting up" );
    INFO( classname(), "Done Setting up" );
}


void TofCalibrationMaker::fillTot(  ){

    if ( !ds ){
        ERROR( classname(), "Invalid DataSource " );
        return;
    }


    /**
     * Make histos
     */
    
    book->cd( "totStep_"+ts(iteration) );
    HistoBins dtBins( config, "b.dt" );
    for ( int i = 0; i < nElements; i++ ){
        HistoBins * hb = corrections[ i ]->getTotBins();
        TH2D * h = new TH2D( ("tot_" + ts( i )).c_str(), "dt vs. Tot", hb->nBins(), hb->getBins().data(), dtBins.nBins(), dtBins.getBins().data() );
        TH2D * hCorr = new TH2D( ("corrTot_" + ts( i )).c_str(), "dt vs. Tot", hb->nBins(), hb->getBins().data(), dtBins.nBins(), dtBins.getBins().data() );
        
        book->add( "tot_" + ts( i ), h );
        book->add( "corrTot_" + ts( i ), hCorr );
    }


    TF1 * piBetaCut = new TF1( "piBetaCut", "TMath::Sqrt([0]*[0]+x*x)/x",0,5);
    piBetaCut->SetParameter(0, 0.4 );

    TaskTimer t;
    t.start();

    Int_t nEvents = (Int_t)chain->GetEntries();
    
    nEventsToProcess = config.getInt( nodePath + ".DataSource:maxEvents", nEvents );

    if ( nEventsToProcess > nEvents )
        nEventsToProcess = nEvents;
    
    INFO( classname(), "Loaded: " << nEventsToProcess << " events " );
    
    TaskProgress tp( "Plotting Tot", nEventsToProcess );
    INFO( classname(), "Iteration : " << iteration );
    // loop over all events
    for(Int_t i=0; i<nEventsToProcess; i++) {
        // ds->getEntry(i);
        chain->GetEntry(i);

        tp.showProgress( i );

        /**
         * Select good events (should already be done in splitter )
         */
        if ( !keepEvent() )
            continue;

        int nTofHits = ds->get( "nTofHits" );
        for ( int iHit = 0; iHit < nTofHits; iHit++ ){

            if ( !keepTrack( iHit ) )
                continue;

            TRACE( classname(), "good track" );

            int tray   = ds->get( "tray", iHit );
            int module = ds->get( "module", iHit );
            int cell   = ds->get( "cell", iHit );
            int id     = relIndex( tray, module, cell );

            if ( id < 0 )
                continue;

            TRACE( classname(), "tray=" << tray << ", mod=" << module << ", cell=" << cell );
            TRACE( classname(), "id = " << id );

            const double tLength = ds->get( "length", iHit );
            const double p       = ds->get( "pt", iHit ) * TMath::CosH( ds->get( "eta", iHit ) );
            const double tofExp  = expectedTof( tLength, p );
            const double tot     = ds->get( "tot", iHit );
            
            double rawTof        = ds->get( "leTime", iHit ) - ds->get( "tStart" );
            double corrTof       = corrections[id]->tof( rawTof, tot, 0 ); // all corrections
            double tof           = corrections[id]->tofForTot( rawTof, 0 ); // aplly all corrections except tot
            
            double dt = tof - tofExp;
            
            double iBeta = (corrTof / tLength ) * cLight;
            // cut on the inverse beta curve
            // if ( 0 < iteration && (iBeta > piBetaCut->Eval( p ) ) ) continue;   

            book->fill( "tot_" + ts(id), ds->get( "tot", iHit ), dt );
            book->fill( "corrTot_" + ts(id), ds->get( "tot", iHit ), (corrTof - tofExp) );

        } // loop on tofHits
        
    } // end loop on events
    INFO( classname(), "Completed in " << t.elapsed() );
}

void TofCalibrationMaker::fillZLocal(  ){

    if ( !ds ){
        INFO( classname(), "Invalid DataSource " );
        return;
    }


    /**
     * Make histos
     */
    
    book->cd( "zLocalStep_"+ts(iteration) );
    HistoBins dtBins( config, "b.dtZ" );
    for ( int i = 0; i < nElements; i++ ){
        HistoBins * hb = corrections[ i ]->getZBins();
        TH2D * h = new TH2D( ("zLocal_" + ts( i )).c_str(), "dt vs. zLocal", hb->nBins(), hb->getBins().data(), dtBins.nBins(), dtBins.getBins().data() );
        book->add( "zLocal_" + ts( i ), h );
        TH2D * hh = new TH2D( ("corrZLocal_" + ts( i )).c_str(), "dt vs. zLocal", hb->nBins(), hb->getBins().data(), dtBins.nBins(), dtBins.getBins().data() );
        book->add( "corrZLocal_" + ts( i ), hh );
    }

    
    TF1 * piBetaCut = new TF1( "piBetaCut", "TMath::Sqrt([0]*[0]+x*x)/x",0,5);
    piBetaCut->SetParameter(0, 0.4 );

    TaskTimer t;
    t.start();

    Int_t nEvents = (Int_t)chain->GetEntries();
    
    nEventsToProcess = config.getInt( nodePath+".DataSource:maxEvents", nEvents );

    if ( nEventsToProcess > nEvents )
        nEventsToProcess = nEvents;
    
    INFO( classname(), "Loaded: " << nEventsToProcess << " events " );
    
    TaskProgress tp( "Event Loop", nEventsToProcess );
    INFO( classname(), "Iteration : " << iteration );
    // loop over all events
    for(Int_t i=0; i<nEventsToProcess; i++) {
        ds->getEntry(i);

        tp.showProgress( i );

        /**
         * Select good events (should already be done in splitter )
         */
        if ( !keepEvent() )
            continue;

        int nTofHits = ds->get( "nTofHits" );
        for ( int iHit = 0; iHit < nTofHits; iHit++ ){

            if ( !keepTrack( iHit ) )
                continue;

            int tray = ds->get( "tray", iHit );
            int module = ds->get( "module", iHit );
            int cell = ds->get( "cell", iHit );
            int id = relIndex( tray, module, cell );

            if ( id < 0 )
                continue;

            const double zLocal = ds->get( "zLocal", iHit  );
            const double tLength = ds->get( "length", iHit );
            const double p = ds->get( "pt", iHit ) * TMath::CosH( ds->get( "eta", iHit ) );
            const double tofExp = expectedTof( tLength, p );
            const double tot = ds->get( "tot", iHit );
            
            const double rawTof   = ds->get( "leTime", iHit ) - ds->get( "tStart" );
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
    INFO( classname(), "Completed in " << t.elapsed() );


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

        double x1 = config.getDouble( nodePath + ".Bins.zLocal:min" );
        double x2 = config.getDouble( nodePath + ".Bins.zLocal:max" );

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
            if ( !config.getBool( nodePath+".Spline:zLocalPol1", false ) )
                corrections[ i ]->setZ( j, profile->GetBinContent( j+1 ) );
            else {
                zThere = x1 + j * (( x2 - x1 ) / (double)l);
                cout << " zThere " << zThere << endl;
                corrections[ i ]->setZ( j, pol1->Eval( zThere ) );
            }
        }
        corrections[ i ]->updateZSpline();
    } 
    INFO( classname(), "Completed in " << t.elapsed() );


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
    INFO( classname(), "Completed in " << t.elapsed() );
}


void TofCalibrationMaker::inverseBeta(  ){

    if ( !ds ){
        ERROR( classname(), "Invalid DataSource " );
        return;
    }

    /**
     * Make histos
     */
    book->cd( "iBetaStep_"+ts(iteration) );
    HistoBins pBins( config, "b.p" );
    if ( 0 == iteration ){
            HistoBins ibBins( config, "b.iBetaFirst" ); 
            TH2D * h = new TH2D( "inverseBeta" , "1/beta", pBins.nBins(), pBins.getBins().data(), ibBins.nBins(), ibBins.getBins().data() );
            book->add( "inverseBeta", h );
    } else {
            HistoBins ibBins( config, "b.iBeta" ); 
            TH2D * h = new TH2D( "inverseBeta" , "1/beta", pBins.nBins(), pBins.getBins().data(), ibBins.nBins(), ibBins.getBins().data() );
            book->add( "inverseBeta", h );

    }
    
    TaskTimer t;
    t.start();

    Int_t nEvents = (Int_t)chain->GetEntries();
    nEventsToProcess = config.getInt( nodePath+".DataSource:maxEvents", nEvents );

    if ( nEventsToProcess > nEvents )
        nEventsToProcess = nEvents;
    
    INFO( classname(), "Loaded: " << nEventsToProcess << " events " );
    
    TaskProgress tp( "Plotting 1/beta", nEventsToProcess );
    INFO( classname(), "Iteration : " << iteration );
    // loop over all events
    for(Int_t i=0; i<nEventsToProcess; i++) {
        ds->getEntry(i);

        tp.showProgress( i );

        /**
         * Select good events (should already be done in splitter )
         */
        if ( !keepEvent() )
            continue;

        int nTofHits = ds->get( "nTofHits" );
        for ( int iHit = 0; iHit < nTofHits; iHit++ ){

            int tray   = ds->get( "tray", iHit );
            int module = ds->get( "module", iHit );
            int cell   = ds->get( "cell", iHit );
            int id     = relIndex( tray, module, cell );

            if ( id < 0 )
                continue;

            const double zLocal  = ds->get( "zLocal", iHit );
            const double tLength = ds->get( "length", iHit );
            const double p       = ds->get( "pt", iHit ) * TMath::CosH( ds->get( "eta", iHit ) );
            const double tofExp  = expectedTof( tLength, p );
            const double tot     = ds->get( "tot", iHit );
            
            double rawTof        = ds->get( "leTime", iHit ) - ds->get( "tStart" );
            double corrTof       = corrections[id]->tof( rawTof, tot, zLocal ); // all corrections
            
            double iBeta         = (corrTof / tLength )*cLight;

            book->fill( "inverseBeta", p, iBeta );
            
        } // loop on tofHits
    } // end loop on events
    INFO( classname(), "Completed in " << t.elapsed() );

    // reporter->newPage( 1, 1 );
    // book->style( "inverseBeta" )->
    // set( "title", "#beta^{-1} : " + nameFor( 0 ) + " -> " + nameFor( nElements-1 ) )->
    // set( "x", "p [GeV]" )->set( "y", "#beta^{-1}")->
    // set( "draw", "colz" )->set( "logz", 1 )->draw();
    // reporter->savePage();
}


void TofCalibrationMaker::alignT0(){

    if ( !ds ){
        ERROR( classname(), "Invalid DataSource " );
        return;
    }

    /**
     * Make histos
     */
    book->cd( "t0Step_"+ts(iteration) );

    if ( 0 == iteration ){
    	TH2D * h = new TH2D( "elementT0", "t0", nElements, -0.5, nElements-0.5, 200, -100, 100 );
       	book->add( "elementT0", h );
    } else if ( 1 == iteration ){
        TH2D * h = new TH2D( "elementT0", "t0", nElements, -0.5, nElements-0.5, 200, -10, 10 );
       	book->add( "elementT0", h );
    } else {
        TH2D * h = new TH2D( "elementT0", "t0", nElements, -0.5, nElements-0.5, 2000, -5, 5 );
       	book->add( "elementT0", h );
    }
    



    TF1 * piBetaCut = new TF1( "piBetaCut", "TMath::Sqrt([0]*[0]+x*x)/x",0,5);
    piBetaCut->SetParameter(0, 0.4 );

    TaskTimer t;
    t.start();

    Int_t nEvents = (Int_t)chain->GetEntries();
    
    nEventsToProcess = config.getInt( nodePath+".DataSource:maxEvents", nEvents );

    if ( nEventsToProcess > nEvents )
        nEventsToProcess = nEvents;
    
    INFO( classname(), "Loaded: " << nEventsToProcess << " events " );
    
    TaskProgress tp( "T0 Alignment", nEventsToProcess );
    INFO( classname(), "Iteration : " << iteration );
    // loop over all events
    for(Int_t i=0; i<nEventsToProcess; i++) {
        ds->getEntry(i);

        tp.showProgress( i );


        /**
         * Select good events
         */
        if ( !keepEvent() )
            continue;

        int nTofHits = ds->get( "nTofHits" );
        for ( int iHit = 0; iHit < nTofHits; iHit++ ){

            if ( !keepTrack( iHit ) )
                continue;

            
            int tray = ds->get( "tray", iHit );
            int module = ds->get( "module", iHit );
            int cell = ds->get( "cell", iHit );
            int id = relIndex( tray, module, cell );

            if ( id < 0 )
                continue;



            const double zLocal     = ds->get( "zLocal", iHit );
            const double tLength    = ds->get( "length", iHit );
            const double p          = ds->get( "pt", iHit ) * TMath::CosH( ds->get( "eta", iHit ) );
            const double tofExp     = expectedTof( tLength, p );
            const double tot        = ds->get( "tot", iHit );
            
            const double rawTof = ds->get( "leTime", iHit ) - ds->get( "tStart" );
            const double corrTof = corrections[id]->tof( rawTof, tot, 0 ); // all corrections
            const double iBeta = ( corrTof / tLength ) * cLight;
    
            double dt = corrTof - tofExp;

            // cut on the inverse beta curve
            if ( 0 < iteration && (iBeta > piBetaCut->Eval( p ) ) ) continue;
            
            //book->fill( "iBeta", p, iBeta );    
            book->fill( "elementT0", id, dt );

        } // loop on tofHits
    } // end loop on events
    INFO( classname(), "Completed in " << t.elapsed() );

    
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


    // reporter->newPage();

    // book->style( "elementT0" )->
    // set( "title", "T0 : " + nameFor( 0 ) + " -> " + nameFor( nElements-1 ) )->
    // set( "y", "#Delta TOF_{m} - TOF_{exp}")->
    // set( "draw", "colz" )->set( "logz", 1 )->draw();

    // reporter->savePage();

}




bool TofCalibrationMaker::keepEvent( ){
    if ( 0 >= ds->get( "numberOfVpdEast" ) || 0 >= ds->get( "numberOfVpdWest" ) ) 
        return false;
    if ( ds->get( "vR" ) > 1.0 ) 
        return false;
    if ( TMath::Abs( ds->get( "vpdVz" )  - ds->get("vertexZ") ) > 6.0 ) 
        return false;
    return true;
}
bool TofCalibrationMaker::keepTrack( int iHit ){
    double p = ds->get( "pt", iHit ) * TMath::CosH( ds->get( "eta", iHit ) );
    double nSigPi = ds->get("nSigPi", iHit );

    if ( 0.3 > p || 0.6 < p ) 
        return false;
    if ( nSigPi > 2.0 ) 
        return false;
    if ( 25 > ds->get( "nHitsFit", iHit ) ) 
        return false;
    return true;
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
        ERROR( classname(), "" );
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
        ERROR( classname(), "" );
    }
    DEBUG( classname(), "( " << tray << ", " << module << ", " << cell << " ) " );
    DEBUG( classname(), trayRange->min << " -> " << trayRange->max );

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


void TofCalibrationMaker::makeBins( string var, int nBins, double min, double max ){


    vector<double> * vals = new vector<double>[ nElements ];

    TaskTimer t;
    t.start();

    Int_t nEvents = (Int_t)chain->GetEntries();
    
    nEventsToProcess = config.getInt( nodePath+".DataSource:maxEvents", nEvents );

    if ( nEventsToProcess > nEvents )
        nEventsToProcess = nEvents;
    
    INFO( classname(), "Loaded: " << nEventsToProcess << " events " );
    
    TaskProgress tp( "Binning " + var , nEventsToProcess );

    // loop over all events
    for(Int_t i=0; i<nEventsToProcess; i++) {
        ds->getEntry(i);

        tp.showProgress( i );


        /**
         * Select good events
         */
        if ( !keepEvent() )
            continue;

        int nTofHits = ds->get( "nTofHits" );
        for ( int iHit = 0; iHit < nTofHits; iHit++ ){

            if ( !keepTrack( iHit ) )
                continue;

            //cout << " ( " << ds->get( "tray", iHit ) << ", " << ds->get( "module", iHit ) << ", " << ds->get( "cell", iHit ) << " ) " << endl;
            int tray   = ds->get( "tray", iHit );
            int module = ds->get( "module", iHit );
            int cell   = ds->get( "cell", iHit );
            int id     = relIndex( tray, module, cell );

            if ( id < 0 )
                continue;

            double val = ds->get( var, iHit );

            if ( val < min || val > max )
                continue;

            vals[ id ].push_back( val );

        } // loop on tracks
    } // loop on events


    for ( int i = 0; i < nElements; i++  ){

        INFO( classname(), "Making Quantiled " << var << "bins from " << vals[i].size() << " datapoints in " << nBins << " from " << min << " to " << max  );
        vector<double> binEdges = HistoBins::makeQuantileBins( vals[i], nBins, min, max );
        INFO( classname(), var << " bins : " << vts( binEdges ) );

        // use the bin edges here
        if ( "tot" == var )
            corrections[ i ]->makeTotBins( binEdges );

    }
    

    delete[] vals;
    INFO( classname(), "Completed in " << t.elapsed() );

    return;

}

void TofCalibrationMaker::reportTot(){

    // gStyle->SetOptStat( 1111 );
    // book->cd( "totStep_"+ts(iteration) );
    // reporter->newPage( 4, 4);
    // for ( int i = 0; i < nElements; i++ ){
    //     if ( !corrections[ i ]->getTotBins() )
    //         continue;
    //     //logger->info(__FUNCTION__) << corrections[ i ]->getTotBins()->nBins() << endl;
    //     double x1 = corrections[ i ]->getTotBins()->getBins()[ 1 ];
    //     double x2 = corrections[ i ]->getTotBins()->getBins()[ corrections[ i ]->getTotBins()->nBins() - 1 ];
    //     INFO( classname(), "Graphing Spline From ( " << x1 << ", " << x2 << " ) " );
    //     book->style( "tot_"+ts(i) )->
    //         set( "draw", "colz" )->set( "title", nameFor( i ) )->
    //         set( "x", "ToT [ns]")->set( "y", "#Delta TOF_{meas} - TOF_{exp}" )->
    //         draw();
    //     INFO( classname(), "Updating Spline" );
    //     corrections[ i ]->updateTotSpline();
    //     TGraph * g = corrections[ i ]->getTotSpline()->graph( x1, x2, .5 );
    //     g->SetLineColor( kRed );
    //     g->Draw("same");

    //     reporter->next();
    // }
    // reporter->savePage();

    // /**
    //  * Draw the corrected ones
    //  */
    // if ( 1 <= iteration ){
    //     reporter->newPage( 4, 4);
    //     for ( int i = 0; i < nElements; i++ ){

    //         book->style( "corrTot_"+ts(i) )->
    //             set( "draw", "colz" )->set( "title",  nameFor( i ) + " : Corrected ToT" )->
    //             set( "x", "ToT [ns]")->set( "y", "#Delta TOF_{meas} - TOF_{exp}" )->
    //             draw();

    //         reporter->next();
    //     }
    //     reporter->savePage();
    // }
} // report Tot

void TofCalibrationMaker::reportZLocal(){

    // book->cd( "zLocalStep_"+ts(iteration) );
    // reporter->newPage( 3, 3);
    // for ( int i = 0; i < nElements; i++ ){

    //     double x1 = corrections[ i ]->getZBins()->getBins()[ 0 ];
    //     double x2 = corrections[ i ]->getZBins()->getBins()[ corrections[ i ]->getZBins()->nBins() ];
    //     INFO( classname(), "Graphing Spline From ( " << x1 << ", " << x2 << " ) " );
        
        
    //     book->style( "zLocal_"+ts(i) )->
    //         set( "draw", "colz" )->set( "title", nameFor( i ) )->
    //         set( "x", "zLocal [cm]")->set( "y", "#Delta TOF_{meas} - TOF_{exp}" )->
    //         draw();
    //     corrections[ i ]->updateTotSpline();
    //     TGraph * g = corrections[ i ]->getZSpline()->graph( x1, x2, .5 );
    //     g->SetLineColor( kRed );
    //     g->Draw("same");
        
    //     reporter->next();
    // }
    // reporter->savePage();

    // /**
    //  * Draw the corrected ones
    //  */
    // if ( 1 <= iteration ){
    //     reporter->newPage( 4, 4);
    //     for ( int i = 0; i < nElements; i++ ){

    //         book->style( "corrZLocal_"+ts(i) )->
    //             set( "draw", "colz" )->set( "title",  nameFor( i ) + " : Corrected zLocal" )->
    //             set( "x", "zLocal [cm]")->set( "y", "#Delta TOF_{meas} - TOF_{exp}" )->
    //             draw();

    //         reporter->next();
    //     }
    //     reporter->savePage();
    // }
}// report zLocal




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

    // header
    if ( "cell" == splitMode )
       params << (nTrays * nModules * nCells) << endl;
    else if ( "module" == splitMode )
       params << (nTrays * nModules ) << endl;
    else if ( "board" == splitMode )
       params << (nTrays * 8) << endl;
    

    for ( int iTray = 1; iTray <= nTrays; iTray++ ){
        for ( int iMod = 1; iMod <= nModules; iMod++ ){
            for ( int iCell = 1; iCell <= nCells; iCell++ ){
                int id = relIndex( iTray, iMod, iCell );
                if ( id < 0 )
                    continue;

                params << iTray << "\t" << iMod << "\t" << iCell << endl;
                int nBins = corrections[ id ]->getTotBins()->nBins();
                params << nBins << endl;
                for ( int i = 0; i <= nBins; i++ )
                    params << corrections[ id ]->getTotBins()->getBins()[ i ] << " ";
                params << endl;
                for ( int i = 0; i < nBins; i++ )
                    params << corrections[ id ]->getTot( i, true ) << " ";
                params << corrections[ id ]->getTot( nBins-1, true ); // add the last one again for the interpolation
                params << endl;

            }
        }
    }

    params.close();

}

void TofCalibrationMaker::exportZParams( string pFile ) {

    ofstream params( pFile.c_str() );

    // header
    if ( "cell" == splitMode )
       params << (nTrays * nModules * nCells) << endl;
    else if ( "module" == splitMode )
       params << (nTrays * nModules ) << endl;
    else if ( "board" == splitMode )
       params << (nTrays * 8) << endl;
    for ( int iTray = 1; iTray <= nTrays; iTray++ ){
        for ( int iMod = 1; iMod <= nModules; iMod++ ){
            for ( int iCell = 1; iCell <= nCells; iCell++ ){
                int id = relIndex( iTray, iMod, iCell );
                if ( id < 0 )
                    continue;

                params << iTray << "\t" << iMod << "\t" << iCell << endl;
                int nBins = corrections[ id ]->getZBins()->nBins();
                params << nBins << endl;
                for ( int i = 0; i <= nBins; i++ )
                    params << corrections[ id ]->getZBins()->getBins()[ i ] << " ";
                params << endl;
                for ( int i = 0; i < nBins; i++ )
                    params << corrections[ id ]->getZ( i, true ) << " ";
                params << corrections[ id ]->getZ( nBins-1, true );
                params << endl;

            }
        }
    }

    params.close();

}