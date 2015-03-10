

#include "TreeAnalyzer.h"
#include "TofCalibrationMaker.h"
#include "PicoSplitter.h"
#include "SplineMaker.h"
using namespace jdb;

#include <iostream>
#include <exception>
#include <cstdlib>

int main( int argc, char* argv[] ) {

	if ( argc >= 2 ){

		string fileList = "";
		string jobPrefix = "";
  
		// try block to load the XML config if given
		try{
			XmlConfig config( argv[ 1 ] );
			
			// if more arguments are given they will be used for 
			// parallel job running with examples of
			// argv[ 2 ] = fileList.lis
			// argv[ 3 ] = jobPrefix_ 
			int jobId = -1;
			if ( argc >= 5){
				jobId = atoi(argv[ 2 ]);
				fileList = (string) argv[ 3 ];
				jobPrefix = (string) argv[ 4 ];
			}
			
			if ( "DataSource" == config[ "job" ] ){
				TreeAnalyzer *ta = new TreeAnalyzer( &config, "" );
				ta->make();
				ta = new TreeAnalyzer( &config, "" );
				ta->make();
			} else if ( "calibrateTof" == config[ "job" ] ) {
				TofCalibrationMaker tcm( &config, "", jobId, fileList, jobPrefix );
				tcm.make();

			} else if ( "splitPicos" == config[ "job" ] ) {
				PicoSplitter ps( &config, "" );
				ps.make();
			}

		} catch ( exception &e ){
			cout << e.what() << endl;
		}

	}

/*
	vector<double> bins;
	bins.push_back( 0 );
	bins.push_back( 5 );
	bins.push_back( 10 );
	vector<double> vals;
	vals.push_back( 1 );
	vals.push_back( 2 );
	vals.push_back( 3 );

	Reporter rp( "test.pdf" );

	rp.newPage();
	TH1D * h1 = new TH1D( "h", "h", 10, 0, 10 );
	h1->Draw();
	h1->GetYaxis()->SetRangeUser( 0, 5 );
	SplineMaker sm( bins, vals );
	TGraph * graph = sm.graph( -10, 20, 1 );
	graph->Draw();
	rp.savePage();
*/

	return 0;
}
