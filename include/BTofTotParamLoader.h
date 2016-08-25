#ifndef BTOF_TOT_PARAM_LOADER_H
#define BTOF_TOT_PARAM_LOADER_H

class BTofTotParamLoader
{
public:
	BTofTotParamLoader() {}
	~BTofTotParamLoader() {}


	void loadParams( string filename, float * _binEdges, float * _binContents ){

		ifstram infile;
		infile.open( filename.c_str() );

		int calibType = -1;
		infile >> calibType


		if ( 23040 == calibType ){
			loadParams_Cell_by_Cell( infile, _binEdges, _binContents );
		} else if ( 840 == calibType ){

		} else if ( 960 == calibType ){

		}
		return ;
	}

	void loadParams_Cell_by_Cell( ifstream &infile, float * _binEdges, float * _binContents ){

		int mNModule = 32;      // 32 for tofr5++ 
		int mNCell = 6;         // 6 cells per module
		int mNBinMax = 60;      // 60 bins for T-Tot, T-Z correction


		for(int i=0;i<mNTray;i++) {
			for(int j=0;j<mNModule;j++) {
				for(int l=0;l<mNCell;l++){
					
					int nbin = 0;
					int trayId = -1, moduleId = -1, cellId = -1;
					infile >> trayId >> moduleId >> cellId;
					infile >> nbin;
					
					if ( trayId < 1 || moduleId < 1 || cellId < 1 ) continue;

					for(int k = 0; k <= nbin; k++ ) 
						infile >> _binEdges[trayId-1][moduleId-1][cellId-1][k];
					
					for(int k = 0; k <= nbin; k++ ) {
						infile >> _binContents[trayId-1][moduleId-1][cellId-1][k];
					}
				}//cell
			}//module
		}//tray
	}
	
};




#endif