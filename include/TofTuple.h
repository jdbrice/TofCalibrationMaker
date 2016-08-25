//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar 10 16:28:04 2015 by ROOT version 5.34/25
// from TChain tof/
//////////////////////////////////////////////////////////

#ifndef TOF_TUPLE_H
#define TOF_TUPLE_H



#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"

#include "TofCellData.h"
// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class TofTuple {
public:
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         vertexX;
   Float_t         vertexY;
   Float_t         vertexZ;
   Int_t           numberOfVpdEast;
   Int_t           numberOfVpdWest;
   Double_t        tStart;
   Float_t         vpdVz;
   Int_t           nTofHits;
   Int_t           tray[kMaxHits];   //[nTofHits]
   Int_t           module[kMaxHits];   //[nTofHits]
   Int_t           cell[kMaxHits];   //[nTofHits]
   Double_t        leTime[kMaxHits];   //[nTofHits]
   Double_t        tot[kMaxHits];   //[nTofHits]
   Int_t           matchFlag[kMaxHits];   //[nTofHits]
   Float_t         zLocal[kMaxHits];   //[nTofHits]
   Int_t           charge[kMaxHits];   //[nTofHits]
   Float_t         pt[kMaxHits];   //[nTofHits]
   Float_t         eta[kMaxHits];   //[nTofHits]
   Float_t         length[kMaxHits];   //[nTofHits]
   Int_t           nHitsFit[kMaxHits];   //[nTofHits]
   Int_t           nHitsDedx[kMaxHits];   //[nTofHits]
   Float_t         dedx[kMaxHits];   //[nTofHits]
   Float_t         nSigE[kMaxHits];   //[nTofHits]
   Float_t         nSigPi[kMaxHits];   //[nTofHits]
   Float_t         nSigK[kMaxHits];   //[nTofHits]
   Float_t         nSigP[kMaxHits];   //[nTofHits]

   // List of branches
   TBranch        *b_vertexX;   //!
   TBranch        *b_vertexY;   //!
   TBranch        *b_vertexZ;   //!
   TBranch        *b_numberOfVpdEast;   //!
   TBranch        *b_numberOfVpdWest;   //!
   TBranch        *b_tStart;   //!
   TBranch        *b_vpdVz;   //!
   TBranch        *b_nTofHits;   //!
   TBranch        *b_tray;   //!
   TBranch        *b_module;   //!
   TBranch        *b_cell;   //!
   TBranch        *b_leTime;   //!
   TBranch        *b_tot;   //!
   TBranch        *b_matchFlag;   //!
   TBranch        *b_zLocal;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_length;   //!
   TBranch        *b_nHitsFit;   //!
   TBranch        *b_nHitsDedx;   //!
   TBranch        *b_dedx;   //!
   TBranch        *b_nSigE;   //!
   TBranch        *b_nSigPi;   //!
   TBranch        *b_nSigK;   //!
   TBranch        *b_nSigP;   //!

   TofTuple(TTree *tree=0);
   virtual ~TofTuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

