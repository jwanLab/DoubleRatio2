/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* AliAnaysisTaskMyTask
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */




#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <vector>
#include <algorithm>
#include <random> // std::default_random_engine  
#include <chrono> // std::chrono::system_clock 
// ROOT classes
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TList.h"
#include "TH1I.h"
#include "TH3I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TExMap.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TComplex.h"
#include "AliAnalysisTask.h"
// Alice analysis base class
#include "AliAnalysisTaskSE.h"
// Alice analysis additional classes
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
// Alice AOD classes
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODVZERO.h"
#include "AliAODZDC.h"

#include "AliAODv0.h"
#include "AliAODTrack.h"
#include "AliAODHeader.h"

// Alice classes
#include "AliCentrality.h"
#include "AliEventplane.h"
//#include "AliEventCuts.h"
#include "AliAnalysisUtils.h"
// Alice MC classes
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliAODMCParticle.h"
// Alice "V" classes
#include "AliVParticle.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliVVZERO.h"
// Alice PID classes
#include "AliAODPid.h"
#include "AliAODpidUtil.h"
#include "AliPID.h"
#include "AliPIDCombined.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"


#include "AliAnalysisTaskMyTask.h"
class AliAnalysisTaskMyTask;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskMyTask) // classimp: necessary for root

AliAnalysisTaskMyTask::AliAnalysisTaskMyTask() : AliAnalysisTaskSE(), 
    fAOD(nullptr), 
    fOutputList(nullptr), 
    mHarmonic(2.),
    fFilterBit(768),
    fPtMin(0.2),
    fPtMax(5.),
    fEtaMax(0.8),
    fNhitsMin(80),
    fChi2Max(4.0),
    fDeDxMin(10)
    // fDcaXyMax(3.),
    // fDcaZMax(3.)
{
 runNum       = -999;
 oldRunNum    = -999;
 runNumBin    = -999;
 for (int i=0; i<3; ++i) vtx[i] = -999;
 vzBin        = -999;
 cent         = -999;
 centBin      = -999;
 hEvtCount  = nullptr;
 hRunNumBin = nullptr;
 hCent      = nullptr;
 for (int i=0; i<3; ++i) hCentCorr[i] = nullptr;
 for (int i=0; i<2; ++i) hVxy[i]      = nullptr;
 for (int i=0; i<2; ++i) hVz[i]       = nullptr;
 for (int i = 0; i < 2; ++i) fHist2DMultCentQA[i] = nullptr;


  // Track-wise
 for (int i=0; i<2; ++i) hPt[i]    = nullptr;
 for (int i=0; i<2; ++i) hEta[i]   = nullptr;
 for (int i=0; i<8; ++i) hBeforePhi[i]   = nullptr;
 for (int i=0; i<8; ++i) hAfterPhi[i]    = nullptr;
 for (int i=0; i<2; ++i) hNhits[i] = nullptr;
 // for (int i=0; i<2; ++i) hDcaXy[i] = nullptr;
 // for (int i=0; i<2; ++i) hDcaZ[i]  = nullptr;
 hPDedx = nullptr;

  // pile up
 fSPDCutPU               = nullptr;
 fV0CutPU                = nullptr;
 fCenCutLowPU            = nullptr;
 fCenCutHighPU           = nullptr;
 fMultCutPU              = nullptr;

 //NUE
 IsDoNUE               = true;
 fListNUE              = nullptr;
 hNUEweightPlus        = nullptr;
 hNUEweightMinus       = nullptr;

 //NUA
 IsDoNUA              = true;
 fListNUA             = nullptr;

 hCorrectNUAPos       = nullptr;
 hCorrectNUANeg       = nullptr;

 for (int i=0; i<8; ++i) 
 {
  pos2Plane[i]         = nullptr;
  neg2Plane[i]         = nullptr;
  pos3Plane[i]         = nullptr;
  neg3Plane[i]         = nullptr;
 }
 for (int i=0; i<8; ++i) Res2Square[i]  = nullptr;
 for (int i=0; i<8; ++i) Res3Square[i]  = nullptr;
 for (int i=0; i<8; ++i) 
 {
  arcPos2Plane[i]       = nullptr;
  arcNeg2Plane[i]       = nullptr;
 }
 

  //CME
  //N(S)
  for (int i=0; i<8; ++i) hNDeltaSinSPsi2[i] = nullptr;
  //N(S_sf)
  for (int i=0; i<8; ++i) hNDeltaSinSPsi2_sf[i] = nullptr;
  //N(SVert)
  for (int i=0; i<8; ++i) hNDeltaSinSVertPsi2[i] = nullptr;
  //N(SVert_sf)
  for (int i=0; i<8; ++i) hNDeltaSinSVertPsi2_sf[i] = nullptr;
  //N(S)
  for (int i=0; i<8; ++i) hNDeltaSinSPsi3[i] = nullptr;
  //N(S_sf)
  for (int i=0; i<8; ++i) hNDeltaSinSPsi3_sf[i] = nullptr;
  //N(SVert)
  for (int i=0; i<8; ++i) hNDeltaSinSVertPsi3[i] = nullptr;
  //N(SVert_sf)
  for (int i=0; i<8; ++i) hNDeltaSinSVertPsi3_sf[i] = nullptr;

  //CMW
  //N(S)
  for (int i=0; i<8; ++i) hNDeltaCosSPsi2[i] = nullptr;
  //N(S_sf)
  for (int i=0; i<8; ++i) hNDeltaCosSPsi2_sf[i] = nullptr;
  //N(SVert)
  for (int i=0; i<8; ++i) hNDeltaCosSVertPsi2[i] = nullptr;
  //N(SVert_sf)
  for (int i=0; i<8; ++i) hNDeltaCosSVertPsi2_sf[i] = nullptr;
  //N(S)
  for (int i=0; i<8; ++i) hNDeltaCosSPsi3[i] = nullptr;
  //N(S_sf)
  for (int i=0; i<8; ++i) hNDeltaCosSPsi3_sf[i] = nullptr;
  //N(SVert)
  for (int i=0; i<8; ++i) hNDeltaCosSVertPsi3[i] = nullptr;
  //N(SVert_sf)
  for (int i=0; i<8; ++i) hNDeltaCosSVertPsi3_sf[i] = nullptr;

}
//_____________________________________________________________________________
AliAnalysisTaskMyTask::AliAnalysisTaskMyTask(const char* name) : AliAnalysisTaskSE(name),
    fAOD(nullptr), 
    fOutputList(nullptr), 
    mHarmonic(2.),
    fFilterBit(768),
    fPtMin(0.2),
    fPtMax(5.),
    fEtaMax(0.8),
    fNhitsMin(80),
    fChi2Max(4.),
    fDeDxMin(10)
   // fDcaXyMax(3.),
   // fDcaZMax(3.)
{
 runNum       = -999;
 oldRunNum    = -999;
 runNumBin    = -999;
 for (int i=0; i<3; ++i) vtx[i] = -999;
 vzBin        = -999;
 cent         = -999;
 centBin      = -999;
 hEvtCount  = nullptr;
 hRunNumBin = nullptr;
 hCent      = nullptr;
 for (int i=0; i<3; ++i) hCentCorr[i] = nullptr;
 for (int i=0; i<2; ++i) hVxy[i]      = nullptr;
 for (int i=0; i<2; ++i) hVz[i]       = nullptr;
 for (int i = 0; i < 2; ++i) fHist2DMultCentQA[i] = nullptr;



  // Track-wise
 for (int i=0; i<2; ++i) hPt[i]    = nullptr;
 for (int i=0; i<2; ++i) hEta[i]   = nullptr;
 for (int i=0; i<8; ++i) hBeforePhi[i]   = nullptr;
 for (int i=0; i<8; ++i) hAfterPhi[i]    = nullptr;
 for (int i=0; i<2; ++i) hNhits[i] = nullptr;
 // for (int i=0; i<2; ++i) hDcaXy[i] = nullptr;
 // for (int i=0; i<2; ++i) hDcaZ[i]  = nullptr;
 hPDedx = nullptr;

  // pile up
 fSPDCutPU               = nullptr;
 fV0CutPU                = nullptr;
 fCenCutLowPU            = nullptr;
 fCenCutHighPU           = nullptr;
 fMultCutPU              = nullptr;

  //NUE
 IsDoNUE               = true;
 fListNUE              = nullptr;
 hNUEweightPlus        = nullptr;
 hNUEweightMinus       = nullptr;

 //NUA
 IsDoNUA              = true;
 fListNUA             = nullptr;

 hCorrectNUAPos       = nullptr;
 hCorrectNUANeg       = nullptr;

 for (int i=0; i<8; ++i) 
 {
  pos2Plane[i]         = nullptr;
  neg2Plane[i]         = nullptr;
  pos3Plane[i]         = nullptr;
  neg3Plane[i]         = nullptr;
 }
 for (int i=0; i<8; ++i) Res2Square[i]  = nullptr;
 for (int i=0; i<8; ++i) Res3Square[i]  = nullptr;
 for (int i=0; i<8; ++i) 
 {
  arcPos2Plane[i]       = nullptr;
  arcNeg2Plane[i]       = nullptr;
 }

  //CME
  //N(S)
  for (int i=0; i<8; ++i) hNDeltaSinSPsi2[i] = nullptr;
  //N(S_sf)
  for (int i=0; i<8; ++i) hNDeltaSinSPsi2_sf[i] = nullptr;
  //N(SVert)
  for (int i=0; i<8; ++i) hNDeltaSinSVertPsi2[i] = nullptr;
  //N(SVert_sf)
  for (int i=0; i<8; ++i) hNDeltaSinSVertPsi2_sf[i] = nullptr;
  //N(S)
  for (int i=0; i<8; ++i) hNDeltaSinSPsi3[i] = nullptr;
  //N(S_sf)
  for (int i=0; i<8; ++i) hNDeltaSinSPsi3_sf[i] = nullptr;
  //N(SVert)
  for (int i=0; i<8; ++i) hNDeltaSinSVertPsi3[i] = nullptr;
  //N(SVert_sf)
  for (int i=0; i<8; ++i) hNDeltaSinSVertPsi3_sf[i] = nullptr;

  //CMW
  //N(S)
  for (int i=0; i<8; ++i) hNDeltaCosSPsi2[i] = nullptr;
  //N(S_sf)
  for (int i=0; i<8; ++i) hNDeltaCosSPsi2_sf[i] = nullptr;
  //N(SVert)
  for (int i=0; i<8; ++i) hNDeltaCosSVertPsi2[i] = nullptr;
  //N(SVert_sf)
  for (int i=0; i<8; ++i) hNDeltaCosSVertPsi2_sf[i] = nullptr;
  //N(S)
  for (int i=0; i<8; ++i) hNDeltaCosSPsi3[i] = nullptr;
  //N(S_sf)
  for (int i=0; i<8; ++i) hNDeltaCosSPsi3_sf[i] = nullptr;
  //N(SVert)
  for (int i=0; i<8; ++i) hNDeltaCosSVertPsi3[i] = nullptr;
  //N(SVert_sf)
  for (int i=0; i<8; ++i) hNDeltaCosSVertPsi3_sf[i] = nullptr;

 // constructor
 DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                     // this chain is created by the analysis manager, so no need to worry about it, 
                                     // it does its work automatically
 DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                     // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskMyTask::~AliAnalysisTaskMyTask()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskMyTask::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use 
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file
    //
    fOutputList = new TList();          // this is a list which will contain all of your histograms
    fOutputList->SetName(GetName());                                    // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)
    // event-wise
    hEvtCount = new TH1I("evtCount","",20, 1, 21);
    hEvtCount->GetXaxis()->SetBinLabel(1,"All");
    hEvtCount->GetXaxis()->SetBinLabel(2,"Info");
    hEvtCount->GetXaxis()->SetBinLabel(3,"Run number");
    hEvtCount->GetXaxis()->SetBinLabel(4,"Vertex");
    hEvtCount->GetXaxis()->SetBinLabel(5,"Cent");
    hEvtCount->GetXaxis()->SetBinLabel(6,"Pile up");
  
    hEvtCount->GetXaxis()->SetBinLabel(10,"Manager");
    hEvtCount->GetXaxis()->SetBinLabel(11,"Handler");
    hEvtCount->GetXaxis()->SetBinLabel(12,"AOD");
    hEvtCount->GetXaxis()->SetBinLabel(13,"PID");
    hEvtCount->GetXaxis()->SetBinLabel(14,"Utils");
    // hEvtCount->GetXaxis()->SetBinLabel(15,"MultSel");
  
   
    hEvtCount->GetXaxis()->SetBinLabel(17,"TPC plane");
    hEvtCount->GetXaxis()->SetBinLabel(18,"VZERO plane");
    hEvtCount->GetXaxis()->SetBinLabel(19,"ZDC plane");
    hEvtCount->GetXaxis()->SetBinLabel(20,"loops end");
    fOutputList->Add(hEvtCount);

     // 10h
   TString runNumList[12]={
     "296623","296622","296621","296619","296618","296616","296615","296594","296553","296552","296551","296550"
   };
   // 11h
   // TString runNum[39]={"170387","170040","170268","170228","170207","169838","170159","170204","170311","170084",
   //                     "169835","170088","170593","170203","170270","169846","170163","170388","170155","170083",
   //                     "170572","169837","169855","170306","170269","170089","170309","170091","170081","170230",
   //                     "170085","170315","170027","170193","170312","170313","170308","169858","169859"};
   hRunNumBin = new TH1I("runNumBin","",20,0,20);
   for (int i=0; i<12; ++i) {    
     hRunNumBin->GetXaxis()->SetBinLabel(i+1,runNumList[i].Data());
   }
   fOutputList->Add(hRunNumBin);
 
   hCent = new TH1D("centrality","",100,0,100);
   fOutputList->Add(hCent);
   hCentCorr[0] = new TH2D("centcorr0","",100,0,100,100,0,100);
   hCentCorr[1] = new TH2D("centcorr1","",100,0,100,100,0,100);
   hCentCorr[2] = new TH2D("centcorr2","",100,0,100,100,0,100);
   for (int i=0; i<3; ++i) fOutputList->Add(hCentCorr[i]);
 
   hVxy[0] = new TH2D("vxy0","",100,-0.5,0.5,100,-0.5,0.5);
   hVxy[1] = new TH2D("vxy1","",100,-0.5,0.5,100,-0.5,0.5);
   hVz[0]  = new TH1D("vz0","",200,-50,50);
   hVz[1]  = new TH1D("vz1","",200,-50,50);
   for (int i=0; i<2; ++i) fOutputList->Add(hVxy[i]);
   for (int i=0; i<2; ++i) fOutputList->Add(hVz[i]);

    fHist2DMultCentQA[0] = new TH2D("fHist2DMultCentQA_BfCut", ";centV0M;multFB32", 100, 0, 100, 20, 0, 5000);
    fHist2DMultCentQA[1] = new TH2D("fHist2DMultCentQA_AfCut", ";centV0M;multFB32", 100, 0, 100, 20, 0, 5000);
    fOutputList->Add(fHist2DMultCentQA[0]);
    fOutputList->Add(fHist2DMultCentQA[1]);

 
 
   // track-wise
   hPt[0] = new TH1D("hPtBeforeCut", "", 200, 0., 20.);
   hPt[1] = new TH1D("hPtAfterCut", "", 200, 0., 20.);
   for (int i=0; i<2; ++i) fOutputList->Add(hPt[i]);
   hEta[0] = new TH1D("hEtaBeforeCut", "", 200, -10., 10.);
   hEta[1] = new TH1D("hEtaAfterCut",  "", 200, -10., 10.);
   for (int i=0; i<2; ++i) fOutputList->Add(hEta[i]);
   for (int i=0; i<8; ++i){
   hBeforePhi[i] = new TH2D(Form("hPhiBeforeCut_cent%i",i), "", 400, 0, 2*TMath::Pi(), 100, -3., 3.);
   hAfterPhi[i]  = new TH2D(Form("hPhiAfterCut_cent%i",i), "", 400, 0, 2*TMath::Pi(), 100, -3., 3.);
   fOutputList->Add(hBeforePhi[i]);
   fOutputList->Add(hAfterPhi[i]);
   }
   hNhits[0] = new TH1D("hNhitsBeforeCut", "", 200, 0., 200.);
   hNhits[1] = new TH1D("hNhitsAfterCut",  "", 200, 0., 200.);
   for (int i=0; i<2; ++i) fOutputList->Add(hNhits[i]);
  // hDcaXy[0] = new TH1D("hDcaXyBeforeCut", "", 100, 0., 10.);
  // hDcaXy[1] = new TH1D("hDcaXyAfterCut",  "", 100, 0., 10.);
  // for (int i=0; i<2; ++i) fOutputList->Add(hDcaXy[i]);
  // hDcaZ[0] = new TH1D("hDcaZBeforeCut", "", 100, 0., 10.);
  // hDcaZ[1] = new TH1D("hDcaZAfterCut",  "", 100, 0., 10.);
  // for (int i=0; i<2; ++i) fOutputList->Add(hDcaZ[i]);
   hPDedx = new TH2D("hPDedx", "", 400, -10., 10., 400, 0, 1000);
   fOutputList->Add(hPDedx);

   //pile up
   fSPDCutPU = new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000);
   
   Double_t parV0[8] = {43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
   fV0CutPU  = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
   fV0CutPU->SetParameters(parV0);
   
   Double_t parV0CL0[6] = {0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};
   fCenCutLowPU  = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
   fCenCutLowPU->SetParameters(parV0CL0);
   fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
   fCenCutHighPU->SetParameters(parV0CL0);
   
   Double_t parFB32[8] = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};
   fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90);
   fMultCutPU->SetParameters(parFB32);

   ////////////////////////
   // NUE
   ////////////////////////
   // Load Calibration Files
   // The global read-in lists and hists are loaded here.
   // They do not need to be loaded run by run.
   if (IsDoNUE) {
     if (!fListNUE) {
       std::cout<<("NUE list not found")<<std::endl;
       return;
     }
     
       hNUEweightPlus  = (TH1D*)fListNUE->FindObject("trkEfficiencyChrgPos");
       hNUEweightMinus = (TH1D*)fListNUE->FindObject("trkEfficiencyChrgNeg");
     
   }

    ////////////////////////
    // NUA
    ////////////////////////
    if (IsDoNUA) {
      if (!fListNUA) {
        std::cout<<("NUA list not found")<<std::endl;
        return;
      }
        hCorrectNUAPos = new TH3F();
        hCorrectNUANeg = new TH3F();
    }

    // Plane
    for(int i=0; i<8; ++i){
    pos2Plane[i]         = new TH1D(Form("pos2Plane_cent%i",i),"", 180, -0.1, 2*TMath::Pi()+0.1); 
    fOutputList->Add(pos2Plane[i]);
    neg2Plane[i]         = new TH1D(Form("neg2Plane_cent%i",i),"", 180, -0.1, 2*TMath::Pi()+0.1);
    fOutputList->Add(neg2Plane[i]);
    pos3Plane[i]         = new TH1D(Form("pos3Plane_cent%i",i),"", 180, -0.1, 2*TMath::Pi()+0.1);
    fOutputList->Add(pos3Plane[i]);
    neg3Plane[i]         = new TH1D(Form("neg3Plane_cent%i",i),"", 180, -0.1, 2*TMath::Pi()+0.1);
    fOutputList->Add(neg3Plane[i]);
    }

   for(int i=0; i<8; ++i){
    Res2Square[i]    = new TProfile(Form("Res2Square_cent%i",i),"",1,0,1);
    fOutputList->Add(Res2Square[i]);
    Res3Square[i]       = new TProfile(Form("Res3Square_cent%i",i),"",1,0,1);
    fOutputList->Add(Res3Square[i]);
   }
   for(int i=0; i<8; ++i){
    arcPos2Plane[i]         = new TH1D(Form("arcPos2Plane_cent%i",i),"", 180, -1*TMath::Pi()-0.1, TMath::Pi()+0.1); 
    fOutputList->Add(arcPos2Plane[i]);
    arcNeg2Plane[i]         = new TH1D(Form("arcNeg2Plane_cent%i",i),"", 180, -1*TMath::Pi()-0.1, TMath::Pi()+0.1);
    fOutputList->Add(arcNeg2Plane[i]);
    }
   
   //CME
   //N(S)
   for (int i=0; i<8; ++i) hNDeltaSinSPsi2[i] = new TH1D(Form("hNDeltaSinSPsi2_cent%i",i),"",75,-1.5,1.5);
   //N(S_sf)
   for (int i=0; i<8; ++i) hNDeltaSinSPsi2_sf[i] = new TH1D(Form("hNDeltaSinSPsi2_sf_cent%i",i),"",75,-1.5,1.5);
   //N(SVert)
   for (int i=0; i<8; ++i) hNDeltaSinSVertPsi2[i] = new TH1D(Form("hNDeltaSinSVertPsi2_cent%i",i),"",75,-1.5,1.5);
   //N(SVert_sf)
   for (int i=0; i<8; ++i) hNDeltaSinSVertPsi2_sf[i] = new TH1D(Form("hNDeltaSinSVertPsi2_sf_cent%i",i),"",75,-1.5,1.5);
   
   //N(S)
   for (int i=0; i<8; ++i) hNDeltaSinSPsi3[i] = new TH1D(Form("hNDeltaSinSPsi3_cent%i",i),"",75,-1.5,1.5);
   //N(S_sf)
   for (int i=0; i<8; ++i) hNDeltaSinSPsi3_sf[i] = new TH1D(Form("hNDeltaSinSPsi3_sf_cent%i",i),"",75,-1.5,1.5);
   //N(SVert)
   for (int i=0; i<8; ++i) hNDeltaSinSVertPsi3[i] = new TH1D(Form("hNDeltaSinSVertPsi3_cent%i",i),"",75,-1.5,1.5);
   //N(SVert_sf)
   for (int i=0; i<8; ++i) hNDeltaSinSVertPsi3_sf[i] = new TH1D(Form("hNDeltaSinSvertPsi3_sf_cent%i",i),"",75,-1.5,1.5);
   
   for (int i=0; i<8; ++i) fOutputList->Add(hNDeltaSinSPsi2[i]);
   for (int i=0; i<8; ++i) fOutputList->Add(hNDeltaSinSPsi2_sf[i]);
   for (int i=0; i<8; ++i) fOutputList->Add(hNDeltaSinSVertPsi2[i]);
   for (int i=0; i<8; ++i) fOutputList->Add(hNDeltaSinSVertPsi2_sf[i]);
   for (int i=0; i<8; ++i) fOutputList->Add(hNDeltaSinSPsi3[i]);
   for (int i=0; i<8; ++i) fOutputList->Add(hNDeltaSinSPsi3_sf[i]);
   for (int i=0; i<8; ++i) fOutputList->Add(hNDeltaSinSVertPsi3[i]);
   for (int i=0; i<8; ++i) fOutputList->Add(hNDeltaSinSVertPsi3_sf[i]);
   
   //CMW
   //N(S)
   for (int i=0; i<8; ++i) hNDeltaCosSPsi2[i] = new TH1D(Form("hNDeltaCosSPsi2_cent%i",i),"",75,-1.5,1.5);
   //N(S_sf)
   for (int i=0; i<8; ++i) hNDeltaCosSPsi2_sf[i] = new TH1D(Form("hNDeltaCosSPsi2_sf_cent%i",i),"",75,-1.5,1.5);
   //N(SVert)
   for (int i=0; i<8; ++i) hNDeltaCosSVertPsi2[i] = new TH1D(Form("hNDeltaCosSVertPsi2_cent%i",i),"",75,-1.5,1.5);
   //N(SVert_sf)
   for (int i=0; i<8; ++i) hNDeltaCosSVertPsi2_sf[i] = new TH1D(Form("hNDeltaCosSVertPsi2_sf_cent%i",i),"",75,-1.5,1.5);
   
   //N(S)
   for (int i=0; i<8; ++i) hNDeltaCosSPsi3[i] = new TH1D(Form("hNDeltaCosSPsi3_cent%i",i),"",75,-1.5,1.5);
   //N(S_sf)
   for (int i=0; i<8; ++i) hNDeltaCosSPsi3_sf[i] = new TH1D(Form("hNDeltaCosSPsi3_sf_cent%i",i),"",75,-1.5,1.5);
   //N(SVert)
   for (int i=0; i<8; ++i) hNDeltaCosSVertPsi3[i] = new TH1D(Form("hNDeltaCosSVertPsi3_cent%i",i),"",75,-1.5,1.5);
   //N(SVert_sf)
   for (int i=0; i<8; ++i) hNDeltaCosSVertPsi3_sf[i] = new TH1D(Form("hNDeltaCosSvertPsi3_sf_cent%i",i),"",75,-1.5,1.5);
   
   for (int i=0; i<8; ++i) fOutputList->Add(hNDeltaCosSPsi2[i]);
   for (int i=0; i<8; ++i) fOutputList->Add(hNDeltaCosSPsi2_sf[i]);
   for (int i=0; i<8; ++i) fOutputList->Add(hNDeltaCosSVertPsi2[i]);
   for (int i=0; i<8; ++i) fOutputList->Add(hNDeltaCosSVertPsi2_sf[i]);
   for (int i=0; i<8; ++i) fOutputList->Add(hNDeltaCosSPsi3[i]);
   for (int i=0; i<8; ++i) fOutputList->Add(hNDeltaCosSPsi3_sf[i]);
   for (int i=0; i<8; ++i) fOutputList->Add(hNDeltaCosSVertPsi3[i]);
   for (int i=0; i<8; ++i) fOutputList->Add(hNDeltaCosSVertPsi3_sf[i]);
   
    
 
    PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the 
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliAnalysisTaskMyTask::UserExec(Option_t *)
{
 
    hEvtCount->Fill(1);

    AliAnalysisManager* manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager) {
    AliError(Form("%s: Could not get Analysis Manager", GetName()));
  } else hEvtCount->Fill(10);
  AliAODInputHandler* handler = (AliAODInputHandler*)manager->GetInputEventHandler();
  if (!handler) {
    AliError(Form("%s: Could not get Input Handler", GetName()));
  } else hEvtCount->Fill(11);
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) {
    AliError(Form("%s: Could not get AOD event", GetName()));
    return;
  } else hEvtCount->Fill(12);
  AliPIDResponse* fPID = handler->GetPIDResponse();
  if (!fPID) {
    AliError(Form("%s: Could not get PIDResponse", GetName()));
  } else hEvtCount->Fill(13);
  AliAnalysisUtils* fUtils = new AliAnalysisUtils();
  if (!fUtils) {
    AliError(Form("%s: Could not get AliAnalysisUtils", GetName()));
  } else hEvtCount->Fill(14);
  // AliMultSelection* fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
  // if (!fMultSel) {
  //   AliError(Form("%s: Could not get AliMultSelection", GetName()));
  // } else hEvtCount->Fill(15);
  if (!manager || !handler || !fAOD || !fUtils) return;
  hEvtCount->Fill(2);
    

  //------------------
  // event-wise
  //------------------

  // runNumber
  
  runNum = fAOD->GetRunNumber();
  if (runNum != oldRunNum) {  
  if (!LoadCalibHistForThisRun()) return;
  oldRunNum = runNum;
  }
  runNumBin = GetRunNumBin(runNum);
  if (runNumBin<0) return;
  hRunNumBin->Fill(runNumBin);
  hEvtCount->Fill(3);
   
    // vertex
  // AliAODVertex* fVtx = fAOD->GetPrimaryVertex();
  // vtx[0] = (double)fVtx->GetX();
  // vtx[1] = (double)fVtx->GetY();
  // vtx[2] = (double)fVtx->GetZ();
  // double vzSPD  = fAOD->GetPrimaryVertexSPD()->GetZ();
  // hVxy[0]->Fill(vtx[0], vtx[1]);
  // hVz[0]->Fill(vtx[2]);
  // if (fabs(vtx[0])<1e-6 || fabs(vtx[1])<1e-6 || fabs(vtx[2])<1e-6) return;
// 
  // if (fabs(vtx[2])>10) return;
  // if (fabs(vtx[2]-vzSPD)>0.5) return;
  // hVxy[1]->Fill(vtx[0], vtx[1]);
  // hVz[1]->Fill(vtx[2]);
  // for (int i = 0; i < 20; ++i) {
  //   if (vtx[2] > -10+i*1 && vtx[2] < -10+(i+1)*1) {vzBin = i; break;}
  // }
  // if (vzBin<0) return;



   const AliVVertex* vtTrc = fAOD->GetPrimaryVertex();
   const AliVVertex* vtSPD = fAOD->GetPrimaryVertexSPD();
   vtx[0] = (double)vtTrc->GetX();
   vtx[1] = (double)vtTrc->GetY();
   vtx[2] = (double)vtTrc->GetZ();
   hVxy[0]->Fill(vtx[0], vtx[1]);
   hVz[0]->Fill(vtx[2]);

   double covTrc[6],covSPD[6];
   vtTrc->GetCovarianceMatrix(covTrc);
   vtSPD->GetCovarianceMatrix(covSPD);
   double dz = vtTrc->GetZ()-vtSPD->GetZ();
   double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
   double errTrc = TMath::Sqrt(covTrc[5]);
   double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
   if (TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrc>20) return;

   hVxy[1]->Fill(vtx[0], vtx[1]);
   hVz[1]->Fill(vtx[2]);

   for (int i = 0; i < 20; ++i) {
      if (vtx[2] > -10+i*1 && vtx[2] < -10+(i+1)*1) {vzBin = i; break;}
   } 
   if (vzBin<-990) return;
   hEvtCount->Fill(4);
   

  // centrality   
  // run1   
  AliMultSelection* fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
  cent  = fMultSel->GetMultiplicityPercentile("V0M");
  double cent2  = fMultSel->GetMultiplicityPercentile("CL1");
 // double cent3  = fMultSel->GetMultiplicityPercentile("TRK"); // 0

 // cent = fAOD->GetCentrality()->GetCentralityPercentile("V0M");
 // double cent2 = fAOD->GetCentrality()->GetCentralityPercentile("CL1");
 // double cent3 = fAOD->GetCentrality()->GetCentralityPercentile("TRK"); //"CL0"
  hCentCorr[0]->Fill(cent,cent2);
 // hCentCorr[2]->Fill(cent,cent3);
  if (fabs(cent-cent2)>7.5) return;
  // if (fabs(cent-cent3)>5) return; // ANA-280
  if (cent<0 || cent>=80) return;
  hCentCorr[1]->Fill(cent,cent2);
  centBin = (int)cent/10;  //centbin 
  hCent->Fill(cent);
  hEvtCount->Fill(5);
  

   // pileup 18q
   if (!RejectEvtTFFit()) return;
   hEvtCount->Fill(6);
  
    

  //------------------
  //* loop trk
  //------------------

    vector<double> vecPhi;
    vector<int> vecCharge;
    vector<double> vecEta;
    vector<double> vecWeight;

    //TPC QxQy
    double sumCos2Pos = 0.;
    double sumCos2Neg = 0.;
    double sumSin2Pos = 0.;
    double sumSin2Neg = 0.;
    double sumCos3Pos = 0.;
    double sumCos3Neg = 0.;
    double sumSin3Pos = 0.;
    double sumSin3Neg = 0.;
   int nTrk = fAOD->GetNumberOfTracks();
   // int multESD = fAOD->GetNumberOfESDTrack();
   // if (multESD-3.38*nTrk>15000) return;  // 
   for (int iTrk = 0; iTrk < nTrk; ++iTrk) {
     AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrk));
      if (!track) {
         AliError(Form("%s: Could not get Track", GetName()));
         continue;
      }
     if (!track->TestFilterBit(fFilterBit)) continue;
     double pt     = track->Pt();
     double eta    = track->Eta();
     double phi    = track->Phi();
     int    charge = track->Charge();
     int    nhits  = track->GetTPCNcls();
     double dedx   = track->GetTPCsignal();
     double chi2   = track->Chi2perNDF();

   

     hBeforePhi[centBin]->Fill(phi,eta);
     hPt[0]->Fill(pt);
     hEta[0]->Fill(eta); 
     hNhits[0]->Fill(nhits);
 
     if (pt<fPtMin || pt>fPtMax) continue;
     if (fabs(eta)>fEtaMax) continue;
     if (fabs(nhits)<fNhitsMin) continue;
     if (chi2>fChi2Max) continue;
     if (dedx<fDeDxMin) continue;
     

     //------------------
     // NUE & NUA
     //------------------
      double weight=1;
      if (IsDoNUE) {
      double wEffi = GetNUECor(charge, pt);
      if (wEffi<0) continue;
      else weight *= wEffi;
    }
      if (IsDoNUA) {
      double wAcc = GetNUACor(charge, phi, eta, vtx[2]);
      if (wAcc<0) continue;
      else weight *= wAcc;
        
    }
     hAfterPhi[centBin]->Fill(phi, eta, weight);
     hPt[1]->Fill(pt,weight); 
     hEta[1]->Fill(eta);
     hNhits[1]->Fill(nhits);
     hPDedx->Fill(track->P()*charge, dedx);
 

     vecPhi.push_back(phi);
     vecCharge.push_back(charge);
     vecEta.push_back(eta);
     vecWeight.push_back(weight);

    //Qn
    if(eta>0.){
        sumCos2Pos += weight*cos(2*phi);
        sumSin2Pos += weight*sin(2*phi);
        sumCos3Pos += weight*cos(3*phi);
        sumSin3Pos += weight*sin(3*phi);
      }else if(eta<0.){
        sumCos2Neg += weight*cos(2*phi);
        sumSin2Neg += weight*sin(2*phi);
        sumCos3Neg += weight*cos(3*phi);
        sumSin3Neg += weight*sin(3*phi);
      }

   }


   //arcPlane
  double arcpsi2Pos = 0.5*TMath::ATan2(sumSin2Pos,sumCos2Pos);
  double arcpsi2Neg = 0.5*TMath::ATan2(sumSin2Neg,sumCos2Neg);
  arcPos2Plane[centBin]->Fill(arcpsi2Pos);
  arcNeg2Plane[centBin]->Fill(arcpsi2Neg);
   //TPC Plane
  TVector2 Q2TPCPos;
  Q2TPCPos.Set(sumCos2Pos,sumSin2Pos);
  double psi2Pos = Q2TPCPos.Phi()/2.;
  pos2Plane[centBin]->Fill(psi2Pos);
  TVector2 Q2TPCNeg;
  Q2TPCNeg.Set(sumCos2Neg,sumSin2Neg);
  double psi2Neg = Q2TPCNeg.Phi()/2.;
  neg2Plane[centBin]->Fill(psi2Neg);

  TVector2 Q3TPCPos;
  Q3TPCPos.Set(sumCos3Pos,sumSin3Pos);
  double psi3Pos = Q3TPCPos.Phi()/3.;
  pos3Plane[centBin]->Fill(psi3Pos);
  TVector2 Q3TPCNeg;
  Q3TPCNeg.Set(sumCos3Neg,sumSin3Neg);
  double psi3Neg = Q3TPCNeg.Phi()/3.;
  neg3Plane[centBin]->Fill(psi3Neg);
  hEvtCount->Fill(17);

   Res2Square[centBin]->Fill(0.5,cos(2*(psi2Pos-psi2Neg)));
   Res3Square[centBin]->Fill(0.5,cos(3*(psi3Pos-psi3Neg)));

  vector<int> vecCharge_sf;
  vecCharge_sf.assign(vecCharge.begin(),vecCharge.end());
  unsigned seed = std::chrono::system_clock::now ().time_since_epoch ().count ();  
  std::shuffle (vecCharge_sf.begin(), vecCharge_sf.end(), std::default_random_engine (seed)); 
  //random_shuffle(vecCharge_sf.begin(), vecCharge_sf.end());
 

  //Sin
  double sumSin2DeltaPhiPos = 0.;
  double sumSin2DeltaPhiPos_sf = 0.;
  double sumSin2DeltaPhiNeg = 0.;
  double sumSin2DeltaPhiNeg_sf = 0.;

  double sumSin3DeltaPhiPos = 0.;
  double sumSin3DeltaPhiPos_sf = 0.;
  double sumSin3DeltaPhiNeg = 0.;
  double sumSin3DeltaPhiNeg_sf = 0.;

  double sumSin2DeltaPhiVertPos = 0.;
  double sumSin2DeltaPhiVertPos_sf = 0.;
  double sumSin2DeltaPhiVertNeg = 0.;
  double sumSin2DeltaPhiVertNeg_sf = 0.;

  double sumSin3DeltaPhiVertPos = 0.;
  double sumSin3DeltaPhiVertPos_sf = 0.;
  double sumSin3DeltaPhiVertNeg = 0.;
  double sumSin3DeltaPhiVertNeg_sf = 0.;

  //Cos
  double sumCos2DeltaPhiPos = 0.;
  double sumCos2DeltaPhiPos_sf = 0.;
  double sumCos2DeltaPhiNeg = 0.;
  double sumCos2DeltaPhiNeg_sf = 0.;

  double sumCos3DeltaPhiPos = 0.;
  double sumCos3DeltaPhiPos_sf = 0.;
  double sumCos3DeltaPhiNeg = 0.;
  double sumCos3DeltaPhiNeg_sf = 0.;

  double sumCos2DeltaPhiVertPos = 0.;
  double sumCos2DeltaPhiVertPos_sf = 0.;
  double sumCos2DeltaPhiVertNeg = 0.;
  double sumCos2DeltaPhiVertNeg_sf = 0.;

  double sumCos3DeltaPhiVertPos = 0.;
  double sumCos3DeltaPhiVertPos_sf = 0.;
  double sumCos3DeltaPhiVertNeg = 0.;
  double sumCos3DeltaPhiVertNeg_sf = 0.;

  double nPos = 0;
  double nNeg = 0;
  double nPos_sf = 0;
  double nNeg_sf = 0;


  for (vector<double>::size_type iTrk = 0; iTrk < vecPhi.size(); iTrk++) {
    double phi = vecPhi[iTrk];
    int  charge  = vecCharge[iTrk];
    int  charge_sf = vecCharge_sf[iTrk];
    double eta = vecEta[iTrk];
    double weight = vecWeight[iTrk];

    // TVector2 Q2TPCPos_tmp;
    // TVector2 Q2TPCNeg_tmp;
    // TVector2 Q3TPCPos_tmp;
    // TVector2 Q3TPCNeg_tmp;
    // Q2TPC_tmp.Set(sumCos2Phi - TMath::Cos(2 * phi), sumSin2Phi - TMath::Sin(2 * phi));
    // Q3TPC_tmp.Set(sumCos3Phi - TMath::Cos(3 * phi), sumSin3Phi - TMath::Sin(3 * phi));
    // double psi2TPC_tmp = Q2TPC_tmp.Phi()/2;
    // double psi3TPC_tmp = Q3TPC_tmp.Phi()/3;
     if (eta<0.)
     {
           if(charge > 0) 
           {
             sumSin2DeltaPhiPos += weight*TMath::Sin(phi-psi2Pos);
             sumSin3DeltaPhiPos += weight*TMath::Sin(3/2*(phi-psi3Pos));
             sumSin2DeltaPhiVertPos += weight*TMath::Sin(phi-psi2Pos-TMath::Pi()/2);  
             sumSin3DeltaPhiVertPos += weight*TMath::Sin(3/2*(phi-psi3Pos-TMath::Pi()/3));
       
             sumCos2DeltaPhiPos += weight*TMath::Cos(2*(phi-psi2Pos));
             sumCos3DeltaPhiPos += weight*TMath::Cos(3*(phi-psi3Pos));
             sumCos2DeltaPhiVertPos += weight*TMath::Cos(2*(phi-psi2Pos-TMath::Pi()/2));
             sumCos3DeltaPhiVertPos += weight*TMath::Cos(3*(phi-psi3Pos-TMath::Pi()/3));
             nPos += weight;
           }
           if(charge_sf > 0) 
           {
             sumSin2DeltaPhiPos_sf += weight*TMath::Sin(phi-psi2Pos);
             sumSin3DeltaPhiPos_sf += weight*TMath::Sin(3/2*(phi-psi3Pos));
             sumSin2DeltaPhiVertPos_sf += weight*TMath::Sin(phi-psi2Pos-TMath::Pi()/2);
             sumSin3DeltaPhiVertPos_sf += weight*TMath::Sin(3/2*(phi-psi3Pos-TMath::Pi()/3));
       
             sumCos2DeltaPhiPos_sf += weight*TMath::Cos(2*(phi-psi2Pos));
             sumCos3DeltaPhiPos_sf += weight*TMath::Cos(3*(phi-psi3Pos));
             sumCos2DeltaPhiVertPos_sf += weight*TMath::Cos(2*(phi-psi2Pos-TMath::Pi()/2));
             sumCos3DeltaPhiVertPos_sf += weight*TMath::Cos(3*(phi-psi3Pos-TMath::Pi()/3));
             nPos_sf += weight;
           }
           if(charge < 0) 
           {
             sumSin2DeltaPhiNeg += weight*TMath::Sin(phi-psi2Pos);
             sumSin3DeltaPhiNeg += weight*TMath::Sin(3/2*(phi-psi3Pos));
             sumSin2DeltaPhiVertNeg += weight*TMath::Sin(phi-psi2Pos-TMath::Pi()/2);
             sumSin3DeltaPhiVertNeg += weight*TMath::Sin(3/2*(phi-psi3Pos-TMath::Pi()/3));
       
             sumCos2DeltaPhiNeg += weight*TMath::Cos(2*(phi-psi2Pos));
             sumCos3DeltaPhiNeg += weight*TMath::Cos(3*(phi-psi3Pos));
             sumCos2DeltaPhiVertNeg += weight*TMath::Cos(2*(phi-psi2Pos-TMath::Pi()/2));
             sumCos3DeltaPhiVertNeg += weight*TMath::Cos(3*(phi-psi3Pos-TMath::Pi()/3));
             nNeg += weight;
       
           }
           if(charge_sf < 0) 
           {
             sumSin2DeltaPhiNeg_sf += weight*TMath::Sin(phi-psi2Pos);
             sumSin3DeltaPhiNeg_sf += weight*TMath::Sin(3/2*(phi-psi3Pos));
             sumSin2DeltaPhiVertNeg_sf += weight*TMath::Sin(3/2*(phi-psi3Pos-TMath::Pi()/3));
       
             sumCos2DeltaPhiNeg_sf += weight*TMath::Cos(2*(phi-psi2Pos));
             sumCos3DeltaPhiNeg_sf += weight*TMath::Cos(3*(phi-psi3Pos));
             sumCos2DeltaPhiVertNeg_sf += weight*TMath::Cos(2*(phi-psi2Pos-TMath::Pi()/2));
             sumCos3DeltaPhiVertNeg_sf += weight*TMath::Cos(3*(phi-psi3Pos-TMath::Pi()/3));
             nNeg_sf += weight;
           }
     } else if(eta>0.){
      if(charge > 0) 
           {
             sumSin2DeltaPhiPos += weight*TMath::Sin(phi-psi2Neg);
             sumSin3DeltaPhiPos += weight*TMath::Sin(3/2*(phi-psi3Neg));
             sumSin2DeltaPhiVertPos += weight*TMath::Sin(phi-psi2Neg-TMath::Pi()/2);
             sumSin3DeltaPhiVertPos += weight*TMath::Sin(3/2*(phi-psi3Neg-TMath::Pi()/3));
       
             sumCos2DeltaPhiPos += weight*TMath::Cos(2*(phi-psi2Neg));
             sumCos3DeltaPhiPos += weight*TMath::Cos(3*(phi-psi3Neg));
             sumCos2DeltaPhiVertPos += weight*TMath::Cos(2*(phi-psi2Neg-TMath::Pi()/2));
             sumCos3DeltaPhiVertPos += weight*TMath::Cos(3*(phi-psi3Neg-TMath::Pi()/3));
             nPos += weight;
           }
           if(charge_sf > 0) 
           {
             sumSin2DeltaPhiPos_sf += weight*TMath::Sin(phi-psi2Neg);
             sumSin3DeltaPhiPos_sf += weight*TMath::Sin(3/2*(phi-psi3Neg));
             sumSin2DeltaPhiVertPos_sf += weight*TMath::Sin(phi-psi2Neg-TMath::Pi()/2);
             sumSin3DeltaPhiVertPos_sf += weight*TMath::Sin(3/2*(phi-psi3Neg-TMath::Pi()/3));
       
             sumCos2DeltaPhiPos_sf += weight*TMath::Cos(2*(phi-psi2Neg));
             sumCos3DeltaPhiPos_sf += weight*TMath::Cos(3*(phi-psi3Neg));
             sumCos2DeltaPhiVertPos_sf += weight*TMath::Cos(2*(phi-psi2Neg-TMath::Pi()/2));
             sumCos3DeltaPhiVertPos_sf += weight*TMath::Cos(3*(phi-psi3Neg-TMath::Pi()/3));
             nPos_sf += weight; 
           }
           if(charge < 0) 
           {
             sumSin2DeltaPhiNeg += weight*TMath::Sin(phi-psi2Neg);
             sumSin3DeltaPhiNeg += weight*TMath::Sin(3/2*(phi-psi3Neg));
             sumSin2DeltaPhiVertNeg += weight*TMath::Sin(phi-psi2Neg-TMath::Pi()/2);
             sumSin3DeltaPhiVertNeg += weight*TMath::Sin(3/2*(phi-psi3Neg-TMath::Pi()/3));
       
             sumCos2DeltaPhiNeg += weight*TMath::Cos(2*(phi-psi2Neg));
             sumCos3DeltaPhiNeg += weight*TMath::Cos(3*(phi-psi3Neg));
             sumCos2DeltaPhiVertNeg += weight*TMath::Cos(2*(phi-psi2Neg-TMath::Pi()/2));
             sumCos3DeltaPhiVertNeg += weight*TMath::Cos(3*(phi-psi3Neg-TMath::Pi()/3));
             nNeg += weight;
       
           }
           if(charge_sf < 0) 
           {
             sumSin2DeltaPhiNeg_sf += weight*TMath::Sin(phi-psi2Neg);
             sumSin3DeltaPhiNeg_sf += weight*TMath::Sin(3/2*(phi-psi3Neg));
             sumSin2DeltaPhiVertNeg_sf += weight*TMath::Sin(phi-psi2Neg-TMath::Pi()/2);
             sumSin3DeltaPhiVertNeg_sf += weight*TMath::Sin(3/2*(phi-psi3Neg-TMath::Pi()/3));
       
             sumCos2DeltaPhiNeg_sf += weight*TMath::Cos(2*(phi-psi2Neg));
             sumCos3DeltaPhiNeg_sf += weight*TMath::Cos(3*(phi-psi3Neg));
             sumCos2DeltaPhiVertNeg_sf += weight*TMath::Cos(2*(phi-psi2Neg-TMath::Pi()/2));
             sumCos3DeltaPhiVertNeg_sf += weight*TMath::Cos(3*(phi-psi3Neg-TMath::Pi()/3));
             nNeg_sf += weight;
           }

     }

 }
    if(abs(nPos * nNeg)<1e-6 || abs(nPos_sf*nNeg_sf)<1e-6) return;
  
    //Sin
    sumSin2DeltaPhiPos /= (double)nPos;
    sumSin2DeltaPhiPos_sf /= (double)nPos_sf;
    sumSin2DeltaPhiVertPos /= (double)nPos;
    sumSin2DeltaPhiVertPos_sf /= (double)nPos_sf;

    sumSin2DeltaPhiNeg /= (double)nNeg;
    sumSin2DeltaPhiNeg_sf /= (double)nNeg_sf;
    sumSin2DeltaPhiVertNeg /= (double)nNeg;
    sumSin2DeltaPhiVertNeg_sf /= (double)nNeg_sf;

    sumSin3DeltaPhiPos /= (double)nPos;
    sumSin3DeltaPhiPos_sf /= (double)nPos_sf;
    sumSin3DeltaPhiVertPos /= (double)nPos;
    sumSin3DeltaPhiVertPos_sf /= (double)nPos_sf;

    sumSin3DeltaPhiNeg /= (double)nNeg;
    sumSin3DeltaPhiNeg_sf /= (double)nNeg_sf;
    sumSin3DeltaPhiVertNeg /= (double)nNeg;
    sumSin3DeltaPhiVertNeg_sf /= (double)nNeg_sf;

    //Cos
    sumCos2DeltaPhiPos /= (double)nPos;
    sumCos2DeltaPhiPos_sf /= (double)nPos_sf;
    sumCos2DeltaPhiVertPos /= (double)nPos;
    sumCos2DeltaPhiVertPos_sf /= (double)nPos_sf;

    sumCos2DeltaPhiNeg /= (double)nNeg;
    sumCos2DeltaPhiNeg_sf /= (double)nNeg_sf;
    sumCos2DeltaPhiVertNeg /= (double)nNeg;
    sumCos2DeltaPhiVertNeg_sf /= (double)nNeg_sf;

    sumCos3DeltaPhiPos /= (double)nPos;
    sumCos3DeltaPhiPos_sf /= (double)nPos_sf;
    sumCos3DeltaPhiVertPos /= (double)nPos;
    sumCos3DeltaPhiVertPos_sf /= (double)nPos_sf;

    sumCos3DeltaPhiNeg /= (double)nNeg;
    sumCos3DeltaPhiNeg_sf /= (double)nNeg_sf;
    sumCos3DeltaPhiVertNeg /= (double)nNeg;
    sumCos3DeltaPhiVertNeg_sf /= (double)nNeg_sf;

    //Sin
    double deltaSinS2 = sumSin2DeltaPhiPos - sumSin2DeltaPhiNeg;
    double deltaSinS2_sf = sumSin2DeltaPhiPos_sf - sumSin2DeltaPhiNeg_sf;
    double deltaSinS2Vert = sumSin2DeltaPhiVertPos - sumSin2DeltaPhiVertNeg;
    double deltaSinS2Vert_sf = sumSin2DeltaPhiVertPos_sf - sumSin2DeltaPhiVertNeg_sf;

    double deltaSinS3 = sumSin3DeltaPhiPos - sumSin3DeltaPhiNeg;
    double deltaSinS3_sf = sumSin3DeltaPhiPos_sf - sumSin3DeltaPhiNeg_sf;
    double deltaSinS3Vert = sumSin3DeltaPhiVertPos - sumSin3DeltaPhiVertNeg;
    double deltaSinS3Vert_sf = sumSin3DeltaPhiVertPos_sf - sumSin3DeltaPhiVertNeg_sf;

    //Cos
    double deltaCosS2 = sumCos2DeltaPhiPos - sumCos2DeltaPhiNeg;
    double deltaCosS2_sf = sumCos2DeltaPhiPos_sf - sumCos2DeltaPhiNeg_sf;
    double deltaCosS2Vert = sumCos2DeltaPhiVertPos - sumCos2DeltaPhiVertNeg;
    double deltaCosS2Vert_sf = sumCos2DeltaPhiVertPos_sf - sumCos2DeltaPhiVertNeg_sf;

    double deltaCosS3 = sumCos3DeltaPhiPos - sumCos3DeltaPhiNeg;
    double deltaCosS3_sf = sumCos3DeltaPhiPos_sf - sumCos3DeltaPhiNeg_sf;
    double deltaCosS3Vert = sumCos3DeltaPhiVertPos - sumCos3DeltaPhiVertNeg;
    double deltaCosS3Vert_sf = sumCos3DeltaPhiVertPos_sf - sumCos3DeltaPhiVertNeg_sf;

    //N(deltaSinS)
    hNDeltaSinSPsi2[centBin]->Fill(deltaSinS2);
    //N(deltaSinS_sf)
    hNDeltaSinSPsi2_sf[centBin]->Fill(deltaSinS2_sf);
    //N(deltaSinSVert)
    hNDeltaSinSVertPsi2[centBin]->Fill(deltaSinS2Vert);
    //N(deltaSinSVert_sf)
    hNDeltaSinSVertPsi2_sf[centBin]->Fill(deltaSinS2Vert_sf);
    //N(deltaSinS)
    hNDeltaSinSPsi3[centBin]->Fill(deltaSinS3);
    //N(deltaSinS_sf)
    hNDeltaSinSPsi3_sf[centBin]->Fill(deltaSinS3_sf);
    //N(deltaSinSVert)
    hNDeltaSinSVertPsi3[centBin]->Fill(deltaSinS3Vert);
    //N(deltaSinSVert_sf)
    hNDeltaSinSVertPsi3_sf[centBin]->Fill(deltaSinS3Vert_sf);

    //N(deltaSinS)
    hNDeltaCosSPsi2[centBin]->Fill(deltaCosS2);
    //N(deltaSinS_sf)
    hNDeltaCosSPsi2_sf[centBin]->Fill(deltaCosS2_sf);
    //N(deltaSinSVert)
    hNDeltaCosSVertPsi2[centBin]->Fill(deltaCosS2Vert);
    //N(deltaSinSVert_sf)
    hNDeltaCosSVertPsi2_sf[centBin]->Fill(deltaCosS2Vert_sf);
    //N(deltaSinS)
    hNDeltaCosSPsi3[centBin]->Fill(deltaCosS3);
    //N(deltaSinS_sf)
    hNDeltaCosSPsi3_sf[centBin]->Fill(deltaCosS3_sf);
    //N(deltaSinSVert)
    hNDeltaCosSVertPsi3[centBin]->Fill(deltaCosS3Vert);
    //N(deltaSinSVert_sf)
    hNDeltaCosSVertPsi3_sf[centBin]->Fill(deltaCosS3Vert_sf);

    hEvtCount->Fill(20);
    PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file
}




//---------------------------------------------------
int AliAnalysisTaskMyTask::GetRunNumBin(int runNum)
{
  int runNumBin=-1;
  // 10h
  int runNumList[12]={
    296623,296622,296621,296619,296618,296616,296615,296594,296553,296552,296551,296550
  };
  // 11h
  // int runNumList[39]={170387,170040,170268,170228,170207,169838,170159,170204,170311,170084,
  //                     169835,170088,170593,170203,170270,169846,170163,170388,170155,170083,
  //                     170572,169837,169855,170306,170269,170089,170309,170091,170081,170230,
  //                     170085,170315,170027,170193,170312,170313,170308,169858,169859};
  for (int i = 0; i < 12; ++i) {
    if (runNum==runNumList[i]) {runNumBin=i; break;}
    else continue;
  }
  return runNumBin;
}


//---------------------------------------------------

bool AliAnalysisTaskMyTask::RejectEvtTFFit()
{
  Float_t centV0M = -999;
  Float_t centCL1 = -999;
  Float_t centCL0 = -999;

  AliMultSelection* fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
  if (!fMultSelection) {
    printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
    exit(1);
  }
  centV0M = (Float_t) cent;
  centCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
  centCL0 = fMultSelection->GetMultiplicityPercentile("CL0");
 

  Int_t nITSClsLy0 = fAOD->GetNumberOfITSClusters(0);
  Int_t nITSClsLy1 = fAOD->GetNumberOfITSClusters(1);
  Int_t nITSCls = nITSClsLy0 + nITSClsLy1;

  AliAODTracklets* aodTrkl = (AliAODTracklets*)fAOD->GetTracklets();
  Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();

  const Int_t nTracks = fAOD->GetNumberOfTracks();
  Int_t multTrk = 0;

  for (Int_t it = 0; it < nTracks; it++) {
      AliAODTrack* aodTrk = (AliAODTrack*)fAOD->GetTrack(it);
      if (!aodTrk) {
          delete aodTrk;
          continue;
      }
      if (aodTrk->TestFilterBit(32)) {
        if ((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2))
        multTrk++;
      }
  }

  fHist2DMultCentQA[0]->Fill(centV0M, multTrk); //  Mult(FB32) Vs Cent(V0M) 
  AliAODVZERO* aodV0 = fAOD->GetVZEROData();
  Float_t  multV0a = aodV0->GetMTotV0A();
  Float_t  multV0c = aodV0->GetMTotV0C();
  Float_t  multV0Tot = multV0a + multV0c;
  UShort_t multV0aOn = aodV0->GetTriggerChargeA(); 
  UShort_t multV0cOn = aodV0->GetTriggerChargeC();
  UShort_t multV0On = multV0aOn + multV0cOn;
  
  
  // //
  // Int_t tpcClsTot = fAOD->GetNumberOfTPCClusters();
  // Float_t nclsDif = Float_t(tpcClsTot) - (53182.6 + 113.326*multV0Tot - 0.000831275*multV0Tot*multV0Tot);

  // pile-up cuts
  if (centCL0 < fCenCutLowPU->Eval(centV0M)) return false;
  if (centCL0 > fCenCutHighPU->Eval(centV0M)) return false;
  if (Float_t(nITSCls) > fSPDCutPU->Eval(nITSTrkls)) return false;
  if (multV0On < fV0CutPU->Eval(multV0Tot)) return false;
  if (Float_t(multTrk) < fMultCutPU->Eval(centV0M)) return false;
  if (((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) return false;
  if (fAOD->IsIncompleteDAQ()) return false;


  fHist2DMultCentQA[1]->Fill(centV0M, multTrk); //  Mult(FB32) Vs Cent(V0M)
  return true;
}


//_____________________________________________________________________________
void AliAnalysisTaskMyTask::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________


bool AliAnalysisTaskMyTask::LoadCalibHistForThisRun()
{
  //18q/r NUA
    if (IsDoNUA) {
      hCorrectNUAPos -> Reset();
      hCorrectNUANeg -> Reset();
      hCorrectNUAPos = (TH3F*) fListNUA->FindObject(Form("fHist_NUA_VzPhiEta_kPID%dPos_Run%d",0,runNum));
      hCorrectNUANeg = (TH3F*) fListNUA->FindObject(Form("fHist_NUA_VzPhiEta_kPID%dNeg_Run%d",0,runNum));
      if (!hCorrectNUAPos) return false;
      if (!hCorrectNUANeg) return false;
    }
    return true;
}


//---------------------------------------------------


double AliAnalysisTaskMyTask::GetNUACor(int charge, double phi, double eta, double vz)
{
  double weightNUA = 1;
  
  if (vzBin<0 || centBin<0 || runNum<0) return -1;
 
      // Rihan and Protty 's NUA Results
    if (charge>0) {
      if (!hCorrectNUAPos) return -1;
      int iBinNUA = hCorrectNUAPos->FindBin(vz,phi,eta);
      if (hCorrectNUAPos->GetBinContent(iBinNUA)>0) weightNUA = (double)hCorrectNUAPos->GetBinContent(iBinNUA);
      return  weightNUA;
    } else if (charge<0) {
      if (!hCorrectNUANeg) return -1;
      int iBinNUA = hCorrectNUANeg->FindBin(vz,phi,eta);
      if (hCorrectNUANeg->GetBinContent(iBinNUA)>0) weightNUA = (double)hCorrectNUANeg->GetBinContent(iBinNUA);
      return weightNUA;
    }
    // In Rihan and Protty 's NUA results, the phi distribution is independent on centrality and particle charge
  
  return weightNUA;
}


//---------------------------------------------------

double AliAnalysisTaskMyTask::GetNUECor(int charge, double pt)
{
  double weightNUE = 1;

 
    if (charge>0) {
   //   hNUEweightPlus = (TH1D*)fListNUE->FindObject("trkEfficiencyChrgPos");
      if (!hNUEweightPlus) return -1;
      int ptBin = hNUEweightPlus->GetXaxis()->FindBin(pt);
      if (hNUEweightPlus->GetBinContent(ptBin)>0) {
        weightNUE = 1./ hNUEweightPlus->GetBinContent(ptBin);
      }
      else return -1;
    }
    if (charge<0) {
   //   hNUEweightMinus = (TH1D*)fListNUE->FindObject("trkEfficiencyChrgNeg");
      if (!hNUEweightMinus) return -1;
      int ptBin = hNUEweightMinus->GetXaxis()->FindBin(pt);
      if (hNUEweightMinus->GetBinContent(ptBin)>0) {
        weightNUE = 1./ hNUEweightMinus->GetBinContent(ptBin);
      }
      else return -1;
    }

  return weightNUE;

}














