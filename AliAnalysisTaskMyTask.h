/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskMyTask_H
#define AliAnalysisTaskMyTask_H

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskMyTask : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTaskMyTask();
                                AliAnalysisTaskMyTask(const char *name);
        virtual                 ~AliAnalysisTaskMyTask();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
         
        void SetListForNUE(TList *flist)         {this->fListNUE        = (TList*)flist->Clone();}
        void SetListForNUA(TList *flist)         {this->fListNUA        = (TList*)flist->Clone();}

        void SetNUEOn(bool doNUE)                {this->IsDoNUE      =       doNUE;}
        void SetNUAOn(bool doNUA)                {this->IsDoNUA      =       doNUA;}

    private:
     
        double                 GetNUACor(int charge, double phi, double eta, double vz);
        double                 GetNUECor(int charge, double pt);
        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list
        

        double mHarmonic;
        unsigned int fFilterBit;
        double fPtMin;
        double fPtMax;
        double fEtaMax;
        double fNhitsMin;
        double fChi2Max;
        double fDeDxMin;
        //double fDcaXyMax;
        //double fDcaZMax;
 
        int   GetRunNumBin(int runNum);
        bool  LoadCalibHistForThisRun();
        bool  RejectEvtTFFit();
      
        int runNum;
        int oldRunNum;
        int runNumBin;
        double vtx[3];
        int vzBin;
        double cent;
        int centBin;
    

          //Event-wise
        TH1I *hEvtCount;
        TH1I *hRunNumBin;
        TH1D *hCent;
        TH2D *hCentCorr[3];
        TH2D *hVxy[2];
        TH1D *hVz[2];
        TH2D*     fHist2DMultCentQA[2];

        

          // Track-wise
        TH1D *hPt[2];
        TH1D *hEta[2];
        TH2D *hBeforePhi[8];
        TH2D *hAfterPhi[8];
        TH1D *hNhits[2];
       // TH1D *hDcaXy[2];
       // TH1D *hDcaZ[2];
        TH2D *hPDedx;

        

        // pile up function
        TF1*                       fSPDCutPU;
        TF1*                        fV0CutPU;
        TF1*                    fCenCutLowPU;
        TF1*                   fCenCutHighPU;
        TF1*                      fMultCutPU;


        ////////////////////////
        // NUE
        ////////////////////////
        //10h/15o
        bool                        IsDoNUE;
        TList*                      fListNUE; // read list for NUE
        TH1D*                 hNUEweightPlus;
        TH1D*                hNUEweightMinus;

         ////////////////////////
         // NUA
         ////////////////////////
         bool            IsDoNUA;
         TList*          fListNUA; // read lists for NUA
     
 
         TH3F*    hCorrectNUAPos; // Protty
         TH3F*    hCorrectNUANeg; // Protty

        //plane
        TH1D*       pos2Plane[8];
        TH1D*       neg2Plane[8];
        TH1D*       pos3Plane[8];
        TH1D*       neg3Plane[8];
        TProfile*   Res2Square[8];
        TProfile*   Res3Square[8];
        TH1D*     arcPos2Plane[8];
        TH1D*     arcNeg2Plane[8];

        //CME
         //N(S) (150,-1.5,1.5)
         TH1D* hNDeltaSinSPsi2[8];
         //N(S_sf)
         TH1D* hNDeltaSinSPsi2_sf[8];
         //N(SVert)
         TH1D* hNDeltaSinSVertPsi2[8];
         //N(SVert_sf)
         TH1D* hNDeltaSinSVertPsi2_sf[8];
         //N(S) (150,-1.5,1.5)
         TH1D* hNDeltaSinSPsi3[8];
         //N(S_sf)
         TH1D* hNDeltaSinSPsi3_sf[8];
         //N(SVert)
         TH1D* hNDeltaSinSVertPsi3[8];
         //N(SVert_sf)
         TH1D* hNDeltaSinSVertPsi3_sf[8];
       
         //CMW
         //N(S) (150,-1.5,1.5)
         TH1D* hNDeltaCosSPsi2[8];
         //N(S_sf)
         TH1D* hNDeltaCosSPsi2_sf[8];
         //N(SVert)
         TH1D* hNDeltaCosSVertPsi2[8];
         //N(SVert_sf)
         TH1D* hNDeltaCosSVertPsi2_sf[8];
         //N(S) (150,-1.5,1.5)
         TH1D* hNDeltaCosSPsi3[8];
         //N(S_sf)
         TH1D* hNDeltaCosSPsi3_sf[8];
         //N(SVert)
         TH1D* hNDeltaCosSVertPsi3[8];
         //N(SVert_sf)
         TH1D* hNDeltaCosSVertPsi3_sf[8];
        


        AliAnalysisTaskMyTask(const AliAnalysisTaskMyTask&); // not implemented
        AliAnalysisTaskMyTask& operator=(const AliAnalysisTaskMyTask&); // not implemented

        ClassDef(AliAnalysisTaskMyTask, 1);
};

#endif
