#include "parameterisation.h"
#include "TChain.h"
#include "TSystem.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TF1.h"
#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include "TImage.h"
#include <iostream>
#include <iomanip>
#include <fstream>

#include <cstdio>

#include <vector>
#include <string>

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//
// TABLE OF CONTENTS:
// 1. PRINTPSFILE - DNE
// 2. FIF - SQUARE ROOT PARAMETERISATION
// 3. MAKELOG BINNING - DNE
// 4. INITIALISE  (1 and 2)
// 5. PROCESS     (1 and 2)
// 6. WIDTH CALCULATION
// 7. SQUARE ROOT FIT
// 8. PULL
// 9. TERMINATE 
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double round(long double number, int precision) {
  int decimals = std::pow(10, precision);
  return (std::round(number * decimals)) / decimals;
}

void printpsfile(std::string, TFile* f);


/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

Double_t fitf(Double_t *x,Double_t *par) {
  Double_t fitval = TMath::Sqrt(par[0] + par[1]*x[0]*x[0]);

  return fitval;
}

double gaussFunc(double *x, double *par) {
   double s;
   s = par[0]*exp(-(x[0]-par[1])*(x[0]-par[1])/par[3]/par[3]/2   );
   return s;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double* MakeLogBinning(int n_bins, double xmin, double xmax){
   double *binning = new double[n_bins+1];
   double delta_x_log = (log10(xmax)-log10(xmin))/n_bins;
   binning[0]=xmin;
   for(int i=1;i<=n_bins;++i){
      binning[i] = pow(10,log10(binning[i-1])+ delta_x_log);

   }
   return binning;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void Init() {
  // define the boundaries for
  int nParams= 5;
  nbins      = 50;
  nbinspull  = 200; // nbinspull 80 default
  maxntracks = 500;
  d0min      = -2.         ; d0max     = 2.         ;
  z0min      = -100.       ; z0max     = 100.       ;
  phimin     = -TMath::Pi(); phimax    = TMath::Pi();
  etamin     = -2.4        ; etamax    = 2.4        ;
  invpt_min  = -0.5e-3     ; invpt_max = 0.5e-3     ;
  pTmin      = 1000           ; ptmax     = 100000  ;
  abscurvmax = .5e-3;
//  Dd0 = 0.4; Drelpt = .00005; Dcurv = .05; Dphi = .015; Dz0 = 2.8; Deta = .020;
 Dd0 = 0.8;    Drelpt = .00010; Dcurv = .15; Dphi = .05; Dz0 = 4.0; Deta = .040;

  trackParam_range[0] = Dd0;
  trackParam_range[1] = Dz0;
  trackParam_range[2] = Deta;
  trackParam_range[3] = Dphi;
  trackParam_range[4] = Drelpt;
  
  ninvptbins = 26; //26
  stepsize = (invpt_max - invpt_min)/ninvptbins;
  invptval = invpt_min;
 
  while(invptval <= invpt_max+0.5*stepsize){
    invptbinsvec.push_back(invptval);
    invptval = invptval +stepsize;

  }
  double etabinsarray[] = {0.0,0.5,1.0,1.5,2.0,2.5};
//  double etabinsarray[] = {0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5};
  double d0binsarray[]  = {-2.0 ,-1.0 ,0 ,1.0 ,2.0 };
  double z0binsarray[]  = {-110.0 ,-50.0  ,0.0 ,50.0  ,110.0 };
  double ptbinsarray[] = {0,10,20,30,40,50};

  for (auto &bin :etabinsarray){etabins.push_back(bin);}
  for (auto &bin :d0binsarray ){d0bins.push_back(bin); }
  for (auto &bin :z0binsarray ){z0bins.push_back(bin); }
  for (auto &bin :ptbinsarray ){ptbinsvec.push_back(bin); }
  
  for(int itp = 0; itp < nParams;itp++){
    string trackParam = trackParam_list[itp];
    TString TtrackParam(trackParam);
    TString TtrackParam_title(";" + trackParam + "(reco-truth)/#sigma;N Tracks");

    double range      = trackParam_range[itp];
    for (int iibl = 0; iibl < 2;iibl++){
      string iblname = iblnames[iibl];
      for( int ieta = 0; ieta < etabins.size();ieta++){
	       for( int iipt = 0; iipt < invptbinsvec.size();iipt++){
	          string hist_res_name_string = "hist_res" + trackParam + "_" + iblname + "_eta" + to_string(ieta) + "_ipt" + to_string(iipt);
	          TString hist_res_name(hist_res_name_string);
	          TString hist_res_title(";#Delta" + trackParam + "(mm);N Tracks");
	       	  if(iibl == 0){
		          nbins = 500;
		          hist_res[itp][iibl][ieta][iipt] = new TH1F( hist_res_name,hist_res_title,nbins, -range,range);
		        }else{
		          nbins = 500;
		          hist_res[itp][iibl][ieta][iipt] = new TH1F( hist_res_name,hist_res_title,nbins, -range,range);
		        }
	       }
      }
    }
    pull_res[itp] = new TH1F( TtrackParam,TtrackParam_title,nbinspull,-5,5);

  }

  ientry2=0;
}

void Init2() {
  // define the boundaries for

  for(int itp = 0; itp < nParams;itp++){
    string trackParam = trackParam_list[itp];
    TString TtrackParam(trackParam);
    TString TtrackParam_title(";" + trackParam + "(reco-truth)/#sigma;N Tracks");

    double range      = trackParam_range[itp];
    for (int iibl = 0; iibl < 2;iibl++){
      string iblname = iblnames[iibl];
      for( int ieta = 0; ieta < etabins.size();ieta++){
	       for( int iipt = 0; iipt < invptbinsvec.size();iipt++){
	          string hist_res_name_string = "hist_rebinned_res" + trackParam + "_" + iblname + "_eta" + to_string(ieta) + "_ipt" + to_string(iipt);
	          TString hist_res_name(hist_res_name_string);
	          TString hist_res_title(";#Delta" + trackParam + "(mm);N Tracks");
	       	  if(iibl == 0){
		          nbins = 25;
              hist_std[itp][iibl][ieta][iipt] = hist_res[itp][iibl][ieta][iipt]->GetStdDev();
              delete hist_res[itp][iibl][ieta][iipt];
		          hist_res[itp][iibl][ieta][iipt] = new TH1F( hist_res_name,hist_res_title,nbins, -4.5*hist_std[itp][iibl][ieta][iipt],4.5*hist_std[itp][iibl][ieta][iipt]);
		        }else{
		          nbins = 25;
              hist_std[itp][iibl][ieta][iipt] = hist_res[itp][iibl][ieta][iipt]->GetStdDev();
              delete hist_res[itp][iibl][ieta][iipt];
		          hist_res[itp][iibl][ieta][iipt] = new TH1F( hist_res_name,hist_res_title,nbins, -4.5*hist_std[itp][iibl][ieta][iipt],4.5*hist_std[itp][iibl][ieta][iipt]);
		        }
	       }
      }
    }
  }

  ientry2=0;
}



/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void Process(Long64_t ientry) {
  t_ftkdata->GetEntry(ientry);

  // collect information from the run-evt from the truth
  Int_t evtnumber_ftk = tracks->eventNumber();
  Int_t runnumber_ftk = tracks->runNumber();

  // first try to start from last number
  for (; ientry2 < t_truth->GetEntries(); ientry2++) {
    t_evtinfo->GetEntry(ientry2);
    if (EventNumber==evtnumber_ftk && RunNumber==runnumber_ftk) break;
  }

  // only enter this loop if we didn't find it above
  if (ientry2 == t_truth->GetEntries()) {
    for (ientry2 = 0; ientry2 < t_truth->GetEntries(); ientry2++) {
      t_evtinfo->GetEntry(ientry2);
      if (EventNumber==evtnumber_ftk && RunNumber==runnumber_ftk) break;
    }
  }

  if (ientry2 == t_truth->GetEntries()) {
    cerr << Form("Mismatch between entries: (%d,%d)",runnumber_ftk,evtnumber_ftk) << endl;
    ientry2=0;
    return;
  }

  t_truth->GetEntry(ientry2);
  if (ientry%10000==0) { Info("Process","Event %lld, (Run,Evt) = (%d,%d)",ientry, runnumber_ftk, evtnumber_ftk);}

  Int_t ntracks = tracks->getNTracks();
  if (Use1stStage == 1)
  ntracks = tracks->getNTracksI();

  FTKBarcodeMM ftkmatchinfo;
  const FTKTrack *curtrack;
  if (Use1stStage != -1) {

    for (Int_t itrk=0;itrk!=ntracks;++itrk) { // loop over the FTK tracks
      curtrack = (Use1stStage ? tracks->getTrackI(itrk) : tracks->getTrack(itrk));
      for ( Int_t ii=0; ii < 12; ++ii){
      	FTKHit hit = curtrack->getFTKHit(ii);
      }
      if (curtrack->getBitmask()&(1<<0) || !(curtrack->getBitmask()&(1<<0)) ) {
	       FTKHit hit = curtrack->getFTKHit(0);
	       bool isIBL = hit.getPlane() == 0;
         if (curtrack->getEventIndex()==0) {
	           // Information on FTK tracks relative to the hard-scattering events
	           // are collected in a vector and later used to build matching maps
	         if (curtrack->getBarcodeFrac()>.5) {
	           ftkmatchinfo.insert(pair<MatchInfo, const FTKTrack*>(MatchInfo(curtrack->getBarcode(),curtrack->getEventIndex()),curtrack));
	         }
	       }
	
      	// apply some basic quality selection
      	double d0 = curtrack->getIP();
        double z0 = curtrack->getZ0();
        double phi = curtrack->getPhi();
        double curv = curtrack->getHalfInvPt();
        double eta = curtrack->getEta();
	      if (d0<d0min || d0>d0max) continue;
	      if (z0<z0min || z0>z0max) continue;
      	if (phi<phimin || phi > phimax) continue;
      	if (curv<-abscurvmax || curv>abscurvmax) continue;
      	if (eta<etamin || eta>etamax) continue;
	
      } // end loop over the FTK tracks
    } // end if not using roads  
  }

  Int_t ntruth = truthTracks->size();
  Int_t ntruth_good(0);

  vector<FTKTruthTrack>::const_iterator itr = truthTracks->begin();
  vector<FTKTruthTrack>::const_iterator itrE = truthTracks->end();
  for (;itr!=itrE;++itr) { // loop over the truth tracks
    const FTKTruthTrack &curtruth = (*itr);

    int barcode = curtruth.getBarcode();
    double px = curtruth.getPX();
    double py = curtruth.getPY();
    double pt = TMath::Sqrt(px*px+py*py);
    double invpt = 1./(2.0*pt);
    double d0 = curtruth.getD0();
    double z0 = curtruth.getZ();
    double curv = curtruth.getQ()*invpt;
    double phi = curtruth.getPhi();
    double eta = curtruth.getEta();
    double qOverp = curv;
    int pdgcode = curtruth.getPDGCode();

    if (barcode>100000 || barcode==0) continue;
    if ( pt < ptmincut ) continue;
    if (d0<d0min || d0>d0max) continue;
    if (z0<z0min || z0>z0max) continue;
    if (curv<-abscurvmax || curv>abscurvmax) continue;
    if (phi<phimin || phi>phimax) continue;
    if (eta<etamin || eta>etamax) continue;
    if (curtruth.getEventIndex()!=0 && curtruth.getQ()==0) continue;
    //EDIT TRACK CRITERIA

    // match the barcode and event index values
    MatchInfo reftruth(barcode,curtruth.getEventIndex());
    pair<FTKBarcodeMM::const_iterator,FTKBarcodeMM::const_iterator> mrange = ftkmatchinfo.equal_range(reftruth);
    if (mrange.first != mrange.second) {
      const FTKTrack *bestftk(0x0);
      for(FTKBarcodeMM::const_iterator ftkI = mrange.first;ftkI!=mrange.second;++ftkI) {
        if (!bestftk){
          bestftk = (*ftkI).second;
        } else if (bestftk->getBarcodeFrac()<(*ftkI).second->getBarcodeFrac()) {
          bestftk = (*ftkI).second;
        }
      }
      
      if (bestftk) {
	       FTKHit hit = bestftk->getFTKHit(0);
	       bool isIBL = hit.getPlane() == 0;
	       for(int ieta = 1; ieta < etabins.size(); ieta++){
	         eta_min = etabins.at(ieta-1);
	         eta_max = etabins.at(ieta);
	         for(int iipt = 1; iipt < invptbinsvec.size();iipt++){
	           invpt_min = invptbinsvec.at(iipt-1);
	           invpt_max = invptbinsvec.at(iipt);
	           if ( TMath::Abs(eta) < eta_min || TMath::Abs(eta) > eta_max) continue;
  		       if ( qOverp < invpt_min        || qOverp > invpt_max       ) continue;

             double invpt_fill = 1./(2*bestftk->getPt());              
             hist_res[0][isIBL][ieta-1][iipt-1]->Fill(d0-bestftk->getIP());
             hist_res[1][isIBL][ieta-1][iipt-1]->Fill(z0-bestftk->getZ0());
             hist_res[2][isIBL][ieta-1][iipt-1]->Fill(eta-bestftk->getEta());
             hist_res[3][isIBL][ieta-1][iipt-1]->Fill(phi-bestftk->getPhi());
			       hist_res[4][isIBL][ieta-1][iipt-1]->Fill(qOverp -invpt_fill);
        	 } // invpt bins
	       } // eta bins 
      } //best ftk loop
    }// matching
  } // end loop over truth tracks





} // processing 


void Process2(Long64_t ientry) {
  t_ftkdata->GetEntry(ientry);

  // collect information from the run-evt from the truth
  Int_t evtnumber_ftk = tracks->eventNumber();
  Int_t runnumber_ftk = tracks->runNumber();

  // first try to start from last number
  for (; ientry2 < t_truth->GetEntries(); ientry2++) {
    t_evtinfo->GetEntry(ientry2);
    if (EventNumber==evtnumber_ftk && RunNumber==runnumber_ftk) break;
  }

  // only enter this loop if we didn't find it above
  if (ientry2 == t_truth->GetEntries()) {
    for (ientry2 = 0; ientry2 < t_truth->GetEntries(); ientry2++) {
      t_evtinfo->GetEntry(ientry2);
      if (EventNumber==evtnumber_ftk && RunNumber==runnumber_ftk) break;
    }
  }

  if (ientry2 == t_truth->GetEntries()) {
    cerr << Form("Mismatch between entries: (%d,%d)",runnumber_ftk,evtnumber_ftk) << endl;
    ientry2=0;
    return;
  }

  t_truth->GetEntry(ientry2);
  if (ientry%10000==0) { Info("Process","Event %lld, (Run,Evt) = (%d,%d)",ientry, runnumber_ftk, evtnumber_ftk);}

  Int_t ntracks = tracks->getNTracks();
  if (Use1stStage == 1)
  ntracks = tracks->getNTracksI();

  FTKBarcodeMM ftkmatchinfo;
  const FTKTrack *curtrack;
  if (Use1stStage != -1) {

    for (Int_t itrk=0;itrk!=ntracks;++itrk) { // loop over the FTK tracks
      curtrack = (Use1stStage ? tracks->getTrackI(itrk) : tracks->getTrack(itrk));
      for ( Int_t ii=0; ii < 12; ++ii){
      	FTKHit hit = curtrack->getFTKHit(ii);
      }
      if (curtrack->getBitmask()&(1<<0) || !(curtrack->getBitmask()&(1<<0)) ) {
	       FTKHit hit = curtrack->getFTKHit(0);
	       bool isIBL = hit.getPlane() == 0;
         if (curtrack->getEventIndex()==0) {
	           // Information on FTK tracks relative to the hard-scattering events
	           // are collected in a vector and later used to build matching maps
	         if (curtrack->getBarcodeFrac()>.5) {
	           ftkmatchinfo.insert(pair<MatchInfo, const FTKTrack*>(MatchInfo(curtrack->getBarcode(),curtrack->getEventIndex()),curtrack));
	         }
	       }
	
      	// apply some basic quality selection
      	double d0 = curtrack->getIP();
        double z0 = curtrack->getZ0();
        double phi = curtrack->getPhi();
        double curv = curtrack->getHalfInvPt();
        double eta = curtrack->getEta();
	      if (d0<d0min || d0>d0max) continue;
	      if (z0<z0min || z0>z0max) continue;
      	if (phi<phimin || phi > phimax) continue;
      	if (curv<-abscurvmax || curv>abscurvmax) continue;
      	if (eta<etamin || eta>etamax) continue;
	
      } // end loop over the FTK tracks
    } // end if not using roads  
  }

  Int_t ntruth = truthTracks->size();
  Int_t ntruth_good(0);

  vector<FTKTruthTrack>::const_iterator itr = truthTracks->begin();
  vector<FTKTruthTrack>::const_iterator itrE = truthTracks->end();
  for (;itr!=itrE;++itr) { // loop over the truth tracks
    const FTKTruthTrack &curtruth = (*itr);

    int barcode = curtruth.getBarcode();
    double px = curtruth.getPX();
    double py = curtruth.getPY();
    double pt = TMath::Sqrt(px*px+py*py);
    double invpt = 1./(2.0*pt);
    double d0 = curtruth.getD0();
    double z0 = curtruth.getZ();
    double curv = curtruth.getQ()*invpt;
    double phi = curtruth.getPhi();
    double eta = curtruth.getEta();
    double qOverp = curv;
    int pdgcode = curtruth.getPDGCode();

    if (barcode>100000 || barcode==0) continue;
    if ( pt < ptmincut ) continue;
    if (d0<d0min || d0>d0max) continue;
    if (z0<z0min || z0>z0max) continue;
    if (curv<-abscurvmax || curv>abscurvmax) continue;
    if (phi<phimin || phi>phimax) continue;
    if (eta<etamin || eta>etamax) continue;
    if (curtruth.getEventIndex()!=0 && curtruth.getQ()==0) continue;
    //EDIT TRACK CRITERIA

    // match the barcode and event index values
    MatchInfo reftruth(barcode,curtruth.getEventIndex());
    pair<FTKBarcodeMM::const_iterator,FTKBarcodeMM::const_iterator> mrange = ftkmatchinfo.equal_range(reftruth);
    if (mrange.first != mrange.second) {
      const FTKTrack *bestftk(0x0);
      for(FTKBarcodeMM::const_iterator ftkI = mrange.first;ftkI!=mrange.second;++ftkI) {
        if (!bestftk){
          bestftk = (*ftkI).second;
        } else if (bestftk->getBarcodeFrac()<(*ftkI).second->getBarcodeFrac()) {
          bestftk = (*ftkI).second;
        }
      }
      
      if (bestftk) {
	       FTKHit hit = bestftk->getFTKHit(0);
	       bool isIBL = hit.getPlane() == 0;
	       for(int ieta = 1; ieta < etabins.size(); ieta++){
	         eta_min = etabins.at(ieta-1);
	         eta_max = etabins.at(ieta);
	         for(int iipt = 1; iipt < invptbinsvec.size();iipt++){
	           invpt_min = invptbinsvec.at(iipt-1);
	           invpt_max = invptbinsvec.at(iipt);
	           if ( TMath::Abs(eta) < eta_min || TMath::Abs(eta) > eta_max) continue;
  		       if ( qOverp < invpt_min        || qOverp > invpt_max       ) continue;

             double invpt_fill = 1./(2*bestftk->getPt());              
             hist_res[0][isIBL][ieta-1][iipt-1]->Fill(d0-bestftk->getIP());
             hist_res[1][isIBL][ieta-1][iipt-1]->Fill(z0-bestftk->getZ0());
             hist_res[2][isIBL][ieta-1][iipt-1]->Fill(eta-bestftk->getEta());
             hist_res[3][isIBL][ieta-1][iipt-1]->Fill(phi-bestftk->getPhi());
			       hist_res[4][isIBL][ieta-1][iipt-1]->Fill(qOverp -invpt_fill);
        	 } // invpt bins
	       } // eta bins 
      } //best ftk loop
    }// matching
  } // end loop over truth tracks





} // processing 




/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////


void width_calculation(){
  double ipt_arr[200];
  std::copy(invptbinsvec.begin(), invptbinsvec.end(), ipt_arr);
  for( int i = 0; i < invptbinsvec.size();i++){ipt_arr[i] = ipt_arr[i] + 0.5*stepsize;}

  for(int itp = 0;itp < nParams ; itp++){
    string trackParam = trackParam_list[itp];
    double trackRange = trackParam_range[itp];

    for(int iibl =0;iibl < 2;iibl++){
      string iblname = iblnames[iibl];

      for(int ieta =0;ieta < etabins.size()-1;ieta++){
        widths.clear();
        width_errors.clear();

        for(int iipt=0;iipt < invptbinsvec.size()-1;iipt++){
	         Double_t par[6]; /* parameter array */
           /* DEFINE FUNCTIONS WE ARE GOING TO FIT.*/
           /* 1. SINGLE GAUSSIAN /CORE GAUSSIAN ESTIMATION*/
           /* 2. TAIL GAUSSIAN ESTIMATION*/
           /* 3. DOUBLE GAUSSIAN */ 

           //hist_std[itp][iibl][ieta][iipt] = hist_res[itp][iibl][ieta][iipt]->GetStdDev();

           //hist_res[itp][iibl][ieta][iipt]->SetRangeUser(-histrange,histrange)
//           double histrange = 4.0*hist_std[itp][iibl][ieta][iipt];
//           hist_res[itp][iibl][ieta][iipt]->SetAxisRange(-histrange,histrange,"X");
           //hist_res[itp][iibl][ieta][iipt]->SetBins(50,-histrange,histrange);
           //std::cout << hist_res[itp][iibl][ieta][iipt]->GetSize() <<std::endl;
//          string rebinned_name = "hist_res_rebinned_" + trackParam + "_" + iblname + "_eta" + to_string(ieta) + "_ipt" + to_string(iipt);
  //        hist_res[itp][iibl][ieta][iipt] = (TH1F*)hist_res_temp[itp][iibl][ieta][iipt]->Rebin(4);
//          hist_res[itp][iibl][ieta][iipt]->SetBins(50,-histrange,histrange);
 //         hist_res[itp][iibl][ieta][iipt] 
//           hist_res[itp][iibl][ieta][iipt]->SetAxisRange(-histrange, histrange, "X");
  //         hist_res[itp][iibl][ieta][iipt]->Rebin(4);


//           double coreRange = hist_std[itp][iibl][ieta][iipt];
//           double tailRange = 0.0;
//           if(coreRange*3.0 > trackRange){ 
//               tailRange = trackRange;   
//           }
//           else{
//               tailRange = 3.0*coreRange;
//           }

//          coreRange = trackRange;
 //         tailRange = trackRange;
          double  histrange = 4.5*hist_std[itp][iibl][ieta][iipt];
          double  corerange = hist_std[itp][iibl][ieta][iipt];
	         TF1 *g1= new TF1 ("m1","gaus",-corerange,corerange);
	         TF1 *g2= new TF1 ("m2","gaus",-histrange,histrange);
	         TF1 *f1= new TF1("double_gaus","gaus(0)+gaus(3)",-histrange,histrange);

           /*  FIT THE FIRST GAUSSIAN */
           /*  NO INITIAL CONDITIONS  */
           //g1->FixParameter(1,0.0);
           hist_res[itp][iibl][ieta][iipt]->Fit(g1,"QR0");

           //double width       = g1->GetParameter(2);
           //double width_error = g1->GetParError(2);
           // EDIT GAUSSIAN
           g2->SetParLimits(3,0.2*par[0],par[0]); // height of second gaussian.
           // g2->SetParLimits(4,none,none);
           g2->SetParLimits(5,par[2],10*par[2]); //width of second gaussian
           //g2->FixParameter(3,none);
           g2->FixParameter(4,0); // // mean of second gaussian.
           // g2->FixParameter(5,none);

           /*  FIT THE SECOND GAUSSIAN */
           /*  INITIAL CONDITIONS GENERATED FROM FIRST GAUSSIAN  */
	         hist_res[itp][iibl][ieta][iipt]->Fit(g2,"QR0");

        	 /* GET PARAMETERS FROM BOTH FITS. */
           /* MAP THE PARAMETERS TO A 6 ELEMENT PAR ARRAY */
	         g1->GetParameters(&par[0]);
	         g2->GetParameters(&par[3]);
	  
           /* Overwrite the parameters from the initial gaussian estimations */
	         par[3] = 0.2*par[0]; // height
	         par[4] = par[1];  // mean 
	         par[5] = 3.0*par[2]; //sigma (width)

           /* SET THE PARAMETERS OF THE DOUBLE GAUSSIAN GAUSS(0) + GAUSS(3) */
	         f1->SetParameters(par);

           /* SET RANGES OF PARAMETERS */
           //f1->SetParLimits(0,none,none);
           //f1->SetParLimits(1,none,none);
           //f1->SetParLimits(2,none,none);
           f1->SetParLimits(3,0.2*par[0],par[0]);
           //f1->SetParLimits(4,none,none);
	         f1->SetParLimits(5,par[2],10*par[2]);

           /* ===========================================================================================*/
           /* FIX PARAMETERS OF F1 */
           //f1->FixParameters(0,none);
           //f1->FixParameters(1,none);
           //f1->FixParameters(2,none);
           //f1->FixParameters(3,none);
           //f1->FixParameters(4,none);
           //f1->FixParameters(5,none);
           /* ===========================================================================================*/

           /* FIT THE HISTOGRAM TO A DOUBLE GAUSSIAN */
        	 hist_res[itp][iibl][ieta][iipt]->Fit(f1,"QR0");

           /* RETREIVE THE CORE AND THE TAIL GAUSSIANS */
           //TF1 *gaussian = hist_res[itp][iibl][ieta][iipt]->GetFunction("m1");
	         //TF1 *gaussian2 = hist_res[itp][iibl][ieta][iipt]->GetFunction("m2");

           /* RETREIVE THE PARAMETERS OF BOTH GAUSSIANS */
           /* CORE GAUSSIAN */
           double n1           = f1->GetParameter(0); double error_n1     = f1->GetParError(0);
           double mu1          = f1->GetParameter(1); double error_mu1    = f1->GetParError(1);
           double sigma1       = f1->GetParameter(2); double error_sigma1 = f1->GetParError(2);
           double var1         = sigma1*sigma1;

           /* TAILS GAUSSIAN */
           double n2           = f1->GetParameter(3);  double error_n2     = f1->GetParError(3);
           double mu2          = f1->GetParameter(4);  double error_mu2    = f1->GetParError(4);
           double sigma2       = f1->GetParameter(5);  double error_sigma2 = f1->GetParError(5);
           double var2         = sigma2*sigma2;


           /* SINGLE GAUSSIAN */
           //double width = g1->GetParameter(2);
           //double width_error = g1->GetParError(2);

          /* DOUBLE GAUSSIAN */
           /* ===================================================================================================================== */
          /* WANT TO RETREIVE THE WIDTH AND THE ERROR IN THAT WIDTH. */
          /* WIDTH FOR A SUM OF RANDOMLY DISTRIBUTED VARIABLES IS A WELL DEFINED CONCEPT. */
          /* FOR TWO GAUSSIANS WE HAVE: */
          //double mubar = (n1*mu1+n2*mu2)/(n1+n2);

          //double width = sqrt(((var1 +mu1*mu1)*n1 + (var2+mu2*mu2)*n2 )/(n1+n2) - mubar)
          /* WHICH FOR VERY SMALL MEANS EQUALS */
          /* double width = sqrt((var1*n1 + var2*n2 )/(n1+n2) ); */
          /* double width = sqrt( weight1*variance1 + weight2*variance2); */
          /* ===================================================================================================================== */

          /* WEIGHTS OF EACH GAUSSIAN FIRST */
          double weight1 = n1/(n1+n2);
          double weight2 = n2/(n1+n2);
          
          double width = sqrt(weight1*var1 +weight2*var2); // APPROXIMATION
          /*CALCULATE THE ERROR OF THE VARIANCE */
          double a_n1     = 1/(n1+n2);
	        double a_n2     = -pow(n1+n2,-2);
  	      double sigma_w1 =  sqrt( a_n1*a_n1*error_n1*error_n1 + a_n2*a_n2*error_n2*error_n2);
 
        	double b_n1 = -pow(n1+n2,-2);
	        double b_n2 = 1/(n1+n2);
	        double sigma_w2 =  sqrt( b_n1*b_n1*error_n1*error_n1 + b_n2*b_n2*error_n2*error_n2);
          /* error^2 = S(1)^2*E(W1)^2 + S(2)^2*E(W2)^2 + W1^2*E(W1)^2 + W2^2*E(W2)^2 */
          double width_error = sqrt(var1*sigma_w1*sigma_w1+ var2*sigma_w2*sigma_w2 +weight1*weight1*sigma_w1*sigma_w1 + weight2*weight2*sigma_w2*sigma_w2);

          widths.push_back(width);
          width_errors.push_back(width_error);
          A_width[itp][iibl][ieta][iipt] = width;
          A_width_errors[itp][iibl][ieta][iipt] = width_error;
        } //invptbins
      } //eta bins
    } // ibl
  } // parameters 
} // void 


/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////



void sqrt_fit(){
  double ipt_arr[200];
  std::copy(invptbinsvec.begin(), invptbinsvec.end(), ipt_arr);
  for( int i = 0; i < invptbinsvec.size();i++){ipt_arr[i] = ipt_arr[i] + 0.5*stepsize;}

  for(int itp = 0;itp < nParams ; itp++){
    string trackParam = trackParam_list[itp];
    double trackRange = trackParam_range[itp];
    for(int iibl =0;iibl < 2;iibl++){
      string iblname = iblnames[iibl];
      for(int ieta =0;ieta < etabins.size()-1;ieta++){
        widths.clear();
        width_errors.clear();

          double width_arr[200];
          double width_errors_arr[200]; 
          for (int i = 0; i < ninvptbins; i++) { width_arr[i]        = A_width[itp][iibl][ieta][i];}
          for (int i = 0; i < ninvptbins; i++) { width_errors_arr[i] = A_width_errors[itp][iibl][ieta][i];}

          Int_t n = ninvptbins;        
          TGraphErrors *graph_linear = new TGraphErrors(n,ipt_arr,width_arr,0,width_errors_arr);
          // EDIT FIT 
          graph_linear->Fit("pol1","q"," ",0.0,invpt_max);
          TF1 *linear_fit = graph_linear->GetFunction("pol1");
          double linear_par0 = linear_fit->GetParameter(0); 
          double linear_par1 = linear_fit->GetParameter(1);

          /* SQUARE ROOT FIT THE ERROR VS CURVATURE */
          TGraphErrors *graph_sqrt = new TGraphErrors(n,ipt_arr,width_arr,0,width_errors_arr);
          TF1 *func = new TF1("sqrtfit",fitf,-invpt_max,invpt_max,2);
          func->SetParameters(width_arr[0],linear_par1*linear_par1);
          double minimumfit = width_arr[13];//std::Min(width_arr[13],width_arr[14]);
          func->SetParLimits(0,0.,1.22*minimumfit*minimumfit);
          func->SetParNames ("a","b");
          graph_sqrt->Fit("sqrtfit","q","",-invpt_max,invpt_max);
          TF1 *sqrt_fit = graph_sqrt->GetFunction("sqrtfit");

          double sqroot_par0 = sqrt_fit->GetParameter(0);
          double sqroot_par1 = sqrt_fit->GetParameter(1);

          A_sq_par0[itp][iibl][ieta] = sqroot_par0;
          A_sq_par1[itp][iibl][ieta] = sqroot_par1;
      } //eta bins
    } // ibl
  } // parameters 
} // void 


/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////


void Pull (Long64_t ientry) {
  t_ftkdata->GetEntry(ientry);

  // collect information from the run-evt from the truth
  Int_t evtnumber_ftk = tracks->eventNumber();
  Int_t runnumber_ftk = tracks->runNumber();

  // first try to start from last number
  for (; ientry2 < t_truth->GetEntries(); ientry2++) {
    t_evtinfo->GetEntry(ientry2);
    if (EventNumber==evtnumber_ftk && RunNumber==runnumber_ftk) break;
  }

  // only enter this loop if we didn't find it above
  if (ientry2 == t_truth->GetEntries()) {
    for (ientry2 = 0; ientry2 < t_truth->GetEntries(); ientry2++) {
      t_evtinfo->GetEntry(ientry2);
      if (EventNumber==evtnumber_ftk && RunNumber==runnumber_ftk) break;
    }
  }

  if (ientry2 == t_truth->GetEntries()) {
    cerr << Form("Mismatch between entries: (%d,%d)",runnumber_ftk,evtnumber_ftk) << endl;
    ientry2=0;
    return;
  }

  t_truth->GetEntry(ientry2);

  if (ientry%10000==0) {
    Info("Process","Event %lld, (Run,Evt) = (%d,%d)",ientry, runnumber_ftk, evtnumber_ftk);
  }

  Int_t ntracks = tracks->getNTracks();
  if (Use1stStage == 1)
  ntracks = tracks->getNTracksI();
  //histontracks_ftk->Fill(ntracks);
  FTKBarcodeMM ftkmatchinfo;
  const FTKTrack *curtrack;
  if (Use1stStage != -1) {

    for (Int_t itrk=0;itrk!=ntracks;++itrk) { // loop over the FTK tracks
      curtrack = (Use1stStage ? tracks->getTrackI(itrk) : tracks->getTrack(itrk));


      for ( Int_t ii=0; ii < 12; ++ii){
        FTKHit hit = curtrack->getFTKHit(ii);
      }
      if (curtrack->getBitmask()&(1<<0) || !(curtrack->getBitmask()&(1<<0)) ) {
         FTKHit hit = curtrack->getFTKHit(0);
         bool isIBL = hit.getPlane() == 0;

         if (curtrack->getEventIndex()==0) {
             // Information on FTK tracks relative to the hard-scattering events
             // are collected in a vector and later used to build matching maps
           if (curtrack->getBarcodeFrac()>.5) {
             ftkmatchinfo.insert(pair<MatchInfo, const FTKTrack*>(MatchInfo(curtrack->getBarcode(),curtrack->getEventIndex()),curtrack));
           }
         }
  
        // apply some basic quality selection
        double eta = curtrack->getEta();
        double d0 = curtrack->getIP();
        double z0 = curtrack->getZ0();
        double phi = curtrack->getPhi();
        double curv = curtrack->getHalfInvPt();
        if (d0<d0min || d0>d0max) continue;
        if (z0<z0min || z0>z0max) continue; 
        if (phi<phimin || phi > phimax) continue;
        if (curv<-abscurvmax || curv>abscurvmax) continue;
        if (eta<etamin || eta>etamax) continue;
  
      } // end loop over the FTK tracks
    } // end if not using roads  
  }
 
  /*
   * Study the list of the truth tracks to collect elements useful for the tracks
   * efficiency
   */
  Int_t ntruth = truthTracks->size();
  Int_t ntruth_good(0);
  //  histontracks_truth->Fill(ntruth);
  vector<FTKTruthTrack>::const_iterator itr = truthTracks->begin();
  vector<FTKTruthTrack>::const_iterator itrE = truthTracks->end();
  for (;itr!=itrE;++itr) { // loop over the truth tracks
    const FTKTruthTrack &curtruth = (*itr);

    int barcode   = curtruth.getBarcode();
    double px     = curtruth.getPX();
    double py     = curtruth.getPY();
    double pt     = TMath::Sqrt(px*px+py*py);
    double invpt  = 1./(2.0*pt);
    double d0     = curtruth.getD0();
    double z0     = curtruth.getZ();
    double curv   = curtruth.getQ()*invpt;
    double phi    = curtruth.getPhi();
    double eta    = curtruth.getEta();
    double qOverp = curv;
    int pdgcode   = curtruth.getPDGCode();


    if (barcode>100000   || barcode==0)      continue;
    if ( pt               < ptmincut )       continue;
    if (d0<d0min         || d0>d0max)        continue;
    if (z0<z0min         || z0>z0max)        continue;
    if (curv<-abscurvmax || curv>abscurvmax) continue;
    if (phi<phimin       || phi>phimax)      continue;
    if (eta<etamin       || eta>etamax)      continue;
    //EDIT TRACK CRITERIA
    if (curtruth.getEventIndex()!=0 && curtruth.getQ()==0) continue;

    // match the barcode and event index values
    MatchInfo reftruth(barcode,curtruth.getEventIndex());

    pair<FTKBarcodeMM::const_iterator,FTKBarcodeMM::const_iterator> mrange = ftkmatchinfo.equal_range(reftruth);

    if (mrange.first != mrange.second) {

      const FTKTrack *bestftk(0x0);
      for(FTKBarcodeMM::const_iterator ftkI = mrange.first;ftkI!=mrange.second;++ftkI) {
        if (!bestftk) {
    // std::cout <<"bestftk check: " << (*ftkI).second << std::endl;
          bestftk = (*ftkI).second;
        } else if (bestftk->getBarcodeFrac()<(*ftkI).second->getBarcodeFrac()) {
          bestftk = (*ftkI).second;
        }
      }
      
      if (bestftk) {
	       FTKHit hit = bestftk->getFTKHit(0);
	       bool isIBL = hit.getPlane() == 0;
	       //	       if (isIBL){
         for(int ieta = 1; ieta < etabins.size(); ieta++){
             eta_min = etabins.at(ieta-1);
             eta_max = etabins.at(ieta);

             for(int iipt = 1; iipt < invptbinsvec.size();iipt++){
               invpt_min = invptbinsvec.at(iipt-1);
               invpt_max = invptbinsvec.at(iipt);

               if ( TMath::Abs(eta) < eta_min || TMath::Abs(eta) > eta_max ) continue;
               if (qOverp < invpt_min || qOverp > invpt_max) continue;
               double widthcalc[5];
               
               for(int itp = 0; itp < nParams; itp++){
                 widthcalc[itp] = TMath::Sqrt(A_sq_par0[itp][isIBL][ieta-1] + A_sq_par1[itp][isIBL][ieta-1]*qOverp*qOverp ); 
               }
               double d0pull     = (d0-bestftk->getIP())/widthcalc[0];
               double z0pull     = (z0-bestftk->getZ0())/widthcalc[1];
               double etapull    = (eta-bestftk->getEta())/widthcalc[2];
               double phipull    = (phi-bestftk->getPhi())/widthcalc[3];  
               double invpt_fill = 1./(2.0*bestftk->getPt());
               double iptpull    = (qOverp-invpt_fill)/widthcalc[4];

               pull_res[0]->Fill(d0pull);
               pull_res[1]->Fill(z0pull);
               pull_res[2]->Fill(etapull);
               pull_res[3]->Fill(phipull);
               pull_res[4]->Fill(iptpull);

              } // invpt bins
            } // eta bins 
	 //	  }// ibl bin

      } //best ftk loop
    }// matching
  } // end loop over truth tracks

} // processing 


/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void Terminate(std::string& outputname) {

  Info("Terminate","Adding the histograms to the file: %s", outputname.c_str());
  TFile *ofile = TFile::Open(outputname.c_str(), "recreate");

  double ipt_arr[200];
  std::copy(invptbinsvec.begin(), invptbinsvec.end(), ipt_arr);
  for( int i = 0; i < invptbinsvec.size();i++){ipt_arr[i] = ipt_arr[i] + 0.5*stepsize;}

  ofstream sqrtIBL("sqrtIBL.log");
  ofstream sqrtnoIBL("sqrtnoIBL.log");
  ofstream linearIBL("linearIBL.log");
  ofstream linearnoIBL("linearnoIBL.log");
  for(int itp = 0;itp < nParams ; itp++){
    string trackParam = trackParam_list[itp];
    double trackRange = trackParam_range[itp];
    string trackParam_unit = trackParam_units[itp];

    for(int iibl =0;iibl < 2;iibl++){
      string iblname = iblnames[iibl];

      for(int ieta =0;ieta < etabins.size()-1;ieta++){
    	widths.clear();
	    width_errors.clear();

    	for(int iipt=0;iipt < invptbinsvec.size()-1;iipt++){
        TCanvas *c = new TCanvas("c", "c", 800, 650);

	      Double_t par[6]; /* parameter array */
//        hist_std[itp][iibl][ieta][iipt] = hist_res[itp][iibl][ieta][iipt]->GetStdDev();
      //std::cout << hist_res[itp][iibl][ieta][iipt]->GetSize() <<std::endl;

      // // EDIT GAUSSIAN
      double  histrange = 4.5*hist_std[itp][iibl][ieta][iipt];
      double  corerange = hist_std[itp][iibl][ieta][iipt];
	    TF1 *g1= new TF1 ("m1","gaus",-corerange,corerange);
	    TF1 *g2= new TF1 ("m2","gaus",-histrange,histrange);
	    TF1 *f1= new TF1("double_gaus","gaus(0)+gaus(3)",-histrange,histrange);

      TF1 *core = new TF1 ( "core","gaus",-histrange,histrange);
      TF1 *tails = new TF1 ( "tails","gaus",-histrange,histrange);

      f1->SetLineColor(kRed+1);
      core->SetLineColor(kBlue+1);
      tails->SetLineColor(kGreen+2);

      f1->SetLineStyle(1);
      core->SetLineStyle(2);
      tails->SetLineStyle(4);

          hist_res[itp][iibl][ieta][iipt]->Fit(g1,"QR0");


           g2->SetParLimits(3,0.2*par[0],par[0]); // height of second gaussian.
           // g2->SetParLimits(4,none,none);
           g2->SetParLimits(5,par[2],10*par[2]); //width of second gaussian
           //g2->FixParameter(3,none);
           g2->FixParameter(4,0); // // mean of second gaussian.
           // g2->FixParameter(5,none);

           /*  FIT THE SECOND GAUSSIAN */
           /*  INITIAL CONDITIONS GENERATED FROM FIRST GAUSSIAN  */
	         hist_res[itp][iibl][ieta][iipt]->Fit(g2,"QR0");

        	 /* GET PARAMETERS FROM BOTH FITS. */
           /* MAP THE PARAMETERS TO A 6 ELEMENT PAR ARRAY */
	         g1->GetParameters(&par[0]);
	         g2->GetParameters(&par[3]);
	  
           /* Overwrite the parameters from the initial gaussian estimations */
	         par[3] = 0.2*par[0]; // height
	         par[4] = par[2];  // mean 
	         par[5] = 3.0*par[2]; //sigma (width)

           /* SET THE PARAMETERS OF THE DOUBLE GAUSSIAN GAUSS(0) + GAUSS(3) */
	         f1->SetParameters(par);

           /* SET RANGES OF PARAMETERS */
           f1->SetParLimits(0,0.0,5000.0);
           //f1->SetParLimits(1,none,none);
           //f1->SetParLimits(2,none,none);
           f1->SetParLimits(3,0.2*par[0],par[0]);
           //f1->SetParLimits(4,none,none);
	         f1->SetParLimits(5,par[2],5*par[2]);

           /* ===========================================================================================*/
           /* FIX PARAMETERS OF F1 */
           //f1->FixParameters(0,none);
           f1->FixParameter(1,0.0);
           //f1->FixParameters(2,none);
           //f1->FixParameters(3,none);
           f1->FixParameter(4,0.0);
           //f1->FixParameters(5,none);
           /* ===========================================================================================*/

           /* FIT THE HISTOGRAM TO A DOUBLE GAUSSIAN */
        	 hist_res[itp][iibl][ieta][iipt]->Fit(f1,"QR");
           core->SetParameters(f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2));
           tails->SetParameters(f1->GetParameter(3),f1->GetParameter(4),f1->GetParameter(5));





       //g1->FixParameter(1,0.0);
   //    hist_res[itp][iibl][ieta][iipt]->Fit(g1,"QR");
   //    hist_res[itp][iibl][ieta][iipt]->Fit(g2,"QR+");


//        g1->GetParameters(&par[0]);
	      //	  hist_res[itp][iibl][ieta][iipt]->Draw();
	      //          g1->SetParLimits(5,0.0,5*par[3]);
	      //          double height = par[0];
        /////////g2->SetParLimits(0,5.,0.2*par[0]);
        /////////g2->FixParameter(4,0);
        /////////g2->SetParLimits(2,0.0,3*par[3]);
	      //	  g2->FixParameter(1,0);
        //hist_res[itp][iibl][ieta][iipt]->Fit(g2,"R");
	      //hist_res[itp][iibl][ieta][iipt]->Draw("sames");
        /* get parameters from the fit, it first fits & takes the parameter from there */

        //g2->GetParameters(&par[3]);

        //par[3] = 0.1*par[3];
        //par[4] = 0.0;
        //par[5] = par[5];

        //f1->SetParameters(par);
//	      f1->SetParLimits(5,0.0,5*par[3]);

        //f1->SetParLimits(3,5.,0.2*par[0]);
	      //f1->SetParLimits(5,0.0,3*par[3]);
	      //	  f1->FixParameter(4,0);
//        hist_res[itp][iibl][ieta][iipt]->Fit(f1,"R");

       //TF1 *gaussian = hist_res[itp][iibl][ieta][iipt]->GetFunction("single_gaus");
       // TF1 *gaussian2 = hist_res[itp][iibl][ieta][iipt]->GetFunction("m2");
        //TF1 *dbgaussian = hist_res[itp][iibl][ieta][iipt]->GetFunction("double_gaus");

        //double ncore = f1->GetParameter(0);
        //double ntails = f1->GetParameter(3);
        //double ntotal = ncore + ntails;
        //double sigmacore = f1->GetParameter(2);
        //double sigmatails = f1->GetParameter(5);

        //double wcore = ncore / ntotal;
        //double wtails = ntails / ntotal;
        //double wtotal = wtails + wcore;
        //double width_new = sqrt(wcore*sigmacore*sigmacore + wtails*sigmatails*sigmatails);
        //double width = width_new;
        //double width_error = g1->GetParError(2);

	      //if(iibl == 0){
        //  width = sigmacore;
        //}else{
	      //  double n1 = f1->GetParameter(0);
	      //  double mu1= f1->GetParameter(1);
	      //  double s1 = f1->GetParameter(2);
	      //  double n2 = f1->GetParameter(3);
	      //  double mu2= f1->GetParameter(4);
	      //  double s2 = f1->GetParameter(5);

    	  //  double var1 = s1*s1;
	      //  double var2 = s2*s2;
	      //  double mubar = (n1*mu1+n2*mu2)/(n1+n2);

	        //width = sqrt(((var1 +mu1*mu1)*n1 + (var2+mu2*mu2)*n2 )/(n1+n2) - mubar);
	      //  width = sqrt(wcore*sigmacore*sigmacore + wtails*sigmatails*sigmatails);
	          //	    width = wcore*sigmacore + wtails*sigmatails; //sqrt(wcore*sigmacore*sigmacore + wtails*sigmatails*sigmatails);
       // }

        //widths.push_back(width);
        //width_errors.push_back(width_error);
	      //A_width[itp][iibl][ieta][iipt] = width;

    	  gStyle->SetOptFit(1111);
	      //std::cout << "after gaussian fit" << std::endl;
	      hist_res[itp][iibl][ieta][iipt]->Draw();
        core->Draw("SAME");
        tails->Draw("SAME");

        // currently draw plot 1 and plot 2
        // 

      TLatex* atlas_title = new TLatex();
      atlas_title->SetTextAlign(11);
      atlas_title->SetTextSize(0.04);
      atlas_title->SetTextFont(42);
      atlas_title->SetNDC();
      std::string tempname ="#bf{#it{ATLAS}} Simulation " + iblname;
      const char *atlas_title_name = tempname.c_str();
      atlas_title->DrawLatex(.15, 0.83, atlas_title_name);


      eta_min = etabins.at(ieta);
      eta_max = etabins.at(ieta+1);
      invpt_min = invptbinsvec.at(iipt);
      invpt_max = invptbinsvec.at(iipt+1);


        TLatex* etapt = new TLatex();
      etapt->SetTextAlign(11);
      etapt->SetTextSize(0.04);
      etapt->SetTextFont(42);
      etapt->SetNDC();

      std::string etapttempname ="#eta in [" + to_string(round(eta_min,2)).substr(0,3) + "," + to_string(round(eta_max,2)).substr(0,5) + "] and 1/2p_{T} in [" + to_string(round(invpt_min*1e3,5)).substr(0,4) + "," + to_string(round(invpt_max*1e3,2)).substr(0,4) + "]";
      
      const char *etapt_name = etapttempname.c_str();
      etapt->DrawLatex(.15, 0.93, etapt_name);

//        TF1 *fb2 = new TF1("fa3","gaus",-trackRange,trackRange);
//        fb2->SetParameters(g1->GetParameter(0),g1->GetParameter(1),g1->GetParameter(2));
//        fb2->Draw("L SAME");
//        TF1 *single_gaussian = hist_res[itp][iibl][ieta][iipt]->GetFunction("single_gaus");
 //       single_gaussian->Draw("same");
//        hist_res[itp][iibl][ieta][iipt]->Draw("single_gaus","","SAME")
	      //	  hist_res[itp][iibl][ieta][iipt]->Draw("double_gaus","","SAME");

    	  c->Update();
	      TImage *img = TImage::Create();

	      img->FromPad(c);

    	  string hist_res_name_string = "hist_res" + trackParam + "_" + iblname + "_eta" + to_string(ieta) + "_ipt" + to_string(iipt)+".png";
	      TString hist_res_name("updatefolder/" + hist_res_name_string);

	      img->WriteImage(hist_res_name);
	      ofile->Add(hist_res[itp][iibl][ieta][iipt]);
	    }


	    double width_arr[200];
	    double width_errors_arr[200];

      for (int i = 0; i < ninvptbins; i++) { width_arr[i]        = A_width[itp][iibl][ieta][i];}
      for (int i = 0; i < ninvptbins; i++) { width_errors_arr[i] = A_width_errors[itp][iibl][ieta][i];}

	    Int_t n = ninvptbins;

	    TCanvas *clinear = new TCanvas;
	    TGraphErrors *graph_linear = new TGraphErrors(n,ipt_arr,width_arr,0,width_errors_arr);
	    graph_linear->Fit("pol1","q"," ",0.0,invpt_max);
	    TF1 *linear_fit = graph_linear->GetFunction("pol1");
	    double linear_par0 = linear_fit->GetParameter(0);
	    double linear_par1 = linear_fit->GetParameter(1);
	    std::string linear_name = "linear_" + trackParam + "_" + iblname + "_eta" + std::to_string(ieta);
	    const char *graph_name = linear_name.c_str();

	    TString graph_title(";Q/2p_{T} (MeV^{-1});#sigma(" + trackParam + ")" + trackParam_unit );
	    graph_linear->SetNameTitle(graph_name,graph_title);
	    graph_linear->SetFillColor(2);
	    graph_linear->Draw("ep");
  
	    TString Yaxisname("#sigma(" + trackParam + ")"+ trackParam_unit );
	    TString Xaxisname("Q/2p_{T} [MeV]^{-1}");
	    graph_linear->GetYaxis()->SetTitle(Yaxisname);
	    graph_linear->GetXaxis()->SetTitle(Xaxisname);
	    //        graph_linear->SetLineWidth(1);
	    graph_linear->SetMarkerStyle(20);
	    graph_linear->SetMarkerColor(2);
	    //gStyle->SetOptFit(1111);
	    TLatex* TitleATLAS = new TLatex();
	    //TString atlastitle("ATLAS Simulation, " + iblname); 
	    //	TitleATLAS->DrawLatex(0,0.75*width_arr[0],atlastitle);


    	TImage *imglinear = TImage::Create();
	    graph_linear->Draw("ep");
	    clinear->Update();

    	TString linear_save("updatefolder/" +linear_name + ".png");
	    imglinear->FromPad(clinear);
	    imglinear->WriteImage(linear_save);

	    ofile->Add(graph_linear);

	    TCanvas *csqroot = new TCanvas("cc", "cc", 800, 650);
	    TGraphErrors *graph_sqrt = new TGraphErrors(n,ipt_arr,width_arr,0,width_errors_arr);
	    TImage *imgsqroot = TImage::Create();
	    csqroot->cd();
	    //	graph_sqrt->Draw("ep");
	    //gStyle->SetOptFit(1111);

    	std::string sqroot_name = "sqroot_" + trackParam + "_" + iblname + "_eta" + std::to_string(ieta);
      TString sqroot_save("updatefolder/" + sqroot_name + ".png");

	    TF1 *func = new TF1("sqrtfit",fitf,-invpt_max,invpt_max,2);
      func->SetParameters(width_arr[0],linear_par1*linear_par1);
      double minimumfit = width_arr[13];//std::Min(width_arr[13],width_arr[14]);
       func->SetParLimits(0,0,1.22*minimumfit*minimumfit);
      //	func->SetParameters(width_arr[0],linear_par1*linear_par1);
	    func->SetParNames ("a","b");

	    graph_sqrt->Fit("sqrtfit","q", "",-invpt_max,invpt_max);
	    TF1 *sqrt_fit = graph_sqrt->GetFunction("sqrtfit");
	    double sqroot_par0 = sqrt_fit->GetParameter(0);
	    double sqroot_par1 = sqrt_fit->GetParameter(1);

      const char *graph_sqroot_name = sqroot_name.c_str();
      TString sqroot_title(";Q/2p_{T} (MeV^{-1});#sigma(" + trackParam + ")" + trackParam_unit );
      graph_sqrt->SetNameTitle(graph_sqroot_name,sqroot_title);
	    //	graph_sqrt->Draw("ep");
	    graph_sqrt->SetFillColor(1);
	    graph_sqrt->SetLineWidth(1);
      graph_sqrt->SetMarkerStyle(20);
      graph_sqrt->SetMarkerColor(2);
	    //	gStyle->SetOptFit(1111);
	    //	graph_sqrt->Draw("ep");

	    TString YaxisSQname("#sigma(" + trackParam + ")"+ trackParam_unit );
	    TString XaxisSQname("Q/2p_{T} [MeV]^{-1}");

	    graph_sqrt->GetYaxis()->SetTitle(YaxisSQname);
	    graph_sqrt->GetXaxis()->SetTitle(XaxisSQname);
      graph_sqrt->SetTitle("");
      graph_sqrt->Draw("AEP");
      TLatex* tl = new TLatex();
      tl->SetTextAlign(11);
      tl->SetTextSize(0.04);
      tl->SetTextFont(42);
      tl->SetNDC();
      std::string tempname2 ="#bf{#it{ATLAS}} Simulation " + iblname;
      const char *tl_name = tempname2.c_str();
      tl->DrawLatex(.15, 0.83, tl_name);


      eta_min = etabins.at(ieta);
      eta_max = etabins.at(ieta+1);



        TLatex* etapt = new TLatex();
      etapt->SetTextAlign(11);
      etapt->SetTextSize(0.04);
      etapt->SetTextFont(42);
      etapt->SetNDC();

      std::string etapttempname ="#eta in [" + to_string(round(eta_min,2)).substr(0,3) + "," + to_string(round(eta_max,2)).substr(0,3) + "]" ;
      
      const char *etapt_name = etapttempname.c_str();
      etapt->DrawLatex(.15, 0.93, etapt_name);




      csqroot->Update();
	    imgsqroot->FromPad(csqroot);
	    imgsqroot->WriteImage(sqroot_save);


	    gStyle->SetOptFit(1111);
	    csqroot->Update();
	    imgsqroot->FromPad(csqroot);
	    //	imgsqroot->Update();
	    imgsqroot->WriteImage(sqroot_save);

	    ofile->Add(graph_sqrt);



	    if ( iibl == 0 ) {
	      sqrtnoIBL << "nBLConsts.set(FTKTrackParam::" << std::left << std::setw(3) << trackParam << ", " << ieta << ", sqroot, " << std::scientific << sqroot_par0 <<  ",  " << std::scientific <<  TMath::Abs(sqroot_par1) << " );" << endl;
	    }else if (iibl == 1){
	      sqrtIBL << "nomConsts.set(FTKTrackParam::" <<  std::left << std::setw(3) << trackParam << ", " << ieta << ", sqroot, " << std::scientific << sqroot_par0 <<  ",  " << std::scientific <<  TMath::Abs(sqroot_par1) << " );" << endl;
	    }

      if ( iibl == 0 ) {
        linearnoIBL << "nBLConsts.set(FTKTrackParam::"<<  std::left << std::setw(3)  <<  trackParam << ", " << ieta << ", linear, " << std::scientific << linear_par0 <<  ",  " << std::scientific << linear_par1 << " );" << endl;
      }else if (iibl == 1){
        linearIBL << "nomConsts.set(FTKTrackParam::" <<  std::left << std::setw(3)  <<  trackParam << ", " << ieta << ", linear, " << std::scientific <<linear_par0 <<  ",  " <<  std::scientific << linear_par1 << " );" << endl;
      }

      }
    }
  }



   for( int itp = 0; itp < nParams; itp++){
     pull_res[itp]->Fit("gaus","q","",-5.,5.);
     TF1 *gaussian = pull_res[itp]->GetFunction("gaus");

     TCanvas *c = new TCanvas;
     c->SetLogy();
     pull_res[itp]->Draw();
     c->Update();
     TImage *img = TImage::Create();
     img->FromPad(c);
     string trackParam = trackParam_list[itp];
     string hist_res_name_string = "updatefolder/pull_res" + trackParam +".png";
     TString hist_res_name(hist_res_name_string);
     gStyle->SetOptFit();
     c->Update();
     img->WriteImage(hist_res_name);

     
     ofile->Add(pull_res[itp]);
   }

  ofile->ls();
  ofile->Write();
  ofile->Close();

  if ( psfile!="" ) printpsfile(psfile, ofile);

}


/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////


void printpsfile(std::string psfile, TFile* f) {
   // print the relevant efficiency curves
   // into one PS file
   gROOT->SetBatch();
   gStyle->SetEndErrorSize(0);

   f = TFile::Open(f->GetName(),"READ");

   TCanvas c("Efficiency","Efficiency",800,600);
   
   std::vector<TGraphAsymmErrors*> graphs{
      (TGraphAsymmErrors*)f->Get("greffeff"), 
	 (TGraphAsymmErrors*)f->Get("greffd0"), 
	 (TGraphAsymmErrors*)f->Get("greffz0"),
	 (TGraphAsymmErrors*)f->Get("greffcurv"), 
	 (TGraphAsymmErrors*)f->Get("greffeta"),
	 (TGraphAsymmErrors*)f->Get("greffetaabs"),
	 (TGraphAsymmErrors*)f->Get("greffphi"), 
	 (TGraphAsymmErrors*)f->Get("greffpt"),
	 (TGraphAsymmErrors*)f->Get("greffpt_lg"),
	 (TGraphAsymmErrors*)f->Get("greffptlo_lg"),
	 };
   
   unsigned int n = 0;
   for ( auto gr : graphs ){
      gPad->SetLogx(TString(gr->GetName()).Contains("_lg"));
      gr->SetMinimum(0.5); // TDR: 0.5
      gr->SetMaximum(1.0); // TDR: 1.0
      gr->SetMarkerStyle(22);
      gr->SetMarkerColor(1);
      gr->GetXaxis()->SetMoreLogLabels();
      gr->Draw("AEP"); // TDR: AEP

      if (n==0 ) c.Print(TString(psfile)+"(");
      else if (n==graphs.size()-1 ) c.Print(TString(psfile)+")");
      else c.Print(psfile.c_str());
      n++;
   }
   cout<<" Wrote .ps-file: "<<psfile<<endl;

}


/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////


int main(int argc, char **argv) {
  int events;
  std::string output;
  std::vector<std::string> files;
  ptmincut = 1000;
  dx = -0.5;
  dy = -0.5;

  try {
    std::string appName = boost::filesystem::basename(argv[0]);
    // Define and parse the program options
    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()
      ("help,h", "Print this help message")
      ("output,o", po::value<std::string>(&output)->default_value("ftk_efficiency_test2.root"), "The name out of the output file")
      ("events,e", po::value<int>(&events)->default_value(-1), "The number of events to run over. Set to -1 to use all events")
      ("tower,t", po::value<int>(&towerNumber)->default_value(-1), "the number of the tower in ftkdata output data. Use -1 if you have fully merged output")
      // ("use-first-stage", po::bool_switch(&Use1stStage)->default_value(false), "Use events 1st stage tracks instead of 2nd stage tracks")
      ("use-first-stage", po::value<int>(&Use1stStage)->default_value(0), "-1: Use roads, 1: Use 1st stage tracks, 0(default): Use 2nd stage tracks")
      ("files", po::value<std::vector<std::string> >(&files)->multitoken(), "FTK NTUP files")
       ("ptmincut", po::value<double>(&ptmincut), "min pt cut on truth tracks")
       ("dx", po::value<double>(&dx)->default_value(-0.5), "dx")
       ("dy", po::value<double>(&dy)->default_value(-0.5), "dx")
       ("psfile", "Produce postscript file with efficieny plots");

    po::variables_map vm;
    try
    {
      po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
      // --help option
      if (vm.count("help"))
      {
        cout << desc << "\n";
        return 1;
      }
      po::notify(vm);
    }
    catch(std::exception& e)
    {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      return 1;
    }

    if ( vm.count("psfile") ) 
       psfile = TString(output).ReplaceAll(".root",".ps");

    Init();
    t_ftkdata = new TChain("ftkdata");
    t_truth = new TChain("truthtracks");
    t_evtinfo = new TChain("evtinfo");

    // add the input files
    std::vector<std::string>::const_iterator files_it = files.begin();
    for (; files_it != files.end(); ++files_it) {
      t_ftkdata->Add((*files_it).c_str());
      t_truth->Add((*files_it).c_str());
      t_evtinfo->Add((*files_it).c_str());
    }

    t_ftkdata->AddFriend(t_truth);
    t_ftkdata->AddFriend(t_evtinfo);
    t_evtinfo->SetBranchAddress("RunNumber",&RunNumber);
    t_evtinfo->SetBranchAddress("EventNumber",&EventNumber);
    //  t_ftkdata->SetBranchAddress("FTKMergedTracksStream",&tracks);
    if (towerNumber != -1) {
      t_ftkdata->SetBranchAddress(Form("FTKMergedTracksStream%d",towerNumber),&tracks);
    } else {
      t_ftkdata->SetBranchAddress(Form("FTKMergedTracksStream"),&tracks);
    }

    if (Use1stStage == -1) {
      if (towerNumber != -1)
    	t_ftkdata->SetBranchAddress(Form("FTKMergedRoadsStream%d",towerNumber),&roads);
      else
    	t_ftkdata->SetBranchAddress(Form("FTKMergedRoadsStream"),&roads);
    }

    t_truth->SetBranchAddress("TruthTracks",&truthTracks);

    Int_t nloop = t_ftkdata->GetEntries();
    if (events > 0) {
      nloop = events;
    }
    if (nloop > t_ftkdata->GetEntries()) {
      nloop = t_ftkdata->GetEntries();
    }
    for (int ientry = 0; ientry < nloop; ++ientry) {
      Process(ientry);
    }
    Init2();
    for (int ientry = 0; ientry < nloop; ++ientry) {
      Process2(ientry);
    }
    width_calculation();
    sqrt_fit();
    for (int ientry = 0; ientry < nloop; ++ientry){
      Pull(ientry);
    }
    Terminate(output);
  } // end try
  catch(std::exception& e)
  {
    std::cerr << "Unhandled Exception reached the top of main: "
              << e.what() << ", application will now exit" << std::endl;
    return 1;
  }
  return 0;
}



         // if(0){
        //   std::cout << "==================================" <<std::endl; 
        //   std::cout << "=============  TRUTH  ============" <<std::endl; 
        //   std::cout << "==================================" <<std::endl;
        //   std::cout << "PX: " << px <<std::endl;
        //   std::cout << "PY: " << py <<std::endl;
        //   std::cout << "PT: " << TMath::Abs(pt) <<std::endl;
        //   std::cout << "INVPT: " << invpt <<std::endl;
        //   std::cout << "qOverp: " << qOverp <<std::endl;  
        //   std::cout << "d0: " << d0 <<std::endl;
        //   std::cout << "z0: " << z0 <<std::endl;
        //   std::cout << "curv: " << curv <<std::endl;
        //   std::cout << "phi: " << phi <<std::endl;
        //   std::cout << "eta: " << eta <<std::endl;
        //   std::cout << "pdgcode: " << pdgcode <<std::endl;
        //   std::cout << "==================================" <<std::endl; 
        //   std::cout << "=============   MATCHED   ========" <<std::endl; 
        //   std::cout << "==================================" <<std::endl;
        //   std::cout << "PT: "    << TMath::Abs(bestftk->getPt()) <<std::endl;
        //   std::cout << "INVPT: " << 1/(2.0*bestftk->getPt()) <<std::endl;
        //   std::cout << "d0: "    << bestftk->getIP() <<std::endl;
        //   std::cout << "z0: "    << bestftk->getZ0() <<std::endl;
        //   std::cout << "phi: "   << bestftk->getPhi() <<std::endl;
        //   std::cout << "eta: "   << bestftk->getEta()<<std::endl;
        //   std::cout << "==================================" <<std::endl; 
        //   std::cout << "===========  MAXIMUM EXTENT ======" <<std::endl; 
        //   std::cout << "==================================" <<std::endl; 
        //   std::cout << "abscurvmax: " << abscurvmax  <<std::endl; 
        //   std::cout << "d0min     : " << d0min       <<std::endl; 
        //   std::cout << "z0min     : " << z0min       <<std::endl; 
        //   std::cout << "phimin    : " << phimin      <<std::endl; 
        //   std::cout << "etamin    : " << etamin      <<std::endl; 
        //   std::cout << "==================================" <<std::endl; 
        //   std::cout << "===========  MAXIMUM EXTENT ======" <<std::endl; 
        //   std::cout << "==================================" <<std::endl; 
        //   std::cout << "abscurvmax: [0.5e-3]"  <<std::endl; 
        //   std::cout << "d0min     : [MM]    "  <<std::endl; 
        //   std::cout << "z0min     : [MM]    "  <<std::endl; 
        //   std::cout << "phimin    : [Rad]   "  <<std::endl; 
        //   std::cout << "etamin    : [Eta]   "  <<std::endl; 
        //   std::cout << "==================================" <<std::endl; 
        //   std::cout << "======= TRUTH - MATCHED   ========" <<std::endl; 
        //   std::cout << "==================================" <<std::endl;
        //   std::cout << "PT - matched: "    << TMath::Abs(pt) - TMath::Abs(bestftk->getPt()) <<std::endl;
        //   std::cout << "qOverp - INVPT: "  << qOverp - 1/(2.0*bestftk->getPt()) <<std::endl;
        //   std::cout << "curv - INVPT: "    << curv - 1/(2.0*bestftk->getPt()) <<std::endl;
        //   std::cout << "d0 - matched: "    << d0  - bestftk->getIP() <<std::endl;
        //   std::cout << "z0 - matched: "    << z0  - bestftk->getZ0() <<std::endl;
        //   std::cout << "phi- matched: "    << phi - bestftk->getPhi() <<std::endl;
        //   std::cout << "eta- matched: "    << eta - bestftk->getEta() <<std::endl;
        //   std::cout << "==================================" <<std::endl;
        //   std::cout << "isIBL: "    << isIBL <<std::endl;
        //   std::cout << "==================================" <<std::endl;
         // }


/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
