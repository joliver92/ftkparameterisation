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
#include "TRatioPlot.h"
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

////////////////////////////////// <SQUARE ROOT FIT DEFINITION>
Double_t fitf(Double_t *x,Double_t *par) {
  Double_t fitval = TMath::Sqrt(par[0] + par[1]*x[0]*x[0]);

  return fitval;
}
////////////////////////////////// </SQUARE ROOT FIT DEFINITION>


////////////////////////////////// <EXPLICIT GAUSSIAN FUNCTION>
double gaussFunc(double *x, double *par) {
   double s;
   s = par[0]*exp(-(x[0]-par[1])*(x[0]-par[1])/par[3]/par[3]/2   );
   return s;
}
////////////////////////////////// </EXPLICIT GAUSSIAN FUNCTION>

////////////////////////////////// <EXPLICIT GAUSSIAN FUNCTION>
double doublegaussfunc(double *x, double *par) {
   double s;
   Double_t argcore = 0;
   if (par[2] != 0) argcore = (x[0] - par[1])/par[2];
   Double_t argtail = 0;
   if (par[4] != 0) argtail = (x[0] - par[3])/par[4];
   //core = par[0]*exp(-0.5*argcore*argcore);
   //tail = 0.11*par[0]*exp(-0.5*argtail*argtail);
   s = par[0]*exp(-0.5*argcore*argcore) + 0.11*par[0]*exp(-0.5*argtail*argtail);
//   s = par[0]*exp(-(x[0]-par[1])*(x[0]-par[1])/par[2]/par[2]/2   ) + 0.11*par[0]*exp(-(x[0]-par[3])*(x[0]-par[3])/par[4]/par[4]/2   )    ;
   return s;
}
////////////////////////////////// </EXPLICIT GAUSSIAN FUNCTION>



////////////////////////////////// <INITIALISE>
void Init() {
  // define the boundaries for
  int nParams= 5;
  nbins      = 500;
  nbinspull  = 80; // nbinspull 80 default
  maxntracks = 500;
  d0min      = -2.         ; d0max     = 2.         ;
  z0min      = -100.       ; z0max     = 100.       ;
  phimin     = -TMath::Pi(); phimax    = TMath::Pi();
  etamin     = -2.4        ; etamax    = 2.4        ;
  invpt_min  = -0.5e-3     ; invpt_max = 0.5e-3     ;
  double invptmin = -0.5e-3; double invptmax = 0.5e-3;
  pTmin      = 1000           ; ptmax     = 100000  ;
  abscurvmax = 0.5e-3;

 Dd0 = 0.8;    Drelpt = .00010; Dcurv = .15; Dphi = .05; Dz0 = 4.0; Deta = .040;

  TP_ranges_max.assign({d0max,z0max,etamax,phimax,invptmax});
  TP_ranges_min.assign({d0min,z0min,etamin,phimin,invptmin});
  TP_ranges_res.assign({Dd0  ,Dz0  ,Deta  ,Dphi  ,Drelpt});
  etabin_min = 0;   etabin_max = 0;
  invptbin_min = 0; invptbin_max = 0;


  trackParam_range[0] = Dd0;
  trackParam_range[1] = Dz0;
  trackParam_range[2] = Deta;
  trackParam_range[3] = Dphi;
  trackParam_range[4] = Drelpt;
  
  ninvptbins = 26; //26 //must be even.

  middlebin1 = std::floor((ninvptbins-1)/2.0);
  middlebin2= std::ceil((ninvptbins-1)/2.0);
  stepsize = (invpt_max - invpt_min)/ninvptbins;
  invptval = invpt_min;
 
  while(invptval <= invpt_max+0.5*stepsize){
    invptbinsvec.push_back(invptval);
    invptval = invptval +stepsize;
   
  }
  //ipt_arr[200];
  std::copy(invptbinsvec.begin(), invptbinsvec.end(), ipt_arr);
  for( unsigned int i = 0; i < invptbinsvec.size();i++){ipt_arr[i] = ipt_arr[i] + 0.5*stepsize;}
 
  //double etabinsarray[] = {0.0,0.5,1.0,1.5};
  //,2.0,2.5};
 double etabinsarray[] = {0.0,0.5,1.0,1.5,2.0,2.5};
 //double etabinsarray[] = {0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5};
 
  for (auto &bin :etabinsarray){etabins.push_back(bin);}
  iblnames[0] = "noIBL";
  iblnames[1] = "IBL";
  nibl = 2;
  for(int itp = 0; itp < nParams;itp++){
     string trackParam = trackParam_list[itp];
     double range      = trackParam_range[itp];
     TString TtrackParam(trackParam);
     TString hist_res_title(";#Delta" + trackParam + "[ftk-truth](mm);N Tracks");
     TString TtrackParam_title(";" + trackParam + "(ftk-truth)/#sigma(ftk);N Tracks");
        

    pull_res[itp] = new TH1F( TtrackParam,TtrackParam_title,nbinspull,-5,5);
    TString TtrackParamSG(trackParam + "SG");
    pull_res_sg[itp] = new TH1F( TtrackParamSG,TtrackParam_title,nbinspull,-5,5);
    for (int ipulltype =0;ipulltype < 3;ipulltype++){
      for (int itp2 = 0; itp2 < nParams;itp2++){
        string trackParam2 = trackParam_list[itp2];

        TString pull_title(    ";" + trackParam2 + ";Rel Track Frequency");
        TString pulllow_title( ";" + trackParam2 + ";Rel Track Frequency");
        TString pullhigh_title(";" + trackParam2 + ";Rel Track Frequency");

        TString Tpulllow( trackParam + "lowpull_" + trackParam2);
        TString Tpullhigh(trackParam + "highpull_" + trackParam2);
        TString Tpull(    trackParam + "fullpull_" + trackParam2);
    
        TString Tpulltype[]  = {Tpull,Tpulllow,Tpullhigh};
        TString Tpulltitle[] = {pull_title,pulllow_title,pullhigh_title};

        TP_pull[itp][itp2][ipulltype] = new TH1F(Tpulltype[ipulltype],Tpulltitle[ipulltype],30,TP_ranges_min[itp2],TP_ranges_max[itp2]);
      }
    }

    for (int iibl = 0; iibl < nibl;iibl++){
      string iblname = iblnames[iibl];
      for( unsigned int ieta = 0; ieta < etabins.size();ieta++){
	       for( unsigned int iipt = 0; iipt < invptbinsvec.size();iipt++){
	          TString hist_res_name("hist_res" + trackParam + "_" + iblname + "_eta" + to_string(ieta) + "_ipt" + to_string(iipt));
		        hist_res[itp][iibl][ieta][iipt] = new TH1F( hist_res_name,hist_res_title,nbins, -range,range);

	       }
      }
    }
  }
  ientry2=0;
}
////////////////////////////////// <INITIALISE>

////////////////////////////////// <INITIALISE 2 >
void Init2() {
  nbins = 61 ;
  for(int itp = 0; itp < nParams;itp++){
    string trackParam = trackParam_list[itp];
    TString TtrackParam(trackParam);
    TString TtrackParam_title(";" + trackParam + "(reco-truth)/#sigma;N Tracks");

    for (int iibl = 0; iibl < nibl;iibl++){
      string iblname = iblnames[iibl];
      for( unsigned int ieta = 0; ieta < etabins.size();ieta++){
	       for( unsigned int iipt = 0; iipt < invptbinsvec.size();iipt++){
	          TString hist_res_name("hist_rebinned_res" + trackParam + "_" + iblname + "_eta" + to_string(ieta) + "_ipt" + to_string(iipt));
	          TString hist_res_title(";#Delta" + trackParam + "[ftk-truth](mm);N Tracks");
            hist_std[itp][iibl][ieta][iipt] = hist_res[itp][iibl][ieta][iipt]->GetStdDev();
            delete hist_res[itp][iibl][ieta][iipt];
		        hist_res[itp][iibl][ieta][iipt] = new TH1F( hist_res_name,hist_res_title,nbins, -4.5*hist_std[itp][iibl][ieta][iipt],4.5*hist_std[itp][iibl][ieta][iipt]);

	       }
      }
    }
  }

  ientry2=0;
}
////////////////////////////////// </INITIALISE 2 >

////////////////////////////////// <PROCESS>
void Process(Long64_t ientry) {

  //////////////////////////////////////////////////////////////////////////////////////////////
  ///// FTK TRACK DETERMINATION
  t_ftkdata->GetEntry(ientry);

  Int_t evtnumber_ftk = ftktracks->eventNumber();
  Int_t runnumber_ftk = ftktracks->runNumber();

  for (; ientry2 < t_truth->GetEntries(); ientry2++) {
    t_evtinfo->GetEntry(ientry2);
    if (EventNumber==evtnumber_ftk && RunNumber==runnumber_ftk) break;
  }

  if (ientry2 == t_truth->GetEntries()) {
    for (ientry2 = 0; ientry2 < t_truth->GetEntries(); ientry2++) {
      t_evtinfo->GetEntry(ientry2);
    if (EventNumber==evtnumber_ftk && RunNumber==runnumber_ftk) break;
    }
  }
  // if we restarted from the beginning then return an error. Bad luck. 
  if (ientry2 == t_truth->GetEntries()) {
    cerr << Form("Mismatch between entries: (%d,%d)",runnumber_ftk,evtnumber_ftk) << endl;
    ientry2=0;
    return;
  }
  t_truth->GetEntry(ientry2);

  if (ientry%10000==0) { Info("Process","Event %lld, (Run,Evt) = (%d,%d)",ientry, runnumber_ftk, evtnumber_ftk);}

  // tracks is defined in main. Retreive the number of tracks 
  Int_t ntracks = ftktracks->getNTracks();

  FTKBarcodeMM ftkmatchinfo;
  const FTKTrack *curtrack;
 
  for (Int_t itrk=0;itrk!=ntracks;++itrk) { // loop over the FTK tracks
    curtrack =  ftktracks->getTrack(itrk);
    if (curtrack->getEventIndex()==0) {
      // Information on FTK tracks relative to the hard-scattering events
      // are collected in a vector and later used to build matching maps
      if (curtrack->getBarcodeFrac()>.5) {
          ftkmatchinfo.insert(pair<MatchInfo, const FTKTrack*>(MatchInfo(curtrack->getBarcode(),curtrack->getEventIndex()),curtrack));
      }
    } //end eventIndex == 0
  }//loop over ftk tracks 
  
  /////
  //////////////////////////////////////////////////////////////////////////////////////////////
  
  //////////////////////////////////////////////////////////////////////////////////////////////
  ///// LOOP THROUGH TRUTH TRACKS 
  Int_t ntruth = truthtracks->size();

  vector<FTKTruthTrack>::const_iterator itr = truthtracks->begin(); 
  vector<FTKTruthTrack>::const_iterator itrE = truthtracks->end();

  double nthftktrack = 0;
  //double ntruthtracks= 1;
  double ntracks_passed = 0;
  //std::cout << "NEWLOOP" << std::endl;
  for (;itr!=itrE;++itr) { // loop over the truth tracks
    const FTKTruthTrack &curtruth = (*itr);
    

    if (curtruth.getEventIndex()!=0 && curtruth.getQ()==0) continue;
    int    barcode   = curtruth.getBarcode();
    if (barcode>100000 || barcode==0) continue;
    double px = curtruth.getPX(); 
    double py = curtruth.getPY(); 
    double pt = TMath::Sqrt(px*px+py*py);
    if ( pt < ptmincut ) continue;
    double invpt     = 1./(2.0*pt);
    double d0        = curtruth.getD0();
    if (d0<d0min || d0>d0max) continue;
    double z0        = curtruth.getZ();
    if (z0<z0min || z0>z0max) continue;
    double curv      = curtruth.getQ()*invpt;
    if (curv < -abscurvmax || curv>abscurvmax) continue;
    double phi       = curtruth.getPhi();
    if (phi<phimin || phi>phimax) continue;
    double eta       = curtruth.getEta();
    if (eta<etamin || eta>etamax) continue;
    double qOverp    = curv;
    ntracks_passed = ntracks_passed + 1;
 

    MatchInfo  reftruth(barcode,curtruth.getEventIndex());

    pair<FTKBarcodeMM::const_iterator,FTKBarcodeMM::const_iterator> mrange = ftkmatchinfo.equal_range(reftruth);

    if (mrange.first != mrange.second) {
      const FTKTrack *bestftk(0x0);

      for(FTKBarcodeMM::const_iterator ftkI = mrange.first;ftkI!=mrange.second;++ftkI) {
        // for every iterator where the barcode and event index of the first are not equal to the barcode and event index of the second. 
        if (!bestftk){
          bestftk = (*ftkI).second; //ftktrack
        } else if (bestftk->getBarcodeFrac()<(*ftkI).second->getBarcodeFrac()) {
          bestftk = (*ftkI).second; // if the barcodefraction is less than the new one. Replace.
        }
      }
     
      if (bestftk) {   
        //std::cout <<"BARCODES::FTK,TRUTH" << "(" << bestftk->getBarcode() << "," << barcode << ")" << std::endl;  
        //std::cout <<"BARCODES::FRAC::" << bestftk->getBarcodeFrac() << std::endl;  
        //if (bestftk->getEventIndex()!=0) continue;
        //int    barcodeftk   = bestftk->getBarcode();
        //if (barcodeftk>100000 || barcodeftk==0) continue;
        double ptftk = bestftk->getPt();
        //if ( ptftk < ptmincut )continue;
        double invpt_ftk     = 1./(2.0*ptftk);
        //double d0ftk        = bestftk->getIP();
        //if (d0ftk<d0min || d0ftk>d0max) continue;
        //double z0ftk        = bestftk->getZ0();
        //if (z0ftk<z0min || z0ftk>z0max) continue;
        //double curvftk      = invpt_ftk;
        //if (curvftk < -abscurvmax || curvftk>abscurvmax) continue;
        //double phiftk       = bestftk->getPhi();
        //if (phiftk<phimin || phiftk>phimax) continue;
        //double etaftk       = bestftk->getEta();
        //if (etaftk<etamin || etaftk>etamax) continue;
        
        nthftktrack = nthftktrack + 1;
	      FTKHit hit  = bestftk->getFTKHit(0);
	      bool isIBL  = hit.getPlane() == 0;
        //isIBL =0;
        double TP_truth[] = {d0,z0,eta,phi,qOverp};
        double TP_ftk[]   = {bestftk->getIP(),bestftk->getZ0(),bestftk->getEta(),bestftk->getPhi(),invpt_ftk};
	      for(unsigned int ieta = 1; ieta < etabins.size(); ieta++){
	         etabin_min = etabins.at(ieta-1);
	         etabin_max = etabins.at(ieta);
           //std::cout << "(etamin,etamax)::" << "(" << etabin_min << ", " << etabin_max << ")" << std::endl;
           if ( TMath::Abs(eta) < etabin_min || TMath::Abs(eta) > etabin_max) continue;

	         for(unsigned int iipt = 1; iipt < invptbinsvec.size();iipt++){
	           invptbin_min = invptbinsvec.at(iipt-1);
	           invptbin_max = invptbinsvec.at(iipt);
  		       if ( qOverp < invptbin_min || qOverp > invptbin_max        ) continue;

             for(int itp = 0; itp < 5; itp++){
                hist_res[itp][isIBL][ieta-1][iipt-1]->Fill(TP_ftk[itp] - TP_truth[itp]);
             }
        	 } // invpt bins
	       } // eta bins
      } //best ftk loop
    }// matching

  } // end loop over truth tracks
  if(0){
    std::cout << "=========================================================" <<  std::endl;
    std::cout << "=====================END OF TRUTH LOOP ==================" <<  std::endl;
      std::cout << "pre quality control" << std::endl;
      std::cout << "ftk track: " << ntracks << " Truth Tracks: " << ntruth << std::endl;
      std::cout << "post quality control: " << std::endl;
      std::cout << "ftk tracks: " << nthftktrack << " ntracks_passed::" << ntracks_passed << std::endl;
      if(nthftktrack > ntracks_passed){
        std::cout << "WARNING::MORE FTK TRACKS" << std::endl;
      }
      std::cout << "=========================================================" <<  std::endl;
  }

}
////////////////////////////////// </PROCESS>

////////////////////////////////// <WIDTH CALCULATION>
void width_calculation(){
  for(int itp = 0;itp < nParams ; itp++){
    string trackParam = trackParam_list[itp];
    //double trackRange = trackParam_range[itp];

    for(int iibl =0;iibl < nibl;iibl++){
      string iblname  = iblnames[iibl];

      for(unsigned int ieta =0;ieta < etabins.size()-1;ieta++){
        //widths.clear();
        //width_errors.clear();

        for(unsigned int iipt=0;iipt < invptbinsvec.size()-1;iipt++){
	         Double_t par[6]; /* parameter array */

          double  histrange = 4.5*hist_std[itp][iibl][ieta][iipt];
          double  corerange = hist_std[itp][iibl][ieta][iipt];
	         TF1 *g1 = new TF1 ("m1","gaus",-corerange,corerange);
	         TF1 *g2 = new TF1 ("m2","gaus",-histrange,histrange);
	         //TF1 *f1 = new TF1("double_gaus","gaus(0)+gaus(3)",-histrange,histrange);
           TF1 *dg = new TF1("dgtest",doublegaussfunc,-histrange,histrange,5);

          TF1 *singlegaussian = new TF1 ("singlegaussian","gaus",-histrange,histrange);

          hist_res[itp][iibl][ieta][iipt]->Fit(g1,"QR0");
          hist_res[itp][iibl][ieta][iipt]->Fit(singlegaussian,"QR0");

           double width       = singlegaussian->GetParameter(2);
           double width_error = singlegaussian->GetParError(2);

          par[0] = g1->GetParameter(0); //single gaussian height
          par[1] = g1->GetParameter(1); //sg mean 
          par[2] = g1->GetParameter(2); //sg sigma

          par[3] = par[1]    ;//mean
          par[4] = 3.0*par[2]; //sigma
          dg->SetParameters(par);
          dg->SetParLimits(0,0,7000);
          dg->FixParameter(1,0.0);
          dg->FixParameter(3,0.0);
          dg->SetParLimits(4,par[2],10*par[2]);
          hist_res[itp][iibl][ieta][iipt]->Fit(dg,"QR0");
          TF1 *coretest = new TF1 ("coretest","gaus",-histrange,histrange);
          TF1 *tailtest = new TF1 ("tailtest","gaus",-histrange,histrange);
          coretest->SetParameters(dg->GetParameter(0),dg->GetParameter(1),dg->GetParameter(2) );
          tailtest->SetParameters(0.11*dg->GetParameter(0),dg->GetParameter(3),dg->GetParameter(4));

          double corewidth  = dg->GetParameter(2); 
          double tailwidth = dg->GetParameter(4);


          double corewidth_error = dg->GetParError(2); 
          double tailwidth_error = dg->GetParError(4); 


           double n1           = dg->GetParameter(0);// double error_n1     = dg->GetParError(0);
         //  double mu1          = dg->GetParameter(1); double error_mu1    = dg->GetParError(1);
         //  double sigma1       = dg->GetParameter(2); double error_sigma1 = dg->GetParError(2);
         //  double var1         = sigma1*sigma1;

         //  double yield1 = coretest->Integral(-histrange,histrange); 
         //  double yield2 = tailtest->Integral(-histrange,histrange);
           /* TAILS GAUSSIAN */
           double n2           = 0.11*dg->GetParameter(0);//  double error_n2     = dg->GetParameter(0);
         //  double mu2          = dg->GetParameter(3);  double error_mu2    = dg->GetParError(3);
         //  double sigma2       = dg->GetParameter(4);  double error_sigma2 = dg->GetParError(4);
         //  double var2         = sigma2*sigma2;

          //EDIT CHANGE TO INTEGRALS.
         // double weight1 = n1/(n1+n2);
         // double weight2 = n2/(n1+n2);      


          //double width = sqrt(weight1*var1 +weight2*var2); // APPROXIMATION
          /*CALCULATE THE ERROR OF THE VARIANCE */
         // double a_n1     = 1/(n1+n2);
	       // double a_n2     = -pow(n1+n2,-2);
  	     // double sigma_w1 =  sqrt( a_n1*a_n1*error_n1*error_n1 + a_n2*a_n2*error_n2*error_n2);
 
        	//double b_n1 = -pow(n1+n2,-2);
	        //double b_n2 = 1/(n1+n2);
	        //double sigma_w2 =  sqrt( b_n1*b_n1*error_n1*error_n1 + b_n2*b_n2*error_n2*error_n2);


        //  std::cout << "n1::Integral: " << core->Integral(-histrange,histrange)  << std::endl;
         // std::cout << "n2::Integral: " << tails->Integral(-histrange,histrange) << std::endl;
         // std::cout << "w1::Integral: " << core->Integral(-histrange,histrange)/(core->Integral(-histrange,histrange) + tails->Integral(-histrange,histrange) ) << std::endl;
         // std::cout << "w1::Yield   : " << n1/(n1+n2)                            << std::endl;

          

           /* SINGLE GAUSSIAN */
           //double width = g1->GetParameter(2);
           //double width_error = g1->GetParError(2);

          /* DOUBLE GAUSSIAN */
           /* ===================================================================================================================== */
          /* WANT TO RETREIVE THE WIDTH AND THE ERROR IN THAT WIDTH. */
          /* WIDTH FOR A SUM OF RANDOMLY DISTRIBUTED VARIABLES IS A WELL DEFINED CONCEPT. */
          /* FOR TWO GAUSSIANS WE HAVE: */
          //double mubar = (n1*mu1+n2*mu2)/(n1+n2);
          // sqrt(((var1 +mu1*mu1)*n1/(n1+n2) + (var2+mu2*mu2)*n2/(n1+n2) - mubar)
          //double width = sqrt(((var1 +mu1*mu1)*n1 + (var2+mu2*mu2)*n2 )/(n1+n2) - mubar)
          /* WHICH FOR VERY SMALL MEANS EQUALS */
          /* double width = sqrt((var1*n1 + var2*n2 )/(n1+n2) ); */
          /* double width = sqrt( weight1*variance1 + weight2*variance2); */
          /* ===================================================================================================================== */

          /* WEIGHTS OF EACH GAUSSIAN FIRST */
          //n1 = f1->Integral();
          //n2 = 

          /* error^2 = S(1)^2*E(W1)^2 + S(2)^2*E(W2)^2 + W1^2*E(W1)^2 + W2^2*E(W2)^2 */
         // double width_error = sqrt(var1*sigma_w1*sigma_w1+ var2*sigma_w2*sigma_w2 +weight1*weight1*sigma_w1*sigma_w1 + weight2*weight2*sigma_w2*sigma_w2);

          //widths.push_back(width);
          //width_errors.push_back(width_error);

          //corewidths.push_back(corewidth);
          //tailwidths.push_back(tailwidth);

          //corewidth_errors.push_back(corewidth_error);
          //tailwidth_errors.push_back(tailwidth_error);

          //ratios.push_back((n2/n1)*100);
          A_ratio[itp][iibl][ieta][iipt] = (n2/n1)*100;

          A_corewidth[itp][iibl][ieta][iipt] = corewidth;
          A_tailwidth[itp][iibl][ieta][iipt] = tailwidth;

          A_corewidth_errors[itp][iibl][ieta][iipt] = corewidth_error;
          A_tailwidth_errors[itp][iibl][ieta][iipt] = tailwidth_error;

          A_width[itp][iibl][ieta][iipt]        = width;
          A_width_errors[itp][iibl][ieta][iipt] = width_error;

          delete coretest;
          delete tailtest;
          delete singlegaussian;
          delete g1;
          delete g2;
          delete dg;
        } //invptbins
      } //eta bins
    } // ibl
  } // parameters 
} 
////////////////////////////////// </WIDTH CALCULATION>

////////////////////////////////// <SQUARE ROOT FITS>
void sqrt_fit(){
  for(int itp = 0;itp < nParams ; itp++){
    string trackParam = trackParam_list[itp];
    //double trackRange = trackParam_range[itp];
    for(int iibl =0;iibl < nibl;iibl++){
      string iblname = iblnames[iibl];
      for(unsigned int ieta =0;ieta < etabins.size()-1;ieta++){
          double width_arr[200];
          double width_errors_arr[200]; 

          for (unsigned int i = 0; i < invptbinsvec.size()-1; i++) { 
            width_arr[i]        = A_width[itp][iibl][ieta][i];
            width_errors_arr[i] = A_width_errors[itp][iibl][ieta][i];
          }
     
          TGraphErrors *graph_linear = new TGraphErrors(ninvptbins,ipt_arr,width_arr,0,width_errors_arr);
          TF1 *linfit = new TF1("linearfit","pol1",0.0,invpt_max);
         
          //std::cout << "widtharr[end]:: " << width_arr[ninvptbins-1] << std::endl;
          //std::cout << "widtharr[ninvptbins]:: " << width_arr[ninvptbins] << std::endl;
          //std::cout << "widtharr[mid]:: " << width_arr[middlebin1]   << std::endl;
          
          linfit->SetParameters((width_arr[middlebin1] + width_arr[middlebin2])/2.0 , (width_arr[middlebin1]-width_arr[ninvptbins-1])/(invpt_max)  );
          
          //linfit->SetParLimits(0,0.0 , 1.05*(width_arr[middlebin1] + width_arr[middlebin2])/2.0);
          
          //std::cout << "0.95*(width_arr[middlebin1] + width_arr[middlebin2])/2.0::" << 0.95*(width_arr[middlebin1] + width_arr[middlebin2])/2.0 << std::endl;
          //std::cout << "1.05*(width_arr[middlebin1] + width_arr[middlebin2])/2.0::" << 1.05*(width_arr[middlebin1] + width_arr[middlebin2])/2.0 << std::endl;
          graph_linear->Fit(linfit,"QR");
          
          double linpar0 = linfit->GetParameter(0); 
          double linpar1 = linfit->GetParameter(1);
          
          TGraphErrors *graph_sqrt = new TGraphErrors(ninvptbins,ipt_arr,width_arr,0,width_errors_arr);
          
          TF1 *func = new TF1("sqrtfit",fitf,-invpt_max,invpt_max,2);
          
          double a = linpar0*linpar0;
          double rangemax = invpt_max;
          double linfuncmax = linpar1*rangemax + linpar0;
          double b = (linfuncmax*linfuncmax - a*a)/rangemax/rangemax;
          
          func->SetParameters(a,b);
          func->SetParLimits(0,0.95*linpar0*linpar0, 1.05*linpar0*linpar0);
          
          //func->SetParNames ("a","b");
          graph_sqrt->Fit(func,"QR");
          
          A_sq_par0[itp][iibl][ieta] = func->GetParameter(0);
          A_sq_par1[itp][iibl][ieta] = func->GetParameter(1);
          
          delete func;
          delete graph_linear;
          delete graph_sqrt;
      } //eta bins
    } // ibl
  } // parameters 
}

void coresqrt_fit(){
  for(int itp = 0;itp < nParams ; itp++){
    string trackParam = trackParam_list[itp];
    //double trackRange = trackParam_range[itp];
    for(int iibl =0;iibl < nibl;iibl++){
      string iblname = iblnames[iibl];
      for(unsigned int ieta =0;ieta < etabins.size()-1;ieta++){
          double width_arr[200];
          double width_errors_arr[200]; 

          for (unsigned int i = 0; i < invptbinsvec.size()-1; i++) { 
            width_arr[i]        = A_corewidth[itp][iibl][ieta][i];
            width_errors_arr[i] = A_corewidth_errors[itp][iibl][ieta][i];
          }
     
          TGraphErrors *graph_linear = new TGraphErrors(ninvptbins,ipt_arr,width_arr,0,width_errors_arr);
          TF1 *linfit = new TF1("linearfit","pol1",0.0,invpt_max);
          
          linfit->SetParameters((width_arr[middlebin1] + width_arr[middlebin2])/2.0 , (width_arr[middlebin1]-width_arr[ninvptbins-1])/(invpt_max*invpt_max)  );

          //linfit->SetParLimits(0,0.95*(width_arr[middlebin1] + width_arr[middlebin2])/2.0 , 1.05*(width_arr[middlebin1] + width_arr[middlebin2])/2.0);

          graph_linear->Fit(linfit,"QR");

          double linpar0 = linfit->GetParameter(0); 
          double linpar1 = linfit->GetParameter(1);

          TGraphErrors *graph_sqrt = new TGraphErrors(ninvptbins,ipt_arr,width_arr,0,width_errors_arr);

          TF1 *func = new TF1("sqrtfit",fitf,-invpt_max,invpt_max,2);

          double a = linpar0*linpar0;
          double rangemax = invpt_max;
          double linfuncmax = linpar1*rangemax + linpar0;
          double b = (linfuncmax*linfuncmax - a*a)/rangemax/rangemax;

          func->SetParameters(a,b);
          func->SetParLimits(0,0.95*linpar0*linpar0, 1.05*linpar0*linpar0);
          //func->SetParNames ("a","b");
          graph_sqrt->Fit(func,"QR");

          A_coresq_par0[itp][iibl][ieta] = func->GetParameter(0);
          A_coresq_par1[itp][iibl][ieta] = func->GetParameter(1);
          delete func;
          delete graph_linear;
          delete graph_sqrt;
      } //eta bins
    } // ibl
  } // parameters 
}

void tailsqrt_fit(){
  for(int itp = 0;itp < nParams ; itp++){
    string trackParam = trackParam_list[itp];
    //double trackRange = trackParam_range[itp];
    for(int iibl =0;iibl < nibl;iibl++){
      string iblname = iblnames[iibl];
      for(unsigned int ieta =0;ieta < etabins.size()-1;ieta++){
          double width_arr[200];
          double width_errors_arr[200]; 

          for (unsigned int i = 0; i < invptbinsvec.size()-1; i++) { 
            width_arr[i]        = A_tailwidth[itp][iibl][ieta][i];
            width_errors_arr[i] = A_tailwidth_errors[itp][iibl][ieta][i];
          }
     
          TGraphErrors *graph_linear = new TGraphErrors(ninvptbins,ipt_arr,width_arr,0,width_errors_arr);
          TF1 *linfit = new TF1("linearfit","pol1",0.0,invpt_max);
          
          linfit->SetParameters((width_arr[middlebin1] + width_arr[middlebin2])/2.0 , (width_arr[middlebin1]-width_arr[ninvptbins-1])/(invpt_max*invpt_max)  );

         // linfit->SetParLimits(0,0.95*(width_arr[middlebin1] + width_arr[middlebin2])/2.0 , 1.05*(width_arr[middlebin1] + width_arr[middlebin2])/2.0);

          graph_linear->Fit(linfit,"QR");

          double linpar0 = linfit->GetParameter(0); 
          double linpar1 = linfit->GetParameter(1);

          TGraphErrors *graph_sqrt = new TGraphErrors(ninvptbins,ipt_arr,width_arr,0,width_errors_arr);

          TF1 *func = new TF1("sqrtfit",fitf,-invpt_max,invpt_max,2);

          double a = linpar0*linpar0;
          double rangemax = invpt_max;
          double linfuncmax = linpar1*rangemax + linpar0;
          double b = (linfuncmax*linfuncmax - a*a)/rangemax/rangemax;

          func->SetParameters(a,b);
          func->SetParLimits(0,0.95*linpar0*linpar0, 1.05*linpar0*linpar0);
          //func->SetParNames ("a","b");
          graph_sqrt->Fit(func,"QR");

          A_tailsq_par0[itp][iibl][ieta] = func->GetParameter(0);
          A_tailsq_par1[itp][iibl][ieta] = func->GetParameter(1);
          delete func;
          delete graph_linear;
          delete graph_sqrt;
      } //eta bins
    } // ibl
  } // parameters 
}



/////////////////////////////////  <PULL>
void Pull (Long64_t ientry) {
  //////////////////////////////////////////////////////////////////////////////////////////////
  ///// SET UP FTK TRACKS 

  t_ftkdata->GetEntry(ientry);

  Int_t evtnumber_ftk = ftktracks->eventNumber();
  Int_t runnumber_ftk = ftktracks->runNumber();

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

  if (ientry%10000==0) { Info("Pull","Event %lld, (Run,Evt) = (%d,%d)",ientry, runnumber_ftk, evtnumber_ftk); }

  Int_t ntracks = ftktracks->getNTracks();

  FTKBarcodeMM ftkmatchinfo;
  const FTKTrack *curtrack;

  for (Int_t itrk=0;itrk!=ntracks;++itrk) { // loop over the FTK tracks
    curtrack = ftktracks->getTrack(itrk); 
    if (curtrack->getEventIndex()==0) {
            // Information on FTK tracks relative to the hard-scattering events
            // are collected in a vector and later used to build matching maps
      if (curtrack->getBarcodeFrac()>.5) {
        ftkmatchinfo.insert(pair<MatchInfo, const FTKTrack*>(MatchInfo(curtrack->getBarcode(),curtrack->getEventIndex()),curtrack));
        // 
      }
    }
  }
 
   //////////////////////////////////////////////////////////////////////////////////////////////
   ///// LOOP THROUGH TRUTH TRACKS

  double ntracks_passed = 0;
  vector<FTKTruthTrack>::const_iterator itr = truthtracks->begin();
  vector<FTKTruthTrack>::const_iterator itrE = truthtracks->end();
  for (;itr!=itrE;++itr) { 
    const FTKTruthTrack &curtruth = (*itr);

    if (curtruth.getEventIndex()!=0 && curtruth.getQ()==0) continue;
    int barcode   = curtruth.getBarcode();
    if (barcode>100000   || barcode==0)      continue;
    double px     = curtruth.getPX();
    double py     = curtruth.getPY();
    double pt     = TMath::Sqrt(px*px+py*py);
    if ( pt  < ptmincut )       continue;
    double invpt  = 1./(2.0*pt);
    double d0     = curtruth.getD0();
    if (d0<d0min || d0>d0max)        continue;
    double z0     = curtruth.getZ();
    if (z0<z0min || z0>z0max)        continue;
    double curv   = curtruth.getQ()*invpt;
    if (curv<-abscurvmax || curv>abscurvmax) continue;
    double phi    = curtruth.getPhi();
    
    if (phi<phimin       || phi>phimax)      continue;
    double eta    = curtruth.getEta();
    if (eta<etamin       || eta>etamax)      continue;
    double qOverp = curv;

    ntracks_passed =  ntracks_passed + 1;

    // match the barcode and event index values
    MatchInfo reftruth(barcode,curtruth.getEventIndex());

    pair<FTKBarcodeMM::const_iterator,FTKBarcodeMM::const_iterator> mrange = ftkmatchinfo.equal_range(reftruth);

    if (mrange.first != mrange.second) {
      const FTKTrack *bestftk(0x0);
      for(FTKBarcodeMM::const_iterator ftkI = mrange.first;ftkI!=mrange.second;++ftkI) {
        if (!bestftk) {
          bestftk = (*ftkI).second;
        } else if (bestftk->getBarcodeFrac()<(*ftkI).second->getBarcodeFrac()) {
          bestftk = (*ftkI).second;
        }
      }



        if (bestftk) {
          

             double ptftk = bestftk->getPt();
             double invpt_ftk     = 1./(2.0*ptftk);

          FTKHit hit = bestftk->getFTKHit(0);
	        bool isIBL = hit.getPlane() == 0;
          //isIBL =0;
          double TP_truth[] = {d0,z0,eta,phi,qOverp};
          double TP_ftk[]   = {bestftk->getIP(),bestftk->getZ0(),bestftk->getEta(),bestftk->getPhi(),invpt_ftk}; 
          TF1 *coregaussian = new TF1 ("coregaussian","gaus",-15,15.);
          TF1 *tailgaussian = new TF1 ("tailgaussian","gaus",-15.,15.);     

               double widthcalc;  
               double tailwidthcalc;
               double corewidthcalc;    
               double corrected_width;  
         for(unsigned int ieta = 1; ieta < etabins.size(); ieta++){
             etabin_min = etabins.at(ieta-1);
             etabin_max = etabins.at(ieta);
            // std::cout << "(min,max)::" << etabin_min << "," << etabin_max << std::endl;
             if ( TMath::Abs(eta) < etabin_min || TMath::Abs(eta) > etabin_max ) continue;



                //FIX FIX FIX FIX
               for(int itp = 0; itp < 5; itp++){
              
                  widthcalc     = TMath::Sqrt(A_sq_par0[itp][isIBL][ieta-1]     + A_sq_par1[itp][isIBL][ieta-1]*qOverp*qOverp ); 
                  corewidthcalc = TMath::Sqrt(A_coresq_par0[itp][isIBL][ieta-1] + A_coresq_par1[itp][isIBL][ieta-1]*qOverp*qOverp );
                  tailwidthcalc = TMath::Sqrt(A_tailsq_par0[itp][isIBL][ieta-1] + A_tailsq_par1[itp][isIBL][ieta-1]*qOverp*qOverp );

                  coregaussian->SetParameters(1.0,0.0,corewidthcalc );
                  tailgaussian->SetParameters(1.0,0.0,tailwidthcalc);

                  double yield1 = coregaussian->Integral(-15.,15.);
                  double yield2 = tailgaussian->Integral(-15.,15.);

                  double weight1 = yield1 / (yield1 + 0.11*yield2); //fixed ratio 
                  double weight2 = 0.11*yield2 / (yield1 + 0.11*yield2); // fixed ratio

                  corrected_width = TMath::Sqrt(  weight1*corewidthcalc*corewidthcalc + weight2*tailwidthcalc*tailwidthcalc);

                  double pullTP    = (TP_ftk[itp] - TP_truth[itp])/corrected_width;
                  double pullTP_sg = (TP_ftk[itp] - TP_truth[itp])/widthcalc;

                  for( int itp2 = 0;itp2 < 5; itp2++){
                     TP_pull[itp][itp2][0]->Fill(TP_truth[itp2]);
                     if (pullTP < 5.0){
                        TP_pull[itp][itp2][1]->Fill(TP_truth[itp2]);
                     }
                     else if (pullTP > 5.0) {
                       TP_pull[itp][itp2][2]->Fill(TP_truth[itp2]);
                     }
                  }
                  //pull_IBL[itp][isIBL]->Fill(pullTP);

                  pull_res[itp]->Fill(pullTP);
                  pull_res_sg[itp]->Fill(pullTP_sg);
                  
              }
              delete coregaussian;
              delete tailgaussian;
             // } // invpt bins
            } // eta bins 
	 //	  }// ibl bin

      } //best ftk loop
    }// matching
   // continue;

  } // end loop over truth tracks

} // processing 
///////////////////////////////// </PULL>


void Terminate(std::string& outputname,std::string& outputfolder) {


  Info("Terminate","Adding the histograms to the file: %s", outputname.c_str());
  TFile *ofile = TFile::Open(outputname.c_str(), "recreate");
 
  //double ipt_arr[200];
  //std::copy(invptbinsvec.begin(), invptbinsvec.end(), ipt_arr);
  //for( unsigned int i = 0; i < invptbinsvec.size();i++){ipt_arr[i] = ipt_arr[i] + 0.5*stepsize;}
 
  ofstream sqrtIBL("sqrtIBL.log");
  ofstream sqrtnoIBL("sqrtnoIBL.log");
  ofstream linearIBL("linearIBL.log");
  ofstream linearnoIBL("linearnoIBL.log");


  //Info("Plots","Producing Histograms");
  for(int itp = 0;itp < nParams ; itp++){
    string trackParam = trackParam_list[itp];
    // std::cout << "pull4" << std::endl;
    //double trackRange = trackParam_range[itp];
    string trackParam_unit = trackParam_units[itp];
 
    for(int iibl =0;iibl < nibl;iibl++){
      string iblname = iblnames[iibl];

        TLatex* atlas_title = new TLatex();
  
        atlas_title->SetTextAlign(11);
        atlas_title->SetTextSize(0.04);
        atlas_title->SetTextFont(42);
        atlas_title->SetNDC();
        std::string tempname ="#bf{#it{ATLAS}} Internal " + iblname ;
        const char *atlas_title_name = tempname.c_str();

      for(unsigned int ieta =0;ieta < etabins.size()-1;ieta++){
    	  widths.clear();
	      width_errors.clear();
        etabin_min = etabins.at(ieta);
        etabin_max = etabins.at(ieta+1);

        TLatex* eta_latex = new TLatex();
        eta_latex->SetTextAlign(11);
        eta_latex->SetTextSize(0.04);
        eta_latex->SetTextFont(42);
        eta_latex->SetNDC();

        std::string etatempname_linear ="#eta in [" + to_string(round(etabin_min,2)).substr(0,3) + "," + to_string(round(etabin_max,2)).substr(0,3) + "]" ;
      
        const char *eta_label = etatempname_linear.c_str();


     // Info("Plots","Eta bin: %d", ieta);
    	for(unsigned int iipt=0;iipt < invptbinsvec.size()-1;iipt++){
        
        //////////////////////////////////////////////////////////////////////////////////////
        ///// SET UP FOR ALL PLOTS 
        invptbin_min = invptbinsvec.at(iipt);
        invptbin_max = invptbinsvec.at(iipt+1);

        double nprecision = 4;
        if(invptbin_min < 0.0  || invptbin_max < 0.0){
          nprecision = 5;
        }  
        std::string iptmin_string = to_string(round(invptbin_min*1e3,5)).substr(0,nprecision);
        std::string iptmax_string = to_string(round(invptbin_max*1e3,5)).substr(0,nprecision);

        TLatex* etapt = new TLatex();
        etapt->SetTextAlign(11);
        etapt->SetTextSize(0.04);
        etapt->SetTextFont(42);
        etapt->SetNDC();
      
        std::string etapttempname ="#eta in [" + to_string(round(etabin_min,2)).substr(0,3) + "," + to_string(round(etabin_max,2)).substr(0,3) + "] and 1/2p_{T} in [" + iptmin_string + "," + iptmax_string + "]";
        
        const char *etapt_name = etapttempname.c_str();



      //////////  END OF COMMON USE DEFINITIONS
      /////////////////////////////////////////////////////////////////////////////





        ///////////////////////////////////////////////////////////////////
        //////// HISTOGRAMS 

        ///// CANVAS CREATION
        TCanvas *c = new TCanvas("c", "c", 800, 650);
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
        pad1->SetLogy();
        pad1->SetBottomMargin(0.05); // Upper and lower plot are joined
        pad1->SetGridx();         // Vertical grid
        pad1->Draw();             // Draw the upper pad: pad1
        pad1->cd();               // pad1 becomes the current pad

	      Double_t par[6]; /* parameter array */


        ////// FIT GENERATION 
        double  histrange = 4.5*hist_std[itp][iibl][ieta][iipt];
        double  corerange = hist_std[itp][iibl][ieta][iipt];

	      TF1 *g1 = new TF1 ("m1","gaus",-corerange,corerange);
        TF1 *dg = new TF1("dgtest",doublegaussfunc,-histrange,histrange,5);

        hist_res[itp][iibl][ieta][iipt]->Fit(g1,"QR0");

        par[0] = g1->GetParameter(0); //single gaussian height
        par[1] = g1->GetParameter(1); //sg mean 
        par[2] = g1->GetParameter(2); //sg sigma

        par[3] = par[1]    ;//mean
        par[4] = 3.0*par[2]; //sigma
        dg->SetParameters(par);
        dg->SetParLimits(0,0,7000);
        dg->FixParameter(1,0.0);
        dg->FixParameter(3,0.0);
        dg->SetParLimits(4,par[2],10*par[2]);
        hist_res[itp][iibl][ieta][iipt]->Fit(dg,"QR0");
        TF1 *coretest = new TF1 ("coretest","gaus",-histrange,histrange);
        TF1 *tailtest = new TF1 ("tailtest","gaus",-histrange,histrange);
        coretest->SetParameters(dg->GetParameter(0),dg->GetParameter(1),dg->GetParameter(2) );
        tailtest->SetParameters(0.11*dg->GetParameter(0),dg->GetParameter(3),dg->GetParameter(4));

        dg->SetLineColor(kRed+1);
        coretest->SetLineColor(kBlue+1);
        tailtest->SetLineColor(kGreen+2);

        dg->SetLineStyle(1);
        coretest->SetLineStyle(2);
        tailtest->SetLineStyle(4);

        gStyle->SetOptFit(1111);

        //// DRAW HISTOGRAM
        hist_res[itp][iibl][ieta][iipt]->Draw();
        dg->Draw("SAME");
        coretest->Draw("SAME");
        tailtest->Draw("SAME");

        c->cd();

        atlas_title->DrawLatex(.12, 0.85, atlas_title_name);
        etapt->DrawLatex(.12, 0.95, etapt_name);
        
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.2);
        pad2->SetGridx(); // vertical grid
        pad2->SetGridy();
        pad2->Draw();
        pad2->cd();       // pad2 becomes the current pad
      
        TH1F *h3 = (TH1F*)hist_res[itp][iibl][ieta][iipt]->Clone("h3");
        h3->GetFunction("dgtest")->SetBit(TF1::kNotDraw);
        h3->SetLineColor(kBlack);
        h3->SetMinimum(0);  // Define Y ..
        h3->SetMaximum(3); // .. range
        h3->Sumw2();
        h3->SetStats(0);      // No statistics on lower plot
        h3->Divide(dg);
       
        //         std::cout << "s" << std::endl;
        h3->SetMarkerStyle(7);
        //        std::cout << "t" << std::endl;
        h3->Draw("ep");       // Draw the ratio plot
 
        h3->SetMinimum(0);  // Define Y ..
        h3->SetMaximum(3); // .. range
   
        h3->GetYaxis()->SetTitle("MC/Fit");
        h3->GetYaxis()->SetNdivisions(505);
        h3->GetYaxis()->SetTitleSize(20);
        h3->GetYaxis()->SetTitleFont(43);
        h3->GetYaxis()->SetTitleOffset(1.55);
        h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        h3->GetYaxis()->SetLabelSize(15);

        // X axis ratio plot settings
        h3->GetXaxis()->SetTitleSize(20);
        h3->GetXaxis()->SetTitleFont(43);
        h3->GetXaxis()->SetTitleOffset(4.);
        h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        h3->GetXaxis()->SetLabelSize(15);

    	  c->Update();


	      TImage *img = TImage::Create();

	      img->FromPad(c);

    	  string hist_res_name_string = "hist_res" + trackParam + "_" + iblname + "_eta" + to_string(ieta) + "_ipt" + to_string(iipt)+".png";
	      TString hist_res_name("output/" + hist_res_name_string);

	      img->WriteImage(hist_res_name);
	      ofile->Add(hist_res[itp][iibl][ieta][iipt]);
        delete c;
	    }
    //// END HISTOGRAM PRODUCTION 
    ///////////////////////////////////////////////////////////////////////////////

      double width_arr[200];
	    double width_errors_arr[200];
      double widthcore_arr[200];
	    double widthcore_errors_arr[200];
      double widthtail_arr[200];
	    double widthtail_errors_arr[200];
      double ratio_arr[200];
      double ratio_errors_arr[200];
      for (int i = 0; i < ninvptbins; i++) { 
        width_arr[i]            = A_width[itp][iibl][ieta][i];
        width_errors_arr[i]     = A_width_errors[itp][iibl][ieta][i];
        widthcore_arr[i]        = A_corewidth[itp][iibl][ieta][i];
        widthcore_errors_arr[i] = A_corewidth_errors[itp][iibl][ieta][i];
        widthtail_arr[i]        = A_tailwidth[itp][iibl][ieta][i];
        widthtail_errors_arr[i] = A_tailwidth_errors[itp][iibl][ieta][i];
        ratio_arr[i]            = A_ratio[itp][iibl][ieta][i];
        ratio_errors_arr[i]     = 0.0;  
      }


    //////////////////////////////////////
    /// LINEAR PLOTS 

	    Int_t n = ninvptbins;

	    TCanvas *clinear = new TCanvas("clinear", "clinear", 800, 650);
      clinear->cd();

      TGraphErrors *graph_linear = new TGraphErrors(n,ipt_arr,width_arr,0,width_errors_arr);
      TF1 *linfit = new TF1("linearfit","pol1",0.0,invpt_max);
      graph_linear->Fit(linfit,"QR");

      double linpar0 = linfit->GetParameter(0); 
      double linpar1 = linfit->GetParameter(1);

      std::string linear_name = "linear_" + trackParam + "_" + iblname + "_eta" + std::to_string(ieta);
      const char *graph_linear_name = linear_name.c_str();
      TString linear_title(";Q/2p_{T} (MeV^{-1});#sigma(" + trackParam + ")" + trackParam_unit );
      TString YaxisLinname("#sigma(" + trackParam + ")"+ trackParam_unit );
	    TString XaxisLinname("Q/2p_{T} [MeV]^{-1}");

      graph_linear->SetNameTitle(graph_linear_name,linear_title);
	    graph_linear->SetFillColor(1);
	    graph_linear->SetLineWidth(1);
      graph_linear->SetMarkerStyle(20);
      graph_linear->SetMarkerColor(2);

	    graph_linear->GetYaxis()->SetTitle(YaxisLinname);
	    graph_linear->GetXaxis()->SetTitle(XaxisLinname);
      graph_linear->SetTitle("");
      graph_linear->Draw("AEP");

      atlas_title->DrawLatex(.12, 0.85, atlas_title_name);
      eta_latex->DrawLatex(.12,0.95,eta_label);

      clinear->Update();
      TString linear_save("output/" + linear_name + ".png");
      TImage *imglinear = TImage::Create();
	    imglinear->FromPad(clinear);
	    imglinear->WriteImage(linear_save);
      delete clinear;
      //delete tlinear;
      delete graph_linear;
      delete linfit;


//////////////////////////////////////////////////////////////////
 
        TCanvas *csqroot = new TCanvas("cc", "cc", 800, 650);
        TGraphErrors *graph_sqrt = new TGraphErrors(n,ipt_arr,width_arr,0,width_errors_arr);

        TF1 *func = new TF1("sqrtfit",fitf,-invpt_max,invpt_max,2);
        double a = linpar0*linpar0;
        double rangemax = invpt_max;
        double linfuncmax = linpar1*rangemax + linpar0;
        double b = (linfuncmax*linfuncmax - a*a)/rangemax/rangemax;
        
        csqroot->cd();
        func->SetParameters(a,b);

        graph_sqrt->Fit(func,"QR");
        double sqroot_par0 = func->GetParameter(0);
        double sqroot_par1 = func->GetParameter(1);

        std::string sqroot_name ="sqroot_" + trackParam + "_" + iblname + "_eta" + std::to_string(ieta);
        const char *graph_sqroot_name = sqroot_name.c_str();
        TString sqroot_title(";Q/2p_{T} (MeV^{-1});#sigma(" + trackParam + ")" + trackParam_unit );
        TString YaxisSQname("#sigma(" + trackParam + ")"+ trackParam_unit );
        TString XaxisSQname("Q/2p_{T} [MeV]^{-1}");

        graph_sqrt->SetNameTitle(graph_sqroot_name,sqroot_title);
        graph_sqrt->SetFillColor(1);
        graph_sqrt->SetLineWidth(1);
        graph_sqrt->SetMarkerStyle(20);
        graph_sqrt->SetMarkerColor(2);

        graph_sqrt->GetYaxis()->SetTitle(YaxisSQname);
        graph_sqrt->GetXaxis()->SetTitle(XaxisSQname);
        graph_sqrt->SetTitle("");
        graph_sqrt->Draw("AEP");

        atlas_title->DrawLatex(.12, 0.85, atlas_title_name);
        eta_latex->DrawLatex(.12,0.95,eta_label);

        gStyle->SetOptFit(1111);

        csqroot->Update();
        TImage *imgsqroot = TImage::Create();
        TString sqroot_save( "output/" + sqroot_name + ".png");
        imgsqroot->FromPad(csqroot);
        imgsqroot->WriteImage(sqroot_save);

        ofile->Add(graph_sqrt);
        delete csqroot;
        delete func;
        delete graph_sqrt;
        if ( iibl == 0 ) {
	        sqrtnoIBL   << "nBLConsts.set(FTKTrackParam::" << std::left << std::setw(3) << trackParam << ", " << ieta << ", sqroot, " << std::scientific << sqroot_par0 <<  ",  " << std::scientific <<  TMath::Abs(sqroot_par1) << " );" << endl;
          linearnoIBL << "nBLConsts.set(FTKTrackParam::" << std::left << std::setw(3) << trackParam << ", " << ieta << ", linear, " << std::scientific << linpar0     <<  ",  " << std::scientific <<  linpar1                 << " );" << endl;
	      }else if (iibl == 1){
	        sqrtIBL     << "nomConsts.set(FTKTrackParam::" << std::left << std::setw(3) << trackParam << ", " << ieta << ", sqroot, " << std::scientific << sqroot_par0 <<  ",  " << std::scientific <<  TMath::Abs(sqroot_par1) << " );" << endl;
	        linearIBL   << "nomConsts.set(FTKTrackParam::" << std::left << std::setw(3) << trackParam << ", " << ieta << ", linear, " << std::scientific << linpar0     <<  ",  " << std::scientific <<  linpar1                 << " );" << endl;
        } 

      
      //delete graph_linear;
/////////////////////

//////////////////////////////////////////////////////////////////
////// CORE

      if(1){
        TCanvas *csqroot = new TCanvas("sqrtcore", "sqrtcore", 800, 650);

        TGraphErrors *graph_linear = new TGraphErrors(n,ipt_arr,widthcore_arr,0,widthcore_errors_arr);
        TF1 *linfit = new TF1("linearfit","pol1",0.0,invpt_max);
        graph_linear->Fit(linfit,"QR");
        linpar0 = linfit->GetParameter(0); 
        linpar1 = linfit->GetParameter(1);

        TGraphErrors *graph_sqrt = new TGraphErrors(n,ipt_arr,widthcore_arr,0,widthcore_errors_arr);

        TF1 *func = new TF1("sqrtfit",fitf,-invpt_max,invpt_max,2);
        double a = linpar0*linpar0;
        double rangemax = invpt_max;
        double linfuncmax = linpar1*rangemax + linpar0;
        double b = (linfuncmax*linfuncmax - a*a)/rangemax/rangemax;
        
        csqroot->cd();
        func->SetParameters(a,b);
        graph_sqrt->Fit(func,"QR");

        graph_sqrt->SetNameTitle(graph_sqroot_name,sqroot_title);
        graph_sqrt->SetFillColor(1);
        graph_sqrt->SetLineWidth(1);
        graph_sqrt->SetMarkerStyle(20);
        graph_sqrt->SetMarkerColor(2);

        graph_sqrt->GetYaxis()->SetTitle(YaxisSQname);
        graph_sqrt->GetXaxis()->SetTitle(XaxisSQname);
        graph_sqrt->SetTitle("");
        graph_sqrt->Draw("AEP");

        atlas_title->DrawLatex(.12, 0.85, atlas_title_name);
        eta_latex->DrawLatex(.12,0.95,eta_label);

        gStyle->SetOptFit(1111);

        csqroot->Update();
        TImage *imgsqroot = TImage::Create();
        std::string sqroot_name = "sqroot_core_" + trackParam + "_" + iblname + "_eta" + std::to_string(ieta);
        TString sqroot_save("output/" + sqroot_name + ".png");
        imgsqroot->FromPad(csqroot);
        imgsqroot->WriteImage(sqroot_save);

        ofile->Add(graph_sqrt);
        delete csqroot;
        
        delete func;
        delete graph_sqrt;
        delete graph_linear;
      }
/////////////////////


//////////////////////////////////////////////////////////////////
////// TAILS

      if(1){
        TCanvas *csqroot = new TCanvas("sqrttail", "sqrttail", 800, 650);

        TGraphErrors *graph_linear = new TGraphErrors(n,ipt_arr,widthtail_arr,0,widthtail_errors_arr);
        TF1 *linfit = new TF1("linearfit","pol1",0.0,invpt_max);
        graph_linear->Fit(linfit,"QR");
        linpar0 = linfit->GetParameter(0); 
        linpar1 = linfit->GetParameter(1);

        TGraphErrors *graph_sqrt = new TGraphErrors(n,ipt_arr,widthtail_arr,0,widthtail_errors_arr);

        TF1 *func = new TF1("sqrtfit",fitf,-invpt_max,invpt_max,2);
        double a = linpar0*linpar0;
        double rangemax = invpt_max;
        double linfuncmax = linpar1*rangemax + linpar0;
        double b = (linfuncmax*linfuncmax - a*a)/rangemax/rangemax;
        
        csqroot->cd();
        func->SetParameters(a,b);
        graph_sqrt->Fit(func,"QR");

        graph_sqrt->SetNameTitle(graph_sqroot_name,sqroot_title);
        graph_sqrt->SetFillColor(1);
        graph_sqrt->SetLineWidth(1);
        graph_sqrt->SetMarkerStyle(20);
        graph_sqrt->SetMarkerColor(2);

        graph_sqrt->GetYaxis()->SetTitle(YaxisSQname);
        graph_sqrt->GetXaxis()->SetTitle(XaxisSQname);
        graph_sqrt->SetTitle("");
        graph_sqrt->Draw("AEP");

        atlas_title->DrawLatex(.12, 0.85, atlas_title_name);
        eta_latex->DrawLatex(.12,0.95,eta_label);

        gStyle->SetOptFit(1111);

        csqroot->Update();
        TImage *imgsqroot = TImage::Create();
        sqroot_name = "sqroot_tail_" + trackParam + "_" + iblname + "_eta" + std::to_string(ieta);
        TString sqroot_save("output/" + sqroot_name + ".png");
        imgsqroot->FromPad(csqroot);
        imgsqroot->WriteImage(sqroot_save);

        ofile->Add(graph_sqrt);
        delete csqroot;
        delete func;
        delete graph_sqrt;
      }
/////////
//////////////////////////////////////

//////////////////////////////////////
//// RATIO PLOTS 

	    TCanvas *cratio = new TCanvas("ratio", "ratio", 800, 650);
	    TGraphErrors *graph_ratio = new TGraphErrors(n,ipt_arr,ratio_arr,0,ratio_errors_arr);
	    cratio->cd();

      std::string ratio_name = "ratio_" + trackParam + "_" + iblname + "_eta" + std::to_string(ieta);
      const char *graph_ratio_name = ratio_name.c_str();
      TString ratio_title(";Q/2p_{T} (MeV^{-1});n2/n1(" + trackParam + ")");
      TString Yaxisrationame("n2/n1(" + trackParam + ")"+ trackParam_unit );
	    TString Xaxisrationame("Q/2p_{T} [MeV]^{-1}");

      graph_ratio->SetNameTitle(graph_ratio_name,ratio_title);
	    graph_ratio->SetFillColor(1);
	    graph_ratio->SetLineWidth(1);
      graph_ratio->SetMarkerStyle(20);
      graph_ratio->SetMarkerColor(2);

	    graph_ratio->GetYaxis()->SetTitle(Yaxisrationame);
	    graph_ratio->GetXaxis()->SetTitle(Xaxisrationame);
      graph_ratio->SetTitle("");
      graph_ratio->Draw("AEP");
      
      atlas_title->DrawLatex(.12, 0.85, atlas_title_name);
      eta_latex->DrawLatex(.12,0.95,eta_label);

      gStyle->SetOptFit(1111);
      cratio->Update();
	    TImage *imgratio = TImage::Create();
      
      TString ratio_save("output/" + ratio_name + ".png");

	    imgratio->FromPad(cratio);
	    imgratio->WriteImage(ratio_save);

	    ofile->Add(graph_ratio);
      delete cratio;

      }
    }
  }

          //      std::cout << "pull37" << std::endl;

   for( int itp = 0; itp < nParams; itp++){
     for (int itp2 =0; itp2 < nParams;itp2++){
       TCanvas *canvas_highlow = new TCanvas;
       //TP_pull[itp][itp2][0]->Scale(2.0/TP_pull[itp][itp2][0]->Integral());
       TP_pull[itp][itp2][1]->Scale(1.0/TP_pull[itp][itp2][1]->Integral());
       TP_pull[itp][itp2][2]->Scale(1.0/TP_pull[itp][itp2][2]->Integral());

       TP_pull[itp][itp2][0]->SetLineColor(kRed);
       TP_pull[itp][itp2][1]->SetLineColor(kGreen+1);
       TP_pull[itp][itp2][2]->SetLineColor(kBlue+1);

       TP_pull[itp][itp2][0]->SetMinimum(0.0);
       TP_pull[itp][itp2][1]->SetMinimum(0.0);
       TP_pull[itp][itp2][2]->SetMinimum(0.0);

       TP_pull[itp][itp2][0]->SetMaximum(0.15);
       TP_pull[itp][itp2][1]->SetMaximum(0.15);
       TP_pull[itp][itp2][2]->SetMaximum(0.15);
//       TP_pull[itp][itp2][0]->Draw("HIST");
       TP_pull[itp][itp2][1]->Draw("HIST");
       TP_pull[itp][itp2][2]->Draw("HIST SAME");

       canvas_highlow->Update();
       TImage *imghighlow = TImage::Create();
       imghighlow->FromPad(canvas_highlow);
       string trackParam1 = trackParam_list[itp];
       string trackParam2 = trackParam_list[itp2];
       string hist_res_name_string = "output/TPpull" +trackParam1 +"_" + trackParam2 +".png";
       TString hist_res_name(hist_res_name_string);
       gStyle->SetOptFit();
       canvas_highlow->Update();
       imghighlow->WriteImage(hist_res_name);
       delete canvas_highlow;
     }
   }

   

   for( int itp = 0; itp < nParams; itp++){
     pull_res[itp]->Fit("gaus","Q","",-5.,5.);

     TCanvas *c = new TCanvas;

     TPad *padpull1 = new TPad("padpull1", "padpull1", 0, 0.3, 1, 1.0);
     padpull1->SetBottomMargin(0.05); // Upper and lower plot are joined
     padpull1->SetGridx();         // Vertical grid
     padpull1->Draw();             // Draw the upper pad: pad1
     padpull1->cd();               // pad1 becomes the current pad


     padpull1->SetLogy();
     pull_res_sg[itp]->SetLineColor(kRed+1);
     pull_res[itp]->SetLineColor(kBlue+1);
     pull_res_sg[itp]->SetLineStyle(2);
     pull_res[itp]->SetLineStyle(1);

 
     pull_res[itp]->Draw();
     pull_res_sg[itp]->Draw("SAME");


      c->cd();
      TPad *padpull2 = new TPad("padpull2", "padpull2", 0, 0.05, 1, 0.3);
      padpull2->SetTopMargin(0);
      padpull2->SetBottomMargin(0.2);
      padpull2->SetGridx(); // vertical grid
      padpull2->SetGridy();
      padpull2->Draw();
      padpull2->cd();       // pad2 becomes the current pad

       TH1F *hpull3 = (TH1F*)pull_res[itp]->Clone("hpull3");
       hpull3->GetFunction("gaus")->SetBit(TF1::kNotDraw);
       hpull3->SetLineColor(kBlack);
       hpull3->SetMinimum(0);  // Define Y ..
       hpull3->SetMaximum(2); // .. range
       hpull3->Sumw2();
       hpull3->SetStats(0);      // No statistics on lower plot
       //hpull3->Divide(gaussian);
       hpull3->Divide(pull_res_sg[itp]);
       hpull3->SetMarkerStyle(7);
       hpull3->Draw("ep");       // Draw the ratio plot
       hpull3->SetMinimum(0);  // Define Y ..
       hpull3->SetMaximum(2); // .. range

       hpull3->GetYaxis()->SetTitle("DG/SG");
       hpull3->GetYaxis()->SetNdivisions(505);
       hpull3->GetYaxis()->SetTitleSize(20);
       hpull3->GetYaxis()->SetTitleFont(43);
       hpull3->GetYaxis()->SetTitleOffset(1.0);
       hpull3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
       hpull3->GetYaxis()->SetLabelSize(15);

       // X axis ratio plot settings
       hpull3->GetXaxis()->SetTitleSize(20);
       hpull3->GetXaxis()->SetTitleFont(43);
       hpull3->GetXaxis()->SetTitleOffset(4.);
       hpull3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
       hpull3->GetXaxis()->SetLabelSize(15);
       
       padpull2->Update();


        
       c->cd();

       c->Update();

     padpull1->Update();
     TImage *img = TImage::Create();
     img->FromPad(c);
     string trackParam = trackParam_list[itp];
     string hist_res_name_string = "output/pull_res" + trackParam +".png";
     TString hist_res_name(hist_res_name_string);
     gStyle->SetOptFit();
     c->Update();
     img->WriteImage(hist_res_name);

     ofile->Add(pull_res[itp]);
     ofile->Add(pull_res_sg[itp]);
     delete c;

   }

  ofile->ls();
  ofile->Write();
  ofile->Close();


}


/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////





/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////


int main(int argc, char **argv) {
  int events;
  std::string output;
  std::string outputfolder;
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
      ("outputfolder,f",po::value<std::string>(&outputfolder)->default_value("output"),"The name of the outputfolder")
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
    std::cout << "TEST::INIT1" << std::endl;
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


    t_ftkdata->SetBranchAddress(Form("FTKMergedTracksStream"),&ftktracks);
    t_truth->SetBranchAddress("TruthTracks",&truthtracks);

    Int_t nloop = t_ftkdata->GetEntries();
    if (events > 0) {
      nloop = events;
    }
    if (nloop > t_ftkdata->GetEntries()) {
      nloop = t_ftkdata->GetEntries();
    }
    std::cout << "TEST::PROCESS" << std::endl;
    for (int ientry = 0; ientry < nloop; ++ientry) {
      Process(ientry);
    }
    reprocessing = 1;
    std::cout << "TEST::INIT2" << std::endl;
    Init2();
    std::cout << "TEST::PROCESS2" << std::endl;
    for (int ientry = 0; ientry < nloop; ++ientry) {
      Process(ientry);
    }
    std::cout << "TEST::WIDTH CALCULATION" << std::endl;
    width_calculation();
    std::cout << "TEST::SQROOT" << std::endl;
    sqrt_fit();
      std::cout << "TEST::CORE" << std::endl;
    coresqrt_fit();
        std::cout << "TEST::TAIL" << std::endl;

    tailsqrt_fit();
    std::cout << "TEST::PULL" << std::endl;

    for (int ientry = 0; ientry < nloop; ++ientry){
      Pull(ientry);
    }
    std::cout << "TEST::TERMINATE" << std::endl;

    Terminate(output,outputfolder);
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


          // SINGLE GAUSSIAN 

          // fit gaussian 
          // retreive width from gaussian 
          // plot 

          //width 
          //width error

          //DOUBLE GAUSSIAN 

          // fit the core gaussian 
          // retreive this as an initial parameter 
          // fit f1 with the core gaussian as initial parameters. 
          // set g2 parameter limits. 
          // retreive width from f1 
          // plot core  
          // plot tail 

          //width 
          //width error

          //core width 
          //core error

          //tail width 
          //tail error

