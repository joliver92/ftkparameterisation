#include <TMath.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <map>
#include <iostream>
#include "../TrigFTKSim/FTKTrackStream.h"
#include "../TrigFTKSim/FTKRoadStream.h"
#include "../TrigFTKSim/FTKTruthTrack.h"
#include "TChain.h"
#include "TFile.h"
#include "TProfile.h"

using namespace std;

class MatchInfo {
private:
  int m_barcode;
  int m_evtindex;

public:
  MatchInfo() : m_barcode(0), m_evtindex(-1) {;}
  MatchInfo(int v1, int v2) : m_barcode(v1), m_evtindex(v2) {;}
  bool operator==(const MatchInfo& o) const { return (m_barcode==o.m_barcode)&&(m_evtindex==o.m_evtindex); }
  bool operator<(const MatchInfo& o) const { if (m_evtindex!=o.m_evtindex) return (m_evtindex<o.m_evtindex); else return m_barcode<o.m_barcode; }
};

typedef multimap<MatchInfo,const FTKTrack*> FTKBarcodeMM;
unsigned nbins;
unsigned nbinspull;
double maxntracks;
double d0min, d0max;
double z0min, z0max;
double z0_min, z0_max;
double eta_min,eta_max;
double d0_min, d0_max;
double pt_min,pt_max;
double invpt_min,invpt_max;
double etabin_min;
double etabin_max;
double invptbin_min;
double invptbin_max;
double ipt_arr[200];
int nibl;



double stepsize;
int ninvptbins;
int middlebin1;
int middlebin2;
double invptval;
int nParams=5;
//int z0bins;
//int d0bins;
//int etabins;
//double z0step;
//double d0step;
//double etastep;

double phimin, phimax;
double etamin, etamax;
double abscurvmax;
double ptmax;
double Dd0;
double Dphi;
double Drelpt;
double Dcurv;
double Deta;
double Dz0;
double ptmincut;
double dx, dy;
Double_t ptbins[20];
Double_t etabinsarray[20];
Double_t d0binsarray[20];
Double_t z0binsarray[20];
Double_t ptbinsarray[20];
int reprocessing = 0;
//TString Tpulltype[3];
//TString Tpulltitle[3];

TH1F *pull_res[5];
TH1F *pull_res_sg[5];
TH1F *TP_pull[5][5][3];

// block of generic control histograms for the FTK tracks
TH1F *histoz0res_ptgt30;
TH1F *histoz0res_ptlt30;

TH2F *histocoordmasketa_ftk;
TH2F *histocoordmaskz0_ftk;
TH2F *histocoordmaskphi_ftk;
TH2F *histonmisseta_ftk;
TH2F *histochisqndfeta_ftk;

TH1F *histontracks_ftk;
TH1F *histod0_ftk;
TH1F *histod0_ftks[20];

TH1F *historesz0_hlt[20];
TH1F *historesd0_hlt[20];
TH1F *historeseta_hlt[20];

TH1F *hist_resz0_IBL[20];
TH1F *hist_resd0_IBL[20];
TH1F *hist_reseta_IBL[20];

TH1F *hist_resz0_noIBL[20];
TH1F *hist_resd0_noIBL[20];
TH1F *hist_reseta_noIBL[20];

TH1F *hist_resd0_IBL_test[30][30];
TH1F *hist_resd0_IBL_test_invpt[30][30];
TH1F *hist_resd0_noIBL_test_invpt[30][30];

TH1F *hist_res[10][2][20][200];
TH1F *hist_res_temp[10][2][20][200];

//corewidths.push_back(corewidth);
//tailwidths.push_back(tailwidth);

//corewidth_errors.push_back(corewidth_error);
//tailwidth_errors.push_back(tailwidth_error);

double A_corewidth[10][2][20][200];
double A_tailwidth[10][2][20][200];

double A_ratio[10][2][20][200];

double A_corewidth_errors[10][2][20][200];
double A_tailwidth_errors[10][2][20][200];



double A_width[10][2][20][200];
double A_width_errors[10][2][20][200];
double A_sq_par0[10][2][20];
double A_sq_par1[10][2][20];

double A_coresq_par0[10][2][20];
double A_coresq_par1[10][2][20];
double A_tailsq_par0[10][2][20];
double A_tailsq_par1[10][2][20];

double hist_std[10][2][20][200];

TProfile *hist_res_vs[10][2][20][200][10];

//string iblnames[2];
//string trackParam_list[10];
string iblnames[] = {"noIBL","IBL"};
string trackParam_list[] = {"d0","z0","eta","phi","Ipt"};
double trackParam_range[5];
std::vector<double>  TP_ranges_max;
std::vector<double>  TP_ranges_min;
std::vector<double>  TP_ranges_res;
string trackParam_units[] = {"[mm]","[mm]"," "," ","[MeV]^{-1}"};
//TGraphErrors *g[20];

TH1F *pthighpull;
TH1F *ptlowpull;


TH1F *histod0_ftk_poseta;
TH1F *histod0_ftk_negeta;
TH1F *histod0_ftk_pt2;
TH1F *histod0_ftk_pt5;
TH1F *histoz0_ftk;
TH1F *histophi_ftk;
TH2F *histophid0_ftk;
TH2F *histophid0_ftk_pt2;
TH2F *histophid0_ftk_pt5;
TH1F *histocurv_ftk;
TH1F *histoeta_ftk;
TH2F *histoetaphi_ftk;


TH2F *histoetaphi_ftk_IBL;
TH2F *histoetaphi_ftk_PixL0;
TH2F *histoetaphi_ftk_PixL1;
TH2F *histoetaphi_ftk_PixL2;
TH2F *histoetaz0_ftk;
TH2F *histoetaz0_ftk_IBL;
TH2F *histoetaz0_ftk_PixL0;
TH2F *histoetaz0_ftk_PixL1;
TH2F *histoetaz0_ftk_PixL2;
TH2F *histoetaphi_truth;
TH2F *histoetaphi_truthM;
TH2F *histoetaz0_truth;
TH2F *histoetaz0_truthM;
TH1F *histopt_ftk;

TH1F *histopt_ftk_lg;
TH1F *histopt_ftklo_lg;

TH2F *histoetaz0det_ftk_IBL;
TH2F *histoetaz0det_ftk_PixL0;
TH2F *histoetaz0det_ftk_PixL1;
TH2F *histoetaz0det_ftk_PixL2;

TH2F *histophiz0_ftk;
TH2F *histophiz0_ftk_IBL;
TH2F *histophiz0_ftk_PixL0;
TH2F *histophiz0_ftk_PixL1;
TH2F *histophiz0_ftk_PixL2;


TH1F *histopt_ftkzoom;
TH1F *histocurv_ftkzoom;

// FTK for fakes
TH1F *histontracks_goodftk;


TH1F *histopt_goodftk_lg;
TH1F *histopt_goodftklo_lg;

TH1F *histopt_goodftkUlo_lg;
TH1F *histopt_goodftkU_lg;

TH1F *histoetaabs_truth;
TH1F *histoeff_truth;
TH1F *histopt_truthlo_lg;
TH1F *histopt_truth_lg;
TH1F *histoetaabs_truthM;
TH1F *histoeff_truthM;

TH1F *histopt_truthMlo_lg;
TH1F *histopt_truthM_lg;


TH1F *histod0_goodftk;
TH1F *histoz0_goodftk;
TH1F *histocurv_goodftk;
TH1F *histoeta_goodftk;
TH1F *histophi_goodftk;
TH1F *histopt_goodftk;

TH1F *histontracks_goodftkU;
TH1F *histod0_goodftkU;
TH1F *histoz0_goodftkU;
TH1F *histocurv_goodftkU;
TH1F *histoeta_goodftkU;
TH1F *histophi_goodftkU;
TH1F *histopt_goodftkU;

// block of distribution related to truth tracks
TH1F *histontracks_truth;
TH1F *histod0_truth;
TH1F *histoz0_truth;
TH1F *histocurv_truth;
TH1F *histoeta_truth;
TH1F *histophi_truth;
TH1F *histopt_truth;

TH1F *histontracks_truthM;
TH1F *histod0_truthM;
TH1F *histoz0_truthM;
TH1F *histocurv_truthM;
TH1F *histoeta_truthM;
TH1F *histophi_truthM;
TH1F *histopt_truthM;
TH1F *histod0res;
TH1F *histoz0res;

TProfile *histod0res_veta;
TProfile *histoz0res_veta;
TProfile *histod0res_vphi;
TProfile *histoz0res_vphi;
TProfile *histod0res_vz0;
TProfile *histoz0res_vz0;

TH1F *histocurvres;
TH1F *histoetares;
TH1F *histophires;
TH1F *historelptres;
TH2F *historelptrespt;

// block of distribution related to truth muon tracks
TH1F *histontracks_truth_muon;
TH1F *histod0_truth_muon;
TH1F *histoz0_truth_muon;
TH1F *histocurv_truth_muon;
TH1F *histoeta_truth_muon;
TH1F *histophi_truth_muon;
TH1F *histopt_truth_muon;
TH1F *histopt_truth_muonlo_lg;
TH1F *histopt_truth_muon_lg;

TH1F *histontracks_truthM_muon;
TH1F *histod0_truthM_muon;
TH1F *histoz0_truthM_muon;
TH1F *histocurv_truthM_muon;
TH1F *histoeta_truthM_muon;
TH1F *histophi_truthM_muon;
TH1F *histopt_truthM_muon;
TH1F *histopt_truthM_muonlo_lg;
TH1F *histopt_truthM_muon_lg;

// Things to access variables!
FTKTrackStream *ftktracks(0);
FTKRoadStream *roads(0);
std::vector<FTKTruthTrack> *truthtracks(0);
//std::vector<double> z0array;
//std::vector<double> d0array;
//std::vector<double> etaarray;

std::vector<double> d0bins;
std::vector<double> z0bins;
std::vector<double> etabins;
std::vector<double> ptbinsvec;
std::vector<double> invptbinsvec;
std::vector<double> widthvec;
std::vector<double> widtherr;

std::vector<double> ratios;

std::vector<double> corewidths;
std::vector<double> tailwidths;
std::vector<double> corewidth_errors;
std::vector<double> tailwidth_errors;
//tailwidths.push_back(tailwidth);

//corewidth_errors.push_back(corewidth_error);
//tailwidth_errors.push_back(tailwidth_error);



std::vector<double> widths;
std::vector<double> width_errors;
Int_t RunNumber, EventNumber;
TChain *t_ftkdata;
TChain *t_truth;
TChain *t_evtinfo;

int towerNumber;
/* TString outputname; */
std::string outputname;
std::string psfile;
int ientry2;
Int_t Use1stStage;
Double_t pTmin;
Double_t pTmax;
