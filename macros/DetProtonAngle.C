//this script is used to find angle acceptance such that one can 
//determeine the spectrometer's angle
#include "stdlib.h"
#include <iostream>
#include "math.h"
using namespace std;

#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TQObject.h>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TString.h"
#include "TCut.h"
#include "TCutG.h"
#include "TPaveText.h"
#include "TText.h"
#include "TPad.h"
#include "TF1.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

double GetBeamEnergy()
{
  double pBeam=0.0;
  TTree *config = (TTree*) gROOT->FindObject("config");
  config->SetBranchAddress("Beam",&pBeam);
  config->GetEntry(0);
  cout<< "The beam energy for this run is "<< pBeam <<" GeV"<<endl;
  return pBeam;
}

double GetConfigLeaf(const char *leaf="VDAngle")
{
  double val=0.0;
  TTree *config = (TTree*) gROOT->FindObject("config");
  config->SetBranchAddress(leaf,&val);
  config->GetEntry(0);
  //cout<< "In this run,  "<< leaf <<" = "<<val<<endl;
  return val;
}


void DetPhotonAngle(double P0_HMS=0.0)
{
  gROOT->Reset();
  gStyle->SetOptStat(100);
  gStyle->SetStatH(0.03);
  gStyle->SetStatX(1.0-gStyle->GetPadRightMargin());


  TTree *track0 = (TTree*) gROOT->FindObject("track0");
  track0->AddFriend("track1");

  //double Beam = GetBeamEnergy();
  const double deg = atan(1.0)/45.;

  double PhotonAngle = GetConfigLeaf("VDAngle")/deg;
  bool bIsSBS = (PhotonAngle>180.)?true:false;

  double ProtonAngle = 0;
  if(bIsSBS) ProtonAngle = GetConfigLeaf("SuperBigBiteAngle")/deg;
  else ProtonAngle = GetConfigLeaf("HMSAngle")/deg;

  TCut pcut = "1000*ElasXS*P0*P0/3.14159*(track1.Pvb>0)";
  if (P0_HMS > 0) pcut=Form("1000*ElasXS*P0*P0/3.14159*(abs(track1.Pvb-%.3f)<%.3f)",P0_HMS,0.09*P0_HMS);

  double x1,y1,x2,y2;

  TCanvas *c61 = new TCanvas("c61","Determine Pronton arm angle",800,500);
  c61->Divide(2,1);
  c61->cd(1);
  track0->Draw("track0.Theta0*57.3>>hTg",pcut);
  TH1F* hTg = (TH1F*) gROOT->FindObject("hTg");
  hTg->SetTitle("Weigted by RCS XS, Pr detected;  #theta_{#gamma} (deg)");
  cout<<"Best Gamma angle is "<< hTg->GetMean()<<" deg\n";

  x1=x2=hTg->GetMean();
  y1=0;y2=hTg->GetMaximum()*1.03;
  TLine *v_w = new TLine(x1,y1,x2,y2);
  v_w->SetLineWidth(2);
  v_w->SetLineColor(2);
  v_w->Draw("same");


  c61->cd(2);
  track0->Draw("track0.P0>>hPg",pcut);
  TH1F* hPg = (TH1F*) gROOT->FindObject("hPg");
  hPg->SetTitle("Weigted by RCS XS, Pr detected;  E_{#gamma} (GeV)");

  c61->SaveAs(Form("Graph/BestGammaAngle_for_Pr%.1f.png",ProtonAngle));
}

void DetProtonAngle(double P0_HMS=0.0)
{
  gROOT->Reset();
  gStyle->SetOptStat(100);
  gStyle->SetStatH(0.03);
  gStyle->SetStatX(1.0-gStyle->GetPadRightMargin());

  TTree *track0 = (TTree*) gROOT->FindObject("track0");
  track0->AddFriend("track1");

  //double Beam = GetBeamEnergy();
  const double deg = atan(1.0)/45.;

  double PhotonAngle = GetConfigLeaf("VDAngle")/deg;
  bool bIsSBS = (PhotonAngle>180.)?true:false;

  const char* Detector = (bIsSBS)?"SBS":"HMS";

  TCut pcut = "(track1.Pvb>0)";
  if (P0_HMS > 0) pcut=Form("(abs(track1.Pvb-%.3f)<%.3f)",P0_HMS,0.09*P0_HMS);

  TCut pcut_weighted = "1000*ElasXS*P0*P0/3.14159*(track1.Pvb>0)";
  if (P0_HMS > 0) pcut=Form("1000*ElasXS*P0*P0/3.14159*(abs(track1.Pvb-%.3f)<%.3f)",P0_HMS,0.09*P0_HMS);

  double ProtonAngle = 0;
  if(bIsSBS) ProtonAngle = GetConfigLeaf("SuperBigBiteAngle")/deg;
  else ProtonAngle = GetConfigLeaf("HMSAngle")/deg;

  double x1,y1,x2,y2;

  TCanvas *c61 = new TCanvas("c61","Determine Pronton arm angle",800,450);
  c61->Divide(2,1);
  c61->cd(1);
  track0->Draw("track1.Theta0*57.3>>hTp_w",pcut_weighted);
  TH1F* hTp_w = (TH1F*) gROOT->FindObject("hTp_w");
  hTp_w->SetTitle("Weigted by RCS XS, #gamma detected;  #theta_{pr} (deg)");
  cout<<"Best weighted proton angle is "<< hTp_w->GetMean()<<" deg\n";
  
  x1=x2=hTp_w->GetMean();
  y1=0;y2=hTp_w->GetMaximum()*1.03;
  TLine *v_w = new TLine(x1,y1,x2,y2);
  v_w->SetLineWidth(2);
  v_w->SetLineColor(2);
  v_w->Draw("same");
 
  c61->cd(2);
  track0->Draw("track1.Theta0*57.3>>hTp",pcut);
  TH1F* hTp = (TH1F*) gROOT->FindObject("hTp");
  hTp->SetTitle("Not Weigted, #gamma detected;  #theta_{pr} (deg)");
  cout<<"Best not weighted proton angle is "<< hTp->GetMean()<<" deg\n";
  
  x1=x2=hTp->GetMean();
  y1=0;y2=hTp->GetMaximum()*1.03;
  TLine *v = new TLine(x1,y1,x2,y2);
  v->SetLineWidth(2);
  v->SetLineColor(2);
  v->Draw("same");

 
  c61->SaveAs(Form("Graph/Best%sAngle_for_Gamma%.1f.png",Detector,PhotonAngle));
}
