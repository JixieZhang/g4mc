//This script is used to plot wacs events
//It requires that proton is the primary 1 particle and gamma is the primary 2 particle
void draw_HMSCut(double p0=3.85)
{
  TTree *proton = (TTree*)gDirectory->Get("track0");
  TTree *gamma = (TTree*)gDirectory->Get("track1");
  proton->AddFriend(gamma,"gamma");


  char pcut[500],gpcut[500];
  sprintf(pcut,"abs(track0.Pvb-%.3f)<%.3f",p0,p0*0.09);
  sprintf(gpcut,"track1.Pvb>0 && abs(track0.Pvb-%.3f)<%.3f",p0,p0*0.09);


  track1->AddFriend("track0");
  track0->SetAlias("myphi","(Phi0<0)?Phi0+2*3.14159:Phi0");
  track1->SetAlias("track0.myphi","(track0.Phi0<0)?track0.Phi0+2*3.14159:track0.Phi0");

  track1->Draw("track1.Theta0*57.3:track0.Theta0*57.3",pcut,"");
  TGraph *gr_t1vst2=(TGraph*) gROOT->FindObject("Graph")->Clone("gr_t1vst2");
  gr_t1vst2->SetMarkerColor(1);
  gr_t1vst2->SetTitle("Pr_detected(black) and #gamma-p Coincident(color);#theta_{pr};#theta_{#gamma}");
  gr_t1vst2->Draw("pA");
  track1->Draw("track1.Theta0*57.3:track0.Theta0*57.3",gpcut,"colrsame");
  c1->SaveAs("Graph/Theta_gamma_pr.png"); 


  track1->Draw("track0.Theta0*57.3>>ht1",pcut,"");
  ht1->SetLineColor(1); ht1->SetTitle("Pr_detected(black) and #gamma-p Coincident(red);#theta_{pr}");
  track1->Draw("track0.Theta0*57.3>>ht2",gpcut,"same");
  ht2->SetLineColor(2);
  c1->SaveAs("Graph/Theta_pr.png"); 

  track1->Draw("track0.myphi*57.3>>hp1",pcut,"");
  hp1->SetLineColor(1); hp1->SetTitle("Pr_detected(black) and #gamma-p Coincident(red);#phi_{pr}");
  track1->Draw("track0.myphi*57.3>>hp2",gpcut,"same");
  hp2->SetLineColor(2);
  c1->SaveAs("Graph/Phi_pr.png"); 

  track1->Draw("track1.Theta0*57.3>>ht1_g",pcut,"");
  ht1_g->SetLineColor(1); ht1_g->SetTitle("Pr_detected(black) and #gamma-p Coincident(red);#theta_{#gamma}");
  track1->Draw("track1.Theta0*57.3>>ht2_g",gpcut,"same");
  ht2_g->SetLineColor(2);
  c1->SaveAs("Graph/Theta_gamma.png"); 

  track1->Draw("track1.Phi0*57.3>>hp1_g",pcut,"");
  hp1_g->SetLineColor(1); hp1_g->SetTitle("Pr_detected(black) and #gamma-p Coincident(red);#phi_{#gamma}");
  track1->Draw("track1.Phi0*57.3>>hp2_g",gpcut,"same");
  hp2_g->SetLineColor(2);
  c1->SaveAs("Graph/Phi_gamma.png"); 


  track1->Draw("Ei>>hei0","","");
  hei0->SetLineColor(3); hei0->SetTitle("Thrown(green), Pr_detected(black) and #gamma-p Coincident(red); E_{#gamma}");
  track1->Draw("Ei>>hei1",pcut,"same"); hei1->SetLineColor(1);
  track1->Draw("Ei>>hei2",gpcut,"same");
  hei2->SetLineColor(2);
  c1->SaveAs("Graph/Ei_gamma.png"); 

  track1->Draw("track1.P0:track1.Theta0*57.3",pcut,"");
  TGraph *gr_pvst_g=(TGraph*) gROOT->FindObject("Graph")->Clone("gr_pvst_g");
  gr_pvst_g->SetMarkerColor(1);
  gr_pvst_g->SetTitle("Pr_detected(black) and #gamma-p Coincident(color);#theta_{#gamma};E'_{#gamma}");
  gr_pvst_g->Draw("pA");
  track1->Draw("track1.P0:track1.Theta0*57.3",gpcut,"colrsame");
  c1->SaveAs("Graph/PTheta_gamma.png"); 

  track1->Draw("track0.P0:track0.Ei",pcut,"");
  TGraph *gr_pvsei=(TGraph*) gROOT->FindObject("Graph")->Clone("gr_pvsei");
  gr_pvsei->SetMarkerColor(1);
  gr_pvsei->SetTitle("Pr_detected(black) and #gamma-p Coincident(color);Ei_{#gamma};P_{pr}");
  gr_pvsei->Draw("pa");
  track1->Draw("track0.P0:track0.Ei",gpcut,"colrsame");
  c1->SaveAs("Graph/PprVsEi.png"); 


  track1->Draw("-Xvb_tr>>hxvb_tr","Pvb>0","");
  hxvb_tr->SetLineColor(1); hxvb_tr->SetTitle("Photon detected(black) and #gamma-p Coincident(red); -Xvb_tr(mm)");  
  track1->Draw("-Xvb_tr>>hxvb_tr2",gpcut,"same");
  hxvb_tr2->SetLineColor(2);
  c1->SaveAs("Graph/Xvb_tr_g.png"); 

  track1->Draw("Yvb_tr>>hyvb_tr","Pvb>0","");
  hyvb_tr->SetLineColor(1); hyvb_tr->SetTitle("Photon detected(black) and #gamma-p Coincident(red); Yvb_tr(mm)");  
  track1->Draw("Yvb_tr>>hyvb_tr2",gpcut,"same"); 
  hyvb_tr2->SetLineColor(2);
  c1->SaveAs("Graph/Yvb_tr_g.png"); 

  TCut mycut = pcut;
  mycut += "track1.Theta0_tr>-1";
  track1->Draw("track1.Theta0_tr:track1.Phi0_tr>>h2tp(100,-0.5,0.5,100,-0.5,0.5)",mycut,"");
  TGraph *gr_thetavsphi_g=(TGraph*) gROOT->FindObject("Graph")->Clone("gr_thetavsphi_g");
  gr_thetavsphi_g->SetMarkerColor(1);
  gr_thetavsphi_g->SetTitle("Photon detected(black) and #gamma-p Coincident(color);#phi_{#gamma}^{tr};#theta_{#gamma}^{tr}");
  gr_thetavsphi_g->Draw("pA");
  track1->Draw("track1.Theta0_tr:track1.Phi0_tr",gpcut,"colrsame");
  c1->SaveAs("Graph/ThetaPhi_tr_gamma.png"); 
}
