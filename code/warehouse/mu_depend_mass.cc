#include<iostream>
#include<fstream>
#include<string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCut.h"
#include "TEventList.h"

#include "/nfs_home/home1/myamatan/RootUtils/AtlasStyle.C"
#include "func.hh"
#include "sigList.hh"

//TString filename = "../input/merged/signal"+tagger+"/signal"+tagger+"."+ch+".root";
TString mass="0700";
TString filename = "../input/merged/signalZ/RS_G_ZZm"+mass+".2Lep.root";

void mu_depend_mass(){

	SetAtlasStyle();

	TH1F *h_hm; // mu > 30
	TH1F *h_lm; // mu < 30

	TFile *f;
	TTree *tr;

	// Fill

	f = new TFile(filename);
	tr = (TTree*)f->Get("Nominal");

	h_lm= new TH1F("lowmu","lowmu",25,50,150);
	h_hm= new TH1F("highmu","highmu",25,50,150);
	
	SetBranch(tr,"2Lep");
	
	TCut precut = returnTCut(preSele, "2Lep");
	tr->Draw(">>elist",precut);
	TEventList *elist = (TEventList*)gROOT->FindObject("elist");
	long nEvent = elist->GetN();
	for(long l=0; l<nEvent; l++){
		tr->GetEntry(elist->GetEntry(l));
		if( mu <= 30 ) h_lm->Fill( m );
		if( mu > 30 ) h_hm->Fill( m );
	}

	// Plot
	TCanvas *cv;
	TLegend *leg;
	TLatex *lat;
	TLatex *latsub;
	TLatex *latsub2;
	cv = new TCanvas("c","c2",200,10,600,500);
	leg= new TLegend(0.65,0.6,0.88,0.8);
	leg->SetBorderSize(0);
	lat = new TLatex(55, 0.16, "#it{#bf{ATLAS}} #bf{#bf{Simulation Internal}}");
	latsub = new TLatex(55, 0.145, "#scale[0.8]{#bf{#bf{after baseline selection}}}");
	latsub2 = new TLatex(55, 0.13, "#scale[0.8]{#bf{#bf{RS_G "+mass+" GeV}}}");

	h_hm->Scale(1./h_hm->Integral());
	h_lm->Scale(1./h_lm->Integral());
	float k_ymax = h_lm->GetBinContent( h_lm->GetMaximumBin() );

	cv->cd();

	h_lm->SetStats(0);
	h_lm->GetYaxis()->SetRangeUser(0,k_ymax*1.4);
	h_lm->GetYaxis()->SetTitle("A.U");
	h_lm->GetXaxis()->SetTitle("Jet mass [GeV]");
	h_lm->SetLineColor(kBlue);
	h_lm->SetLineWidth(2);
	leg->AddEntry(h_lm, "#mu < 30", "l");

	h_hm->SetStats(0);
	h_hm->GetYaxis()->SetTitle("A.U");
	h_hm->GetXaxis()->SetTitle("Jet mass [GeV]");
	h_hm->SetLineColor(kRed);
	h_hm->SetLineWidth(2);
	leg->AddEntry(h_hm, "#mu >= 30", "l");

	h_lm->Draw("hist");
	h_hm->Draw("hist same");

	leg->Draw();
	lat->Draw();
	latsub->Draw();
	latsub2->Draw();

	cv->Print("../output/mass_mu_dependent."+mass+".pdf");
	//cv[i]->Print("../output/signalEff." + boson.at(i) + "." + nameOption + ".commonWP.png");

}
