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

float mcSum=0.;


void mu_depend_eff(){

	SetAtlasStyle();

	TH1F *hbef[3];
	TH1F *haft[3];
	TH1F *heff[3];

	TF1 *fD2;
	TF1 *fML;
	TF1 *fMU;

	TFile *f;
	TTree *tr;

	TString wpName = "MaxVVSemi_v5";


	std::string ch_tagger = ch.Data();
	ch_tagger += "_";
	ch_tagger += tagger.Data();

	// Fill

	//TString fname = "../input/merged_new/signal"+tagger+"/signal"+tagger+".0Lep.root";
	TString fname = "../input/merged_new/sigDBL/sigDBL_"+ch+"_mc16ad.root";

	f = new TFile(fname);
	tr = (TTree*)f->Get("Nominal");
	for(int i=0; i<3; i++){
	  hbef[i] = new TH1F(Form("hb%i",i),Form("hb%i",i), nPtbins, ptRange);
	  haft[i] = new TH1F(Form("ha%i",i),Form("ha%i",i), nPtbins, ptRange);
	  heff[i] = new TH1F(Form("heff%i",i),Form("heff%i",i), nPtbins, ptRange);
	}
	
	// Smoothed cut parameters
	std::string wpfile = std::string("../output/WPs/SmoothedContained"+tagger+wps[wpName]);
	std::vector<TString> cutWPs = retrieveWPs( wpfile );
	fD2 = new TF1("D2", cutWPs.at(2), 200, 2500);
	fMU = new TF1("MU", cutWPs.at(1), 200, 2500);
	fML = new TF1("ML", cutWPs.at(0), 200, 2500);

	SetBranch(tr, ch);

	float D2cut=0.;
	float MUcut=0.;
	float MLcut=0.;
	TCut precut = returnTCut(preSele, ch);
	tr->Draw(">>elist",precut);
	TEventList *elist = (TEventList*)gROOT->FindObject("elist");
	long nEvent = elist->GetN();
	for(long l=0; l<nEvent; l++){
		tr->GetEntry(elist->GetEntry(l));
		if( std::find(sigDSID[ch_tagger].begin(), sigDSID[ch_tagger].end(), dsid) != sigDSID[ch_tagger].end() ){
			// D2 & Mass window
			D2cut = pt > 2500 ? fD2->Eval(2500) : fD2->Eval(pt);
			MUcut = pt > 2500 ? fMU->Eval(2500) : fMU->Eval(pt);
			MLcut = pt > 2500 ? fML->Eval(2500) : fML->Eval(pt);
			// Denominator
			if( mu < 25 ){
				hbef[0]->Fill( pt );
			}else if( mu < 40 ){
				hbef[1]->Fill( pt );
			}else{
				hbef[2]->Fill( pt );
			}
			if( d2<D2cut && m>MLcut && m<MUcut ){
				if( mu < 25 ){
			  		haft[0]->Fill( pt );
			  	}else if( mu < 40 ){
			  		haft[1]->Fill( pt );
			 	}else{
			  		haft[2]->Fill( pt );
			  	}
			}
		}
	}

	// Efficiency calculation
	for(UInt_t i=0; i<3; i++){
		setEfficiency(hbef[i], haft[i], heff[i]);
	}


	// Plot
	TCanvas *cv;
	TLegend *leg;
	TLatex *lat;
	TLatex *latsub;
	cv = new TCanvas("C","C2",200,10,600,500);
	leg= new TLegend(0.4,0.65,0.78,0.85);
	leg->SetBorderSize(0);
	lat = new TLatex(300, 1.1, "#it{#bf{ATLAS}} #bf{#bf{Simulation Internal}}");
	latsub = new TLatex(350, 0.95, "#bf{"+tagger+" tagger}");
	cv->cd();

	heff[0]->SetStats(0);
	heff[0]->GetYaxis()->SetRangeUser(0,1.2);
	heff[0]->GetYaxis()->SetTitle("Signal efficiency");
	heff[0]->GetXaxis()->SetTitle("Jet pt [GeV]");
	heff[0]->SetMarkerColor(kGreen-3);
	heff[0]->SetMarkerSize(1);
	heff[0]->SetMarkerStyle(8);
	heff[0]->SetLineColor(kGreen-3);
	heff[0]->SetLineWidth(2);
	leg->AddEntry(heff[0], "#mu < 25", "lp");

	heff[1]->SetStats(0);
	heff[1]->GetYaxis()->SetRangeUser(0,1.2);
	heff[1]->GetYaxis()->SetTitle("Signal efficiency");
	heff[1]->GetXaxis()->SetTitle("Jet pt [GeV]");
	heff[1]->SetMarkerColor(kBlue-3);
	heff[1]->SetMarkerSize(1);
	heff[1]->SetMarkerStyle(8);
	heff[1]->SetLineColor(kBlue-3);
	heff[1]->SetLineWidth(2);
	leg->AddEntry(heff[1], "25 < #mu < 40", "lp");

	heff[2]->SetStats(0);
	heff[2]->GetYaxis()->SetRangeUser(0,1.2);
	heff[2]->GetYaxis()->SetTitle("Signal efficiency");
	heff[2]->GetXaxis()->SetTitle("Jet pt [GeV]");
	heff[2]->SetMarkerColor(kRed-3);
	heff[2]->SetMarkerSize(1);
	heff[2]->SetMarkerStyle(8);
	heff[2]->SetLineColor(kRed-3);
	heff[2]->SetLineWidth(2);
	leg->AddEntry(heff[2], "40 < #mu", "lp");

	heff[0]->Draw("");
	heff[1]->Draw("same");
	heff[2]->Draw("same");
	leg->Draw();
	lat->Draw();
	latsub->Draw();

	cv->Print("../output/signalEff_mu_depend." + tagger + ".pdf");
	//cv[i]->Print("../output/signalEff." + boson.at(i) + "." + nameOption + ".commonWP.png");

}
