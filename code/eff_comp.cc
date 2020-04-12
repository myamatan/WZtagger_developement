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
#include "sigList.hh"
#include "func.hh"
//#include "func_old.hh"

std::vector<std::vector<TString> >  jet = {
					    {"tcc","MaxVVSemi_v5","TCC MaxSignificance"},
					    {"lct","Fix50VVJJ_LCT","LCTopo 50\% efficiency"}
};


TString  boson = "W";
const UInt_t jetsize = jet.size();

float mcSum=0.;

void eff_comp(){
	
	ch = "0Lep";

	SetAtlasStyle();

	TH1F *hbef[jetsize];
	TH1F *haft[jetsize];
	TH1F *heff[jetsize];

	TF1 *fD2;
	TF1 *fML;
	TF1 *fMU;

	TFile * f[jetsize];
	TTree *tr[jetsize];

	TString wpD2Name;
	TString wpMUName;
	TString wpMLName;

	// Fill
	for(UInt_t i=0; i<jetsize; i++){

		wpD2Name = jet.at(i).at(1);
		wpMUName = jet.at(i).at(1);
		wpMLName = jet.at(i).at(1);

		TString fname = "../input/merged_"+jet.at(i).at(0)+"_new/signal"+ boson + "/signal" + boson + "." + ch + ".root";
		std::cout << "++++++ " << fname << " ++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
		f[i] = new TFile(fname);
		tr[i] = (TTree*)f[i]->Get("Nominal");
		hbef[i]= new TH1F(Form("hb%i",i),Form("hb%i",i), nPtbins, ptRange);
		haft[i]= new TH1F(Form("ha%i",i),Form("ha%i",i), nPtbins, ptRange);
		heff[i]= new TH1F(Form("heff%i",i),Form("heff%i",i), nPtbins, ptRange);
		
		// D2
		std::string wpfile = std::string("../output/WPs/SmoothedContained"+boson+wps[wpD2Name]);
		std::vector<TString> cutWPs = retrieveWPs( wpfile );
		fD2 = new TF1("D2", cutWPs.at(2), 200, 2500);
		cutWPs.clear();
		// MU
		wpfile = std::string("../output/WPs/SmoothedContained"+boson+wps[wpMUName]);
		cutWPs = retrieveWPs( wpfile );
		fMU = new TF1("MU", cutWPs.at(1), 200, 2500);
		cutWPs.clear();
		// ML
		wpfile = std::string("../output/WPs/SmoothedContained"+boson+wps[wpMLName]);
		cutWPs = retrieveWPs( wpfile );
		fML = new TF1("ML", cutWPs.at(0), 200, 2500);
		cutWPs.clear();

		SetBranch(tr[i], ch);
		
		float D2cut=0.;
		float MUcut=0.;
		float MLcut=0.;
		TCut precut = returnTCut(preSele, ch);
		tr[i]->Draw(">>elist",precut);
		TEventList *elist = (TEventList*)gROOT->FindObject("elist");
		long nEvent = elist->GetN();
		for(long l=0; l<nEvent; l++){
			tr[i]->GetEntry(elist->GetEntry(l));
			// D2 & Mass window
			D2cut = pt > 2500 ? fD2->Eval(2500) : fD2->Eval(pt);
			MUcut = pt > 2500 ? fMU->Eval(2500) : fMU->Eval(pt);
			MLcut = pt > 2500 ? fML->Eval(2500) : fML->Eval(pt);
			hbef[i]->Fill( pt );
			if( d2<D2cut && m>MLcut && m<MUcut ){
				haft[i]->Fill( pt );
			}
			if( l%(int(nEvent/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEvent, float(float(l)/nEvent)*100) << std::endl;
		}
		delete fD2;
		delete fMU;
		delete fML;

	}// bosonSize

	// Efficiency calculation
	for(UInt_t i=0; i<jetsize; i++){
			setEfficiency(hbef[i], haft[i], heff[i]);
	}


	// Plot
	TCanvas *cv;
	TLegend *leg;
	TLatex  *lat;
	TLatex  *latsub;
	cv = new TCanvas("a","a",200,10,600,500);
	//cv = new TCanvas();
	leg= new TLegend(0.4,0.2,0.88,0.4);
	leg->SetBorderSize(0);
	lat = new TLatex(300, 1.1, "#it{#bf{ATLAS}} #bf{#bf{Simulation Work in progress}}");
	latsub = new TLatex(350, 0.95, "#bf{"+boson+" tagger}");
	cv->cd();
	for(UInt_t i=0; i<jetsize; i++){
		if( i==0 ) heff[i]->Draw();
		else heff[i]->Draw("same");
		heff[i]->SetStats(0);
		heff[i]->GetYaxis()->SetRangeUser(0,1.2);
		heff[i]->GetYaxis()->SetTitle("Signal efficiency");
		heff[i]->GetXaxis()->SetTitle("Jet pt [GeV]");
		heff[i]->SetMarkerColor(colorList.at(i));
		heff[i]->SetMarkerSize(1);
		heff[i]->SetMarkerStyle(8);
		heff[i]->SetLineColor(colorList.at(i));
		heff[i]->SetLineWidth(2);
		leg->AddEntry(heff[i], jet.at(i).at(2), "lep");

	}
	leg->Draw();
	lat->Draw();
	latsub->Draw();
	cv->Print("../output/signalEff_comp." + ch + "_" + boson + "_" + nameOption + ".eps");

}
