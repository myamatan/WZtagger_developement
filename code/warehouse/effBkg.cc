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

std::vector<TString>  boson = {"Z", "W"};
//std::vector<TString>  chs = {"0Lep","1Lep","2Lep"};
std::vector<TString>  chs = {"0Lep"};
const UInt_t bosonSize = boson.size();
const UInt_t chsSize = chs.size();

float mcSum=0.;

void effBkg(){

	SetAtlasStyle();

	TH1F *hbef[bosonSize][chsSize];
	TH1F *haft[bosonSize][chsSize];
	TH1F *heff[bosonSize][chsSize];

	TF1 *fD2[bosonSize];
	TF1 *fML[bosonSize];
	TF1 *fMU[bosonSize];

	TFile *f[chsSize];
	TTree *tr[chsSize];

	TString wpName = "MaxVVSemi_v5";
	//TString wpName = "Fix50VVJJ_LCT";

	// Fill
	for(UInt_t j=0; j<chsSize; j++){

		//TString fname = "../input/merged_new/allMC/allMC." + chs.at(j) + ".root";
		//TString fname = "../input/merged_lct_new/allMC/allMC." + chs.at(j) + ".root";
		TString fname = "../input/merged_tcc_new/allMC/allMC." + chs.at(j) + ".root";

		std::cout << "++++++ " << fname << " ++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

		f[j] = new TFile(fname);
		tr[j] = (TTree*)f[j]->Get("Nominal");
		
		SetBranch(tr[j], chs.at(j));

		for(UInt_t i=0; i<bosonSize; i++){
			hbef[i][j]= new TH1F(Form("hb%i-%i",i,j),Form("hb%i-%i",i,j), nPtbins, ptRange);
			haft[i][j]= new TH1F(Form("ha%i-%i",i,j),Form("ha%i-%i",i,j), nPtbins, ptRange);
			heff[i][j]= new TH1F(Form("heff%i-%i",i,j),Form("heff%i-%i",i,j), nPtbins, ptRange);
			
			// Smoothed cut parameters
			//std::string wpfile = std::string("../output/WPs/SmoothedContained"+boson.at(i)+"Tagger_AntiKt10TrackCaloClusterTrimmed_"+chs.at(0)+"_"+sigType+"_"+nameOption+".txt");
			std::string wpfile = std::string("../output/WPs/SmoothedContained"+boson.at(i)+wps[wpName]);
			//std::string wpfile = std::string("../output/WPs/Smoothed"+boson.at(i)+wps[wpName]);
			std::vector<TString> cutWPs = retrieveWPs( wpfile );

			if( j==0 ){
				fD2[i] = new TF1("D2"+boson.at(i), cutWPs.at(2), 200, 2500);
				fMU[i] = new TF1("MU"+boson.at(i), cutWPs.at(1), 200, 2500);
				fML[i] = new TF1("ML"+boson.at(i), cutWPs.at(0), 200, 2500);
			}
		}

		float D2cut=0.;
		float MUcut=0.;
		float MLcut=0.;
		TCut precut = returnTCut(preSele, chs.at(j));
        	tr[j]->Draw(">>elist",precut);
        	TEventList *elist = (TEventList*)gROOT->FindObject("elist");
        	long nEvent = elist->GetN();
		for(long l=0; l<nEvent; l++){
			tr[j]->GetEntry(elist->GetEntry(l));
			for(UInt_t i=0; i<bosonSize; i++){
				// D2 & Mass window 
				D2cut = pt > 2500 ? fD2[i]->Eval(2500) : fD2[i]->Eval(pt);
				MUcut = pt > 2500 ? fMU[i]->Eval(2500) : fMU[i]->Eval(pt);
				MLcut = pt > 2500 ? fML[i]->Eval(2500) : fML[i]->Eval(pt);
				hbef[i][j]->Fill( pt );
				if( d2<D2cut && m>MLcut && m<MUcut ){
					haft[i][j]->Fill( pt );
				}
			}
			if( l%(int(nEvent/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEvent, float(float(l)/nEvent)*100) << std::endl;
		}

	}// chSize

	// Efficiency calculation
	for(UInt_t i=0; i<bosonSize; i++){
		for(UInt_t j=0; j<chsSize; j++){
			setEfficiency(hbef[i][j], haft[i][j], heff[i][j], "rejection");
		}
	}


	// Plot
	TCanvas *cv[bosonSize];
	TLegend *leg[bosonSize];
	TLatex *lat[bosonSize];
	TLatex *latsub[bosonSize];
	for(UInt_t i=0; i<bosonSize; i++){
		cv[i] = new TCanvas();
		leg[i]= new TLegend(0.7,0.66,0.88,0.92);
		leg[i]->SetBorderSize(0);
		lat[i] = new TLatex(250, 100, "#it{ATLAS} #bf{Simulation Work in progress}");
		latsub[i] = new TLatex(400, 90, "#bf{"+boson.at(i)+" tagger}");
		cv[i]->cd();
		for(UInt_t j=0; j<chsSize; j++){
			if( j==0 ) heff[i][j]->Draw();
			else heff[i][j]->Draw("same");
			heff[i][j]->SetStats(0);
			heff[i][j]->GetYaxis()->SetRangeUser(0,110);
			heff[i][j]->GetYaxis()->SetTitle("Background rejection (1/eff)");
			heff[i][j]->GetXaxis()->SetTitle("fatjet pt [GeV]");
			//heff[i][j]->SetMarkerColor(colorList.at(j));
			heff[i][j]->SetMarkerColor(kBlack);
			heff[i][j]->SetMarkerSize(1);
			heff[i][j]->SetMarkerStyle(8);
			//heff[i][j]->SetLineColor(colorList.at(j));
			heff[i][j]->SetLineColor(kBlack);
			heff[i][j]->SetLineWidth(2);
			leg[i]->AddEntry(heff[i][j], chs.at(j), "lep");
		}
		leg[i]->Draw();
		lat[i]->Draw();
		latsub[i]->Draw();

		cv[i]->Print("../output/bkgReje." + boson.at(i) + "." + nameOption + "." + wpName + ".eps");
	}

}
