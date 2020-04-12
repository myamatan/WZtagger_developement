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

//std::vector<TString>  boson = {"Z","W"};
std::vector<TString>  boson = {"Z","W"};
std::vector<TString>  chs = {"1Lep"};
//std::vector<TString>  chs = {"0Lep","1Lep","2Lep"};
const UInt_t bosonSize = boson.size();
const UInt_t chsSize = chs.size();

float mcSum=0.;

void eff(){

	SetAtlasStyle();

	TH1F *hbef[bosonSize][chsSize];
	TH1F *haft[bosonSize][chsSize];
	TH1F *heff[bosonSize][chsSize];

	TF1 *fD2;
	TF1 *fML;
	TF1 *fMU;

	TFile *f[bosonSize][chsSize];
	TTree *tr[bosonSize][chsSize];

	//TString wpName = "MaxVVSemi_v5";
	TString wpName = "Fix50VVJJ_LCT";

	TString wpD2Name = wpName;
	TString wpMUName = wpName;
	TString wpMLName = wpName;

	// Fill
	for(UInt_t i=0; i<bosonSize; i++){
		for(UInt_t j=0; j<chsSize; j++){


			//TString fname = "../input/merged_new/signal"+ boson.at(i) + "/signal" + boson.at(i) + "." + chs.at(j) + ".root";
			//TString fname = "../input/merged_tcc_new/signal"+ boson.at(i) + "/signal" + boson.at(i) + "." + chs.at(j) + ".root";
			TString fname = "../input/merged_lct_new/signal"+ boson.at(i) + "/signal" + boson.at(i) + "." + chs.at(j) + ".root";
			//TString fname = "../input/merged/signal"+ boson.at(i) + "/signal" + boson.at(i) + "." + chs.at(j) + ".root";

			std::cout << "++++++ " << fname << " ++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

			f[i][j] = new TFile(fname);
			tr[i][j] = (TTree*)f[i][j]->Get("Nominal");

			hbef[i][j]= new TH1F(Form("hb%i-%i",i,j),Form("hb%i-%i",i,j), nPtbins, ptRange);
			haft[i][j]= new TH1F(Form("ha%i-%i",i,j),Form("ha%i-%i",i,j), nPtbins, ptRange);
			heff[i][j]= new TH1F(Form("heff%i-%i",i,j),Form("heff%i-%i",i,j), nPtbins, ptRange);
			
			// Smoothed cut parameters
			//std::string wpfile = std::string("../output/WPs/SmoothedContained"+boson.at(i)+"Tagger_AntiKt10TrackCaloClusterTrimmed_"+chs.at(0)+"_"+sigType+"_"+nameOption+".txt");
			//std::string wpfile = std::string("../output/WPs/Smoothed"+boson.at(i)+wps[wpName]);

			// D2
			std::string wpfile = std::string("../output/WPs/SmoothedContained"+boson.at(i)+wps[wpD2Name]);
			std::vector<TString> cutWPs = retrieveWPs( wpfile );
			fD2 = new TF1("D2", cutWPs.at(2), 200, 2500);
			cutWPs.clear();
			// MU
			wpfile = std::string("../output/WPs/SmoothedContained"+boson.at(i)+wps[wpMUName]);
			cutWPs = retrieveWPs( wpfile );
			fMU = new TF1("MU", cutWPs.at(1), 200, 2500);
			cutWPs.clear();
			// ML
			wpfile = std::string("../output/WPs/SmoothedContained"+boson.at(i)+wps[wpMLName]);
			cutWPs = retrieveWPs( wpfile );
			fML = new TF1("ML", cutWPs.at(0), 200, 2500);
			cutWPs.clear();

			SetBranch(tr[i][j], chs.at(j));
			
			float D2cut=0.;
			float MUcut=0.;
			float MLcut=0.;
			TCut precut = returnTCut(preSele, chs.at(j));
        		tr[i][j]->Draw(">>elist",precut);
        		TEventList *elist = (TEventList*)gROOT->FindObject("elist");
        		long nEvent = elist->GetN();
			for(long l=0; l<nEvent; l++){
				tr[i][j]->GetEntry(elist->GetEntry(l));
				// D2 & Mass window
				D2cut = pt > 2500 ? fD2->Eval(2500) : fD2->Eval(pt);
				MUcut = pt > 2500 ? fMU->Eval(2500) : fMU->Eval(pt);
				MLcut = pt > 2500 ? fML->Eval(2500) : fML->Eval(pt);
				hbef[i][j]->Fill( pt );
				if( d2<D2cut && m>MLcut && m<MUcut ){
					haft[i][j]->Fill( pt );
				}
				if( l%(int(nEvent/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEvent, float(float(l)/nEvent)*100) << std::endl;
			}
			delete fD2;
			delete fMU;
			delete fML;

		}// chsSize
	}// bosonSize

	// Efficiency calculation
	for(UInt_t i=0; i<bosonSize; i++){
		for(UInt_t j=0; j<chsSize; j++){
			setEfficiency(hbef[i][j], haft[i][j], heff[i][j]);
		}
	}


	// Plot
	TCanvas *cv[bosonSize];
	TLegend *leg[bosonSize];
	TLatex *lat[bosonSize];
	TLatex *latsub[bosonSize];
	for(UInt_t i=0; i<bosonSize; i++){
		//cv[i] = new TCanvas(Form("%i",i),Form("%i",i),200,10,600,500);
		cv[i] = new TCanvas(Form("%i",i),Form("%i",i));
		leg[i]= new TLegend(0.7,0.2,0.88,0.4);
		leg[i]->SetBorderSize(0);
		lat[i] = new TLatex(300, 1.1, "#it{#bf{ATLAS}} #bf{#bf{Simulation Work in progress}}");
		latsub[i] = new TLatex(350, 0.95, "#bf{"+boson.at(i)+" tagger}");
		cv[i]->cd();
		for(UInt_t j=0; j<chsSize; j++){
			if( j==0 ) heff[i][j]->Draw();
			else heff[i][j]->Draw("same");
			heff[i][j]->SetStats(0);
			heff[i][j]->GetYaxis()->SetRangeUser(0,1.2);
			heff[i][j]->GetYaxis()->SetTitle("Signal efficiency");
			heff[i][j]->GetXaxis()->SetTitle("Jet pt [GeV]");
			//heff[i][j]->SetMarkerColor(colorList.at(j));
			heff[i][j]->SetMarkerColor(kBlue);
			heff[i][j]->SetMarkerSize(1);
			heff[i][j]->SetMarkerStyle(8);
			//heff[i][j]->SetLineColor(colorList.at(j));
			heff[i][j]->SetLineColor(kBlue);
			heff[i][j]->SetLineWidth(2);
			leg[i]->AddEntry(heff[i][j], chs.at(j), "lep");
		}
		leg[i]->Draw();
		lat[i]->Draw();
		latsub[i]->Draw();

		//cv[i]->Print("../output/signalEff." + boson.at(i) + "." + nameOption + "." + wpName + ".eps");
		cv[i]->Print("../output/signalEff." + boson.at(i) + "." + nameOption + "." + wpD2Name + "." + wpMUName + "." + wpMLName + ".eps");
		//cv[i]->Print("../output/signalEff." + boson.at(i) + "." + nameOption + ".commonWP.png");
	}

}
