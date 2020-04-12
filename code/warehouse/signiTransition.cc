#include<iostream>
#include<fstream>
#include<iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TCut.h"
#include "TEventList.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"

#include "/nfs_home/home1/myamatan/RootUtils/AtlasStyle.C"
#include "func.hh"
#include "sigList.hh"

TString filenameSig = "../input/merged/signal"+tagger+"/signal"+tagger+"."+ch+".root";
TString filenameBkg = "../input/merged/allMC/allMC."+ch+".root";

TString outNameOpt = "wZmass.w2tag";
float d2cut_min = 0.5;
float d2cut_max = 2.0;
float d2cut_d = 0.05;

float xMin = d2cut_min; // 450
float xMax = d2cut_max;// 5450
float yScaleMin = 0.8;   // 0.7
float yScaleMax = 1.6;   // 1.6


void signiTransition(){
	
	TFile *fS = new TFile(filenameSig);
	TTree *trS = (TTree*)fS->Get("Nominal");
	TFile *fB= new TFile(filenameBkg);
	TTree *trB = (TTree*)fB->Get("Nominal");
	
	// D2 range
	int d2_range = UInt_t( (d2cut_max-d2cut_min)/d2cut_d );
	
	// Hist set
	TH1F *hb[d2_range];
	TH1F *hs[d2_range][NmassPoints];
	TH1F *hsig[NmassPoints];
	for(int i=0; i<int(d2_range); i++){
		hb[i]= new TH1F(Form("hb_%i",i), Form("hb_%i",i), nbins, rebin21);
		for(UInt_t l=0; l<NmassPoints; l++){
			hs[i][l]= new TH1F(Form("hs_%i_%i",i,l), Form("hs_%i_%i",i,l), nbins, rebin21);
		}
	}
	for(UInt_t l=0; l<NmassPoints; l++){
		hsig[l]= new TH1F(Form("hsig_%i",l), Form("hsig_%i",l), d2_range, d2cut_min, d2cut_max);
	}

	// Background
	std::cout << "Background filling ....." << std::endl;
	SetBranch(trB, ch);
	TCut precut = returnTCut("beforeBosonTag");
        trB->Draw(">>elist",precut);
        TEventList *elistB = (TEventList*)gROOT->FindObject("elist");
	long nEventB = elistB->GetN(); //long nEventB = trB->GetEntries();
	for(long l=0; l<nEventB; l++){
        	trB->GetEntry(elistB->GetEntry(l)); //trB->GetEntry(l);
		for(int i=0; i<d2_range; i++){
			if( d2 < d2cut_min+d2cut_d*i
	                   && nbtag==2
			   && ztag>=2 
			){
				hb[i]->Fill(vvm, weight*lumi/40.);
				//hb[i]->Fill(vvm, 1.);
			}
		}
		if( l%(int(nEventB/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEventB, float(float(l)/nEventB)*100) << std::endl;
	}

	// Signal
	std::cout << "Signal filling ....." << std::endl;
	std::map<int, int> sigMap;
	UInt_t sigCount = 0;
	UInt_t sigArr[10];
	SetBranch(trS, ch);
	precut = returnTCut("beforeBosonTag");
        trS->Draw(">>elist",precut);
        TEventList *elistS = (TEventList*)gROOT->FindObject("elist");
	long nEventS = elistS->GetN(); //long nEventS = trS->GetEntries();
	for(long l=0; l<nEventS; l++){
        	trS->GetEntry(elistS->GetEntry(l)); //trS->GetEntry(l);
		float weightS = weight;
		if( sigWeight.find(mc) != sigWeight.end() ){
			weightS = weight * lumi * sigWeight[mc];
			if( sigMap.find(mc) == sigMap.end() ){
				sigMap[mc] = sigCount; 
				sigArr[sigCount] = mc;
				sigCount++;
			}
		}
		for(int i=0; i<d2_range; i++){
			if( d2 < d2cut_min+d2cut_d*i
			   && nbtag==2
			   && ztag>=2
			){
				hs[i][sigMap[mc]]->Fill(vvm, weightS);
				//hs[i][sigMap[mc]]->Fill(vvm, 1.);
			}
		}
		if( l%(int(nEventS/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEventS, float(float(l)/nEventS)*100) << std::endl;
	}

	// Significance calculation
	std::cout << "Significance calculation....." << std::endl;
	// D2 & Mass window
	float maxSig=0;
	float minSig=1E+5;
	std::vector<int> massList;
	for(int i=0; i<d2_range; i++){
		for(int j=0; j<int(sigCount); j++){
			float signi = significance(hs[i][j], hs[i][j], hb[i], "mass");
			float error = significance(hs[i][j], hs[i][j], hb[i], "mass", "error");
			hsig[j]->SetBinContent(i, signi);
			hsig[j]->SetBinError(i, error);
			if( i==0 ) massList.push_back(stoi(sigName[sigArr[j]]));
			if( signi > maxSig ) maxSig = signi;
			if( signi < minSig ) minSig = signi;
		}
	}
	std::sort(massList.begin(), massList.end());

	std::map<int, int> orderMap;
	std::map<int, int> op_orderMap;
	for(UInt_t m=0; m<sigCount; m++){
		for(UInt_t l=0; l<massList.size(); l++){
			if( stoi(sigName[sigArr[m]]) == massList.at(l) ){
				orderMap[m] = l; // massList --> mass order
				op_orderMap[l] = m; // mass order --> massList
			}
		}
	}

	// Plottintg .......................................
	TCanvas *c1 = new TCanvas("c1","c1",600,600); //cv = new TCanvas("c1","Canvas",200,10,600,480);//600,480
        TLegend *leg= new TLegend(0.5,0.65,0.88,0.88);                                                    
        leg->SetBorderSize(0);
	leg->SetFillStyle(0);

	TLatex *lat;
	TLatex *latsub;
	lat = new TLatex(0.55, maxSig*yScaleMax*0.9, "#scale[0.65]{#it{#bf{ATLAS}} #bf{#bf{Simulation Internal}}}");
	latsub = new TLatex(0.55, maxSig*yScaleMax*0.8, "#scale[0.65]{#bf{Merged Tagged}}");

	//gStyle->SetErrorX(0.0);
	//gStyle->SetEndErrorSize(5);//5

	//c1->SetGridy();
	//for(auto n = combs.begin(); n != combs.end(); ++n){
	for(int m=0; m<int(sigCount); m++){
		int i = op_orderMap[m];
		hsig[i]->SetStats(0);
		hsig[i]->SetTitle("");
                hsig[i]->GetYaxis()->SetTitleSize(0.04);
		hsig[i]->GetXaxis()->SetTitle("D2 upper cut");
		hsig[i]->GetXaxis()->SetRangeUser(xMin,xMax);
		hsig[i]->GetYaxis()->SetTitle("Significance");
		hsig[i]->GetYaxis()->SetRangeUser(minSig*yScaleMin,maxSig*yScaleMax);
		hsig[i]->SetLineColor(colorListSig.at(i));
		hsig[i]->SetLineWidth(2);
		hsig[i]->SetMarkerColor(colorListSig.at(i));
		hsig[i]->SetMarkerStyle(8);
		hsig[i]->SetMarkerSize(0.8);
		leg->AddEntry(hsig[i], "Signal : "+TString(sigName[sigArr[i]])+" [GeV]", "lp");
		if( m!=0 ) hsig[i]->Draw("E0 same");
		else hsig[i]->Draw("E0");
		std::cout << TString(sigName[sigArr[i]]) << " : ";
		for(int k=1; k<=hsig[i]->GetNbinsX(); k++){
			std::cout << hsig[i]->GetBinContent(k) << " ";
		}
		std::cout << std::endl;
	}
	leg->Draw();
	lat->Draw();
	latsub->Draw();

	c1->Print("../output/sigTrans."+ch+"."+tagger+"."+sigType+"."+nameOption+"."+outNameOpt+".eps");
	

}


