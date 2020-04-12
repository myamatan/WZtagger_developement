#include<iostream>
#include<fstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TEventList.h"

#include "/nfs_home/home1/myamatan/RootUtils/AtlasStyle.C"
#include "func.hh"
#include "sigList.hh"

TString filenameSig = "../input/merged/signal"+tagger+"/signal"+tagger+"."+ch+".root";
//TString filenameSig = "../input/merged/signal/RS_G_ZZm1000."+ch+".root";
TString filenameBkg = "../input/merged/allMC/allMC."+ch+".root";

float d2cut_min = 0.5;//0.5
float d2cut_max = 0.6;//1.6
float d2cut_d = 0.1;//0.05

float mcutu_min = 100.;//100
float mcutu_max = 140.;//120
float mcutu_d = 1;//2

float mcutl_min = 52.;//60
float mcutl_max = 80.;//80
float mcutl_d = 1;

void dist(){

	TFile *fS = new TFile(filenameSig);
	TTree *trS = (TTree*)fS->Get("Nominal");
	TFile *fB= new TFile(filenameBkg);
	TTree *trB = (TTree*)fB->Get("Nominal");

	UInt_t d2_range = UInt_t( (d2cut_max-d2cut_min)/d2cut_d );
	UInt_t mu_range = UInt_t( (mcutu_max-mcutu_min)/mcutu_d );
	UInt_t ml_range = UInt_t( (mcutl_max-mcutl_min)/mcutl_d );

	TH1F *hb[d2_range][mu_range][ml_range];
	TH1F *hs[d2_range][mu_range][ml_range][NmassPoints][2];
	for(UInt_t i=0; i<d2_range; i++){
		for(UInt_t j=0; j<mu_range; j++){
			for(UInt_t k=0; k<ml_range; k++){
				TString nameID = TString(std::to_string(i) + "-" + std::to_string(j) + "-" + std::to_string(k) + "-");
				hb[i][j][k] = new TH1F(nameID+"hb", nameID+"hb",nbins, rebin21);
				for(UInt_t l=0; l<NmassPoints; l++){
					nameID = TString(std::to_string(i) + "-" + std::to_string(j) + "-" + std::to_string(k) + "-" + std::to_string(l) + "-");
					hs[i][j][k][l][0] = new TH1F(nameID+"hsBe", nameID+"hsBe", nbins, rebin21);
					hs[i][j][k][l][1] = new TH1F(nameID+"hsAf", nameID+"hsAf", nbins, rebin21);
				}
			}
		}
	}

	// Cut function
	TF1 *fD2;
	TF1 *fML;
	TF1 *fMU;
	float D2cut=0.;
	float MUcut=0.;
	float MLcut=0.;
	std::string wpfile = std::string("../output/SmoothedContained"+tagger+"Tagger_AntiKt10TrackCaloClusterTrimmed_"+ch+"_"+sigType+"_"+nameOption+".txt");
	std::vector<TString> cutWPs = retrieveWPs( wpfile );
	fD2 = new TF1("D2", cutWPs.at(2), 200, 2500);
	fMU = new TF1("MU", cutWPs.at(1), 200, 2500);
	fML = new TF1("ML", cutWPs.at(0), 200, 2500);

	
	ofstream outf(std::string("../output/SmoothedContained")+tagger+"Tagger_AntiKt10TrackCaloClusterTrimmed_"+ch.Data()+"_"+sigType+"_"+nameOption+std::string("_ConstantMassWindow.txt"));

	// Background
	std::cout << "Background filling ....." << std::endl;
	SetBranch(trB, ch);
        TCut precut = returnTCut(preSele);
        trB->Draw(">>elist",precut);
        TEventList *elistB = (TEventList*)gROOT->FindObject("elist");
	long nEventB = elistB->GetN(); //long nEventB = trB->GetEntries();
	for(long l=0; l<nEventB; l++){
        	trB->GetEntry(elistB->GetEntry(l)); //trB->GetEntry(l);
		D2cut = pt > 2500 ? fD2->Eval(2500) : fD2->Eval(pt);
		MUcut = pt > 2500 ? fMU->Eval(2500) : fMU->Eval(pt);
		MLcut = pt > 2500 ? fML->Eval(2500) : fML->Eval(pt);
		for(UInt_t i=0; i<d2_range; i++){
			for(UInt_t j=0; j<mu_range; j++){
				for(UInt_t k=0; k<ml_range; k++){
					if( 
					    //d2 < d2cut_min+d2cut_d*i &&
					    d2 < D2cut &&
					    m  < mcutu_max-mcutu_d*j &&
					    m  > mcutl_min+mcutl_d*k
					){
						hb[i][j][k]->Fill(vvm, weight);
					}
				}
			}
		}
		if( l%200000==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEventB, float(float(l)/nEventB)*100) << std::endl;
	}

	// Signal
	std::cout << "Signal filling ....." << std::endl;
	std::map<int, int> sigMap;
	UInt_t sigArr[10];
	UInt_t sigCount = 0;
	SetBranch(trS, ch);
        precut = returnTCut(preSele);
        trS->Draw(">>elist",precut);
        TEventList *elistS = (TEventList*)gROOT->FindObject("elist");
	long nEventS = elistS->GetN(); //long nEventS = trS->GetEntries();
	for(long l=0; l<nEventS; l++){
        	trS->GetEntry(elistS->GetEntry(l)); //trS->GetEntry(l);
		float weightS = 1.0;
		if( sigWeight[mc]!=0 ) weightS = weight * 40 * sigWeight[mc];
		if( sigMap.find(mc) == sigMap.end() ){
			sigMap[mc] = sigCount; 
			sigArr[sigCount] = mc;
			sigCount++;
		}
		D2cut = pt > 2500 ? fD2->Eval(2500) : fD2->Eval(pt);
		MUcut = pt > 2500 ? fMU->Eval(2500) : fMU->Eval(pt);
		MLcut = pt > 2500 ? fML->Eval(2500) : fML->Eval(pt);
		for(UInt_t i=0; i<d2_range; i++){
		  for(UInt_t j=0; j<mu_range; j++){
		    for(UInt_t k=0; k<ml_range; k++){

		      hs[i][j][k][sigMap[mc]][0]->Fill(vvm, weightS);
		      if( 
			  //d2 < d2cut_min+d2cut_d*i &&
		          d2 < D2cut &&
			  m  < mcutu_max-mcutu_d*j &&
			  m  > mcutl_min+mcutl_d*k
		      ){
		        hs[i][j][k][sigMap[mc]][1]->Fill(vvm, weightS);
		      }

		    }
		  }
		}
		if( l%5000==0 ) std::cout << Form("[%li/%li] : %0.0f %%  done .....", l, nEventS, float(float(l)/nEventS)*100) << std::endl;
	}

	// Significance calculation
	std::cout << "Optimized cut searching  ....." << std::endl;
	TH1F *h1Sig[10];
	for(UInt_t i=0; i<10; i++){
		h1Sig[i] = new TH1F(Form("h1Sig%i",i),Form("h1Sig%i",i),d2_range,d2cut_min,d2cut_max);
	}
	TH2F *hSig = new TH2F("Sig","Sig",sigCount,0,sigCount,d2_range,d2cut_min,d2cut_max);
	float d2cut[10]={0.0};
	float mcutu[10]={0.0};
	float mcutl[10]={0.0};
	float maxSigni[10]={0.0};
	float signiSum=0.0;
	float d2cutAll=0.;
	float mcutuAll=0.;
	float mcutlAll=0.;
	float maxSigniAll=0.;
	for(UInt_t i=0; i<d2_range; i++){
		for(UInt_t j=0; j<mu_range; j++){
			for(UInt_t k=0; k<ml_range; k++){

				signiSum = 0.0;
		        	for(UInt_t l=0; l<sigCount; l++){
			        	float perSig = significance(hs[i][j][k][l][0], hs[i][j][k][l][1], hb[i][j][k], sigType);
					signiSum += perSig;
					if( perSig > maxSigni[l] ){
						d2cut[l] = d2cut_min+d2cut_d*i;
						mcutu[l] = mcutu_max-mcutu_d*j;
						mcutl[l] = mcutl_min+mcutl_d*k;
						maxSigni[l] = perSig;
					}
					hSig->SetBinContent(l, i, perSig);
					h1Sig[l]->SetBinContent(i, perSig);
					h1Sig[l]->SetBinError(i, 0);
				}
				if( signiSum > maxSigniAll ){
					d2cutAll = d2cut_min+d2cut_d*i;
					mcutuAll = mcutu_max-mcutu_d*j;
					mcutlAll = mcutl_min+mcutl_d*k;
					maxSigniAll = signiSum;
				}

			}
		}
	}

	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	std::cout << "+++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	std::cout << "Max Significance : " << maxSigniAll << std::endl;
	std::cout << "MC : All" << std::endl;
	std::cout << Form("D2 cut : %0.2f", d2cutAll) << std::endl;
	std::cout << Form("Upper mass cut : %0.2f", mcutuAll) << std::endl;
	std::cout << Form("Lower mass cut : %0.2f", mcutlAll) << std::endl;
	/*
	for(UInt_t l=0; l<sigCount; l++){
		std::cout << "+++++++++++++++++++++++++++++++++++++++++++" << std::endl;
		std::cout << "Max Significance : " << maxSigni[l] << std::endl;
		std::cout << "MC : " << TString(sigName[sigArr[l]]) << std::endl;
		std::cout << Form("D2 cut : %0.2f", d2cut[l]) << std::endl;
		std::cout << Form("Upper mass cut : %0.2f", mcutu[l]) << std::endl;
		std::cout << Form("Lower mass cut : %0.2f", mcutl[l]) << std::endl;
	}
	*/

	TLegend *leg= new TLegend(0.7,0.6,0.88,0.88);
	leg->SetBorderSize(0);
	TCanvas *c1 = new TCanvas();
	//c1->SetLogy();

	for(UInt_t i=0; i<sigCount; i++){
		if( i==0 ) h1Sig[i]->Draw("HIST");
		else h1Sig[i]->Draw("HIST same");
		h1Sig[i]->SetLineColor(colorList.at(i));
		h1Sig[i]->Scale(h1Sig[0]->Integral()/h1Sig[i]->Integral());
		h1Sig[i]->SetLineWidth(2);
		h1Sig[i]->SetMarkerColor(colorList.at(i));
		h1Sig[i]->SetStats(0);
		h1Sig[i]->GetXaxis()->SetTitle("fatjet D2");
		h1Sig[i]->GetXaxis()->SetRangeUser(0.5, 1.5);
		h1Sig[i]->GetYaxis()->SetTitle("Significance");
		h1Sig[i]->GetYaxis()->SetRangeUser(3E-3, 9E-3);
		leg->AddEntry(h1Sig[i], TString(sigName[sigArr[i]]), "lep");
	}
	leg->Draw("same");


	//c1->Print("../output/dist." + ch + ".png");
		
	outf << "File Name: SmoothedContained" << tagger << "Tagger_AntiKt10TrackCaloClusterTrimmed_" << ch << sigType << "_" << nameOption << "_ConstantMassWindow" << std::endl;
	outf << "DecorationName: Smoothe" << tagger << "Contained50" << std::endl;
	outf << Form("MassCutLow: %e", mcutlAll) << std::endl;;
	outf << Form("MassCutHigh: %e", mcutuAll) << std::endl;;
	outf << Form("D2Cut: x > 2500 ? %e:", fD2->Eval(2500))+cutWPs.at(2) << std::endl;;
	outf.close();

	//hSig->GetXaxis()->SetTitle("Signal mass [TeV]");
	//hSig->GetYaxis()->SetTitle("Fatjet D2");
	//hSig->Draw("colz");



}
