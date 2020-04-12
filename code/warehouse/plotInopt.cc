#include<iostream>
#include<fstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TCut.h"
#include "TEventList.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"

#include "/nfs_home/home1/myamatan/RootUtils/AtlasStyle.C"
#include "func.hh"
#include "sigList.hh"

TString filenameSig = "../input/merged/signal"+tagger+"/signal"+tagger+"."+ch+".root";
TString filenameBkg = "../input/merged/allMC/allMC."+ch+".root";
	
float sum = 0.0;

void plotInopt(){

	TFile *fS = new TFile(filenameSig);
	TTree *trS = (TTree*)fS->Get("Nominal");
	TFile *fB= new TFile(filenameBkg);
	TTree *trB = (TTree*)fB->Get("Nominal");

	// Fill & Pt flat
	TH1F *hb= new TH1F("hb", "hb", nPtbins, ptRange);
	TH1F *hs[NmassPoints];
	TH1F *hpt[NmassPoints];
	for(UInt_t l=0; l<NmassPoints; l++){
		TString nameID = TString(std::to_string(l));
		hs[l] = new TH1F(nameID+"hs", nameID+"hs", nPtbins, ptRange);
		hpt[l] = new TH1F(nameID+"hpt", nameID+"hpt", nPtbins, ptRange);
	}
	TH1F *hptAll= new TH1F("hptAll", "hptAll", nPtbins, ptRange);

	/////////////////////////////
	// Background
	/////////////////////////////
	std::cout << "Background filling ....." << std::endl;
	SetBranch(trB, ch);
	// Filling
	TCut precut = returnTCut(preSele);
        trB->Draw(">>elist",precut);
        TEventList *elistB = (TEventList*)gROOT->FindObject("elist");
	long nEventB = elistB->GetN(); //long nEventB = trB->GetEntries();
	for(long l=0; l<nEventB; l++){
        	trB->GetEntry(elistB->GetEntry(l)); //trB->GetEntry(l);
		hb->Fill(pt, weight*lumi/40.);
		if( l%(int(nEventB/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEventB, float(float(l)/nEventB)*100) << std::endl;
	}

	/////////////////////////////
	// Signal
	/////////////////////////////
	std::cout << "Signal filling ....." << std::endl;
	std::map<int, int> sigMap;
	UInt_t sigCount = 0;
	UInt_t sigArr[10];
	SetBranch(trS, ch);
	precut = returnTCut(preSele);
        trS->Draw(">>elist",precut);
        TEventList *elistS = (TEventList*)gROOT->FindObject("elist");
	long nEventS = elistS->GetN(); //long nEventS = trS->GetEntries();
	// For pt flat
	for(long l=0; l<nEventS; l++){
        	trS->GetEntry(elistS->GetEntry(l)); //trS->GetEntry(l);
		float weightS = weight;
		if( sigWeight.find(mc) != sigWeight.end() ) weightS = weight * lumi * sigWeight[mc];
		if( sigMap.find(mc) == sigMap.end() ){
			sigMap[mc] = sigCount; 
			sigArr[sigCount] = mc;
			sigCount++;
		}
		hpt[sigMap[mc]]->Fill(pt, weightS);
		hptAll->Fill(pt, weightS);
		sum += weightS;
	}
	// Filling
	for(long l=0; l<nEventS; l++){
        	trS->GetEntry(elistS->GetEntry(l)); //trS->GetEntry(l);
		float weightS = weight;
		if( sigWeight.find(mc) != sigWeight.end() ) weightS = weight * lumi * sigWeight[mc];
		if( isPtFlat ) weightS *= ptFlat(hptAll, pt);
		hs[sigMap[mc]]->Fill(pt, weightS);
		if( l%(int(nEventS/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEventS, float(float(l)/nEventS)*100) << std::endl;
	}

	std::cout << "Plot .............." << std::endl;
	TCanvas *c1 = new TCanvas();
	TLegend *leg= new TLegend(0.6,0.6,0.88,0.88);
	leg->SetBorderSize(0);

	c1->cd();
	hb->SetStats(0);
	hb->SetLineColor(kBlack);
	hb->SetLineWidth(2);
	hb->SetFillStyle(3002);
	hb->SetFillColor(kBlack);

	THStack *hstack = new THStack();
	for(UInt_t m=0; m<sigCount; m++){
		hs[m]->SetFillColor(colorList.at(m));
		hs[m]->SetFillStyle(3002);
		hstack->Add(hs[m]);
		leg->AddEntry(hs[m], TString(sigName[sigArr[m]])+Form(": %0.2f", hs[m]->Integral()), "f");
	}
	//hb->Draw("HIST");	
	hstack->Draw("HIST");
	hstack->SetMaximum(sum/hs[0]->GetNbinsX()*1.6); // 1.6 // 4.6
	hstack->GetXaxis()->SetTitle("fatjet pt [GeV]");
	leg->Draw();

	
		
	
}
