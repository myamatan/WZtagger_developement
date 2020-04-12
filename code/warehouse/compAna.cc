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

TString filenameSig = "../input/merged_Max/signal"+tagger+"/signal"+tagger+"."+ch+".root";
TString filenameBkg = "../input/merged_Max/allMC/allMC."+ch+".root";

TString outNameOpt = "wVmass";
float xMin = 0;
float xMax = 5;
int nBin = 20;
float yScaleMin = 0.7;   // 0.7
float yScaleMax = 1.2;   // 1.6
float yRatioMin = 0.;   // 0.7
float yRatioMax = 2.;  // 1.3 

int ID = 0;
TString xName = "D2";
TString h0Name = "Signal-all";
TString h1Name = "Signal-2tag";
TString pad2Name = "Ratio";
std::map<int, std::vector<TString>> nameMap = {
		{0, {"Background","Signal"}},
		{1, {"Signal-all","Signal-2tag"}},
		{2, {"Background-all","Background-2tag"}},
		{3, {"Background-2tag","Signal-2tag"}}
};

void compAna(){

	if( xName=="D2" ){
		xMin = 0;
		xMax = 5;
		nBin = 20;
		yScaleMax = 1.2;
	}
	else if( xName=="m" ){
		xMin = 50;
		xMax = 150;
		nBin = 20;
		yScaleMax = 1.7;
	}

	h0Name = nameMap[ID].at(0);
	h1Name = nameMap[ID].at(1);
	
	TFile *fS = new TFile(filenameSig);
	TTree *trS = (TTree*)fS->Get("Nominal");
	TFile *fB= new TFile(filenameBkg);
	TTree *trB = (TTree*)fB->Get("Nominal");
	
	TH1F *h0 = new TH1F("h0","h0",nBin,xMin,xMax);
	TH1F *h1 = new TH1F("h1","h1",nBin,xMin,xMax);
	TH1F *hr = new TH1F("hr","hr",nBin,xMin,xMax);


	float vtofill = -99;

	// Background
	std::cout << "Background filling ....." << std::endl;
	SetBranch(trB, ch);
	TCut precut = returnTCut("beforeBosonTag");
        trB->Draw(">>elist",precut);
        TEventList *elistB = (TEventList*)gROOT->FindObject("elist");
	long nEventB = elistB->GetN(); //long nEventB = trB->GetEntries();
	if( h0Name.Index("Background")==-1 && h1Name.Index("Background")==-1 ) nEventB = 0;
	for(long l=0; l<nEventB; l++){
        	trB->GetEntry(elistB->GetEntry(l)); //trB->GetEntry(l);
		// variable to fill
		if( xName=="D2" ) vtofill = d2;
		else if( xName=="m" ) vtofill = m;
		// h0 fill
		if( h0Name.Index("Background")!=-1 && ztag>=2 ){
			if( h0Name.Index("2tag")!=-1 ){
				if( nbtag==2 ) h0->Fill(vtofill, weight*lumi/40.);
			}else{
				h0->Fill(vtofill, weight*lumi/40.);
			}
		}
		// h1 fill
		if( h1Name.Index("Background")!=-1 && ztag>=2 ){
			if( h1Name.Index("2tag")!=-1 ){
				if( nbtag==2 ) h1->Fill(vtofill, weight*lumi/40.);
			}else{
				h1->Fill(vtofill, weight*lumi/40.);
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
	if( h0Name.Index("Signal")==-1 && h1Name.Index("Signal")==-1 ) nEventS = 0;
	for(long l=0; l<nEventS; l++){
        	trS->GetEntry(elistS->GetEntry(l)); //trS->GetEntry(l);
		// variable to fill
		if( xName=="D2" ) vtofill = d2;
		else if( xName=="m" ) vtofill = m;
		float weightS = weight;
		if( sigWeight.find(mc) != sigWeight.end() ){
			weightS = weight * lumi * sigWeight[mc];
			if( sigMap.find(mc) == sigMap.end() ){
				sigMap[mc] = sigCount; 
				sigArr[sigCount] = mc;
				sigCount++;
			}
		}
		// h0 fill
		if( h0Name.Index("Signal")!=-1 && ztag>=2 ){
			if( h0Name.Index("2tag")!=-1 ){
				if( nbtag==2 ) h0->Fill(vtofill, weightS);
			}else{
				h0->Fill(vtofill, weightS);
			}
		}
		// h1 fill
		if( h1Name.Index("Signal")!=-1 && ztag>=2 ){
			if( h1Name.Index("2tag")!=-1 ){
				if( nbtag==2 ) h1->Fill(vtofill, weightS);
			}else{
				h1->Fill(vtofill, weightS);
			}
		}
		if( l%(int(nEventS/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEventS, float(float(l)/nEventS)*100) << std::endl;
	}

	//TCanvas *c1 = new TCanvas("c1","c1",600,600); //cv = new TCanvas("c1","Canvas",200,10,600,480);//600,480
	TCanvas *c1 = new TCanvas("c1","c1",200,10,560,480); //cv = new TCanvas("c1","Canvas",200,10,600,480);//600,480
	TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.18,0.95,0.97);//0.23
        TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.05,0.95,0.305);//0.05
        pad1->SetFillColor(kWhite);
        pad2->SetFillColor(kWhite);
        pad1->Draw();
        pad2->Draw();
        TLegend *leg= new TLegend(0.48,0.6,0.88,0.72); 
        leg->SetBorderSize(0);
	leg->SetFillStyle(0);


	gStyle->SetErrorX(0.0);
	gStyle->SetEndErrorSize(5);

	float h0Inte = h0->Integral();
	float h1Inte = h1->Integral();
	h0->Scale(1./h0->Integral());
	h1->Scale(1./h1->Integral());
	float maxSig=h1->GetBinContent( h1->GetMaximumBin() );
	float minSig=h1->GetBinContent( h1->GetMinimumBin() );

	TLatex *lat;
	TLatex *latsub;
	lat = new TLatex(xMin+(xMax-xMin)*0.4, maxSig*yScaleMax*0.9, "#it{#bf{ATLAS}} #bf{#bf{Simulation Internal}}");
	latsub = new TLatex(xMin+(xMax-xMin)*0.4, maxSig*yScaleMax*0.8, "#bf{Merged Untagged+Tagged}");

	// Ratio calculation
	for(int i=0; i<h0->GetNbinsX(); i++){
		float v0  = h0->GetBinContent(i);
		float er0 = h0->GetBinError(i);
		float v1  = h1->GetBinContent(i);
		float er1 = h1->GetBinError(i);

		if( pad2Name=="Ratio" ){
			float r = v1/v0;
			float rer = sqrt( pow(er1/v0,2) + pow(v1/v0/v0*er0,2) );
			if( isnan(r) || isinf(r) ) r = -100.;
			if( isnan(rer) || isinf(rer) ) rer = 0.;
			hr->SetBinContent(i, r);
			hr->SetBinError(i, rer);
		}
		else if( pad2Name=="Significance" ){
			float sig = sqrt( 2*((v1+v0)*log(1+v1/v0)-v1) );
			float com_par = 0.5 * pow(2*((v1+v0)*log(1+v1/v0)-v1), -1/2.);
                        float s_part = 2 * log(1+v1/v0);
                        float b_part = 2 * (log(1+v1/v0)-v1/v0);
                        float siger = sqrt( pow(com_par*s_part*er1, 2)+pow(com_par*b_part*er0,2) );
			if( isnan(sig) || isinf(sig) ) sig = -100.;
			if( isnan(siger) || isinf(siger) ) siger = 0.;
			hr->SetBinContent(i, sig);
			hr->SetBinError(i, siger);
		}
	}

	pad1->cd();
	//pad1->SetGridy();
	pad1->SetBottomMargin(0.16);
	pad1->SetTopMargin(0.05);
	//h0->GetYaxis()->SetNdivisions(505);
	h0->SetStats(0);
	h0->SetTitle("");
        h0->GetYaxis()->SetTitleSize(0.05);
        h0->GetYaxis()->SetTitleOffset(0.84);
	h0->GetXaxis()->SetRangeUser(xMin,xMax);
	h0->GetYaxis()->SetTitle("normalized #event");
	h0->GetYaxis()->SetRangeUser(minSig*yScaleMin,maxSig*yScaleMax);
	h0->SetLineColor(kBlack);
	h0->SetLineWidth(2);

	h1->SetStats(0);
	h1->SetTitle("");
        h1->GetYaxis()->SetTitleSize(0.14);
	h1->GetXaxis()->SetRangeUser(xMin,xMax);
	h1->GetYaxis()->SetTitle("normalized #event");
	h1->GetYaxis()->SetRangeUser(minSig*yScaleMin,maxSig*yScaleMax);
	h1->SetLineColor(kRed);
	h1->SetLineWidth(2);

	//leg->AddEntry(h0, h0Name+Form("  (%0.1f)",h0Inte), "lp");
	//leg->AddEntry(h1, h1Name+Form("  (%0.1f)",h1Inte), "lp");
	leg->AddEntry(h0, "Background", "lp");
	leg->AddEntry(h1, "Signal", "lp");

	h0->Draw("E0 norm hist");
	h1->Draw("E0 norm hist same");
	leg->Draw();
	lat->Draw();
	latsub->Draw();

	pad2->cd();
	pad2->SetGridy();
	pad2->SetBottomMargin(0.25);
	//pad2->SetTopMargin(0.1);
	hr->SetStats(0);
	hr->SetTitle("");
	hr->GetYaxis()->SetNdivisions(505);
	hr->GetXaxis()->SetTitle(xName);
	if( xName=="m" ) hr->GetXaxis()->SetTitle(xName+" [GeV]");
	hr->GetXaxis()->SetRangeUser(xMin,xMax);
	hr->GetYaxis()->SetTitle("Ratio");
	hr->GetYaxis()->SetRangeUser(yRatioMin, yRatioMax);
	hr->SetLineColor(kRed);
	hr->SetLineWidth(2);
	hr->GetXaxis()->SetTitleSize(0.14);
        hr->GetYaxis()->SetTitleSize(0.14);
        hr->GetXaxis()->SetLabelSize(0.14);
        hr->GetYaxis()->SetLabelSize(0.14);
        hr->GetXaxis()->SetTitleOffset(0.84);
        hr->GetYaxis()->SetTitleOffset(0.27);
	hr->Draw("E0 hist");

	
	c1->Print("../output/sigificance."+ch+"."+tagger+"."+Form("MapID-%d",ID)+"."+xName+"."+outNameOpt+".compAna.eps");
	

}


