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
#include "TCut.h"
#include "TEventList.h"

#include "/nfs_home/home1/myamatan/RootUtils/AtlasStyle.C"
#include "func.hh"
#include "sigList.hh"

TString filename = "../input/merged_nv/signal"+tagger+"/signal"+tagger+"."+ch+".root";
//TString filename = "../input/merged/allMC/allMC."+ch+".root";

float d2cut[nPtbins]={0.};
float mUcut[nPtbins]={0.};
float mLcut[nPtbins]={0.};
float ktdrcut[nPtbins]={0.};
float pTbins[nPtbins]={0.};

TH1F *h_original;

float markerSize = 1.2;

void fit3D(){
	
	//SetAtlasStyle();

	TFile *f = new TFile(filename);
	TTree *tr = (TTree*)f->Get("Nominal");

	TLatex *latD2 = new TLatex(400, 2.7, "#color[1]{#it{ATLAS} #bf{Simulation Preliminary}}}");
	TLatex *latD2sub = new TLatex(400, 2.5, "#color[1]{#bf{"+ch+"/"+tagger+" tagger}}");
	TLatex *latMa = new TLatex(400, 140, "#color[1]{#it{ATLAS} #bf{Simulation Preliminary}}");
	TLatex *latMasub = new TLatex(400, 132, "#color[1]{#bf{"+ch+"/"+tagger+" tagger}}");
	TLatex *latKtDR = new TLatex(400, 0.9, "#color[1]{#it{ATLAS} #bf{Simulation Preliminary}}");
	TLatex *latKtDRsub = new TLatex(400, 0.8, "#color[1]{#bf{"+ch+"/"+tagger+" tagger}}");


	std::vector<TString> cutWPs;
	std::string wpfile;
	TF1 *fD2C;
	TF1 *fMLC;
	TF1 *fMUC;
	wpfile = std::string("../output/WPs/SmoothedContained"+tagger+wps["MaxVVJJ"]);
	cutWPs = retrieveWPs( wpfile );
	fD2C = new TF1(Form("D2_comb"), cutWPs.at(2), 200, 1100);
	fMUC = new TF1(Form("MU_comb"), cutWPs.at(1), 200, 1100);
	fMLC = new TF1(Form("ML_comb"), cutWPs.at(0), 200, 1100);

	// Smoothed cut parameters output
	ofstream outf(std::string("../output/WPs/SmoothedContained")+tagger+"Tagger_AntiKt10TrackCaloClusterTrimmed_3D_"+ch.Data()+"_"+sigType+"_"+nameOption+std::string(".txt"));

	// Optiized value input
	string line;
	ifstream ifs(std::string("../output/cutValue.3D.") + ch.Data() + "." + tagger + "." + sigType + "." + nameOption + std::string(".txt"));

	int count=0;
	int countall=0;
    	while (getline(ifs, line, ','))
	{
		if( count==4 ){
			ktdrcut[countall] = stof(line);
			count=0;
			countall++;
		}
		else if( count==3 ){
			mLcut[countall] = stof(line);
			count++;
		}
		else if( count==2 ){
			mUcut[countall] = stof(line);
			count++;
		}
		else if( count==1 ){
			d2cut[countall] = stof(line);
			count++;
		}
		else if( count==0 ){
			pTbins[countall] = stof(line);
			count++;
		}
	}
    	ifs.close();

	for(UInt_t i=0; i<nPtbins; i++){
		if( int(i) >= connectBin && isCombineFit ){
			d2cut[i] = fD2C->Eval(pTbins[i]);
			mUcut[i] = fMUC->Eval(pTbins[i]);
			mLcut[i] = fMLC->Eval(pTbins[i]);
		}
		std::cout << pTbins[i] << "," << d2cut[i] << "," << mUcut[i] << "," << mLcut[i] << "," << ktdrcut[i] << std::endl;
	}

	TGraph *gD2 = new TGraph(nPtbins, pTbins, d2cut);
	TGraph *gMU = new TGraph(nPtbins, pTbins, mUcut);
	TGraph *gML = new TGraph(nPtbins, pTbins, mLcut);
	TGraph *gKtDR = new TGraph(nPtbins, pTbins, ktdrcut);
	TF1 *fD2 = new TF1("fD2", "[0] + [1]*x + [2]*x*x + [3]*x*x*x", 200, 1100);
	TF1 *fMU = new TF1("fMU", "sqrt( pow([0]/x+[1],2) + pow([2]*x+[3],2) )", 200, 1100);
	TF1 *fML = new TF1("fML", "sqrt( pow([0]/x-[1],2) + pow([2]*x+[3],2) )", 200, 1100);
	TF1 *fKtDR = new TF1("fKtDR", "[0] + [1]*x + [2]*x*x + [3]*x*x*x", 200, 1100);
	fD2->SetParameter(0, 0.7);
	fD2->SetParameter(1, 2E-4);
	fD2->SetParameter(2, -1.5E-7);
	fD2->SetParameter(3, 8E-11);
	fMU->SetParameter(0, 180);
	fMU->SetParameter(1, 93);
	fMU->SetParameter(2, 2E-2);
	fMU->SetParameter(3, -28);
	fML->SetParameter(0, 13E+3);
	fML->SetParameter(1, 65);
	fML->SetParameter(2, -1.8E-2);
	fML->SetParameter(3, 68);
	fKtDR->SetParameter(0, 0.7);
	fKtDR->SetParameter(1, 2E-4);
	fKtDR->SetParameter(2, -1.5E-7);
	fKtDR->SetParameter(3, 8E-11);
	gD2->Fit(fD2);
	gMU->Fit(fMU);
	gML->Fit(fML);
	gKtDR->Fit(fKtDR);

	// Reference WP
	TF1 *fD2Max;
	TF1 *fMLMax;
	TF1 *fMUMax;
	TF1 *fD2Fix;
	TF1 *fMLFix;
	TF1 *fMUFix;
	// MaxSignificance
	wpfile = std::string("../input/WPs/SmoothedContained"+tagger+"Tagger_AntiKt10TrackCaloClusterTrimmed_MaxSignificance_MC16a_20171201.dat");
	cutWPs = retrieveWPs( wpfile );
	fD2Max = new TF1("D2Max", cutWPs.at(2), 200, 1100);
	fMUMax = new TF1("MUMax", cutWPs.at(1), 200, 1100);
	fMLMax = new TF1("MLMax", cutWPs.at(0), 200, 1100);
	cutWPs.clear();
	// FixedSignalEfficienncy
	wpfile = std::string("../input/WPs/SmoothedContained"+tagger+"Tagger_AntiKt10TrackCaloClusterTrimmed_FixedSignalEfficiency50_MC16a_20171201.dat");
	cutWPs = retrieveWPs( wpfile );
	fD2Fix = new TF1("D2Fix", cutWPs.at(2), 200, 1100);
	fMUFix = new TF1("MUFix", cutWPs.at(1), 200, 1100);
	fMLFix = new TF1("MLFix", cutWPs.at(0), 200, 1100);
	cutWPs.clear();

	SetBranch(tr, ch);	
	
	// D2 
	TH2F *hD2;
	// Mass Upper
	TH2F *hMa;
	// KtDR
	TH2F *hKtDR;

	// Histo initialization (x:VV mass, y:fatjet pt)
	hD2 = new TH2F("D2", "D2", nPtbins, ptRange, 80, 0, 5);
	hMa = new TH2F("Ma", "Ma", nPtbins, ptRange, 40, 50, 150);
	hKtDR = new TH2F("KtDR", "", nPtbins, ptRange, 80, 0, 3);
	hD2->SetStats(0);
	hMa->SetStats(0);
	hKtDR->SetStats(0);

	// Filling
	float D2cut=0.;
	float MUcut=0.;
	float MLcut=0.;
	float KtDRcut=0.;
	TCut precut = returnTCut("beforeBosonTag");
        tr->Draw(">>elist",precut);
        TEventList *elist = (TEventList*)gROOT->FindObject("elist");
        long nEvent = elist->GetN();//long nEvent = tr->GetEntries();
	for(long l=0; l<nEvent; l++){
		tr->GetEntry(elist->GetEntry(l));//tr->GetEntry(l);
		float weightC = 1.0;
		//if( sigWeight[mc]!=0 ) weightC = weight * 40 * sigWeight[mc];
		// D2 
		D2cut = pt > 1100 ? fD2->Eval(1100) : fD2->Eval(pt);
		if( d2 < D2cut &&
		    1
		){
			hD2->Fill(pt, d2, weightC);
		}
		// Mass cut
		MUcut = pt > 1100 ? fMU->Eval(1100) : fMU->Eval(pt);
		MLcut = pt > 1100 ? fML->Eval(1100) : fML->Eval(pt);
		if( m < MUcut &&
		    m > MLcut &&
		    1
		){
			hMa->Fill(pt, m, weightC);
		}
		// KtDR
		KtDRcut = pt > 1100 ? fKtDR->Eval(1100) : fKtDR->Eval(pt);
		if( ktdr < KtDRcut &&
		    1
                ){
			hKtDR->Fill(pt, ktdr, weightC);
		}
		if( l%100000==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEvent, float(float(l)/nEvent)*100) << std::endl;
	}

	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	std::cout << "File Name: SmoothedContained" << tagger << "Tagger_AntiKt10TrackCaloClusterTrimmed_" << ch << sigType << "_" << nameOption << std::endl;
	std::cout << "DecorationName: Smoothe" << tagger << "Contained50" << std::endl;
	std::cout << Form("MassCutLow: x > 1100 ? %e : sqrt( pow((%e)/x-(%e),2) + pow((%e)*x+(%e),2) )", fML->Eval(1100), fML->GetParameter(0), fML->GetParameter(1), fML->GetParameter(2), fML->GetParameter(3)) << std::endl;;
	std::cout << Form("MassCutHigh: x > 1100 ? %e : sqrt( pow((%e)/x+(%e),2) + pow((%e)*x+(%e),2) )", fMU->Eval(1100), fMU->GetParameter(0), fMU->GetParameter(1), fMU->GetParameter(2), fMU->GetParameter(3)) << std::endl;;
	std::cout << Form("D2Cut: x > 1100 ? %e : (%e)+(%e)*x+(%e)*pow(x,2)+(%e)*pow(x,3)", fD2->Eval(1100), fD2->GetParameter(0), fD2->GetParameter(1), fD2->GetParameter(2), fD2->GetParameter(3)) << std::endl;;
	std::cout << Form("KtDRCut: x > 1100 ? %e : (%e)+(%e)*x+(%e)*pow(x,2)+(%e)*pow(x,3)", fKtDR->Eval(1100), fKtDR->GetParameter(0), fKtDR->GetParameter(1), fKtDR->GetParameter(2), fKtDR->GetParameter(3)) << std::endl;;
	outf << "File Name: SmoothedContained" << tagger << "Tagger_AntiKt10TrackCaloClusterTrimmed_" << ch << sigType << "_" << nameOption << std::endl;
	outf << "DecorationName: Smoothe" << tagger << "Contained50" << std::endl;
	outf << Form("MassCutLow: x > 1100 ? %e : sqrt( pow((%e)/x-(%e),2) + pow((%e)*x+(%e),2) )", fML->Eval(1100), fML->GetParameter(0), fML->GetParameter(1), fML->GetParameter(2), fML->GetParameter(3)) << std::endl;;
	outf << Form("MassCutHigh: x > 1100 ? %e : sqrt( pow((%e)/x+(%e),2) + pow((%e)*x+(%e),2) )", fMU->Eval(1100), fMU->GetParameter(0), fMU->GetParameter(1), fMU->GetParameter(2), fMU->GetParameter(3)) << std::endl;;
	outf << Form("D2Cut: x > 1100 ? %e : (%e)+(%e)*x+(%e)*pow(x,2)+(%e)*pow(x,3)", fD2->Eval(1100), fD2->GetParameter(0), fD2->GetParameter(1), fD2->GetParameter(2), fD2->GetParameter(3)) << std::endl;;
	outf << Form("KtDRCut: x > 1100 ? %e : (%e)+(%e)*x+(%e)*pow(x,2)+(%e)*pow(x,3)", fKtDR->Eval(1100), fKtDR->GetParameter(0), fKtDR->GetParameter(1), fKtDR->GetParameter(2), fKtDR->GetParameter(3)) << std::endl;;
	outf.close();

	TCanvas *c1 = new TCanvas();
	TCanvas *c2 = new TCanvas();
	TCanvas *c3 = new TCanvas();

	if( !isPlotMarker ) markerSize = 0.0;

	c1->cd();
	hD2->GetXaxis()->SetRangeUser(200,1100);
	hD2->GetXaxis()->SetTitle("fatJet pt [GeV]");
	hD2->GetYaxis()->SetRangeUser(0,3);
	hD2->GetYaxis()->SetTitle("fatJet D2");
	hD2->Draw("colz");
	gD2->SetMarkerSize(markerSize);//1.2
	gD2->SetMarkerStyle(8);
	gD2->SetMarkerColor(kRed+1);
	gD2->Draw("P same");
	latD2->Draw("same");
	latD2sub->Draw("same");

	fD2Max->SetLineColor(kBlack);
	fD2Max->SetLineStyle(2);
	fD2Max->SetLineWidth(3);
	fD2Max->Draw("same");
	fD2Fix->SetLineColor(kGray);
	fD2Fix->SetLineStyle(2);
	fD2Fix->SetLineWidth(3);
	fD2Fix->Draw("same");

	c2->cd();
	hMa->GetXaxis()->SetRangeUser(200,1100);
	hMa->GetXaxis()->SetTitle("fatJet pt [GeV]");
	hMa->GetYaxis()->SetRangeUser(50,150);
	hMa->GetYaxis()->SetTitle("fatJet mass [GeV]");
	hMa->Draw("colz");

	gMU->SetMarkerSize(markerSize);//1.2
	gMU->SetMarkerStyle(8);
	gMU->SetMarkerColor(kGreen+2);
	fMU->SetLineColor(kGreen+2);
	gMU->Draw("P same");
	fMU->Draw("same");

	gML->SetMarkerSize(markerSize);//1.2
	gML->SetMarkerStyle(8);
	gML->SetMarkerColor(kRed+1);
	gML->Draw("P same");

	latMa->Draw("same");
	latMasub->Draw("same");

	fMUMax->SetLineColor(kBlack);
	fMUMax->SetLineStyle(2);
	fMUMax->SetLineWidth(3);
	fMUMax->Draw("same");
	fMUFix->SetLineColor(kGray);
	fMUFix->SetLineStyle(2);
	fMUFix->SetLineWidth(3);
	fMUFix->Draw("same");

	fMLMax->SetLineColor(kBlack);
	fMLMax->SetLineStyle(2);
	fMLMax->SetLineWidth(3);
	fMLMax->Draw("same");
	fMLFix->SetLineColor(kGray);
	fMLFix->SetLineStyle(2);
	fMLFix->SetLineWidth(3);
	fMLFix->Draw("same");

	c3->cd();
	hKtDR->GetXaxis()->SetRangeUser(200,1100);
	hKtDR->GetXaxis()->SetTitle("fatJet pt [GeV]");
	hKtDR->GetYaxis()->SetRangeUser(0,1.0);
	hKtDR->GetYaxis()->SetTitle("fatJet KtDR");
	hKtDR->Draw("colz");
	gKtDR->SetMarkerSize(markerSize);//1.2
	gKtDR->SetMarkerStyle(8);
	gKtDR->SetMarkerColor(kRed+1);
	gKtDR->Draw("P same");
	latKtDR->Draw("same");
	latKtDRsub->Draw("same");

	c1->Print("../output/D2.Pt.3D."+ch+"."+tagger+"."+sigType+"."+nameOption+".nofit.png");
        c2->Print("../output/Mass.Pt.3D."+ch+"."+tagger+"."+sigType+"."+nameOption+".nofit.png");
        c3->Print("../output/KtDR.Pt.3D."+ch+"."+tagger+"."+sigType+"."+nameOption+".nofit.png");

}
