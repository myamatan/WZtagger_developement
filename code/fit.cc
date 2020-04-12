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
#include "TLegend.h"

#include "/nfs_home/home1/myamatan/RootUtils/AtlasStyle.C"
#include "sigList.hh"
#include "func.hh"

TString filename = "../input/merged_new/signal"+tagger+"/signal"+tagger+"."+ch+".root";

float d2cut[nPtbins]={0.};
float mUcut[nPtbins]={0.};
float mLcut[nPtbins]={0.};
float pTbins[nPtbins]={0.};
float pTbinsBand[nPtbins]={0.};
float yError[nPtbins]={0.};

TH1F *h_original;

float markerSize = 1.2;

void fit(){
	
	SetAtlasStyle();

	TFile *f = new TFile(filename);
	TTree *tr = (TTree*)f->Get("Nominal");

	TLatex *latD2 = new TLatex(400, 2.7, "#it{#bf{ATLAS}} #bf{#bf{Simulation Internal}}");
	//TLatex *latD2sub = new TLatex(400, 2.5, "#bf{"+ch+"/"+tagger+" tagger}");
	TLatex *latD2sub = new TLatex(400, 2.5, "#scale[0.8]{#bf{#bf{#sqrt{s}=13 TeV, "+tagger+"-Tagger}}}");
	//TLatex *latMa = new TLatex(400, 140, "#it{ATLAS} #bf{Simulation Preliminary}");
	TLatex *latMa = new TLatex(400, 140, "#it{#bf{ATLAS}} #bf{#bf{Simulation Internal}}");
	TLatex *latMasub = new TLatex(400, 134, "#scale[0.8]{#bf{#bf{#sqrt{s}=13 TeV, "+tagger+"-Tagger}}}");

	// Optiized value output
	string line;
	ifstream ifs(std::string("../output/cutValue.") + ch.Data() + "." + tagger + "." + sigType + "." + nameOption + std::string(".txt"));
	//ifstream ifs(std::string("../output/cutValue.") + ch.Data() + "." + tagger + "." + sigType + ".20181210_reOpt" + std::string(".txt"));

	// Connectin
	std::vector<TString> cutWPs;
	std::string wpfile;
	TF1 *fD2C;
	TF1 *fMLC;
	TF1 *fMUC;
	wpfile = std::string("../output/WPs/SmoothedContained"+tagger+wps["MaxVVJJ"]);
	cutWPs = retrieveWPs( wpfile );
	fD2C = new TF1(Form("D2_comb"), cutWPs.at(2), 200, func_xMax);
	fMUC = new TF1(Form("MU_comb"), cutWPs.at(1), 200, func_xMax);
	fMLC = new TF1(Form("ML_comb"), cutWPs.at(0), 200, func_xMax);

	// Smoothed cut parameters
	ofstream outf(std::string("../output/WPs/SmoothedContained")+tagger+"Tagger_AntiKt10TrackCaloClusterTrimmed_"+ch.Data()+"_"+sigType+"_"+nameOption+std::string(".txt"));

	int count=0;
	int countall=0;
    	while (getline(ifs, line, ','))
	{
		if( count==4 ){
			mLcut[countall] = stof(line);
			count=0;
			countall++;
		}
		else if( count==3 ){
			mUcut[countall] = stof(line);
			count++;
		}
		else if( count==2 ){
			d2cut[countall] = stof(line);
			count++;
		}
		else if( count==1 ){
			pTbinsBand[countall] = int(stof(line)/2.);
			yError[countall] = 0.0;
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
		std::cout << pTbins[i] << "," << d2cut[i] << "," << mUcut[i] << "," << mLcut[i] << std::endl;
	}

	TGraph *gD2 = new TGraphErrors(nPtbins, pTbins, d2cut, pTbinsBand, yError);
	TGraph *gMU = new TGraphErrors(nPtbins, pTbins, mUcut, pTbinsBand, yError);
	TGraph *gML = new TGraphErrors(nPtbins, pTbins, mLcut, pTbinsBand, yError);
	TF1 *fD2 = new TF1("fD2", "[0] + [1]*x + [2]*x*x + [3]*x*x*x", 200, func_xMax);
	TF1 *fMU = new TF1("fMU", "sqrt( pow([0]/x+[1],2) + pow([2]*x+[3],2) )", 200, func_xMax);
	TF1 *fML = new TF1("fML", "sqrt( pow([0]/x-[1],2) + pow([2]*x+[3],2) )", 200, func_xMax);

	fD2->SetParameter(0, 0.8);
        fD2->SetParameter(1, 2E-2);
        fD2->SetParameter(2, 1.5E-12);
        fD2->SetParameter(3, 8E-12);
        //fD2->SetParLimits(1, 0,1E+3);

	if( tagger=="W" ){
		// W-tagger
        	fMU->SetParameter(0, 1800);
        	fMU->SetParameter(1, 93); 
        	fMU->SetParameter(2, 4E-1);
        	fMU->SetParameter(3, -1E+1);
        	fMU->SetParLimits(0, 1E+2,1E+8);
	}else{
		// Z-tagger
        	fMU->SetParameter(0, 2.84369E+3);
        	fMU->SetParameter(1, 1.03296E+2); 
        	fMU->SetParameter(2, 3.76696E-2);
        	fMU->SetParameter(3, -2.46148E-1);
        	fMU->SetParLimits(0, 1E+2,1E+8);
	}


        fML->SetParameter(0, 13E+3);
        fML->SetParameter(1, 65);
        fML->SetParameter(2, -1.8E-2);
        fML->SetParameter(3, 68);

	gD2->Fit(fD2);
	gMU->Fit(fMU);
	gML->Fit(fML);

	// Reference WP
	TF1 *fD2Max;
	TF1 *fMLMax;
	TF1 *fMUMax;
	// MaxSignificance
	//wpfile = std::string("../output/WPs/SmoothedContained"+tagger+wps["MaxVVJJ_3D"]);
	wpfile = std::string("../output/WPs/SmoothedContained"+tagger+wps["MaxVVJJ"]);
	cutWPs = retrieveWPs( wpfile );
	fD2Max = new TF1("D2Max", cutWPs.at(2), 200, func_xMax);
	fMUMax = new TF1("MUMax", cutWPs.at(1), 200, func_xMax);
	fMLMax = new TF1("MLMax", cutWPs.at(0), 200, func_xMax);
	//cutWPs.clear();

	SetBranch(tr, ch);	
	
	// D2 
	TH2F *hD2;
	// Mass Upper
	TH2F *hMa;

	// Histo initialization (x:VV mass, y:fatjet pt)
	hD2 = new TH2F("D2", "D2", nPtbinsDraw, ptRangeDraw, 20, 0, 5);
	hMa = new TH2F("Ma", "Ma", nPtbinsDraw, ptRangeDraw, 20, 50, 150);
	//hD2 = new TH2F("D2", "D2", 80, 200, 2500, 80, 0, 5);
	//hMa = new TH2F("Ma", "Ma", 80, 200, 2500, 40, 50, 150);
	hD2->SetStats(0);
	hMa->SetStats(0);

	// Filling
	float D2cut=0.;
	float MUcut=0.;
	float MLcut=0.;
	TCut precut = returnTCut(preSele);
        tr->Draw(">>elist",precut);
        TEventList *elist = (TEventList*)gROOT->FindObject("elist");
        long nEvent = elist->GetN();//long nEvent = tr->GetEntries();
	for(long l=0; l<nEvent; l++){
		tr->GetEntry(elist->GetEntry(l));//tr->GetEntry(l);
		float weightC = 1;
		//if( sigWeight.find(dsid) != sigWeight.end() ) weightC = weight * lumi * sigWeight[dsid];
		// D2 
		D2cut = pt > 2500 ? fD2->Eval(2500) : fD2->Eval(pt);
		hD2->Fill(pt, d2, weightC);
		//hD2->Fill(pt, d2);
		// Mass cut
		MUcut = pt > 2500 ? fMU->Eval(2500) : fMU->Eval(pt);
		MLcut = pt > 2500 ? fML->Eval(2500) : fML->Eval(pt);
		hMa->Fill(pt, m, weightC);
		//hMa->Fill(pt, m);
		if( l%100000==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEvent, float(float(l)/nEvent)*100) << std::endl;
	}


	// Print WP and save to txt file
	//TString MLform = Form("min( sqrt( pow((%e)/x-(%e),2) + pow((%e)*x+(%e),2)),  "+cutWPs.at(0)+")", fML->GetParameter(0), fML->GetParameter(1), fML->GetParameter(2), fML->GetParameter(3));
	//TString MUform = Form("max( sqrt( pow((%e)/x+(%e),2) + pow((%e)*x+(%e),2)), "+cutWPs.at(1)+")", fMU->GetParameter(0), fMU->GetParameter(1), fMU->GetParameter(2), fMU->GetParameter(3));
	//TString D2form = Form("max( (%e)+(%e)*x+(%e)*pow(x,2)+(%e)*pow(x,3), "+ cutWPs.at(2)+")", fD2->GetParameter(0), fD2->GetParameter(1), fD2->GetParameter(2), fD2->GetParameter(3));
	TString MLform = Form("sqrt( pow((%e)/x-(%e),2) + pow((%e)*x+(%e),2))", fML->GetParameter(0), fML->GetParameter(1), fML->GetParameter(2), fML->GetParameter(3));
	TString MUform = Form("sqrt( pow((%e)/x+(%e),2) + pow((%e)*x+(%e),2))", fMU->GetParameter(0), fMU->GetParameter(1), fMU->GetParameter(2), fMU->GetParameter(3));
	TString D2form = Form("(%e)+(%e)*x+(%e)*pow(x,2)+(%e)*pow(x,3)", fD2->GetParameter(0), fD2->GetParameter(1), fD2->GetParameter(2), fD2->GetParameter(3));
	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	std::cout << "File Name: SmoothedContained" << tagger << "Tagger_AntiKt10TrackCaloClusterTrimmed_" << ch << sigType << "_" << nameOption << std::endl;
	std::cout << "DecorationName: Smoothe" << tagger << "Contained50" << std::endl;
	std::cout << "MassCutLow: " << MLform << std::endl;
	std::cout << "MassCutHigh: " << MUform << std::endl;
	std::cout << "D2Cut: " << D2form << std::endl;
	outf << "File Name: SmoothedContained" << tagger << "Tagger_AntiKt10TrackCaloClusterTrimmed_" << ch << sigType << "_" << nameOption << std::endl;
	outf << "DecorationName: Smoothe" << tagger << "Contained50" << std::endl;
	outf << "MassCutLow: " << MLform << std::endl;
	outf << "MassCutHigh: " << MUform << std::endl;
	outf << "D2Cut: " << D2form << std::endl;
	outf.close();

	// Final WP 
	TF1 *fD2F = new TF1("fD2F", D2form, 200, func_xMax);
	TF1 *fMUF = new TF1("fMUF", MUform, 200, func_xMax);
	TF1 *fMLF = new TF1("fMLF", MLform, 200, func_xMax);

	TCanvas *c1 = new TCanvas("c1","cav1",200,10,600,500);
	TCanvas *c2 = new TCanvas("c2","cav2",200,10,600,500);

	TLegend *leg1 = new TLegend(0.2, 0.65,0.75,0.8);
	TLegend *leg2 = new TLegend(0.2, 0.65,0.75,0.8);
	leg1->SetBorderSize(0);
	leg2->SetBorderSize(0);
	leg1->SetFillStyle(0);
	leg2->SetFillStyle(0);
	

	if( !isPlotMarker ) markerSize = 0.0;

	gStyle->SetPalette(kInvertedDarkBodyRadiator);
	//gStyle->SetPalette(kAquamarine);
	
	c1->cd();

	hD2->GetXaxis()->SetRangeUser(200,func_xMax);
	hD2->GetXaxis()->SetTitle("Jet pt [GeV]");
	hD2->GetXaxis()->SetTitleOffset(1.2);
	hD2->GetYaxis()->SetRangeUser(0,3);
	hD2->GetYaxis()->SetTitle("Jet D_{2}");
	hD2->GetYaxis()->SetTitleOffset(1.2);
	hD2->SetMarkerSize(0.0);
	hD2->SetMarkerColor(kWhite);
	//hD2->SetMarkerColor(kOrange-3);
	hD2->Draw("");

	gD2->SetMarkerSize(markerSize);//1.2
	gD2->SetMarkerStyle(8);
	gD2->SetMarkerColor(kAzure+2);
	gD2->SetLineColor(kAzure+2);
	gD2->SetLineWidth(2);
	gStyle->SetEndErrorSize(3);
	gD2->Draw("P same");
	latD2->Draw("same");
	latD2sub->Draw("same");
	fD2F->SetLineColor(kAzure+2);
	fD2F->Draw("same");

	fD2Max->SetLineColor(kBlack);
	fD2Max->SetLineStyle(2);
	fD2Max->SetLineWidth(3);
	fD2Max->Draw("same");
	
	//leg1->AddEntry(gD2, "Combined Max significance", "lp");
	leg1->AddEntry(gD2, "VVSemileptonic Max significance", "lp");
	//leg1->AddEntry(fD2Max, "VVJJ 3D Max significance", "l");
	leg1->AddEntry(fD2Max, "VVJJ Max significance", "l");
	leg1->Draw("same");

	c2->cd();
	hMa->GetXaxis()->SetRangeUser(200,func_xMax);
	hMa->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
	hMa->GetXaxis()->SetTitleOffset(1.2);
	hMa->GetYaxis()->SetRangeUser(50,150);
	hMa->GetYaxis()->SetTitle("Jet mass [GeV]");
	hMa->GetYaxis()->SetTitleOffset(1.2);
	hMa->SetMarkerSize(0.);
	hMa->SetMarkerColor(kWhite);
	hMa->Draw("");

	gMU->SetMarkerSize(markerSize);//1.2
	gMU->SetMarkerStyle(8);
	gMU->SetMarkerColor(kGreen+2);
	gMU->SetLineWidth(2);
	gStyle->SetEndErrorSize(3);
	gMU->SetLineColor(kGreen+2);
	gMU->SetLineWidth(2);
	gMU->Draw("P same");
	fMUF->SetLineColor(kGreen+2);
	fMUF->Draw("same");

	gML->SetMarkerSize(markerSize);//1.2
	gML->SetMarkerStyle(8);
	gML->SetMarkerColor(kGreen+2);
	gML->SetLineWidth(2);
	gStyle->SetEndErrorSize(3);
	gML->SetLineColor(kGreen+2);
	gML->SetLineWidth(2);
	gML->Draw("P same");
	fMLF->SetLineColor(kGreen+2);
	fMLF->Draw("same");

	latMa->Draw("same");
	latMasub->Draw("same");

	fMUMax->SetLineColor(kBlack);
	fMUMax->SetLineStyle(2);
	fMUMax->SetLineWidth(3);
	fMUMax->Draw("same");

	fMLMax->SetLineColor(kBlack);
	fMLMax->SetLineStyle(2);
	fMLMax->SetLineWidth(3);
	fMLMax->Draw("same");

	//leg2->AddEntry(gMU, "Combined Max significance", "lp");
	leg2->AddEntry(gMU, "VVSemileptonic Max significance", "lp");
	//leg2->AddEntry(fMLMax, "VVJJ 3D Max significance", "l");
	leg2->AddEntry(fMLMax, "VVJJ Max significance", "l");
	leg2->Draw("same");

	nameOption += ".fullrange.oldVVJJ";
	c1->Print("../output/D2.Pt."+ch+"."+tagger+"."+sigType+"."+nameOption+".eps");
        c2->Print("../output/Mass.Pt."+ch+"."+tagger+"."+sigType+"."+nameOption+".eps");

}
