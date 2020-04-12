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
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"

#include "/nfs_home/home1/myamatan/RootUtils/AtlasStyle.C"
#include "sigList.hh"
#include "func.hh"

//ch="0Lep";

TString filenameSig = "../input/merged_tcc_new/signal"+tagger+"/signal"+tagger+"."+ch+".root";
TString filenameBkg = "../input/merged_tcc_new/allMC/allMC."+ch+".root";

TString filenameSig_lct = "../input/merged_lct_new/signal"+tagger+"/signal"+tagger+"."+ch+".root";
TString filenameBkg_lct = "../input/merged_lct_new/allMC/allMC."+ch+".root";

TString outNameOpt = "comb";
float xMin = 250; // 450
float xMax = 5450;// 5450
float yScaleMin = 0.7;   // 0.7
float yScaleMax = 1.6;   // 1.6
float yRatioMin = 0.7;   // 0.7
float yRatioMax = 1.3;  // 1.3 

int ratio_deno = 1;

const int Nwps = wps.size();
const int Ncombs = combs.size();
float d2cut = 0;
float mcutu = 0;
float mcutl = 0;

// Compared to LCTopo
void signiComp_LCT(){

	////////////////////////////////////////////////////
	// LCTopo first
	////////////////////////////////////////////////////
	TFile *fS_lct  = new TFile(filenameSig_lct);
	TTree *trS_lct = (TTree*)fS_lct->Get("Nominal");
	TFile *fB_lct  = new TFile(filenameBkg_lct);
	TTree *trB_lct = (TTree*)fB_lct->Get("Nominal");

	TF1 *fD2_lct;
	TF1 *fML_lct;
	TF1 *fMU_lct;
	std::vector<TString> cutWPs_lct = retrieveWPs( std::string("../output/WPs/SmoothedContained"+tagger+wps["Fix50VVJJ_LCT"]) );
	fD2_lct = new TF1(Form("D2_LCT"), cutWPs_lct.at(2), 200, 2500);
	fMU_lct = new TF1(Form("MU_LCT"), cutWPs_lct.at(1), 200, 2500);
	fML_lct = new TF1(Form("ML_LCT"), cutWPs_lct.at(0), 200, 2500);

	TH1F *hb_lct  = new TH1F("hb_lct","hb_lct", nbins, rebin21);
	TH1F *hg_lct  = new TH1F("hg_lct","hg_lct", 50,450,5450);
	TH1F *hgr_lct = new TH1F("hgrlct","hgr_lct", 50,450,5450);
	TH1F *hs_lct[NmassPoints];
	for(UInt_t l=0; l<NmassPoints; l++){
		hs_lct[l]= new TH1F(Form("hs_lct_%i",l), Form("hs_lct_%i",l), nbins, rebin21);
	}

	// Background
	std::cout << "LCT Background filling ....." << std::endl;
	SetBranch(trB_lct, ch);
	TCut precut_lct = returnTCut(preSele);
        trB_lct->Draw(">>elist",precut_lct);
        TEventList *elistB_lct = (TEventList*)gROOT->FindObject("elist");
	long nEventB_lct = elistB_lct->GetN(); //long nEventB = trB->GetEntries();
	for(long l=0; l<nEventB_lct; l++){
        	trB_lct->GetEntry(elistB_lct->GetEntry(l)); //trB->GetEntry(l);
		// D2 &&  Mass window
		d2cut = pt > 2500 ? fD2_lct->Eval(2500) : fD2_lct->Eval(pt);
		mcutu = pt > 2500 ? fMU_lct->Eval(2500) : fMU_lct->Eval(pt);
		mcutl = pt > 2500 ? fML_lct->Eval(2500) : fML_lct->Eval(pt);
		if( d2 < d2cut && m < mcutu && m > mcutl ){
			hb_lct->Fill(vvm, weight*lumi/40.);
		}
		if( l%(int(nEventB_lct/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEventB_lct, float(float(l)/nEventB_lct)*100) << std::endl;
	}
	// Signal
	std::cout << "LCT Signal filling ....." << std::endl;
	std::map<int, int> sigMap_lct;
	UInt_t sigCount_lct = 0;
	UInt_t sigArr_lct[10];
	SetBranch(trS_lct, ch);
	precut_lct = returnTCut(preSele);
        trS_lct->Draw(">>elist",precut_lct);
        TEventList *elistS_lct = (TEventList*)gROOT->FindObject("elist");
	long nEventS_lct = elistS_lct->GetN(); //long nEventS = trS->GetEntries();
	for(long l=0; l<nEventS_lct; l++){
	//for(long l=0; l<1000; l++){
        	trS_lct->GetEntry(elistS_lct->GetEntry(l)); //trS->GetEntry(l);
		float weightS = weight;
		if( sigWeight.find(dsid) != sigWeight.end() ){
			if( ch!="1Lep" ) weightS = weight * lumi * sigWeight[dsid];
			//if( ch=="1Lep" && sigType=="eff" ) weightS = weight / sigWeight[dsid];
			if( sigMap_lct.find(dsid) == sigMap_lct.end() ){
				sigMap_lct[dsid] = sigCount_lct; 
				sigArr_lct[sigCount_lct] = dsid;
				sigCount_lct++;
			}
		}
		// D2 &&  Mass window
		d2cut = pt > 2500 ? fD2_lct->Eval(2500) : fD2_lct->Eval(pt);
		mcutu = pt > 2500 ? fMU_lct->Eval(2500) : fMU_lct->Eval(pt);
		mcutl = pt > 2500 ? fML_lct->Eval(2500) : fML_lct->Eval(pt);
		if( d2 < d2cut && m < mcutu && m > mcutl ){
			hs_lct[sigMap_lct[dsid]]->Fill(vvm, weightS);
		}
		if( l%(int(nEventS_lct/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEventS_lct, float(float(l)/nEventS_lct)*100) << std::endl;
	}

	// Significance calculation
	std::cout << "Significance calculation....." << std::endl;
	float signi_lct[NmassPoints];
	float error_lct[NmassPoints];
	float ratio_lct[NmassPoints];
	float ratioE_lct[NmassPoints];
	std::vector<int> massList_lct;
	// D2 & Mass window
	for(UInt_t m=0; m<sigCount_lct; m++){
		signi_lct[m] = significance(hs_lct[m], hs_lct[m], hb_lct, "mass");
		error_lct[m] = significance(hs_lct[m], hs_lct[m], hb_lct, "mass", "error");
		hg_lct->SetBinContent(hg_lct->FindBin(stoi(sigName[sigArr_lct[m]])), signi_lct[m]);
		hg_lct->SetBinError(hg_lct->FindBin(stoi(sigName[sigArr_lct[m]])), 0.00001);
		hgr_lct->SetBinContent(hgr_lct->FindBin(stoi(sigName[sigArr_lct[m]])), 1.0);
	}

	std::cout << "LCTopo done ----------------" << std::endl;
	////////////////////////////////////////////////////
	
	
	TFile *fS = new TFile(filenameSig);
	TTree *trS = (TTree*)fS->Get("Nominal");
	TFile *fB= new TFile(filenameBkg);
	TTree *trB = (TTree*)fB->Get("Nominal");
	
	// Cut function
	std::vector<TString> cutWPs;
	std::string wpfile;

	TF1 *fD2[Ncombs];
	TF1 *fML[Ncombs];
	TF1 *fMU[Ncombs];
	TH1F *hb[Ncombs];
	TH1F *hs[Ncombs][NmassPoints];
	TH1F *hg[Ncombs];
	TH1F *hgr[Ncombs];
	int n = 0;
	for(auto i = combs.begin(); i != combs.end(); ++i){
		// D2 set
		wpfile = std::string("../output/WPs/SmoothedContained"+tagger+wps[i->second.at(0)]);
		cutWPs = retrieveWPs( wpfile );
		fD2[n] = new TF1(Form("D2_%i",n), cutWPs.at(2), 200, 2500);
		cutWPs.clear();
		// Upper M set
		wpfile = std::string("../output/WPs/SmoothedContained"+tagger+wps[i->second.at(1)]);
		cutWPs = retrieveWPs( wpfile );
		fMU[n] = new TF1(Form("MU_%i",n), cutWPs.at(1), 200, 2500);
		cutWPs.clear();
		// Lower M set
		wpfile = std::string("../output/WPs/SmoothedContained"+tagger+wps[i->second.at(2)]);
		cutWPs = retrieveWPs( wpfile );
		fML[n] = new TF1(Form("ML_%i",n), cutWPs.at(0), 200, 2500);
		cutWPs.clear();
		// Hist set
		hb[n]= new TH1F(Form("hb_%i",n), Form("hb_%i",n), nbins, rebin21);
		for(UInt_t l=0; l<NmassPoints; l++){
			hs[n][l]= new TH1F(Form("hs_%i_%i",n,l), Form("hs_%i_%i",n,l), nbins, rebin21);
		}
		hg[n]= new TH1F(Form("hg_%i",n), Form("hg_%i",n), 50,450,5450);//50
		hgr[n]= new TH1F(Form("hgr_%i",n), Form("hgr_%i",n), 50,450,5450);//50
		n++;
	}

	// Background
	std::cout << "Background filling ....." << std::endl;
	SetBranch(trB, ch);
	TCut precut = returnTCut(preSele);
        trB->Draw(">>elist",precut);
        TEventList *elistB = (TEventList*)gROOT->FindObject("elist");
	long nEventB = elistB->GetN(); //long nEventB = trB->GetEntries();
	for(long l=0; l<nEventB; l++){
        	trB->GetEntry(elistB->GetEntry(l)); //trB->GetEntry(l);
		for(int i=0; i<Ncombs; i++){
			// D2 &&  Mass window
			d2cut = pt > 2500 ? fD2[i]->Eval(2500) : fD2[i]->Eval(pt);
			mcutu = pt > 2500 ? fMU[i]->Eval(2500) : fMU[i]->Eval(pt);
			mcutl = pt > 2500 ? fML[i]->Eval(2500) : fML[i]->Eval(pt);
			if( d2 < d2cut && m < mcutu && m > mcutl ){
				hb[i]->Fill(vvm, weight*lumi/40.);
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
	precut = returnTCut(preSele);
        trS->Draw(">>elist",precut);
        TEventList *elistS = (TEventList*)gROOT->FindObject("elist");
	long nEventS = elistS->GetN(); //long nEventS = trS->GetEntries();
	for(long l=0; l<nEventS; l++){
	//for(long l=0; l<1000; l++){
        	trS->GetEntry(elistS->GetEntry(l)); //trS->GetEntry(l);
		float weightS = weight;
		if( sigWeight.find(dsid) != sigWeight.end() ){
			if( ch!="1Lep" ) weightS = weight * lumi * sigWeight[dsid];
			//if( ch=="1Lep" && sigType=="eff" ) weightS = weight / sigWeight[dsid];
			if( sigMap.find(dsid) == sigMap.end() ){
				sigMap[dsid] = sigCount; 
				sigArr[sigCount] = dsid;
				sigCount++;
			}
		}
		for(int i=0; i<Ncombs; i++){
			// D2 &&  Mass window
			d2cut = pt > 2500 ? fD2[i]->Eval(2500) : fD2[i]->Eval(pt);
			mcutu = pt > 2500 ? fMU[i]->Eval(2500) : fMU[i]->Eval(pt);
			mcutl = pt > 2500 ? fML[i]->Eval(2500) : fML[i]->Eval(pt);
			if( d2 < d2cut && m < mcutu && m > mcutl ){
				hs[i][sigMap[dsid]]->Fill(vvm, weightS);
			}
		}
		if( l%(int(nEventS/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEventS, float(float(l)/nEventS)*100) << std::endl;
	}

	// Significance calculation
	std::cout << "Significance calculation....." << std::endl;
	float signi[Ncombs][NmassPoints];
	float error[Ncombs][NmassPoints];
	float ratio[Ncombs][NmassPoints];
	float ratioE[Ncombs][NmassPoints];
	float maxSig=0.;
	float minSig=1E+10;
	TLine *ls[Ncombs][NmassPoints];
	std::vector<int> massList;
	// D2 & Mass window
	for(int i=0; i<Ncombs; i++){
		for(UInt_t m=0; m<sigCount; m++){
			signi[i][m] = significance(hs[i][m], hs[i][m], hb[i], "mass");
			error[i][m] = significance(hs[i][m], hs[i][m], hb[i], "mass", "error");
			hg[i]->SetBinContent(hg[i]->FindBin(stoi(sigName[sigArr[m]])), signi[i][m]);
			hg[i]->SetBinError(hg[i]->FindBin(stoi(sigName[sigArr[m]])), 0.00001);
			if( i==0 ) massList.push_back(stoi(sigName[sigArr[m]]));
			if( signi[i][m] > maxSig ) maxSig = signi[i][m];
			if( signi[i][m] < minSig ) minSig = signi[i][m];
		}
	}
	sort(massList.begin(), massList.end());


	// Plottintg .......................................
	//TCanvas *c1 = new TCanvas("c1","c1",600,600); //cv = new TCanvas("c1","Canvas",200,10,600,480);//600,480
	TCanvas *c1 = new TCanvas("c1","c1",200,10,560,480); //cv = new TCanvas("c1","Canvas",200,10,600,480);//600,480
	TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.18,0.95,0.97);//0.23
        TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.05,0.95,0.305);//0.05
        pad1->SetFillColor(kWhite);
        pad2->SetFillColor(kWhite);
        pad1->Draw();
        pad2->Draw();
        TLegend *leg= new TLegend(0.3,0.65,0.88,0.92);                                                    
        leg->SetBorderSize(0);

	gStyle->SetErrorX(0.0);
	gStyle->SetEndErrorSize(5);


	std::map<int, int> orderMap;
	std::map<int, int> op_orderMap;
	for(UInt_t m=0; m<sigCount; m++){
		for(UInt_t l=0; l<massList.size(); l++){
			if( stoi(sigName[sigArr[m]]) == massList.at(l) ){
				orderMap[m] = l;
				op_orderMap[l] = m;
			}
		}
	}

	std::map<int, int> orderMap_lct;
	std::map<int, int> op_orderMap_lct;
	for(UInt_t m=0; m<sigCount; m++){
		for(UInt_t l=0; l<massList.size(); l++){
			if( stoi(sigName[sigArr_lct[m]]) == massList.at(l) ){
				orderMap_lct[m] = l;
				op_orderMap_lct[l] = m;
			}
		}
	}

	// Ratio calculation
	int i=0;
	for(auto n = combs.begin(); n != combs.end(); ++n){
		for(UInt_t m=0; m<sigCount; m++){
			float v0  = hg[i]->GetBinContent(hg[i]->FindBin(stoi(sigName[sigArr[m]])));
			float er0 = hg[i]->GetBinError(hg[i]->FindBin(stoi(sigName[sigArr[m]])));
			//float v1  = hg[ratio_deno]->GetBinContent(hg[ratio_deno]->FindBin(stoi(sigName[sigArr[m]])));
			//float er1 = hg[ratio_deno]->GetBinError(hg[ratio_deno]->FindBin(stoi(sigName[sigArr[m]])));
			float v1  = hg_lct->GetBinContent(hg_lct->FindBin(stoi(sigName[sigArr[m]])));
			float er1 = hg_lct->GetBinError(hg_lct->FindBin(stoi(sigName[sigArr[m]])));

			float r = v0/v1;
			float rer = sqrt( pow(er0/v1,2) + pow(v0/v1/v1*er1,2) );
			if( isnan(r) || isinf(r) ) r = -100.;
			if( isnan(rer) || isinf(rer) ) rer = 0.;
			hgr[i]->SetBinContent(hgr[i]->FindBin(stoi(sigName[sigArr[m]])), r);
			hgr[i]->SetBinError(hgr[i]->FindBin(stoi(sigName[sigArr[m]])), rer);
			ratio[i][m] = r;
			error[i][m] = r;
		}
		i++;
	}

	pad1->cd();
	//pad1->SetGridy();
	pad1->SetBottomMargin(0.16);
	pad1->SetTopMargin(0.05);
	i=0;
	for(auto n = combs.begin(); n != combs.end(); ++n){
		hg[i]->SetStats(0);
		hg[i]->SetTitle("");
		hg[i]->GetXaxis()->SetLabelSize(0);
                hg[i]->GetYaxis()->SetTitleSize(0.14);
		hg[i]->GetXaxis()->SetTitle("Signal mass [GeV]");
		hg[i]->GetXaxis()->SetRangeUser(xMin,xMax);
		hg[i]->GetYaxis()->SetTitle("Significance");
		hg[i]->GetYaxis()->SetRangeUser(minSig*yScaleMin,maxSig*yScaleMax);
		hg[i]->SetLineColor(colorListSig.at(i));
		hg[i]->SetLineWidth(2);
		hg[i]->SetMarkerColor(colorListSig.at(i));
		hg[i]->SetMarkerStyle(8);
		hg[i]->SetMarkerSize(1);
		leg->AddEntry(hg[i], n->first, "lp");
		if( i!=0 ) hg[i]->Draw("E1 same");
		else hg[i]->Draw("E1");
		// For line connection
		for(UInt_t m=0; m<sigCount; m++){
			if( orderMap[m]+1 != int(sigCount) ){
				int next_m = op_orderMap[orderMap[m]+1];
				ls[i][m] = new TLine(stoi(sigName[sigArr[m]]), signi[i][m], stoi(sigName[sigArr[next_m]]), signi[i][next_m]);
				ls[i][m]->SetLineColor(colorListSig.at(i));
				ls[i][m]->SetLineStyle(1);
				ls[i][m]->SetLineWidth(2);
				ls[i][m]->Draw("same");
			}
		}
		i++;
	}
	// LCTopo plot
	hg_lct->SetStats(0);
	hg_lct->SetLineColor(kBlack);
	hg_lct->SetLineWidth(2);
	hg_lct->SetLineStyle(2);
	hg_lct->SetMarkerColor(kBlack);
	hg_lct->SetMarkerStyle(8);
	hg_lct->SetMarkerSize(0);
	leg->AddEntry(hg_lct, "LCTopo 50% fixed efficiency", "lp");
	hg_lct->Draw("E1 same");
	// For line connection
	TLine *ls_lct[NmassPoints];
	for(UInt_t m=0; m<sigCount; m++){
		if( orderMap_lct[m]+1 != int(sigCount) ){
			int next_m = op_orderMap_lct[orderMap_lct[m]+1];
			ls_lct[m] = new TLine(stoi(sigName[sigArr_lct[m]]), signi_lct[m], stoi(sigName[sigArr_lct[next_m]]), signi_lct[next_m]);
			ls_lct[m]->SetLineColor(kBlack);
			ls_lct[m]->SetLineStyle(2);
			ls_lct[m]->SetLineWidth(2);
			ls_lct[m]->Draw("same");
		}
	}
	leg->Draw();

	pad2->cd();
	pad2->SetGridy();
	pad2->SetBottomMargin(0.25);
	//pad2->SetTopMargin(0.1);
	i=0;
	TLine *lsr[Ncombs][NmassPoints];
	for(auto n = combs.begin(); n != combs.end(); ++n){
		hgr[i]->SetStats(0);
		hgr[i]->SetTitle("");
		hgr[i]->GetYaxis()->SetNdivisions(505);
		hgr[i]->GetXaxis()->SetTitle("Signal mass [GeV]");
		hgr[i]->GetXaxis()->SetRangeUser(xMin,xMax);
		hgr[i]->GetYaxis()->SetTitle("Ratio");
		hgr[i]->GetYaxis()->SetRangeUser(yRatioMin, yRatioMax);
		hgr[i]->SetLineColor(colorListSig.at(i));
		hgr[i]->SetLineWidth(2);
		hgr[i]->SetMarkerColor(colorListSig.at(i));
		hgr[i]->SetMarkerStyle(8);
		hgr[i]->SetMarkerSize(1);
		hgr[i]->GetXaxis()->SetTitleSize(0.14);
                hgr[i]->GetYaxis()->SetTitleSize(0.14);
                hgr[i]->GetXaxis()->SetLabelSize(0.14);
                hgr[i]->GetYaxis()->SetLabelSize(0.14);
                hgr[i]->GetXaxis()->SetTitleOffset(0.88);
                hgr[i]->GetYaxis()->SetTitleOffset(0.27);
		if( i!=0 ) hgr[i]->Draw("E1 same");
		else hgr[i]->Draw("E1");
		// For line connection
		for(UInt_t m=0; m<sigCount; m++){
			if( orderMap[m]+1 != int(sigCount) ){
				int next_m = op_orderMap[orderMap[m]+1];
				lsr[i][m] = new TLine(stoi(sigName[sigArr[m]]), ratio[i][m], stoi(sigName[sigArr[next_m]]), ratio[i][next_m]);
				lsr[i][m]->SetLineColor(colorListSig.at(i));
				lsr[i][m]->SetLineStyle(1);
				lsr[i][m]->SetLineWidth(2);
				lsr[i][m]->Draw("same");
			}
		}
		i++;
	}
	// LCTopo plot
	// For line connection
	TLine *lsr_lct = new TLine(xMin, 1.0, xMax, 1.0);
	lsr_lct->SetLineColor(kBlack);
	lsr_lct->SetLineStyle(2);
	lsr_lct->SetLineWidth(2);
	lsr_lct->Draw("same");
	
	// Significance summary
	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	std::cout << "Signal:	";	
	for(auto i = combs.begin(); i != combs.end(); ++i){
		std::cout << i->first << "	";
	}
	std::cout << std::endl;
	for(UInt_t m=0; m<sigCount; m++){
		std::cout << sigName[sigArr[m]] << "	";
		for(int i=0; i<Ncombs; i++){
			std::cout << Form("%0.3f", signi[i][m]) << "	";
		}
		std::cout << std::endl;
	}
	
	c1->Print("../output/sigificance."+ch+"."+tagger+"."+sigType+"."+nameOption+"."+Form("%d",Ncombs)+"."+outNameOpt+".png");
	

}


