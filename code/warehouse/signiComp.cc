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

TString filenameSig = "../input/merged_tcc/signal"+tagger+"/signal"+tagger+"."+ch+".root";
TString filenameBkg = "../input/merged_tcc/allMC/allMC."+ch+".root";

TString outNameOpt = "comb";
float xMin = 250; // 450
float xMax = 5450;// 5450
float yScaleMin = 0.7;   // 0.7
float yScaleMax = 1.6;   // 1.6
float yRatioMin = 0.7;   // 0.7
float yRatioMax = 1.3;  // 1.3 

int ratio_deno = 1;

void signiComp(){
	
	TFile *fS = new TFile(filenameSig);
	TTree *trS = (TTree*)fS->Get("Nominal");
	TFile *fB= new TFile(filenameBkg);
	TTree *trB = (TTree*)fB->Get("Nominal");
	
	// Cut function
	const int Nwps = wps.size();
	const int Ncombs = combs.size();
	float d2cut = 0;
	float mcutu = 0;
	float mcutl = 0;
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
			if( 
			    //d2 < d2cut && 
			    d2 > d2cut && 
 			    m < mcutu  &&
			    m > mcutl  &&
			    1
			){
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
			if( 
			    //d2 < d2cut && 
			    d2 > d2cut && 
 			    m < mcutu  &&
			    m > mcutl  &&
			    1
			){
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
			//hg[i]->SetBinError(hg[i]->FindBin(stoi(sigName[sigArr[m]])), error[i][m]);
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

	// Ratio calculation
	int i=0;
	for(auto n = combs.begin(); n != combs.end(); ++n){
		for(UInt_t m=0; m<sigCount; m++){
			float v0  = hg[i]->GetBinContent(hg[i]->FindBin(stoi(sigName[sigArr[m]])));
			float er0 = hg[i]->GetBinError(hg[i]->FindBin(stoi(sigName[sigArr[m]])));
			float v1  = hg[ratio_deno]->GetBinContent(hg[ratio_deno]->FindBin(stoi(sigName[sigArr[m]])));
			float er1 = hg[ratio_deno]->GetBinError(hg[ratio_deno]->FindBin(stoi(sigName[sigArr[m]])));

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


