#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TCut.h"
#include "TEventList.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "THStack.h"
#include "TF1.h"
#include "TLatex.h"

#include "/nfs_home/home1/myamatan/RootUtils/AtlasStyle.C"
#include "func.hh"
#include "sigList.hh"

//TString filenameSig = "../input/merged_truthStu/signal"+tagger+"/signal"+tagger+"."+ch+".root";
//TString filenameBkg = "../input/merged_truthStu/Zjets/Zjets."+ch+".root";
TString filenameSig = "../input/merged_new/signal"+tagger+"/signal"+tagger+"."+ch+".root";
TString filenameBkg = "../input/merged_new/Zjets/Zjets."+ch+".root";

// QG_ratio
float n_quark = 20007.3;
float n_gluon = 3062.39;
float qratio_max = 1.0;
float qratio_min = 0.0;
float qratio_d = 0.1;

TString value = "D2";
int nBin = 20;
float xMin = 0.0; 
float xMax = 5.0;
float yScaleMin = 0.8;
float yScaleMax = 2.0; // 1.5
float sigScale = 1.0;

void qgRatio(){


	SetAtlasStyle();

	if( value=="D2" ){
		nBin = 20;
		xMin = 0.0; 
		xMax = 5.0;
		sigScale = 1.0;
	}
	else if( value=="m" ){
		nBin = 12;
		xMin = 40.; //40
		xMax = 160.;//160
		sigScale = 0.4;
	}

	TFile *fS = new TFile(filenameSig);
	TTree *trS = (TTree*)fS->Get("Nominal");
	TFile *fB= new TFile(filenameBkg);
	TTree *trB = (TTree*)fB->Get("Nominal");

	// quarl ratio range
	qratio_max += qratio_d;
	UInt_t q_range = UInt_t( (qratio_max-qratio_min)/qratio_d );
	float q_scale = (n_gluon+n_quark) / n_quark;
	float g_scale = (n_gluon+n_quark) / n_gluon;

	// Fcuntion
	TF1 *fD2;
	TF1 *fML;
	TF1 *fMU;
	float d2cut = 0;
	float mcutu = 0;
	float mcutl = 0;
	std::string wpfile = std::string("../output/WPs/SmoothedContained"+tagger+wps["MaxVVJJ"]);
	std::vector<TString> cutWPs = retrieveWPs( wpfile );
	fD2 = new TF1(Form("D2"), cutWPs.at(2), 200, 2500);
	fMU = new TF1(Form("MU"), cutWPs.at(1), 200, 2500);
	fML = new TF1(Form("ML"), cutWPs.at(0), 200, 2500);

	// Histo initialization (x:VV mass, y:fatjet pt)
	TString nameID;
	TH1F *hbkg[q_range+1];
	TH1F *hbkgQ[q_range+1];
	TH1F *hbkgG[q_range+1];
	THStack *hbkgS[q_range+1];
	for(UInt_t i=0; i<q_range+1; i++){
		nameID = TString(std::to_string(i));
		hbkg[i] = new TH1F(nameID+"hbkg", nameID+"hbkg", nBin, xMin, xMax);
		hbkgQ[i] = new TH1F(nameID+"hbkgQ", nameID+"hbkgQ", nBin, xMin, xMax);
		hbkgG[i] = new TH1F(nameID+"hbkgG", nameID+"hbkgG", nBin, xMin, xMax);
		hbkgS[i] = new THStack();
	}
	TH1F *hsig = new TH1F("hsig", "hsig", nBin, xMin, xMax);


	float vtofill;
	/////////////////////////////
	// Background
	/////////////////////////////
	std::cout << "Background filling ....." << std::endl;
	SetBranch(trB, ch);
	float qsum  = 0.;
	float gsum  = 0.;
	// Filling
	TCut precut = returnTCut("beforeBosonTag");
        trB->Draw(">>elist",precut);
        TEventList *elistB = (TEventList*)gROOT->FindObject("elist");
	long nEventB = elistB->GetN();
	for(long l=0; l<nEventB; l++){
        	trB->GetEntry(elistB->GetEntry(l));
		// D2 &&  Mass window
		d2cut = pt > 2500 ? fD2->Eval(2500) : fD2->Eval(pt);
		mcutu = pt > 2500 ? fMU->Eval(2500) : fMU->Eval(pt);
		mcutl = pt > 2500 ? fML->Eval(2500) : fML->Eval(pt);
		vtofill = -99;
		if( value=="D2" ) vtofill = d2;
		else if( value=="m" ) vtofill = m;
		if(
		    //d2 < d2cut &&
		    m < mcutu &&
		    m > mcutl &&
		    //pt < 2500 && pt >= 650 && 
		    pt < 650 && pt >= 200 && 
		    1
		){
			// Raw filling
			if( pdgid==21 ) hbkgG[q_range]->Fill(vtofill, weight*lumi/40.);
			else hbkgQ[q_range]->Fill(vtofill, weight*lumi/40.);
			hbkg[q_range]->Fill(vtofill, weight*lumi/40.);
			// Ratio scan
			for(UInt_t i=0; i<q_range; i++){
				if( pdgid==21 ){
					hbkgG[i]->Fill(vtofill, weight*lumi/40.*g_scale*(1.-(qratio_min+qratio_d*i)));
					hbkg[i]->Fill(vtofill, weight*lumi/40.*g_scale*(1.-(qratio_min+qratio_d*i)));
				}else{
					hbkgQ[i]->Fill(vtofill, weight*lumi/40.*q_scale*(qratio_min+qratio_d*i));
					hbkg[i]->Fill(vtofill, weight*lumi/40.*q_scale*(qratio_min+qratio_d*i));
				}
			}
		}
		if( l%(int(nEventB/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEventB, float(float(l)/nEventB)*100) << std::endl;
	}

	std::cout << "Bkg q/g : " <<  hbkgQ[q_range]->Integral(0, nBin+1) << " / " << hbkgG[q_range]->Integral(0, nBin+1) << std::endl;
	std::cout << "q(g)_ratio	q	g	q+g" << std::endl;
	for(UInt_t i=0; i<q_range; i++){
		std::cout << qratio_min+qratio_d*i << "		" << hbkgQ[i]->Integral(0, nBin+1) << "	" << hbkgG[i]->Integral(0, nBin+1) << "	" << hbkgQ[i]->Integral(0, nBin+1)+hbkgG[i]->Integral(0, nBin+1)  << std::endl;
	}

	/////////////////////////////
	// Signal
	/////////////////////////////
	std::cout << "Signal filling ....." << std::endl;
	std::map<int, int> sigMap;
	UInt_t sigCount = 0;
	UInt_t sigArr[10];
	SetBranch(trS, ch);
	precut = returnTCut("beforeBosonTag");
        trS->Draw(">>elist",precut);
        TEventList *elistS = (TEventList*)gROOT->FindObject("elist");
	long nEventS = elistS->GetN();
	// Filling
	for(long l=0; l<nEventS; l++){
        	trS->GetEntry(elistS->GetEntry(l));
		// D2 &&  Mass window
		d2cut = pt > 2500 ? fD2->Eval(2500) : fD2->Eval(pt);
		mcutu = pt > 2500 ? fMU->Eval(2500) : fMU->Eval(pt);
		mcutl = pt > 2500 ? fML->Eval(2500) : fML->Eval(pt);
		vtofill = -99;
                if( value=="D2" ) vtofill = d2;
                else if( value=="m" ) vtofill = m;
		float weightS = weight;
		if( sigWeight.find(mc) != sigWeight.end() ){
			weightS = weight * lumi * sigWeight[mc];
			if( sigMap.find(mc) == sigMap.end() ){
				sigMap[mc] = sigCount; 
				sigArr[sigCount] = mc;
				sigCount++;
			}
		}
		if(
		    //d2 < d2cut &&
		    m < mcutu &&
		    m > mcutl &&
		    //pt < 2500 && pt >= 650 && 
		    pt < 650 && pt >= 200 && 
		    1
		){
			hsig->Fill(vtofill, weightS);
		}
		if( l%(int(nEventS/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEventS, float(float(l)/nEventS)*100) << std::endl;
	}

	// Plottintg .......................................
	TCanvas *c;
        TLegend *leg;   

/*
	// Ratio scan
	for(UInt_t i=0; i<q_range; i++){
		c = new TCanvas("cv","cv",600,600);
		c->cd();
		c->SetGridy();
        	leg = new TLegend(0.5,0.65,0.88,0.88);                    
        	leg->SetBorderSize(0);
		hsig->Scale(sigScale*hbkg[i]->Integral(0,nBin+1)/hsig->Integral(0,nBin+1));
		float k_ymax = hsig->GetBinContent( hsig->GetMaximumBin() );
		setHisto(hbkg[i], kBlack, kBlack, xMin, xMax, 0.0, k_ymax*yScaleMax, value, 1, false);
		setHisto(hbkgQ[i], kGreen-9, kBlack, xMin, xMax, 0.0, k_ymax*yScaleMax, value, 1, true);
		setHisto(hbkgG[i], kRed-9, kBlack, xMin, xMax, 0.0, k_ymax*yScaleMax, value, 1, true);
		hbkgS[i]->Add(hbkgQ[i]);
		hbkgS[i]->Add(hbkgG[i]);
		setHisto(hsig, kBlue-7, kBlue-7, xMin, xMax, 0.0, k_ymax*yScaleMax, value, 1, false);
		leg->AddEntry(hbkg[i], Form("Z+jets"), "l");
		leg->AddEntry(hbkgQ[i], Form("Quark jet (%0.1f)",(qratio_min+qratio_d*i)), "f");
		leg->AddEntry(hbkgG[i], Form("Gluon jet (%0.1f)",((1-(qratio_min+qratio_d*i)))), "f");
		leg->AddEntry(hsig, "Signal", "l");
		
		hsig->Draw("hist");
		//hbkgS[i]->Draw("hist same");
		hbkg[i]->Draw("hist same");
		hsig->Draw("hist same");
		c->RedrawAxis();
		leg->Draw();
		//c->Print("../output/qgRatio."+ch+"."+tagger+Form(".qRatio%d",i)+"."+value+".png");
		//c->Print("../output/qgRatio."+ch+"."+tagger+Form(".qRatio%d",i)+"."+value+".mcut.png");
		c->Print("../output/qgRatio."+ch+"."+tagger+Form(".qRatio%d",i)+"."+value+"200-650.mcut.png");
		//c->Print("../output/qgRatio."+ch+"."+tagger+Form(".qRatio%d",i)+"."+value+"650-2500.D2cut.png");
		delete c;
		delete leg;	hsig->Scale(hsig->Integral(0,nBin+1)/hbkg[i]->Integral(0,nBin+1)/sigScale);
		hsig->Scale(hsig->Integral(0,nBin+1)/hbkg[i]->Integral(0,nBin+1)/sigScale);
	}
*/

	// Normal Plot
	std::cout << "Plotting ......................." << std::endl;
	c = new TCanvas("cv","cv",600,500);
	TLatex *lat = new TLatex(0.25, 215E-3, "#it{#bf{ATLAS}} #bf{#bf{Simulation Preliminary}}");;
	TLatex *latsub = new TLatex(0.25, 200E-3, "#scale[0.8]{#bf{#bf{#sqrt{s}=13 TeV, after Jet mass window cut}}}");;
	//TLatex *latsub2 = new TLatex(0.25, 240E-3, "#scale[0.8]{#bf{#bf{200 < p_{T} < 650 Gev}}}");;
	TLatex *latsub2 = new TLatex(0.25, 185E-3, "#scale[0.8]{#bf{#bf{650 < p_{T} < 2500 Gev}}}");;
	c->cd();
	//c->SetGridy();
        leg = new TLegend(0.5,0.55,0.88,0.78);                    
        leg->SetBorderSize(0);
	//hsig->Scale(sigScale*hbkg[q_range]->Integral(0,nBin+1)/hsig->Integral(0,nBin+1));
	hsig->Scale(1./hsig->Integral(0,nBin+1));
	//hbkgQ[q_range]->Scale(1./hbkgQ[q_range]->Integral(0,nBin+1));
	//hbkgG[q_range]->Scale(1./hbkgG[q_range]->Integral(0,nBin+1));
	float k_ymax = hsig->GetBinContent( hsig->GetMaximumBin() );
	//float k_ymax = hsig->GetBinContent( hbkgQ[q_range]->GetMaximumBin() );
	//float q_ratio = hbkgQ[q_range]->Integral(0, nBin+1) / hbkg[q_range]->Integral(0, nBin+1);
	//float g_ratio = hbkgG[q_range]->Integral(0, nBin+1) / hbkg[q_range]->Integral(0, nBin+1);
	setHisto(hbkg[q_range], kBlack, kBlack, xMin, xMax, 0.0, k_ymax*yScaleMax, value, 1, false);
	setHisto(hbkgQ[q_range], kGreen-9, kGreen-3, xMin, xMax, 0.0, k_ymax*yScaleMax, value, 1, false);
	setHisto(hbkgG[q_range], kRed-9, kRed-3, xMin, xMax, 0.0, k_ymax*yScaleMax, value, 1, false);
	hbkgS[q_range]->Add(hbkgQ[q_range]);
	hbkgS[q_range]->Add(hbkgG[q_range]);
	setHisto(hsig, kBlue-7, kBlue-7, xMin, xMax, 0.0, k_ymax*yScaleMax, value, 1, false);
	//leg->AddEntry(hbkg[q_range], Form("Z+jets"), "l");
	//leg->AddEntry(hbkgQ[q_range], Form("Quark jet (%0.2f)",q_ratio), "f");
	//leg->AddEntry(hbkgG[q_range], Form("Gluon jet (%0.2f)",g_ratio), "f");
	leg->AddEntry(hbkgQ[q_range], Form("Quark-induced jet"), "l");
	leg->AddEntry(hbkgG[q_range], Form("Gluon-induced jet"), "l");
	//leg->AddEntry(hsig, "Signal", "l");
	
	hsig->Draw("hist");
	//hbkgS[q_range]->Draw("hist same");
	//hbkg[q_range]->Draw("hist same");
	//hbkgG[q_range]->Draw("hist");
	//hbkgQ[q_range]->Draw("hist same");
	//hsig->Draw("hist same");

	c->RedrawAxis();
	//leg->Draw();
	//lat->Draw();
	//latsub->Draw();
	//latsub2->Draw();

	//c->Print("../output/qgRatio."+ch+"."+tagger+Form(".qRatioRaw")+"."+value+".png");
	//c->Print("../output/qgRatio."+ch+"."+tagger+Form(".qRatioRaw")+"."+value+"mcut.png");
	c->Print("../output/qgRatio."+ch+"."+tagger+Form(".qRatioRaw")+"."+value+"200-650.mcut.pdf");
	//c->Print("../output/qgRatio."+ch+"."+tagger+Form(".qRatioRaw")+"."+value+"650-2500.mcut.pdf");


}
