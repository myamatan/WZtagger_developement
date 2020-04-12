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
#include "TLatex.h"

#include "/nfs_home/home1/myamatan/RootUtils/AtlasStyle.C"
#include "sigList.hh"
#include "func.hh"

TString filenameSig = "../input/merged_tcc_new/signal"+tagger+"/signal"+tagger+"."+ch+".root";
TString filenameBkg = "../input/merged_tcc_new/allMC/allMC."+ch+".root";
TString filenameSig_lct = "../input/merged_lct_new/signal"+tagger+"/signal"+tagger+"."+ch+".root";
TString filenameBkg_lct = "../input/merged_lct_new/allMC/allMC."+ch+".root";

TString outNameOpt = "comb";
float xMin = 250; // 450
float xMax = 5450;// 5450
float yScaleMin = 0.7;   // 0.7
float yScaleMax = 2.3;   // 1.6
float yRatioMin = 0.7;   // 0.7
float yRatioMax = 1.5;  // 1.3 

int ratio_deno = 1;

const int Nwps = wps.size();
const int Ncombs = combs.size();
float d2cut = 0;
float mcutu = 0;
float mcutl = 0;
float d2cut_hp = 0;
float mcutu_hp = 0;
float mcutl_hp = 0;
float d2cut_lp = 0;
float mcutu_lp = 0;
float mcutl_lp = 0;

// compared to LCTopo with HP + LP
void signiComp_LCT_HPLP(){

	SetAtlasStyle();

	////////////////////////////////////////////////////
	// LCTopo
	////////////////////////////////////////////////////
	TFile *fS_lct  = new TFile(filenameSig_lct);
	TTree *trS_lct = (TTree*)fS_lct->Get("Nominal");
	TFile *fB_lct  = new TFile(filenameBkg_lct);
	TTree *trB_lct = (TTree*)fB_lct->Get("Nominal");

	// HP cut
	std::vector<TString> cutWPs_lct = retrieveWPs( std::string("../output/WPs/SmoothedContained"+tagger+wps["Fix50VVJJ_LCT"]) );
	TF1 *fD2_lct_hp;
	TF1 *fML_lct_hp;
	TF1 *fMU_lct_hp;
	fD2_lct_hp = new TF1(Form("D2_LCT_HP"), cutWPs_lct.at(2), 200, 2500);
	fMU_lct_hp = new TF1(Form("MU_LCT_HP"), cutWPs_lct.at(1), 200, 2500);
	fML_lct_hp = new TF1(Form("ML_LCT_HP"), cutWPs_lct.at(0), 200, 2500);
	cutWPs_lct.clear();
	// LP cut
	cutWPs_lct = retrieveWPs( std::string("../output/WPs/SmoothedContained"+tagger+wps["Fix80VVJJ_LCT"]) );
	TF1 *fD2_lct_lp;
	TF1 *fML_lct_lp;
	TF1 *fMU_lct_lp;
	fD2_lct_lp = new TF1(Form("D2_LCT_LP"), cutWPs_lct.at(2), 200, 2500);
	fMU_lct_lp = new TF1(Form("MU_LCT_LP"), cutWPs_lct.at(1), 200, 2500);
	fML_lct_lp = new TF1(Form("ML_LCT_LP"), cutWPs_lct.at(0), 200, 2500);

	// HP hist
	TH1F *hb_lct_hp  = new TH1F("hb_lct_hp","hb_lct_hp", nbins, rebin21);
	TH1F *hb_lct_hp_af  = new TH1F("hb_lct_hp_af","hb_lct_hp_af", nbins, rebin21); // af -> After fit
	TH1F *hs_lct_hp[NmassPoints];
	for(UInt_t l=0; l<NmassPoints; l++){
		hs_lct_hp[l]= new TH1F(Form("hs_lct_hp_%i",l), Form("hs_lct_hp_%i",l), nbins, rebin21);
	}
	// LP hist
	TH1F *hb_lct_lp  = new TH1F("hb_lct_lp","hb_lct_lp", nbins, rebin21);
	TH1F *hb_lct_lp_af  = new TH1F("hb_lct_lp_af","hb_lct_lp_af", nbins, rebin21);
	TH1F *hs_lct_lp[NmassPoints];
	for(UInt_t l=0; l<NmassPoints; l++){
		hs_lct_lp[l]= new TH1F(Form("hs_lct_lp_%i",l), Form("hs_lct_lp_%i",l), nbins, rebin21);
	}

	TH1F *hg_lct  = new TH1F("hg_lct","hg_lct", 50,450,5450);
	TH1F *hgr_lct = new TH1F("hgr_lct","hgr_lct", 50,450,5450);

	// Background
	std::cout << "LCT Background filling ....." << std::endl;
	SetBranch(trB_lct, ch);
	TCut precut_lct = returnTCut(preSele);
        trB_lct->Draw(">>elist",precut_lct);
        TEventList *elistB_lct = (TEventList*)gROOT->FindObject("elist");
	long nEventB_lct = elistB_lct->GetN(); //long nEventB = trB->GetEntries();
	for(long l=0; l<nEventB_lct; l++){
        	trB_lct->GetEntry(elistB_lct->GetEntry(l)); //trB->GetEntry(l);
		// D2 &&  Mass window (HP)
		d2cut_hp = pt > 2500 ? fD2_lct_hp->Eval(2500) : fD2_lct_hp->Eval(pt);
		mcutu_hp = pt > 2500 ? fMU_lct_hp->Eval(2500) : fMU_lct_hp->Eval(pt);
		mcutl_hp = pt > 2500 ? fML_lct_hp->Eval(2500) : fML_lct_hp->Eval(pt);
		// D2 &&  Mass window (LP)
		d2cut_lp = pt > 2500 ? fD2_lct_lp->Eval(2500) : fD2_lct_lp->Eval(pt);
		mcutu_lp = pt > 2500 ? fMU_lct_lp->Eval(2500) : fMU_lct_lp->Eval(pt);
		mcutl_lp = pt > 2500 ? fML_lct_lp->Eval(2500) : fML_lct_lp->Eval(pt);
		if( d2 < d2cut_hp && m < mcutu_hp && m > mcutl_hp ){
			hb_lct_hp->Fill(vvm, weight*lumi/40.);
		}else if( d2 < d2cut_lp && m < mcutu_lp && m > mcutl_lp ){
			hb_lct_lp->Fill(vvm, weight*lumi/40.);
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
		// D2 &&  Mass window (HP)
		d2cut_hp = pt > 2500 ? fD2_lct_hp->Eval(2500) : fD2_lct_hp->Eval(pt);
		mcutu_hp = pt > 2500 ? fMU_lct_hp->Eval(2500) : fMU_lct_hp->Eval(pt);
		mcutl_hp = pt > 2500 ? fML_lct_hp->Eval(2500) : fML_lct_hp->Eval(pt);
		// D2 &&  Mass window (LP)
		d2cut_lp = pt > 2500 ? fD2_lct_lp->Eval(2500) : fD2_lct_lp->Eval(pt);
		mcutu_lp = pt > 2500 ? fMU_lct_lp->Eval(2500) : fMU_lct_lp->Eval(pt);
		mcutl_lp = pt > 2500 ? fML_lct_lp->Eval(2500) : fML_lct_lp->Eval(pt);
		if( d2 < d2cut_hp && m < mcutu_hp && m > mcutl_hp ){
			hs_lct_hp[sigMap_lct[dsid]]->Fill(vvm, weightS);
		}else if( d2 < d2cut_lp && m < mcutu_lp && m > mcutl_lp ){
			hs_lct_lp[sigMap_lct[dsid]]->Fill(vvm, weightS);
		}
		if( l%(int(nEventS_lct/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEventS_lct, float(float(l)/nEventS_lct)*100) << std::endl;
	}

	// Convert to "per GeV"
	for(int i=1; i<nbins+1; i++){
		float vwidth = rebin21[i] - rebin21[i-1];
		hb_lct_hp->SetBinContent(i, hb_lct_hp->GetBinContent(i)/vwidth);
		hb_lct_lp->SetBinContent(i, hb_lct_lp->GetBinContent(i)/vwidth);
		hb_lct_hp->SetBinError(i, hb_lct_hp->GetBinError(i)/vwidth);
		hb_lct_lp->SetBinError(i, hb_lct_lp->GetBinError(i)/vwidth);
		for(UInt_t m=0; m<sigCount_lct; m++){
			hs_lct_hp[m]->SetBinContent(i, hs_lct_hp[m]->GetBinContent(i)/vwidth);
			hs_lct_lp[m]->SetBinContent(i, hs_lct_hp[m]->GetBinContent(i)/vwidth);
			hs_lct_hp[m]->SetBinError(i, hs_lct_hp[m]->GetBinError(i)/vwidth);
			hs_lct_lp[m]->SetBinError(i, hs_lct_hp[m]->GetBinError(i)/vwidth);
		}
	}

	///////////////////////
	// Bkg m(llJ) smoothing
	///////////////////////
	// Fit
	TF1 *f_lct_hp = new TF1("f_lct_hp","[0]*(1-x/sqrt(13e+3))/( (x/sqrt(13e+3))^([1]+[2]*log(x/sqrt(13e+3))) )", fitRangeMin, 6000);
	TF1 *f_lct_lp = new TF1("f_lct_lp","[0]*(1-x/sqrt(13e+3))/( (x/sqrt(13e+3))^([1]+[2]*log(x/sqrt(13e+3))) )", fitRangeMin, 6000);

	f_lct_hp->SetParameter(0, -7e+0);
        f_lct_hp->SetParameter(1, 8.1);
        f_lct_hp->SetParameter(2, 0.14);
        f_lct_hp->SetParLimits(0, -1E+8, 0);
        f_lct_hp->SetParLimits(1, 0, 10);
        f_lct_hp->SetParLimits(2, 0, 10);

	f_lct_lp->SetParameter(0, -7e+0);
        f_lct_lp->SetParameter(1, 8.1);
        f_lct_lp->SetParameter(2, 0.14);
        f_lct_lp->SetParLimits(0, -1E+8, 0);
        f_lct_lp->SetParLimits(1, 0, 10);
        f_lct_lp->SetParLimits(2, 0, 10);


	// Fit
	hb_lct_hp->Fit("f_lct_hp","N","",fitRangeMin,6000);	
	hb_lct_lp->Fit("f_lct_lp","N","",fitRangeMin,6000);	

	float chi2_lct_hp = f_lct_hp->GetChisquare() / f_lct_hp->GetNDF();
	float chi2_lct_lp = f_lct_lp->GetChisquare() / f_lct_lp->GetNDF();

	// Smoothing
	for(int i=1; i<nbins+1; i++){
		if( rebin21[i] > smoothXmin ){
			float vwidth = rebin21[i] - rebin21[i-1];
			hb_lct_hp_af->SetBinContent(i, f_lct_hp->Integral(rebin21[i-1],rebin21[i])/vwidth);
			hb_lct_lp_af->SetBinContent(i, f_lct_lp->Integral(rebin21[i-1],rebin21[i])/vwidth);
			hb_lct_hp_af->SetBinError(i, f_lct_hp->IntegralError(rebin21[i-1],rebin21[i])/sqrt(vwidth));
			hb_lct_lp_af->SetBinError(i, f_lct_lp->IntegralError(rebin21[i-1],rebin21[i])/sqrt(vwidth));
		}else{
			hb_lct_hp_af->SetBinContent(i, hb_lct_hp->GetBinContent(i));
			hb_lct_lp_af->SetBinContent(i, hb_lct_lp->GetBinContent(i));
			hb_lct_hp_af->SetBinError(i, hb_lct_hp->GetBinError(i));
			hb_lct_lp_af->SetBinError(i, hb_lct_lp->GetBinError(i));
		}
	}
	

	// Plot
	TCanvas *cv = new TCanvas("cv","cv",200,10,560,480);
        TLegend *leg2;                                                   
	TLatex *lat2 = new TLatex(500, 1E+5*0.1, "#it{#bf{ATLAS}} #bf{#bf{Simulation Work in progress}}");
	TLatex *lat_lct_hp = new TLatex(2500, 1E-1, Form("#scale[0.8]{#bf{#bf{LCTopo,HP (chi2/NDF=%3.2f)}}}",chi2_lct_hp));
	TLatex *lat_lct_lp = new TLatex(2500, 1E-1, Form("#scale[0.8]{#bf{#bf{LCTopo,LP (chi2/NDF=%3.2f)}}}",chi2_lct_lp));
	cv->SetLogy();
	cv->cd();
        leg2= new TLegend(0.45,0.61,0.90,0.85);                                                    
        leg2->SetBorderSize(0);
	leg2->SetFillStyle(0);
	hb_lct_hp_af->GetYaxis()->SetRangeUser(1E-7,1E+5);
	hb_lct_hp_af->GetXaxis()->SetTitle("m_{T}(vvJ)[GeV]");
	hb_lct_hp_af->SetLineColor(kBlue);
	hb_lct_hp_af->Draw("hist");
	f_lct_hp->SetLineColor(kRed);
	f_lct_hp->DrawCopy("same");
	f_lct_hp->SetFillColor(kRed);
	f_lct_hp->SetFillStyle(3018);
	f_lct_hp->SetLineWidth(2);
	f_lct_hp->SetMarkerSize(0);
	f_lct_hp->Draw("E3 same");
	hb_lct_hp->Draw("E0 same");
	leg2->AddEntry(hb_lct_hp,   "Pre-fit  Bkg","lep");
	leg2->AddEntry(f_lct_hp,    "Post-fit Bkg function","l");
	leg2->AddEntry(hb_lct_hp_af,"Post-fit Bkg (Binned)","l");
	if( tagger=="W" ){
		hs_lct_hp[sigMap_lct[302265]]->SetLineColor(kGreen);
		hs_lct_hp[sigMap_lct[302265]]->Draw("hist same");
		leg2->AddEntry(hs_lct_hp[sigMap_lct[302265]],"HVTWZ 5TeV","l");
	}else{
		hs_lct_hp[sigMap_lct[303326]]->SetLineColor(kGreen);
		hs_lct_hp[sigMap_lct[303326]]->Draw("hist same");
		leg2->AddEntry(hs_lct_hp[sigMap_lct[303326]],"RSG_ZZ 5TeV","l");
	}
	leg2->Draw();
	lat2->Draw();
	lat_lct_hp->Draw();
	cv->Print("../output/smoothing.lct.hp."+ch+"."+tagger+"."+sigType+"."+nameOption+"."+Form("%d",Ncombs)+"."+outNameOpt+".png");
	delete leg2;

        leg2= new TLegend(0.45,0.61,0.90,0.85);                       
        leg2->SetBorderSize(0);
	leg2->SetFillStyle(0);
	hb_lct_lp_af->GetYaxis()->SetRangeUser(1E-7,1E+5);
	hb_lct_lp_af->GetXaxis()->SetTitle("m_{T}(vvJ)[GeV]");
	hb_lct_lp_af->SetLineColor(kBlue);
	hb_lct_lp_af->Draw("hist");
	f_lct_lp->SetLineColor(kRed);
	f_lct_lp->DrawCopy("same");
	f_lct_lp->SetFillColor(kRed);
	f_lct_lp->SetFillStyle(3018);
	f_lct_lp->SetLineWidth(2);
	f_lct_lp->SetMarkerSize(0);
	f_lct_lp->Draw("E3 same");
	hb_lct_lp->Draw("E0 same");
	leg2->AddEntry(hb_lct_lp,   "Pre-fit  Bkg","lep");
	leg2->AddEntry(f_lct_lp,    "Post-fit Bkg function","l");
	leg2->AddEntry(hb_lct_lp_af,"Post-fit Bkg (Binned)","l");
	if( tagger=="W" ){
		hs_lct_lp[sigMap_lct[302265]]->SetLineColor(kGreen);
		hs_lct_lp[sigMap_lct[302265]]->Draw("hist same");
		leg2->AddEntry(hs_lct_lp[sigMap_lct[302265]],"HVTWZ 5TeV","l");
	}else{
		hs_lct_lp[sigMap_lct[303326]]->SetLineColor(kGreen);
		hs_lct_lp[sigMap_lct[303326]]->Draw("hist same");
		leg2->AddEntry(hs_lct_lp[sigMap_lct[303326]],"RSG_ZZ 5TeV","l");
	}
	leg2->Draw();
	lat2->Draw();
	lat_lct_lp->Draw();
	cv->Print("../output/smoothing.lct.lp."+ch+"."+tagger+"."+sigType+"."+nameOption+"."+Form("%d",Ncombs)+"."+outNameOpt+".png");
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        delete leg2;                                                   
	delete lat2;
	delete lat_lct_hp;
	delete lat_lct_lp;


	// Significance calculation
	std::cout << "Significance calculation....." << std::endl;

	float signi_lct_hp[NmassPoints];
	float error_lct_hp[NmassPoints];
	float ratio_lct_hp[NmassPoints];
	float ratioE_lct_hp[NmassPoints];

	float signi_lct_lp[NmassPoints];
	float error_lct_lp[NmassPoints];
	float ratio_lct_lp[NmassPoints];
	float ratioE_lct_lp[NmassPoints];

	float signi_lct[NmassPoints];
	float error_lct[NmassPoints];
	std::vector<int> massList_lct;
	// D2 & Mass window
	for(UInt_t m=0; m<sigCount_lct; m++){
		if( isMCSmooth ){
			signi_lct_hp[m] = significance(hs_lct_hp[m], hs_lct_hp[m], hb_lct_hp_af, "mass");
			error_lct_hp[m] = significance(hs_lct_hp[m], hs_lct_hp[m], hb_lct_hp_af, "mass", "error");
			signi_lct_lp[m] = significance(hs_lct_lp[m], hs_lct_lp[m], hb_lct_lp_af, "mass");
			error_lct_lp[m] = significance(hs_lct_lp[m], hs_lct_lp[m], hb_lct_lp_af, "mass", "error");
		}else{
			signi_lct_hp[m] = significance(hs_lct_hp[m], hs_lct_hp[m], hb_lct_hp, "mass");
			error_lct_hp[m] = significance(hs_lct_hp[m], hs_lct_hp[m], hb_lct_hp, "mass", "error");
			signi_lct_lp[m] = significance(hs_lct_lp[m], hs_lct_lp[m], hb_lct_lp, "mass");
			error_lct_lp[m] = significance(hs_lct_lp[m], hs_lct_lp[m], hb_lct_lp, "mass", "error");
		}

		signi_lct[m] = sqrt( pow(signi_lct_hp[m],2) + pow(signi_lct_lp[m],2) );
		//signi_lct[m] = sqrt( pow(signi_lct_hp[m],2) );

		error_lct[m] = sqrt( pow(error_lct_hp[m],2) + pow(error_lct_lp[m],2) );
		hg_lct->SetBinContent(hg_lct->FindBin(stoi(sigName[sigArr_lct[m]])), signi_lct[m]);
		hg_lct->SetBinError(hg_lct->FindBin(stoi(sigName[sigArr_lct[m]])), 0.00001);
		hgr_lct->SetBinContent(hgr_lct->FindBin(stoi(sigName[sigArr_lct[m]])), 1.0);
	}

	std::cout << "LCTopo done ----------------" << std::endl;


	////////////////////////////////////////////////////
	// TCC
	////////////////////////////////////////////////////
	
	TFile *fS = new TFile(filenameSig);
	TTree *trS = (TTree*)fS->Get("Nominal");
	TFile *fB= new TFile(filenameBkg);
	TTree *trB = (TTree*)fB->Get("Nominal");
	
	// Cut function
	std::vector<TString> cutWPs;
	std::string wpfile;

	TH1F *hg[Ncombs];
	TH1F *hgr[Ncombs];

	// HP cut function & histogram
	TF1 *fD2_hp[Ncombs];
	TF1 *fML_hp[Ncombs];
	TF1 *fMU_hp[Ncombs];
	TH1F *hb_hp[Ncombs];
	TH1F *hb_hp_af[Ncombs];
	TH1F *hs_hp[Ncombs][NmassPoints];
	int n = 0;
	for(auto i = combs.begin(); i != combs.end(); ++i){
		// D2 set
		wpfile = std::string("../output/WPs/SmoothedContained"+tagger+wps[i->second.at(0)]);
		cutWPs = retrieveWPs( wpfile );
		fD2_hp[n] = new TF1(Form("D2_hp_%i",n), cutWPs.at(2), 200, 2500);
		cutWPs.clear();
		// Upper M set
		wpfile = std::string("../output/WPs/SmoothedContained"+tagger+wps[i->second.at(1)]);
		cutWPs = retrieveWPs( wpfile );
		fMU_hp[n] = new TF1(Form("MU_hp_%i",n), cutWPs.at(1), 200, 2500);
		cutWPs.clear();
		// Lower M set
		wpfile = std::string("../output/WPs/SmoothedContained"+tagger+wps[i->second.at(2)]);
		cutWPs = retrieveWPs( wpfile );
		fML_hp[n] = new TF1(Form("ML_hp_%i",n), cutWPs.at(0), 200, 2500);
		cutWPs.clear();
		// Hist set
		hb_hp[n]= new TH1F(Form("hb_hp_%i",n), Form("hb_hp_%i",n), nbins, rebin21);
		hb_hp_af[n]= new TH1F(Form("hb_hp_af_%i",n), Form("hb_hp_af_%i",n), nbins, rebin21);
		for(UInt_t l=0; l<NmassPoints; l++){
			hs_hp[n][l]= new TH1F(Form("hs_hp_%i_%i",n,l), Form("hs_hp_%i_%i",n,l), nbins, rebin21);
		}
		hg[n]= new TH1F(Form("hg_%i",n), Form("hg_%i",n), 50,450,5450);//50
		hgr[n]= new TH1F(Form("hgr_%i",n), Form("hgr_%i",n), 50,450,5450);//50
		n++;
	}
	// LP cut function & histogram
	TF1 *fD2_lp[Ncombs];
	TF1 *fML_lp[Ncombs];
	TF1 *fMU_lp[Ncombs];
	TH1F *hb_lp[Ncombs];
	TH1F *hb_lp_af[Ncombs];
	TH1F *hs_lp[Ncombs][NmassPoints];
	n = 0;
	for(auto i = combs_LP.begin(); i != combs_LP.end(); ++i){
		// D2 set
		wpfile = std::string("../output/WPs/SmoothedContained"+tagger+wps[i->second.at(0)]);
		cutWPs = retrieveWPs( wpfile );
		fD2_lp[n] = new TF1(Form("D2_lp_%i",n), cutWPs.at(2), 200, 2500);
		cutWPs.clear();
		// Upper M set
		wpfile = std::string("../output/WPs/SmoothedContained"+tagger+wps[i->second.at(1)]);
		cutWPs = retrieveWPs( wpfile );
		fMU_lp[n] = new TF1(Form("MU_lp_%i",n), cutWPs.at(1), 200, 2500);
		cutWPs.clear();
		// Lower M set
		wpfile = std::string("../output/WPs/SmoothedContained"+tagger+wps[i->second.at(2)]);
		cutWPs = retrieveWPs( wpfile );
		fML_lp[n] = new TF1(Form("ML_lp_%i",n), cutWPs.at(0), 200, 2500);
		cutWPs.clear();
		// Hist set
		hb_lp[n]= new TH1F(Form("hb_lp_%i",n), Form("hb_lp_%i",n), nbins, rebin21);
		hb_lp_af[n]= new TH1F(Form("hb_lp_af_%i",n), Form("hb_lp_af_%i",n), nbins, rebin21);
		for(UInt_t l=0; l<NmassPoints; l++){
			hs_lp[n][l]= new TH1F(Form("hs_lp_%i_%i",n,l), Form("hs_lp_%i_%i",n,l), nbins, rebin21);
		}
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
			// D2 &&  Mass window (HP)
			d2cut_hp = fD2_hp[i]->Eval(pt);
			mcutu_hp = fMU_hp[i]->Eval(pt);
			mcutl_hp = fML_hp[i]->Eval(pt);
			// D2 &&  Mass window (LP)
			d2cut_lp = fD2_lp[i]->Eval(pt);
			mcutu_lp = fMU_lp[i]->Eval(pt);
			mcutl_lp = fML_lp[i]->Eval(pt);
			if(
				d2 < d2cut_hp &&
				m < mcutu_hp  &&
				m > mcutl_hp
			 ){
				hb_hp[i]->Fill(vvm, weight*lumi/40.);
			}else if(
				//d2 < d2cut_lp &&
				d2 >= d2cut_hp &&
				m < mcutu_hp  &&
				m > mcutl_hp
			){
				hb_lp[i]->Fill(vvm, weight*lumi/40.);
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
			// D2 &&  Mass window (HP)
			d2cut_hp = fD2_hp[i]->Eval(pt);
			mcutu_hp = fMU_hp[i]->Eval(pt);
			mcutl_hp = fML_hp[i]->Eval(pt);
			// D2 &&  Mass window (LP)
			d2cut_lp = fD2_lp[i]->Eval(pt);
			mcutu_lp = fMU_lp[i]->Eval(pt);
			mcutl_lp = fML_lp[i]->Eval(pt);
			if(
				d2 < d2cut_hp &&
				m < mcutu_hp  &&
				m > mcutl_hp
			){
				hs_hp[i][sigMap[dsid]]->Fill(vvm, weightS);
			}else if(	
				//d2 < d2cut_lp &&
				d2 >= d2cut_hp &&
				m < mcutu_hp  &&
				m > mcutl_hp
			){
				hs_lp[i][sigMap[dsid]]->Fill(vvm, weightS);
			}
		}
		if( l%(int(nEventS/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEventS, float(float(l)/nEventS)*100) << std::endl;
	}


	// Convert to "per GeV"
	for(int s=0; s<Ncombs; s++){
		for(int i=1; i<nbins+1; i++){
			float vwidth = rebin21[i] - rebin21[i-1];
			hb_hp[s]->SetBinContent(i, hb_hp[s]->GetBinContent(i)/vwidth);
			hb_lp[s]->SetBinContent(i, hb_lp[s]->GetBinContent(i)/vwidth);
			hb_hp[s]->SetBinError(i, hb_hp[s]->GetBinError(i)/vwidth);
			hb_lp[s]->SetBinError(i, hb_lp[s]->GetBinError(i)/vwidth);
			for(UInt_t m=0; m<sigCount; m++){
				hs_hp[s][m]->SetBinContent(i, hs_hp[s][m]->GetBinContent(i)/vwidth);
				hs_lp[s][m]->SetBinContent(i, hs_hp[s][m]->GetBinContent(i)/vwidth);
				hs_hp[s][m]->SetBinError(i, hs_hp[s][m]->GetBinError(i)/vwidth);
				hs_lp[s][m]->SetBinError(i, hs_hp[s][m]->GetBinError(i)/vwidth);
			}
		}
	}


	///////////////////////
	// Bkg m(llJ) smoothing
	///////////////////////
	// Fit
	TF1 *f_tcc_hp[Ncombs];
	TF1 *f_tcc_lp[Ncombs];
		
	TLatex *lat_tcc_hp;
	TLatex *lat_tcc_lp;

	for(int i=0; i<Ncombs; i++){

		f_tcc_hp[i] = new TF1(Form("f_tcc_hp_%d",i),"[0]*(1-x/sqrt(13e+3))/( (x/sqrt(13e+3))^([1]+[2]*log(x/sqrt(13e+3))) )", fitRangeMin, 6000);
		f_tcc_lp[i] = new TF1(Form("f_tcc_lp_%d",i),"[0]*(1-x/sqrt(13e+3))/( (x/sqrt(13e+3))^([1]+[2]*log(x/sqrt(13e+3))) )", fitRangeMin, 6000);
	
		f_tcc_hp[i]->SetParameter(0, -7e+0);
	        f_tcc_hp[i]->SetParameter(1, 8.1);
	        f_tcc_hp[i]->SetParameter(2, 0.14);
	        f_tcc_hp[i]->SetParLimits(0, -1E+8, 0);
	        f_tcc_hp[i]->SetParLimits(1, 0, 10);
	        f_tcc_hp[i]->SetParLimits(2, 0, 10);
	
		f_tcc_lp[i]->SetParameter(0, -7e+0);
	        f_tcc_lp[i]->SetParameter(1, 8.1);
	        f_tcc_lp[i]->SetParameter(2, 0.14);
	        f_tcc_lp[i]->SetParLimits(0, -1E+8, 0);
	        f_tcc_lp[i]->SetParLimits(1, 0, 10);
	        f_tcc_lp[i]->SetParLimits(2, 0, 10);

		hb_hp[i]->Fit(f_tcc_hp[i],"N","",fitRangeMin,6000);	
		hb_lp[i]->Fit(f_tcc_lp[i],"N","",fitRangeMin,6000);	

		float chi2_tcc_hp = f_tcc_hp[i]->GetChisquare() / f_tcc_hp[i]->GetNDF();
		float chi2_tcc_lp = f_tcc_lp[i]->GetChisquare() / f_tcc_lp[i]->GetNDF();

		// Smoothing
		for(int j=1; j<nbins+1; j++){
			if( rebin21[j] > smoothXmin ){
				float vwidth = rebin21[j] - rebin21[j-1];
				hb_hp_af[i]->SetBinContent(j, f_tcc_hp[i]->Integral(rebin21[j-1],rebin21[j])/vwidth);
				hb_lp_af[i]->SetBinContent(j, f_tcc_lp[i]->Integral(rebin21[j-1],rebin21[j])/vwidth);
				hb_hp_af[i]->SetBinError(j, f_tcc_hp[i]->IntegralError(rebin21[j-1],rebin21[j])/sqrt(vwidth));
				hb_lp_af[i]->SetBinError(j, f_tcc_lp[i]->IntegralError(rebin21[j-1],rebin21[j])/sqrt(vwidth));
				//std::cout << f_tcc_lp[i]->IntegralError(rebin21[j-1],rebin21[j]) << std::endl;
			}else{
				hb_hp_af[i]->SetBinContent(j, hb_hp[i]->GetBinContent(j));
				hb_lp_af[i]->SetBinContent(j, hb_lp[i]->GetBinContent(j));
				hb_hp_af[i]->SetBinError(j, hb_hp[i]->GetBinError(j));
				hb_lp_af[i]->SetBinError(j, hb_lp[i]->GetBinError(j));
			}
		}
	
		// Plot
		lat2 = new TLatex(500, 1E+5*0.1, "#it{#bf{ATLAS}} #bf{#bf{Simulation Work in progress}}");
		lat_tcc_hp = new TLatex(2500, 1E-1, Form("#scale[0.8]{#bf{#bf{TCC,HP (chi2/NDF=%3.2f)}}}",chi2_tcc_hp));
		lat_tcc_lp = new TLatex(2500, 1E-1, Form("#scale[0.8]{#bf{#bf{TCC,LP (chi2/NDF=%3.2f)}}}",chi2_tcc_lp));
		cv->SetLogy();
		cv->cd();
        	leg2= new TLegend(0.45,0.61,0.90,0.85);                                                    
        	leg2->SetBorderSize(0);
		leg2->SetFillStyle(0);
		hb_hp_af[i]->GetYaxis()->SetRangeUser(1E-7,1E+5);
		hb_hp_af[i]->GetXaxis()->SetTitle("m_{T}(vvJ)[GeV]");
		hb_hp_af[i]->SetLineColor(kBlue);
		hb_hp_af[i]->Draw("hist");
		f_tcc_hp[i]->SetLineColor(kRed);
		f_tcc_hp[i]->DrawCopy("same");
		f_tcc_hp[i]->SetFillColor(kRed);
		f_tcc_hp[i]->SetFillStyle(3018);
		f_tcc_hp[i]->SetLineWidth(2);
		f_tcc_hp[i]->SetMarkerSize(0);
		f_tcc_hp[i]->Draw("E3 same");
		hb_hp[i]->Draw("E0 same");
		leg2->AddEntry(hb_hp[i],   "Pre-fit  Bkg","lep");
		leg2->AddEntry(f_tcc_hp[i],    "Post-fit Bkg function","l");
		leg2->AddEntry(hb_hp_af[i],"Post-fit Bkg (Binned)","l");
		if( tagger=="W" ){
			hs_hp[i][sigMap[302265]]->SetLineColor(kGreen);
			hs_hp[i][sigMap[302265]]->Draw("hist same");
			leg2->AddEntry(hs_hp[i][sigMap[302265]],"HVTWZ 5TeV","l");
		}else{
			hs_hp[i][sigMap[303326]]->SetLineColor(kGreen);
			hs_hp[i][sigMap[303326]]->Draw("hist same");
			leg2->AddEntry(hs_hp[i][sigMap[303326]],"RSG_ZZ 5TeV","l");
		}
		leg2->Draw();
		lat2->Draw();
		lat_tcc_hp->Draw();
		cv->Print(Form("../output/smoothing.tcc.hp.%dthWP.",i)+ch+"."+tagger+"."+sigType+"."+nameOption+"."+Form("%d",Ncombs)+"."+outNameOpt+".png");
		delete leg2;

        	leg2= new TLegend(0.45,0.61,0.90,0.85);                       
        	leg2->SetBorderSize(0);
		leg2->SetFillStyle(0);
		hb_lp_af[i]->GetYaxis()->SetRangeUser(1E-7,1E+5);
		hb_lp_af[i]->GetXaxis()->SetTitle("m_{T}(vvJ)[GeV]");
		hb_lp_af[i]->SetLineColor(kBlue);
		hb_lp_af[i]->Draw("hist");
		f_tcc_lp[i]->SetLineColor(kRed);
		f_tcc_lp[i]->DrawCopy("same");
		f_tcc_lp[i]->SetFillColor(kRed);
		f_tcc_lp[i]->SetFillStyle(3018);
		f_tcc_lp[i]->SetLineWidth(2);
		f_tcc_lp[i]->SetMarkerSize(0);
		f_tcc_lp[i]->Draw("E3 same");
		hb_lp[i]->Draw("E0 same");
		leg2->AddEntry(hb_lp[i],   "Pre-fit  Bkg","lep");
		leg2->AddEntry(f_tcc_lp[i],    "Post-fit Bkg function","l");
		leg2->AddEntry(hb_lp_af[i],"Post-fit Bkg (Binned)","l");
		if( tagger=="W" ){
			hs_lp[i][sigMap[302265]]->SetLineColor(kGreen);
			hs_lp[i][sigMap[302265]]->Draw("hist same");
			leg2->AddEntry(hs_lp[i][sigMap[302265]],"HVTWZ 5TeV","l");
		}else{
			hs_lp[i][sigMap[303326]]->SetLineColor(kGreen);
			hs_lp[i][sigMap[303326]]->Draw("hist same");
			leg2->AddEntry(hs_lp[i][sigMap[303326]],"RSG_ZZ 5TeV","l");
		}
		leg2->Draw();
		lat2->Draw();
		lat_tcc_lp->Draw();
		cv->Print(Form("../output/smoothing.tcc.lp.%dthWP.",i)+ch+"."+tagger+"."+sigType+"."+nameOption+"."+Form("%d",Ncombs)+"."+outNameOpt+".png");
	}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Significance calculation
	std::cout << "Significance calculation....." << std::endl;
	float signi_hp[Ncombs][NmassPoints];
	float error_hp[Ncombs][NmassPoints];
	float ratio_hp[Ncombs][NmassPoints];
	float ratioE_hp[Ncombs][NmassPoints];
	float signi_lp[Ncombs][NmassPoints];
	float error_lp[Ncombs][NmassPoints];
	float ratio_lp[Ncombs][NmassPoints];
	float ratioE_lp[Ncombs][NmassPoints];

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
			if( isMCSmooth ){
				signi_hp[i][m] = significance(hs_hp[i][m], hs_hp[i][m], hb_hp_af[i], "mass");
				error_hp[i][m] = significance(hs_hp[i][m], hs_hp[i][m], hb_hp_af[i], "mass", "error");
				signi_lp[i][m] = significance(hs_lp[i][m], hs_lp[i][m], hb_lp_af[i], "mass");
				error_lp[i][m] = significance(hs_lp[i][m], hs_lp[i][m], hb_lp_af[i], "mass", "error");
			}else{
				signi_hp[i][m] = significance(hs_hp[i][m], hs_hp[i][m], hb_hp[i], "mass");
				error_hp[i][m] = significance(hs_hp[i][m], hs_hp[i][m], hb_hp[i], "mass", "error");
				signi_lp[i][m] = significance(hs_lp[i][m], hs_lp[i][m], hb_lp[i], "mass");
				error_lp[i][m] = significance(hs_lp[i][m], hs_lp[i][m], hb_lp[i], "mass", "error");
			}

			signi[i][m] = sqrt( pow(signi_hp[i][m],2) + pow(signi_lp[i][m],2) );
			//signi[i][m] = signi_hp[i][m];

			error[i][m] = sqrt( pow(error_hp[i][m],2) + pow(error_lp[i][m],2) );
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
        TLegend *leg= new TLegend(0.25,0.46,0.90,0.84);                                                    
        leg->SetBorderSize(0);
	leg->SetFillStyle(0);

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
			float v1  = hg_lct->GetBinContent(hg_lct->FindBin(stoi(sigName[sigArr[m]])));
			float er1 = hg_lct->GetBinError(hg_lct->FindBin(stoi(sigName[sigArr[m]])));

			float r = v0/v1;
			float rer = sqrt( pow(er0/v1,2) + pow(v0/v1/v1*er1,2) );
			if( isnan(r) || isinf(r) ) r = -100.;
			if( isnan(rer) || isinf(rer) ) rer = 0.;
			hgr[i]->SetBinContent(hgr[i]->FindBin(stoi(sigName[sigArr[m]])), r);
			hgr[i]->SetBinError(hgr[i]->FindBin(stoi(sigName[sigArr[m]])), rer);
			ratio[i][m] = r;
			error[i][m] = rer;
		}
		i++;
	}

	pad1->cd();
	//pad1->SetGridy();
	TLatex *lat = new TLatex(500, maxSig*yScaleMax*0.92, "#it{#bf{ATLAS}} #bf{#bf{Simulation Work in progress}}");
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
	//leg->AddEntry(hg_lct, "LCTopo HP(50%)+LP(80%)", "lp");
	leg->AddEntry(hg_lct, "LCTopo", "lp");
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
	lat->Draw();

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


