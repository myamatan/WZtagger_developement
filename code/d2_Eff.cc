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
#include "TLegend.h"
#include "TCut.h"
#include "TEventList.h"
#include "TGaxis.h"

#include "/nfs_home/home1/myamatan/RootUtils/AtlasStyle.C"
#include "sigList.hh"
#include "func.hh"

TString  boson = "W";
TString  chs = "1Lep";
TString  memo = "roughBin";
int nbin = 15;

std::vector<std::vector<float>> ptRanges = {
	{500,600},
	{600,700},
	{700,3000}
};
const UInt_t NptRange = ptRanges.size();

float mcSum=0.;

void d2_Eff(TString cutDef="loose", TString WP="MaxVVJJ_3D"){

	SetAtlasStyle();

	TH1F *hbef[NptRange];
	TH1F *haft[NptRange];
	TH1F *heff[NptRange];
	TH1F *hr[NptRange];

	TF1 *fD2;
	TF1 *fML;
	TF1 *fMU;

	TFile *f;
	TTree *tr;

	//TString wpName = "MaxVVSemi_v5";
	//TString wpName = "MaxVVJJ_3D";
	//TString wpName = "MaxVVJJ";
	TString wpName = WP;

	// Fill
	//TString fname = "../input/merged_new/allMC/allMC." + chs + ".root";
	nameOption = "Wjets";
	TString fname = "../input/merged_new/"+nameOption+"/"+nameOption+"." + chs + ".root";


	std::cout << "++++++ " << fname << " ++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

	f = new TFile(fname);
	tr = (TTree*)f->Get("Nominal");
		
	SetBranch(tr, chs);

	for(UInt_t i=0; i<NptRange; i++){
		hbef[i] = new TH1F(Form("hb%i",i),  Form("hb%i",i),     nbin,50,200);
		haft[i] = new TH1F(Form("ha%i",i),  Form("ha%i",i),     nbin,50,200);
		heff[i] = new TH1F(Form("heff%i",i),Form("heff%i",i),   nbin,50,200);
		hr[i]   = new TH1F(Form("hr%i",i),Form("hr%i",i),   nbin,50,200);
	}
			
	// Smoothed cut parameters
	std::string wpfile = std::string("../output/WPs/SmoothedContained"+boson+wps[wpName]);
	std::vector<TString> cutWPs = retrieveWPs( wpfile );
	fD2 = new TF1("D2"+boson, cutWPs.at(2), 200, 2500);
	fMU = new TF1("MU"+boson, cutWPs.at(1), 200, 2500);
	fML = new TF1("ML"+boson, cutWPs.at(0), 200, 2500);

	float D2cut=0.;
	float MUcut=0.;
	float MLcut=0.;
	TCut precut;
	if( cutDef!="loose") precut = "fatjet_D2>0 && fatjet_m >=50 &&" + preCutfor1Lep;
	else precut = "fatjet_D2>0 && fatjet_m >=50 &&" + preCutfor1Lep_loose;
	tr->Draw(">>elist",precut);
	TEventList *elist = (TEventList*)gROOT->FindObject("elist");
	long nEvent = elist->GetN();
	for(long l=0; l<nEvent; l++){
		tr->GetEntry(elist->GetEntry(l));
		// D2 & Mass window 
		D2cut = pt > 2500 ? fD2->Eval(2500) : fD2->Eval(pt);
		MUcut = pt > 2500 ? fMU->Eval(2500) : fMU->Eval(pt);
		MLcut = pt > 2500 ? fML->Eval(2500) : fML->Eval(pt);
		for(UInt_t i=0; i<NptRange; i++){
		  if( pt > ptRanges.at(i).at(0) && pt <= ptRanges.at(i).at(1) ){
		    hbef[i]->Fill( m, weight );
		    if( d2<D2cut ) haft[i]->Fill( m, weight );
		  }
		}
		if( l%(int(nEvent/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEvent, float(float(l)/nEvent)*100) << std::endl;
	}

	// Efficiency calculation
	for(UInt_t i=0; i<NptRange; i++){
		setEfficiency(hbef[i], haft[i], heff[i]);
		setRatio(heff[i], heff[0], hr[i]);
	}


	// TCanvas
        TCanvas *cv;
        TPad *pad1;
        TPad *pad2;
        TLegend *leg;
	TLatex    *lat;
	TLatex *latsub;
        std::cout << "Preparing canvas................";
        cv = new TCanvas("c1","Canvas",200,10,500,480);//600,480
        pad1 = new TPad("pad1","This is pad1",0.01,0.19,0.95,0.97);//0.2
        pad2 = new TPad("pad2","This is pad2",0.01,0.05,0.95,0.313);//0.05
        pad1->SetFillColor(kWhite);
        pad2->SetFillColor(kWhite);
        pad1->Draw();
        pad2->Draw();
        std::cout << "done" << std::endl;
	gStyle->SetEndErrorSize(1.8);

	lat = new TLatex(55, 1.1, "#it{#bf{ATLAS}} #scale[0.8]{#it{Simulation Work in progress}}");
	latsub = new TLatex(55, 0.93, "#scale[0.8]{#splitline{"+boson+" tagger, D2 cut efficiency}{}}");

	std::cout << "Plotting Pad1...................";
        pad1->cd();
        //pad1->SetLogy();
        leg = new TLegend(0.6,0.56,0.88,0.85);
        float k_Nbin  = float(heff[0]->GetNbinsX()+1);
        float k_range = float(heff[0]->GetXaxis()->GetXmax() - heff[0]->GetXaxis()->GetXmin());
        float k_xmin  = float(heff[0]->GetXaxis()->GetXmin()) + k_range/k_Nbin*float(heff[0]->FindFirstBinAbove());
        float k_xmax  = float(heff[0]->GetXaxis()->GetXmin()) + k_range/k_Nbin*float(heff[0]->FindLastBinAbove());
        float XrangeMax = k_xmax;
        float XrangeMin = k_xmin;
        float k_ymax = heff[0]->GetBinContent( heff[0]->GetMaximumBin() );
        float YrangeMax = k_ymax * 1.1;//40.1
        float YrangeMin = k_ymax * 0;//1000.

	leg->SetBorderSize(0);
        for(UInt_t i=0; i<NptRange; i++){
                heff[i]->SetStats(0);
                heff[i]->SetLabelSize(0);
                //heff[i]->GetYaxis()->SetRangeUser(YrangeMin, YrangeMax);
                heff[i]->GetYaxis()->SetRangeUser(0., 1.2);
                heff[i]->GetYaxis()->SetTitle("Efficiency");
                ((TGaxis*)heff[i]->GetYaxis())->SetMaxDigits(4);
                heff[i]->SetLineColor(colorList.at(i));
                heff[i]->SetMarkerColor(colorList.at(i));
                if( i==0 ) heff[i]->Draw("E1");
                else heff[i]->Draw("E1 same");
		leg->AddEntry(heff[i], Form("%d < p_{T}(J) < %d GeV",int(ptRanges.at(i).at(0)),int(ptRanges.at(i).at(1))), "lep");
        }

	// Draw leg and latex
        leg->Draw();
	lat->Draw();
	latsub->Draw();
        pad1->RedrawAxis();
        std::cout << "done" << std::endl;

        std::cout << "Plotting Pad2...................";
        pad2->cd();
        pad2->SetBottomMargin(0.25);
        pad2->SetGridy();
        for(UInt_t j=1; j<NptRange; j++){
                hr[j]->SetStats(0);
                hr[j]->GetYaxis()->SetTitle("Ratio");
                hr[j]->GetXaxis()->SetTitle("m_{J} [GeV]");
                hr[j]->SetLineColor(colorList.at(j));
                hr[j]->SetMarkerColor(colorList.at(j));
                hr[j]->SetMarkerSize(0.8);
                hr[j]->GetYaxis()->SetRangeUser(0.4,1.6);
                hr[j]->GetYaxis()->SetNdivisions(505);
                hr[j]->GetXaxis()->SetTitleOffset(0.88);
                hr[j]->GetYaxis()->SetTitleOffset(0.47);
                hr[j]->GetXaxis()->SetTitleSize(0.13);
                hr[j]->GetYaxis()->SetTitleSize(0.13);
                hr[j]->GetXaxis()->SetLabelSize(0.15);
                hr[j]->GetYaxis()->SetLabelSize(0.15);
                if( j==0 ) hr[j]->Draw("E1");
                else hr[j]->Draw("same E1");
        }
        std::cout << "done" << std::endl;

	cv->Print("../output/d2EffOnMass_ratio."+cutDef+"."+boson + "." + nameOption + "." + wpName + "." + memo + ".eps");

/*

	// Plot
	TCanvas    *cv;
	TLegend   *leg;
	TLatex    *lat;
	TLatex *latsub;
	cv = new TCanvas("c1","Canvas",200,10,500,480);
	leg= new TLegend(0.6,0.56,0.88,0.85);
	leg->SetBorderSize(0);
	lat = new TLatex(55, 1.1, "#it{#bf{ATLAS}} #scale[0.8]{#it{Simulation Work in progress}}");
	latsub = new TLatex(55, 0.93, "#scale[0.8]{#splitline{"+boson+" tagger, D2 cut efficiency}{}}");
	cv->cd();
	for(UInt_t i=0; i<NptRange; i++){
	  if( i==0 ) heff[i]->Draw();
	  else heff[i]->Draw("same");
	  heff[i]->SetStats(0);
	  heff[i]->GetYaxis()->SetRangeUser(0,1.2);
	  heff[i]->GetYaxis()->SetTitle("Efficiency");
	  heff[i]->GetXaxis()->SetTitle("m_{J} [GeV]");
	  heff[i]->SetMarkerColor(colorList.at(i));
	  heff[i]->SetMarkerSize(1);
	  heff[i]->SetMarkerStyle(8);
	  heff[i]->SetLineColor(colorList.at(i));
	  heff[i]->SetLineWidth(2);
	  leg->AddEntry(heff[i], Form("%d < p_{T}(J) < %d GeV",int(ptRanges.at(i).at(0)),int(ptRanges.at(i).at(1))), "lep");
	}
	leg->Draw();
	lat->Draw();
	latsub->Draw();
	cv->Print("../output/d2EffOnMass." + boson + "." + nameOption + "." + wpName + ".eps");
*/

}
