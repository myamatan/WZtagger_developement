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

std::vector<std::vector<float>> ptRanges = {
	{500,600},
	{600,700},
	{700,3000}
};
const UInt_t NptRange = ptRanges.size();

float mcSum=0.;

void mShape_bypT(TString cutDef="loose", TString WP="MaxVVSemi_v5"){

	SetAtlasStyle();

	TH1F *h[NptRange];
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
		h[i] = new TH1F(Form("h%i",i),  Form("h%i",i),     30,50,200);
		hr[i] = new TH1F(Form("hr%i",i),  Form("hr%i",i),  30,50,200);
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
	//TCut precut = returnTCut("beforeBosonTag", chs);
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
		for(int i=0; i<NptRange; i++){
		  if( pt > ptRanges.at(i).at(0) && pt <= ptRanges.at(i).at(1) ){
		    if( d2<D2cut ) h[i]->Fill( m, weight );
		  }
		}
		if( l%(int(nEvent/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEvent, float(float(l)/nEvent)*100) << std::endl;
	}

	// Ratio calculation
	for(UInt_t i=0; i<NptRange; i++){
		h[i]->Scale(1./h[i]->Integral());
	}
	for(UInt_t i=0; i<NptRange; i++){
		setRatio(h[0], h[i], hr[i]);
	}

	// TCanvas
        TCanvas *cv;
        TPad *pad1;
        TPad *pad2;
        TLegend *leg;
        TLatex *latex;
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

	std::cout << "Plotting Pad1...................";
        pad1->cd();
        //pad1->SetLogy();
        leg = new TLegend(0.5 ,0.2 ,0.92,0.48);
        float k_Nbin  = float(h[0]->GetNbinsX()+1);
        float k_range = float(h[0]->GetXaxis()->GetXmax() - h[0]->GetXaxis()->GetXmin());
        float k_xmin  = float(h[0]->GetXaxis()->GetXmin()) + k_range/k_Nbin*float(h[0]->FindFirstBinAbove());
        float k_xmax  = float(h[0]->GetXaxis()->GetXmin()) + k_range/k_Nbin*float(h[0]->FindLastBinAbove());
        float XrangeMax = k_xmax;
        float XrangeMin = k_xmin;
        float k_ymax = h[0]->GetBinContent( h[0]->GetMaximumBin() );
        float YrangeMax = k_ymax * 1.1;//40.1
        float YrangeMin = k_ymax * 0;//1000.
        latex = new TLatex(XrangeMin+(XrangeMax-XrangeMin)*0.35,
                        YrangeMin+(YrangeMax-YrangeMin)*0.45,
                        "#splitline{ #it{#bf{ATLAS}} #it{#scale[0.8]{Simulation Internal}} }{ #scale[0.8]{#sqrt{s} = 13 TeV, "+boson+" tagger} }");

	leg->SetBorderSize(0);
        for(int i=0; i<NptRange; i++){
                h[i]->SetStats(0);
                h[i]->SetLabelSize(0);
                h[i]->GetYaxis()->SetRangeUser(YrangeMin, YrangeMax);
                h[i]->GetYaxis()->SetTitle("normalized #event");
                ((TGaxis*)h[i]->GetYaxis())->SetMaxDigits(4);
                h[i]->SetLineColor(colorList.at(i));
                h[i]->SetMarkerColor(colorList.at(i));
                if( i==0 ) h[i]->Draw("E1");
                else h[i]->Draw("E1 same");
		leg->AddEntry(h[i], Form("%d < p_{T}(J) < %d GeV",int(ptRanges.at(i).at(0)),int(ptRanges.at(i).at(1))), "lep");
        }

	// Draw leg and latex
        leg->Draw();
        latex->Draw();
        pad1->RedrawAxis();
        std::cout << "done" << std::endl;

        std::cout << "Plotting Pad2...................";
        pad2->cd();
        pad2->SetBottomMargin(0.25);
        pad2->SetGridy();
        for(int j=1; j<NptRange; j++){
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

	cv->Print("../output/mShape_byPt." + cutDef + "." + boson + "." + nameOption + "." + wpName + ".eps");

}
