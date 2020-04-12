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
#include "THStack.h"

#include "/nfs_home/home1/myamatan/RootUtils/AtlasStyle.C"
#include "func.hh"
#include "sigList.hh"

std::vector<TString>  filename = {

				   "allMC",
				   "Zjets",
				   "Wjets",
				   "ttbar",
				   "diboson",
				   //"signal",
				   "RS_G_ZZm0700",
				   "RS_G_ZZm3000",

				   //"RS_G_ZZm0700",
				   //"allMC",
				   //"RS_G_ZZm3000",
				   //"RS_G_ZZm5000",


				   //"HVTm0700",
				   //"allMC",
				   //"HVTm3000",
				   //"HVTm5000",

				   //"HVTWWm0700",
				   //"allMC",
				   //"HVTWWm3000",
				   //"HVTWWm5000",
};
TString stackList = "Zjets.Wjets.diboson.ttbar";

float mcSum=0.;

void compare(){

	UInt_t nfile = filename.size();
	TH1F *hD2[nfile];
	TH1F *hM[nfile];
	THStack *hD2s = new THStack();
	THStack *hMs = new THStack();

	TFile *f;
	TTree *tr;

	TLatex *latD2 = new TLatex(400, 2.7, "#it{ATLAS} #bf{Simulation Preliminary}");
	TLatex *latMa = new TLatex(400, 140, "#it{ATLAS} #bf{Simulation Preliminary}");

	for(UInt_t i=0; i<nfile; i++){

		std::cout << filename.at(i) << " +++++++++++++++++++++++++++++++++++++++" << std::endl;

		TString sample = filename.at(i);
		TString fname;
		if( nameOption == "Model" ){
			if( sample.Index("00") != -1 ) fname = "../input/merged/signal"+ tagger + "/" + sample + ".root";
			else fname = "../input/merged/" + sample(0,sample.Index(".")) + tagger  + "/" + sample + tagger + ".root";
		}else{
			if( sample.Index("00") != -1 ) fname = "../input/merged/signal"+ tagger + "/" + sample + "." + ch+ ".root";
			else if( sample=="signal" ) fname = "../input/merged/signal"+ tagger + "/signal" + tagger + "." + ch+ ".root";
			else fname = "../input/merged/" + sample + "/" + sample + "." + ch+ ".root";
		}

		f = new TFile(fname);
		tr = (TTree*)f->Get("Nominal");

		hD2[i]= new TH1F(Form("hD2%i",i),Form("hD2%i",i),40, 0,5);
		hM[i]= new TH1F(Form("hM%i",i),Form("hM%i",i),40,50,150);

		SetBranch(tr, ch);

		float D2cut;
		float MUcut;
		float MLcut;

		TCut precut = returnTCut("beforeBosonTag");
        	tr->Draw(">>elist",precut);
        	TEventList *elist = (TEventList*)gROOT->FindObject("elist");
        	long nEvent = elist->GetN();//long nEvent = tr->GetEntries();
		for(long l=0; l<nEvent; l++){
			tr->GetEntry(elist->GetEntry(l));//tr->GetEntry(l);
			if( 1 /* ztag>=0*/ ){

				if( sample=="allMC" ) mcSum += weight;

				float weightD2 = weight;
				float weightM = weight;
				if( sigWeight[mc]!=0 ){
					weightD2 = weight * 40 * sigWeight[mc] * 15; //30 //60
					weightM = weight * 40 * sigWeight[mc] * 15;  //15 //30
				}
				// D2 
				//D2cut = pt > 2500 ? fD2->Eval(2500) : fD2->Eval(pt);
				//if( d2 < D2cut ){
					hD2[i]->Fill(d2, weightD2);
				//}

				// Mass cut
				//MUcut = pt > 2500 ? fMU->Eval(2500) : fMU->Eval(pt);
				//MLcut = pt > 2500 ? fML->Eval(2500) : fML->Eval(pt);
				//if( m < MUcut && m > MLcut ){
					hM[i]->Fill(m, weightM);
				//}

			}
			if( l%(int(nEvent/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEvent, float(float(l)/nEvent)*100) << std::endl;
		}
	}
	
	TCanvas *c1 = new TCanvas();
        TCanvas *c2 = new TCanvas();
	
	TLegend *leg= new TLegend(0.6,0.6,0.88,0.88);
	leg->SetBorderSize(0);

	c1->cd();
	k_ymax = 5E+3;
	//c1->SetLogy();
	c1->SetGridx();
	bool isFirstDraw = true;
	std::vector<int> signalID;
	for(UInt_t i=0; i<nfile; i++){
		hD2[i]->SetStats(0);
		hD2[i]->SetTitle("");
		hD2[i]->GetXaxis()->SetRangeUser(0,5);
		hD2[i]->GetXaxis()->SetTitle("fatJet D2");
 		// Signal scaling
		if( TString(filename.at(i)).Index("00")!=-1 || TString(filename.at(i)).Index("signal")!=-1 ){
			hD2[i]->Scale(mcSum/2/hD2[i]->Integral());
			k_ymax = hD2[i]->GetBinContent( hD2[i]->GetMaximumBin() );
			signalID.push_back(i);
		}
		hD2[i]->GetYaxis()->SetRangeUser(1E-2, k_ymax*3.);
		hD2[i]->SetLineColor(colorList.at(i));
		hD2[i]->SetLineWidth(2);
		if( stackList.Index(filename.at(i)) != -1 ){
			hD2[i]->SetLineColor(kBlack);
			hD2[i]->SetFillColor(colorList.at(i));
			hD2s->Add(hD2[i]);
		};
		if( isFirstDraw && stackList.Index(filename.at(i)) == -1 ){
			hD2[i]->Draw("HIST");
			isFirstDraw = false;
		}
		else if( stackList.Index(filename.at(i)) == -1 ){
			hD2[i]->Draw("HIST same");
		}
		if( filename.at(i) != "allMC" ) leg->AddEntry(hD2[i], filename.at(i), "f");
	}
	hD2s->Draw("HIST same");
	for(UInt_t i=0; i<signalID.size(); i++){
		hD2[signalID.at(i)]->Draw("HIST same");
	}
	leg->Draw("same");
	latD2->Draw("same");
	c1->RedrawAxis();

	c2->cd();
	isFirstDraw = true;
	//c2->SetLogy();
	k_ymax = 5E+3;
	c2->SetGridx();
	for(UInt_t i=0; i<nfile; i++){
		hM[i]->SetStats(0);
		hM[i]->SetTitle("");
		hM[i]->GetXaxis()->SetRangeUser(50,150);
		hM[i]->GetXaxis()->SetTitle("fatJet M [GeV]");
		if( TString(filename.at(i)).Index("00")!=-1 || TString(filename.at(i)).Index("signal")!=-1 ){
			hM[i]->Scale(mcSum/4/hM[i]->Integral());
			k_ymax = hM[i]->GetBinContent( hM[i]->GetMaximumBin() );
		}
		hM[i]->GetYaxis()->SetRangeUser(1E-2, k_ymax*1.2);
		hM[i]->SetLineColor(colorList.at(i));
		hM[i]->SetLineWidth(2);
		if( stackList.Index(filename.at(i)) != -1 ){
			hM[i]->SetLineColor(kBlack);
			hM[i]->SetFillColor(colorList.at(i));
			hMs->Add(hM[i]);
		};
		if( isFirstDraw && stackList.Index(filename.at(i)) == -1 ){
			hM[i]->Draw("HIST");
			isFirstDraw = false;
		}
		else if( stackList.Index(filename.at(i)) == -1 ){
			hM[i]->Draw("HIST same");
		}
	}
	hMs->Draw("HIST same");
	for(UInt_t i=0; i<signalID.size(); i++){
		hM[signalID.at(i)]->Draw("HIST same");
	}
	leg->Draw("same");
	latMa->Draw("same");
	c2->RedrawAxis();

	c1->Print("../output/D2."+ch+"."+tagger+"."+nameOption+".png");
	c2->Print("../output/Mass."+ch+"."+tagger+"."+nameOption+".png");

}
