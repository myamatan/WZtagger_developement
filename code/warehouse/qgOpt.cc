#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TString.h"
#include "TCut.h"
#include "TEventList.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "THStack.h"

#include "/nfs_home/home1/myamatan/RootUtils/AtlasStyle.C"
#include "func.hh"
#include "sigList.hh"

TString filenameSig = "../input/merged_truthStu/signal"+tagger+"/signal"+tagger+"."+ch+".root";
TString filenameBkg = "../input/merged_truthStu/Zjets/Zjets."+ch+".root";

// QG_ratio
float n_quark = 20007.3;
float n_gluon = 3062.39;
float qratio_max = 1.0;
float qratio_min = 0.0;
float qratio_d = 0.2;

// D2 & mass cut scan range
float d2cut_max = 2.0;
float d2cut_min = 0.5;
float d2cut_d = 0.002;

float mcutu_max = 111;//120
float mcutu_min = 110;//90
float mcutu_d = 1;//2

float mcutl_max = 71;//90
float mcutl_min = 70;//50
float mcutl_d = 1;

float xMin = qratio_min; // 450
float xMax = qratio_max;// 5450
float yScaleMin[3] = {0.7, 0.7, 0.7};
float yScaleMax[3] = {1.6, 1.6, 1.6};


void qgOpt(){

	SetAtlasStyle();

	TFile *fS = new TFile(filenameSig);
	TTree *trS = (TTree*)fS->Get("Nominal");
	TFile *fB= new TFile(filenameBkg);
	TTree *trB = (TTree*)fB->Get("Nominal");

	// D2 & MU & ML scan range
	UInt_t d2_range = UInt_t( (d2cut_max-d2cut_min)/d2cut_d );
	UInt_t mu_range = UInt_t( (mcutu_max-mcutu_min)/mcutu_d );
	UInt_t ml_range = UInt_t( (mcutl_max-mcutl_min)/mcutl_d );

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
	//std::string wpfile = std::string("../output/WPs/SmoothedContained"+tagger+wps["FixVVJJ"]);
	std::vector<TString> cutWPs = retrieveWPs( wpfile );
	fD2 = new TF1(Form("D2"), cutWPs.at(2), 200, 2500);
	fMU = new TF1(Form("MU"), cutWPs.at(1), 200, 2500);
	fML = new TF1(Form("ML"), cutWPs.at(0), 200, 2500);

	// Histo initialization (x:VV mass, y:fatjet pt)
	TString nameID;
	TH1F *hbkg[d2_range][mu_range][ml_range][q_range];
	TH1F *hsig[d2_range][mu_range][ml_range][NmassPoints];
	for(UInt_t i=0; i<d2_range; i++){
		for(UInt_t j=0; j<mu_range; j++){
			for(UInt_t k=0; k<ml_range; k++){
				for(UInt_t q=0; q<q_range; q++){
					nameID = TString(std::to_string(i)+"-"+std::to_string(j)+"-"+std::to_string(k)+"-"+std::to_string(q)+"-");
					hbkg[i][j][k][q] = new TH1F(nameID+"hbkg", nameID+"hbkg", nbins, rebin21);
				}
				for(UInt_t q=0; q<NmassPoints; q++){
					nameID = TString(std::to_string(i)+"-"+std::to_string(j)+"-"+std::to_string(k)+"-"+std::to_string(q)+"-");
					hsig[i][j][k][q] = new TH1F(nameID+"hsig", nameID+"hsig", nbins, rebin21);
				}
			}
		}
	}

	TH1F *h[NmassPoints][3];
	for(UInt_t i=0; i<NmassPoints; i++){
		for(UInt_t j=0; j<3; j++){
			nameID = TString(std::to_string(i)+"-"+std::to_string(j));
			h[i][j] = new TH1F(nameID+"hsig", nameID+"hsig", q_range, qratio_min, qratio_max);
		}
	}

	/////////////////////////////
	// Background
	/////////////////////////////
	std::cout << "Background filling ....." << std::endl;
	SetBranch(trB, ch);
	float qsum  = 0.;
	float gsum  = 0.;
	// Filling
	TCut precut = returnTCut(preSele);
        trB->Draw(">>elist",precut);
        TEventList *elistB = (TEventList*)gROOT->FindObject("elist");
	long nEventB = elistB->GetN();
	for(long l=0; l<nEventB; l++){
        	trB->GetEntry(elistB->GetEntry(l));
		// D2 &&  Mass window
		d2cut = pt > 2500 ? fD2->Eval(2500) : fD2->Eval(pt);
		mcutu = pt > 2500 ? fMU->Eval(2500) : fMU->Eval(pt);
		mcutl = pt > 2500 ? fML->Eval(2500) : fML->Eval(pt);
		// Raw filling
		if( pdgid==21 ) qsum += weight*lumi/40.;
		else gsum += weight*lumi/40.;
		// D2 &&  Mass window
		for(UInt_t i=0; i<d2_range; i++){
		  for(UInt_t j=0; j<mu_range; j++){
		    for(UInt_t k=0; k<ml_range; k++){

		    	if( d2 < d2cut_min+d2cut_d*i &&
			    m < mcutu_max-mcutu_d*j &&
			    m > mcutl_min+mcutl_d*k &&
			    //d2 < d2cut &&
			    //m < mcutu &&
			    //m > mcutl &&
			    1 
			){
				for(UInt_t m=0; m<q_range; m++){
					if( pdgid==21 ) hbkg[i][j][k][m]->Fill(vvm, weight*lumi/40. * g_scale * (1.-(qratio_min+qratio_d*m+qratio_d/2)));
					else hbkg[i][j][k][m]->Fill(vvm, weight*lumi/40. * q_scale * (qratio_min+qratio_d*m+qratio_d/2));
				}
			}

		    }
		  } 
		}
		if( l%(int(nEventB/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEventB, float(float(l)/nEventB)*100) << std::endl;
	}

	std::cout << "Bkg q/g : " <<  qsum << " / " << gsum << std::endl;

	/////////////////////////////
	// Signal
	/////////////////////////////
	std::cout << "Signal filling ....." << std::endl;
	std::map<int, int> sigMap;
	UInt_t sigCount = 0;
	UInt_t sigArr[10];
	SetBranch(trS, ch);
	precut = returnTCut(preSele);
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
		float weightS = weight;
		if( sigWeight.find(mc) != sigWeight.end() ){
			weightS = weight * lumi * sigWeight[mc];
			if( sigMap.find(mc) == sigMap.end() ){
				sigMap[mc] = sigCount; 
				sigArr[sigCount] = mc;
				sigCount++;
			}
		}
		// D2 &&  Mass window
		for(UInt_t i=0; i<d2_range; i++){
		  for(UInt_t j=0; j<mu_range; j++){
		    for(UInt_t k=0; k<ml_range; k++){

		    	if( d2 < d2cut_min+d2cut_d*i &&
			    m < mcutu_max-mcutu_d*j &&
			    m > mcutl_min+mcutl_d*k &&
			    //d2 < d2cut &&
			    //m < mcutu &&
			    //m > mcutl &&
			    1
			){
				hsig[i][j][k][sigMap[mc]]->Fill(vvm, weightS);
			}

		    }
		  } 
		}
		if( l%(int(nEventS/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEventS, float(float(l)/nEventS)*100) << std::endl;
	}

	// Significance calculation
	std::cout << "Optimized cut searching  ....." << std::endl;
	float cut[q_range][NmassPoints][3];
	float maxSigni=0.0;
	float signi=0.0;
	std::vector<int> massList;
	for(UInt_t i=0; i<q_range; i++){
		for(UInt_t m=0; m<sigCount; m++){
			// D2 & Mass window
			for(UInt_t j=0; j<d2_range; j++){
				for(UInt_t k=0; k<mu_range; k++){
					for(UInt_t l=0; l<ml_range; l++){
						//float s = hsig[j][k][l][m]->Integral(0, nbins+1);
						//float b = hbkg[j][k][l][i]->Integral(0, nbins+1);
						signi = significance(hsig[j][k][l][m], hsig[j][k][l][m], hbkg[j][k][l][i], "mass");
					    	//signi = 2*( (s+b)*log(1+s/b)-s );
						if( isnan(signi) || isinf(signi) || signi<0 ) signi = 0.;
						if( signi > maxSigni){
							cut[i][m][0] = d2cut_min+d2cut_d*j;
							cut[i][m][1] = mcutu_max-mcutu_d*k;
							cut[i][m][2] = mcutl_min+mcutl_d*l;
							maxSigni = signi;
						}
					}
				}
			}
			maxSigni=0.;
			signi=0.;
			if( i==0 ) massList.push_back(stoi(sigName[sigArr[m]]));
		}
	}
	std::sort(massList.begin(), massList.end());
	std::map<int, int> orderMap;
	std::map<int, int> op_orderMap;
	for(UInt_t m=0; m<sigCount; m++){
		for(UInt_t l=0; l<massList.size(); l++){
			if( stoi(sigName[sigArr[m]]) == massList.at(l) ){
				orderMap[m] = l; // massList --> mass order
				op_orderMap[l] = m; // mass order --> massList
			}
		}
	}

	// Filling
	float maxCon[3]={0};
	float minCon[3]={1E+9,1E+9,1E+9};
	std::cout << "Mass	D2 cut	Mucut	Mlcut" << std::endl;
	for(UInt_t i=0; i<q_range; i++){
		for(UInt_t m=0; m<sigCount; m++){
			for(UInt_t s=0; s<3; s++){
				float vtofill = cut[i][m][s];
				float ver = d2cut_d;
				if( s==1 ) ver = mcutu_d;
				else if( s==2 ) ver = mcutl_d;
				h[m][s]->SetBinContent(i,vtofill);
				h[m][s]->SetBinError(i,ver);
				if( vtofill > maxCon[s] ) maxCon[s] = vtofill;
				if( vtofill < minCon[s] ) minCon[s] = vtofill;
			}
			if( m==0 ) std::cout << stoi(sigName[sigArr[m]]) << "	" << cut[i][m][0] << "	" << cut[i][m][1] << "	" << cut[i][m][2] << std::endl;
		}
	}
	
	std::cout << "Range ++++++++++++++++++" << std::endl;
	for(UInt_t s=0; s<3; s++){
		std::cout << maxCon[s] << " : " << minCon[s] << std::endl;
	}

	// Plottintg .......................................
	std::cout << "Plotting +++++++++++++++++" << std::endl;
	std::map<int, bool> plotList ={
			{700, false},
			{1000, true},
			{2000, true},
			{3000, true},
			{4000, false},
			{5000, false}
	};
	bool isfirstDraw = true; 
	
	TCanvas *c;
        TLegend *leg;
	TLatex *lat;
	TLatex *latsub;
	for(int i=0; i<1; i++){
	//for(int i=0; i<3; i++){
		c = new TCanvas(Form("c%d",i),Form("c%d",i),600,500); //cv = new TCanvas("c1","Canvas",200,10,600,480);//600,480
		lat = new TLatex(0.05, 2.4, "#it{#bf{ATLAS}} #bf{#bf{Simulation Preliminary}}");
		latsub = new TLatex(0.06, 2.25, "#scale[0.8]{#bf{#bf{#sqrt{s}=13 TeV, 70 < m(J) < 110}}}");
		c->cd();
		//c->SetGridy();
        	leg= new TLegend(0.18,0.58,0.58,0.76);                    
        	leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		for(int m=0; m<int(sigCount); m++){
			int j = op_orderMap[m];
			h[j][i]->SetStats(0);
			h[j][i]->SetTitle("");
        	        h[j][i]->GetYaxis()->SetTitleSize(0.04);
        	        h[j][i]->GetYaxis()->SetTitleOffset(1.22);
			h[j][i]->GetXaxis()->SetTitle("Fraction of quark-induced jet");
			h[j][i]->GetXaxis()->SetRangeUser(xMin,xMax);
			h[j][i]->GetYaxis()->SetTitle("Optimal D2 threshold");
			h[j][i]->GetYaxis()->SetRangeUser(minCon[i]*yScaleMin[i],maxCon[i]*yScaleMax[i]);
			h[j][i]->SetLineColor(colorListSig.at(j));
			h[j][i]->SetLineWidth(2);
			h[j][i]->SetMarkerColor(colorListSig.at(j));
			h[j][i]->SetMarkerStyle(8);
			h[j][i]->SetMarkerSize(1.1);
			if( plotList[stoi(sigName[sigArr[j]])] ){
				leg->AddEntry(h[j][i], "RS_G "+TString(sigName[sigArr[j]])+" GeV", "lp");
				if( !isfirstDraw ) h[j][i]->Draw("E0 same");
				else h[j][i]->Draw("E0");
				isfirstDraw = false;
			}
		}
		leg->Draw();
		lat->Draw();
		latsub->Draw();
		c->Print("../output/qgOpt."+ch+"."+tagger+Form(".%d",i)+".pdf");

		//delete c;
		//delete leg;
	}

}

