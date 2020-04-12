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

#include "/nfs_home/home1/myamatan/RootUtils/AtlasStyle.C"
#include "func.hh"
#include "sigList.hh"

TString filenameSig = "../input/merged/signal"+tagger+"/signal"+tagger+"."+ch+".root";
TString filenameBkg = "../input/merged/Zjets/Zjets."+ch+".root";

float d2cut_min = 0.5;
float d2cut_max = 2.0;
float d2cut_d = 0.02;

float mcutu_min = 92.;//92
float mcutu_max = 126.;//126
float mcutu_d = 2;//2

float mcutl_min = 54.;//54
float mcutl_max = 86.;//86
float mcutl_d = 2;

float n_quark = 2619.2;
float n_gluon = 1612.11;
float qratio_max = 1.0;
float qratio_min = 0.0;
float qratio_d = 0.05;

float xMin = qratio_min; // 450
float xMax = qratio_max;// 5450
float yScaleMin = 0.8;   // 0.7
float yScaleMax = 1.5;   // 1.6

float yMin[3] = {0.7, 0.8, 0.8};
float yMax[3] = {1.6, 2., 2.};

void simpleOpt(){

	TFile *fS = new TFile(filenameSig);
	TTree *trS = (TTree*)fS->Get("Nominal");
	TFile *fB= new TFile(filenameBkg);
	TTree *trB = (TTree*)fB->Get("Nominal");

	// Optiized value output
	ofstream outf(std::string("../output/cutValue.") + ch.Data() + "." + tagger + "." + sigType + "." + nameOption + std::string(".txt"));
	
	// D2 & MU & ML scan range
	UInt_t d2_range = UInt_t( (d2cut_max-d2cut_min)/d2cut_d );
	UInt_t mu_range = UInt_t( (mcutu_max-mcutu_min)/mcutu_d );
	UInt_t ml_range = UInt_t( (mcutl_max-mcutl_min)/mcutl_d );

	// quarl ratio range
	qratio_max += qratio_d;
	UInt_t q_range = UInt_t( (qratio_max-qratio_min)/qratio_d );
	float q_scale = (n_gluon+n_quark) / n_quark;
	float g_scale = (n_gluon+n_quark) / n_gluon;

	// D2 & Mass window
	float bkg[d2_range][mu_range][ml_range][q_range];
	float sig[d2_range][mu_range][ml_range][NmassPoints];
	
	// Histo initialization (x:VV mass, y:fatjet pt)
	TString nameID;
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
	float qsumAll = 0.;
	float gsumAll = 0.;
	float qsum[q_range];
	float gsum [q_range];
	std::cout << "Background filling ....." << std::endl;
	SetBranch(trB, ch);
	// Filling
	TCut precut = returnTCut(preSele);
        trB->Draw(">>elist",precut);
        TEventList *elistB = (TEventList*)gROOT->FindObject("elist");
	long nEventB = elistB->GetN(); //long nEventB = trB->GetEntries();
	for(long l=0; l<nEventB; l++){
        	trB->GetEntry(elistB->GetEntry(l)); //trB->GetEntry(l);
		if( pdgid==21 ) gsumAll += weight;
		else qsumAll += weight;
		// D2 &&  Mass window
		for(UInt_t i=0; i<d2_range; i++){
			for(UInt_t j=0; j<mu_range; j++){
				for(UInt_t k=0; k<ml_range; k++){

					if( d2 < d2cut_min+d2cut_d*i && m < mcutu_max-mcutu_d*j && m > mcutl_min+mcutl_d*k ){
						for(UInt_t m=0; m<q_range; m++){
							if( pdgid==21 ) bkg[i][j][k][m] += weight*lumi/40. * g_scale * (1.-(qratio_min+qratio_d*m));
							else bkg[i][j][k][m] += weight*lumi/40. * q_scale * (qratio_min+qratio_d*m);
						}
					}

				}
			}
		}
		for(UInt_t m=0; m<q_range; m++){
			if( pdgid==21 ) gsum[m] += weight*lumi/40. * g_scale * (1.-(qratio_min+qratio_d*m));
			else qsum[m] += weight*lumi/40. * q_scale * (qratio_min+qratio_d*m);;
		}
		if( l%(int(nEventB/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEventB, float(float(l)/nEventB)*100) << std::endl;
	}
	std::cout << "Bkg q/g : " << qsumAll << " / " << gsumAll << std::endl;
	
	std::cout << "q(g)_ratio	q	g	q+g" << std::endl;
	for(UInt_t m=0; m<q_range; m++){
		std::cout << qratio_min+qratio_d*m << "		" << qsum[m] << "	" << gsum[m] << "	" << qsum[m]+gsum[m]<< std::endl;
	}

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
	long nEventS = elistS->GetN(); //long nEventS = trS->GetEntries();
	// Filling
	for(long l=0; l<nEventS; l++){
        	trS->GetEntry(elistS->GetEntry(l)); //trS->GetEntry(l);
		float weightS = weight;
		if( sigWeight.find(mc) != sigWeight.end() ){
			//weightS = weight * lumi * sigWeight[mc];
			if( sigMap.find(mc) == sigMap.end() ){
				sigMap[mc] = sigCount; 
				sigArr[sigCount] = mc;
				sigCount++;
			}
		}
		// D2 & Mass window
		for(UInt_t i=0; i<d2_range; i++){
			for(UInt_t j=0; j<mu_range; j++){
				for(UInt_t k=0; k<ml_range; k++){
					if( d2 < d2cut_min+d2cut_d*i && m < mcutu_max-mcutu_d*j && m > mcutl_min+mcutl_d*k){
						sig[i][j][k][sigMap[mc]] += weightS;
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
						signi = 0.;
						float s = sig[j][k][l][m];
						float b = bkg[j][k][l][i];
					    	signi = 2*( (s+b)*log(1+s/b)-s );
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
			if( i==0 ) massList.push_back(stoi(sigName[sigArr[m]]));
			maxSigni=0.;
			signi=0.;
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
	float maxSig[3]={0};
	float minSig[3]={1E+9,1E+9,1E+9};
	for(UInt_t i=0; i<q_range; i++){
		for(UInt_t m=0; m<sigCount; m++){
			for(UInt_t s=0; s<3; s++){
				float vtofill = cut[i][m][s];
				float ver = d2cut_d;
				if( s==1 ) ver = mcutu_d;
				else if( s==2 ) ver = mcutl_d;
				h[m][s]->SetBinContent(i,vtofill);
				h[m][s]->SetBinError(i,ver);
				if( vtofill > maxSig[s] ) maxSig[s] = vtofill;
				if( vtofill < minSig[s] ) minSig[s] = vtofill;
			}
		}
	}

	// Plottintg .......................................
	TCanvas *c[3];
        TLegend *leg[3];           
	for(int i=0; i<3; i++){
		c[i] = new TCanvas(Form("c%d",i),Form("c%d",i),600,600); //cv = new TCanvas("c1","Canvas",200,10,600,480);//600,480
		c[i]->cd();
		c[i]->SetGridy();
        	leg[i]= new TLegend(0.5,0.65,0.88,0.88);                    
        	leg[i]->SetBorderSize(0);
		for(int m=0; m<int(sigCount); m++){
			int j = op_orderMap[m];
			h[j][i]->SetStats(0);
			h[j][i]->SetTitle("");
        	        h[j][i]->GetYaxis()->SetTitleSize(0.04);
        	        h[j][i]->GetYaxis()->SetTitleOffset(1.22);
			h[j][i]->GetXaxis()->SetTitle("quark ratio (1 - gluon ratio)");
			h[j][i]->GetXaxis()->SetRangeUser(xMin,xMax);
			h[j][i]->GetYaxis()->SetTitle("Optimized cut value");
			//h[j][i]->GetYaxis()->SetRangeUser(minSig[i]*yScaleMin,maxSig[i]*yScaleMax);
			h[j][i]->GetYaxis()->SetRangeUser(minSig[i]*yMin[i],maxSig[i]*yMax[i]);
			h[j][i]->SetLineColor(colorListSig.at(j));
			h[j][i]->SetLineWidth(2);
			h[j][i]->SetMarkerColor(colorListSig.at(j));
			h[j][i]->SetMarkerStyle(8);
			h[j][i]->SetMarkerSize(0.8);
			leg[i]->AddEntry(h[j][i], "Signal : "+TString(sigName[sigArr[j]])+" [GeV]", "lp");
			if( m!=0 ) h[j][i]->Draw("E0 same");
			else h[j][i]->Draw("E0");
		}
		leg[i]->Draw();
		c[i]->Print("../output/simpleComp."+ch+"."+tagger+Form(".%d",i)+".png");
	}

}
