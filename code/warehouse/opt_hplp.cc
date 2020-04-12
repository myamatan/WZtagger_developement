#include<iostream>
#include<fstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TString.h"
#include "TCut.h"
#include "TEventList.h"
#include "TF1.h"

#include "/nfs_home/home1/myamatan/RootUtils/AtlasStyle.C"
#include "sigList.hh"
#include "func.hh"

TString filenameSig = "../input/merged_new/signal"+tagger+"/signal"+tagger+"."+ch+".root";
TString filenameBkg = "../input/merged_new/allMC/allMC."+ch+".root";

float d2cut_max = 1.7;//1.8
float d2cut_min = 1.0;//1.0
float d2cut_d = 0.05;//0.05

float mcutu_max = 117;//117, 125
float mcutu_min = 95;//95 , 100 
float mcutu_d = 2;//2

float mcutl_max = 76;//76  , 81
float mcutl_min = 58;//58  , 67
float mcutl_d = 2;//2

void opt_hplp(){

	TFile *fS = new TFile(filenameSig);
	TTree *trS = (TTree*)fS->Get("Nominal");
	TFile *fB= new TFile(filenameBkg);
	TTree *trB = (TTree*)fB->Get("Nominal");

	// Optiized value output
	ofstream outf(std::string("../output/cutValue.") + ch.Data() + "." + tagger + "." + sigType + "." + nameOption + std::string(".txt"));

	TFile *f = new TFile("../output/cutValue." + ch + "." + tagger + "." + sigType + "." + nameOption + ".root","recreate");
	TTree *tr = new TTree("col","col");
	Float_t D2, MU, ML, SIG, PT, S, B, eB, eS;
	tr->Branch("PT",&PT);
	tr->Branch("D2",&D2);
	tr->Branch("MU",&MU);
	tr->Branch("ML",&ML);
	tr->Branch("SIG",&SIG);
	tr->Branch("S",&S);
	tr->Branch("B",&B);
	tr->Branch("eS",&eS);
	tr->Branch("eB",&eB);
	

	// Cut Fcuntion
	TF1 *fD2;
	TF1 *fML;
	TF1 *fMU;
	std::string wpfile = std::string("../output/WPs/SmoothedContained"+tagger+wps["MaxVVSemi"]);
	std::vector<TString> cutWPs = retrieveWPs( wpfile );
	fD2 = new TF1(Form("D2"), cutWPs.at(2), 200, 2500);
	fMU = new TF1(Form("MU"), cutWPs.at(1), 200, 2500);
	fML = new TF1(Form("ML"), cutWPs.at(0), 200, 2500);
	
	UInt_t d2_range = UInt_t( (d2cut_max-d2cut_min)/d2cut_d );
	UInt_t mu_range = UInt_t( (mcutu_max-mcutu_min)/mcutu_d );
	UInt_t ml_range = UInt_t( (mcutl_max-mcutl_min)/mcutl_d );

	// D2 & Mass window
	TH2F *hb[d2_range][mu_range][ml_range][2];
	TH2F *hs[d2_range][mu_range][ml_range][NmassPoints][2];
	TH2F *hsAll[d2_range][mu_range][ml_range][2];
	TH2F *hb_LP[d2_range][mu_range][ml_range][2];
	TH2F *hs_LP[d2_range][mu_range][ml_range][NmassPoints][2];
	TH2F *hsAll_LP[d2_range][mu_range][ml_range][2];


	// Pt flat
	TH1F *hpt[NmassPoints];
	for(UInt_t l=0; l<NmassPoints; l++){
		TString nameID = TString(std::to_string(l));
		hpt[l] = new TH1F(nameID+"h", nameID+"h", nPtbins, ptRange);
	}
	TH1F *hptAll= new TH1F("hptAll", "hptAll", nPtbins, ptRange);
	TH1F *hptb= new TH1F("hptb", "hptb", nPtbins, ptRange);

	// Significance fill
	TH3F *h3D[nPtbins];
	for(int i=0; i<nPtbins; i++){
		h3D[i] = new TH3F(Form("signi_%d",i),Form("signi_%d",i),d2_range,d2cut_min,d2cut_max,mu_range,mcutu_min,mcutu_max,ml_range,mcutl_min,mcutl_max);
	}

	// Histo initialization (x:VV mass, y:fatjet pt)
	for(UInt_t i=0; i<d2_range; i++){
		for(UInt_t j=0; j<mu_range; j++){
			for(UInt_t k=0; k<ml_range; k++){
				TString nameID = TString(std::to_string(i) + "-" + std::to_string(j) + "-" + std::to_string(k));
				hb[i][j][k][0] = new TH2F(nameID+"hbb", nameID+"hbb",nbins, rebin21, nPtbins, ptRange);
				hb[i][j][k][1] = new TH2F(nameID+"hba", nameID+"hba",nbins, rebin21, nPtbins, ptRange);
				hb_LP[i][j][k][0] = new TH2F(nameID+"hbb_LP", nameID+"hbb_LP",nbins, rebin21, nPtbins, ptRange);
				hb_LP[i][j][k][1] = new TH2F(nameID+"hba_LP", nameID+"hba_LP",nbins, rebin21, nPtbins, ptRange);
				for(UInt_t l=0; l<NmassPoints; l++){
					nameID = TString(std::to_string(i) + "-" + std::to_string(j) + "-" + std::to_string(k) + "-" + std::to_string(l));
					hs[i][j][k][l][0] = new TH2F(nameID+"hsb", nameID+"hsb", nbins, rebin21, nPtbins, ptRange);
					hs[i][j][k][l][1] = new TH2F(nameID+"hsa", nameID+"hsa", nbins, rebin21, nPtbins, ptRange);
					hs_LP[i][j][k][l][0] = new TH2F(nameID+"hsb_LP", nameID+"hsb_LP", nbins, rebin21, nPtbins, ptRange);
					hs_LP[i][j][k][l][1] = new TH2F(nameID+"hsa_LP", nameID+"hsa_LP", nbins, rebin21, nPtbins, ptRange);
				}
				hsAll[i][j][k][0] = new TH2F(nameID+"hsAllb", nameID+"hsAllb", nbins, rebin21, nPtbins, ptRange);
				hsAll[i][j][k][1] = new TH2F(nameID+"hsAlla", nameID+"hsAlla", nbins, rebin21, nPtbins, ptRange);
				hsAll_LP[i][j][k][0] = new TH2F(nameID+"hsAllb_LP", nameID+"hsAllb_LP", nbins, rebin21, nPtbins, ptRange);
				hsAll_LP[i][j][k][1] = new TH2F(nameID+"hsAlla_LP", nameID+"hsAlla_LP", nbins, rebin21, nPtbins, ptRange);
			}
		}
	}

	/////////////////////////////
	// Background
	/////////////////////////////
	std::cout << "Background filling ....." << std::endl;
	SetBranch(trB, ch);
	// Filling
	TCut precut = returnTCut(preSele);
        trB->Draw(">>elist",precut);
        TEventList *elistB = (TEventList*)gROOT->FindObject("elist");
	long nEventB = elistB->GetN(); //long nEventB = trB->GetEntries();
	for(long l=0; l<nEventB; l++){
        	trB->GetEntry(elistB->GetEntry(l)); //trB->GetEntry(l);
		// D2 &&  Mass window
		for(UInt_t i=0; i<d2_range; i++){
			for(UInt_t j=0; j<mu_range; j++){
				for(UInt_t k=0; k<ml_range; k++){
					hb[i][j][k][0]->Fill(vvm, pt, weight*lumi/40.);
					hb_LP[i][j][k][0]->Fill(vvm, pt, weight*lumi/40.);
					if( d2 < d2cut_max-d2cut_d*i &&
					    m  < mcutu_max-mcutu_d*j &&
			 		    m  > mcutl_min+mcutl_d*k &&
					    //m  < fMU->Eval(pt) &&
					    //m  > fML->Eval(pt) &&
					    1
					){
						hb[i][j][k][1]->Fill(vvm, pt, weight*lumi/40.);
					}else if( //d2 >= d2cut_max-d2cut_d*i &&
					          m  < mcutu_max-mcutu_d*j &&
			 		          m  > mcutl_min+mcutl_d*k &&
					          //m  < fMU->Eval(pt) &&
					          //m  > fML->Eval(pt) &&
					          1
					){
						hb_LP[i][j][k][1]->Fill(vvm, pt, weight*lumi/40.);
					}
				}
			}
		}
		if( l%(int(nEventB/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEventB, float(float(l)/nEventB)*100) << std::endl;
	}

	/////////////////////////////
	// Signal
	/////////////////////////////
	std::cout << "Signal filling ....." << std::endl;
	std::map<int, int> sigMap;
	UInt_t sigCount = 0;
	UInt_t sigArr[30];
	SetBranch(trS, ch);
	precut = returnTCut(preSele);
        trS->Draw(">>elist",precut);
        TEventList *elistS = (TEventList*)gROOT->FindObject("elist");
	long nEventS = elistS->GetN(); //long nEventS = trS->GetEntries();
	// For pt flat
	for(long l=0; l<nEventS; l++){
        	trS->GetEntry(elistS->GetEntry(l)); //trS->GetEntry(l);
		float weightS = weight;
		if( sigWeight.find(dsid) != sigWeight.end() ){
			if( ch!="1Lep" ) weightS = weight * lumi * sigWeight[dsid];
			//if( ch=="1Lep" && sigType=="eff" ) weightS = weight / sigWeight[dsid];
		}
		if( sigMap.find(dsid) == sigMap.end() ){
			sigMap[dsid] = sigCount; 
			sigArr[sigCount] = dsid;
			sigCount++;
		}
		hpt[sigMap[dsid]]->Fill(pt, weightS);
		hptAll->Fill(pt, weightS);
	}
	// Filling
	for(long l=0; l<nEventS; l++){
        	trS->GetEntry(elistS->GetEntry(l)); //trS->GetEntry(l);
		float weightS = weight;
		if( sigWeight.find(dsid) != sigWeight.end() ){
			if( ch!="1Lep" ) weightS = weight * lumi * sigWeight[dsid];
			//if( ch=="1Lep" && sigType=="eff" ) weightS = weight / sigWeight[dsid];
		}
		if( isPtFlat ) weightS *= ptFlat(hptAll, pt);
		// D2 & Mass window
		for(UInt_t i=0; i<d2_range; i++){
			for(UInt_t j=0; j<mu_range; j++){
				for(UInt_t k=0; k<ml_range; k++){
					hs[i][j][k][sigMap[dsid]][0]->Fill(vvm, pt, weightS);
					hsAll[i][j][k][0]->Fill(vvm, pt, weightS);
					hs_LP[i][j][k][sigMap[dsid]][0]->Fill(vvm, pt, weightS);
					hsAll_LP[i][j][k][0]->Fill(vvm, pt, weightS);
					if( d2 < d2cut_max-d2cut_d*i &&
					    m  < mcutu_max-mcutu_d*j &&
			 		    m  > mcutl_min+mcutl_d*k &&
					    //m < fMU->Eval(pt) &&
					    //m > fML->Eval(pt) &&
					    1
					){  
						hs[i][j][k][sigMap[dsid]][1]->Fill(vvm, pt, weightS);
						hsAll[i][j][k][1]->Fill(vvm, pt, weightS);
					}else if( //d2 >= d2cut_max-d2cut_d*i &&
					          m  < mcutu_max-mcutu_d*j &&
			 		          m  > mcutl_min+mcutl_d*k &&
					          //m < fMU->Eval(pt) &&
					          //m > fML->Eval(pt) &&
					          1
					){  
						hs_LP[i][j][k][sigMap[dsid]][1]->Fill(vvm, pt, weightS);
						hsAll_LP[i][j][k][1]->Fill(vvm, pt, weightS);
					}
				}
			}
		}
		if( l%(int(nEventS/10.))==0 ) std::cout << Form("[%li/%li] : %0.0f %% done .....", l, nEventS, float(float(l)/nEventS)*100) << std::endl;
	}

	// Significance calculation
	std::cout << "Optimized cut searching  ....." << std::endl;
	float d2cut[nPtbins]={0};
	float mcutu[nPtbins]={0};
	float mcutl[nPtbins]={0};
	float maxSigni=-999;
	float signi=0.0;
	float pre_signi=0.0;
	bool  isFill=false;
	for(UInt_t i=0; i<nPtbins; i++){
		// D2 & Mass window
		std::cout << "Pt : " << ptRange[i] << "+++++" << std::endl;
		for(UInt_t j=0; j<d2_range; j++){
			for(UInt_t k=0; k<mu_range; k++){
				for(UInt_t l=0; l<ml_range; l++){
					signi = 0.;
				    	//signi = significance(hsAll[j][k][l][0], hsAll[j][k][l][1], hb[j][k][l][1], i+1, sigType);
					//if( isnan(signi) || signi<=0 ) signi = 0.0;
					float sumS=0.0;
					float sumB=0.0;
					float sumSe=0.0;
					float sumBe=0.0;
		       		 	for(UInt_t m=0; m<sigCount; m++){
						float pre_sig = pow(significance(hs[j][k][l][m][0], hs[j][k][l][m][1], hb[j][k][l][1], i+1, sigType), 2);
						float pre_sig_LP = pow(significance(hs_LP[j][k][l][m][0], hs_LP[j][k][l][m][1], hb_LP[j][k][l][1], i+1, sigType), 2);
						if( !isnan(pre_sig) && !isnan(pre_sig_LP) && pre_sig>0 && pre_sig_LP>0 ) signi += sqrt(pow(pre_sig,2)+pow(pre_sig_LP,2));
						for(int s=1; s<=hs[j][k][l][m][1]->GetNbinsX(); s++){
	        				        sumS += hs[j][k][l][m][1]->GetBinContent(s, i+1);
	        				        sumSe += pow(hs[j][k][l][m][1]->GetBinError(s, i+1), 2);
	        				        if( m==0 ){
								sumB += hb[j][k][l][1]->GetBinContent(s, i+1);
								sumBe += pow(hb[j][k][l][1]->GetBinError(s, i+1), 2);
							}
						}
					}
					h3D[i]->SetBinContent(j,k,l,signi);// x:D2,y:mu,z;ml
					PT  = ptRange[i];
					D2  = d2cut_max-d2cut_d*j;
					MU  = mcutu_max-mcutu_d*k;
					ML  = mcutl_min+mcutl_d*l;
					SIG = signi;
					S   = sumS;
					B   = sumB;
					eS  = sqrt(sumSe);
					eB  = sqrt(sumBe);
					tr->Fill();
					if( (!isFill && signi>maxSigni) || (isFill && signi/maxSigni>1.0005 && signi/pre_signi<1.2) ){ // 1.005
						d2cut[i] = d2cut_max-d2cut_d*j;
						mcutu[i] = mcutu_max-mcutu_d*k;
						mcutl[i] = mcutl_min+mcutl_d*l;
						maxSigni = signi;
						isFill   = true;
					}
					pre_signi = signi;
				}
			}
		}
		h3D[i]->Write();
		maxSigni=-999.;
		signi=0.;
		isFill=false;
	}
	tr->Write();

	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	std::cout << "Pt	:bin width	:D2	:Upper Mass	: Lower Mass" << std::endl;	
	for(UInt_t i=0; i<nPtbins; i++){
		std::cout << ptRange[i] << "	" << int(ptRange[i+1]-ptRange[i]) << "	" << d2cut[i] << "	" << mcutu[i] << "	" << mcutl[i] << std::endl;
		if( i==nPtbins-1 ) outf << int((ptRange[i]+ptRange[i+1])/2.) << "," << int(ptRange[i+1]-ptRange[i]) << "," << d2cut[i] << "," << mcutu[i] << "," << mcutl[i];
		else outf << int((ptRange[i]+ptRange[i+1])/2.) << "," << int(ptRange[i+1]-ptRange[i]) << "," << d2cut[i] << "," << mcutu[i] << "," << mcutl[i] << ",";
		//pyoutf << int((ptRange[i]+ptRange[i+1])/2.) << "," << d2cut[i] << "," << mcutu[i] << "," << mcutl[i] << std::endl;
	}
	outf.close();

}
