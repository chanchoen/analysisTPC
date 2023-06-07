#include <iostream>
#include <sstream>
#include <string>

#include "GETDecoder.h"
#include "GETPad.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TTree.h"
#include "TPaveText.h"
#include "TLegend.h"

using namespace TMath;
using namespace std;

const double par3 = 2.685;
// const double par4 = 16.01;
const double par4 = 21.19;
const double par1ToADC = TMath::Power(par3 * par4 / TMath::E(), par3);

const double fDriftvelocity = 55; //[mm/us]
const double fPadHeight = 11.9; //[mm]
const double fPadWeith = 1.9; //[mm]
const double fPadGap = 0.1; //[mm]
const double fDistance = 550; //[mm]
const double fhalfH = 89; //[mm]
const double fTimeoffsetRW = 172.629; //[mm]
const double fTimeoffsetLW = -171.667; //[mm]

// double transform_x(double xPos[8]){
//     return xPos[8]-(fPadWeith+fPadGap)*16;
// }

// double transform_y(double yPos[8]){
//     return -yPos[8]+(fDistance+(fPadHeight+fPadGap)*3.5);
// }

// double transform_z(double zPos[8]){
//     return zPos[8]-fhalfH;
// }

double sigFcn(double* x, double* par) {
    return x[0] < par[2] ? par[0] : par[1] * TMath::Power(x[0] - par[2], par3) * TMath::Exp(-(x[0] - par[2]) / par4);
}
double fitFcn(double* x, double* par) {
    return sigFcn(x, par) < par[3] ? sigFcn(x, par) : par[3];
}

int arrSigFcnTR[350], arrSigFcnTL[350];
double sigFcnTR(double* x, double* par) {
    int idx = x[0] - (par[0] - 149);
    if (idx < 0) {
        return par[1] * arrSigFcnTR[0];
    } else if (idx >= 350) {
        return par[1] * arrSigFcnTR[349];
    } else {
        return par[1] * arrSigFcnTR[idx];
    }
}
double sigFcnTL(double* x, double* par) {
    int idx = x[0] - (par[0] - 149);
    if (idx < 0) {
        return par[1] * arrSigFcnTL[0];
    } else if (idx >= 350) {
        return par[1] * arrSigFcnTL[349];
    } else {
        return par[1] * arrSigFcnTR[idx];
    }
}

std::string getLabel(int asadId, int agetId, int chanId);

void createSignalTemplate();
void fitSignalTemplate();
void resolutionTest();
void resolutionTestL();
void readSigFcn();
void ADC_to_Hits();
void resolutionTestLandR();

int main() {
    // readSigFcn();
    // ADC_to_Hits();
    // createSignalTemplate();
    // fitSignalTemplate();
    // resolutionTest();
    // resolutionTestL();
    resolutionTestLandR();
}

// void ADC_to_Hits() {
//     GETDecoder decoder;
//     GETPad pad;

//     auto hSignal = new TH1D("hSignal", "", 500, 0.5, 500.5);
//     auto hPedestal = new TH1D("hPedestal", "", 1000, -499.5, 500.5);
//     auto fFitFcnTR = new TF1("fFitFcnTR", fitFcnTR, 1, 500, 4);
//     auto fFitFcTL = new TF1("fFitFcnTL", sigFcnTL, 1, 500, 4);
//     //
//     int runID;
//     int eventID;
//     bool isBroken;
//     double eventTime;
//     double eventDiffTime;
//     int ADCs[4][4][68][512];  // [nAsAds][nAGETs][nChannels][nTimebuckets]
//     //
//     TFile* fileIn = new TFile("/home/shlee/workspace/treeOfADC_all_2.root", "read");
//     TTree* treeIn = (TTree*)fileIn->Get("treeOfADC");
//     treeIn->SetBranchAddress("RunID", &runID);
//     treeIn->SetBranchAddress("EventID", &eventID);
//     treeIn->SetBranchAddress("IsBroken", &isBroken);
//     treeIn->SetBranchAddress("EventTime", &eventTime);
//     treeIn->SetBranchAddress("EventDiffTime", &eventDiffTime);
//     treeIn->SetBranchAddress("ADC", ADCs);
//     //
//     int E_CsI_RD;
//     int E_CsI_RU;
//     int E_CsI_LD;
//     int E_CsI_LU;
//     //
//     int nHits_TR;
//     int xId_TR[256];
//     int yId_TR[256];
//     double driftTime_TR[256];
//     double ADC_TR[256];
//     //
//     int nHits_TL;
//     int xId_TL[256];
//     int yId_TL[256];
//     double driftTime_TL[256];
//     double ADC_TL[256];
//     //
//     TFile* fout = new TFile("./treeOfHits_newSig.root", "recreate");
//     TTree* treeOut = new TTree("treeOfHits", "");
//     treeOut->Branch("RunID", &runID, "RunID/I");
//     treeOut->Branch("EventID", &eventID, "EventID/I");
//     treeOut->Branch("IsBroken", &isBroken, "IsBroken/O");
//     treeOut->Branch("EventTime", &eventTime, "EventTime/D");
//     treeOut->Branch("EventDiffTime", &eventDiffTime, "EventDiffTime/D");
//     // CsI
//     treeOut->Branch("E_CsI_RD", &E_CsI_RD, "E_CsI_RD/I");
//     treeOut->Branch("E_CsI_RU", &E_CsI_RU, "E_CsI_RU/I");
//     treeOut->Branch("E_CsI_LD", &E_CsI_LD, "E_CsI_LD/I");
//     treeOut->Branch("E_CsI_LU", &E_CsI_LU, "E_CsI_LU/I");
//     // TPC Right
//     treeOut->Branch("nHits_TR", &nHits_TR, "nHits_TR/I");
//     treeOut->Branch("xId_TR", xId_TR, "xId_TR[nHits_TR]/I");
//     treeOut->Branch("yId_TR", yId_TR, "yId_TR[nHits_TR]/I");
//     treeOut->Branch("driftTime_TR", driftTime_TR, "driftTime_TR[nHits_TR]/D");
//     treeOut->Branch("ADC_TR", ADC_TR, "ADC_TR[nHits_TR]/D");
//     // TPC Left
//     treeOut->Branch("nHits_TL", &nHits_TL, "nHits_TL/I");
//     treeOut->Branch("xId_TL", xId_TL, "xId_TL[nHits_TL]/I");
//     treeOut->Branch("yId_TL", yId_TL, "yId_TL[nHits_TL]/I");
//     treeOut->Branch("driftTime_TL", driftTime_TL, "driftTime_TL[nHits_TL]/D");
//     treeOut->Branch("ADC_TL", ADC_TL, "ADC_TL[nHits_TL]/D");

//     for (int idx = 0; idx < treeIn->GetEntries(); idx++) {
//         // for (int idx = 0; idx < treeIn->GetEntries(); idx++) {
//         treeIn->GetEntry(idx);
//         if (idx % 1000 == 0) std::cout << idx << std::endl;
//         // // Init
//         E_CsI_RD = 0;
//         E_CsI_RU = 0;
//         E_CsI_LD = 0;
//         E_CsI_LU = 0;
//         nHits_TR = 0;
//         nHits_TL = 0;
//         memset(xId_TR, 0, sizeof(xId_TR));
//         memset(yId_TR, 0, sizeof(yId_TR));
//         memset(driftTime_TR, 0, sizeof(driftTime_TR));
//         memset(ADC_TR, 0, sizeof(ADC_TR));
//         memset(ADC_TR, 0, sizeof(ADC_TR));

//         memset(xId_TL, 0, sizeof(xId_TL));
//         memset(yId_TL, 0, sizeof(yId_TL));
//         memset(driftTime_TL, 0, sizeof(driftTime_TL));
//         memset(ADC_TL, 0, sizeof(ADC_TL));
//         memset(ADC_TL, 0, sizeof(ADC_TL));
//         if (isBroken) {
//             treeOut->Fill();
//             continue;
//         }
//         // CsI
//         for (int asadId = 1; asadId <= 3; asadId += 2) {
//             // rd, ru, ld, lu
//             int agetId = 3;
//             int chanIds[2][2] = {{49, 15}, {66, 32}};
//             for (int chanId : chanIds[(asadId - 1) / 2]) {
//                 std::string label = getLabel(asadId, agetId, chanId);
//                 int FPNChanId = decoder.FPNChanId(chanId);
//                 for (int buckId = 1; buckId <= 500; buckId++) {
//                     int ADC = ADCs[asadId][agetId][chanId][buckId];
//                     int FPN = ADCs[asadId][agetId][FPNChanId][buckId];
//                     if (buckId <= 140) hPedestal->Fill(ADC - FPN);
//                     hSignal->SetBinContent(buckId, ADC - FPN);
//                 }
//                 int ped = hPedestal->GetBinCenter(hPedestal->GetMaximumBin());
//                 int maxBin = hSignal->GetMaximumBin();
//                 int maxADC = hSignal->GetMaximum() - ped;
//                 hPedestal->Reset("ICESM");
//                 hSignal->Reset("ICESM");
//                 if (maxBin < 185 || maxBin > 200) continue;
//                 if (label == "crd" && maxADC > E_CsI_RD) E_CsI_RD = maxADC;
//                 if (label == "cru" && maxADC > E_CsI_RU) E_CsI_RU = maxADC;
//                 if (label == "cld" && maxADC > E_CsI_LD) E_CsI_LD = maxADC;
//                 if (label == "clu" && maxADC > E_CsI_LU) E_CsI_LU = maxADC;
//             }
//         }
//         // TPC
//         for (int asadId = 0; asadId <= 2; asadId += 2) {
//             for (int agetId = 0; agetId < 4; agetId++) {
//                 for (int chanId = 0; chanId < 68; chanId++) {
//                     if (decoder.IsFPN(chanId) || decoder.IsDead(chanId)) continue;
//                     int FPNChanId = decoder.FPNChanId(chanId);
//                     for (int buckId = 1; buckId <= 400; buckId++) {
//                         int ADC = ADCs[asadId][agetId][chanId][buckId] - ADCs[asadId][agetId][FPNChanId][buckId];
//                         if (buckId <= 150) {
//                             hPedestal->Fill(ADC);
//                         } else {
//                             hSignal->SetBinContent(buckId, ADC);
//                         }
//                     }
//                     int pedADC = hPedestal->GetBinCenter(hPedestal->GetMaximumBin());
//                     int maxBin = hSignal->GetMaximumBin();
//                     int maxADC = hSignal->GetMaximum() - pedADC;
//                     hSignal->Reset("ICESM");
//                     hPedestal->Reset("ICESM");
//                     if (150 < maxBin && maxBin <= 400 && 50 < maxADC) {
//                         for (int buckId = 1; buckId <= 500; buckId++) {
//                             int ADC = ADCs[asadId][agetId][chanId][buckId] - ADCs[asadId][agetId][FPNChanId][buckId];
//                             hSignal->SetBinContent(buckId, ADC - pedADC);
//                         }
//                         if (maxADC > 3400) {
//                             fFitFcn->SetParameters(0, maxADC / par1ToADC, maxBin - par3 * par4, maxADC + 50);
//                             fFitFcn->SetParLimits(1, (maxADC - 5) / par1ToADC, (maxADC + 5000) / par1ToADC);
//                             // fFitFcn->SetParLimits(2, maxBin - 100, maxBin);
//                         } else {
//                             fFitFcn->SetParameters(0, maxADC / par1ToADC, maxBin - par3 * par4, maxADC + 50);
//                             fFitFcn->SetParLimits(1, (maxADC - 5) / par1ToADC, (maxADC + 5) / par1ToADC);
//                             // fFitFcn->SetParLimits(2, maxBin - par3 * par4 - 10, maxBin - par3 * par4 + 10);
//                         }
//                         fFitFcn->SetParLimits(0, -10, 10);
//                         fFitFcn->SetParLimits(3, maxADC - 50, maxADC + 50);
//                         fFitFcn->SetLineColor(kRed);
//                         fFitFcn->SetLineWidth(2);
//                         fFitFcn->SetRange(maxBin - 150, maxBin + 100);
//                         hSignal->Fit(fFitFcn, "QNR");
//                         double pars[4];
//                         fFitFcn->GetParameters(pars);
//                         if (asadId == 0) {
//                             xId_TR[nHits_TR] = pad.GetXId(agetId, chanId);
//                             yId_TR[nHits_TR] = pad.GetYId(agetId, chanId);
//                             driftTime_TR[nHits_TR] = pars[2];
//                             ADC_TR[nHits_TR] = pars[1] * par1ToADC;
//                             nHits_TR++;
//                         } else {
//                             xId_TL[nHits_TL] = pad.GetXId(agetId, chanId);
//                             yId_TL[nHits_TL] = pad.GetYId(agetId, chanId);
//                             driftTime_TL[nHits_TL] = pars[2];
//                             ADC_TL[nHits_TL] = pars[1] * par1ToADC;
//                             nHits_TL++;
//                         }
//                     }
//                     hSignal->Reset("ICESM");
//                 }
//             }
//         }
//         treeOut->Fill();
        
//     }
//     treeOut->Write(0, TObject::kOverwrite);
//     fout->Close();
//     fileIn->Close();
// }

void resolutionTest() {
    GETPad pad;

    TCanvas* cPad;
    gStyle->SetPadTopMargin(0.02);
    gStyle->SetPadBottomMargin(0.085);
    gStyle->SetPadLeftMargin(0.095);
    gStyle->SetPadRightMargin(0.02);
    cPad = new TCanvas("cPad", "", 1500, 750);
    cPad->Divide(2,1);
    bool isPaletteAxisOn = false;

    auto c1 = new TCanvas("c1", "", 1500, 750);
    auto c2 = new TCanvas("c2", "", 1500, 750);
    c1->Divide(2,1);
    c2->Divide(2,1);
    auto c3 = new TCanvas("c3", "", 1500, 1500);
    c3->Divide(2,2);
    auto c4 = new TCanvas("c4", "", 1500, 750);
    c4->Divide(2,1);

    auto fGaus = new TF1("fGaus", "gaus", 0, 31);               // The cluster of x_pad direction
    auto fTrackYX = new TF1("fTrackYX", "pol1", -6.25, 93.75);  // x_pad vs y_pad
    auto fTrackYZ = new TF1("fTrackYZ", "pol1", -6.25, 93.75);  // L_drift(z) vs y_pad
    auto fTrack_thetaYX = new TF1("fTrack_thetaYX", "gaus", -10, 10);
    auto hYZ_transform_dx = new TH1D("hYZ_transform_dx-DT", ";#it{#Delta X}_{transformed} [mm]", 400, -400, 1000);
    auto hYZ_transform_dz = new TH1D("hYZ_transform_dz-DT", ";#it{#Delta Z}_{transformed} [mm]", 400, -400, 1000);

    auto hYX_transform_dy_pad = new TH1D("hYX_transform_dy-pad", ";#it{#Delta Y}_{transformed} [mm]", 200, -100, 100);
    auto hYX_transform_dz_pad = new TH1D("hYX_transform_dz-pad", ";#it{#Delta Z}_{transformed} [mm]", 500, -1000, 1000);
    auto hTarget_dxdy = new TH2D("hTarget_dxdy", "target_(Z=0 plane); #it{#Delta X} [mm]; #it{#Delta Y} [mm]", 400, -100, 300, 200, -100, 100);
    auto hdeltaZ_diff = new TH1D("hdeltaZ_diff", "#it{#Delta Z_{drift}- #Delta Z_{pad}}; diff [mm]", 500, -1000, 1000);

    auto hResoX_Pad = new TH1D("hResoX_Pad", ";#it{#Deltax}_{Pad} [mm];", 1000, -1, 1);   // x_pad resoution
    auto hResoL_Drift = new TH1D("hResoL_Drift", ";#it{#DeltaL}_{Drift} [mm];", 100, -4, 4);  // L_drift(z) resoution
    auto hThetaYX_Pad = new TH1D("hThetaYX_Pad", ";#it{#theta}_{YX} [#circ];", 100, -5, 5);  // 
    auto hThetaYL_Drift = new TH1D("hThetaYL_Drift", ";#it{#theta}_{YL} [#circ];", 30, -15, 15);  // 
    auto hConstantYX_Pad = new TH1D("hConstantYX_Pad", ";#it{#Delta C}_{Pad} [mm];", 100, -100, 100);
    auto hConstantYZ_Drift = new TH1D("hConstantYX_Drift", ";#it{#Delta C}_{Drift} [mm];", 200, 0, 800);
    auto fConstantYX_Pad = new TF1("fConstantYX_Pad", "gaus", -100, 100);
    auto fConstantYZ_Drift = new TF1("fConstantYZ_Drift", "gaus", 300, 400);
    auto height_Drift= new TH2D("height_Drift", "height_drift", 32, 0, 32, 8, 0, 8);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Event Info.
    int eventID;
    bool isBroken;
    // Right TPC
    int nHits_TR;
    int xId_TR[256];
    int yId_TR[256];
    double driftTime_TR[256];
    double ADC_TR[256];
    // int driftTime_TR[256];
    // int ADC_TR[256];
    // Read file
    TFile* fileIn = new TFile("./treeOfHits.root", "read");
    // TFile* fileIn = new TFile("./treeOfHitsOrg.root", "read");
    TTree* treeIn = (TTree*)fileIn->Get("treeOfHits");
    // Event Info.
    treeIn->SetBranchAddress("EventID", &eventID);
    treeIn->SetBranchAddress("IsBroken", &isBroken);
    // TPC Right
    treeIn->SetBranchAddress("nHits_TR", &nHits_TR);
    treeIn->SetBranchAddress("xId_TR", xId_TR);
    treeIn->SetBranchAddress("yId_TR", yId_TR);
    treeIn->SetBranchAddress("driftTime_TR", driftTime_TR);
    treeIn->SetBranchAddress("ADC_TR", ADC_TR);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (int idx = 0;idx < treeIn->GetEntries(); idx++) {
    // for (int idx = 0; idx < treeIn->GetEntries(); idx++) {
        treeIn->GetEntry(idx);
        if (isBroken) continue;
        if (idx % 1000 == 0) std::cout << idx << std::endl;
        // TPC Right
        int adcPad_TR[32][8] = {0};
        int timePad_TR[32][8] = {0};
        bool isSatdPad_TR[32][8] = {0};
        for (int hitIdx = 0; hitIdx < nHits_TR; hitIdx++) {
            int xId = xId_TR[hitIdx];
            int yId = yId_TR[hitIdx];
            adcPad_TR[xId][yId] = ADC_TR[hitIdx];
            timePad_TR[xId][yId] = driftTime_TR[hitIdx];
        }
        //
        int nClusters = 0;
        double xPos[8], yPos[8], zPos[8];
        // double xPos_TF[8], yPos_TF[8], zPos_TF[8];
        auto trackYX = new TGraph();
        auto trackYZ = new TGraph();
        auto trackXY = new TGraph();  // For pad drawing

        // auto trackYX_TF= new TGraph();
        // auto trackYZ_TF= new TGraph();

        for (int yId = 1; yId < 7; yId++) {
            auto hCluster = new TH1D("hCluster", "", 32, -0.5, 31.5);
            for (int xId = 0; xId < 32; xId++) {
                hCluster->SetBinContent(xId + 1, adcPad_TR[xId][yId]);
                pad.ADC->Fill(2. * xId, 12. * yId, adcPad_TR[xId][yId]);  // For pad drawing
                height_Drift->Fill(xId, yId, timePad_TR[xId][yId]);
                // std::cout<< xId<< " "<<yId<<" "<<timePad_TR[xId][yId]<<std::endl;
            }
            int maxBin = hCluster->GetMaximumBin();
            int maxADC = hCluster->GetMaximum();
            if (maxADC > 3500) {
                hCluster->SetBinContent(maxBin - 1, 0);
                hCluster->SetBinContent(maxBin, 0);
                hCluster->SetBinContent(maxBin + 1, 0);
                maxBin = hCluster->GetMaximumBin();
                maxADC = hCluster->GetMaximum();
            }
            int max_xId = maxBin - 1;
            // Only if cluster exists
            if (maxADC > 50 && 2 <= max_xId && max_xId <= 30) {
                if (hCluster->GetBinContent(maxBin - 1) > 50 && hCluster->GetBinContent(maxBin + 1) > 50) {
                    fGaus->SetRange(max_xId - 2, max_xId + 2);
                    hCluster->Fit("fGaus", "QNR");
                    double pars[3];
                    fGaus->GetParameters(pars);
                    xPos[nClusters] = 2. * pars[1];
                    yPos[nClusters] = 12. * yId;
                    zPos[nClusters] = 1.1 * timePad_TR[max_xId][yId]; //1.1=tb * driftvelocity
                    // xPos_TF = transform_x(xPos[nClusters]);
                    // yPos_TF = transform_x(yPos);
                    // zPos_TF = transform_x(zPos);
                    trackYX->SetPoint(trackYX->GetN(), yPos[nClusters], xPos[nClusters]);
                    trackYZ->SetPoint(trackYZ->GetN(), yPos[nClusters], zPos[nClusters]);
                    trackXY->SetPoint(trackXY->GetN(), xPos[nClusters], yPos[nClusters]);  // For pad drawing

                    // trackYX_TF->SetPoint(trackYX_TF->GetN(), yPos_TF[nClusters], xPos_TF[nClusters]);
                    // trackYZ_TF->SetPoint(trackYZ_TF->GetN(), yPos_TF[nClusters], zPos_TF[nClusters]);
                    nClusters++;
                }
            }
            hCluster->Delete();
        }
        TH1D* hProjDrifttime = height_Drift->ProjectionY();


        if (nClusters >= 3) { //nClusters 3 or 4 
            trackYX->Fit("fTrackYX", "QN");
            trackYZ->Fit("fTrackYZ", "QN");
            double parsYX[2], parsYZ[2];
            fTrackYX->GetParameters(parsYX);
            fTrackYZ->GetParameters(parsYZ);
            
            //Right Wing
            double genCoord_drift_NR= -(fDistance+(fPadHeight+fPadGap)*3.5)*parsYZ[1]+parsYZ[0]-fhalfH;
            double dx_genCoord_drift= (-(fDistance+(fPadHeight+fPadGap)*3.5)*parsYZ[1]+parsYZ[0]-fhalfH-fTimeoffsetRW)/(Cos(40*Pi()/180)-parsYZ[1]*Sin(40*Pi()/180));
            double dz_genCoord_drift= (-(fDistance+(fPadHeight+fPadGap)*3.5)*parsYZ[1]+parsYZ[0]-fhalfH-fTimeoffsetRW)/(Sin(40*Pi()/180)+parsYZ[1]*Cos(40*Pi()/180));
            // double dy_genCoord_pad=-(fDistance-(fPadHeight+fPadGap)*3.5)*parsYX[1]+parsYX[0]-(fPadWeith+fPadGap)*16;
            // double dz_genCoord_pad=-(dy_genCoord_pad)/parsYX[1]*Cos(40*Pi()/180);
            double dy_genCoord_pad=-(fDistance-(fPadHeight+fPadGap)*3.5)*parsYX[1]+parsYX[0]-(fPadWeith+fPadGap)*15.5;
            double dz_genCoord_pad=-(dy_genCoord_pad)/parsYX[1]*Cos(40*Pi()/180);
            double dz_diff= dz_genCoord_drift-dz_genCoord_pad;

            hYZ_transform_dx->Fill(dx_genCoord_drift);
            hYZ_transform_dz->Fill(dz_genCoord_drift);
            hYX_transform_dy_pad->Fill(dy_genCoord_pad);
            hYX_transform_dz_pad->Fill(dz_genCoord_pad);
            hTarget_dxdy->Fill(dx_genCoord_drift, dy_genCoord_pad);
            hdeltaZ_diff->Fill(dz_diff);
            hThetaYX_Pad->Fill(ATan(parsYX[1]) * RadToDeg());
            hThetaYL_Drift->Fill(ATan(parsYZ[1]) * RadToDeg());
            if (genCoord_drift_NR>150){
                hConstantYX_Pad->Fill(dy_genCoord_pad);
                hConstantYZ_Drift->Fill(genCoord_drift_NR);
            }
            for (int ii = 0; ii < nClusters; ii++) {
                double dx = xPos[ii] - parsYX[1] * yPos[ii] - parsYX[0];
                double dz = zPos[ii] - parsYZ[1] * yPos[ii] - parsYZ[0];
                hResoX_Pad->Fill(dx);
                hResoL_Drift->Fill(dz);
            }
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // if (dz_genCoord_drift>200){
                // hYZ_transform_dx->Fill(dx_genCoord_drift);
                // hYZ_transform_dz->Fill(dz_genCoord_drift);
                // hYX_transform_dy_pad->Fill(dy_genCoord_pad);
                // hYX_transform_dz_pad->Fill(dz_genCoord_pad);
            cPad->cd(1);
            pad.ADC->SetStats(0);
            pad.ADC->SetMinimum(1);
            pad.ADC->SetMaximum(3500);
            // if (!isPaletteAxisOn) {
            //     auto palette = new TPaletteAxis(66., -5., 71., 89., pad.ADC);
            //     pad.ADC->GetListOfFunctions()->Add(palette);
            //     isPaletteAxisOn = true;
            // }
            pad.ADC->Draw("colz");
            pad.Frame->Draw("same");
            fTrackYX->SetLineWidth(3.);
            fTrackYX->SetLineColor(kRed);
            fTrackYX->SetParameters(-parsYX[0] / parsYX[1], 1. / parsYX[1]);
            fTrackYX->Draw("same");
            trackXY->SetMarkerStyle(20);
            trackXY->SetMarkerSize(1.);
            trackXY->SetMarkerColor(kBlack);
            trackXY->Draw("same p");

            cPad->cd(2);
            // height_Drift->SetStats(0);
            // height_Drift->SetMinimum(0);
            // height_Drift->SetMaximum(500);
            // height_Drift->Draw("LEGO");
            // pad.Frame->Draw("same");
            // pad.DriftTime->SetMinimum(1);
            hProjDrifttime->SetStats(0);
            hProjDrifttime->SetMinimum(0);
            // hProjDrifttime->Fit("pol1", "QN");
            // hProjDrifttime->Draw();
                // cPad->SaveAs(Form("./track_dz200/track_%d.png", idx));
            
            
            // cPad->cd(1);
            // trackYX_TF->Draw();
            // cPad->cd(2);
            // trackYZ_TF->Draw();
            // cPad->SaveAs(Form("./track_delta/del_track_%d.png", idx));
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        }
        trackYX->Delete();
        trackYZ->Delete();
        trackXY->Delete();  // For pad drawing
        // trackYX_TF->Delete();
        // trackYZ_TF->Delete();
        height_Drift->Reset();
        pad.Clear();        // For pad drawing
        // double gparsYX[3], gparYZ[3];
        // fConstantYX_Pad->GetParameters(gparsYX);
        // fConstantYZ_Drift->GetParameters(gparYZ);
    }
    hConstantYX_Pad->Fit("fConstantYX_Pad");
    hConstantYZ_Drift ->Fit("fConstantYZ_Drift");

    fileIn->Close();

    c1->cd(1);
    hThetaYX_Pad->Draw();
    c1->cd(2);
    // hResoX_Pad->Draw();
    hConstantYX_Pad->Draw();
    // c1->SaveAs("./thetaYX_resolution.png");

    c2->cd(1);
    hThetaYL_Drift->Draw();
    c2->cd(2);
    // hResoL_Drift->Draw();
    hConstantYZ_Drift->Draw();
    // c2->SaveAs("./thetaYL_Drift_resolution.png");

    c3->cd(1);
    // hTarget_dxdy->Draw("Colz");
    // hdeltaZ_diff->Draw();
    hYZ_transform_dx->Draw();
    c3->cd(2);
    hYZ_transform_dz->Draw();
    c3->cd(3);
    hYX_transform_dy_pad->Draw();
    c3->cd(4);
    hYX_transform_dz_pad->Draw();
    c3->SaveAs("./hYZ_Rtransform.png");
    // c3->SaveAs("./hdeltaZdiff.png");

    c4->cd(1);
    hConstantYX_Pad->Draw();
    fConstantYX_Pad->Draw("same");

    c4->cd(2);
    hConstantYZ_Drift->Draw();
    fConstantYZ_Drift->Draw("same");

    // c4->SaveAs("./hXYZ_constant.png");
}

void resolutionTestL() {
    GETPad pad;

    TCanvas* cPad;
    gStyle->SetPadTopMargin(0.02);
    gStyle->SetPadBottomMargin(0.085);
    gStyle->SetPadLeftMargin(0.095);
    gStyle->SetPadRightMargin(0.02);
    cPad = new TCanvas("cPad", "", 1500, 750);
    cPad->Divide(2,1);
    bool isPaletteAxisOn = false;

    auto c1 = new TCanvas("c1", "", 2100, 750);
    auto c2 = new TCanvas("c2", "", 2100, 750);
    c1->Divide(3,1);
    c2->Divide(3,1);
    auto c3 = new TCanvas("c3", "", 1500, 1500);
    c3->Divide(2,2);

    auto fGaus = new TF1("fGaus", "gaus", 0, 31);               // The cluster of x_pad direction
    auto fTrackYX = new TF1("fTrackYX", "pol1", -6.25, 93.75);  // x_pad vs y_pad
    auto fTrackYZ = new TF1("fTrackYZ", "pol1", -6.25, 93.75);  // L_drift(z) vs y_pad

    auto hYZ_transform_dx = new TH1D("hYZ_transform_dx-DT", ";#it{#Delta X}_{transformed} [mm]", 200, -100, 400);
    auto hYZ_transform_dz = new TH1D("hYZ_transform_dz-DT", ";#it{#Delta Z}_{transformed} [mm]", 200, -400, 100);
    auto hYX_transform_dy_pad = new TH1D("hYX_transform_dy-pad", ";#it{#Delta Y}_{transformed} [mm]", 200, -100, 100);
    auto hYX_transform_dz_pad = new TH1D("hYX_transform_dz-pad", ";#it{#Delta Z}_{transformed} [mm]", 500, -1000, 1000);
    auto hTarget_dxdy = new TH2D("hTarget_dxdy", "target_(Z=0 plane); #it{#Delta X} [mm]; #it{#Delta Y} [mm]", 400, -100, 300, 200, -100, 100);

    auto hResoX_Pad = new TH1D("hResoX_Pad", ";#it{#Deltax}_{Pad} [mm];", 1000, -1, 1);   // x_pad resoution
    auto hResoL_Drift = new TH1D("hResoL_Drift", ";#it{#DeltaL}_{Drift} [mm];", 100, -4, 4);  // L_drift(z) resoution
    auto hThetaYX_Pad = new TH1D("hThetaYX_Pad", ";#it{#theta}_{YX} [#circ];", 100, -5, 5);  // 
    auto hThetaYL_Drift = new TH1D("hThetaYL_Drift", ";#it{#theta}_{YL} [#circ];", 50, -10, 10);  // 
    auto hConstantYX_Pad = new TH1D("hConstantYX_Pad", ";Constant_{Pad} [mm];", 80, -40, 40);
    auto hConstantYL_Drift = new TH1D("hConstantYL_Drift", ";Constant_{Drift} [mm];", 600, -300, 300);
    auto hDistanceYX_Pad = new TH1D("hDistanceYX_Pad", ";#frac{Constant}{theta} Distance_{Pad} [mm];", 500, -1000, 1000);
    auto hDistanceYL_Drift = new TH1D("hDistanceYL_Drift", ";#frac{Constant}{theta} Distance_{Drift} [mm];", 500, -1000, 1000);
    auto fConstantYL_Drift = new TF1("fConstantYL_Drift", "gaus", -300, -100);
    auto height_Drift= new TH2D("height_Drift", "height_drift", 32, 0, 32, 8, 0, 8);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Event Info.
    int eventID;
    bool isBroken;
    // Left TPC
    int nHits_TL;
    int xId_TL[256];
    int yId_TL[256];
    double driftTime_TL[256];
    double ADC_TL[256];
    // int driftTime_TR[256];
    // int ADC_TR[256];
    // Read file
    TFile* fileIn = new TFile("./treeOfHits.root", "read");
    // TFile* fileIn = new TFile("./treeOfHitsOrg.root", "read");
    TTree* treeIn = (TTree*)fileIn->Get("treeOfHits");
    // Event Info.
    treeIn->SetBranchAddress("EventID", &eventID);
    treeIn->SetBranchAddress("IsBroken", &isBroken);
    // TPC Left
    treeIn->SetBranchAddress("nHits_TL", &nHits_TL);
    treeIn->SetBranchAddress("xId_TL", xId_TL);
    treeIn->SetBranchAddress("yId_TL", yId_TL);
    treeIn->SetBranchAddress("driftTime_TL", driftTime_TL);
    treeIn->SetBranchAddress("ADC_TL", ADC_TL);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (int idx = 0;idx < treeIn->GetEntries(); idx++) {
    // for (int idx = 0; idx < treeIn->GetEntries(); idx++) {
        treeIn->GetEntry(idx);
        if (isBroken) continue;
        if (idx % 1000 == 0) std::cout << idx << std::endl;
        // TPC Left
        int adcPad_TL[32][8] = {0};
        int timePad_TL[32][8] = {0};
        bool isSatdPad_TL[32][8] = {0};
        for (int hitIdx = 0; hitIdx < nHits_TL; hitIdx++) {
            int xId = xId_TL[hitIdx];
            int yId = yId_TL[hitIdx];
            adcPad_TL[xId][yId] = ADC_TL[hitIdx];
            timePad_TL[xId][yId] = driftTime_TL[hitIdx];
        }
        //
        int nClusters = 0;
        double xPos[8], yPos[8], zPos[8];
        auto trackYX = new TGraph();
        auto trackYZ = new TGraph();
        auto trackXY = new TGraph();  // For pad drawing
        for (int yId = 1; yId < 7; yId++) {
            auto hCluster = new TH1D("hCluster", "", 32, -0.5, 31.5);
            for (int xId = 0; xId < 32; xId++) {
                hCluster->SetBinContent(xId + 1, adcPad_TL[xId][yId]);
                pad.ADC->Fill(2. * xId, 12. * yId, adcPad_TL[xId][yId]);  // For pad drawing
                height_Drift->Fill(xId, yId, timePad_TL[xId][yId]);
                // std::cout<< xId<< " "<<yId<<" "<<timePad_TR[xId][yId]<<std::endl;
            }
            int maxBin = hCluster->GetMaximumBin();
            int maxADC = hCluster->GetMaximum();
            if (maxADC > 3500) {
                hCluster->SetBinContent(maxBin - 1, 0);
                hCluster->SetBinContent(maxBin, 0);
                hCluster->SetBinContent(maxBin + 1, 0);
                maxBin = hCluster->GetMaximumBin();
                maxADC = hCluster->GetMaximum();
            }
            int max_xId = maxBin - 1;
            // Only if cluster exists
            if (maxADC > 50 && 2 <= max_xId && max_xId <= 30) {
                if (hCluster->GetBinContent(maxBin - 1) > 50 && hCluster->GetBinContent(maxBin + 1) > 50) {
                    fGaus->SetRange(max_xId - 2, max_xId + 2);
                    hCluster->Fit("fGaus", "QNR");
                    double pars[3];
                    fGaus->GetParameters(pars);
                    xPos[nClusters] = 2. * pars[1];
                    yPos[nClusters] = 12. * yId;
                    zPos[nClusters] = 1.1 * timePad_TL[max_xId][yId];
                    trackYX->SetPoint(trackYX->GetN(), yPos[nClusters], xPos[nClusters]);
                    trackYZ->SetPoint(trackYZ->GetN(), yPos[nClusters], zPos[nClusters]);
                    trackXY->SetPoint(trackXY->GetN(), xPos[nClusters], yPos[nClusters]);  // For pad drawing
                    nClusters++;
                }
            }
            hCluster->Delete();
        }
        if (nClusters >= 3) { //nClusters 3 or 4 
            trackYX->Fit("fTrackYX", "QN");
            trackYZ->Fit("fTrackYZ", "QN");
            double parsYX[2], parsYZ[2];
            fTrackYX->GetParameters(parsYX);
            fTrackYZ->GetParameters(parsYZ);
            //Left wing
            double genCoord_drift_NR= -(fDistance+(fPadHeight+fPadGap)*3.5)*parsYZ[1]-parsYZ[0]+fhalfH;
            double dx_genCoord_drift= (-(fDistance+(fPadHeight+fPadGap)*3.5)*parsYZ[1]-parsYZ[0]+fhalfH-fTimeoffsetLW)/(Cos(40*Pi()/180)+parsYZ[1]*Sin(40*Pi()/180));
            double dz_genCoord_drift= (-(fDistance+(fPadHeight+fPadGap)*3.5)*parsYZ[1]-parsYZ[0]+fhalfH-fTimeoffsetLW)/(-Sin(40*Pi()/180)+parsYZ[1]*Cos(40*Pi()/180));
            // double dy_genCoord_pad=-(fDistance-(fPadHeight+fPadGap)*3.5)*parsYX[1]+parsYX[0]-(fPadWeith+fPadGap)*16;
            // double dz_genCoord_pad=-(dy_genCoord_pad)/parsYX[1]*Cos(40*Pi()/180);
            double dy_genCoord_pad=(fDistance+(fPadHeight+fPadGap)*3.5)*parsYX[1]+parsYX[0]-(fPadWeith+fPadGap)*15.5;
            double dz_genCoord_pad=-(dy_genCoord_pad)/parsYX[1]*Cos(40*Pi()/180);
            double dz_diff= dz_genCoord_drift-dz_genCoord_pad;
            

            hYZ_transform_dx->Fill(dx_genCoord_drift);
            hYZ_transform_dz->Fill(dz_genCoord_drift);
            hYX_transform_dy_pad->Fill(dy_genCoord_pad);
            hYX_transform_dz_pad->Fill(dz_genCoord_pad);
            hTarget_dxdy->Fill(dx_genCoord_drift, dy_genCoord_pad);

            hThetaYX_Pad->Fill(ATan(parsYX[1]) * RadToDeg());
            hThetaYL_Drift->Fill(ATan(parsYZ[1]) * RadToDeg());
            hConstantYX_Pad->Fill(parsYX[0]-(fPadWeith+fPadGap)*16);
            if(genCoord_drift_NR<-150){
                hConstantYL_Drift->Fill(genCoord_drift_NR);
            }
            hDistanceYX_Pad->Fill((parsYX[0]-(fPadWeith+fPadGap)*16)/parsYX[1]);
            hDistanceYL_Drift->Fill((parsYZ[0]-fhalfH)/parsYZ[1]);

            for (int ii = 0; ii < nClusters; ii++) {
                double dx = xPos[ii] - parsYX[1] * yPos[ii] - parsYX[0];
                double dz = zPos[ii] - parsYZ[1] * yPos[ii] - parsYZ[0];
                hResoX_Pad->Fill(dx);
                hResoL_Drift->Fill(dz);
            }
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            cPad->cd(1);
            pad.ADC->SetStats(0);
            pad.ADC->SetMinimum(1);
            pad.ADC->SetMaximum(3500);
            // if (!isPaletteAxisOn) {
            //     auto palette = new TPaletteAxis(66., -5., 71., 89., pad.ADC);
            //     pad.ADC->GetListOfFunctions()->Add(palette);
            //     isPaletteAxisOn = true;
            // }
            pad.ADC->Draw("colz");
            pad.Frame->Draw("same");
            fTrackYX->SetLineWidth(3.);
            fTrackYX->SetLineColor(kRed);
            fTrackYX->SetParameters(-parsYX[0] / parsYX[1], 1. / parsYX[1]);
            fTrackYX->Draw("same");
            trackXY->SetMarkerStyle(20);
            trackXY->SetMarkerSize(1.);
            trackXY->SetMarkerColor(kBlack);
            trackXY->Draw("same p");

            cPad->cd(2);
            height_Drift->SetStats(0);
            height_Drift->SetMinimum(0);
            height_Drift->SetMaximum(500);
            height_Drift->Draw("LEGO");
            // pad.Frame->Draw("same");
            // pad.DriftTime->SetMinimum(1);

            // cPad->SaveAs(Form("./track/track_%d.png", idx));
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        }
        trackYX->Delete();
        trackYZ->Delete();
        trackXY->Delete();  // For pad drawing
        height_Drift->Reset();
        pad.Clear();        // For pad drawing
        
    }
    hConstantYL_Drift ->Fit("fConstantYL_Drift");

    fileIn->Close();
    c1->cd(1);
    hThetaYX_Pad->Draw();
    c1->cd(2);
    // hResoX_Pad->Draw();
    hConstantYX_Pad->Draw();
    c1->cd(3);
    hDistanceYX_Pad->Draw();
    // c1->SaveAs("./thetaYX_resolution.png");

    c2->cd(1);
    hThetaYL_Drift->Draw();
    c2->cd(2);
    // hResoL_Drift->Draw();
    hConstantYL_Drift->Draw();
    fConstantYL_Drift->Draw("same");
    c2->cd(3);
    hDistanceYL_Drift->Draw();
    c2->SaveAs("./thetaYL_Drift_resolution.png");

    c3->cd(1);
    // hTarget_dxdy->Draw("Colz");
    // hdeltaZ_diff->Draw();
    hYZ_transform_dx->Draw();
    c3->cd(2);
    hYZ_transform_dz->Draw();
    c3->cd(3);
    hYX_transform_dy_pad->Draw();
    c3->cd(4);
    hYX_transform_dz_pad->Draw();
    c3->SaveAs("./hYZ_Ltransform.png");
    // c3->SaveAs("./hdeltaZdiff.png");
}

void resolutionTestLandR(){
    GETPad padR;
    GETPad padL;

    TCanvas* cPad;
    gStyle->SetPadTopMargin(0.02);
    gStyle->SetPadBottomMargin(0.085);
    gStyle->SetPadLeftMargin(0.095);
    gStyle->SetPadRightMargin(0.02);
    // gStyle->SetStatX(0.3);
    cPad = new TCanvas("cPad", "", 1500, 750);
    cPad->Divide(2,1);
    bool isPaletteAxisOn = false;

    auto c1 = new TCanvas("c1", "", 750, 750);
    auto c2 = new TCanvas("c2", "", 1500, 1500);
    auto c3 = new TCanvas("c3", "", 1500, 1500);
    c2->Divide(2,2);
    c3->Divide(2,2);

    auto fGausR = new TF1("fGausR", "gaus", 0, 31);               // The cluster of x_pad direction
    auto fGausL = new TF1("fGausL", "gaus", 0, 31);               // The cluster of x_pad direction

    auto fTrackYXR = new TF1("fTrackYXR", "pol1", -6.25, 93.75);  // x_pad vs y_pad
    auto fTrackYZR = new TF1("fTrackYZR", "pol1", -6.25, 93.75);  // L_drift(z) vs y_pad
    auto fTrackYXL = new TF1("fTrackYXL", "pol1", -6.25, 93.75);  // x_pad vs y_pad
    auto fTrackYZL = new TF1("fTrackYZL", "pol1", -6.25, 93.75);  // L_drift(z) vs y_pad

    auto hYZ_Rtransform_dx = new TH1D("hYZ_Rtransform_dx-DT", "Rightwing #it{#Delta X};#it{#Delta X}_{transformed} [mm]", 100, -400, 100);
    auto hYZ_Rtransform_dz = new TH1D("hYZ_Rtransform_dz-DT", "Rightwing #it{#Delta Z};#it{#Delta Z}_{transformed} [mm]", 30, -75, 75);
    auto hYX_Rtransform_dy_pad = new TH1D("hYX_Rtransform_dy-pad", "Rightwing #it{#Delta Y};#it{#Delta Y}_{transformed} [mm]", 100, -100, 100);
    auto hYX_Rtransform_dz_pad = new TH1D("hYX_Rtransform_dz-pad", "Rightwing #it{#Delta Z};#it{#Delta Z}_{transformed} [mm]", 200, -1000, 1000);
    
    auto hYZ_Ltransform_dx = new TH1D("hYZ_Ltransform_dx-DT", "Leftwing #it{#Delta X};#it{#Delta X}_{transformed} [mm]", 100, -100, 400);
    auto hYZ_Ltransform_dz = new TH1D("hYZ_Ltransform_dz-DT", "Leftwing #it{#Delta Z};#it{#Delta Z}_{transformed} [mm]", 30, -75, 75);
    auto hYX_Ltransform_dy_pad = new TH1D("hYX_Ltransform_dy-pad", "Leftwing #it{#Delta Y};#it{#Delta Y}_{transformed} [mm]", 100, -100, 100);
    auto hYX_Ltransform_dz_pad = new TH1D("hYX_Ltransform_dz-pad", "Leftwing #it{#Delta Z};#it{#Delta Z}_{transformed} [mm]", 200, -1000, 1000);
    auto hZ_LandR = new TH2D("hZ_LandR", "dZ_position; #it{#Delta Z}_{R} [mm]; #it{#Delta Y}_{L} [mm]", 30, -75, 75, 30, 75, 75);
    auto fYZ_Rfit= new TF1("fYZ_Rfit", "gaus", -50, 50);
    auto fYZ_Lfit= new TF1("fYZ_Lfit", "gaus", -50, 50);
    // auto hDiff_RandL_dz = new TH1D("hDiff_RandL_dz", ";#it{#Delta Z}_{diff} [mm];", 100, -100, 100);
    // auto hDiff_RandL_dy = new TH1D("hDiff_RandL_dy", ";#it{#Delta Y}_{diff} [mm];", 100, -100, 100);
    // auto hDiff_RandL_dx = new TH1D("hDiff_RandL_dx", ";#it{#Delta X}_{diff} [mm];", 100, -100, 100);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Event Info.
    int eventID;
    bool isBroken;
    // Left TPC
    int nHits_TL;
    int xId_TL[256];
    int yId_TL[256];
    double driftTime_TL[256];
    double ADC_TL[256];
    // Right TPC
    int nHits_TR;
    int xId_TR[256];
    int yId_TR[256];
    double driftTime_TR[256];
    double ADC_TR[256];
    // Read file
    TFile* fileIn = new TFile("./treeOfHits.root", "read");
    TTree* treeIn = (TTree*)fileIn->Get("treeOfHits");
    // Event Info.
    treeIn->SetBranchAddress("EventID", &eventID);
    treeIn->SetBranchAddress("IsBroken", &isBroken);
    // TPC Left
    treeIn->SetBranchAddress("nHits_TL", &nHits_TL);
    treeIn->SetBranchAddress("xId_TL", xId_TL);
    treeIn->SetBranchAddress("yId_TL", yId_TL);
    treeIn->SetBranchAddress("driftTime_TL", driftTime_TL);
    treeIn->SetBranchAddress("ADC_TL", ADC_TL);
    // TPC Right
    treeIn->SetBranchAddress("nHits_TR", &nHits_TR);
    treeIn->SetBranchAddress("xId_TR", xId_TR);
    treeIn->SetBranchAddress("yId_TR", yId_TR);
    treeIn->SetBranchAddress("driftTime_TR", driftTime_TR);
    treeIn->SetBranchAddress("ADC_TR", ADC_TR);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (int idx = 0;idx < treeIn->GetEntries(); idx++) {
    // for (int idx = 0; idx < treeIn->GetEntries(); idx++) {
        treeIn->GetEntry(idx);
        if (isBroken) continue;
        if (idx % 1000 == 0) std::cout << idx << std::endl;
        // TPC Left
        int adcPad_TL[32][8] = {0};
        int timePad_TL[32][8] = {0};
        bool isSatdPad_TL[32][8] = {0};
        for (int hitIdx = 0; hitIdx < nHits_TL; hitIdx++) {
            int xId = xId_TL[hitIdx];
            int yId = yId_TL[hitIdx];
            adcPad_TL[xId][yId] = ADC_TL[hitIdx];
            timePad_TL[xId][yId] = driftTime_TL[hitIdx];
        }
        // TPC Right
        int adcPad_TR[32][8] = {0};
        int timePad_TR[32][8] = {0};
        bool isSatdPad_TR[32][8] = {0};
        for (int hitIdx = 0; hitIdx < nHits_TR; hitIdx++) {
            int xId = xId_TR[hitIdx];
            int yId = yId_TR[hitIdx];
            adcPad_TR[xId][yId] = ADC_TR[hitIdx];
            timePad_TR[xId][yId] = driftTime_TR[hitIdx];
        }
        int nClustersR = 0;
        double xPosR[8], yPosR[8], zPosR[8];
        int nClustersL = 0;
        double xPosL[8], yPosL[8], zPosL[8];
        auto trackYXR = new TGraph();
        auto trackYZR = new TGraph();
        // auto trackXYR = new TGraph();
        auto trackYXL = new TGraph();
        auto trackYZL = new TGraph();
        // auto trackXYL = new TGraph();
        
        for (int yId = 1; yId < 7; yId++) {
            auto hClusterR = new TH1D("hClusterR", "", 32, -0.5, 31.5);
            auto hClusterL = new TH1D("hClusterL", "", 32, -0.5, 31.5);
            for (int xId = 0; xId < 32; xId++) {
                hClusterR->SetBinContent(xId + 1, adcPad_TR[xId][yId]);
                padR.ADC->Fill(2. * xId, 12. * yId, adcPad_TR[xId][yId]);  // For pad drawing
                hClusterL->SetBinContent(xId + 1, adcPad_TL[xId][yId]);
                padL.ADC->Fill(2. * xId, 12. * yId, adcPad_TL[xId][yId]);  // For pad drawing
                
            }
            int maxBinR = hClusterR->GetMaximumBin();
            int maxADCR = hClusterR->GetMaximum();
            int maxBinL = hClusterL->GetMaximumBin();
            int maxADCL = hClusterL->GetMaximum();
            if (maxADCR > 3500) {
                hClusterR->SetBinContent(maxBinR - 1, 0);
                hClusterR->SetBinContent(maxBinR, 0);
                hClusterR->SetBinContent(maxBinR + 1, 0);
                maxBinR = hClusterR->GetMaximumBin();
                maxADCR = hClusterR->GetMaximum();
            }
            if (maxADCL > 3500) {
                hClusterL->SetBinContent(maxBinL - 1, 0);
                hClusterL->SetBinContent(maxBinL, 0);
                hClusterL->SetBinContent(maxBinL + 1, 0);
                maxBinL = hClusterL->GetMaximumBin();
                maxADCL = hClusterL->GetMaximum();
            }
            int max_xIdR = maxBinR - 1;
            int max_xIdL = maxBinL - 1;



            // Only if cluster exists
            if (maxADCR > 50 && 2 <= max_xIdR && max_xIdR <= 30) {
                if (hClusterR->GetBinContent(maxBinR - 1) > 50 && hClusterR->GetBinContent(maxBinR + 1) > 50) {
                    fGausR->SetRange(max_xIdR - 2, max_xIdR + 2);
                    hClusterR->Fit("fGausR", "QNR");
                    double parsR[3];
                    fGausR->GetParameters(parsR);
                    xPosR[nClustersR] = 2. * parsR[1];
                    yPosR[nClustersR] = 12. * yId;
                    zPosR[nClustersR] = 1.1 * timePad_TR[max_xIdR][yId];
                    trackYXR->SetPoint(trackYXR->GetN(), yPosR[nClustersR], xPosR[nClustersR]);
                    trackYZR->SetPoint(trackYZR->GetN(), yPosR[nClustersR], zPosR[nClustersR]);
                    // trackXYR->SetPoint(trackXYR->GetN(), xPosR[nClustersR], yPosR[nClustersR]);  // For pad drawing
                    nClustersR++;
                }
            }
            hClusterR->Delete();
            if (maxADCL > 50 && 2 <= max_xIdL && max_xIdL <= 30) {
                if (hClusterL->GetBinContent(maxBinL - 1) > 50 && hClusterL->GetBinContent(maxBinL + 1) > 50) {
                    fGausL->SetRange(max_xIdL - 2, max_xIdL + 2);
                    hClusterL->Fit("fGausL", "QNR");
                    double parsL[3];
                    fGausL->GetParameters(parsL);
                    xPosL[nClustersL] = 2. * parsL[1];
                    yPosL[nClustersL] = 12. * yId;
                    zPosL[nClustersL] = 1.1 * timePad_TL[max_xIdL][yId];
                    trackYXL->SetPoint(trackYXL->GetN(), yPosL[nClustersL], xPosL[nClustersL]);
                    trackYZL->SetPoint(trackYZL->GetN(), yPosL[nClustersL], zPosL[nClustersL]);
                    // trackXYL->SetPoint(trackXYL->GetN(), xPosL[nClustersL], yPosL[nClustersL]);  // For pad drawing
                    nClustersL++;
                }
            }
            hClusterL->Delete();
        }
        if (16<=nHits_TL && nHits_TL<=45 && 18<=nHits_TR && nHits_TR <=49){ //3sigma range
            if (nClustersR >= 3 && nClustersL >= 3) { //nClusters 3 or 4
                //Right Wing
                trackYXR->Fit("fTrackYXR", "QN");
                trackYZR->Fit("fTrackYZR", "QN");
                double parsYXR[2], parsYZR[2];
                fTrackYXR->GetParameters(parsYXR);
                fTrackYZR->GetParameters(parsYZR); 
                double dx_genCoord_driftR= (-(fDistance+(fPadHeight+fPadGap)*3.5)*parsYZR[1]+parsYZR[0]-fhalfH-fTimeoffsetRW)/(Cos(40*Pi()/180)-parsYZR[1]*Sin(40*Pi()/180));
                double dz_genCoord_driftR= (-(fDistance+(fPadHeight+fPadGap)*3.5)*parsYZR[1]+parsYZR[0]-fhalfH-fTimeoffsetRW)/(Sin(40*Pi()/180)+parsYZR[1]*Cos(40*Pi()/180));
                double dy_genCoord_padR=(fDistance-(fPadHeight+fPadGap)*3.5)*parsYXR[1]-parsYXR[0]+(fPadWeith+fPadGap)*15.5;
                double dz_genCoord_padR=-(dy_genCoord_padR)/parsYXR[1]*Cos(40*Pi()/180);
                
                //Left wing
                trackYXL->Fit("fTrackYXL", "QN");
                trackYZL->Fit("fTrackYZL", "QN");
                double parsYXL[2], parsYZL[2];
                fTrackYXL->GetParameters(parsYXL);
                fTrackYZL->GetParameters(parsYZL);
                double dx_genCoord_driftL= (-(fDistance+(fPadHeight+fPadGap)*3.5)*parsYZL[1]-parsYZL[0]+fhalfH-fTimeoffsetLW)/(Cos(40*Pi()/180)+parsYZL[1]*Sin(40*Pi()/180));
                double dz_genCoord_driftL= (-(fDistance+(fPadHeight+fPadGap)*3.5)*parsYZL[1]-parsYZL[0]+fhalfH-fTimeoffsetLW)/(-Sin(40*Pi()/180)+parsYZL[1]*Cos(40*Pi()/180));
                double dy_genCoord_padL=(fDistance+(fPadHeight+fPadGap)*3.5)*parsYXL[1]+parsYXL[0]-(fPadWeith+fPadGap)*15.5;
                double dz_genCoord_padL=-(dy_genCoord_padL)/parsYXL[1]*Cos(40*Pi()/180);
                if(-50<=dz_genCoord_driftR && dz_genCoord_driftR<=50 && -50<=dz_genCoord_driftL && dz_genCoord_driftL<= 50){
                    //Right Fill
                    hYZ_Rtransform_dx->Fill(dx_genCoord_driftR);
                    hYZ_Rtransform_dz->Fill(dz_genCoord_driftR);
                    hYX_Rtransform_dy_pad->Fill(dy_genCoord_padR);
                    hYX_Rtransform_dz_pad->Fill(dz_genCoord_padR);
                    //Left Fill
                    hYZ_Ltransform_dx->Fill(dx_genCoord_driftL);
                    hYZ_Ltransform_dz->Fill(dz_genCoord_driftL);
                    hYX_Ltransform_dy_pad->Fill(dy_genCoord_padL);
                    hYX_Ltransform_dz_pad->Fill(dz_genCoord_padL);

                    hZ_LandR->Fill(dz_genCoord_driftR, dz_genCoord_driftL);
                    //RandL
                    double diff_RandL_dz= dz_genCoord_driftR-dz_genCoord_driftL;
                    double diff_RandL_dy= dy_genCoord_padR-dy_genCoord_padL;
                    double diff_RandL_dx= dx_genCoord_driftR-dx_genCoord_driftL;
                    // hDiff_RandL_dz->Fill(diff_RandL_dz);
                    // hDiff_RandL_dy->Fill(diff_RandL_dy);
                    // hDiff_RandL_dx->Fill(diff_RandL_dx);
                }
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // cPad->cd(1);
                // padR.ADC->SetStats(0);
                // padR.ADC->SetMinimum(1);
                // padR.ADC->SetMaximum(3500);
                // padR.ADC->Draw("colz");
                // padR.Frame->Draw("same");
                // fTrackYXR->SetLineWidth(3.);
                // fTrackYXR->SetLineColor(kRed);
                // fTrackYXR->SetParameters(-parsYXR[0] / parsYXR[1], 1. / parsYXR[1]);
                // fTrackYXR->Draw("same");
                // // trackXYR->SetMarkerStyle(20);
                // // trackXYR->SetMarkerSize(1.);
                // // trackXYR->SetMarkerColor(kBlack);
                // // trackXYR->Draw("same p");

                // cPad->cd(2);
                // padL.ADC->SetStats(0);
                // padL.ADC->SetMinimum(1);
                // padL.ADC->SetMaximum(3500);
                // padL.ADC->Draw("colz");
                // padL.Frame->Draw("same");
                // fTrackYXL->SetLineWidth(3.);
                // fTrackYXL->SetLineColor(kRed);
                // fTrackYXL->SetParameters(-parsYXL[0] / parsYXL[1], 1. / parsYXL[1]);
                // fTrackYXL->Draw("same");
                // // trackXYL->SetMarkerStyle(20);
                // // trackXYL->SetMarkerSize(1.);
                // // trackXYL->SetMarkerColor(kBlack);
                // // trackXYL->Draw("same p");

                // cPad->SaveAs(Form("./trackRandL/track_%d.png", idx));
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            }
        }
        trackYXR->Delete();
        trackYZR->Delete();
        // trackXYR->Delete();
        trackYXL->Delete();
        trackYZL->Delete();
        // trackXYL->Delete();
        padR.Clear();
        padL.Clear();
    }
    fileIn->Close();
    hYZ_Rtransform_dz->Fit("fYZ_Rfit");
    hYZ_Ltransform_dz->Fit("fYZ_Lfit");
    double parafitR[3], parafitL[3];
    fYZ_Rfit->GetParameters(parafitR);
    fYZ_Lfit->GetParameters(parafitL);
    TPaveText* pave = new TPaveText(0.4,0.7,0.1,0.9);
    pave->AddText(Form("Xmean %.1f", parafitR[1]));
    pave->AddText(Form("Xsigma %.1f", parafitR[2]));
    pave->AddText(Form("Ymean %.1f", parafitL[1]));
    pave->AddText(Form("Ysigma %.1f", parafitL[2]));
    

    c1->cd();
    hZ_LandR->Draw("colz");
    // pave->Draw();
    // c1->cd(1);
    // hDiff_RandL_dx->Draw();
    // c1->cd(2);
    // hDiff_RandL_dy->Draw();
    // c1->cd(3);
    // hDiff_RandL_dz->Draw();
    c1->SaveAs("./dz_RandL.png");

    c2->cd(1);
    hYZ_Rtransform_dx->Draw();
    c2->cd(2);
    hYZ_Rtransform_dz->Draw();
    c2->cd(3);
    hYX_Rtransform_dy_pad->Draw();
    c2->cd(4);
    hYX_Rtransform_dz_pad->Draw();
    c2->SaveAs("./hYZ_Rtransform_RandL_.png");

    c3->cd(1);
    hYZ_Ltransform_dx->Draw();
    c3->cd(2);
    fYZ_Lfit->Draw();
    hYZ_Ltransform_dz->Draw("same");
    c3->cd(3);
    hYX_Ltransform_dy_pad->Draw();
    c3->cd(4);
    hYX_Ltransform_dz_pad->Draw();
    c3->SaveAs("./hYZ_Ltransform_RandL_.png");
}

void readSigFcn() {
    TFile* fileIn = new TFile("./hSignalTemplate.root", "read");
    auto hSignal1DTemplateTR = (TH1D*)fileIn->Get("hSignal1DTemplateTR");
    auto hSignal1DTemplateTL = (TH1D*)fileIn->Get("hSignal1DTemplateTL");
    for (int buckId = 1; buckId <= 350; buckId++) {
        arrSigFcnTR[buckId - 1] = hSignal1DTemplateTR->GetBinContent(buckId);
        arrSigFcnTL[buckId - 1] = hSignal1DTemplateTL->GetBinContent(buckId);
    }
    fileIn->Close();
}

void fitSignalTemplate() {
    GETDecoder decoder;
    auto c1 = new TCanvas("c1", "", 1000, 1000);
    // auto c2 = new TCanvas("c2", "", 1000, 1000);

    auto fFitFcn = new TF1("fFitFcn", fitFcn, 1, 500, 4);
    auto fsigFcn = new TF1("fsigFcn", sigFcnTR, 1, 350, 3);
    auto hSignal1DFactory = new TH1D("hSignal1DFactory", "", 500, 0.5, 500.5);

    TFile* fileIn = new TFile("/home/shlee/workspace/treeOfADC_all_2.root", "read");
    TTree* treeIn = (TTree*)fileIn->Get("treeOfADC");
    int eventID;
    bool isBroken;
    int ADCs[4][4][68][512];  // [nAsAds][nAGETs][nChannels][nTimebuckets]
    treeIn->SetBranchAddress("IsBroken", &isBroken);
    treeIn->SetBranchAddress("ADC", ADCs);
    //
    for (int idx = 0; idx < treeIn->GetEntries(); idx++) {
        treeIn->GetEntry(idx);
        if (idx % 1000 == 0) std::cout << idx << std::endl;
        if (isBroken) continue;
        // TPC
        for (int asadId = 0; asadId <= 0; asadId += 2) {
            for (int agetId = 0; agetId < 4; agetId++) {
                for (int chanId = 0; chanId < 68; chanId++) {
                    if (decoder.IsFPN(chanId) || decoder.IsDead(chanId)) continue;
                    int FPNChanId = decoder.FPNChanId(chanId);
                    auto hSignal = (TH1D*)hSignal1DFactory->Clone("hSignal");
                    auto hPedestal = new TH1D("hPedestal", "", 1000, -499.5, 500.5);
                    for (int buckId = 1; buckId <= 400; buckId++) {
                        int ADC = ADCs[asadId][agetId][chanId][buckId] - ADCs[asadId][agetId][FPNChanId][buckId];
                        if (buckId <= 150) {
                            hPedestal->Fill(ADC);
                        } else {
                            hSignal->SetBinContent(buckId, ADC);
                        }
                    }
                    int pedADC = hPedestal->GetBinCenter(hPedestal->GetMaximumBin());
                    int maxBin = hSignal->GetMaximumBin();
                    int maxADC = hSignal->GetMaximum() - pedADC;
                    hSignal->Delete();
                    hPedestal->Delete();
                    hSignal = (TH1D*)hSignal1DFactory->Clone("hSignal");
                    if (150 < maxBin && maxBin <= 400 && 50 < maxADC && maxADC < 3500) {
                        for (int buckId = 1; buckId <= 500; buckId++) {
                            int ADC = ADCs[asadId][agetId][chanId][buckId] - ADCs[asadId][agetId][FPNChanId][buckId];
                            hSignal->SetBinContent(buckId, ADC - pedADC);
                        }
                        c1->cd();
                        hSignal->Draw();
                        fsigFcn->SetRange(maxBin - 50, maxBin + 50);
                        if (asadId == 0) {
                            // if (maxADC > 3400) {
                            //     fFitFcn->SetParameters(0, maxADC / par1ToADC, maxBin - par3 * par4, maxADC + 50);
                            //     fFitFcn->SetParLimits(1, (maxADC - 5) / par1ToADC, (maxADC + 5000) / par1ToADC);
                            //     fFitFcn->SetParLimits(2, maxBin - 100, maxBin);
                            // } else {
                            //     fFitFcn->SetParameters(0, maxADC / par1ToADC, maxBin - par3 * par4, maxADC + 50);
                            //     fFitFcn->SetParLimits(1, (maxADC - 5) / par1ToADC, (maxADC + 5) / par1ToADC);
                            //     fFitFcn->SetParLimits(2, maxBin - par3 * par4 - 10, maxBin - par3 * par4 + 10);
                            // }
                            // fFitFcn->SetParLimits(0, -10, 10);
                            // fFitFcn->SetParLimits(3, maxADC - 50, maxADC + 50);
                            // fFitFcn->SetLineColor(kRed);
                            // fFitFcn->SetLineWidth(2);
                            // fFitFcn->SetRange(maxBin - 150, maxBin + 100);
                            // hSignal->Fit(fFitFcn, "QNR");
                            // c1->cd();
                            // hSignal->SetAxisRange(-100, 4000, "Y");
                            // hSignal->SetLineColor(kBlue);
                            // hSignal->Draw("hist");
                            // double pars[4];
                            // fFitFcn->GetParameters(pars);
                            // hSignal->SetTitle(Form("%.1lf + %.3lf x (x[0] - %.1lf)^2.685 x exp(-(x[0] - %.1lf) / %.2lf);TimeBucket (10 ns);ADC", pars[0], pars[1], pars[2], pars[2], par4));
                            // fFitFcn->SetRange(1, 500);
                            // fFitFcn->SetNpx(500);
                            // fFitFcn->Draw("same");
                            // c1->Update();
                            // c1->SaveAs(Form("./signal/event%05d_%d_%d.png", idx, agetId, chanId));
                            fsigFcn->SetParameters(maxBin, maxADC / 1000., 0.5);
                            fsigFcn->SetParLimits(1, maxADC * 0.9 / 1000., maxADC * 1.1 / 1000.);
                            fsigFcn->SetParLimits(2, 0, 1.);
                        } else {
                            fsigFcn->SetParameters(maxBin, maxADC / 1000., -0.5);
                            fsigFcn->SetParLimits(1, maxADC * 0.9 / 1000., maxADC * 1.1 / 1000.);
                            fsigFcn->SetParLimits(2, -1., 0);
                        }
                        hSignal->Fit(fsigFcn, "INQR");
                        double pars[2];
                        fsigFcn->GetParameters(pars);
                        hSignal->SetXTitle("Timebuckets [20 ns]");
                        hSignal->SetYTitle("ADC");
                        // hSignal->SetTitle(Form("%d %.1lf %d %.1lf", maxBin, pars[0], maxADC, pars[1] * 1000.));
                        fsigFcn->SetRange(1, 500);
                        fsigFcn->Draw("same");
                        c1->SaveAs(Form("./signal/signal_%d_%d_%d_%d.png", idx, asadId, agetId, chanId));
                    }
                    hSignal->Delete();
                }
            }
        }
    }
    fileIn->Close();
}

void createSignalTemplate() {
    GETDecoder decoder;
    auto c1 = new TCanvas("c1", "", 1000, 1000);
    auto c2 = new TCanvas("c2", "", 1000, 1000);
    auto hSignal1DFactory = new TH1D("hSignal1DFactory", "", 500, 0.5, 500.5);
    auto hSignal2DFactory = new TH2D("hSignal2DFactory", "", 500, 0.5, 500.5, 4096, -500.5, 3595.5);
    auto hSignal1DTemplateTR = new TH1D("hSignal1DTemplateTR", "", 350, 0.5, 350.5);
    auto hSignal1DTemplateTL = new TH1D("hSignal1DTemplateTL", "", 350, 0.5, 350.5);
    auto hSignal2DTemplateTR = new TH2D("hSignal2DTemplateTR", "", 350, 0.5, 350.5, 1500, -49.5, 1450.5);
    auto hSignal2DTemplateTL = new TH2D("hSignal2DTemplateTL", "", 350, 0.5, 350.5, 1500, -49.5, 1450.5);

    TFile* fileIn = new TFile("/home/shlee/workspace/treeOfADC_all_2.root", "read");
    TTree* treeIn = (TTree*)fileIn->Get("treeOfADC");
    int eventID;
    bool isBroken;
    int ADCs[4][4][68][512];  // [nAsAds][nAGETs][nChannels][nTimebuckets]
    treeIn->SetBranchAddress("IsBroken", &isBroken);
    treeIn->SetBranchAddress("ADC", ADCs);
    //
    TFile* fileOut = new TFile("./hSignalTemplate.root", "recreate");
    int signalTemplateTR[250];
    int signalTemplateTL[250];
    //
    // for (int idx = 0; idx < 1000; idx++) {
    for (int idx = 0; idx < treeIn->GetEntries(); idx++) {
        treeIn->GetEntry(idx);
        if (idx % 1000 == 0) std::cout << idx << std::endl;
        if (isBroken) continue;
        // TPC
        for (int asadId = 0; asadId <= 2; asadId += 2) {
            for (int agetId = 0; agetId < 4; agetId++) {
                for (int chanId = 0; chanId < 68; chanId++) {
                    if (decoder.IsFPN(chanId) || decoder.IsDead(chanId)) continue;
                    int FPNChanId = decoder.FPNChanId(chanId);
                    auto hSignal = (TH1D*)hSignal1DFactory->Clone("hSignal");
                    auto hPedestal = new TH1D("hPedestal", "", 100, -49.5, 50.5);
                    for (int buckId = 1; buckId <= 400; buckId++) {
                        int ADC = ADCs[asadId][agetId][chanId][buckId] - ADCs[asadId][agetId][FPNChanId][buckId];
                        if (buckId <= 150) {
                            if (-50 < ADC && ADC <= 50) {
                                hPedestal->Fill(ADC);
                            }
                        } else {
                            hSignal->SetBinContent(buckId, ADC);
                        }
                    }
                    int pedADC = hPedestal->GetBinCenter(hPedestal->GetMaximumBin());
                    int maxBin = hSignal->GetMaximumBin();
                    int maxADC = hSignal->GetMaximum() - pedADC;
                    hSignal->Delete();
                    hPedestal->Delete();
                    if (150 < maxBin && maxBin <= 400 && 500 < maxADC && maxADC < 3500) {
                        int ADC0 = ADCs[asadId][agetId][chanId][maxBin - 150] - ADCs[asadId][agetId][FPNChanId][maxBin - 150] - pedADC;
                        for (int buckId = maxBin - 150; buckId <= maxBin + 200; buckId++) {
                            int ADC = ADCs[asadId][agetId][chanId][buckId] - ADCs[asadId][agetId][FPNChanId][buckId] - pedADC;
                            if (asadId == 0) hSignal2DTemplateTR->Fill(buckId - maxBin + 150, (1000. / maxADC) * (ADC - ADC0));
                            if (asadId == 2) hSignal2DTemplateTL->Fill(buckId - maxBin + 150, (1000. / maxADC) * (ADC - ADC0));
                        }
                    }
                }
            }
        }
    }
    fileIn->Close();
    for (int buckId = 1; buckId <= 350; buckId++) {
        auto hTempTR = hSignal2DTemplateTR->ProjectionY("hTempTR", buckId, buckId);
        auto hTempTL = hSignal2DTemplateTL->ProjectionY("hTempTL", buckId, buckId);
        int modeADC_TR = hTempTR->GetBinCenter(hTempTR->GetMaximumBin());
        int modeADC_TL = hTempTL->GetBinCenter(hTempTL->GetMaximumBin());
        hSignal1DTemplateTR->SetBinContent(buckId, modeADC_TR);
        hSignal1DTemplateTL->SetBinContent(buckId, modeADC_TL);
        hTempTR->Delete();
        hTempTL->Delete();
    }
    c1->cd();
    hSignal2DTemplateTR->Draw("colz");
    hSignal1DTemplateTR->SetLineColor(kRed);
    hSignal1DTemplateTR->Draw("same");
    c2->cd();
    hSignal2DTemplateTL->Draw("colz");
    hSignal1DTemplateTL->SetLineColor(kRed);
    hSignal1DTemplateTL->Draw("same");
    //
    hSignal1DTemplateTR->Write();
    hSignal2DTemplateTR->Write();
    hSignal1DTemplateTL->Write();
    hSignal2DTemplateTL->Write();
    fileOut->Close();
}

std::string getLabel(int asadId, int agetId, int chanId) {
    std::stringstream label;
    std::string detector;
    std::string rightOrLeft = (asadId == 0 || asadId == 1) ? "r" : "l";
    std::string downOrUp;
    std::string ohmicOrJunction;
    if (chanId == 11 || chanId == 22 || chanId == 45 || chanId == 56) {
        label << "FPN";
        return label.str();
    }
    if (asadId == 0 || asadId == 2) {
        detector = "t";
    } else if (agetId == 2) {
        label << "none";
        return label.str();
    } else if (asadId == 1 && agetId == 3 && chanId == 49) {
        detector = "c";
        downOrUp = "d";
    } else if (asadId == 1 && agetId == 3 && chanId == 15) {
        detector = "c";
        downOrUp = "u";
    } else if (asadId == 3 && agetId == 3 && chanId == 66) {
        detector = "c";
        downOrUp = "d";
    } else if (asadId == 3 && agetId == 3 && chanId == 32) {
        detector = "c";
        downOrUp = "u";
    } else if (agetId == 0 && chanId <= 12) {
        detector = "s";
        ohmicOrJunction = "o";
        downOrUp = (chanId <= 3)    ? "d"
                   : (chanId <= 7)  ? "u"
                   : (chanId <= 12) ? "b"
                                    : "";
    } else if (agetId == 1 && chanId <= 50) {
        detector = "s";
        ohmicOrJunction = "j";
        downOrUp = (chanId <= 16)   ? "d"
                   : (chanId <= 33) ? "u"
                   : (chanId <= 50) ? "b"
                                    : "";
    } else {
        label << "none";
        return label.str();
    }

    label << detector << rightOrLeft << downOrUp << ohmicOrJunction;
    return label.str();
}
