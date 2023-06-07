#include "GETPad.h"
//
GETPad::GETPad() {
    // int xId = 0;
    // int yId = 0;
    // for (int agetId = 0; agetId < 4; agetId++) {
    //     switch (agetId) {
    //         case 0:
    //             xId = 16;
    //             yId = 3;
    //             break;
    //         case 1:
    //             xId = 16;
    //             yId = 7;
    //             break;
    //         case 2:
    //             xId = 8;
    //             yId = 7;
    //             break;
    //         case 3:
    //             xId = 0;
    //             yId = 7;
    //             break;
    //     }
    //     for (int chanId = 0; chanId < 68; chanId++) {
    //         if (chanId == 11 || chanId == 22 || chanId == 45 || chanId == 56) {
    //             this->xId_[agetId][chanId] = -1;
    //             this->yId_[agetId][chanId] = -1;
    //             continue;
    //         }
    //         if (agetId == 0 || agetId == 1) {
    //             this->agetId_[xId][yId] = agetId;
    //             this->chanId_[xId][yId] = chanId;
    //             this->xId_[agetId][chanId] = xId++;
    //             this->yId_[agetId][chanId] = yId;
    //             if (xId == 32) {
    //                 xId = 16;
    //                 yId--;
    //             }
    //         } else if (agetId == 2 || agetId == 3) {
    //             this->agetId_[xId][yId] = agetId;
    //             this->chanId_[xId][yId] = chanId;
    //             this->xId_[agetId][chanId] = xId;
    //             this->yId_[agetId][chanId] = yId--;
    //             if (yId == -1) {
    //                 xId++;
    //                 yId = 7;
    //             }
    //         }
    //     }
    // }
    int xId = 0;
    int yId = 0;

    for (int agetId = 0; agetId < 4; agetId++) {
        xId = 16;
        yId = 6 - (2 * agetId);

        int nChanId = 0;
        int padSorterIdx = 0;
        for (int chanId = 0; chanId < 68; chanId++) {
            if (chanId == 11 || chanId == 22 || chanId == 45 || chanId == 56) continue;
            nChanId++;
            int idxFPN = (nChanId - 1) / 16;
            if (chanId == 33 || chanId == 34 || chanId == 35) {
                xId = 0;
                if (chanId == 34) {
                    xId = 1;
                    yId++;
                    padSorterIdx = 14;
                }
                this->agetId_[xId][yId] = agetId;
                this->chanId_[xId][yId] = chanId;
                this->xId_[agetId][chanId] = xId;
                this->yId_[agetId][chanId] = yId;
                continue;
            }

            int sign = -1;
            if (nChanId % 2 == 0) {
                sign = 1;
                padSorterIdx++;
            } else if (chanId >= 38) {
                padSorterIdx -= 2;
            }
            xId = 16 + sign * (padSorterIdx);
            this->agetId_[xId][yId] = agetId;
            this->chanId_[xId][yId] = chanId;
            this->xId_[agetId][chanId] = xId;
            this->yId_[agetId][chanId] = yId;
        }
    }
    this->AddPad();
}
void GETPad::AddPad() {
    double xMin = 15.5 * (this->dx_ + this->gap_) - 52;
    double xMax = 52 + 15.5 * (this->dx_ + this->gap_);
    double yMin = 3.5 * (this->dy_ + this->gap_) - 52;
    double yMax = 52 + 3.5 * (this->dy_ + this->gap_);
    int padId = 0;
    while (gROOT->FindObject(Form("cPad_%d", padId)) != nullptr) {
        padId++;
    }
    this->ADC = new TH2Poly(Form("hPad_ADC_%d", padId), ";#it{x}_{Pad} [mm];#it{y}_{Pad} [mm]", xMin, xMax, yMin, yMax);
    this->DriftTime = new TH2Poly(Form("hPad_driftTime_%d", padId), ";[mm];[mm]", xMin, xMax, yMin, yMax);
    this->Frame = new TH2Poly(Form("Frame_%d", padId), "", xMin, xMax, yMin, yMax);
    this->DUMMY = new TH1D(Form("DUMMY_%d", padId), "", 1, 0, 1);
    // Add the bins
    double x[4], y[4];
    // Set Rectangle pattern parameters
    double xCenter = 0., yCenter = 0.;
    for (int xId = 0; xId < 32; xId++) {
        for (int yId = 0; yId < 8; yId++) {
            x[0] = xCenter - 0.5 * this->dx_;
            y[0] = yCenter - 0.5 * this->dy_;
            x[1] = x[0];
            y[1] = yCenter + 0.5 * this->dy_;
            x[2] = xCenter + 0.5 * this->dx_;
            y[2] = y[1];
            x[3] = x[2];
            y[3] = y[0];
            this->ADC->AddBin(4, x, y);
            this->DriftTime->AddBin(4, x, y);
            this->Frame->AddBin(4, x, y);
            // Go up
            yCenter += this->dy_ + this->gap_;
        }
        xCenter += this->dx_ + this->gap_;
        yCenter = 0.;
    }
}
void GETPad::Clear() {
    this->ADC->ClearBinContents();
    this->DriftTime->ClearBinContents();
}
int GETPad::GetAgetId(int xId, int yId) {
    return this->agetId_[xId][yId];
}
int GETPad::GetChanId(int xId, int yId) {
    return this->chanId_[xId][yId];
}
int GETPad::GetXId(int agetId, int chanId) {
    return this->xId_[agetId][chanId];
}
int GETPad::GetYId(int agetId, int chanId) {
    return this->yId_[agetId][chanId];
}
double GETPad::GetX(int agetId, int chanId) {
    return (this->dx_ + this->gap_) * this->GetXId(agetId, chanId);
}
double GETPad::GetY(int agetId, int chanId) {
    return (this->dy_ + this->gap_) * this->GetYId(agetId, chanId);
}
