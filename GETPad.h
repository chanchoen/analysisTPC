#ifndef __GETPad_h
#define __GETPad_h

#include <iostream>
#include <string>

#include "TROOT.h"
#include "TH1D.h"
#include "TH2Poly.h"
#include "TCanvas.h"

class GETPad {
private:
    const double dx_ = 1.7;
    const double dy_ = 11.7;
    const double gap_ = 0.3;
    int xId_[4][68];
    int yId_[4][68];
    int agetId_[32][8];
    int chanId_[32][8];
    TH1D* DUMMY;

public:
    GETPad();

    TH2Poly* Frame;
    TH2Poly* ADC;
    TH2Poly* DriftTime;
    void AddPad();
    void Clear();

    int GetAgetId(int xId, int yId);
    int GetChanId(int xId, int yId);
    int GetXId(int agetId, int chanId);
    int GetYId(int agetId, int chanId);
    double GetX(int agetId, int chanId);
    double GetY(int agetId, int chanId);
};
#endif
