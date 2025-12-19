#include "TROOT.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TASImage.h"
#include "TApplication.h"
#include "TSystem.h"

#include <cassert>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <chrono>
#include <string>
#include <iomanip>

using namespace std;

typedef struct { uint8_t red, green, blue; } COLOR;

static inline int color_distance2(int r1,int g1,int b1,int r2,int g2,int b2){
    int dr=r1-r2, dg=g1-g2, db=b1-b2;
    return 2*dr*dr + 4*dg*dg + 3*db*db;
}

static inline int R(UInt_t p){ return (p>>16)&0xFF; }
static inline int G(UInt_t p){ return (p>>8 )&0xFF; }
static inline int B(UInt_t p){ return (p    )&0xFF; }

static inline int clampi(int x,int lo,int hi){
    return (x<lo)?lo:((x>hi)?hi:x);
}

static inline int irand(int lo, int hi){
    return lo + (rand() % (hi - lo + 1));
}

static inline double urand01(){
    return (double)rand() / (double)RAND_MAX;
}

// my: deltaE only (no swap here)
static inline int deltaE_swap(const vector<COLOR> &sourcepic,
                              const vector<COLOR> &targetpic,
                              int a, int b)
{
    const COLOR &Sa = sourcepic[a];
    const COLOR &Sb = sourcepic[b];
    const COLOR &Ta = targetpic[a];
    const COLOR &Tb = targetpic[b];

    int oldE = color_distance2((int)Sa.red,(int)Sa.green,(int)Sa.blue, (int)Ta.red,(int)Ta.green,(int)Ta.blue)
             + color_distance2((int)Sb.red,(int)Sb.green,(int)Sb.blue, (int)Tb.red,(int)Tb.green,(int)Tb.blue);

    int newE = color_distance2((int)Sb.red,(int)Sb.green,(int)Sb.blue, (int)Ta.red,(int)Ta.green,(int)Ta.blue)
             + color_distance2((int)Sa.red,(int)Sa.green,(int)Sa.blue, (int)Tb.red,(int)Tb.green,(int)Tb.blue);

    return (newE - oldE);
}

static long long total_energy(const vector<COLOR> &sourcepic, const vector<COLOR> &targetpic){
    long long E = 0;
    int N = (int)sourcepic.size();
    for(int i=0;i<N;i++){
        E += color_distance2((int)sourcepic[i].red,(int)sourcepic[i].green,(int)sourcepic[i].blue,
                             (int)targetpic[i].red,(int)targetpic[i].green,(int)targetpic[i].blue);
    }
    return E;
}

static double EstimateT0_image(const vector<COLOR> &sourcepic,
                               const vector<COLOR> &targetpic,
                               int numPix, int W, int H,
                               int nProbe, double factor,
                               double pGlobal, int radius0, int radiusMin)
{
    vector<int> dEs;
    dEs.reserve(nProbe);

    for(int k=0;k<nProbe;k++){
        int a = rand() % numPix;
        int b = a;

        if (urand01() < pGlobal){
            b = rand() % numPix;
        } else {
            int ax = a % W;
            int ay = a / W;

            int rad = radius0;
            if (rad < radiusMin) rad = radiusMin;

            int dx = irand(-rad, rad);
            int dy = irand(-rad, rad);

            int bx = clampi(ax + dx, 0, W-1);
            int by = clampi(ay + dy, 0, H-1);
            b = by*W + bx;
        }

        if(a==b){ k--; continue; }

        int dE = deltaE_swap(sourcepic, targetpic, a, b);
        dEs.push_back(std::abs(dE)); // my: stable T0 estimate
    }

    if(dEs.empty()) return 1.0;

    std::sort(dEs.begin(), dEs.end());
    int idx = (int)(0.90 * (dEs.size() - 1));
    return factor * (double)dEs[idx];
}

int main(int argc, char **argv) {

    if (argc < 3) {
        cout << "Usage: simpix image1 image2 <output=out.png> [--batch]" << endl;
        return 0;
    }

    TString fsrc = argv[1];
    TString ftgt = argv[2];

    bool batchMode=false;
    TString fout = "newout.png";

    for(int i=3;i<argc;i++){
        string a = argv[i];
        if(a=="--batch") batchMode=true;
        else if(!a.empty() && a[0] != '-') fout = argv[i];
    }

    cout << "Reading images: source= " << fsrc << " target= " << ftgt << endl;
    cout << "Output= " << fout << endl;

    TApplication theApp("App", &argc, argv);
    if(batchMode) gROOT->SetBatch(kTRUE);

    TASImage *src = new TASImage(fsrc.Data());
    TASImage *tgt = new TASImage(ftgt.Data());
    TASImage *out = new TASImage(*src);

    assert(src->GetWidth() == tgt->GetWidth() && src->GetHeight() == tgt->GetHeight());
    cout << "Pixel Geometry: " << src->GetWidth() << " x " << src->GetHeight() << endl;

    Long_t numPixL = (Long_t)src->GetWidth() * (Long_t)src->GetHeight();
    int numPix = (int)numPixL;

    UInt_t *srcPix = src->GetArgbArray();
    UInt_t *tgtPix = tgt->GetArgbArray();
    UInt_t *outPix = out->GetArgbArray();

    vector<COLOR> sourcepic(numPix);
    vector<COLOR> targetpic(numPix);
    vector<COLOR> bestpic(numPix);

    for(int i=0;i<numPix;i++){
        sourcepic[i].red   = (uint8_t)R(srcPix[i]);
        sourcepic[i].green = (uint8_t)G(srcPix[i]);
        sourcepic[i].blue  = (uint8_t)B(srcPix[i]);

        targetpic[i].red   = (uint8_t)R(tgtPix[i]);
        targetpic[i].green = (uint8_t)G(tgtPix[i]);
        targetpic[i].blue  = (uint8_t)B(tgtPix[i]);
    }

    srand(12345);

    int W = src->GetWidth();
    int H = src->GetHeight();

    // my: scale steps with image size (640x480 ~ 15M), cap at 60M
    long long runsLL = 50LL * (long long)numPix;
    if (runsLL < 8000000LL)  runsLL = 8000000LL;
    if (runsLL > 60000000LL) runsLL = 60000000LL;
    int runs = (int)runsLL;

    double pGlobalMin = 0.03; // my: rare long swaps late
    double pGlobalMax = 0.20; // my: more mixing early

    int radiusMin = 2;
    int radius0 = std::max(16, std::min(256, std::max(W, H) / 4));

    long long curE  = total_energy(sourcepic, targetpic);

    // my: bestpic should match bestE we print
    long long bestE = curE;
    long long bestCandE = curE;
    std::memcpy(bestpic.data(), sourcepic.data(), (size_t)numPix * sizeof(COLOR));

    // my: avoid copying bestpic on every tiny improvement
    int improveSinceSave = 0;
    const int SAVE_EVERY_IMPROVES = 200;
    const long long SAVE_MIN_DROP = 3000000;

    cout << "start E = " << curE << endl;

    int nProbe = std::min(50000, numPix);
    double T0  = EstimateT0_image(sourcepic, targetpic, numPix, W, H,
                                  nProbe, 1.2, pGlobalMax, radius0, radiusMin);

    // my: not too cold for short runs
    double Tf = 50.0;
    if (Tf >= T0) Tf = 0.01 * T0;

    double alpha = exp(log(Tf / T0) / (double)runs);
    double T = T0;

    cout << "Auto T0 = " << T0 << " (probe=" << nProbe << ")\n";
    cout.setf(std::ios::fixed);
    cout << "Tf=" << Tf << " alpha=" << std::setprecision(9) << alpha << "\n";

    auto t_start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < runs; i++) {

        int a = rand() % numPix;
        int b = a;

        int radiusNow = (int)(radius0 * (T / T0)); // my: big early, small late
        if (radiusNow < radiusMin) radiusNow = radiusMin;

        double frac = T / T0;
        if (frac < 0) frac = 0;
        if (frac > 1) frac = 1;
        double pGlobalNow = pGlobalMin + (pGlobalMax - pGlobalMin) * frac;

        if (urand01() < pGlobalNow){
            b = rand() % numPix;
        } else {
            int ax = a % W;
            int ay = a / W;
            int dx = irand(-radiusNow, radiusNow);
            int dy = irand(-radiusNow, radiusNow);
            int bx = clampi(ax + dx, 0, W-1);
            int by = clampi(ay + dy, 0, H-1);
            b = by*W + bx;
        }

        if (a == b) { T *= alpha; continue; }

        int dE = deltaE_swap(sourcepic, targetpic, a, b);
        long long newE = curE + (long long)dE;

        bool accept = false;
        if (dE <= 0) accept = true;
        else{
            double u = urand01();
            if (u < exp(-(double)dE / T)) accept = true;
        }

        if (accept){
            std::swap(sourcepic[a], sourcepic[b]);
            curE = newE;

            if (curE < bestCandE) {
                bestCandE = curE;
                improveSinceSave++;

                if (improveSinceSave >= SAVE_EVERY_IMPROVES || (bestE - bestCandE) >= SAVE_MIN_DROP) {
                    bestE = bestCandE;
                    std::memcpy(bestpic.data(), sourcepic.data(), (size_t)numPix * sizeof(COLOR));
                    improveSinceSave = 0;
                }
            }
        }

        T *= alpha;

        if (i % 100000 == 0) {
            cout << i << " T=" << T << " curE=" << curE << " bestE=" << bestE << endl;
        }
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dt = t_end - t_start;

    cout << "runtime_seconds = " << dt.count() << endl;
    cout << "runtime_minutes = " << dt.count()/60.0 << endl;

    unsigned a255 = 255;
    for (int i = 0; i < numPix; i++) {
        int rr = (int)bestpic[i].red;
        int gg = (int)bestpic[i].green;
        int bb = (int)bestpic[i].blue;
        outPix[i] = (a255 << 24) | (rr << 16) | (gg << 8) | bb;
    }

    out->WriteImage(fout.Data());

    TCanvas *c1 = new TCanvas("c1", "images", 1024, 768);
    c1->Divide(2,2);

    c1->cd(1); src->Draw("X");
    c1->cd(2); tgt->Draw("X");
    c1->cd(3); out->Draw("X");

    auto base = [](TString s){
      s = gSystem->BaseName(s);
      s.ReplaceAll(".png","");
      s.ReplaceAll(".jpg","");
      s.ReplaceAll(".jpeg","");
      return s;
    };

    TString tag = Form("%s_to_%s", base(fsrc).Data(), base(ftgt).Data());
    c1->Print(Form("collage_%s.png", tag.Data()));

    if (!batchMode) {
      cout << "Press ^c to exit" << endl;
      theApp.SetIdleTimer(30,".q");
      theApp.Run();
    }

    return 0;
}
