#include "TROOT.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TASImage.h"
#include "TApplication.h"
#include "TSystem.h"

#include "assert.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <chrono>


using namespace std;

typedef struct {
    int red, green, blue;
} COLOR;

// fast distance (no sqrt / no pow)
static inline int color_distance2(int r1,int g1,int b1,int r2,int g2,int b2){
    int dr=r1-r2, dg=g1-g2, db=b1-b2;
    // weights help a bit
    return 2*dr*dr + 4*dg*dg + 3*db*db;
}

// pull rgb from UInt_t ARGB
static inline int R(UInt_t p){ return (p>>16)&0xFF; }
static inline int G(UInt_t p){ return (p>>8 )&0xFF; }
static inline int B(UInt_t p){ return (p    )&0xFF; }

static inline int clampi(int x,int lo,int hi){
    return (x<lo)?lo:((x>hi)?hi:x);
}

double exponentialCooling(int iteration, int maxIterations) {
    double T0 = 200000.0;
    double Tend = 50.0;
    double t = (double)iteration / (double)(maxIterations-1);
    return T0 * pow(Tend/T0, t);
}


// propose a swap and compute dE (O(1))
// also does the swap, so caller must undo if rejected
static inline int swap_deltaE(vector<COLOR> &sourcepic,
                             const vector<COLOR> &targetpic,
                             int a, int b)
{
    const COLOR &Sa = sourcepic[a];
    const COLOR &Sb = sourcepic[b];
    const COLOR &Ta = targetpic[a];
    const COLOR &Tb = targetpic[b];

    int oldE = color_distance2(Sa.red,Sa.green,Sa.blue, Ta.red,Ta.green,Ta.blue)
             + color_distance2(Sb.red,Sb.green,Sb.blue, Tb.red,Tb.green,Tb.blue);

    int newE = color_distance2(Sb.red,Sb.green,Sb.blue, Ta.red,Ta.green,Ta.blue)
             + color_distance2(Sa.red,Sa.green,Sa.blue, Tb.red,Tb.green,Tb.blue);

    // do the swap
    std::swap(sourcepic[a], sourcepic[b]);

    return (newE - oldE);
}

// compute full energy once (O(N))
long long total_energy(const vector<COLOR> &sourcepic, const vector<COLOR> &targetpic){
    long long E = 0;
    int N = (int)sourcepic.size();
    for(int i=0;i<N;i++){
        E += color_distance2(sourcepic[i].red,sourcepic[i].green,sourcepic[i].blue,
                             targetpic[i].red,targetpic[i].green,targetpic[i].blue);
    }
    return E;
}

int main(int argc, char **argv) {

    if (argc < 3) {
        cout << "Usage: simpix image1 image2 <output=out.png> [--batch]" << endl;
        return 0;
    }

    TString fsrc = argv[1];
    TString ftgt = argv[2];
    TString fout;
    argc > 3 ? fout = argv[3] : fout = "newout.png";

    bool batchMode=false;
    for(int i=1;i<argc;i++){
        if(string(argv[i])=="--batch") batchMode=true;
    }

    cout << "Reading images: source= " << fsrc << " target= " << ftgt << endl;
    cout << "Output= " << fout << endl;

    TApplication theApp("App", &argc, argv);
    if(batchMode) gROOT->SetBatch(kTRUE);

    // create image objects
    TASImage *src = new TASImage(fsrc.Data());
    TASImage *tgt = new TASImage(ftgt.Data());
    TASImage *out = new TASImage(*src); // start with a copy of the source

    // Test image geometry, exit if they are not the same dimensions
    assert(src->GetWidth() == tgt->GetWidth() && src->GetHeight() == tgt->GetHeight());
    cout << "Pixel Geometry: " << src->GetWidth() << " x " << src->GetHeight() << endl;

    Long_t numPixL = (Long_t)src->GetWidth() * (Long_t)src->GetHeight();
    int numPix = (int)numPixL;

    // *** The work happens here

    UInt_t *srcPix = src->GetArgbArray();
    UInt_t *tgtPix = tgt->GetArgbArray();
    UInt_t *outPix = out->GetArgbArray();

    // use vectors (stack will die for 1080p)
    vector<COLOR> sourcepic(numPix);
    vector<COLOR> targetpic(numPix);
    vector<COLOR> bestpic(numPix);

    // init from images (this was missing)
    for(int i=0;i<numPix;i++){
        sourcepic[i].red   = R(srcPix[i]);
        sourcepic[i].green = G(srcPix[i]);
        sourcepic[i].blue  = B(srcPix[i]);

        targetpic[i].red   = R(tgtPix[i]);
        targetpic[i].green = G(tgtPix[i]);
        targetpic[i].blue  = B(tgtPix[i]);
    }

    // SA state
    srand(12345);

    long long curE  = total_energy(sourcepic, targetpic);
    long long bestE = curE;
    bestpic = sourcepic;

    cout << "start E = " << curE << endl;

    
    int runs = 8000000; // tune
    int W = src->GetWidth();
    int H = src->GetHeight();
    double pLocal = 0.97;
    int radius = 16;
    
    auto t_start = std::chrono::high_resolution_clock::now();


    for (int i = 0; i < runs; i++) {
        double T = exponentialCooling(i, runs);

        // pick two positions (mostly local swaps)
        int a = rand() % numPix;
        int b;

        if (((double)rand()/RAND_MAX) < pLocal){
            int ax = a % W;
            int ay = a / W;
            int dx = (rand() % (2*radius+1)) - radius;
            int dy = (rand() % (2*radius+1)) - radius;
            int bx = clampi(ax + dx, 0, W-1);
            int by = clampi(ay + dy, 0, H-1);
            b = by*W + bx;
        }else{
            b = rand() % numPix;
        }

        if (a == b) continue;

        // do swap + get dE
        int dE = swap_deltaE(sourcepic, targetpic, a, b);
        long long newE = curE + dE;

        // accept / reject
        bool accept = false;
        if (dE <= 0) accept = true;
        else{
            double u = (double)rand() / (double)RAND_MAX;
            if (u < exp(-(double)dE / T)) accept = true;
        }

        if (accept){
            curE = newE;

            if (curE < bestE){
                bestE = curE;
                bestpic = sourcepic;
            }
        } else {
            // undo swap if rejected
            std::swap(sourcepic[a], sourcepic[b]);
        }

        if (i % 100000 == 0) {
            cout << i << " curE=" << curE << " bestE=" << bestE << endl;
        }
    }
    
    auto t_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dt = t_end - t_start;

    cout << "runtime_seconds = " << dt.count() << endl;
    cout << "runtime_minutes = " << dt.count()/60.0 << endl;


    // write best image to outPix
    unsigned a255 = 255;
    for (int i = 0; i < numPix; i++) {
        int rr = clampi(bestpic[i].red,   0, 255);
        int gg = clampi(bestpic[i].green, 0, 255);
        int bb = clampi(bestpic[i].blue,  0, 255);
        outPix[i] = (a255 << 24) | (rr << 16) | (gg << 8) | bb;
    }

    // *************************

    if(!batchMode){
        TCanvas *c1 = new TCanvas("c1", "images", 1024, 768);
        c1->Divide(2, 2);

        c1->cd(1); c1->Draw(); src->Draw("X");
        c1->cd(2); tgt->Draw("X");
        c1->cd(3); out->Draw("X");
        c1->Print("newcollage.png");
    }

    out->WriteImage(fout.Data());

    if(batchMode) return 0;

    cout << "Press ^c to exit" << endl;
    theApp.SetIdleTimer(30, ".q");
    theApp.Run();
    return 0;
}

