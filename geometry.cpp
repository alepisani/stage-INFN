#include <TApplication.h>
#include <TCanvas.h>
#include <TView.h>
#include <TList.h>
#include <TPolyLine3D.h>
#include "TH1F.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TGraph.h"
#include "TFile.h"
#include "TPad.h"
#include "TPolyMarker3D.h"
#include "TMarker3DBox.h"
#include <iostream>
using namespace std;

int events = 10000;    //max 100000 otherwise my laptop explode!!
double pi = TMath::Pi();

//trigger 1 (bottom trigger)
static constexpr double TR1CenterZ = 0.;
static constexpr double TR1Thickness = 2;
static constexpr std::array<double, 3> TR1Size = {154.6,32.5,TR1Thickness};
static constexpr double TR1GapY = 1.9;
//trigger 2 (top trigger)
static constexpr double TR2CenterZ = TR1CenterZ+60.5;
static constexpr double TR2Thickness = 8;
static constexpr std::array<double, 3> TR2Size = {36,150,TR2Thickness};
static constexpr double TR2GapX = 2; 
//layers
static constexpr int PixelNCols = 1024;
static constexpr int PixelNRows = 512;
static constexpr double PixelSizeCols = 0.02924;
static constexpr double PixelSizeRows = 0.02688;
static constexpr double ChipSizeX = PixelSizeCols*PixelNCols;
static constexpr double ChipSizeY = PixelSizeRows*PixelNRows;
static constexpr double ChipDistanceX = 0.150;
static constexpr double ChipDistanceY = 0.150;
static constexpr double ChipStaveDistanceY = 7.22312;
static constexpr std::array<double, 3> StaveZ = {17.825+TR1CenterZ,17.825+8.5+TR1CenterZ,17.825+17+TR1CenterZ};
//useful variables
int hmgt = 0;           //how many generated tracks
int hmgthTR1 = 0;       //how many generated tracks hitted TR1
int hmgthL2 = 0;        //how many generated tracks hitted L2
int hmgthL1 = 0;        //how many generated tracks hitted L1
int hmgthL0 = 0;        //how many generated tracks hitted L0
double err_cl = 1;      //errore cluster fisso atm but have to change later on
bool track_generation = true;    //true: generate tracks joining pTR2-qTR1 //false: generate tracks pTR2,theta,phi with appropriate distribution




//---------------------------------------------------------------------------------------------------------------------------

void geometry() {

    TCanvas* geom = new TCanvas("geom", "3D View", 800, 600);
    TView* display = TView::CreateView(1);
    display->SetRange(-100, -100, 0, 100, 100, 70);
    display->ShowAxis();

    //draw trigger 1
    TMarker3DBox *TR1_0 = new TMarker3DBox(0,0,0,TR1Size[0]/2,TR1Size[1]/2,TR1Size[2]/2,0,0);
    TR1_0->SetLineColor(kRed);
    TR1_0->SetLineWidth(3);
    TR1_0->Draw();
    TMarker3DBox *TR1_1 = new TMarker3DBox(0,TR1Size[1]+TR1GapY,0,TR1Size[0]/2,TR1Size[1]/2,TR1Size[2]/2,0,0);
    TR1_1->SetLineColor(kRed);
    TR1_1->SetLineWidth(3);
    TR1_1->Draw();
    TMarker3DBox *TR1_2 = new TMarker3DBox(0,2*(TR1Size[1]+TR1GapY),0,TR1Size[0]/2,TR1Size[1]/2,TR1Size[2]/2,0,0);
    TR1_2->SetLineColor(kRed);
    TR1_2->SetLineWidth(3);
    TR1_2->Draw();
    TMarker3DBox *TR1_3 = new TMarker3DBox(0,-(TR1Size[1]+TR1GapY),0,TR1Size[0]/2,TR1Size[1]/2,TR1Size[2]/2,0,0);
    TR1_3->SetLineColor(kRed);
    TR1_3->SetLineWidth(3);
    TR1_3->Draw();
    TMarker3DBox *TR1_4 = new TMarker3DBox(0,-2*(TR1Size[1]+TR1GapY),0,TR1Size[0]/2,TR1Size[1]/2,TR1Size[2]/2,0,0);
    TR1_4->SetLineColor(kRed);
    TR1_4->SetLineWidth(3);
    TR1_4->Draw();

    //draw trigger 2
    TMarker3DBox *TR2_0 = new TMarker3DBox(-1.5*(TR2Size[0]+TR2GapX),0,TR2CenterZ, TR2Size[0]/2, TR2Size[1]/2, TR2Size[2]/2, 0, 0);
    TR2_0->SetLineColor(kRed);
    TR2_0->SetLineWidth(3);
    TR2_0->Draw();
    TMarker3DBox *TR2_1 = new TMarker3DBox(-0.5*(TR2Size[0]+TR2GapX),0,TR2CenterZ, TR2Size[0]/2, TR2Size[1]/2, TR2Size[2]/2, 0, 0);
    TR2_1->SetLineColor(kRed);
    TR2_1->SetLineWidth(3);
    TR2_1->Draw();
    TMarker3DBox *TR2_2 = new TMarker3DBox(0.5*(TR2Size[0]+TR2GapX),0,TR2CenterZ, TR2Size[0]/2, TR2Size[1]/2, TR2Size[2]/2, 0, 0);
    TR2_2->SetLineColor(kRed);
    TR2_2->SetLineWidth(3);
    TR2_2->Draw();
    TMarker3DBox *TR2_3 = new TMarker3DBox(1.5*(TR2Size[0]+TR2GapX),0,TR2CenterZ, TR2Size[0]/2, TR2Size[1]/2, TR2Size[2]/2, 0, 0);
    TR2_3->SetLineColor(kRed);
    TR2_3->SetLineWidth(3);
    TR2_3->Draw();

    //not sure about the thickness

    //draw layer 1
    //row, cols (00) 
    //row0
    TMarker3DBox *stave1_00 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_00->Draw();
    TMarker3DBox *stave1_01 = new TMarker3DBox(-1*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_01->Draw();
    TMarker3DBox *stave1_02 = new TMarker3DBox(0,-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_02->Draw();
    TMarker3DBox *stave1_03 = new TMarker3DBox(1*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_03->Draw();
    TMarker3DBox *stave1_04 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_04->Draw();
    //row1
    TMarker3DBox *stave1_10 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_10->Draw();
    TMarker3DBox *stave1_11 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_11->Draw();
    TMarker3DBox *stave1_12 = new TMarker3DBox(0,-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_12->Draw();
    TMarker3DBox *stave1_13 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_13->Draw();
    TMarker3DBox *stave1_14 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_14->Draw();
    //row2
    TMarker3DBox *stave1_20 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_20->Draw();
    TMarker3DBox *stave1_21 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_21->Draw();
    TMarker3DBox *stave1_22 = new TMarker3DBox(0,-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_22->Draw();
    TMarker3DBox *stave1_23 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_23->Draw();
    TMarker3DBox *stave1_24 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_24->Draw();
    //row3
    TMarker3DBox *stave1_30 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_30->Draw();
    TMarker3DBox *stave1_31 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_31->Draw();
    TMarker3DBox *stave1_32 = new TMarker3DBox(0,-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_32->Draw();
    TMarker3DBox *stave1_33 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_33->Draw();
    TMarker3DBox *stave1_34 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_34->Draw();
    //row4
    TMarker3DBox *stave1_40 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_40->Draw();
    TMarker3DBox *stave1_41 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_41->Draw();
    TMarker3DBox *stave1_42 = new TMarker3DBox(0,-0.5*(ChipDistanceY+ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_42->Draw();
    TMarker3DBox *stave1_43 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_43->Draw();
    TMarker3DBox *stave1_44 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_44->Draw();
    //row5
    TMarker3DBox *stave1_50 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_50->Draw();
    TMarker3DBox *stave1_51 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_51->Draw();
    TMarker3DBox *stave1_52 = new TMarker3DBox(0,+0.5*(ChipDistanceY+ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_52->Draw();
    TMarker3DBox *stave1_53 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_53->Draw();
    TMarker3DBox *stave1_54 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_54->Draw();
    //row6
    TMarker3DBox *stave1_60 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_60->Draw();
    TMarker3DBox *stave1_61 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_61->Draw();
    TMarker3DBox *stave1_62 = new TMarker3DBox(0,+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_62->Draw();
    TMarker3DBox *stave1_63 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_63->Draw();
    TMarker3DBox *stave1_64 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_64->Draw();
    //row7
    TMarker3DBox *stave1_70 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_70->Draw();
    TMarker3DBox *stave1_71 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_71->Draw();
    TMarker3DBox *stave1_72 = new TMarker3DBox(0,+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_72->Draw();
    TMarker3DBox *stave1_73 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_73->Draw();
    TMarker3DBox *stave1_74 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_74->Draw();
    //row8
    TMarker3DBox *stave1_80 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_80->Draw();
    TMarker3DBox *stave1_81 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_81->Draw();
    TMarker3DBox *stave1_82 = new TMarker3DBox(0,+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_82->Draw();
    TMarker3DBox *stave1_83 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_83->Draw();
    TMarker3DBox *stave1_84 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_84->Draw();
    //row9
    TMarker3DBox *stave1_90 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_90->Draw();
    TMarker3DBox *stave1_91 = new TMarker3DBox(-1*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_91->Draw();
    TMarker3DBox *stave1_92 = new TMarker3DBox(0,+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_92->Draw();
    TMarker3DBox *stave1_93 = new TMarker3DBox(1*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_93->Draw();
    TMarker3DBox *stave1_94 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_94->Draw();



    //layer2
    //row0
    TMarker3DBox *stave2_00 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_00->Draw();
    TMarker3DBox *stave2_01 = new TMarker3DBox(-1*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_01->Draw();
    TMarker3DBox *stave2_02 = new TMarker3DBox(0,-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_02->Draw();
    TMarker3DBox *stave2_03 = new TMarker3DBox(1*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_03->Draw();
    TMarker3DBox *stave2_04 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_04->Draw();
    //row1
    TMarker3DBox *stave2_10 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_10->Draw();
    TMarker3DBox *stave2_11 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_11->Draw();
    TMarker3DBox *stave2_12 = new TMarker3DBox(0,-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_12->Draw();
    TMarker3DBox *stave2_13 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_13->Draw();
    TMarker3DBox *stave2_14 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_14->Draw();
    //row2
    TMarker3DBox *stave2_20 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_20->Draw();
    TMarker3DBox *stave2_21 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_21->Draw();
    TMarker3DBox *stave2_22 = new TMarker3DBox(0,-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_22->Draw();
    TMarker3DBox *stave2_23 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_23->Draw();
    TMarker3DBox *stave2_24 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_24->Draw();
    //row3
    TMarker3DBox *stave2_30 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_30->Draw();
    TMarker3DBox *stave2_31 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_31->Draw();
    TMarker3DBox *stave2_32 = new TMarker3DBox(0,-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_32->Draw();
    TMarker3DBox *stave2_33 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_33->Draw();
    TMarker3DBox *stave2_34 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_34->Draw();
    //row4
    TMarker3DBox *stave2_40 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_40->Draw();
    TMarker3DBox *stave2_41 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_41->Draw();
    TMarker3DBox *stave2_42 = new TMarker3DBox(0,-0.5*(ChipDistanceY+ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_42->Draw();
    TMarker3DBox *stave2_43 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_43->Draw();
    TMarker3DBox *stave2_44 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_44->Draw();
    //row5
    TMarker3DBox *stave2_50 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_50->Draw();
    TMarker3DBox *stave2_51 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_51->Draw();
    TMarker3DBox *stave2_52 = new TMarker3DBox(0,+0.5*(ChipDistanceY+ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_52->Draw();
    TMarker3DBox *stave2_53 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_53->Draw();
    TMarker3DBox *stave2_54 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_54->Draw();
    //row6
    TMarker3DBox *stave2_60 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_60->Draw();
    TMarker3DBox *stave2_61 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_61->Draw();
    TMarker3DBox *stave2_62 = new TMarker3DBox(0,+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_62->Draw();
    TMarker3DBox *stave2_63 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_63->Draw();
    TMarker3DBox *stave2_64 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_64->Draw();
    //row7
    TMarker3DBox *stave2_70 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_70->Draw();
    TMarker3DBox *stave2_71 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_71->Draw();
    TMarker3DBox *stave2_72 = new TMarker3DBox(0,+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_72->Draw();
    TMarker3DBox *stave2_73 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_73->Draw();
    TMarker3DBox *stave2_74 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_74->Draw();
    //row8
    TMarker3DBox *stave2_80 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_80->Draw();
    TMarker3DBox *stave2_81 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_81->Draw();
    TMarker3DBox *stave2_82 = new TMarker3DBox(0,+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_82->Draw();
    TMarker3DBox *stave2_83 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_83->Draw();
    TMarker3DBox *stave2_84 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_84->Draw();
    //row9
    TMarker3DBox *stave2_90 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_90->Draw();
    TMarker3DBox *stave2_91 = new TMarker3DBox(-1*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_91->Draw();
    TMarker3DBox *stave2_92 = new TMarker3DBox(0,+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_92->Draw();
    TMarker3DBox *stave2_93 = new TMarker3DBox(1*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_93->Draw();
    TMarker3DBox *stave2_94 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_94->Draw();

    

    //layer3
    //row0
    TMarker3DBox *stave3_00 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_00->Draw();
    TMarker3DBox *stave3_01 = new TMarker3DBox(-1*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_01->Draw();
    TMarker3DBox *stave3_02 = new TMarker3DBox(0,-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_02->Draw();
    TMarker3DBox *stave3_03 = new TMarker3DBox(1*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_03->Draw();
    TMarker3DBox *stave3_04 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_04->Draw();
    //row1
    TMarker3DBox *stave3_10 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_10->Draw();
    TMarker3DBox *stave3_11 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_11->Draw();
    TMarker3DBox *stave3_12 = new TMarker3DBox(0,-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_12->Draw();
    TMarker3DBox *stave3_13 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_13->Draw();
    TMarker3DBox *stave3_14 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_14->Draw();
    //row2
    TMarker3DBox *stave3_20 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_20->Draw();
    TMarker3DBox *stave3_21 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_21->Draw();
    TMarker3DBox *stave3_22 = new TMarker3DBox(0,-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_22->Draw();
    TMarker3DBox *stave3_23 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_23->Draw();
    TMarker3DBox *stave3_24 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_24->Draw();
    //row3
    TMarker3DBox *stave3_30 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_30->Draw();
    TMarker3DBox *stave3_31 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_31->Draw();
    TMarker3DBox *stave3_32 = new TMarker3DBox(0,-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_32->Draw();
    TMarker3DBox *stave3_33 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_33->Draw();
    TMarker3DBox *stave3_34 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_34->Draw();
    //row4
    TMarker3DBox *stave3_40 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_40->Draw();
    TMarker3DBox *stave3_41 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_41->Draw();
    TMarker3DBox *stave3_42 = new TMarker3DBox(0,-0.5*(ChipDistanceY+ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_42->Draw();
    TMarker3DBox *stave3_43 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_43->Draw();
    TMarker3DBox *stave3_44 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_44->Draw();
    //row5
    TMarker3DBox *stave3_50 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_50->Draw();
    TMarker3DBox *stave3_51 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_51->Draw();
    TMarker3DBox *stave3_52 = new TMarker3DBox(0,+0.5*(ChipDistanceY+ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_52->Draw();
    TMarker3DBox *stave3_53 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_53->Draw();
    TMarker3DBox *stave3_54 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_54->Draw();
    //row6
    TMarker3DBox *stave3_60 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_60->Draw();
    TMarker3DBox *stave3_61 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_61->Draw();
    TMarker3DBox *stave3_62 = new TMarker3DBox(0,+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_62->Draw();
    TMarker3DBox *stave3_63 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_63->Draw();
    TMarker3DBox *stave3_64 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_64->Draw();
    //row7
    TMarker3DBox *stave3_70 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_70->Draw();
    TMarker3DBox *stave3_71 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_71->Draw();
    TMarker3DBox *stave3_72 = new TMarker3DBox(0,+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_72->Draw();
    TMarker3DBox *stave3_73 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_73->Draw();
    TMarker3DBox *stave3_74 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_74->Draw();
    //row8
    TMarker3DBox *stave3_80 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_80->Draw();
    TMarker3DBox *stave3_81 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_81->Draw();
    TMarker3DBox *stave3_82 = new TMarker3DBox(0,+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_82->Draw();
    TMarker3DBox *stave3_83 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_83->Draw();
    TMarker3DBox *stave3_84 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_84->Draw();
    //row9
    TMarker3DBox *stave3_90 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_90->Draw();
    TMarker3DBox *stave3_91 = new TMarker3DBox(-1*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_91->Draw();
    TMarker3DBox *stave3_92 = new TMarker3DBox(0,+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_92->Draw();
    TMarker3DBox *stave3_93 = new TMarker3DBox(1*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_93->Draw();
    TMarker3DBox *stave3_94 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_94->Draw();



    //TH1F(name, title, nbins, xlow, xup)
    TH1F* hxTR2 = new TH1F("hxTR2", "x_trigger2", events, -TR2Size[0]*2.5, TR2Size[0]*2.5);
    TH1F* hyTR2 = new TH1F("hyTR2", "y_trigger2 ", events, -TR2Size[1]/2, TR2Size[1]/2);
    TH1F* hzTR2 = new TH1F("hzTR2", "z_trigger2 ", events, -TR2Size[2]/2+TR2CenterZ, TR2Size[2]/2+TR2CenterZ);
    TH1F* hxTR1 = new TH1F("hxTR1", "x_trigger1", events, -TR1Size[0]/2, TR1Size[0]/2);
    TH1F* hyTR1 = new TH1F("hyTR1", "y_trigger1", events, -TR1Size[1]*3, TR1Size[1]*3);
    TH1F* hzTR1 = new TH1F("hzTR1", "z_trigger1 ", events, -TR1Size[2]/2+TR1CenterZ, TR1Size[2]/2+TR1CenterZ);
    TH1F* hphi = new TH1F("hphi", "phi", events, -pi, pi);
    TH1F* htheta = new TH1F("htheta", "theta", events, -pi/2, pi/2);

    //puoi settare il seed for reproducibility
    TRandom3 *rnd = new TRandom3(10); 

    //MC
    for (int i=0; i < events; i++,hmgt++){
        double xTR2, yTR2, zTR2, xTR1, yTR1, zTR1, phi, theta;

        //i layer acquisiscono il segnale quando arriva un segnale AND dagli scintillatori
        //genero traccie misurabili come traccie che passano nei due scintillatori TR2, TR1
        if(track_generation){
            xTR2 = -100;
            double xTR2_fake = rnd->Uniform(-TR2Size[0]*2,TR2Size[0]*2);     
            if(xTR2_fake>0 && xTR2_fake<TR2Size[0]){xTR2 = xTR2_fake + 0.5*TR2GapX;}
            if(xTR2_fake<0 && xTR2_fake>-TR2Size[0]){xTR2 = xTR2_fake - 0.5*TR2GapX;}
            if(xTR2_fake<2*TR2Size[0] && xTR2_fake>TR2Size[0]){xTR2 = xTR2_fake + 1.5*TR2GapX;}
            if(xTR2_fake>-2*TR2Size[0] && xTR2_fake<-TR2Size[0]){xTR2 = xTR2_fake - 1.5*TR2GapX;}
            yTR2 = rnd->Uniform(-TR2Size[1]/2,TR2Size[1]/2);
            zTR2 = rnd->Uniform(TR2CenterZ-TR2Thickness/2,TR2CenterZ+TR2Thickness/2);

            yTR1 = -100;
            xTR1 = rnd->Uniform(-TR1Size[0]/2,TR1Size[0]/2);
            double yTR1_fake = rnd->Uniform(-TR1Size[1]*2.5,TR1Size[1]*2.5);     
            if(yTR1_fake<TR1Size[1]/2 && yTR1_fake>-TR1Size[1]/2){yTR1 = yTR1_fake;}
            if(yTR1_fake<1.5*TR1Size[1] && yTR1_fake>0.5*TR1Size[1]){yTR1 = yTR1_fake + TR1GapY;}
            if(yTR1_fake<2.5*TR1Size[1] && yTR1_fake>1.5*TR1Size[1]){yTR1 = yTR1_fake + 2*TR1GapY;}
            if(yTR1_fake>-1.5*TR1Size[1] && yTR1_fake<-0.5*TR1Size[1]){yTR1 = yTR1_fake - TR1GapY;}
            if(yTR1_fake>-2.5*TR1Size[1] && yTR1_fake<-1.5*TR1Size[1]){yTR1 = yTR1_fake - 2*TR1GapY;}
            zTR1 = rnd->Uniform(TR1CenterZ-TR1Thickness/2,TR1CenterZ+TR1Thickness/2);

            phi = TMath::ATan((yTR2-yTR1)/(xTR2-xTR1));
            theta = TMath::ATan((xTR1-xTR2)/((TMath::Cos(phi))*(zTR2-zTR1)));   

            hmgthTR1++;
        }  
        if(!track_generation){
            xTR2 = -100;
            double xTR2_fake = rnd->Uniform(-TR2Size[0]*2,TR2Size[0]*2);     
            if(xTR2_fake>0 && xTR2_fake<TR2Size[0]){xTR2 = xTR2_fake + 0.5*TR2GapX;}
            if(xTR2_fake<0 && xTR2_fake>-TR2Size[0]){xTR2 = xTR2_fake - 0.5*TR2GapX;}
            if(xTR2_fake<2*TR2Size[0] && xTR2_fake>TR2Size[0]){xTR2 = xTR2_fake + 1.5*TR2GapX;}
            if(xTR2_fake>-2*TR2Size[0] && xTR2_fake<-TR2Size[0]){xTR2 = xTR2_fake - 1.5*TR2GapX;}
            yTR2 = rnd->Uniform(-TR2Size[1]/2,TR2Size[1]/2);
            zTR2 = rnd->Uniform(TR2CenterZ-TR2Thickness/2,TR2CenterZ+TR2Thickness/2);

            phi = rnd->Uniform(pi)-pi/2;
            double THETA = rnd->Uniform(pi)-pi/2;
            theta = pow(TMath::Cos(THETA),2); 

            xTR1 = xTR2 + (zTR2-TR1CenterZ)*(TMath::Sin(theta))*(TMath::Cos(phi))*(1/(TMath::Cos(theta)));
            yTR1 = yTR2 + (zTR2-TR1CenterZ)*(TMath::Sin(theta))*(TMath::Sin(phi)*(1/(TMath::Cos(theta))));
            //check how many tracks hitted TR1
            if(xTR1 < TR1Size[0]/2 && xTR1 > -TR1Size[0]/2 &&
                (
                    (yTR1 < (2.5*TR1Size[1]+2*TR1GapY) && yTR1 > (1.5*TR1Size[1]+2*TR1GapY)) ||
                    (yTR1 < (1.5*TR1Size[1]+1*TR1GapY) && yTR1 > (0.5*TR1Size[1]+1*TR1GapY)) ||
                    (yTR1 < (0.5*TR1Size[1]+0*TR1GapY) && yTR1 > -(0.5*TR1Size[1]+0*TR1GapY)) ||
                    (yTR1 < -(0.5*TR1Size[1]+1*TR1GapY) && yTR1 > -(1.5*TR1Size[1]+1*TR1GapY)) ||
                    (yTR1 < -(1.5*TR1Size[1]+2*TR1GapY) && yTR1 > -(2.5*TR1Size[1]+2*TR1GapY)) 
                )
              ){
                hmgthTR1++;
            }
        }
        hxTR2->Fill(xTR2);
        hyTR2->Fill(yTR2);
        hzTR2->Fill(zTR2);
        hxTR1->Fill(xTR1);
        hyTR1->Fill(yTR1);
        hzTR1->Fill(zTR1);
        hphi->Fill(phi);
        htheta->Fill(theta);

        //plotting tracks
        Double_t x_line[2] = {xTR2, xTR1};
        Double_t y_line[2] = {yTR2, yTR1};
        Double_t z_line[2] = {zTR2, zTR1};
        TPolyLine3D* line_track = new TPolyLine3D(2, x_line, y_line, z_line);
        line_track->SetLineColor(kBlue);
        line_track->SetLineWidth(2);
        line_track->Draw();

        //calcolimo intersezione con i piani
        double xL2 = xTR2 + (zTR2-StaveZ[2])*(TMath::Sin(theta))*(TMath::Cos(phi))*(1/(TMath::Cos(theta)));
        double yL2 = yTR2 + (zTR2-StaveZ[2])*(TMath::Sin(theta))*(TMath::Sin(phi))*(1/(TMath::Cos(theta)));
        double zL2 = StaveZ[2];
        double xL1 = xTR2 + (zTR2-StaveZ[1])*(TMath::Sin(theta))*(TMath::Cos(phi))*(1/(TMath::Cos(theta)));
        double yL1 = yTR2 + (zTR2-StaveZ[1])*(TMath::Sin(theta))*(TMath::Sin(phi))*(1/(TMath::Cos(theta)));
        double zL1 = StaveZ[1];
        double xL0 = xTR2 + (zTR2-StaveZ[0])*(TMath::Sin(theta))*(TMath::Cos(phi))*(1/(TMath::Cos(theta)));
        double yL0 = yTR2 + (zTR2-StaveZ[0])*(TMath::Sin(theta))*(TMath::Sin(phi))*(1/(TMath::Cos(theta)));
        double zL0 = StaveZ[0];

        //check if the track hitted the staves in layer 2
        if((xL2 < ChipSizeX*2.5 + ChipDistanceX && xL2 > -(ChipSizeX*2.5 + ChipDistanceX)) &&
           ((yL2 < ChipSizeY*5 + ChipStaveDistanceY*2 + ChipDistanceY*2.5 && yL2 > ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*2.5) ||
            (yL2 < ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*1.5 && yL2 > ChipSizeY*3 + ChipStaveDistanceY*2 + ChipDistanceY*1.5) ||
            (yL2 < ChipSizeY*3 + ChipStaveDistanceY*1 + ChipDistanceY*1.5 && yL2 > ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*1.5) ||
            (yL2 < ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*0.5 && yL2 > ChipSizeY*1 + ChipStaveDistanceY*1 + ChipDistanceY*0.5) ||
            (yL2 < ChipSizeY*1 + ChipStaveDistanceY*0 + ChipDistanceY*0.5 && yL2 > ChipSizeY*0 + ChipStaveDistanceY*0 + ChipDistanceY*0.5) ||
            (yL2 > -(ChipSizeY*1 + ChipStaveDistanceY*0 + ChipDistanceY*0.5) && yL2 < -(ChipSizeY*0 + ChipStaveDistanceY*0 + ChipDistanceY*0.5)) ||
            (yL2 > -(ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*0.5) && yL2 < -(ChipSizeY*1 + ChipStaveDistanceY*1 + ChipDistanceY*0.5)) ||
            (yL2 > -(ChipSizeY*3 + ChipStaveDistanceY*1 + ChipDistanceY*1.5) && yL2 < -(ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*1.5)) ||
            (yL2 > -(ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*1.5) && yL2 < -(ChipSizeY*3 + ChipStaveDistanceY*2 + ChipDistanceY*1.5)) ||
            (yL2 > -(ChipSizeY*5 + ChipStaveDistanceY*2 + ChipDistanceY*2.5) && yL2 < -(ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*2.5)))
        ){
            hmgthL2++;
        }
        if((xL1 < ChipSizeX*2.5 + ChipDistanceX && xL1 > -(ChipSizeX*2.5 + ChipDistanceX)) &&
           ((yL1 < ChipSizeY*5 + ChipStaveDistanceY*2 + ChipDistanceY*2.5 && yL1 > ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*2.5) ||
            (yL1 < ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*1.5 && yL1 > ChipSizeY*3 + ChipStaveDistanceY*2 + ChipDistanceY*1.5) ||
            (yL1 < ChipSizeY*3 + ChipStaveDistanceY*1 + ChipDistanceY*1.5 && yL1 > ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*1.5) ||
            (yL1 < ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*0.5 && yL1 > ChipSizeY*1 + ChipStaveDistanceY*1 + ChipDistanceY*0.5) ||
            (yL1 < ChipSizeY*1 + ChipStaveDistanceY*0 + ChipDistanceY*0.5 && yL1 > ChipSizeY*0 + ChipStaveDistanceY*0 + ChipDistanceY*0.5) ||
            (yL1 > -(ChipSizeY*1 + ChipStaveDistanceY*0 + ChipDistanceY*0.5) && yL1 < -(ChipSizeY*0 + ChipStaveDistanceY*0 + ChipDistanceY*0.5)) ||
            (yL1 > -(ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*0.5) && yL1 < -(ChipSizeY*1 + ChipStaveDistanceY*1 + ChipDistanceY*0.5)) ||
            (yL1 > -(ChipSizeY*3 + ChipStaveDistanceY*1 + ChipDistanceY*1.5) && yL1 < -(ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*1.5)) ||
            (yL1 > -(ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*1.5) && yL1 < -(ChipSizeY*3 + ChipStaveDistanceY*2 + ChipDistanceY*1.5)) ||
            (yL1 > -(ChipSizeY*5 + ChipStaveDistanceY*2 + ChipDistanceY*2.5) && yL1 < -(ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*2.5)))
        ){
            hmgthL1++;
        }
        if((xL0 < ChipSizeX*2.5 + ChipDistanceX && xL0 > -(ChipSizeX*2.5 + ChipDistanceX)) &&
           ((yL0 < ChipSizeY*5 + ChipStaveDistanceY*2 + ChipDistanceY*2.5 && yL0 > ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*2.5) ||
            (yL0 < ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*1.5 && yL0 > ChipSizeY*3 + ChipStaveDistanceY*2 + ChipDistanceY*1.5) ||
            (yL0 < ChipSizeY*3 + ChipStaveDistanceY*1 + ChipDistanceY*1.5 && yL0 > ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*1.5) ||
            (yL0 < ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*0.5 && yL0 > ChipSizeY*1 + ChipStaveDistanceY*1 + ChipDistanceY*0.5) ||
            (yL0 < ChipSizeY*1 + ChipStaveDistanceY*0 + ChipDistanceY*0.5 && yL0 > ChipSizeY*0 + ChipStaveDistanceY*0 + ChipDistanceY*0.5) ||
            (yL0 > -(ChipSizeY*1 + ChipStaveDistanceY*0 + ChipDistanceY*0.5) && yL0 < -(ChipSizeY*0 + ChipStaveDistanceY*0 + ChipDistanceY*0.5)) ||
            (yL0 > -(ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*0.5) && yL0 < -(ChipSizeY*1 + ChipStaveDistanceY*1 + ChipDistanceY*0.5)) ||
            (yL0 > -(ChipSizeY*3 + ChipStaveDistanceY*1 + ChipDistanceY*1.5) && yL0 < -(ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*1.5)) ||
            (yL0 > -(ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*1.5) && yL0 < -(ChipSizeY*3 + ChipStaveDistanceY*2 + ChipDistanceY*1.5)) ||
            (yL0 > -(ChipSizeY*5 + ChipStaveDistanceY*2 + ChipDistanceY*2.5) && yL0 < -(ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*2.5)))
        ){
            hmgthL0++;
        }

        /*
        TMarker3DBox *err_L2 = new TMarker3DBox(xL2,yL2,zL2,err_cl,err_cl,0,0,0);
        err_L2->SetLineWidth(2);
        err_L2->SetLineColor(kGreen);
        err_L2->Draw();

        TMarker3DBox *err_L1 = new TMarker3DBox(xL1,yL1,zL1,err_cl,err_cl,0,0,0);
        err_L1->SetLineWidth(2);
        err_L1->SetLineColor(kBlue);
        err_L1->Draw();

        TMarker3DBox *err_L0 = new TMarker3DBox(xL0,yL0,zL0,err_cl,err_cl,0,0,0);
        err_L0->SetLineWidth(2);
        err_L0->SetLineColor(kRed);
        err_L0->Draw();
        */

        
    }


    char file[200];
    sprintf(file,"geom");
    geom->SaveAs(file);
    //output mc distributions
    sprintf(file,"geomHist.root");
    TFile f(file,"RECREATE");
    hxTR2->Write();
    hyTR2->Write();
    hzTR2->Write();
    hxTR1->Write();
    hyTR1->Write();
    hzTR1->Write();
    hphi->Write();
    htheta->Write();
    f.Close();

    //double phi_max = pi/2;
    //double phi_thetamax = TMath::ATan((-TR2Size[1]-(2.5*TR1Size[1]+2*TR1GapY))/(-(2*TR2Size[0]+1.5*TR2GapX)-TR1Size[0]/2));
    //double theta_max = TMath::ATan((-(-(2*TR2Size[0]+1.5*TR2GapX)-TR1Size[0]/2))/(TMath::Cos(phi_thetamax)*(TR2CenterZ+TR2Size[2]/2)-TR1CenterZ-TR1Size[2]/2));
    if(track_generation){
        cout << "GENERIAMO LE TRACK UNENDO pTR2-qTR1" << endl;
    }
    if(!track_generation){
        cout << "GENERIAMO LE TRACK DA pTR2 CON DISTRIBUZIONE DEGLI ANGOLI" << endl;
    }
    //cout << "angolo phi max: " << phi_max << endl;
    //cout << "angolo theta max: " << theta_max << endl;
    cout << "how many generated tracks: " << hmgt << endl;
    cout << "how many generated tracks hitted TR1: " << hmgthTR1 << endl;
    cout << "how many generated tracks hitted L2: " << hmgthL2 << endl;
    cout << "how many generated tracks hitted L1: " << hmgthL1 << endl;
    cout << "how many generated tracks hitted L0: " << hmgthL0 << endl;


}




