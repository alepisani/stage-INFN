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

double pi = 3.141592;
double h = 8.5; //mm
float shift_zB1 = 17.825;                  // [mm] distance between trigger layer and firts tracking layer
static constexpr double TR1CenterZ = 0.;
static constexpr double TR2CenterZ = TR1CenterZ+60.5;
float shift_z2T = TR2CenterZ - shift_zB1 - 2*h;                  // [mm] distance between trigger layer and firts tracking layer
std::array<double, 3> StaveZ = { shift_zB1, shift_zB1 + h, shift_zB1 + 2*h};
float x_lay = 13.76256*10; //mm
float y_lay = 29.94176*5; //mm
float t_lay = 0.; //mm
static constexpr double TR1Thickness = 2;
static constexpr double TR2Thickness = 8;
static constexpr std::array<double, 3> TR1Size = {154.6,32.5*5,TR1Thickness};
static constexpr std::array<double, 3> TR2Size = {36*5,150,TR2Thickness};


float offset = 0.;
vector<float> c_lay0 = { offset, offset, 0. };
int iteration = 100;
TRandom3 *rnd = new TRandom3(0);            //puoi settare il seed for reproducibility
double err_cl = 1. ;
vector<vector<double>> clusters_lay2 = {};
vector<vector<double>> clusters_lay1 = {};
vector<vector<double>> clusters_lay0 = {};
vector<vector<vector<double>>> real_tracks = {};
int howmanyrecotracks = 0;
int howmanyrealtracks = 0;

//turret composto da 3 stave sovrapposti
//rivelatore composto da 5 turret affiancati


//--------------------------------------------------------------------------------

void tracks() {

    cout << "offset: " << offset << endl;
    cout << "MC tracks: " << iteration << endl;
   
    TCanvas* canvas = new TCanvas("canvas", "3D View", 800, 600);
    TView* layer_0 = TView::CreateView(1);
    layer_0->SetRange(0, 0, 0, 200, 200, 70);
    layer_0->ShowAxis();
    
    TCanvas* hyp = new TCanvas("hyp", "3D View", 800, 600);
    TView* hp_tracks = TView::CreateView(1);
    hp_tracks->SetRange(0, 0, 0, 200, 200, 70);
    hp_tracks->ShowAxis();


    canvas->cd();
    //nb TMarker3DBox le cordinate date sono il centro del box disegnato e le dimensioni la meta dei lati reali

    //trigger levels
    TMarker3DBox *TR1 = new TMarker3DBox(TR1Size[0]/2, TR1Size[1]/2, TR1CenterZ, TR1Size[0]/2, TR1Size[1]/2, 0, 0, 0);
    TMarker3DBox *TR2 = new TMarker3DBox(TR2Size[0]/2, TR2Size[1]/2, TR2CenterZ, TR2Size[0]/2, TR2Size[1]/2, 0, 0, 0);   
                    //magari disegna ogni singolo trigger, quelli disegnati ora sono tutti e 5 insieme, draw edges
    //layers  
    TMarker3DBox *layer0 = new TMarker3DBox(x_lay/2, y_lay/2, StaveZ[2], x_lay/2, y_lay/2, 0, 0, 0); //x,y,z coordinates //x,y,z dimension //theta,phi    
    TMarker3DBox *layer1 = new TMarker3DBox(x_lay/2, y_lay/2, StaveZ[1], x_lay/2, y_lay/2, 0, 0, 0); //x,y,z coordinates //x,y,z dimension //theta,phi    
    TMarker3DBox *layer2 = new TMarker3DBox(x_lay/2, y_lay/2, StaveZ[0], x_lay/2, y_lay/2, 0, 0, 0); //x,y,z coordinates //x,y,z dimension //theta,phi  
    TR1->SetLineColor(kRed);
    TR2->SetLineColor(kRed);
    layer0->SetLineColor(kRed);
    layer1->SetLineColor(kRed);
    layer2->SetLineColor(kRed);
    TR1->SetLineWidth(4);
    TR2->SetLineWidth(4);
    layer0->SetLineWidth(4);
    layer1->SetLineWidth(4);
    layer2->SetLineWidth(4);  
    TR1->Draw();
    TR2->Draw();
    layer0->Draw(); 
    layer1->Draw();      
    layer2->Draw();
    hyp->cd();
    TR1->Draw();
    TR2->Draw();
    layer0->Draw(); 
    layer1->Draw();
    layer2->Draw();
    canvas->cd();

    //TH1F(name, title, nbins, xlow, xup)
    TH1F* hxTR2 = new TH1F("hxTR2", "x_trigger2", iteration, 0, TR2Size[0]);
    TH1F* hyTR2 = new TH1F("hyTR2", "y_trigger2 ", iteration, 0, TR2Size[1]);
    TH1F* hxTR1 = new TH1F("hxTR1", "x_trigger1", iteration, 0, TR1Size[0]);
    TH1F* hyTR1 = new TH1F("hyTR1", "y_trigger1", iteration, 0, TR1Size[1]);
    TH1F* hphi = new TH1F("hphi", "phi", iteration, -pi, pi);
    TH1F* htheta = new TH1F("htheta", "theta", iteration, 0, pi/2);
    TH1F* hTheta = new TH1F("hTheta", "Theta", iteration, 0, pi/2);

    //puoi settare il seed for reproducibility
    TRandom3 *rnd = new TRandom3(132517689); 

    //MC
    for (int i=0; i < iteration; i++){

        //genero le traccie a partire dal TR1
        double xTR2 = rnd->Uniform(TR2Size[0]);
        double yTR2 = rnd->Uniform(TR2Size[1]);
        double phi = rnd->Uniform(2.*pi)-pi;
        double THETA = rnd->Uniform(pi/2);
        double theta = pow(TMath::Cos(THETA),2); 
        double xTR1 = xTR2 + (TR2CenterZ-TR1CenterZ)*(TMath::Sin(theta))*(TMath::Cos(phi));
        double yTR1 = yTR2 + (TR2CenterZ-TR1CenterZ)*(TMath::Sin(theta))*(TMath::Sin(phi));
        hxTR2->Fill(xTR2);
        hyTR2->Fill(yTR2);
        hxTR1->Fill(xTR1);
        hyTR1->Fill(yTR1);
        hphi->Fill(phi);
        htheta->Fill(theta);
        hTheta->Fill(THETA);

        //genera punti sui different layer
        double xL2 = xTR2 + (TR2CenterZ-2*h-shift_zB1)*(TMath::Sin(theta))*(TMath::Cos(phi));
        double yL2 = yTR2 + (TR2CenterZ-2*h-shift_zB1)*(TMath::Sin(theta))*(TMath::Sin(phi));
        double zL2 = shift_zB1 + 2*h;
        double xL1 = xTR2 + (TR2CenterZ-h-shift_zB1)*(TMath::Sin(theta))*(TMath::Cos(phi));
        double yL1 = yTR2 + (TR2CenterZ-h-shift_zB1)*(TMath::Sin(theta))*(TMath::Sin(phi));
        double zL1 = shift_zB1 + h;
        double xL0 = xTR2 + (TR2CenterZ-shift_zB1)*(TMath::Sin(theta))*(TMath::Cos(phi));
        double yL0 = yTR2 + (TR2CenterZ-shift_zB1)*(TMath::Sin(theta))*(TMath::Sin(phi));
        double zL0 = shift_zB1;
         


        Double_t x_line[2] = {xTR2, xTR1};
        Double_t y_line[2] = {yTR2, yTR1};
        Double_t z_line[2] = {TR2CenterZ, TR1CenterZ};
        TPolyLine3D* line_track = new TPolyLine3D(2, x_line, y_line, z_line);
        line_track->SetLineColor(kBlue);
        line_track->SetLineWidth(1);
        line_track->Draw();
        howmanyrealtracks++;
        

        hyp->cd();
        TMarker3DBox *pb = new TMarker3DBox(xL2, yL2, shift_zB1 + 2*h, err_cl, err_cl, 0, 0, 0); 
        pb->SetLineColor(kBlack);
        if(xL2<=x_lay && xL2>0 && yL2<=y_lay && yL2>0){
            pb->Draw();
            canvas->cd();
            pb->Draw();
            hyp->cd();
        }
        TMarker3DBox *mb = new TMarker3DBox(xL1, yL1, shift_zB1 + h, err_cl, err_cl, 0, 0, 0);
        mb->SetLineColor(kBlack);
        if(xL1<=x_lay && xL1>0 && yL1<=y_lay && yL1>0){
            mb->Draw();
            canvas->cd();
            mb->Draw();
            hyp->cd();
        }
        TMarker3DBox *qb = new TMarker3DBox(xL0, yL0, shift_zB1, err_cl, err_cl, 0, 0, 0); 
        qb->SetLineColor(kBlack);
        if(xL0<=x_lay && xL0>0 && yL0<=y_lay && yL0>0){
            qb->Draw();
            canvas->cd();
            qb->Draw();
            hyp->cd();
        }
    canvas->cd();        

        //store point hitting the layer in vector
        //in ogni clusters_layN ci sono tutti i "pixel" attivi per quel dato layer
        vector<double> cl_lay2 = {xL2, yL2, shift_zB1 + 2*h};
        vector<double> cl_lay1 = {xL1, yL1, shift_zB1 + h};
        vector<double> cl_lay0 = {xL0, yL0, shift_zB1};
        if(xL2<=x_lay && xL2>0 && yL2<=y_lay && yL2>0){
            clusters_lay2.push_back(cl_lay2);
        }
        if(xL1<=x_lay && xL1>0 && yL1<=y_lay && yL1>0){
            clusters_lay1.push_back(cl_lay1);
        }
        if(xL0<=x_lay && xL0>0 && yL0<=y_lay && yL0>0){
            clusters_lay0.push_back(cl_lay0);
        }        
    }
    cout << "dim cl2:" << clusters_lay2.size() << endl; 
    cout << "dim cl1:" << clusters_lay1.size() << endl;
    cout << "dim cl0:" << clusters_lay0.size() << endl;


    
    //now, compute segment clusters_lay02 for each point and look for if the segment touches a cluster_lay1
    //nb con clusters_lay hai delle matrici nx3 con n quanti clusters ci sono su quel layer
    int size_cl_lay2 = clusters_lay2.size();
    int size_cl_lay1 = clusters_lay1.size();
    int size_cl_lay0 = clusters_lay0.size();

    if(size_cl_lay2 >= size_cl_lay0){
        for(int m=0; m<clusters_lay2.size(); m++){
            for (int n=0; n<clusters_lay0.size(); n++){
                double hp_xL1 = (clusters_lay2[m][0]+clusters_lay0[n][0])/2;
                double hp_yL1 = (clusters_lay2[m][1]+clusters_lay0[n][1])/2;
                for(int l=0; l<clusters_lay1.size(); l++){
                    if(hp_xL1 == clusters_lay1[l][0] && hp_yL1 == clusters_lay1[l][1]){
                        howmanyrecotracks++;
                        hyp->cd();
                        Double_t reco_xline[3] = {clusters_lay2[m][0], hp_xL1, clusters_lay0[n][0]};
                        Double_t reco_yline[3] = {clusters_lay2[m][1], hp_yL1, clusters_lay0[n][1]};
                        Double_t reco_zline[3] = {shift_zB1 + 2*h, shift_zB1 + h, shift_zB1};
                        TPolyLine3D* reco = new TPolyLine3D(3, reco_xline, reco_yline, reco_zline);
                        reco->SetLineColor(kBlue);  
                        reco->SetLineWidth(2);
                        reco->Draw();
                        TMarker3DBox *err0 = new TMarker3DBox(reco_xline[0], reco_yline[0], reco_zline[0], 0.5, 0.5, 0, 0, 0); 
                        err0->SetLineColor(kBlack);
                        err0->Draw();
                        TMarker3DBox *err1 = new TMarker3DBox(reco_xline[1], reco_yline[1], reco_zline[1], 0.5, 0.5, 0, 0, 0); 
                        err1->SetLineColor(kBlack);
                        err1->Draw();
                        TMarker3DBox *err2 = new TMarker3DBox(reco_xline[2], reco_yline[2], reco_zline[2], 0.5, 0.5, 0, 0, 0);
                        err2->SetLineColor(kBlack);
                        err2->Draw();
                        canvas->cd();
                    }
                }
            }
        }    
    }
    if(size_cl_lay0 > size_cl_lay2){
        for(int n=0; n<clusters_lay0.size(); n++){
            for (int m=0; m<clusters_lay2.size(); m++){
                double hp_xL1 = (clusters_lay2[m][0]+clusters_lay0[n][0])/2;
                double hp_yL1 = (clusters_lay2[m][1]+clusters_lay0[n][1])/2;
                for(int l=0; l<clusters_lay1.size(); l++){
                    if(hp_xL1 == clusters_lay1[l][0] && hp_yL1 == clusters_lay1[l][1]){
                        howmanyrecotracks++;
                        hyp->cd();
                        Double_t reco_xline[3] = {clusters_lay2[m][0], hp_xL1, clusters_lay0[n][0]};
                        Double_t reco_yline[3] = {clusters_lay2[m][1], hp_yL1, clusters_lay0[n][1]};
                        Double_t reco_zline[3] = {shift_zB1 + 2*h, shift_zB1 + h, shift_zB1};
                        TPolyLine3D* reco = new TPolyLine3D(3, reco_xline, reco_yline, reco_zline);
                        reco->SetLineColor(kBlue);  
                        reco->SetLineWidth(2);
                        reco->Draw();
                        TMarker3DBox *err0 = new TMarker3DBox(reco_xline[0], reco_yline[0], reco_zline[0], err_cl, err_cl, 0, 0, 0); 
                        err0->SetLineColor(kBlack);
                        err0->Draw();
                        TMarker3DBox *err1 = new TMarker3DBox(reco_xline[1], reco_yline[1], reco_zline[1], err_cl, err_cl, 0, 0, 0); 
                        err1->SetLineColor(kBlack);
                        err1->Draw();
                        TMarker3DBox *err2 = new TMarker3DBox(reco_xline[2], reco_yline[2], reco_zline[2], err_cl, err_cl, 0, 0, 0); 
                        err2->SetLineColor(kBlack);
                        err2->Draw();
                        canvas->cd();
                    }
                }
            }
        }
    }
    
    cout << "reco tracks:" << howmanyrecotracks << endl;
    cout << "real tracks:" << howmanyrealtracks << endl;

    //double acceptance = (double)npassedAcc/(double)events;
    //calcola accettanza geometrica
    

    //writing data into ttree to view them in tbrowser
    char file[200];
    sprintf(file,"tracks");
    canvas->SaveAs(file);
    hyp->SaveAs(file);
    //output mc distributions
    sprintf(file,"outputHists.root");
    TFile f(file,"RECREATE");
    hxTR2->Write();
    hyTR2->Write();
    hxTR1->Write();
    hyTR1->Write();
    hphi->Write();
    htheta->Write();
    hTheta->Write();
    f.Close();
}



