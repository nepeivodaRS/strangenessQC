#include "TStyle.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRatioPlot.h"
#include "TLegend.h"

void StyleCanvas(TCanvas *canvas, Float_t LMargin, Float_t RMargin, Float_t TMargin, Float_t BMargin)
{
  canvas->SetFillColor(0);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetLeftMargin(LMargin);
  canvas->SetRightMargin(RMargin);
  canvas->SetTopMargin(TMargin);
  canvas->SetBottomMargin(BMargin);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
}

void DrawHorLine(Double_t x, Double_t y){
  TLine *line = new TLine(0, y, x, y);
  line->SetLineStyle(2); // Set the line style to dashed (2)
  line->SetLineColor(kBlack); // Set the line color (kRed is a ROOT predefined color)
  line->Draw("same"); // Draw the line on the same canvas
}

void DrawVertLine(Double_t x, Double_t yMin, Double_t yMax, Color_t color){
  TLine *line = new TLine(x, yMin, x, yMax);
  line->SetLineStyle(2); // Set the line style to dashed (2)
  line->SetLineWidth(2); // Set the line width
  line->SetLineColor(color); // Set the line color (kRed is a ROOT predefined color)
  line->Draw("same"); // Draw the line on the same canvas
}

const Int_t numParticles = 7;

TString invMassNames[numParticles] = {
  "InvMassK0S", 
  "InvMassLambda", 
  "InvMassAntiLambda",
  "InvMassXiMinus",
  "InvMassXiPlus",
  "InvMassOmegaMinus",
  "InvMassOmegaPlus",
};

Float_t pdgMass[numParticles] = {0.497611, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245};

Int_t color[8] = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};

TString particleNames[] = {"K0S","Lambda","AntiLambda","XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus"};
TString particleSymnbols[] = {"K^{0}_{S}", "#Lambda", "#bar{#Lambda}","#Xi^{-}", "#Xi^{+}", "#Omega^{-}", "#Omega^{+}"};
