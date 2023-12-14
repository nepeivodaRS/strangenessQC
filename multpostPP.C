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

void StylePad(TPad *pad, Float_t LMargin, Float_t RMargin, Float_t TMargin, Float_t BMargin)
{
  pad->SetFillColor(0);
  pad->SetTickx(1);
  pad->SetTicky(1);
  pad->SetLeftMargin(LMargin);
  pad->SetRightMargin(RMargin);
  pad->SetTopMargin(TMargin);
  pad->SetBottomMargin(BMargin);
}

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title, 
                Bool_t XRange, Float_t XLow, Float_t XUp, Float_t xOffset, Float_t yOffset, Float_t mSize, 
                Float_t xTitleSize, Float_t yTitleSize)
{
  histo->GetYaxis()->SetRangeUser(Low, Up);
  if (XRange)
    histo->GetXaxis()->SetRangeUser(XLow, XUp);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(mSize);
  histo->GetXaxis()->SetTitle(titleX);
  //histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleSize(xTitleSize);
  histo->GetXaxis()->SetTitleOffset(xOffset);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(yTitleSize);
  //histo->GetYaxis()->SetLabelSize(0.05);
  histo->GetYaxis()->SetTitleOffset(yOffset);
  histo->SetTitle(title);
}

void DrawHorLine(Double_t x, Double_t y){
  TLine *line = new TLine(0, y, x, y);
  line->SetLineStyle(2); // Set the line style to dashed (2)
  line->SetLineColor(kBlack); // Set the line color (kRed is a ROOT predefined color)
  line->Draw("same"); // Draw the line on the same canvas
}

void multpostPP(TString fileList = "postPPresults/listPP.txt",
                TString PathOut = "multPP.root")
{
  gROOT->SetBatch(kTRUE);
  // Files with histograms
  std::vector<std::string> name;
  std::vector<std::string> nameLegend;
  std::ifstream file(Form("%s", fileList.Data()));

  std::string remove = "/Users/rnepeiv/workLund/PhD_work/run3QCPbPb/qcTaskDev/postPPscripts/postPPresults/";
  std::string remove2 = ".root";
  cout << "List:" << endl;
  cout << fileList.Data() << endl;
  cout << "Files:" << endl;

  if (file.is_open())
  {
    std::string line;
    while (std::getline(file, line))
    {
      name.push_back(line);
      cout << line << endl;
      size_t pos = line.find(remove);
      if (pos != std::string::npos)
      {
        line.erase(pos, remove.length());
      }
      size_t pos2 = line.find(remove2);
      if (pos2 != std::string::npos)
      {
        line.erase(pos2, remove2.length());
      }
      nameLegend.push_back(line);
    }
    file.close();
  }
  else
  {
    std::cerr << "Unable to open fileList!" << std::endl;
  }

  const Int_t numFiles = name.size();
  cout << "Number of files: " << numFiles << endl;

  TFile *fileIn[numFiles];

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

  TH1F *histoRatioMean[numParticles];
  TH1F *histoRatioWidth[numParticles];
  TH1F *histoRatioYield[numParticles];

  TDirectory *fitDirectories[numFiles];

  TH1F *hMeans[numFiles][numParticles];
  TH1F *hSigmas[numFiles][numParticles];
  TH1F *hYields[numFiles][numParticles];

  TH1F *hYieldsSpecies[3];
  TCanvas *canvasYieldSpecies[3];
  TPad *padYieldSpeciesLow[3];
  TPad *padYieldSpeciesUp[3];
  TH1F *hYieldsRatioSpecies[numFiles][3];
  TH1F *hYieldsDenomSpecies[numFiles][3];

  TH1F *hMeansRatio[numFiles][numParticles];
  TH1F *hSigmasRatio[numFiles][numParticles];
  TH1F *hYieldsRatio[numFiles][numParticles];

  TH1F *hMeansDenom[numParticles];
  TH1F *hSigmasDenom[numParticles];
  TH1F *hYieldsDenom[numParticles];

  TCanvas *canvasMean[numParticles];
  TPad *padMeanLow[numParticles];
  TPad *padMeanUp[numParticles];

  TCanvas *canvasSigma[numParticles];
  TPad *padSigmaLow[numParticles];
  TPad *padSigmaUp[numParticles];

  TCanvas *canvasYield[numParticles];
  TPad *padYieldLow[numParticles];
  TPad *padYieldUp[numParticles];

  Float_t pdgMass[numParticles] = {0.497611, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245};
  Float_t pdgMassError[numParticles] = {0.013, 0.006, 0.006, 0.07, 0.07, 0.29, 0.29};

  Float_t meanYLow[numParticles] = {0.493, 1.114, 1.114, 1.315, 1.315, 1.668, 1.668};
  Float_t meanYUp[numParticles] = {0.505, 1.117, 1.117, 1.328, 1.328, 1.677, 1.677};
  Float_t meanRatioLow[numParticles] = {0.995, 0.998, 0.998, 0.996, 0.996, 0.996, 0.996};
  Float_t meanRatioUp[numParticles] = {1.005, 1.002, 1.002, 1.004, 1.004, 1.004, 1.004};

  Float_t sigmaYLow[numParticles] = {0.003, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005};
  Float_t sigmaYUp[numParticles] = {0.015, 0.008, 0.008, 0.004, 0.004, 0.01, 0.01};
  Float_t sigmaRatioLow[numParticles] = {0.7, 0.5, 0.5, 0.8, 0.8, 0.8, 0.8};
  Float_t sigmaRatioUp[numParticles] = {1.3, 1.5, 1.5, 1.2, 1.2, 1.2, 1.2};

  Float_t yieldYLow[numParticles] = {1e-3, 1e-3, 1e-3, 1e-4, 1e-4, 1e-5, 1e-5};
  Float_t yieldYUp[numParticles] = {100, 10, 10, 0.1, 0.1, 0.01, 0.01};
  Float_t yieldRatioLow[numParticles] = {0.8, 0.8, 0.8, 0.8, 0.8, 0.2, 0.2};
  Float_t yieldRatioUp[numParticles] = {1.5, 1.5, 1.5, 1.2, 1.2, 2, 2};

  Int_t color[8] = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};

  TString particleNames[] = {"K0S","Lambda","AntiLambda","XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus"};
  TString particleSymnbols[] = {"K^{0}_{S}", "#Lambda", "#bar{#Lambda}","#Xi^{-}", "#Xi^{+}", "#Omega^{-}", "#Omega^{+}"};

  for (Int_t iFile = 0; iFile < numFiles; iFile++) {
    // Open input file
    fileIn[iFile] = TFile::Open(Form("%s", name[iFile].c_str()));
    if (!fileIn[iFile] || fileIn[iFile]->IsZombie()) {
        std::cerr << "Error opening input file!" << std::endl;
        return;
    }
    // Open `fitParams` directory
    fitDirectories[iFile] = fileIn[iFile]->GetDirectory("fitParams");
    if (!fitDirectories[iFile])
    {
      std::cerr << "Directory `fitParams` is not found!" << std::endl;
      return;
    }

    // Iterate over particles
    for (Int_t iPart = 0; iPart < numParticles; iPart++) {
      if (iFile == 0) {
        hMeansDenom[iPart] = (TH1F *)fitDirectories[iFile]->Get("Mean_" + invMassNames[iPart]);
        hMeansDenom[iPart]->SetName(Form("meanDenomHisto_%i", iPart));

        hSigmasDenom[iPart] = (TH1F *)fitDirectories[iFile]->Get("Sigma_" + invMassNames[iPart]);
        hSigmasDenom[iPart]->SetName(Form("sigmaDenomHisto_%i", iPart));

        hYieldsDenom[iPart] = (TH1F *)fitDirectories[iFile]->Get("Yield_" + invMassNames[iPart]);
        hYieldsDenom[iPart]->SetName(Form("yieldDenomHisto_%i", iPart));
      }
      // Get mean histogram
      hMeans[iFile][iPart] = (TH1F *)fitDirectories[iFile]->Get("Mean_" + invMassNames[iPart]);
      hMeansRatio[iFile][iPart] = (TH1F *)hMeans[iFile][iPart]->Clone("MeanClone_" + invMassNames[iPart]);
      hMeansRatio[iFile][iPart]->Divide(hMeansDenom[iPart]);

      // Get sigma histogram
      hSigmas[iFile][iPart] = (TH1F *)fitDirectories[iFile]->Get("Sigma_" + invMassNames[iPart]);
      hSigmasRatio[iFile][iPart] = (TH1F *)hSigmas[iFile][iPart]->Clone("SigmaClone_" + invMassNames[iPart]);
      hSigmasRatio[iFile][iPart]->Divide(hSigmasDenom[iPart]);

      // Get yield histogram
      hYields[iFile][iPart] = (TH1F *)fitDirectories[iFile]->Get("Yield_" + invMassNames[iPart]);
      hYieldsRatio[iFile][iPart] = (TH1F *)hYields[iFile][iPart]->Clone("YieldClone_" + invMassNames[iPart]);
      hYieldsRatio[iFile][iPart]->Divide(hYieldsDenom[iPart]);

      if (iPart == 1){
        hYieldsDenomSpecies[iFile][0] = (TH1F *)hYields[iFile][iPart]->Clone(Form("DenomComaprisonHistoK0s_%s", name[iFile].c_str()));
      }
      if (iPart == 2) {
        hYieldsRatioSpecies[iFile][0] = (TH1F *)hYields[iFile][iPart]->Clone("YieldCompClone_" + invMassNames[iPart]);
        hYieldsRatioSpecies[iFile][0]->Divide(hYieldsDenomSpecies[iFile][0]); 
      }
    
    }
  }

  gStyle->SetOptStat(0);
  for (Int_t iPart = 0; iPart < numParticles; iPart++) {
    // Mean
    canvasMean[iPart] = new TCanvas("mean_" + particleNames[iPart], particleNames[iPart], 800, 600);
    StyleCanvas(canvasMean[iPart], 0.15, 0.05, 0.05, 0.15);
    padMeanUp[iPart] = new TPad("pad1" + particleNames[iPart], "pad1" + particleNames[iPart], 0, 0.36, 1, 1);
    padMeanLow[iPart] = new TPad("pad2" + particleNames[iPart], "pad2" + particleNames[iPart], 0, 0.01, 1, 0.35);
    StylePad(padMeanUp[iPart], 0.15, 0.05, 0.05, 0.01);
    StylePad(padMeanLow[iPart], 0.15, 0.05, 0.03, 0.2);
    TLegend *legMean = new TLegend(0.65, 0.75, 0.85, 0.9);
    legMean->SetBorderSize(0);
    legMean->SetFillStyle(0);
    canvasMean[iPart]->cd();
    padMeanUp[iPart]->Draw();
    padMeanLow[iPart]->Draw();

    // Sigma
    canvasSigma[iPart] = new TCanvas("sigma_" + particleNames[iPart], particleNames[iPart], 800, 600);
    StyleCanvas(canvasSigma[iPart], 0.15, 0.05, 0.05, 0.15);
    padSigmaUp[iPart] = new TPad("pad1" + particleNames[iPart], "pad1" + particleNames[iPart], 0, 0.36, 1, 1);
    padSigmaLow[iPart] = new TPad("pad2" + particleNames[iPart], "pad2" + particleNames[iPart], 0, 0.01, 1, 0.35);
    StylePad(padSigmaUp[iPart], 0.15, 0.05, 0.05, 0.01);
    StylePad(padSigmaLow[iPart], 0.15, 0.05, 0.03, 0.2);
    TLegend *legSigma = new TLegend(0.65, 0.75, 0.85, 0.9);
    legSigma->SetBorderSize(0);
    legSigma->SetFillStyle(0);
    canvasSigma[iPart]->cd();
    padSigmaUp[iPart]->Draw();
    padSigmaLow[iPart]->Draw();

    // Yield
    canvasYield[iPart] = new TCanvas("yield_" + particleNames[iPart], particleNames[iPart], 800, 600);
    StyleCanvas(canvasYield[iPart], 0.15, 0.05, 0.05, 0.15);
    padYieldUp[iPart] = new TPad("pad1" + particleNames[iPart], "pad1" + particleNames[iPart], 0, 0.36, 1, 1);
    padYieldLow[iPart] = new TPad("pad2" + particleNames[iPart], "pad2" + particleNames[iPart], 0, 0.01, 1, 0.35);
    StylePad(padYieldUp[iPart], 0.15, 0.05, 0.05, 0.01);
    StylePad(padYieldLow[iPart], 0.15, 0.05, 0.03, 0.2);
    TLegend *legYield = new TLegend(0.65, 0.75, 0.85, 0.9);
    legYield->SetBorderSize(0);
    legYield->SetFillStyle(0);
    canvasYield[iPart]->cd();
    padYieldUp[iPart]->Draw();
    padYieldLow[iPart]->Draw();    

    for (Int_t iFile = 0; iFile < numFiles; iFile++) {
      // Mean
      // Up
      padMeanUp[iPart]->cd();
      hMeans[iFile][iPart]->GetYaxis()->SetMaxDigits(4);
      hMeans[iFile][iPart]->GetYaxis()->SetDecimals(kTRUE);
      hMeans[iFile][iPart]->GetXaxis()->SetTitle("");
      hMeans[iFile][iPart]->SetLineColor(color[iFile]);
      hMeans[iFile][iPart]->GetXaxis()->SetLabelSize(0.);
      StyleHisto(hMeans[iFile][iPart], meanYLow[iPart], meanYUp[iPart], color[iFile], 20, "", hMeans[iFile][iPart]->GetYaxis()->GetTitle(), "", 0, 0, 0, 1.5, 1.0, 1, 0.0, 0.05);
      TAxis *axisMean = hMeans[iFile][iPart]->GetYaxis();
      axisMean->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
      hMeans[iFile][iPart]->Draw("same");
      legMean->AddEntry(hMeans[iFile][iPart], nameLegend[iFile].c_str(), "pl");
      // Low
      padMeanLow[iPart]->cd();
      hMeansRatio[iFile][iPart]->GetYaxis()->SetTitle("");
      hMeansRatio[iFile][iPart]->GetXaxis()->SetLabelSize(0.08);
      hMeansRatio[iFile][iPart]->GetXaxis()->SetTitleOffset(1.2);
      hMeansRatio[iFile][iPart]->GetYaxis()->SetTitleOffset(0.7);
      hMeansRatio[iFile][iPart]->GetYaxis()->SetLabelSize(0.08);
      StyleHisto(hMeansRatio[iFile][iPart], meanRatioLow[iPart], meanRatioUp[iPart], color[iFile], 20, "#it{p}_{T} (GeV/#it{c})", Form("Ratio to %s", nameLegend[0].c_str()), "", 0, 0, 0, 1.0, 0.6, 1, 0.08, 0.08);
      if(!(iFile == 0)) {
        hMeansRatio[iFile][iPart]->Draw("same");
      }

      // Sigma
      // Up
      padSigmaUp[iPart]->cd();
      hSigmas[iFile][iPart]->GetYaxis()->SetMaxDigits(4);
      hSigmas[iFile][iPart]->GetYaxis()->SetDecimals(kTRUE);
      hSigmas[iFile][iPart]->GetXaxis()->SetTitle("");
      hSigmas[iFile][iPart]->SetLineColor(color[iFile]);
      hSigmas[iFile][iPart]->GetXaxis()->SetLabelSize(0.);
      StyleHisto(hSigmas[iFile][iPart], sigmaYLow[iPart], sigmaYUp[iPart], color[iFile], 20, "", hSigmas[iFile][iPart]->GetYaxis()->GetTitle(), "", 0, 0, 0, 1.5, 1.0, 1, 0.0, 0.05);
      TAxis *axisSigma = hSigmas[iFile][iPart]->GetYaxis();
      axisSigma->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
      hSigmas[iFile][iPart]->Draw("same");
      legSigma->AddEntry(hSigmas[iFile][iPart], nameLegend[iFile].c_str(), "pl");
      // Low
      padSigmaLow[iPart]->cd();
      hSigmasRatio[iFile][iPart]->GetYaxis()->SetTitle("");
      hSigmasRatio[iFile][iPart]->GetXaxis()->SetLabelSize(0.08);
      hSigmasRatio[iFile][iPart]->GetXaxis()->SetTitleOffset(1.2);
      hSigmasRatio[iFile][iPart]->GetYaxis()->SetTitleOffset(0.7);
      hSigmasRatio[iFile][iPart]->GetYaxis()->SetLabelSize(0.08);
      StyleHisto(hSigmasRatio[iFile][iPart], sigmaRatioLow[iPart], sigmaRatioUp[iPart], color[iFile], 20, "#it{p}_{T} (GeV/#it{c})", Form("Ratio to %s", nameLegend[0].c_str()), "", 0, 0, 0, 1.0, 0.6, 1, 0.08, 0.08);
      if(!(iFile == 0)) {
        hSigmasRatio[iFile][iPart]->Draw("same");
      }

      // Yield
      // Up
      padYieldUp[iPart]->cd();
      hYields[iFile][iPart]->GetYaxis()->SetMaxDigits(4);
      hYields[iFile][iPart]->GetYaxis()->SetDecimals(kTRUE);
      hYields[iFile][iPart]->GetXaxis()->SetTitle("");
      hYields[iFile][iPart]->SetLineColor(color[iFile]);
      hYields[iFile][iPart]->GetXaxis()->SetLabelSize(0.);
      StyleHisto(hYields[iFile][iPart], yieldYLow[iPart], yieldYUp[iPart], color[iFile], 20, "", hYields[iFile][iPart]->GetYaxis()->GetTitle(), "", 0, 0, 0, 1.5, 1.0, 1, 0.0, 0.05);
      padYieldUp[iPart]->SetLogy();
      TAxis *axisYield = hYields[iFile][iPart]->GetYaxis();
      axisYield->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
      hYields[iFile][iPart]->Draw("same");
      legYield->AddEntry(hYields[iFile][iPart], nameLegend[iFile].c_str(), "pl");
      // Low
      padYieldLow[iPart]->cd();
      hYieldsRatio[iFile][iPart]->GetYaxis()->SetTitle("");
      hYieldsRatio[iFile][iPart]->GetXaxis()->SetLabelSize(0.08);
      hYieldsRatio[iFile][iPart]->GetXaxis()->SetTitleOffset(1.2);
      hYieldsRatio[iFile][iPart]->GetYaxis()->SetTitleOffset(0.7);
      hYieldsRatio[iFile][iPart]->GetYaxis()->SetLabelSize(0.08);
      StyleHisto(hYieldsRatio[iFile][iPart], yieldRatioLow[iPart], yieldRatioUp[iPart], color[iFile], 20, "#it{p}_{T} (GeV/#it{c})", Form("Ratio to %s", nameLegend[0].c_str()), "", 0, 0, 0, 1.0, 0.6, 1, 0.08, 0.08);
      if(!(iFile == 0)) {
        hYieldsRatio[iFile][iPart]->Draw("same");
      }
    }

    TLegend *LegendTitle = new TLegend(0.25, 0.7, 0.55, 0.9);
    LegendTitle->SetFillStyle(0);
    LegendTitle->SetTextAlign(33);
    LegendTitle->SetTextSize(0.06);
    LegendTitle->SetTextFont(42);
    LegendTitle->SetLineColorAlpha(0.,0.);
    LegendTitle->SetFillColorAlpha(0.,0.);
    LegendTitle->SetBorderSize(0.);
    LegendTitle->AddEntry("", "#bf{ALICE Work In Progress}", "");
    LegendTitle->AddEntry("", "Pb--Pb, #sqrt{#it{s}} = 5.36 TeV", "");
    LegendTitle->AddEntry("", particleSymnbols[iPart], "");

    // Mean
    padMeanUp[iPart]->cd();
    LegendTitle->Draw();
    legMean->Draw();
    DrawHorLine(6.0, pdgMass[iPart]);
    padMeanLow[iPart]->cd();
    DrawHorLine(6.0, 1.0);
    gPad->Update();
    canvasMean[iPart]->Update();
    canvasMean[iPart]->SaveAs("multResults/pdf/" + particleNames[iPart] + ".pdf(");

    // Sigma
    padSigmaUp[iPart]->cd();
    LegendTitle->Draw();
    legSigma->Draw();
    padSigmaLow[iPart]->cd();
    DrawHorLine(6.0, 1.0);
    gPad->Update();
    canvasSigma[iPart]->Update();
    canvasSigma[iPart]->SaveAs("multResults/pdf/" + particleNames[iPart] + ".pdf");

    // Yield
    padYieldUp[iPart]->cd();
    LegendTitle->Draw();
    legYield->Draw();
    padYieldLow[iPart]->cd();
    DrawHorLine(6.0, 1.0);
    gPad->Update();
    canvasYield[iPart]->Update();
    canvasYield[iPart]->SaveAs("multResults/pdf/" + particleNames[iPart] + ".pdf");

    delete canvasMean[iPart];
    delete canvasSigma[iPart];
    delete canvasYield[iPart];
  }

  // Yield
  canvasYieldSpecies[0] = new TCanvas("yieldComaprison_" + particleNames[1], particleNames[1], 800, 600);
  StyleCanvas(canvasYieldSpecies[0], 0.15, 0.05, 0.05, 0.15);
  padYieldSpeciesUp[0] = new TPad("pad1Comaprison" + particleNames[1], "pad1Comaprison" + particleNames[1], 0, 0.36, 1, 1);
  padYieldSpeciesLow[0] = new TPad("pad2Comparison" + particleNames[1], "pad2Comparison" + particleNames[1], 0, 0.01, 1, 0.35);
  StylePad(padYieldSpeciesUp[0], 0.15, 0.05, 0.05, 0.01);
  StylePad(padYieldSpeciesLow[0], 0.15, 0.05, 0.03, 0.2);
  TLegend *legYieldRato = new TLegend(0.65, 0.75, 0.85, 0.9);
  legYieldRato->SetBorderSize(0);
  legYieldRato->SetFillStyle(0);
  canvasYieldSpecies[0]->cd();
  padYieldSpeciesUp[0]->Draw();
  padYieldSpeciesLow[0]->Draw();

  for (Int_t iFile = 0; iFile < numFiles; iFile++) {
    // Up
    padYieldSpeciesUp[0]->cd();
    hYields[iFile][1]->GetYaxis()->SetMaxDigits(4);
    hYields[iFile][1]->GetYaxis()->SetDecimals(kTRUE);
    hYields[iFile][1]->GetXaxis()->SetTitle("");
    hYields[iFile][1]->SetLineColor(color[iFile]);
    hYields[iFile][1]->GetXaxis()->SetLabelSize(0.);
    StyleHisto(hYields[iFile][1], yieldYLow[1], yieldYUp[1], color[iFile], 20, "", hYields[iFile][1]->GetYaxis()->GetTitle(), "", 0, 0, 0, 1.5, 1.0, 1, 0.0, 0.05);
    padYieldSpeciesUp[0]->SetLogy();
    TAxis *axisYield = hYields[iFile][1]->GetYaxis();
    axisYield->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
    hYields[iFile][1]->Draw("same");
    legYieldRato->AddEntry(hYields[iFile][1], nameLegend[iFile].c_str(), "pl");

    hYields[iFile][2]->GetYaxis()->SetMaxDigits(4);
    hYields[iFile][2]->GetYaxis()->SetDecimals(kTRUE);
    hYields[iFile][2]->GetXaxis()->SetTitle("");
    hYields[iFile][2]->SetLineColor(color[iFile]);
    hYields[iFile][2]->SetLineStyle(2);
    hYields[iFile][2]->GetXaxis()->SetLabelSize(0.);
    StyleHisto(hYields[iFile][2], yieldYLow[2], yieldYUp[2], color[iFile], 20, "", hYields[iFile][2]->GetYaxis()->GetTitle(), "", 0, 0, 0, 1.5, 1.0, 1, 0.0, 0.05);
    padYieldSpeciesUp[0]->SetLogy();
    hYields[iFile][2]->Draw("same");

    // Low
    padYieldSpeciesLow[0]->cd();
    hYieldsRatioSpecies[iFile][0]->GetYaxis()->SetTitle("");
    hYieldsRatioSpecies[iFile][0]->GetXaxis()->SetLabelSize(0.08);
    hYieldsRatioSpecies[iFile][0]->GetXaxis()->SetTitleOffset(1.2);
    hYieldsRatioSpecies[iFile][0]->GetYaxis()->SetTitleOffset(0.7);
    hYieldsRatioSpecies[iFile][0]->GetYaxis()->SetLabelSize(0.08);
    StyleHisto(hYieldsRatioSpecies[iFile][0], 0, 1.0, color[iFile], 20, "#it{p}_{T} (GeV/#it{c})", "Anti-Particle to Particle Ratio", "", 0, 0, 0, 1.0, 0.6, 1, 0.08, 0.08);
    hYieldsRatioSpecies[iFile][0]->Draw("same");
  }

  TLegend *LegendTitle = new TLegend(0.25, 0.7, 0.55, 0.9);
  LegendTitle->SetFillStyle(0);
  LegendTitle->SetTextAlign(33);
  LegendTitle->SetTextSize(0.06);
  LegendTitle->SetTextFont(42);
  LegendTitle->SetLineColorAlpha(0.,0.);
  LegendTitle->SetFillColorAlpha(0.,0.);
  LegendTitle->SetBorderSize(0.);
  LegendTitle->AddEntry("", "#bf{ALICE Work In Progress}", "");
  LegendTitle->AddEntry("", "Pb--Pb, #sqrt{#it{s}} = 5.36 TeV", "");
  LegendTitle->AddEntry("", particleSymnbols[1] + " and " + particleSymnbols[2], "");
  // Yield
  padYieldSpeciesUp[0]->cd();
  legYieldRato->Draw();
  LegendTitle->Draw();
  padYieldSpeciesLow[0]->cd();
  gPad->Update();
  canvasYieldSpecies[0]->Update();
  canvasYieldSpecies[0]->SaveAs("multResults/pdf/" + particleNames[1] + ".pdf)");


  for (Int_t iFile = 0; iFile < numFiles; iFile++) {
    fileIn[iFile]->Close();
  }

  gROOT->SetBatch(kFALSE);
}