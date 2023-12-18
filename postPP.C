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

Bool_t reject = kTRUE;
Double_t fparab(Double_t *x, Double_t *par)
{
  const Int_t numPart = 7;
  // Signal region
  Float_t liminf[numPart] = {0.460, 1.105, 1.105, 1.315, 1.315, 1.668, 1.668};
  Float_t limsup[numPart] = {0.522, 1.125, 1.125, 1.328, 1.328, 1.677, 1.677};
  Int_t part = par[3];
  if (reject && x[0] > liminf[part] && x[0] < limsup[part])
  {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0];
}

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

void DrawVertLine(Double_t x, Double_t yMin, Double_t yMax, Color_t color){
  TLine *line = new TLine(x, yMin, x, yMax);
  line->SetLineStyle(2); // Set the line style to dashed (2)
  line->SetLineWidth(2); // Set the line width
  line->SetLineColor(color); // Set the line color (kRed is a ROOT predefined color)
  line->Draw("same"); // Draw the line on the same canvas
}

void postPP(TString fileList = "listQC.txt", // PP and QC task
            TString PathOut = "postPP.root")
{
  gROOT->SetBatch(kTRUE);
  // Files with histograms
  std::vector<std::string> name;
  std::vector<std::string> nameLegend;
  std::ifstream file(Form("%s", fileList.Data()));

  std::string remove = "/Users/rnepeiv/workLund/PhD_work/run3QCPbPb/qcTaskDev/results/combined/AnalysisResults_";
  std::string remove2 = ".root";

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

  cout << "List:" << endl;
  cout << fileList.Data() << endl;
  cout << "Files:" << endl;

  const Int_t numFiles = name.size();
  cout << "Number of files: " << numFiles << endl;

  TFile *fileIn[numFiles];
  TFile *fileOut[numFiles];

  const Int_t numParticles = 7;

  TDirectory *partDirectories[numFiles][numParticles];

  TString invMassNames[numParticles] = {
    "InvMassK0S", 
    "InvMassLambda", 
    "InvMassAntiLambda",
    "InvMassXiMinus",
    "InvMassXiPlus",
    "InvMassOmegaMinus",
    "InvMassOmegaPlus",
  };

  TString dirNames[numParticles] = {
    "lf-strangenessqcpp/k0S", 
    "lf-strangenessqcpp/lambda", 
    "lf-strangenessqcpp/antiLambda",
    "lf-strangenessqcpp/xi", 
    "lf-strangenessqcpp/antixi", 
    "lf-strangenessqcpp/omega",
    "lf-strangenessqcpp/antiomega"
  };

  TH3F *hInvMass3D[numFiles][numParticles];
  TH2F *hInvMass2D[numFiles][numParticles];
  TH2F *hInvMass2Dclones[numFiles][numParticles];
  TH1F *hInvMass1D[numFiles][numParticles];

  Int_t const numPtBins = 6; // maximum number of bins for a particle (in the following array)
  Int_t const numPtBinsPart[numParticles] = {
    6,
    4,
    4,
    1,
    1,
    1,
    1
  };
  Float_t ptBins[numParticles][numPtBins + 1] = {
    {0., 0.7, 1.0, 1.5, 2., 3., 6.},
    {0., 1., 1.5, 2., 6.},
    {0., 1., 1.5, 2., 6.},
    {0., 6.},
    {0., 6.},
    {0., 6.},
    {0., 6.}
  };

  TH2F *hInvMass2Dpt[numFiles][numParticles][numPtBins];
  TH1F *hInvMass1Dpt[numFiles][numParticles][numPtBins];

  Float_t minRangeSignal[numParticles] = {0.47, 1.112, 1.112, 1.315, 1.315, 1.668, 1.668};
  Float_t maxRangeSignal[numParticles] = {0.52, 1.122, 1.120, 1.328, 1.328, 1.677, 1.677};
  Float_t minRange[numParticles] =       {0.435, 1.095, 1.10, 1.305, 1.305, 1.655, 1.655};
  Float_t maxRange[numParticles] =       {0.560, 1.140, 1.14, 1.340, 1.340, 1.690, 1.690};

  Float_t pdgMass[numParticles] = {0.497611, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245};

  Float_t meanYLow[numParticles] = {0.493, 1.114, 1.114, 1.315, 1.315, 1.668, 1.668};
  Float_t meanYUp[numParticles] = {0.505, 1.117, 1.117, 1.328, 1.328, 1.677, 1.677};

  Float_t sigmaYLow[numParticles] = {0.003, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005};
  Float_t sigmaYUp[numParticles] = {0.015, 0.008, 0.008, 0.004, 0.004, 0.01, 0.01};

  Float_t yieldYLow[numParticles] = {1e-3, 1e-3, 1e-3, 1e-4, 1e-4, 1e-5, 1e-5};
  Float_t yieldYUp[numParticles] = {100, 10, 10, 0.1, 0.1, 0.01, 0.01};

  TH1F *hMeans[numFiles][numParticles];
  TH1F *hSigmas[numFiles][numParticles];
  TH1F *hYields[numFiles][numParticles];

  Float_t mean[numFiles][numParticles][numPtBins];
  Float_t errMean[numFiles][numParticles][numPtBins];

  Float_t sigma[numFiles][numParticles][numPtBins];
  Float_t errSigma[numFiles][numParticles][numPtBins];

  Float_t yieldWithBG[numFiles][numParticles][numPtBins];
  Float_t yieldBG[numFiles][numParticles][numPtBins];
  Float_t errYieldBG[numFiles][numParticles][numPtBins];

  TF1 *total[numFiles][numParticles][numPtBins];
  TF1 *bkgparab[numFiles][numParticles][numPtBins]; // bg function without signal region
  TF1 *bkgparabDraw[numFiles][numParticles][numPtBins]; // bg function in full range

  // fFitResultTotal is not used
  TFitResultPtr fFitResultTotal[numFiles][numParticles][numPtBins];
  // fFitResultParab is used for bg integral error
  TFitResultPtr fFitResultParab[numFiles][numParticles][numPtBins];

  TCanvas *canvasInvMass[numFiles][numParticles][numPtBins];
  cout << ptBins[0][0] << " " << ptBins[0][1] << endl;

  for (Int_t iFile = 0; iFile < numFiles; iFile++)
  {
    // Open input file
    fileIn[iFile] = TFile::Open(Form("%s", name[iFile].c_str()));
    if (!fileIn[iFile] || fileIn[iFile]->IsZombie()) {
        std::cerr << "Error opening input file!" << std::endl;
        return;
    }

    TDirectory* eventSeldir = fileIn[iFile]->GetDirectory("lf-strangenessqc/eventSelection");
    if (!eventSeldir)
    {
      std::cerr << "`lf-strangenessqc/eventSelection` directory is not found!" << std::endl;
      return;
    }

    TH1F* hNevents = (TH1F *)eventSeldir->Get("hVertexZRec");
    if (!hNevents)
    {
      std::cerr << "Histogram `hVertexZRec` is not found!" << std::endl;
      return;
    }
    hNevents->Draw();
    gPad->SaveAs(Form("postPPresults/%s.pdf[", nameLegend[iFile].c_str()));
    Float_t norm = hNevents->GetEntries();

    fileIn[iFile]->cd();
    // Create output file
    fileOut[iFile] = new TFile(Form("postPPresults/%s.root", nameLegend[iFile].c_str()), "recreate");

    // Create dir to store all the fitting parameters
    TDirectory* fitParamsOutDir = fileOut[iFile]->mkdir("fitParams");
    // Create dir to store all the invariant mass histograms
    TDirectory* invMassOutDir = fileOut[iFile]->mkdir("invMassHists");
    // Iterate over particles
    for (Int_t iPart = 0; iPart < numParticles; iPart++) {
      // Open `iPart` directory
      partDirectories[iFile][iPart] = fileIn[iFile]->GetDirectory(dirNames[iPart]);
      if (!partDirectories[iFile][iPart])
      {
        std::cerr << "Directory " << dirNames[iPart] << " not found!" << std::endl;
        return;
      }
      // Get inv mass histogram
      hInvMass3D[iFile][iPart] = (TH3F *)partDirectories[iFile][iPart]->Get(invMassNames[iPart]);
      if (!hInvMass3D[iFile][iPart])
      {
        std::cerr << "Histogram " << invMassNames[iPart] << " not found!" << std::endl;
        return;
      }

      // Make it 2D
      hInvMass2D[iFile][iPart] = static_cast<TH2F*>(hInvMass3D[iFile][iPart]->Project3D("xy"));
      hInvMass2Dclones[iFile][iPart] = (TH2F*)hInvMass2D[iFile][iPart]->Clone(Form("2DHist_%d_%d", iFile, iPart));
      hInvMass1D[iFile][iPart] = (TH1F *)hInvMass2Dclones[iFile][iPart]->ProjectionX();
      hInvMass1D[iFile][iPart]->Sumw2();
      hInvMass1D[iFile][iPart]->SetName("fullInvMass_" + invMassNames[iPart]);

      // Create mean, sigma and yield hsitograms
      hMeans[iFile][iPart] = new TH1F("Mean_" + invMassNames[iPart], "Mean_" + invMassNames[iPart], numPtBinsPart[iPart], ptBins[iPart]);
      hSigmas[iFile][iPart] = new TH1F("Sigma_" + invMassNames[iPart], "Sigma_" + invMassNames[iPart], numPtBinsPart[iPart], ptBins[iPart]);
      hYields[iFile][iPart] = new TH1F("Yield_" + invMassNames[iPart], "Yield_" + invMassNames[iPart], numPtBinsPart[iPart], ptBins[iPart]);
      for (Int_t iPt = 0; iPt < numPtBinsPart[iPart]; iPt++) {
        hInvMass2Dpt[iFile][iPart][iPt] = (TH2F*)hInvMass2D[iFile][iPart]->Clone(Form("2DHistInPtBin_%d_%d_%d", iFile, iPart, iPt));
        hInvMass2Dpt[iFile][iPart][iPt]->GetYaxis()->SetRangeUser(ptBins[iPart][iPt], ptBins[iPart][iPt+1]);
        hInvMass1Dpt[iFile][iPart][iPt] = (TH1F *)hInvMass2Dpt[iFile][iPart][iPt]->ProjectionX();
        hInvMass1Dpt[iFile][iPart][iPt]->Sumw2();
        hInvMass1Dpt[iFile][iPart][iPt]->GetYaxis()->SetTitle("Counts");
        hInvMass1Dpt[iFile][iPart][iPt]->SetTitle(Form("%.2f-%.2f pt bin", ptBins[iPart][iPt], ptBins[iPart][iPt+1]));
        // Fit inv mass histogram
        canvasInvMass[iFile][iPart][iPt] = new TCanvas(invMassNames[iPart] + Form("_%d", iPt), Form("%s_%d_%d_%d", hInvMass3D[iFile][iPart]->GetName(), iFile, iPart, iPt), 1000, 800);
        StyleCanvas(canvasInvMass[iFile][iPart][iPt], 0.14, 0.05, 0.11, 0.15);

        canvasInvMass[iFile][iPart][iPt]->cd();
        bkgparab[iFile][iPart][iPt] = new TF1(Form("parab_%d_%d_%d", iFile, iPart, iPt), fparab, minRange[iPart], maxRange[iPart], 4);
        bkgparab[iFile][iPart][iPt]->SetLineColor(kGreen);
        bkgparab[iFile][iPart][iPt]->FixParameter(3, iPart); // Fix Particle type to reject signal range
        fFitResultParab[iFile][iPart][iPt] = hInvMass1Dpt[iFile][iPart][iPt]->Fit(bkgparab[iFile][iPart][iPt], "R0QS");
        Double_t parBG[4];
        bkgparab[iFile][iPart][iPt]->GetParameters(parBG);
        total[iFile][iPart][iPt] = new TF1(Form("total_%d_%d", iFile, iPart), "gaus(0) + gaus(3) + pol2(6)", minRange[iPart], maxRange[iPart]);
        total[iFile][iPart][iPt]->SetLineColor(kRed);
        total[iFile][iPart][iPt]->SetParName(0, "norm");
        total[iFile][iPart][iPt]->SetParName(1, "mean");
        total[iFile][iPart][iPt]->SetParName(2, "sigma");
        total[iFile][iPart][iPt]->SetParName(3, "norm2");
        total[iFile][iPart][iPt]->SetParName(4, "mean2");
        total[iFile][iPart][iPt]->SetParName(5, "sigma2");
        total[iFile][iPart][iPt]->SetParName(6, "p0");
        total[iFile][iPart][iPt]->SetParName(7, "p1");
        total[iFile][iPart][iPt]->SetParName(8, "p2");
        Float_t peakValue = hInvMass1Dpt[iFile][iPart][iPt]->GetBinContent(hInvMass1Dpt[iFile][iPart][iPt]->GetMaximumBin());
        total[iFile][iPart][iPt]->SetParameters(peakValue, pdgMass[iPart], 0.001, peakValue, pdgMass[iPart], 0.001);
        total[iFile][iPart][iPt]->SetParLimits(0, 0., 1.2 * peakValue);
        total[iFile][iPart][iPt]->SetParLimits(1, minRangeSignal[iPart], maxRangeSignal[iPart]);
        total[iFile][iPart][iPt]->SetParLimits(2, 0.0001, 0.05);
        total[iFile][iPart][iPt]->SetParLimits(3, 0., 1.2 * peakValue);
        total[iFile][iPart][iPt]->SetParLimits(4, minRangeSignal[iPart], maxRangeSignal[iPart]);
        total[iFile][iPart][iPt]->SetParLimits(5, 0.0001, 0.03);
        total[iFile][iPart][iPt]->FixParameter(6, parBG[0]);
        total[iFile][iPart][iPt]->FixParameter(7, parBG[1]);
        total[iFile][iPart][iPt]->FixParameter(8, parBG[2]);
        total[iFile][iPart][iPt]->SetNpx(1e4);
        total[iFile][iPart][iPt]->SetLineWidth(3);
        fFitResultTotal[iFile][iPart][iPt] = hInvMass1Dpt[iFile][iPart][iPt]->Fit(total[iFile][iPart][iPt], "SRB+Q");
        bkgparabDraw[iFile][iPart][iPt] = new TF1(Form("parabToDraw_%d_%d_%d", iFile, iPart, iPt), "pol2", minRange[iPart], maxRange[iPart]);
        bkgparabDraw[iFile][iPart][iPt]->FixParameter(0, parBG[0]);
        bkgparabDraw[iFile][iPart][iPt]->FixParameter(1, parBG[1]);
        bkgparabDraw[iFile][iPart][iPt]->FixParameter(2, parBG[2]);
        bkgparabDraw[iFile][iPart][iPt]->SetLineStyle(2);
        bkgparabDraw[iFile][iPart][iPt]->SetNpx(1e4);
        bkgparabDraw[iFile][iPart][iPt]->SetLineWidth(3);
        bkgparabDraw[iFile][iPart][iPt]->SetLineColor(kCyan+2); //{kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
        bkgparabDraw[iFile][iPart][iPt]->Draw("same");

        mean[iFile][iPart][iPt] = (total[iFile][iPart][iPt]->GetParameter(1) + total[iFile][iPart][iPt]->GetParameter(4))/2.;
        errMean[iFile][iPart][iPt] = sqrt(pow(total[iFile][iPart][iPt]->GetParError(1), 2) + pow(total[iFile][iPart][iPt]->GetParError(4), 2));
        sigma[iFile][iPart][iPt] = (total[iFile][iPart][iPt]->GetParameter(2) + total[iFile][iPart][iPt]->GetParameter(5))/2.;
        errSigma[iFile][iPart][iPt] = sqrt(pow(total[iFile][iPart][iPt]->GetParError(2), 2) + pow(total[iFile][iPart][iPt]->GetParError(5), 2));

        Double_t leftSignal = mean[iFile][iPart][iPt] - 5 * sigma[iFile][iPart][iPt];
        Double_t rightSignal = mean[iFile][iPart][iPt] + 5 * sigma[iFile][iPart][iPt];
        Double_t yaxisMin = hInvMass1Dpt[iFile][iPart][iPt]->GetMinimum();
        DrawVertLine(leftSignal, yaxisMin, peakValue, kBlack);
        DrawVertLine(rightSignal, yaxisMin, peakValue, kBlack);
        DrawVertLine(minRange[iPart], yaxisMin, peakValue, kMagenta+1);
        DrawVertLine(maxRange[iPart], yaxisMin, peakValue, kMagenta+1);

        gPad->Update();
        canvasInvMass[iFile][iPart][iPt]->Update();
        fileOut[iFile]->cd();
        invMassOutDir->cd();
        canvasInvMass[iFile][iPart][iPt]->Write();
        canvasInvMass[iFile][iPart][iPt]->SaveAs(Form("postPPresults/%s.pdf", nameLegend[iFile].c_str()));
        delete canvasInvMass[iFile][iPart][iPt];

        hMeans[iFile][iPart]->SetBinContent(iPt + 1, mean[iFile][iPart][iPt]);
        hMeans[iFile][iPart]->SetBinError(iPt + 1, errMean[iFile][iPart][iPt]);

        hSigmas[iFile][iPart]->SetBinContent(iPt + 1, sigma[iFile][iPart][iPt]);
        hSigmas[iFile][iPart]->SetBinError(iPt + 1, errSigma[iFile][iPart][iPt]);

        yieldBG[iFile][iPart][iPt] = bkgparabDraw[iFile][iPart][iPt]->Integral(leftSignal, rightSignal)/hInvMass1Dpt[iFile][iPart][iPt]->GetBinWidth(1);
        errYieldBG[iFile][iPart][iPt] = bkgparabDraw[iFile][iPart][iPt]->IntegralError(leftSignal, rightSignal, fFitResultParab[iFile][iPart][iPt]->GetParams(), (fFitResultParab[iFile][iPart][iPt]->GetCovarianceMatrix()).GetMatrixArray());

        yieldWithBG[iFile][iPart][iPt] = 0;
        for (Int_t bin = hInvMass1Dpt[iFile][iPart][iPt]->GetXaxis()->FindBin(leftSignal); 
             bin <= hInvMass1Dpt[iFile][iPart][iPt]->GetXaxis()->FindBin(rightSignal); 
             bin++)
        {
          yieldWithBG[iFile][iPart][iPt] += hInvMass1Dpt[iFile][iPart][iPt]->GetBinContent(bin);
        }

        hYields[iFile][iPart]->SetBinContent(iPt + 1, (yieldWithBG[iFile][iPart][iPt] - yieldBG[iFile][iPart][iPt]) / hYields[iFile][iPart]->GetBinWidth(iPt + 1) / norm);
        hYields[iFile][iPart]->SetBinError(iPt + 1, sqrt(yieldWithBG[iFile][iPart][iPt] + pow(errYieldBG[iFile][iPart][iPt], 2)) / hYields[iFile][iPart]->GetBinWidth(iPt + 1) / norm);
      }

      fitParamsOutDir->cd();

      StyleHisto(hMeans[iFile][iPart], meanYLow[iPart], meanYUp[iPart], kBlack, 20, "#it{p}_{T} (GeV/#it{c})", "Mean (GeV/#it{c}^{2})", "", 0, 0, 0, 1.0, 1.25, 1, 0.04, 0.04);
      hMeans[iFile][iPart]->Draw();
      gPad->SaveAs(Form("postPPresults/%s.pdf", nameLegend[iFile].c_str()));
      hMeans[iFile][iPart]->Write();

      StyleHisto(hSigmas[iFile][iPart], sigmaYLow[iPart], sigmaYUp[iPart], kBlack, 20, "#it{p}_{T} (GeV/#it{c})", "Sigma (GeV/#it{c}^{2})", "", 0, 0, 0, 1.0, 1.25, 1, 0.04, 0.04);
      hSigmas[iFile][iPart]->Draw();
      gPad->SaveAs(Form("postPPresults/%s.pdf", nameLegend[iFile].c_str()));
      hSigmas[iFile][iPart]->Write();

      StyleHisto(hYields[iFile][iPart], yieldYLow[iPart], yieldYUp[iPart], kBlack, 20, "#it{p}_{T} (GeV/#it{c})", "1/#it{N}_{evt} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}", "", 0, 0, 0, 1.0, 1.1, 1, 0.04, 0.04);
      hYields[iFile][iPart]->Draw();
      gPad->SetLogy();
      if (iPart == (numParticles - 1)) {
        gPad->SaveAs(Form("postPPresults/%s.pdf]", nameLegend[iFile].c_str()));
      } else {
        gPad->SaveAs(Form("postPPresults/%s.pdf", nameLegend[iFile].c_str()));
      }
      hYields[iFile][iPart]->Write();
      gPad->SetLogy(0);

      invMassOutDir->cd();
      hInvMass1D[iFile][iPart]->Write();
    }
    fileIn[iFile]->Close();
    fileOut[iFile]->Close();
  }
  gROOT->SetBatch(kFALSE);
}