#include "postPP.h"

Bool_t reject = kTRUE;
Double_t fparab(Double_t *x, Double_t *par)
{
  const Int_t numPart = 7;
  // Signal region for the BG estimation
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

void postPP(TString fileList = "listQC.txt") // PP and QC task
{
  gROOT->SetBatch(kTRUE);
  // Files with histograms
  std::vector<std::string> name;
  std::vector<std::string> nameLegend;
  std::ifstream file(Form("%s", fileList.Data()));

  std::string remove = "/Users/rnepeiv/workLund/PhD_work/run3QCPbPb/qcTaskDev/results/combined/AnalysisResults_";
  std::string remove2 = ".root";
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

  cout << "List:" << endl;
  cout << fileList.Data() << endl;

  const Int_t numFiles = name.size();
  cout << "Number of files: " << numFiles << endl;

  TFile *fileIn[numFiles];
  TFile *fileOut[numFiles];

  const Int_t numParticles = 7;

  TDirectory *partDirectories[numFiles][numParticles];
  TDirectory *eventSeldir[numFiles];

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

  Float_t minRangeSignal[numParticles] = {0.48, 1.112, 1.112, 1.315, 1.315, 1.668, 1.668};
  Float_t maxRangeSignal[numParticles] = {0.51, 1.122, 1.120, 1.328, 1.328, 1.677, 1.677};
  Float_t minRange[numParticles] =       {0.435, 1.095, 1.10, 1.305, 1.305, 1.655, 1.655};
  Float_t maxRange[numParticles] =       {0.560, 1.140, 1.14, 1.340, 1.340, 1.690, 1.690};

  Float_t meanYLow[numParticles] = {0.493, 1.114, 1.114, 1.315, 1.315, 1.668, 1.668};
  Float_t meanYUp[numParticles] = {0.505, 1.117, 1.117, 1.328, 1.328, 1.677, 1.677};

  Float_t sigmaYLow[numParticles] = {0.003, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005};
  Float_t sigmaYUp[numParticles] = {0.015, 0.008, 0.008, 0.004, 0.004, 0.01, 0.01};

  Float_t yieldYLow[numParticles] = {1e-3, 1e-3, 1e-3, 1e-4, 1e-4, 1e-5, 1e-5};
  Float_t yieldYUp[numParticles] = {100, 10, 10, 0.1, 0.1, 0.01, 0.01};

  TH1F *hMeans[numFiles][numParticles];
  TH1F *hSigmas[numFiles][numParticles];
  TH1F *hYields[numFiles][numParticles];
  TH1F *hPurity[numFiles][numParticles];

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

  TCanvas *canvasNEvents[numFiles];
  TH1F *hNevents[numFiles];

  gStyle->SetErrorX(0);

  for (Int_t iFile = 0; iFile < numFiles; iFile++)
  {
    // Open input file
    fileIn[iFile] = TFile::Open(Form("%s", name[iFile].c_str()));
    if (!fileIn[iFile] || fileIn[iFile]->IsZombie()) {
        std::cerr << "Error opening input file!" << std::endl;
        return;
    }

    fileIn[iFile]->cd();

    eventSeldir[iFile] = fileIn[iFile]->GetDirectory("lf-strangenessqc/eventSelection");
    if (!eventSeldir[iFile])
    {
      std::cerr << "`lf-strangenessqc/eventSelection` directory is not found!" << std::endl;
      return;
    }

    hNevents[iFile] = (TH1F *)eventSeldir[iFile]->Get("hVertexZRec");
    if (!hNevents[iFile])
    {
      std::cerr << "Histogram `hVertexZRec` is not found!" << std::endl;
      return;
    }

    // Create output file
    fileOut[iFile] = new TFile(Form("postPPresults/%s.root", nameLegend[iFile].c_str()), "recreate");

    // Create dir to store all the fitting parameters
    TDirectory* fitParamsOutDir = fileOut[iFile]->mkdir("fitParams");
    // Create dir to store all the invariant mass histograms
    TDirectory* invMassOutDir = fileOut[iFile]->mkdir("invMassHists");

    canvasNEvents[numFiles] = new TCanvas(Form("%s_%d", hNevents[iFile]->GetName(), iFile), Form("%s_%d", hNevents[iFile]->GetName(), iFile), 1000, 800);
    hNevents[iFile]->Draw();
    canvasNEvents[numFiles]->Write();
    canvasNEvents[numFiles]->SaveAs(Form("postPPresults/%s.pdf(", nameLegend[iFile].c_str()));
    Float_t norm = hNevents[iFile]->GetEntries();

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
      hPurity[iFile][iPart] = new TH1F("Purity_" + invMassNames[iPart], "Purity_" + invMassNames[iPart], numPtBinsPart[iPart], ptBins[iPart]);
      for (Int_t iPt = 0; iPt < numPtBinsPart[iPart]; iPt++) {
        hInvMass2Dpt[iFile][iPart][iPt] = (TH2F*)hInvMass2D[iFile][iPart]->Clone(Form("2DHistInPtBin_%d_%d_%d", iFile, iPart, iPt));
        hInvMass2Dpt[iFile][iPart][iPt]->GetYaxis()->SetRangeUser(ptBins[iPart][iPt]+1e-6, ptBins[iPart][iPt+1]-1e-6);
        hInvMass1Dpt[iFile][iPart][iPt] = (TH1F *)hInvMass2Dpt[iFile][iPart][iPt]->ProjectionX();
        hInvMass1Dpt[iFile][iPart][iPt]->Sumw2();
        hInvMass1Dpt[iFile][iPart][iPt]->GetYaxis()->SetTitle("Counts");
        hInvMass1Dpt[iFile][iPart][iPt]->SetTitle(Form("%.2f-%.2f pt bin", ptBins[iPart][iPt], ptBins[iPart][iPt+1]));
        // Fit inv mass histogram
        canvasInvMass[iFile][iPart][iPt] = new TCanvas(Form("%s_%d_%d_%d", hInvMass3D[iFile][iPart]->GetName(), iFile, iPart, iPt), Form("%s_%d_%d_%d", hInvMass3D[iFile][iPart]->GetName(), iFile, iPart, iPt), 1000, 800);
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
        total[iFile][iPart][iPt]->SetParLimits(2, 0.001, 0.01);
        total[iFile][iPart][iPt]->SetParLimits(3, 0., 1.2 * peakValue);
        total[iFile][iPart][iPt]->SetParLimits(4, minRangeSignal[iPart], maxRangeSignal[iPart]);
        total[iFile][iPart][iPt]->SetParLimits(5, 0.001, 0.01);
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

        // Carolina code
        Double_t N1a = total[iFile][iPart][iPt]->GetParameter(0);      // N1
        Double_t N2a = total[iFile][iPart][iPt]->GetParameter(3);      // N2
        Double_t mu1a = total[iFile][iPart][iPt]->GetParameter(1);     // mu1
        Double_t mu2a = total[iFile][iPart][iPt]->GetParameter(4);     // mu2
        Double_t sigma1a = total[iFile][iPart][iPt]->GetParameter(2);  // sigma1
        Double_t sigma2a = total[iFile][iPart][iPt]->GetParameter(5);  // sigma2
        TMatrixD cova = fFitResultTotal[iFile][iPart][iPt]->GetCovarianceMatrix();

        Double_t mu_wa = (N1a*mu1a + N2a*mu2a)/(N1a+N2a);
        Double_t sigma_wa = (N1a*sigma1a + N2a*sigma2a)/(N1a+N2a);

        Double_t sa = N1a + N2a;
        Double_t wa_mu = N1a*mu1a + N2a*mu2a;
        Double_t wa_sigma = N1a*sigma1a + N2a*sigma2a; 

        Double_t mu_wa_step = pow((mu1a - mu2a),2)*(pow(N1a,2)*cova(3,3) + pow(N2a,2)*cova(0,0))
                    + 2*cova(0,3)*(wa_mu-sa*mu1a)*(wa_mu-sa*mu2a)
                    + pow(sa,2)*(pow(N1a,2)*cova(1,1) + pow(N2a,2)*cova(4,4) + 2*N1a*N2a*cova(1,4))
                    - 2*sa * (N1a * (cova(0,1)*(wa_mu-sa*mu1a) + cova(3,1)*(wa_mu-sa*mu2a)) + 
                          N2a * (cova(0,4)*(wa_mu-sa*mu1a) + cova(3,4)*(wa_mu-sa*mu2a)));

        mean[iFile][iPart][iPt] = mu_wa;
        errMean[iFile][iPart][iPt] = sqrt(mu_wa_step / pow(sa,4));

        Double_t sigma_wa_step = pow((sigma1a - sigma2a),2)*(pow(N1a,2)*cova(3,3) + pow(N2a,2)*cova(0,0))
                    + 2*cova(0,3)*(wa_sigma-sa*sigma1a)*(wa_sigma-sa*sigma2a)
                    + pow(sa,2)*(pow(N1a,2)*cova(2,2) + pow(N2a,2)*cova(5,5) + 2*N1a*N2a*cova(2,5))
                    - 2*sa * (N1a * (cova(0,2)*(wa_sigma-sa*sigma1a) + cova(3,2)*(wa_sigma-sa*sigma2a)) + 
                          N2a * (cova(0,5)*(wa_sigma-sa*sigma1a) + cova(3,5)*(wa_sigma-sa*sigma2a)));

        sigma[iFile][iPart][iPt] = sigma_wa;
        errSigma[iFile][iPart][iPt] = sqrt(sigma_wa_step / pow(sa,4));

        // cout << "sigma1: " << total[iFile][iPart][iPt]->GetParameter(2) << " err1: " << total[iFile][iPart][iPt]->GetParError(2) << " ampl1: " <<  total[iFile][iPart][iPt]->GetParameter(0) << std::endl;
        // cout << "sigma2: " << total[iFile][iPart][iPt]->GetParameter(5) << " err2: " << total[iFile][iPart][iPt]->GetParError(5) << " ampl2: " <<  total[iFile][iPart][iPt]->GetParameter(3) <<  std::endl;
        // cout << "total sigma: " << sigma[iFile][iPart][iPt] << " total err: " << errSigma[iFile][iPart][iPt] <<  std::endl;

        Double_t leftSignal = mean[iFile][iPart][iPt] - 6 * sigma[iFile][iPart][iPt];
        Double_t rightSignal = mean[iFile][iPart][iPt] + 6 * sigma[iFile][iPart][iPt];
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

        hPurity[iFile][iPart]->SetBinContent(iPt + 1, (yieldWithBG[iFile][iPart][iPt] - yieldBG[iFile][iPart][iPt])/(yieldWithBG[iFile][iPart][iPt]));
        hPurity[iFile][iPart]->SetBinError(iPt + 1, 0);
      }

      fitParamsOutDir->cd();

      StyleHisto(hMeans[iFile][iPart], meanYLow[iPart], meanYUp[iPart], kBlack, 20, "#it{p}_{T} (GeV/#it{c})", "Mean (GeV/#it{c}^{2})", "", 0, 0, 0, 1.0, 1.25, 1, 0.04, 0.04);
      hMeans[iFile][iPart]->Draw("P");
      gPad->SaveAs(Form("postPPresults/%s.pdf", nameLegend[iFile].c_str()));
      hMeans[iFile][iPart]->Write();

      StyleHisto(hSigmas[iFile][iPart], sigmaYLow[iPart], sigmaYUp[iPart], kBlack, 20, "#it{p}_{T} (GeV/#it{c})", "Sigma (GeV/#it{c}^{2})", "", 0, 0, 0, 1.0, 1.25, 1, 0.04, 0.04);
      hSigmas[iFile][iPart]->Draw("P");
      gPad->SaveAs(Form("postPPresults/%s.pdf", nameLegend[iFile].c_str()));
      hSigmas[iFile][iPart]->Write();

      StyleHisto(hPurity[iFile][iPart], 0, 1, kBlack, 20, "#it{p}_{T} (GeV/#it{c})", "S/(B+S)", "", 0, 0, 0, 1.0, 1.25, 1, 0.04, 0.04);
      hPurity[iFile][iPart]->Draw("P");
      gPad->SaveAs(Form("postPPresults/%s.pdf", nameLegend[iFile].c_str()));
      hPurity[iFile][iPart]->Write();

      StyleHisto(hYields[iFile][iPart], yieldYLow[iPart], yieldYUp[iPart], kBlack, 20, "#it{p}_{T} (GeV/#it{c})", "1/#it{N}_{evt} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}", "", 0, 0, 0, 1.0, 1.1, 1, 0.04, 0.04);
      hYields[iFile][iPart]->Draw("P");
      gPad->SetLogy();
      if (iPart == (numParticles - 1)) {
        gPad->SaveAs(Form("postPPresults/%s.pdf)", nameLegend[iFile].c_str()));
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