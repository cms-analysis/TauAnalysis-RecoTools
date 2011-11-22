
#include <TFile.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TROOT.h>

#include <iostream>
#include <iomanip>

TH1* getHistogram(TFile* inputFile, const TString& dqmDirectory, const TString& meName)
{  
  //std::cout << "inputFile = " << inputFile->GetName() << std::endl;

  TString histogramName = TString(dqmDirectory);
  if ( !histogramName.EndsWith("/") ) histogramName.Append("/");
  histogramName.Append("beforeAddPUreweight/");
  histogramName.Append(meName);

  TH1* histogram = (TH1*)inputFile->Get(histogramName.Data());
  //std::cout << "histogramName = " << histogramName.Data() << ": histogram = " << histogram;
  //if ( histogram ) std::cout << ", integral = " << histogram->Integral();
  //std::cout << std::endl; 

  if ( !histogram->GetSumw2N() ) histogram->Sumw2();

  return histogram;
}

TGraphErrors* compMEtResolution_vs_PileUp(TFile* inputFile, TF1* fit, 
					  TObjArray& processes, 
					  const TString& data_or_mcType, const TString& runPeriod, const TString& projection)
{
  const int numVertices = 24;

  TGraphErrors* graph = new TGraphErrors(numVertices);

  for ( int iVertex = 1; iVertex <= numVertices; ++iVertex ) {

    TH2* histogram_sum = 0;
  
    int numProcesses = processes.GetEntries();
    for ( int iProcess = 0; iProcess < numProcesses; ++iProcess ) {
      TObjString* process = dynamic_cast<TObjString*>(processes.At(iProcess));
      
      TString histogramName_process = Form("%sVsQtNumVerticesEq%i", projection.Data(), iVertex);
      TH2* histogram_process = dynamic_cast<TH2*>(getHistogram(inputFile, process->GetString(), histogramName_process));
      
      if ( histogram_sum ) {
	histogram_sum->Add(histogram_process);
      } else {
	histogram_sum = (TH2*)histogram_process->Clone(TString(histogram_process->GetName()).Append("sum"));
      }
    }

    TH1* histogram_1d = 0;
    if ( fit ) {
      TString histogramName_1d = TString(histogram_sum->GetName()).Append("_1d");
      histogram_1d = new TH1D(histogramName_1d.Data(), histogramName_1d.Data(), 75, -75., +75.);
      if ( !histogram_1d->GetSumw2N() ) histogram_1d->Sumw2();

      int numBinsX = histogram_sum->GetNbinsX();
      for ( int iBinX = 1; iBinX <= numBinsX; ++iBinX ) {
	double x = histogram_sum->GetXaxis()->GetBinCenter(iBinX);

	int numBinsY = histogram_sum->GetNbinsY();
	for ( int iBinY = 1; iBinY <= numBinsY; ++iBinY ) {
	  double y = histogram_sum->GetYaxis()->GetBinCenter(iBinY);

	  double diff_y = y - fit->Eval(x);

	  int bin_1d = histogram_1d->FindBin(diff_y);

	  double binContent_1d = histogram_1d->GetBinContent(bin_1d);
	  double binError_1d = histogram_1d->GetBinError(bin_1d);

	  double addBinContent = histogram_sum->GetBinContent(iBinX, iBinY);
          double addBinError = histogram_sum->GetBinError(iBinX, iBinY);

	  histogram_1d->SetBinContent(bin_1d, binContent_1d + addBinContent);
	  histogram_1d->SetBinError(bin_1d, TMath::Sqrt(binError_1d*binError_1d + addBinError*addBinError));
	}
      }
    } else {
      histogram_1d = histogram_sum->ProjectionY();
    }

    double x = iVertex;
    double y = histogram_1d->GetRMS();
    double yErr = histogram_1d->GetRMSError();

    graph->SetPoint(iVertex - 1, x, y);
    graph->SetPointError(iVertex - 1, 0., yErr);

    delete histogram_1d;
  }

  TString graphName = Form("graph_%s_%s_%s", data_or_mcType.Data(), runPeriod.Data(), projection.Data());
  graph->SetName(graphName.Data());

  return graph;
}

TGraphErrors* compMEtResolution_vs_PileUp_data(TFile* inputFile, TF1* fit, const TString& runPeriod, const TString& projection)
{
  TObjArray processes;
  processes.Add(new TObjString("Data"));
  return compMEtResolution_vs_PileUp(inputFile, fit, processes, "data", runPeriod, projection);
}

TGraphErrors* compMEtResolution_vs_PileUp_mc(TFile* inputFile, TF1* fit, const TString& runPeriod, const TString& projection)
{
  TObjArray processes;
  processes.Add(new TObjString("simDYtoMuMu"));
  processes.Add(new TObjString("simTTplusJets"));
  processes.Add(new TObjString("simWW"));
  processes.Add(new TObjString("simWZ"));
  processes.Add(new TObjString("simZZ"));
  processes.Add(new TObjString("simQCD"));
  return compMEtResolution_vs_PileUp(inputFile, fit, processes, "mc", runPeriod, projection);
}

TF1* fitMEtResolution_vs_PileUp(TGraph* graph)
{
  std::cout << "fitting " << graph->GetName() << ":" << std::endl;
  TString fitName = TString(graph->GetName()).Append("_fit");
  TF1* fit = new TF1(fitName.Data(), "TMath::Sqrt([0]*[0] + TMath::Power(x, [1])*[2]*[2])", 0.5, 24.5);
  fit->SetParameter(0, 8.);
  fit->SetParameter(1, 1.);
  //fit->FixParameter(1, 1.);
  fit->SetParameter(2, 4.);
  graph->Fit(fit, "0");
  return fit;
}

void makeMEtResolution_vs_PileUpPlot(const TString& metType, 
				     const TString& projection, const TString& yAxisTitle,
				     TFile* inputFile_2011runA, 
				     TF1* fit_2011runA_data, TF1* fit_2011runA_mc,
				     TFile* inputFile_2011runB,
				     TF1* fit_2011runB_data, TF1* fit_2011runB_mc)
{
  TGraphErrors* graph_data_2011runA = compMEtResolution_vs_PileUp_data(inputFile_2011runA, fit_2011runA_data, "2011runA", projection);
  TF1* fit_data_2011runA = fitMEtResolution_vs_PileUp(graph_data_2011runA);
  
  TGraphErrors* graph_mc_2011runA = compMEtResolution_vs_PileUp_mc(inputFile_2011runA, fit_2011runA_mc, "2011runA", projection);
  TF1* fit_mc_2011runA = fitMEtResolution_vs_PileUp(graph_mc_2011runA);
  
  TGraphErrors* graph_data_2011runB = compMEtResolution_vs_PileUp_data(inputFile_2011runB, fit_2011runB_data, "2011runB", projection);
  TF1* fit_data_2011runB = fitMEtResolution_vs_PileUp(graph_data_2011runB);
  
  TGraphErrors* graph_mc_2011runB = compMEtResolution_vs_PileUp_mc(inputFile_2011runB, fit_2011runB_mc, "2011runB", projection);
  TF1* fit_mc_2011runB = fitMEtResolution_vs_PileUp(graph_mc_2011runB);

  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 720);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.12);
  canvas->SetBottomMargin(0.12);

  TH1* dummyHistogram = new TH1D("dummyHistogram", "dummyHistogram", 24, 0.5, 24.5);
  dummyHistogram->SetTitle("");
  dummyHistogram->SetStats(false);
  dummyHistogram->SetMaximum(25.);
  dummyHistogram->SetMinimum(0.);
  
  TAxis* xAxis = dummyHistogram->GetXaxis();
  xAxis->SetTitle("Num. rec. Vertices");
  xAxis->SetTitleOffset(1.15);

  TAxis* yAxis_top = dummyHistogram->GetYaxis();
  yAxis_top->SetTitle(yAxisTitle.Data());
  yAxis_top->SetTitleOffset(1.2);

  dummyHistogram->Draw("axis");

  graph_data_2011runA->SetLineColor(1);
  graph_data_2011runA->SetMarkerColor(1);
  graph_data_2011runA->SetMarkerStyle(20);
  graph_data_2011runA->Draw("p");

  fit_data_2011runA->SetLineColor(graph_data_2011runA->GetLineColor());
  fit_data_2011runA->SetLineWidth(2);
  fit_data_2011runA->Draw("same");

  graph_mc_2011runA->SetLineColor(1);
  graph_mc_2011runA->SetMarkerColor(1);
  graph_mc_2011runA->SetMarkerStyle(24);
  graph_mc_2011runA->Draw("p");

  fit_mc_2011runA->SetLineColor(graph_mc_2011runA->GetLineColor());
  fit_mc_2011runA->SetLineStyle(2);
  fit_mc_2011runA->SetLineWidth(2);
  fit_mc_2011runA->Draw("same");

  graph_data_2011runB->SetLineColor(2);
  graph_data_2011runB->SetMarkerColor(2);
  graph_data_2011runB->SetMarkerStyle(21);
  graph_data_2011runB->Draw("p");

  fit_data_2011runB->SetLineColor(graph_data_2011runB->GetLineColor());
  fit_data_2011runB->SetLineWidth(2);
  fit_data_2011runB->Draw("same");

  graph_mc_2011runB->SetLineColor(2);
  graph_mc_2011runB->SetMarkerColor(2);
  graph_mc_2011runB->SetMarkerStyle(25);
  graph_mc_2011runB->Draw("p");

  fit_mc_2011runB->SetLineColor(graph_mc_2011runB->GetLineColor());
  fit_mc_2011runB->SetLineStyle(2);
  fit_mc_2011runB->SetLineWidth(2);
  fit_mc_2011runB->Draw("same");

  TLegend* legend = new TLegend(0.145, 0.665, 0.48, 0.87, "", "brNDC"); 
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->AddEntry(graph_data_2011runA, "Run A Data", "p");
  legend->AddEntry(graph_mc_2011runA,   "sim. Run A", "p");
  legend->AddEntry(graph_data_2011runB, "Run B Data", "p");
  legend->AddEntry(graph_mc_2011runB,   "sim. Run B", "p");
  legend->Draw();	

  canvas->Update();
  TString outputFileName_plot = Form("metResolution_vs_PileUp_%s_%s.pdf", metType.Data(), projection.Data());
  canvas->Print(outputFileName_plot.Data());

  delete dummyHistogram;
  delete legend;
  delete canvas;
}

TString getFileName_full(const TString& path, const TString& fileName)
{
  TString fileName_full = path;
  if ( !fileName_full.EndsWith("/") ) fileName_full.Append("/");
  fileName_full.Append(fileName);
  return fileName_full;
}

void setFitParameter(TF1* fit, double r_had, double alpha, double beta)
{
  fit->SetParameter(0, r_had);
  fit->SetParameter(1, alpha);
  fit->SetParameter(2, beta);
}

void makeMEtResolution_vs_PileUpPlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  TString inputFilePath_2011runA = "/data1/veelken/tmp/ZllRecoilCorrection/v4_3_2/2011RunA/";

  TString inputFileName_pfMEt_2011runA = "analyzeZllRecoilCorrectionHistograms_all_pfMEtSmeared.root";
  TString inputFileName_pfMEtType1corrected_2011runA = "analyzeZllRecoilCorrectionHistograms_all_pfMEtTypeIcorrectedSmeared.root";

  TString inputFilePath_2011runB = "/data1/veelken/tmp/ZllRecoilCorrection/v4_3_1/2011RunB/";

  TString inputFileName_pfMEt_2011runB = "analyzeZllRecoilCorrectionHistograms_all_pfMEtSmeared.root";
  TString inputFileName_pfMEtType1corrected_2011runB = "analyzeZllRecoilCorrectionHistograms_all_pfMEtTypeIcorrectedSmeared.root";

  TFile* inputFile_pfMEt_2011runA = 
    new TFile(getFileName_full(inputFilePath_2011runA, 
			       inputFileName_pfMEt_2011runA));
  TFile* inputFile_pfMEt_2011runB = 
    new TFile(getFileName_full(inputFilePath_2011runB, 
			       inputFileName_pfMEt_2011runB));
  TFile* inputFile_pfMEtType1corrected_2011runA = 
    new TFile(getFileName_full(inputFilePath_2011runA, 
			       inputFileName_pfMEtType1corrected_2011runA));
  TFile* inputFile_pfMEtType1corrected_2011runB = 
    new TFile(getFileName_full(inputFilePath_2011runB, 
			       inputFileName_pfMEtType1corrected_2011runB));
  
  TString fit_uParl_formula = "-[0]*x*0.5*(1.0 - TMath::Erf(-[1]*TMath::Power(x, [2])))";
  double fit_uParl_xMin = 0.;
  double fit_uParl_xMax = 500.;

  TF1* fit_uParl_pfMEt_2011runA_data = 
    new TF1("fit_uParl_pfMEt_2011runA_data", 
	    fit_uParl_formula.Data(), fit_uParl_xMin, fit_uParl_xMax);
  setFitParameter(fit_uParl_pfMEt_2011runA_data, 0.938, 0.160, 0.542);
  TF1* fit_uParl_pfMEt_2011runA_mc = 
    new TF1("fit_uParl_pfMEt_2011runA_data", 
	    fit_uParl_formula.Data(), fit_uParl_xMin, fit_uParl_xMax);
  setFitParameter(fit_uParl_pfMEt_2011runA_mc, 0.952, 0.166, 0.560);

  TF1* fit_uParl_pfMEt_2011runB_data = 
    new TF1("fit_uParl_pfMEt_2011runB_data", 
	    fit_uParl_formula.Data(), fit_uParl_xMin, fit_uParl_xMax);
  setFitParameter(fit_uParl_pfMEt_2011runB_data, 0.941, 0.165, 0.523);
  TF1* fit_uParl_pfMEt_2011runB_mc = 
    new TF1("fit_uParl_pfMEt_2011runB_data", 
	    fit_uParl_formula.Data(), fit_uParl_xMin, fit_uParl_xMax);
  setFitParameter(fit_uParl_pfMEt_2011runB_mc, 0.953, 0.155, 0.568);

  TF1* fit_uParl_pfMEtType1corrected_2011runA_data = 
    new TF1("fit_uParl_pfMEtType1corrected_2011runA_data", 
	    fit_uParl_formula.Data(), fit_uParl_xMin, fit_uParl_xMax);
  setFitParameter(fit_uParl_pfMEtType1corrected_2011runA_data, 1.010, 0.114, 0.688);
  TF1* fit_uParl_pfMEtType1corrected_2011runA_mc = 
    new TF1("fit_uParl_pfMEtType1corrected_2011runA_data", 
	    fit_uParl_formula.Data(), fit_uParl_xMin, fit_uParl_xMax);
  setFitParameter(fit_uParl_pfMEtType1corrected_2011runA_mc, 1.004, 0.123, 0.688);

  TF1* fit_uParl_pfMEtType1corrected_2011runB_data = 
    new TF1("fit_uParl_pfMEtType1corrected_2011runB_data", 
	    fit_uParl_formula.Data(), fit_uParl_xMin, fit_uParl_xMax);
  setFitParameter(fit_uParl_pfMEtType1corrected_2011runB_data, 1.010, 0.123, 0.664);
  TF1* fit_uParl_pfMEtType1corrected_2011runB_mc = 
    new TF1("fit_uParl_pfMEtType1corrected_2011runB_data", 
	    fit_uParl_formula.Data(), fit_uParl_xMin, fit_uParl_xMax);
  setFitParameter(fit_uParl_pfMEtType1corrected_2011runB_mc, 1.004, 0.125, 0.673);

  makeMEtResolution_vs_PileUpPlot("rawPFMEt", 
				  "uParl", "rms(u_{#parallel} ) / GeV", 
				  inputFile_pfMEt_2011runA, 
				  fit_uParl_pfMEt_2011runA_data, fit_uParl_pfMEt_2011runA_mc,
				  inputFile_pfMEt_2011runB, 
				  fit_uParl_pfMEt_2011runB_data, fit_uParl_pfMEt_2011runB_mc);
  makeMEtResolution_vs_PileUpPlot("rawPFMEt", 
				  "uPerp", "rms(u_{#perp}  ) / GeV", 
				  inputFile_pfMEt_2011runA, 
				  0, 0, 
				  inputFile_pfMEt_2011runB, 
				  0, 0);
  makeMEtResolution_vs_PileUpPlot("Type1correctedPFMEt", 
				  "uParl", "rms(u_{#parallel} ) / GeV", 
				  inputFile_pfMEtType1corrected_2011runA, 
				  fit_uParl_pfMEtType1corrected_2011runA_data, fit_uParl_pfMEtType1corrected_2011runA_mc,
				  inputFile_pfMEtType1corrected_2011runB, 
				  fit_uParl_pfMEtType1corrected_2011runA_data, fit_uParl_pfMEtType1corrected_2011runA_mc);
  makeMEtResolution_vs_PileUpPlot("Type1correctedPFMEt", 
				  "uPerp", "rms(u_{#perp}  ) / GeV", 
				  inputFile_pfMEtType1corrected_2011runA, 
				  0, 0,
				  inputFile_pfMEtType1corrected_2011runB, 
				  0, 0);

  delete inputFile_pfMEt_2011runA;
  delete inputFile_pfMEt_2011runB;
  delete inputFile_pfMEtType1corrected_2011runA;
  delete inputFile_pfMEtType1corrected_2011runB;
}
