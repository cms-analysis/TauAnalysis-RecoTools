
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

TGraphErrors* compMEtResolution_vs_PileUp(TFile* inputFile, 
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
      
      TString histogramName = Form("%sVsQtEq%i", projection.Data(), iVertex);
      TH2* histogram_process = dynamic_cast<TH2*>(getHistogram(inputFile, process->GetString(), histogramName));
      
      if ( histogram_sum ) {
	histogram_sum->Add(histogram_process);
      } else {
	histogram_sum = (TH2*)histogram_process->Clone(TString(histogram_process->GetName()).Append("sum"));
      }
    }

    TH1* histogram_1d = histogram_sum->ProjectionY();

    double x = iVertex;
    double y = histogram_1d->GetRMS();
    double yErr = histogram_1d->GetRMSError();

    graph->SetPoint(iVertex - 1, x, y);
    graph->SetPointError(iVertex - 1, 0., yErr);
  }

  TString graphName = Form("graph_%s_%s_%s", data_or_mcType.Data(), runPeriod.Data(), projection.Data());
  graph->SetName(graphName.Data());

  return graph;
}

TGraphErrors* compMEtResolution_vs_PileUp_data(TFile* inputFile, const TString& runPeriod, const TString& projection)
{
  TObjArray processes;
  processes.Add(new TObjString("Data"));
  return compMEtResolution_vs_PileUp(inputFile, processes, "data", runPeriod, projection);
}

TGraphErrors* compMEtResolution_vs_PileUp_mc(TFile* inputFile, const TString& runPeriod, const TString& projection)
{
  TObjArray processes;
  processes.Add(new TObjString("simDYtoMuMu"));
  processes.Add(new TObjString("simTTplusJets"));
  processes.Add(new TObjString("simWW"));
  processes.Add(new TObjString("simWZ"));
  processes.Add(new TObjString("simZZ"));
  processes.Add(new TObjString("simQCD"));
  return compMEtResolution_vs_PileUp(inputFile, processes, "mc", runPeriod, projection);
}

TF1* fitMEtResolution_vs_PileUp(TGraph* graph)
{
  std::cout << "fitting " << graph->GetName() << ":" << std::endl;
  TString fitName = TString(graph->GetName()).Append("_fit");
  TF1* fit = new TF1(fitName.Data(), "TMath::Sqrt([0]*[0] + TMath::Power(x, [1])*[2]*[2])", 0.5, 24.5);
  fit->SetParameter(0, 8.);
  fit->SetParameter(1, 1.);
  fit->SetParameter(2, 4.);
  graph->Fit(fit, "0");
  return fit;
}

void makeMEtResolution_vs_PileUpPlot(const TString& metType, 
				     const TString& projection, const TString& yAxisTitle,
				     TFile* inputFile_2011runA, 
				     TFile* inputFile_2011runB)
{
  TGraphErrors* graph_data_2011runA = compMEtResolution_vs_PileUp_data(inputFile_2011runA, "2011runA", projection);
  TF1* fit_data_2011runA = fitMEtResolution_vs_PileUp(graph_data_2011runA);
  
  TGraphErrors* graph_mc_2011runA = compMEtResolution_vs_PileUp_mc(inputFile_2011runA, "2011runA", projection);
  TF1* fit_mc_2011runA = fitMEtResolution_vs_PileUp(graph_mc_2011runA);
  
  TGraphErrors* graph_data_2011runB = compMEtResolution_vs_PileUp_data(inputFile_2011runB, "2011runB", projection);
  TF1* fit_data_2011runB = fitMEtResolution_vs_PileUp(graph_data_2011runB);
  
  TGraphErrors* graph_mc_2011runB = compMEtResolution_vs_PileUp_mc(inputFile_2011runB, "2011runB", projection);
  TF1* fit_mc_2011runB = fitMEtResolution_vs_PileUp(graph_mc_2011runB);

  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 900);
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

  TLegend* legend = new TLegend(0.16, 0.71, 0.42, 0.89, "", "brNDC"); 
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

void makeMEtResolution_vs_PileUpPlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  TString inputFilePath_2011runA = "/data1/veelken/tmp/ZllRecoilCorrection/v4_3_1/2011RunA/";

  TString inputFileName_pfMEt_2011runA = "analyzeZllRecoilCorrectionHistograms_all_pfMEtSmeared.root";
  TString inputFileName_pfMEtType1corrected_2011runA = "analyzeZllRecoilCorrectionHistograms_all_pfMEtTypeIcorrectedSmeared.root";

  TString inputFilePath_2011runB = "/data1/veelken/tmp/ZllRecoilCorrection/v4_3/2011RunB/";

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
  
  makeMEtResolution_vs_PileUpPlot("rawPFMEt", 
				  "uParl", "rms(u_{#parallel} ) / GeV", 
				  inputFile_pfMEt_2011runA,
				  inputFile_pfMEt_2011runB);
  makeMEtResolution_vs_PileUpPlot("rawPFMEt", 
				  "uPerp", "rms(u_{#perp}  ) / GeV", 
				  inputFile_pfMEt_2011runA,
				  inputFile_pfMEt_2011runB);
  makeMEtResolution_vs_PileUpPlot("Type1correctedPFMEt", 
				  "uParl", "rms(u_{#parallel} ) / GeV", 
				  inputFile_pfMEtType1corrected_2011runA,
				  inputFile_pfMEtType1corrected_2011runB);
  makeMEtResolution_vs_PileUpPlot("Type1correctedPFMEt", 
				  "uPerp", "rms(u_{#perp}  ) / GeV", 
				  inputFile_pfMEtType1corrected_2011runA,
				  inputFile_pfMEtType1corrected_2011runB);

  delete inputFile_pfMEt_2011runA;
  delete inputFile_pfMEt_2011runB;
  delete inputFile_pfMEtType1corrected_2011runA;
  delete inputFile_pfMEtType1corrected_2011runB;
}
