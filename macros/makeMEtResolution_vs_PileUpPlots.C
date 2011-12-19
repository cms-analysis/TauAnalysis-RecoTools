
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

#include <string>
#include <map>
#include <iostream>
#include <iomanip>
#include <fstream>

const double p0_data        =  0.8208;
const double p0ErrUp_data   =  0.0358;
const double p0ErrDown_data =  0.0358;
const double p1_data        =  0.7189;
const double p1ErrUp_data   =  0.0074;
const double p1ErrDown_data =  0.0074;
const double p2_data        = -0.0010;
const double p2ErrUp_data   =  0.0003;
const double p2ErrDown_data =  0.0003;

const double p0_mc          =  1.0002;
const double p0ErrUp_mc     =  0.0005;
const double p0ErrDown_mc   =  0.0005;
const double p1_mc          =  0.7067;
const double p1ErrUp_mc     =  0.0008;
const double p1ErrDown_mc   =  0.0008;
const double p2_mc          = -0.0018;
const double p2ErrUp_mc     =  0.0001;
const double p2ErrDown_mc   =  0.0001;

typedef std::map<std::string, double> valueMap1;
typedef std::map<std::string, valueMap1> valueMap2;
typedef std::map<std::string, valueMap2> valueMap3;

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
					  TObjArray& processes, const std::string& central_or_shift,
					  const TString& data_or_mcType, const TString& runPeriod, const TString& projection,
					  double p0, double p1, double p2)
{
  const int numVertices = 24;

  TGraphErrors* graph = new TGraphErrors(numVertices);

  for ( int iVertex = 1; iVertex <= numVertices; ++iVertex ) {

    TH2* histogram_sum = 0;
  
    int numProcesses = processes.GetEntries();
    for ( int iProcess = 0; iProcess < numProcesses; ++iProcess ) {
      TObjString* process = dynamic_cast<TObjString*>(processes.At(iProcess));
      
      TString dqmDirectory = process->GetString();
      if ( central_or_shift != "central" ) dqmDirectory.Append("/").Append(central_or_shift.data());
      TString histogramName_process = Form("%sVsQtNumVerticesEq%i", projection.Data(), iVertex);
      TH2* histogram_process = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectory, histogramName_process));
      
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

    //double x = iVertex;
    if ( (p1*p1 + 4*iVertex*p2 - 4*p0*p2) < 0. ) continue;
    double x = (-p1 + TMath::Sqrt(p1*p1 + 4*iVertex*p2 - 4*p0*p2))/(2*p2);
    if ( x < 0. ) continue;
    //std::cout << "p0 = " << p0 << ", p1 = " << p1 << ", p2 = " << p2 << ":" 
    //	        << " numVtx = " << iVertex << " --> numPU = " << x << std::endl;
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

TGraphErrors* compMEtResolution_vs_PileUp_data(TFile* inputFile, TF1* fit, 
					       const std::string& central_or_shift, const TString& runPeriod, const TString& projection)
{
  TObjArray processes;
  processes.Add(new TObjString("Data"));
  double p0 = p0_data;
  if ( central_or_shift == "vertexRecoEff_p0Up"   ) p0 += p0ErrUp_data;
  if ( central_or_shift == "vertexRecoEff_p0Down" ) p0 -= p0ErrDown_data;
  double p1 = p1_data;
  if ( central_or_shift == "vertexRecoEff_p1Up"   ) p1 += p1ErrUp_data;
  if ( central_or_shift == "vertexRecoEff_p1Down" ) p1 -= p1ErrDown_data;
  double p2 = p2_data;
  if ( central_or_shift == "vertexRecoEff_p2Up"   ) p2 += p2ErrUp_data;
  if ( central_or_shift == "vertexRecoEff_p2Down" ) p2 -= p2ErrDown_data;
  std::string central_or_shift_mod = ( central_or_shift.find("vertexRecoEff") == std::string::npos ) ?
    central_or_shift : "central";
  TGraphErrors* graph = compMEtResolution_vs_PileUp(inputFile, fit, processes, central_or_shift_mod , 
						    "data", runPeriod, projection, p0, p1, p2);
  if ( central_or_shift != "central" ) graph->SetName(Form("%s_%s", graph->GetName(), central_or_shift.data()));
  return graph;
}

TGraphErrors* compMEtResolution_vs_PileUp_mc(TFile* inputFile, TF1* fit, const std::string& central_or_shift,
					     const TString& runPeriod, const TString& projection)
{
  TObjArray processes;
  processes.Add(new TObjString("simDYtoMuMu"));
  processes.Add(new TObjString("simTTplusJets"));
  processes.Add(new TObjString("simWW"));
  processes.Add(new TObjString("simWZ"));
  processes.Add(new TObjString("simZZ"));
  processes.Add(new TObjString("simQCD"));
  double p0 = p0_mc;
  if ( central_or_shift == "vertexRecoEff_p0Up"   ) p0 += p0ErrUp_mc;
  if ( central_or_shift == "vertexRecoEff_p0Down" ) p0 -= p0ErrDown_mc;
  double p1 = p1_mc;
  if ( central_or_shift == "vertexRecoEff_p1Up"   ) p1 += p1ErrUp_mc;
  if ( central_or_shift == "vertexRecoEff_p1Down" ) p1 -= p1ErrDown_mc;
  double p2 = p2_mc;
  if ( central_or_shift == "vertexRecoEff_p2Up"   ) p2 += p2ErrUp_mc;
  if ( central_or_shift == "vertexRecoEff_p2Down" ) p2 -= p2ErrDown_mc;
  std::string central_or_shift_mod = ( central_or_shift.find("vertexRecoEff") == std::string::npos ) ?
    central_or_shift : "central";
  TGraphErrors* graph = compMEtResolution_vs_PileUp(inputFile, fit, processes, central_or_shift_mod , 
						    "mc", runPeriod, projection, p0, p1, p2);
  if ( central_or_shift != "central" ) graph->SetName(Form("%s_%s", graph->GetName(), central_or_shift.data()));
  return graph;
}

TF1* fitMEtResolution_vs_PileUp(TGraph* graph)
{
  std::cout << "fitting " << graph->GetName() << ":" << std::endl;
  TString fitName = TString(graph->GetName()).Append("_fit");
  TF1* fit = new TF1(fitName.Data(), "TMath::Sqrt([0]*[0] + TMath::Power(x, [1])*[2]*[2])", 0., 25.);
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
  TGraphErrors* graph_data_2011runA = 
    compMEtResolution_vs_PileUp_data(inputFile_2011runA, fit_2011runA_data, 
				     "central", "2011runA", projection);
  TF1* fit_data_2011runA = fitMEtResolution_vs_PileUp(graph_data_2011runA);
  
  TGraphErrors* graph_mc_2011runA = 
    compMEtResolution_vs_PileUp_mc(inputFile_2011runA, fit_2011runA_mc, 
				   "central", "2011runA", projection);
  TF1* fit_mc_2011runA = fitMEtResolution_vs_PileUp(graph_mc_2011runA);
  
  TGraphErrors* graph_data_2011runB = 
    compMEtResolution_vs_PileUp_data(inputFile_2011runB, fit_2011runB_data, 
				     "central", "2011runB", projection);
  TF1* fit_data_2011runB = fitMEtResolution_vs_PileUp(graph_data_2011runB);
  
  TGraphErrors* graph_mc_2011runB = 
    compMEtResolution_vs_PileUp_mc(inputFile_2011runB, fit_2011runB_mc, 	
				   "central", "2011runB", projection);
  TF1* fit_mc_2011runB = fitMEtResolution_vs_PileUp(graph_mc_2011runB);

  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 720);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.12);
  canvas->SetBottomMargin(0.12);

  TH1* dummyHistogram = new TH1D("dummyHistogram", "dummyHistogram", 26, -0.5, 25.5);
  dummyHistogram->SetTitle("");
  dummyHistogram->SetStats(false);
  dummyHistogram->SetMaximum(25.);
  dummyHistogram->SetMinimum(0.);
  
  TAxis* xAxis = dummyHistogram->GetXaxis();
  //xAxis->SetTitle("Num. rec. Vertices");
  xAxis->SetTitle("Num. Pile-Up Interactions");
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

double square(double x)
{
  return x*x;
}

void compErr(std::map<std::string, double>& p, double p_central_value, double pErr_central_value, 
	     TObjArray& sysUncertainties, double& errUp, double& errDown)
{
  double errUp2   = 0.;
  double errDown2 = 0.;

  int numSysUncertainties = sysUncertainties.GetEntries();
  assert((numSysUncertainties % 2) == 0);
  for ( int iSysUncertainty = 0; iSysUncertainty < (numSysUncertainties / 2); ++iSysUncertainty ) {
    TObjString* sysUncertaintyUp = dynamic_cast<TObjString*>(sysUncertainties.At(2*iSysUncertainty));
    TObjString* sysUncertaintyDown = dynamic_cast<TObjString*>(sysUncertainties.At(2*iSysUncertainty + 1));

    double pUp   = p[sysUncertaintyUp->GetString().Data()];
    double pDown = p[sysUncertaintyDown->GetString().Data()];

    double pMin  = TMath::Min(pUp, pDown);
    double pMax  = TMath::Max(pUp, pDown);

    if ( pMin < p_central_value ) errDown2 += square(p_central_value - pMin);
    if ( pMax > p_central_value ) errUp2   += square(pMax - p_central_value);
  }

//--- add in quadrature statistical uncertainty of fit
  errUp2   += square(pErr_central_value);
  errDown2 += square(pErr_central_value);

  errUp   = TMath::Sqrt(errUp2);
  errDown = TMath::Sqrt(errDown2);
}

typedef TGraphErrors* (compMEtResolution_function)(TFile*, TF1*, const std::string&, const TString&, const TString&);

void computeSysUncertainties(const TString& projection,
			     TFile* inputFile, TF1* fit, compMEtResolution_function* f, 
			     const TString& runPeriod, TObjArray& sysUncertainties,
			     valueMap3& resolution)
{
  TGraphErrors* graph_central_value = (*f)(inputFile, fit, "central", runPeriod, projection);
  TF1* fit_central_value = fitMEtResolution_vs_PileUp(graph_central_value);	  
						
  double p0_central_value    = fit_central_value->GetParameter(0);
  double p0Err_central_value = fit_central_value->GetParError(0);
  double p1_central_value    = fit_central_value->GetParameter(1);
  double p1Err_central_value = fit_central_value->GetParError(1);
  double p2_central_value    = fit_central_value->GetParameter(2);
  double p2Err_central_value = fit_central_value->GetParError(2);
    
  std::map<std::string, double> p0; // key = central or systematic uncertainty
  std::map<std::string, double> p1; 
  std::map<std::string, double> p2;

  int numSysUncertainties = sysUncertainties.GetEntries();
  for ( int iSysUncertainty = 0; iSysUncertainty < numSysUncertainties; ++iSysUncertainty ) {
    TObjString* sysUncertainty = dynamic_cast<TObjString*>(sysUncertainties.At(iSysUncertainty));
    TGraphErrors* graph_i = (*f)(inputFile, fit, sysUncertainty->GetString().Data(), runPeriod, projection);
    TF1* fit_i = fitMEtResolution_vs_PileUp(graph_i);
    p0[sysUncertainty->GetString().Data()] = fit_i->GetParameter(0);
    p1[sysUncertainty->GetString().Data()] = fit_i->GetParameter(1);
    p2[sysUncertainty->GetString().Data()] = fit_i->GetParameter(2);
  }

  double p0ErrUp, p0ErrDown;
  compErr(p0, p0_central_value, p0Err_central_value, sysUncertainties, p0ErrUp, p0ErrDown);
  resolution[projection.Data()]["sigmaZ"]["value"] = p0_central_value;
  resolution[projection.Data()]["sigmaZ"]["errUp"] = p0ErrUp;
  resolution[projection.Data()]["sigmaZ"]["errDown"] = p0ErrDown;
  std::cout << "sigmaZ = " << p0_central_value << " + " << p0ErrUp << " - " << p0ErrDown << std::endl;

  double p1ErrUp, p1ErrDown;
  compErr(p1, p1_central_value, p1Err_central_value, sysUncertainties, p1ErrUp, p1ErrDown);
  resolution[projection.Data()]["alpha"]["value"] = p1_central_value;
  resolution[projection.Data()]["alpha"]["errUp"] = p1ErrUp;
  resolution[projection.Data()]["alpha"]["errDown"] = p1ErrDown;
  std::cout << "alpha = " << p1_central_value << " + " << p1ErrUp << " - " << p1ErrDown << std::endl;

  double p2ErrUp, p2ErrDown;
  compErr(p2, p2_central_value, p2Err_central_value, sysUncertainties, p2ErrUp, p2ErrDown);
  resolution[projection.Data()]["sigmaMB"]["value"] = p2_central_value;
  resolution[projection.Data()]["sigmaMB"]["errUp"] = p2ErrUp;
  resolution[projection.Data()]["sigmaMB"]["errDown"] = p2ErrDown;
  std::cout << "sigmaMB = " << p2_central_value << " + " << p2ErrUp << " - " << p2ErrDown << std::endl;
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

void printResolutionRow(std::ofstream& table, 
			const std::string& projection,
			const std::string& parameterTitle, const std::string& parameterName,
			valueMap3& resolution_data, valueMap3& resolution_mc)
{
  table << "$" << parameterTitle << "$ & $" 
	<< std::fixed << std::setprecision(3)
	<< resolution_data[projection][parameterName]["value"]  
	<< "^{+" << resolution_data[projection][parameterName]["errUp"] << "}"
    	<< "_{-" << resolution_data[projection][parameterName]["errDown"] << "}$ & $" 
	<< resolution_mc[projection][parameterName]["value"]  
	<< "^{+" << resolution_mc[projection][parameterName]["errUp"] << "}"
    	<< "_{-" << resolution_mc[projection][parameterName]["errDown"] << "}$ \\\\" << std::endl;
}

void printResolutionTable(std::ofstream& table, 
			  const std::string& runPeriod, 
			  valueMap3& resolution_data, valueMap3& resolution_mc)
{
  table << "\\hline" << std::endl;
  table << "\\multicolumn{3}{|c|}{" << runPeriod << "} \\\\" << std::endl;
  table << "\\hline" << std::endl;
  printResolutionRow(table, "uParl", "\\sigma_{Z}^{\\parallel}",  "sigmaZ",  resolution_data, resolution_mc);
  printResolutionRow(table, "uParl", "\\sigma_{MB}^{\\parallel}", "sigmaMB", resolution_data, resolution_mc);
  printResolutionRow(table, "uParl", "\\alpha_{\\parallel}",      "alpha",   resolution_data, resolution_mc);
  table << "\\hline" << std::endl;
  printResolutionRow(table, "uPerp", "\\sigma_{Z}^{\\perp}",      "sigmaZ",  resolution_data, resolution_mc);
  printResolutionRow(table, "uPerp", "\\sigma_{MB}^{\\perp}",     "sigmaMB", resolution_data, resolution_mc);
  printResolutionRow(table, "uPerp", "\\alpha_{\\perp}",          "alpha",   resolution_data, resolution_mc);
  table << "\\hline" << std::endl;
}

void makeMEtResolution_vs_PileUpPlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  TString inputFilePath_2011runA = "/data1/veelken/tmp/ZllRecoilCorrection/v4_5/2011RunA/";

  TString inputFileName_pfMEt_2011runA = "analyzeZllRecoilCorrectionHistograms_all_pfMEtSmeared.root";
  TString inputFileName_pfMEtType1corrected_2011runA = "analyzeZllRecoilCorrectionHistograms_all_pfMEtTypeIcorrectedSmeared.root";

  TString inputFilePath_2011runB = "/data1/veelken/tmp/ZllRecoilCorrection/v4_5/2011RunB/";

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
  double fit_uParl_xMax = 300.;

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
				  "uParl", "RMS(u_{#parallel} ) / GeV", 
				  inputFile_pfMEt_2011runA, 
				  fit_uParl_pfMEt_2011runA_data, fit_uParl_pfMEt_2011runA_mc,
				  inputFile_pfMEt_2011runB, 
				  fit_uParl_pfMEt_2011runB_data, fit_uParl_pfMEt_2011runB_mc);
  makeMEtResolution_vs_PileUpPlot("rawPFMEt", 
				  "uPerp", "RMS(u_{#perp}  ) / GeV", 
				  inputFile_pfMEt_2011runA, 
				  0, 0, 
				  inputFile_pfMEt_2011runB, 
				  0, 0);
  makeMEtResolution_vs_PileUpPlot("Type1correctedPFMEt", 
				  "uParl", "RMS(u_{#parallel} ) / GeV", 
				  inputFile_pfMEtType1corrected_2011runA, 
				  fit_uParl_pfMEtType1corrected_2011runA_data, fit_uParl_pfMEtType1corrected_2011runA_mc,
				  inputFile_pfMEtType1corrected_2011runB, 
				  fit_uParl_pfMEtType1corrected_2011runA_data, fit_uParl_pfMEtType1corrected_2011runA_mc);
  makeMEtResolution_vs_PileUpPlot("Type1correctedPFMEt", 
				  "uPerp", "RMS(u_{#perp}  ) / GeV", 
				  inputFile_pfMEtType1corrected_2011runA, 
				  0, 0,
				  inputFile_pfMEtType1corrected_2011runB, 
				  0, 0);

  TObjArray sysUncertainties_data;
  sysUncertainties_data.Add(new TObjString("vertexRecoEff_p0Up"));
  sysUncertainties_data.Add(new TObjString("vertexRecoEff_p0Down"));
  sysUncertainties_data.Add(new TObjString("vertexRecoEff_p1Up"));
  sysUncertainties_data.Add(new TObjString("vertexRecoEff_p1Down"));
  sysUncertainties_data.Add(new TObjString("vertexRecoEff_p2Up"));
  sysUncertainties_data.Add(new TObjString("vertexRecoEff_p2Down"));

  valueMap3 resolution_pfMEt_2011runA_data; // key = 'uParl'/'uPerp', 'sigmaZ'/'sigmaMB'/'alpha', 'value'/'errUp'/'errDown'

  std::cout << "computing uncertainties on uParl, rawPFMEt, 2011 run A data:" << std::endl;
  computeSysUncertainties("uParl", inputFile_pfMEt_2011runA, 
			  fit_uParl_pfMEt_2011runA_data, &compMEtResolution_vs_PileUp_data,
			  "2011runA", sysUncertainties_data, resolution_pfMEt_2011runA_data);
  std::cout << "computing uncertainties on uPerp, rawPFMEt, 2011 run A data:" << std::endl;
  computeSysUncertainties("uPerp", inputFile_pfMEt_2011runA, 
			  0, &compMEtResolution_vs_PileUp_data,
			  "2011runA", sysUncertainties_data, resolution_pfMEt_2011runA_data);

  valueMap3 resolution_pfMEt_2011runB_data; // key = 'uParl'/'uPerp', 'sigmaZ'/'sigmaMB'/'alpha', 'value'/'errUp'/'errDown'

  std::cout << "computing uncertainties on uParl, rawPFMEt, 2011 run B data:" << std::endl;
  computeSysUncertainties("uParl", inputFile_pfMEt_2011runB, 
			  fit_uParl_pfMEt_2011runB_data, &compMEtResolution_vs_PileUp_data,
			  "2011runB", sysUncertainties_data, resolution_pfMEt_2011runB_data);
  std::cout << "computing uncertainties on uPerp, rawPFMEt, 2011 run B data:" << std::endl;
  computeSysUncertainties("uPerp", inputFile_pfMEt_2011runB, 
			  0, &compMEtResolution_vs_PileUp_data,
			  "2011runB", sysUncertainties_data, resolution_pfMEt_2011runB_data);

  valueMap3 resolution_pfMEtType1corrected_2011runA_data; // key = 'uParl'/'uPerp', 'sigmaZ'/'sigmaMB'/'alpha', 'value'/'errUp'/'errDown'

  std::cout << "computing uncertainties on uParl, Type1correctedPFMEt, 2011 run A data:" << std::endl;
  computeSysUncertainties("uParl", inputFile_pfMEtType1corrected_2011runA, 
			  fit_uParl_pfMEtType1corrected_2011runA_data, &compMEtResolution_vs_PileUp_data,
			  "2011runA", sysUncertainties_data, resolution_pfMEtType1corrected_2011runA_data);
  std::cout << "computing uncertainties on uPerp, Type1correctedPFMEt, 2011 run A data:" << std::endl;
  computeSysUncertainties("uPerp", inputFile_pfMEtType1corrected_2011runA, 
			  0, &compMEtResolution_vs_PileUp_data,
			  "2011runA", sysUncertainties_data, resolution_pfMEtType1corrected_2011runA_data);

  valueMap3 resolution_pfMEtType1corrected_2011runB_data; // key = 'uParl'/'uPerp', 'sigmaZ'/'sigmaMB'/'alpha', 'value'/'errUp'/'errDown'

  std::cout << "computing uncertainties on uParl, Type1correctedPFMEt, 2011 run B data:" << std::endl;
  computeSysUncertainties("uParl", inputFile_pfMEtType1corrected_2011runB, 
			  fit_uParl_pfMEtType1corrected_2011runB_data, &compMEtResolution_vs_PileUp_data,
			  "2011runB", sysUncertainties_data, resolution_pfMEtType1corrected_2011runB_data);
  std::cout << "computing uncertainties on uPerp, Type1correctedPFMEt, 2011 run B data:" << std::endl;
  computeSysUncertainties("uPerp", inputFile_pfMEtType1corrected_2011runB, 
			  0, &compMEtResolution_vs_PileUp_data,
			  "2011runB", sysUncertainties_data, resolution_pfMEtType1corrected_2011runB_data);

  TObjArray sysUncertainties_mc;
  sysUncertainties_mc.Add(new TObjString("vertexRecoEff_p0Up"));
  sysUncertainties_mc.Add(new TObjString("vertexRecoEff_p0Down"));
  sysUncertainties_mc.Add(new TObjString("vertexRecoEff_p1Up"));
  sysUncertainties_mc.Add(new TObjString("vertexRecoEff_p1Down"));
  sysUncertainties_mc.Add(new TObjString("vertexRecoEff_p2Up"));
  sysUncertainties_mc.Add(new TObjString("vertexRecoEff_p2Down"));
  sysUncertainties_mc.Add(new TObjString("jetEnUp"));
  sysUncertainties_mc.Add(new TObjString("jetEnDown"));
  sysUncertainties_mc.Add(new TObjString("jetResUp"));
  sysUncertainties_mc.Add(new TObjString("jetResDown"));
  sysUncertainties_mc.Add(new TObjString("unclEnUp"));
  sysUncertainties_mc.Add(new TObjString("unclEnDown"));
  sysUncertainties_mc.Add(new TObjString("vertexRecoEffUp"));
  sysUncertainties_mc.Add(new TObjString("vertexRecoEffDown"));

  valueMap3 resolution_pfMEt_2011runA_mc; // key = 'uParl'/'uPerp', 'sigmaZ'/'sigmaMB'/'alpha', 'value'/'errUp'/'errDown'

  std::cout << "computing uncertainties on uParl, rawPFMEt, Monte Carlo 2011 run A:" << std::endl;
  computeSysUncertainties("uParl", inputFile_pfMEt_2011runA, 
			  fit_uParl_pfMEt_2011runA_mc, &compMEtResolution_vs_PileUp_mc,
			  "2011runA", sysUncertainties_mc, resolution_pfMEt_2011runA_mc);
  std::cout << "computing uncertainties on uPerp, rawPFMEt, Monte Carlo 2011 run A:" << std::endl;
  computeSysUncertainties("uPerp", inputFile_pfMEt_2011runA, 
			  0, &compMEtResolution_vs_PileUp_mc,
			  "2011runA", sysUncertainties_mc, resolution_pfMEt_2011runA_mc);

  valueMap3 resolution_pfMEt_2011runB_mc; // key = 'uParl'/'uPerp', 'sigmaZ'/'sigmaMB'/'alpha', 'value'/'errUp'/'errDown'

  std::cout << "computing uncertainties on uParl, rawPFMEt, Monte Carlo 2011 run B:" << std::endl;
  computeSysUncertainties("uParl", inputFile_pfMEt_2011runB, 
			  fit_uParl_pfMEt_2011runB_mc, &compMEtResolution_vs_PileUp_mc,
			  "2011runB", sysUncertainties_mc, resolution_pfMEt_2011runB_mc);
  std::cout << "computing uncertainties on uPerp, rawPFMEt, Monte Carlo 2011 run B:" << std::endl;
  computeSysUncertainties("uPerp", inputFile_pfMEt_2011runB, 
			  0, &compMEtResolution_vs_PileUp_mc,
			  "2011runB", sysUncertainties_mc, resolution_pfMEt_2011runB_mc);

  valueMap3 resolution_pfMEtType1corrected_2011runA_mc; // key = 'uParl'/'uPerp', 'sigmaZ'/'sigmaMB'/'alpha', 'value'/'errUp'/'errDown'

  std::cout << "computing uncertainties on uParl, Type1correctedPFMEt, Monte Carlo 2011 run A:" << std::endl;
  computeSysUncertainties("uParl", inputFile_pfMEtType1corrected_2011runA, 
			  fit_uParl_pfMEtType1corrected_2011runA_mc, &compMEtResolution_vs_PileUp_mc,
			  "2011runA", sysUncertainties_mc, resolution_pfMEtType1corrected_2011runA_mc);
  std::cout << "computing uncertainties on uPerp, Type1correctedPFMEt, Monte Carlo 2011 run A:" << std::endl;
  computeSysUncertainties("uPerp", inputFile_pfMEtType1corrected_2011runA, 
			  0, &compMEtResolution_vs_PileUp_mc,
			  "2011runA", sysUncertainties_mc, resolution_pfMEtType1corrected_2011runA_mc);

  valueMap3 resolution_pfMEtType1corrected_2011runB_mc; // key = 'uParl'/'uPerp', 'sigmaZ'/'sigmaMB'/'alpha', 'value'/'errUp'/'errDown'

  std::cout << "computing uncertainties on uParl, Type1correctedPFMEt, Monte Carlo 2011 run B:" << std::endl;
  computeSysUncertainties("uParl", inputFile_pfMEtType1corrected_2011runB, 
			  fit_uParl_pfMEtType1corrected_2011runB_mc, &compMEtResolution_vs_PileUp_mc,
			  "2011runB", sysUncertainties_mc, resolution_pfMEtType1corrected_2011runB_mc);
  std::cout << "computing uncertainties on uPerp, Type1correctedPFMEt, Monte Carlo 2011 run B:" << std::endl;
  computeSysUncertainties("uPerp", inputFile_pfMEtType1corrected_2011runB, 
			  0, &compMEtResolution_vs_PileUp_mc,
			  "2011runB", sysUncertainties_mc, resolution_pfMEtType1corrected_2011runB_mc);

  std::ofstream* resolutionTable_pfMEt = 
    new std::ofstream("makeMEtResolution_vs_PileUpPlots_pfMEt.tex", std::ios::out);
  printResolutionTable(*resolutionTable_pfMEt, "Run period A", 
		       resolution_pfMEt_2011runA_data, resolution_pfMEt_2011runA_mc);
  printResolutionTable(*resolutionTable_pfMEt, "Run period B", 
		       resolution_pfMEt_2011runB_data, resolution_pfMEt_2011runB_mc);
  delete resolutionTable_pfMEt;

  std::ofstream* resolutionTable_pfMEtType1corrected = 
    new std::ofstream("makeMEtResolution_vs_PileUpPlots_pfMEtType1corrected.tex", std::ios::out);
  printResolutionTable(*resolutionTable_pfMEtType1corrected, "Run period A", 
		       resolution_pfMEtType1corrected_2011runA_data, resolution_pfMEtType1corrected_2011runA_mc);
  printResolutionTable(*resolutionTable_pfMEtType1corrected, "Run period B", 
		       resolution_pfMEtType1corrected_2011runB_data, resolution_pfMEtType1corrected_2011runB_mc);
  delete resolutionTable_pfMEtType1corrected;

  delete inputFile_pfMEt_2011runA;
  delete inputFile_pfMEt_2011runB;
  delete inputFile_pfMEtType1corrected_2011runA;
  delete inputFile_pfMEtType1corrected_2011runB;
}
