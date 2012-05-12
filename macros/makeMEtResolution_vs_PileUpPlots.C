
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
#include <assert.h>

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

const bool make_plots            = true;
const bool draw1_data            = true;
const bool draw1_mc              = true;
const bool draw2_data            = true;
const bool draw2_mc              = true;
const bool compute_uncertainties = false;

enum { kGenNumPileUpInteractions, kRecVertices };
const int xAxis_mode = kRecVertices;

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

double square(double x)
{
  return x*x;
}

TGraphErrors* compMEtResolution_vs_PileUp(TFile* inputFile, 
					  TObjArray& processes, const std::string& central_or_shift,
					  const TString& data_or_mcType, const TString& runPeriod, const TString& projection,
					  double p0, double p1, double p2)
{
  int numVertices = 0;
  if      ( std::string(runPeriod.Data()) == "2011runA" ) numVertices = 24;
  else if ( std::string(runPeriod.Data()) == "2011runB" ) numVertices = 24;
  else if ( std::string(runPeriod.Data()) == "2012runA" ) numVertices = 34;
  else assert(0);

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

    TString histogramName_1d = Form("%s_1d", histogram_sum->GetName());
    TH1* histogram_1d = new TH1D(histogramName_1d.Data(), histogramName_1d.Data(), 75, -75., +75.);
    if ( !histogram_1d->GetSumw2N() ) histogram_1d->Sumw2();

    int numBinsX = histogram_sum->GetNbinsX();
    for ( int iBinX = 1; iBinX <= numBinsX; ++iBinX ) {
      TString histogramName_sliceX = Form("%s_sliceX%i", histogram_sum->GetName(), iBinX);
      TH1* histogram_sliceX = histogram_sum->ProjectionY(histogramName_sliceX.Data(), iBinX, iBinX, "e");

      int numBinsY = histogram_sum->GetNbinsY();
      for ( int iBinY = 1; iBinY <= numBinsY; ++iBinY ) {
	double y = histogram_sum->GetYaxis()->GetBinCenter(iBinY);

	double diff_y = y - histogram_sliceX->GetMean();

	int bin_1d = histogram_1d->FindBin(diff_y);

	double binContent_1d = histogram_1d->GetBinContent(bin_1d);
	double binError_1d = histogram_1d->GetBinError(bin_1d);

	double addBinContent = histogram_sum->GetBinContent(iBinX, iBinY);
	double addBinError = histogram_sum->GetBinError(iBinX, iBinY);

	histogram_1d->SetBinContent(bin_1d, binContent_1d + addBinContent);
	histogram_1d->SetBinError(bin_1d, TMath::Sqrt(binError_1d*binError_1d + addBinError*addBinError));
      }

      delete histogram_sliceX;
    }

    double x = -1.;
    if ( xAxis_mode == kGenNumPileUpInteractions ) {
      assert(p1 > 0.);
      if ( p2 != 0. ) {
	if ( (p1*p1 + 4*iVertex*p2 - 4*p0*p2) < 0. ) continue;
	x = (-p1 + TMath::Sqrt(p1*p1 + 4*iVertex*p2 - 4*p0*p2))/(2*p2);
      } else {
	x = (iVertex - p0)/p1;
      }
    } else if ( xAxis_mode == kRecVertices ) {
      x = iVertex;
    } else assert(0);
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

TGraphErrors* compMEtResolution_vs_PileUp_data(TFile* inputFile, 
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
  TGraphErrors* graph = compMEtResolution_vs_PileUp(inputFile, processes, central_or_shift_mod , 
						    "data", runPeriod, projection, p0, p1, p2);
  if ( central_or_shift != "central" ) graph->SetName(Form("%s_%s", graph->GetName(), central_or_shift.data()));
  return graph;
}

TGraphErrors* compMEtResolution_vs_PileUp_mc(TFile* inputFile, const std::string& central_or_shift,
					     const TString& runPeriod, const TString& projection)
{
  TObjArray processes;
  processes.Add(new TObjString("ZplusJets_madgraph"));
  //processes.Add(new TObjString("simTTplusJets"));
  //processes.Add(new TObjString("simWW"));
  //processes.Add(new TObjString("simWZ"));
  //processes.Add(new TObjString("simZZ"));
  //processes.Add(new TObjString("simQCD"));
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
  TGraphErrors* graph = compMEtResolution_vs_PileUp(inputFile, processes, central_or_shift_mod , 
						    "mc", runPeriod, projection, p0, p1, p2);
  if ( central_or_shift != "central" ) graph->SetName(Form("%s_%s", graph->GetName(), central_or_shift.data()));
  return graph;
}

TF1* fitMEtResolution_vs_PileUp(TGraph* graph, const TString& runPeriod)
{
  double xMax = 0;
  if      ( std::string(runPeriod.Data()) == "2011runA" ) xMax = 25.;
  else if ( std::string(runPeriod.Data()) == "2011runB" ) xMax = 25.;
  else if ( std::string(runPeriod.Data()) == "2012runA" ) xMax = 35.;
  else assert(0);

  std::cout << "fitting " << graph->GetName() << ":" << std::endl;
  TString fitName = TString(graph->GetName()).Append("_fit");
  TF1* fit = new TF1(fitName.Data(), "TMath::Sqrt([0]*[0] + TMath::Power(x, [1])*[2]*[2])", 0., xMax);
  fit->SetParameter(0, 8.);
  fit->SetParameter(1, 1.);
  //fit->FixParameter(1, 1.);
  fit->SetParameter(2, 4.);
  graph->Fit(fit, "0");
  return fit;
}

void cloneStyleOptions_graph(TGraph* graph, const TGraph* graph_ref)
{
  graph->SetMarkerStyle(graph_ref->GetMarkerStyle());
  graph->SetMarkerSize(graph_ref->GetMarkerSize());
  graph->SetMarkerColor(graph_ref->GetMarkerColor());
  graph->SetLineColor(graph_ref->GetLineColor());
  graph->SetLineWidth(graph_ref->GetLineWidth());
  graph->SetFillColor(graph_ref->GetFillColor());
  graph->SetFillStyle(graph_ref->GetFillStyle());
}

TGraphErrors* makeGraph_data_div_mc(TGraphErrors* graph_data, TGraphErrors* graph_mc)
{
  assert(graph_data->GetN() == graph_mc->GetN());

  int numPoints = graph_data->GetN();

  TGraphErrors* graph_diff = new TGraphErrors(numPoints);
  
  for ( int iPoint = 0; iPoint < numPoints; ++iPoint ) {
    double x_data, y_data;
    graph_data->GetPoint(iPoint, x_data, y_data);
    double xErr_data = graph_data->GetErrorX(iPoint);
    double yErr_data = graph_data->GetErrorY(iPoint);
    
    double x_mc, y_mc;
    graph_mc->GetPoint(iPoint, x_mc, y_mc);
    assert(x_mc == x_data);
    double xErr_mc = graph_mc->GetErrorX(iPoint);
    double yErr_mc = graph_mc->GetErrorY(iPoint);
    
    double y_diff = (y_data - y_mc)/y_mc;
    //std::cout << "x = " << x_data << ": y(data) = " << y_data << ", y(mc) = " << y_mc << " --> dy = " << y_diff << std::endl;
    double yErr2_diff = 0.;
    if ( y_data > 0. ) yErr2_diff += square(yErr_data/y_data);
    yErr2_diff += square(yErr_mc/y_mc);
    yErr2_diff *= square(y_data/y_mc);

    graph_diff->SetPoint(iPoint, x_data, y_diff);
    graph_diff->SetPointError(iPoint, xErr_data, TMath::Sqrt(yErr2_diff));
  }

  return graph_diff;
}

void cloneStyleOptions_fit(TF1* fit, const TF1* fit_ref)
{
  fit->SetLineColor(fit_ref->GetLineColor());
  fit->SetLineStyle(fit_ref->GetLineStyle());
  fit->SetLineWidth(fit_ref->GetLineWidth());
}

TF1* makeFit_data_div_mc(TF1* fit_data, TF1* fit_mc)
{
  TString fitName_diff = TString(fit_data->GetName()).Append("_div_").Append(fit_mc->GetName());
  TString formula_data = fit_data->GetTitle();
  int numParameter_data = fit_data->GetNpar();
  TString formula_mc = fit_mc->GetTitle();
  int numParameter_mc = fit_mc->GetNpar();
  for ( int iParameter_mc = numParameter_mc - 1; iParameter_mc >= 0; --iParameter_mc ) {
    TString parameterName_old = Form("[%i]", iParameter_mc);
    TString parameterName_new = Form("[%i]", iParameter_mc + numParameter_data);
    formula_mc.ReplaceAll(parameterName_old, parameterName_new);
  }
  TString formula_diff = TString("((").Append(formula_data).Append(")/(").Append(formula_mc).Append(")) - 1.0");
  double xMin_diff = TMath::Max(fit_data->GetXmin(), fit_mc->GetXmin());
  double xMax_diff = TMath::Min(fit_data->GetXmax(), fit_mc->GetXmax());
  TF1* fit_diff = new TF1(fitName_diff.Data(), formula_diff, xMin_diff, xMax_diff);
  for ( int iParameter_data = 0; iParameter_data < numParameter_data; ++iParameter_data ) {
    fit_diff->SetParameter(iParameter_data, fit_data->GetParameter(iParameter_data));
  }
  for ( int iParameter_mc = 0; iParameter_mc < numParameter_mc; ++iParameter_mc ) {
    fit_diff->SetParameter(iParameter_mc + numParameter_data, fit_mc->GetParameter(iParameter_mc));
  }
  return fit_diff;
}

void makeMEtResolution_vs_PileUpPlot(const TString& metType, 
				     const TString& projection, const TString& yAxisTitle,
				     TFile* inputFile1, const TString& runPeriod1, const TString& legendEntry1,
				     TFile* inputFile2, const TString& runPeriod2, const TString& legendEntry2)
{
  TGraphErrors* graph1_data = 0;
  TF1* fit1_data = 0;
  if ( draw1_data ) {
    graph1_data = compMEtResolution_vs_PileUp_data(inputFile1, "central", runPeriod1.Data(), projection);
    fit1_data = fitMEtResolution_vs_PileUp(graph1_data, runPeriod1);
  }

  TGraphErrors* graph1_mc = 0;
  TF1* fit1_mc = 0;
  if ( draw1_mc ) {
    graph1_mc = compMEtResolution_vs_PileUp_mc(inputFile1, "central", runPeriod1.Data(), projection);
    fit1_mc = fitMEtResolution_vs_PileUp(graph1_mc, runPeriod1);
  }

  TGraphErrors* graph2_data = 0;
  TF1* fit2_data = 0;
  if ( draw2_data ) {
    graph2_data = compMEtResolution_vs_PileUp_data(inputFile2, "central", runPeriod2.Data(), projection);
    fit2_data = fitMEtResolution_vs_PileUp(graph2_data, runPeriod2);
  }

  TGraphErrors* graph2_mc = 0;
  TF1* fit2_mc = 0;
  if ( draw2_mc ) {
    graph2_mc = compMEtResolution_vs_PileUp_mc(inputFile2, "central", runPeriod2.Data(), projection);
    fit2_mc = fitMEtResolution_vs_PileUp(graph2_mc, runPeriod2);
  }

  //TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 720);
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 1020);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.12);
  canvas->SetBottomMargin(0.12);

  TPad* topPad = new TPad("topPad", "topPad", 0.00, 0.35, 1.00, 1.00);
  topPad->SetFillColor(10);
  topPad->SetTopMargin(0.04);
  topPad->SetLeftMargin(0.15);
  topPad->SetBottomMargin(0.03);
  topPad->SetRightMargin(0.05);

  TPad* bottomPad = new TPad("bottomPad", "bottomPad", 0.00, 0.00, 1.00, 0.35);
  bottomPad->SetFillColor(10);
  bottomPad->SetTopMargin(0.02);
  bottomPad->SetLeftMargin(0.15);
  bottomPad->SetBottomMargin(0.24);
  bottomPad->SetRightMargin(0.05);

  canvas->cd();
  topPad->Draw();
  topPad->cd();

  double xAxisMax = 0.;
  if ( fit1_data ) xAxisMax = TMath::Max(xAxisMax, fit1_data->GetXmax());
  if ( fit1_mc   ) xAxisMax = TMath::Max(xAxisMax, fit1_mc->GetXmax());
  if ( fit2_data ) xAxisMax = TMath::Max(xAxisMax, fit2_data->GetXmax());
  if ( fit2_mc   ) xAxisMax = TMath::Max(xAxisMax, fit2_mc->GetXmax());  
  TH1* dummyHistogram_top = new TH1D("dummyHistogram_top", "dummyHistogram_top", TMath::Nint(xAxisMax) + 1, -0.5, xAxisMax + 0.5);
  dummyHistogram_top->SetTitle("");
  dummyHistogram_top->SetStats(false);
  dummyHistogram_top->SetMaximum(40.);
  dummyHistogram_top->SetMinimum(0.);
  
  TAxis* xAxis_top = dummyHistogram_top->GetXaxis();
  if      ( xAxis_mode == kGenNumPileUpInteractions ) xAxis_top->SetTitle("Num. Pile-Up Interactions");
  else if ( xAxis_mode == kRecVertices              ) xAxis_top->SetTitle("Num. rec. Vertices");
  else assert(0);

  xAxis_top->SetTitleOffset(1.15);
  xAxis_top->SetLabelColor(10);
  xAxis_top->SetTitleColor(10);

  TAxis* yAxis_top = dummyHistogram_top->GetYaxis();
  yAxis_top->SetTitle(yAxisTitle.Data());
  yAxis_top->SetTitleOffset(1.20);
  yAxis_top->SetTitleSize(0.06);

  dummyHistogram_top->Draw("axis");

  if ( draw1_data ) {
    graph1_data->SetLineColor(1);
    graph1_data->SetMarkerColor(1);
    graph1_data->SetMarkerStyle(20);
    graph1_data->Draw("p");
    
    fit1_data->SetLineColor(graph1_data->GetLineColor());
    fit1_data->SetLineWidth(2);
    fit1_data->Draw("same");
  }

  if ( draw1_mc ) {
    graph1_mc->SetLineColor(1);
    graph1_mc->SetMarkerColor(1);
    graph1_mc->SetMarkerStyle(24);
    graph1_mc->Draw("p");

    fit1_mc->SetLineColor(graph1_mc->GetLineColor());
    fit1_mc->SetLineStyle(2);
    fit1_mc->SetLineWidth(2);
    fit1_mc->Draw("same");
  }

  if ( draw2_data ) {
    graph2_data->SetLineColor(2);
    graph2_data->SetMarkerColor(2);
    graph2_data->SetMarkerStyle(21);
    graph2_data->Draw("p");
    
    fit2_data->SetLineColor(graph2_data->GetLineColor());
    fit2_data->SetLineWidth(2);
    fit2_data->Draw("same");
  }

  if ( draw2_mc ) {
    graph2_mc->SetLineColor(2);
    graph2_mc->SetMarkerColor(2);
    graph2_mc->SetMarkerStyle(25);
    graph2_mc->Draw("p");
    
    fit2_mc->SetLineColor(graph2_mc->GetLineColor());
    fit2_mc->SetLineStyle(2);
    fit2_mc->SetLineWidth(2);
    fit2_mc->Draw("same");
  }

  TLegend* legend = new TLegend(0.165, 0.635, 0.54, 0.87, "", "brNDC"); 
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  if ( draw1_data ) legend->AddEntry(graph1_data, TString(legendEntry1).Append(" Data"), "p");
  if ( draw1_mc   ) legend->AddEntry(graph1_mc,   TString(legendEntry1).Append(" MC"),   "p");
  if ( draw2_data ) legend->AddEntry(graph2_data, TString(legendEntry2).Append(" Data"), "p");
  if ( draw2_mc   ) legend->AddEntry(graph2_mc,   TString(legendEntry2).Append(" MC"),   "p");
  legend->Draw();	

  canvas->cd();
  bottomPad->Draw();
  bottomPad->cd();

  TH1* dummyHistogram_bottom = (TH1*)dummyHistogram_top->Clone("dummyHistogram_bottom");

  TAxis* xAxis_bottom = dummyHistogram_bottom->GetXaxis();
  xAxis_bottom->SetTitle(xAxis_top->GetTitle());
  xAxis_bottom->SetTitleOffset(1.20);
  xAxis_bottom->SetLabelColor(1);
  xAxis_bottom->SetTitleColor(1);
  xAxis_bottom->SetTitleSize(0.08);
  xAxis_bottom->SetLabelOffset(0.02);
  xAxis_bottom->SetLabelSize(0.08);
  xAxis_bottom->SetTickLength(0.055);

  TAxis* yAxis_bottom = dummyHistogram_bottom->GetYaxis();
  yAxis_bottom->SetTitle("#frac{Data - Simulation}{Simulation}");
  yAxis_bottom->SetTitleOffset(0.85);
  yAxis_bottom->SetNdivisions(505);
  yAxis_bottom->CenterTitle();
  yAxis_bottom->SetTitleSize(0.08);
  yAxis_bottom->SetLabelSize(0.08);
  yAxis_bottom->SetTickLength(0.04);

  double yDiffMax = +0.25;

  dummyHistogram_bottom->SetTitle("");
  dummyHistogram_bottom->SetStats(false);
  dummyHistogram_bottom->SetMaximum(+yDiffMax);
  dummyHistogram_bottom->SetMinimum(-yDiffMax);
  dummyHistogram_bottom->SetMarkerColor(1);
  dummyHistogram_bottom->SetMarkerStyle(20);
  dummyHistogram_bottom->Draw("axis");

  double xMin = +1.e+3;
  double xMax = -1.e+3;
  if ( fit1_data ) {
    xMin = TMath::Min(xMin, fit1_data->GetXmin());
    xMax = TMath::Max(xMax, fit1_data->GetXmax());
  } 
  if ( fit2_data ) {
    xMin = TMath::Min(xMin, fit2_data->GetXmin());
    xMax = TMath::Max(xMax, fit2_data->GetXmax());
  } 
  if ( xMax > xMin ) {
    TF1* line_at_zero = new TF1("line_at_zeto", "0.", xMin, xMax);
    line_at_zero->SetLineColor(8);
    line_at_zero->SetLineWidth(1);
    line_at_zero->Draw("same");
  }

  if ( draw1_data && draw1_mc ) {
    TGraphErrors* graph1_data_div_mc = makeGraph_data_div_mc(graph1_data, graph1_mc);
    cloneStyleOptions_graph(graph1_data_div_mc, graph1_data);
    graph1_data_div_mc->Draw("p");

    TF1* fit1_data_div_mc = makeFit_data_div_mc(fit1_data, fit1_mc);
    cloneStyleOptions_fit(fit1_data_div_mc, fit1_data);
    fit1_data_div_mc->Draw("same");
  }

  if ( draw2_data && draw2_mc ) {
    TGraphErrors* graph2_data_div_mc = makeGraph_data_div_mc(graph2_data, graph2_mc);
    cloneStyleOptions_graph(graph2_data_div_mc, graph2_data);
    graph2_data_div_mc->Draw("p");
    
    TF1* fit2_data_div_mc = makeFit_data_div_mc(fit2_data, fit2_mc);
    cloneStyleOptions_fit(fit2_data_div_mc, fit2_data);
    fit2_data_div_mc->Draw("same");
  }

  canvas->Update();
  std::string outputFileName = Form("metResolution_vs_PileUp_%s_%s.eps", metType.Data(), projection.Data());
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  canvas->Print(std::string(outputFileName_plot).append(".pdf").data());

  delete dummyHistogram_top;
  delete legend;
  delete topPad;
  delete dummyHistogram_bottom;
  delete bottomPad;
  delete canvas;
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

typedef TGraphErrors* (compMEtResolution_function)(TFile*, const std::string&, const TString&, const TString&);

void computeSysUncertainties(const TString& projection,
			     TFile* inputFile, compMEtResolution_function* f, 
			     const TString& runPeriod, TObjArray& sysUncertainties,
			     valueMap3& resolution)
{
  TGraphErrors* graph_central_value = (*f)(inputFile, "central", runPeriod, projection);
  TF1* fit_central_value = fitMEtResolution_vs_PileUp(graph_central_value, runPeriod);	  
						
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
    TGraphErrors* graph_i = (*f)(inputFile, sysUncertainty->GetString().Data(), runPeriod, projection);
    TF1* fit_i = fitMEtResolution_vs_PileUp(graph_i, runPeriod);
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

  TString inputFilePath1 = "/data1/veelken/tmp/ZllRecoilCorrection/v5_12_1/2012RunA/";
  TString inputFileName1 = "analyzeZllRecoilCorrectionHistograms_all_pfMEtTypeIcorrectedSmeared.root";
  TString legendEntry1   = "corr. PFMEt";
  TString runPeriod1     = "2012runA";

  TString inputFilePath2 = "/data1/veelken/tmp/ZllRecoilCorrection/v5_12_1/2012RunA/";
  TString inputFileName2 = "analyzeZllRecoilCorrectionHistograms_all_pfMEtByPhilsMVA.root";
  TString legendEntry2   = "MVA MEt";
  TString runPeriod2     = "2012runA";

  TFile* inputFile1 = new TFile(getFileName_full(inputFilePath1, inputFileName1));
  TFile* inputFile2 = new TFile(getFileName_full(inputFilePath2, inputFileName2));
  
  TString fit_uParl_formula = "-[0]*x*0.5*(1.0 - TMath::Erf(-[1]*TMath::Power(x, [2])))";

  double fit_uParl_xMin = 0.;
  double fit_uParl_xMax = 300.;

  makeMEtResolution_vs_PileUpPlot("pfMEtType1corr_vs_mvaMEt", 
				  "uParl", "RMS(u_{#parallel} ) / GeV", 
				  inputFile1, runPeriod1, legendEntry1,
				  inputFile2, runPeriod2, legendEntry2);
  makeMEtResolution_vs_PileUpPlot("pfMEtType1corr_vs_mvaMEt", 
				  "uPerp", "RMS(u_{#perp}  ) / GeV", 
				  inputFile1, runPeriod1, legendEntry1,
				  inputFile2, runPeriod2, legendEntry2);

  if ( compute_uncertainties ) {
    TObjArray sysUncertainties_data;
    sysUncertainties_data.Add(new TObjString("vertexRecoEff_p0Up"));
    sysUncertainties_data.Add(new TObjString("vertexRecoEff_p0Down"));
    sysUncertainties_data.Add(new TObjString("vertexRecoEff_p1Up"));
    sysUncertainties_data.Add(new TObjString("vertexRecoEff_p1Down"));
    sysUncertainties_data.Add(new TObjString("vertexRecoEff_p2Up"));
    sysUncertainties_data.Add(new TObjString("vertexRecoEff_p2Down"));
    
    valueMap3 resolution1_data; // key = 'uParl'/'uPerp', 'sigmaZ'/'sigmaMB'/'alpha', 'value'/'errUp'/'errDown'
    
    std::cout << "computing uncertainties on uParl, data (1):" << std::endl;
    computeSysUncertainties("uParl", inputFile1, 
			    &compMEtResolution_vs_PileUp_data,
			    runPeriod1, sysUncertainties_data, resolution1_data);
    std::cout << "computing uncertainties on uPerp, data (1):" << std::endl;
    computeSysUncertainties("uPerp", inputFile1, 
			    &compMEtResolution_vs_PileUp_data,
			    runPeriod1, sysUncertainties_data, resolution1_data);
    
    valueMap3 resolution2_data; // key = 'uParl'/'uPerp', 'sigmaZ'/'sigmaMB'/'alpha', 'value'/'errUp'/'errDown'
    
    std::cout << "computing uncertainties on uParl, data (2):" << std::endl;
    computeSysUncertainties("uParl", inputFile2, 
			    &compMEtResolution_vs_PileUp_data,
			    runPeriod2, sysUncertainties_data, resolution2_data);
    std::cout << "computing uncertainties on uPerp, data (2):" << std::endl;
    computeSysUncertainties("uPerp", inputFile2, 
			    &compMEtResolution_vs_PileUp_data,
			    runPeriod2, sysUncertainties_data, resolution2_data);
        
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
    
    valueMap3 resolution1_mc; // key = 'uParl'/'uPerp', 'sigmaZ'/'sigmaMB'/'alpha', 'value'/'errUp'/'errDown'
    
    std::cout << "computing uncertainties on uParl, Monte Carlo (1):" << std::endl;
    computeSysUncertainties("uParl", inputFile1, 
			    &compMEtResolution_vs_PileUp_mc,
			    runPeriod1, sysUncertainties_mc, resolution1_mc);
    std::cout << "computing uncertainties on uPerp, Monte Carlo (1):" << std::endl;
    computeSysUncertainties("uPerp", inputFile1, 
			    &compMEtResolution_vs_PileUp_mc,
			    runPeriod1, sysUncertainties_mc, resolution1_mc);
    
    valueMap3 resolution2_mc; // key = 'uParl'/'uPerp', 'sigmaZ'/'sigmaMB'/'alpha', 'value'/'errUp'/'errDown'
    
    std::cout << "computing uncertainties on uParl, Monte Carlo (2):" << std::endl;
    computeSysUncertainties("uParl", inputFile2, 
			    &compMEtResolution_vs_PileUp_mc,
			    runPeriod2, sysUncertainties_mc, resolution2_mc);
    std::cout << "computing uncertainties on uPerp, Monte Carlo (2):" << std::endl;
    computeSysUncertainties("uPerp", inputFile2, 
			    &compMEtResolution_vs_PileUp_mc,
			    runPeriod2, sysUncertainties_mc, resolution2_mc);
        
    std::ofstream* resolutionTable = new std::ofstream("makeMEtResolution_vs_PileUpPlots.tex", std::ios::out);
    printResolutionTable(*resolutionTable, runPeriod1.Data(), resolution1_data, resolution1_mc);
    printResolutionTable(*resolutionTable, runPeriod2.Data(), resolution2_data, resolution2_mc);
    delete resolutionTable;
  }

  delete inputFile1;
  delete inputFile2;
}
