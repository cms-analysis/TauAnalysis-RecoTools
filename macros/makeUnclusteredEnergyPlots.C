
#include <TFile.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TPaveText.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TROOT.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <assert.h>
#include <math.h>
#include <limits>

enum { type_Undefined, type_Calo, type_Track, type_PF };

enum { mode_uParl, mode_uParl_div_qT };

TH1* getHistogram(TFile* inputFile, const std::string& dqmDirectory, const std::string& meName)
{  
  if ( !inputFile ) return 0;

  TString histogramName = "";
  if ( dqmDirectory != "" ) histogramName.Append(Form("/%s", dqmDirectory.data()));
  if ( histogramName.Length() > 0 && !histogramName.EndsWith("/") ) histogramName.Append("/");
  histogramName.Append(meName);

  TH1* histogram = (TH1*)inputFile->Get(histogramName.Data());
  std::cout << "histogramName = " << histogramName.Data() << ": histogram = " << histogram;
  if ( histogram ) std::cout << ", integral = " << histogram->Integral();
  std::cout << std::endl; 

  if ( histogram && !histogram->GetSumw2N() ) histogram->Sumw2();

  //if ( histogram->GetDimension() == 1 ) histogram->Rebin(5);

  //histogram->Scale(1./histogram->Integral());

  return histogram;
}

//-------------------------------------------------------------------------------
double square(double x)
{
  return x*x;
}

TH1* rebinHistogram(const TH1* histogram, unsigned numBinsMin_rebinned, double xMin, double xMax, bool normalize)
{
  TH1* histogram_rebinned = 0;

  if ( histogram ) {
    TAxis* xAxis = histogram->GetXaxis();

    unsigned numBins = xAxis->GetNbins();
    unsigned numBins_withinRange = 0;
    double binWidth = 0.;
    bool isVariableBinning = false;
    for ( unsigned iBin = 1; iBin <= numBins; ++iBin ) {
      double binCenter = xAxis->GetBinCenter(iBin);
      if ( binCenter >= xMin && binCenter <= xMax ) ++numBins_withinRange;
      if ( iBin == 1 ) binWidth = xAxis->GetBinWidth(iBin);
      else if ( TMath::Abs(xAxis->GetBinWidth(iBin) - binWidth) > (1.e-3*binWidth) ) isVariableBinning = true;
    }

    std::cout << "histogram = " << histogram->GetName() << ":" 
              << " numBins(" << xMin << ".." << "xMax) = " << numBins_withinRange << ", integral = " << histogram->Integral() << std::endl;
    
    unsigned numBins_rebinned = numBins_withinRange;
    if ( !isVariableBinning ) {
      for ( int combineNumBins = 5; combineNumBins >= 2; --combineNumBins ) {
	if ( numBins_withinRange >= (combineNumBins*numBinsMin_rebinned) && (numBins % combineNumBins) == 0 ) {
	  numBins_rebinned /= combineNumBins;
	  numBins_withinRange /= combineNumBins;
	}
      }
    }

    std::string histogramName_rebinned = std::string(histogram->GetName()).append("_rebinned");
    if ( !isVariableBinning ) {
      histogram_rebinned = new TH1D(histogramName_rebinned.data(), histogram->GetTitle(), numBins_rebinned, xMin, xMax);
    } else {
      TArrayF binning_withinRange(numBins_withinRange + 1);
      unsigned iBin_withinRange = 1;
      for ( unsigned iBin = 1; iBin <= numBins; ++iBin ) {
	double binCenter = xAxis->GetBinCenter(iBin);
	if ( binCenter >= xMin && binCenter <= xMax ) {
	  if ( iBin_withinRange == 1 ) binning_withinRange[0] = xAxis->GetBinLowEdge(iBin);
	  binning_withinRange[iBin_withinRange] = xAxis->GetBinUpEdge(iBin);
	  ++iBin_withinRange;
	}
      }
      histogram_rebinned = new TH1D(histogramName_rebinned.data(), histogram->GetTitle(), numBins_rebinned, binning_withinRange.GetArray());
    }
    if ( !histogram_rebinned->GetSumw2N() ) histogram_rebinned->Sumw2();

    TAxis* xAxis_rebinned = histogram_rebinned->GetXaxis();
      
    unsigned iBin = 1;
    for ( unsigned iBin_rebinned = 1; iBin_rebinned <= numBins_rebinned; ++iBin_rebinned ) {
      double binContent_rebinned = 0.;
      double binError2_rebinned = 0.;

      double xMin_rebinnedBin = xAxis_rebinned->GetBinLowEdge(iBin_rebinned);
      double xMax_rebinnedBin = xAxis_rebinned->GetBinUpEdge(iBin_rebinned);

      while ( histogram->GetBinCenter(iBin) < xMin_rebinnedBin ) {
	++iBin;
      }

      while ( histogram->GetBinCenter(iBin) >= xMin_rebinnedBin && histogram->GetBinCenter(iBin) < xMax_rebinnedBin ) {
	binContent_rebinned += histogram->GetBinContent(iBin);
	binError2_rebinned += square(histogram->GetBinError(iBin));
	++iBin;
      }

      histogram_rebinned->SetBinContent(iBin_rebinned, binContent_rebinned);
      histogram_rebinned->SetBinError(iBin_rebinned, TMath::Sqrt(binError2_rebinned));
    }

    if ( normalize ) {
      if ( !histogram_rebinned->GetSumw2N() ) histogram_rebinned->Sumw2();
      histogram_rebinned->Scale(1./histogram_rebinned->Integral());
    }

    std::cout << "histogram(rebinned) = " << histogram_rebinned->GetName() << ":" 
              << " numBins = " << histogram_rebinned->GetNbinsX() << ", integral = " << histogram_rebinned->Integral() << std::endl;
  }

  return histogram_rebinned;
}

TH1* compRatioHistogram(const std::string& ratioHistogramName, const TH1* numerator, const TH1* denominator, bool subtractOne = true)
{
  assert(numerator->GetDimension() == denominator->GetDimension());
  assert(numerator->GetNbinsX() == denominator->GetNbinsX());

  TH1* histogramRatio = (TH1*)numerator->Clone(ratioHistogramName.data());
  histogramRatio->Divide(denominator);

  if ( subtractOne ) {
    int nBins = histogramRatio->GetNbinsX();
    for ( int iBin = 1; iBin <= nBins; ++iBin ){
      double binContent = histogramRatio->GetBinContent(iBin);
      histogramRatio->SetBinContent(iBin, binContent - 1.);
    }
  }

  histogramRatio->SetLineColor(numerator->GetLineColor());
  histogramRatio->SetLineWidth(numerator->GetLineWidth());
  histogramRatio->SetMarkerColor(numerator->GetMarkerColor());
  histogramRatio->SetMarkerStyle(numerator->GetMarkerStyle());

  return histogramRatio;
}

void showDistribution(double canvasSizeX, double canvasSizeY,
		      TH1* histogram_ref, const std::string& legendEntry_ref,
		      TH1* histogram2, const std::string& legendEntry2,
		      TH1* histogram3, const std::string& legendEntry3,
		      TH1* histogram4, const std::string& legendEntry4,
		      TH1* histogram5, const std::string& legendEntry5,
		      TH1* histogram6, const std::string& legendEntry6,
		      double xMin, double xMax, unsigned numBinsMin_rebinned, const std::string& xAxisTitle, double xAxisOffset,
		      bool useLogScale, double yMin, double yMax, double yMin_ratio, double yMax_ratio, const std::string& yAxisTitle, double yAxisOffset,
		      double legendX0, double legendY0, 
		      const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
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
  topPad->SetLogy(useLogScale);

  TPad* bottomPad = new TPad("bottomPad", "bottomPad", 0.00, 0.00, 1.00, 0.35);
  bottomPad->SetFillColor(10);
  bottomPad->SetTopMargin(0.02);
  bottomPad->SetLeftMargin(0.15);
  bottomPad->SetBottomMargin(0.24);
  bottomPad->SetRightMargin(0.05);
  bottomPad->SetLogy(false);

  canvas->cd();
  topPad->Draw();
  topPad->cd();

  int colors[6] = { 1, 2, 3, 4, 6, 7 };
  int markerStyles[6] = { 22, 32, 20, 24, 21, 25 };

  TLegend* legend = new TLegend(legendX0, legendY0, legendX0 + 0.42, legendY0 + 0.22, "", "brNDC"); 
  legend->SetBorderSize(0);
  legend->SetFillColor(0);

  TH1* histogram_ref_rebinned = rebinHistogram(histogram_ref, numBinsMin_rebinned, xMin, xMax, true);
  histogram_ref_rebinned->SetTitle("");
  histogram_ref_rebinned->SetStats(false);
  histogram_ref_rebinned->SetMinimum(yMin);
  histogram_ref_rebinned->SetMaximum(yMax);
  histogram_ref_rebinned->SetLineColor(colors[0]);
  histogram_ref_rebinned->SetLineWidth(2);
  histogram_ref_rebinned->SetMarkerColor(colors[0]);
  histogram_ref_rebinned->SetMarkerStyle(markerStyles[0]);
  histogram_ref_rebinned->Draw("e1p");
  legend->AddEntry(histogram_ref_rebinned, legendEntry_ref.data(), "p");

  TAxis* xAxis_top = histogram_ref_rebinned->GetXaxis();
  xAxis_top->SetTitle(xAxisTitle.data());
  xAxis_top->SetTitleOffset(xAxisOffset);
  xAxis_top->SetLabelColor(10);
  xAxis_top->SetTitleColor(10);

  TAxis* yAxis_top = histogram_ref_rebinned->GetYaxis();
  yAxis_top->SetTitle(yAxisTitle.data());
  yAxis_top->SetTitleOffset(yAxisOffset);

  TH1* histogram2_rebinned = 0;
  if ( histogram2 ) {
    histogram2_rebinned = rebinHistogram(histogram2, numBinsMin_rebinned, xMin, xMax, true);
    histogram2_rebinned->SetLineColor(colors[1]);
    histogram2_rebinned->SetLineWidth(2);
    histogram2_rebinned->SetMarkerColor(colors[1]);
    histogram2_rebinned->SetMarkerStyle(markerStyles[1]);
    histogram2_rebinned->Draw("e1psame");
    legend->AddEntry(histogram2_rebinned, legendEntry2.data(), "p");
  }

  TH1* histogram3_rebinned = 0;
  if ( histogram3 ) {
    histogram3_rebinned = rebinHistogram(histogram3, numBinsMin_rebinned, xMin, xMax, true);
    histogram3_rebinned->SetLineColor(colors[2]);
    histogram3_rebinned->SetLineWidth(2);
    histogram3_rebinned->SetMarkerColor(colors[2]);
    histogram3_rebinned->SetMarkerStyle(markerStyles[2]);
    histogram3_rebinned->Draw("e1psame");
    legend->AddEntry(histogram3_rebinned, legendEntry3.data(), "p");
  }

  TH1* histogram4_rebinned = 0;
  if ( histogram4 ) {
    histogram4_rebinned = rebinHistogram(histogram4, numBinsMin_rebinned, xMin, xMax, true);
    histogram4_rebinned->SetLineColor(colors[3]);
    histogram4_rebinned->SetLineWidth(2);
    histogram4_rebinned->SetMarkerColor(colors[3]);
    histogram4_rebinned->SetMarkerStyle(markerStyles[3]);
    histogram4_rebinned->Draw("e1psame");
    legend->AddEntry(histogram4_rebinned, legendEntry4.data(), "p");
  }

  TH1* histogram5_rebinned = 0;
  if ( histogram5 ) {
    histogram5_rebinned = rebinHistogram(histogram5, numBinsMin_rebinned, xMin, xMax, true);
    histogram5_rebinned->SetLineColor(colors[4]);
    histogram5_rebinned->SetLineWidth(2);
    histogram5_rebinned->SetMarkerColor(colors[4]);
    histogram5_rebinned->SetMarkerStyle(markerStyles[4]);
    histogram5_rebinned->Draw("e1psame");
    legend->AddEntry(histogram5_rebinned, legendEntry5.data(), "p");
  }

  TH1* histogram6_rebinned = 0;
  if ( histogram6 ) {
    histogram6_rebinned = rebinHistogram(histogram6, numBinsMin_rebinned, xMin, xMax, true);
    histogram6_rebinned->SetLineColor(colors[5]);
    histogram6_rebinned->SetLineWidth(2);
    histogram6_rebinned->SetMarkerColor(colors[5]);
    histogram6_rebinned->SetMarkerStyle(markerStyles[5]);
    histogram6_rebinned->Draw("e1psame");
    legend->AddEntry(histogram6_rebinned, legendEntry6.data(), "p");
  }

  legend->Draw();

  canvas->cd();
  bottomPad->Draw();
  bottomPad->cd();

  TH1* histogram2_div_ref = 0;
  if ( histogram2 ) {
    std::string histogramName2_div_ref = std::string(histogram2->GetName()).append("_div_").append(histogram_ref->GetName());
    histogram2_div_ref = compRatioHistogram(histogramName2_div_ref, histogram2_rebinned, histogram_ref_rebinned);
    histogram2_div_ref->SetTitle("");
    histogram2_div_ref->SetStats(false);
    histogram2_div_ref->SetMinimum(yMin_ratio);
    histogram2_div_ref->SetMaximum(yMax_ratio);

    TAxis* xAxis_bottom = histogram2_div_ref->GetXaxis();
    xAxis_bottom->SetTitle(xAxis_top->GetTitle());
    xAxis_bottom->SetLabelColor(1);
    xAxis_bottom->SetTitleColor(1);
    xAxis_bottom->SetTitleOffset(1.20);
    xAxis_bottom->SetTitleSize(0.08);
    xAxis_bottom->SetLabelOffset(0.02);
    xAxis_bottom->SetLabelSize(0.08);
    xAxis_bottom->SetTickLength(0.055);
    
    TAxis* yAxis_bottom = histogram2_div_ref->GetYaxis();
    yAxis_bottom->SetTitle("#frac{Data - MC}{Data}");
    yAxis_bottom->SetTitleOffset(0.70);
    yAxis_bottom->SetNdivisions(505);
    yAxis_bottom->CenterTitle();
    yAxis_bottom->SetTitleSize(0.08);
    yAxis_bottom->SetLabelSize(0.08);
    yAxis_bottom->SetTickLength(0.04);  
  
    histogram2_div_ref->Draw("e1p");
  }

  TH1* histogram3_div_ref = 0;
  if ( histogram3 ) {
    std::string histogramName3_div_ref = std::string(histogram3->GetName()).append("_div_").append(histogram_ref->GetName());
    histogram3_div_ref = compRatioHistogram(histogramName3_div_ref, histogram3_rebinned, histogram_ref_rebinned);
    histogram3_div_ref->Draw("e1psame");
  }

  TH1* histogram4_div_ref = 0;
  if ( histogram4 ) {
    std::string histogramName4_div_ref = std::string(histogram4->GetName()).append("_div_").append(histogram_ref->GetName());
    histogram4_div_ref = compRatioHistogram(histogramName4_div_ref, histogram4_rebinned, histogram_ref_rebinned);
    histogram4_div_ref->Draw("e1psame");
  }

  TH1* histogram5_div_ref = 0;
  if ( histogram5 ) {
    std::string histogramName5_div_ref = std::string(histogram5->GetName()).append("_div_").append(histogram_ref->GetName());
    histogram5_div_ref = compRatioHistogram(histogramName5_div_ref, histogram5_rebinned, histogram_ref_rebinned);
    histogram5_div_ref->Draw("e1psame");
  }

  TH1* histogram6_div_ref = 0;
  if ( histogram6 ) {
    std::string histogramName6_div_ref = std::string(histogram6->GetName()).append("_div_").append(histogram_ref->GetName());
    histogram6_div_ref = compRatioHistogram(histogramName6_div_ref, histogram6_rebinned, histogram_ref_rebinned);
    histogram6_div_ref->Draw("e1psame");
  }
  
  canvas->Update();
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  if ( useLogScale ) outputFileName_plot.append("_log");
  else outputFileName_plot.append("_linear");
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  //canvas->Print(std::string(outputFileName_plot).append(".root").data());
  
  delete legend;
  delete histogram2_div_ref;
  delete histogram3_div_ref;
  delete histogram4_div_ref;
  delete histogram5_div_ref;
  delete histogram6_div_ref;
  delete topPad;
  delete bottomPad;
  delete canvas;  
}
//-------------------------------------------------------------------------------

enum { kMean, kRMS };

double getHistogramAvWithinRange(const TH1* histogram, double xMin, double xMax)
{
  int binLowIndex = const_cast<TH1*>(histogram)->FindBin(xMin);
  int binUpIndex  = const_cast<TH1*>(histogram)->FindBin(xMax);
  histogram->GetXaxis()->SetRange(binLowIndex, binUpIndex);

  double xAv = histogram->GetMean();
  
  // reset x-axis range selection 
  histogram->GetXaxis()->SetRange(1., 0.);

  return xAv;
}

TGraphAsymmErrors* makeGraph_mean_or_rms(const std::string& name, const std::string& title, 
					 const TH2* histogram_uParl_vs_eta, const TH1* histogram_qT, double qTmin, double qTmax, int mode, int mean_or_rms, bool divideByBinWidth)
{
  //std::cout << "<makeGraph_mean_or_rms>:" << std::endl;
  //std::cout << " name = " << name << std::endl;

  if ( !(histogram_uParl_vs_eta && histogram_qT) ) return 0;

  double qTav = getHistogramAvWithinRange(histogram_qT, qTmin, qTmax);
  //std::cout << "qT: min = " << qTmin << ", max = " << qTmax << ", av = " << qTav << std::endl;
  
  const TAxis* etaAxis = histogram_uParl_vs_eta->GetXaxis();
  int etaNumBins = etaAxis->GetNbins();

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(etaNumBins);
  graph->SetName(Form("%s_qT%1.1fto%1.1f", name.data(), qTmin, qTmax));
  graph->SetTitle(Form("%s (%1.1f < q_{T} < %1.1f)", title.data(), qTmin, qTmax));  

  for ( int etaBin = 1; etaBin <= etaNumBins; ++etaBin ) {
    double etaMin = etaAxis->GetBinLowEdge(etaBin);
    double etaMax = etaAxis->GetBinUpEdge(etaBin);

    double x        = etaAxis->GetBinCenter(etaBin);
    double xErrUp   = etaMax - x;
    double xErrDown = x - etaMin;

    //std::cout << " eta: min = " << etaMin << ", max = " << etaMax << ", center = " << x << std::endl;

    TString histogramName_uParl_proj = Form("%s_py_etaBin%i", histogram_uParl_vs_eta->GetName(), etaBin);
    TH1D* histogram_uParl_proj = histogram_uParl_vs_eta->ProjectionY(histogramName_uParl_proj.Data(), etaBin, etaBin, "e");
    // CV: skip (qT, eta) bins with limited event statistics
    //if ( !(histogram_uParl_proj->GetEntries() >= 100) ) continue;

    double y = 0.;
    double yErr = 0.;
    if ( mean_or_rms == kMean ) {
      y = -histogram_uParl_proj->GetMean()/qTav;
      yErr = histogram_uParl_proj->GetMeanError()/qTav;
    } else if ( mean_or_rms == kRMS ) {
      y = histogram_uParl_proj->GetRMS()/qTav;
      yErr = histogram_uParl_proj->GetRMSError()/qTav;
    } else assert (0);
    if ( mode == mode_uParl ) {
      y /= qTav;
      yErr /= qTav;
    } else if ( mode != mode_uParl_div_qT ) assert(0);
    //if ( mean_or_rms == kMean ) std::cout << " -uParl/qT = " << y << " +/- " << yErr << " (#entries = " << histogram_uParl_proj->GetEntries() << ")" << std::endl;
    //else if ( mean_or_rms == kRMS ) std::cout << " sigma(-uParl/qT) = " << y << " +/- " << yErr << " (#entries = " << histogram_uParl_proj->GetEntries() << ")" << std::endl;
    if ( divideByBinWidth ) {
      y /= (etaMax - etaMin);
      yErr /= (etaMax - etaMin);
    }

    graph->SetPoint(etaBin - 1, x, y);
    graph->SetPointError(etaBin - 1, xErrDown, xErrUp, yErr, yErr);
  }

  // reset x-axis range selection 
  histogram_qT->GetXaxis()->SetRange(1., 0.);

  return graph;
}

TGraphAsymmErrors* makeGraph_uParl_div_qT(const std::string& name, const std::string& title, 
					  const TH2* histogram_uParl_vs_qT, const TH1* histogram_qT, int mode, bool isCaloMEt)
{
  //std::cout << "<makeGraph_uParl_div_qT>:" << std::endl;
  //std::cout << " name = " << name << std::endl;

  if ( !(histogram_uParl_vs_qT && histogram_qT) ) return 0;
  
  const TAxis* xAxis = histogram_uParl_vs_qT->GetXaxis();
  int numBins = xAxis->GetNbins();

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(numBins);
  graph->SetName(name.data());
  graph->SetTitle(title.data());  

  for ( int iBin = 1; iBin <= numBins; ++iBin ) {
    double qTmin = xAxis->GetBinLowEdge(iBin);
    double qTmax = xAxis->GetBinUpEdge(iBin);

    int binLowIndex = const_cast<TH1*>(histogram_qT)->FindBin(qTmin);
    int binUpIndex  = const_cast<TH1*>(histogram_qT)->FindBin(qTmax);
    histogram_qT->GetXaxis()->SetRange(binLowIndex, binUpIndex);

    //double x        = histogram_qT->GetMean();
    double x        = xAxis->GetBinCenter(iBin);
    double xErrUp   = qTmax - x;
    double xErrDown = x - qTmin;

    TString histogramName_uParl_proj = Form("%s_py_%i", histogram_uParl_vs_qT->GetName(), iBin);
    TH1D* histogram_uParl_proj = histogram_uParl_vs_qT->ProjectionY(histogramName_uParl_proj.Data(), iBin, iBin, "e");
    // CV: skip qT bins with limited event statistics
    //if ( !(histogram_uParl_proj->GetEntries() >= 100) ) continue;

    if ( x > 0. ) {
      double y = -histogram_uParl_proj->GetMean();
      double yErr = histogram_uParl_proj->GetMeanError();      
      if ( mode == mode_uParl ) {
	y /= x;
	yErr /= x;
      } else if ( mode != mode_uParl_div_qT ) assert(0);
      if ( isCaloMEt ) y -= 1.; 
      graph->SetPoint(iBin - 1, x, y);      
      graph->SetPointError(iBin - 1, xErrDown, xErrUp, yErr, yErr);
    }
  }
  
  return graph;
}

//-------------------------------------------------------------------------------
TGraphAsymmErrors* compRatioGraph(const std::string& ratioGraphName, const TGraph* numerator, const TGraph* denominator, bool assertBinningCompatibility)
{
  std::cout << "<compRatioGraph>:" << std::endl;
  std::cout << " ratioGraphName = " << ratioGraphName << std::endl;

  assert(numerator->GetN() == denominator->GetN());
  int nPoints = numerator->GetN();

  TGraphAsymmErrors* graphRatio = new TGraphAsymmErrors(nPoints);
  graphRatio->SetName(ratioGraphName.data());

  for ( int iPoint = 0; iPoint < nPoints; ++iPoint ){
    double x_numerator, y_numerator;
    numerator->GetPoint(iPoint, x_numerator, y_numerator);
    double xErrUp_numerator = 0.;
    double xErrDown_numerator = 0.;
    double yErrUp_numerator = 0.;
    double yErrDown_numerator = 0.;
    if ( dynamic_cast<const TGraphAsymmErrors*>(numerator) ) {
      const TGraphAsymmErrors* numerator_asymmerrors = dynamic_cast<const TGraphAsymmErrors*>(numerator);
      xErrUp_numerator = numerator_asymmerrors->GetErrorXhigh(iPoint);
      xErrDown_numerator = numerator_asymmerrors->GetErrorXlow(iPoint);
      yErrUp_numerator = numerator_asymmerrors->GetErrorYhigh(iPoint);
      yErrDown_numerator = numerator_asymmerrors->GetErrorYlow(iPoint);
    } else if ( dynamic_cast<const TGraphErrors*>(numerator) ) {
      const TGraphErrors* numerator_errors = dynamic_cast<const TGraphErrors*>(numerator);
      xErrUp_numerator = numerator_errors->GetErrorX(iPoint);
      xErrDown_numerator = xErrUp_numerator;
      yErrUp_numerator = numerator_errors->GetErrorY(iPoint);
      yErrDown_numerator = yErrUp_numerator;
    }

    double x_denominator, y_denominator;
    denominator->GetPoint(iPoint, x_denominator, y_denominator);
    //std::cout << "point #" << iPoint << ": x(numerator) = " << x_numerator << ", x(denominator) = " << x_denominator << std::endl;
    if ( assertBinningCompatibility ) {
      if ( TMath::Abs(x_denominator - x_numerator) > 1.e-3 ) {
	std::cerr << "Incompatible binning of numerator and denominator graphs !!" << std::endl;
	std::cout << "point #" << iPoint << ": x(numerator) = " << x_numerator << ", x(denominator) = " << x_denominator << std::endl;
	assert(0);
      } 
    } else {
      x_denominator = x_numerator;
      y_denominator = denominator->Eval(x_denominator);
    }
    double xErrUp_denominator = 0.;
    double xErrDown_denominator = 0.;
    double yErrUp_denominator = 0.;
    double yErrDown_denominator = 0.;
    if ( dynamic_cast<const TGraphAsymmErrors*>(denominator) ) {
      const TGraphAsymmErrors* denominator_asymmerrors = dynamic_cast<const TGraphAsymmErrors*>(denominator);
      xErrUp_denominator = denominator_asymmerrors->GetErrorXhigh(iPoint);
      xErrDown_denominator = denominator_asymmerrors->GetErrorXlow(iPoint);
      yErrUp_denominator = denominator_asymmerrors->GetErrorYhigh(iPoint);
      yErrDown_denominator = denominator_asymmerrors->GetErrorYlow(iPoint);
    } else if ( dynamic_cast<const TGraphErrors*>(denominator) ) {
      const TGraphErrors* denominator_errors = dynamic_cast<const TGraphErrors*>(denominator);
      xErrUp_denominator = denominator_errors->GetErrorX(iPoint);
      xErrDown_denominator = xErrUp_denominator;
      yErrUp_denominator = denominator_errors->GetErrorY(iPoint);
      yErrDown_denominator = yErrUp_denominator;
    }

    double x_ratio = x_numerator;
    double y_ratio = ( y_denominator != 0. ) ? (y_numerator/y_denominator) : 0.;
    double xErrUp_ratio = TMath::Max(xErrUp_numerator, xErrUp_denominator);
    double xErrDown_ratio = TMath::Max(xErrDown_numerator, xErrDown_denominator);
    double yErr2Up_ratio = 0.;
    if ( y_numerator   ) yErr2Up_ratio += square(yErrUp_numerator/y_numerator);
    if ( y_denominator ) yErr2Up_ratio += square(yErrDown_denominator/y_numerator);
    double yErrUp_ratio = TMath::Sqrt(yErr2Up_ratio)*y_ratio;
    double yErr2Down_ratio = 0.;
    if ( y_numerator   ) yErr2Down_ratio += square(yErrDown_numerator/y_numerator);
    if ( y_denominator ) yErr2Down_ratio += square(yErrUp_denominator/y_numerator);
    double yErrDown_ratio = TMath::Sqrt(yErr2Down_ratio)*y_ratio;

    //std::cout << "point #" << iPoint << ": x = " << x_numerator << std::endl;
    //std::cout << " numerator: y = " << y_numerator << " + " << yErrUp_numerator << " - " << yErrDown_numerator << std::endl;
    //std::cout << " denominator: y = " << y_denominator << " + " << yErrUp_denominator << " - " << yErrDown_denominator << std::endl;
    //std::cout << "--> ratio: y = " << y_ratio << " + " << yErrUp_ratio << " - " << yErrDown_ratio << std::endl;

    graphRatio->SetPoint(iPoint, x_ratio, y_ratio);
    graphRatio->SetPointError(iPoint, xErrDown_ratio, xErrUp_ratio, yErrDown_ratio, yErrUp_ratio);
  }
  
  graphRatio->SetLineColor(numerator->GetLineColor());
  graphRatio->SetLineWidth(numerator->GetLineWidth());
  graphRatio->SetMarkerColor(numerator->GetMarkerColor());
  graphRatio->SetMarkerStyle(numerator->GetMarkerStyle());

  return graphRatio;
}

struct graph_with_Fit
{
  graph_with_Fit(TGraphAsymmErrors* graph, TF1* fit, const std::string& legendEntry)
    : graph_(graph),
      fit_(fit),
      legendEntry_(legendEntry)
  {}
  ~graph_with_Fit() {}
  std::string legendEntry_graph() const { return legendEntry_; }
  std::string legendEntry_fit() const { return std::string("Fit").append(" ").append(legendEntry_); }
  TGraphAsymmErrors* graph_;
  TF1* fit_;
  std::string legendEntry_;
};

void showGraphs_with_Fits(double canvasSizeX, double canvasSizeY,
			  const graph_with_Fit* graph_with_Fit_ref,
			  const graph_with_Fit* graph_with_Fit2,
			  const graph_with_Fit* graph_with_Fit3,
			  const graph_with_Fit* graph_with_Fit4,
			  const graph_with_Fit* graph_with_Fit_ratio1,
			  const graph_with_Fit* graph_with_Fit_ratio2,
			  const graph_with_Fit* graph_with_Fit_ratio3,
			  double xMin, double xMax, unsigned numBinsX, const std::string& xAxisTitle, double xAxisOffset,
			  bool useLogScale, double yMin, double yMax, double yMin_ratio, double yMax_ratio, const std::string& yAxisTitle, double yAxisOffset,
			  double legendX0, double legendY0, 
			  const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
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
  topPad->SetLogy(useLogScale);

  TPad* bottomPad = new TPad("bottomPad", "bottomPad", 0.00, 0.00, 1.00, 0.35);
  bottomPad->SetFillColor(10);
  bottomPad->SetTopMargin(0.02);
  bottomPad->SetLeftMargin(0.15);
  bottomPad->SetBottomMargin(0.24);
  bottomPad->SetRightMargin(0.05);
  bottomPad->SetLogy(false);

  canvas->cd();
  topPad->Draw();
  topPad->cd();

  int colors[7] = { 1, 2, 3, 4, 6, 7, 8 };
  int markerStyles[7] = { 20, 21, 22, 23, 25, 26, 32 };

  TLegend* legend = new TLegend(legendX0, legendY0, legendX0 + 0.43, legendY0 + 0.25, "", "brNDC"); 
  legend->SetBorderSize(0);
  legend->SetFillColor(0);

  TH1* dummyHistogram_top = new TH1D("dummyHistogram_top", "dummyHistogram_top", numBinsX, xMin, xMax);
  dummyHistogram_top->SetTitle("");
  dummyHistogram_top->SetStats(false);
  dummyHistogram_top->SetMinimum(yMin);
  dummyHistogram_top->SetMaximum(yMax);

  TAxis* xAxis_top = dummyHistogram_top->GetXaxis();
  xAxis_top->SetTitle(xAxisTitle.data());
  xAxis_top->SetTitleOffset(xAxisOffset);
  xAxis_top->SetLabelColor(10);
  xAxis_top->SetTitleColor(10);

  TAxis* yAxis_top = dummyHistogram_top->GetYaxis();
  yAxis_top->SetTitle(yAxisTitle.data());
  yAxis_top->SetTitleSize(0.045);
  yAxis_top->SetTitleOffset(yAxisOffset);

  dummyHistogram_top->Draw("axis");

  graph_with_Fit_ref->graph_->SetLineColor(colors[0]);
  graph_with_Fit_ref->graph_->SetLineWidth(1);
  graph_with_Fit_ref->graph_->SetMarkerColor(colors[0]);
  graph_with_Fit_ref->graph_->SetMarkerStyle(markerStyles[0]);
  graph_with_Fit_ref->graph_->SetMarkerSize(1);
  graph_with_Fit_ref->graph_->Draw("p");
  legend->AddEntry(graph_with_Fit_ref->graph_, graph_with_Fit_ref->legendEntry_graph().data(), "p");
  if ( graph_with_Fit_ref->fit_ ) {
    graph_with_Fit_ref->fit_->SetLineColor(colors[0]);
    graph_with_Fit_ref->fit_->SetLineWidth(1);
    graph_with_Fit_ref->fit_->Draw("same");
    legend->AddEntry(graph_with_Fit_ref->fit_, graph_with_Fit_ref->legendEntry_fit().data(), "l");
  }

  if ( graph_with_Fit2 && graph_with_Fit2->graph_ ) {
    graph_with_Fit2->graph_->SetLineColor(colors[1]);
    graph_with_Fit2->graph_->SetLineWidth(1);
    graph_with_Fit2->graph_->SetMarkerColor(colors[1]);
    graph_with_Fit2->graph_->SetMarkerStyle(markerStyles[1]);
    graph_with_Fit2->graph_->SetMarkerSize(1);
    graph_with_Fit2->graph_->Draw("p");
    legend->AddEntry(graph_with_Fit2->graph_, graph_with_Fit2->legendEntry_graph().data(), "p");
    if ( graph_with_Fit2->fit_ ) {
      graph_with_Fit2->fit_->SetLineColor(colors[1]);
      graph_with_Fit2->fit_->SetLineWidth(1);
      graph_with_Fit2->fit_->Draw("same");
      legend->AddEntry(graph_with_Fit2->fit_, graph_with_Fit2->legendEntry_fit().data(), "l");
    }
  }

  if ( graph_with_Fit3 && graph_with_Fit3->graph_ ) {
    graph_with_Fit3->graph_->SetLineColor(colors[2]);
    graph_with_Fit3->graph_->SetLineWidth(1);
    graph_with_Fit3->graph_->SetMarkerColor(colors[2]);
    graph_with_Fit3->graph_->SetMarkerStyle(markerStyles[2]);
    graph_with_Fit3->graph_->SetMarkerSize(1);
    graph_with_Fit3->graph_->Draw("p");
    legend->AddEntry(graph_with_Fit3->graph_, graph_with_Fit3->legendEntry_graph().data(), "p");
    if ( graph_with_Fit3->fit_ ) {
      graph_with_Fit3->fit_->SetLineColor(colors[2]);
      graph_with_Fit3->fit_->SetLineWidth(1);
      graph_with_Fit3->fit_->Draw("same");
      legend->AddEntry(graph_with_Fit3->fit_, graph_with_Fit3->legendEntry_fit().data(), "l");
    }
  }

  if ( graph_with_Fit4 && graph_with_Fit4->graph_ ) {
    graph_with_Fit4->graph_->SetLineColor(colors[3]);
    graph_with_Fit4->graph_->SetLineWidth(1);
    graph_with_Fit4->graph_->SetMarkerColor(colors[3]);
    graph_with_Fit4->graph_->SetMarkerStyle(markerStyles[3]);
    graph_with_Fit4->graph_->SetMarkerSize(1);
    graph_with_Fit4->graph_->Draw("p");
    legend->AddEntry(graph_with_Fit4->graph_, graph_with_Fit4->legendEntry_graph().data(), "p");
    if ( graph_with_Fit4->fit_ ) {
      graph_with_Fit4->fit_->SetLineColor(colors[3]);
      graph_with_Fit4->fit_->SetLineWidth(1);
      graph_with_Fit4->fit_->Draw("same");
      legend->AddEntry(graph_with_Fit4->fit_, graph_with_Fit4->legendEntry_fit().data(), "l");
    }
  }
  
  if ( graph_with_Fit_ratio1 && graph_with_Fit_ratio1->graph_ ) {
    graph_with_Fit_ratio1->graph_->SetLineColor(colors[4]);
    graph_with_Fit_ratio1->graph_->SetLineWidth(1);
    graph_with_Fit_ratio1->graph_->SetMarkerColor(colors[4]);
    graph_with_Fit_ratio1->graph_->SetMarkerStyle(markerStyles[4]);
    graph_with_Fit_ratio1->graph_->SetMarkerSize(1);
    legend->AddEntry(graph_with_Fit_ratio1->graph_, graph_with_Fit_ratio1->legendEntry_graph().data(), "p");
    if ( graph_with_Fit_ratio1->fit_ ) {
      graph_with_Fit_ratio1->fit_->SetLineColor(colors[4]);
      graph_with_Fit_ratio1->fit_->SetLineWidth(1);
      legend->AddEntry(graph_with_Fit_ratio1->fit_, graph_with_Fit_ratio1->legendEntry_fit().data(), "l");
    }
  }

  if ( graph_with_Fit_ratio2 && graph_with_Fit_ratio2->graph_ ) {
    graph_with_Fit_ratio2->graph_->SetLineColor(colors[5]);
    graph_with_Fit_ratio2->graph_->SetLineWidth(1);
    graph_with_Fit_ratio2->graph_->SetMarkerColor(colors[5]);
    graph_with_Fit_ratio2->graph_->SetMarkerStyle(markerStyles[5]);
    graph_with_Fit_ratio2->graph_->SetMarkerSize(1);
    legend->AddEntry(graph_with_Fit_ratio2->graph_, graph_with_Fit_ratio2->legendEntry_graph().data(), "p");
    if ( graph_with_Fit_ratio2->fit_ ) {
      graph_with_Fit_ratio2->fit_->SetLineColor(colors[5]);
      graph_with_Fit_ratio2->fit_->SetLineWidth(1);
      legend->AddEntry(graph_with_Fit_ratio2->fit_, graph_with_Fit_ratio2->legendEntry_fit().data(), "l");
    }
  }
  
  if ( graph_with_Fit_ratio3 && graph_with_Fit_ratio3->graph_ ) {
    graph_with_Fit_ratio3->graph_->SetLineColor(colors[6]);
    graph_with_Fit_ratio3->graph_->SetLineWidth(1);
    graph_with_Fit_ratio3->graph_->SetMarkerColor(colors[6]);
    graph_with_Fit_ratio3->graph_->SetMarkerStyle(markerStyles[6]);
    graph_with_Fit_ratio3->graph_->SetMarkerSize(1);
    legend->AddEntry(graph_with_Fit_ratio3->graph_, graph_with_Fit_ratio3->legendEntry_graph().data(), "p");
    if ( graph_with_Fit_ratio3->fit_ ) {
      graph_with_Fit_ratio3->fit_->SetLineColor(colors[6]);
      graph_with_Fit_ratio3->fit_->SetLineWidth(1);
      legend->AddEntry(graph_with_Fit_ratio3->fit_, graph_with_Fit_ratio3->legendEntry_fit().data(), "l");
    }
  }
 
  legend->Draw();

  canvas->cd();
  bottomPad->Draw();
  bottomPad->cd();

  TH1* dummyHistogram_bottom = new TH1D("dummyHistogram_bottom", "dummyHistogram_bottom", numBinsX, xMin, xMax);
  dummyHistogram_bottom->SetTitle("");
  dummyHistogram_bottom->SetStats(false);
  dummyHistogram_bottom->SetMinimum(yMin_ratio);
  dummyHistogram_bottom->SetMaximum(yMax_ratio);

  TAxis* xAxis_bottom = dummyHistogram_bottom->GetXaxis();
  xAxis_bottom->SetTitle(xAxis_top->GetTitle());
  xAxis_bottom->SetLabelColor(1);
  xAxis_bottom->SetTitleColor(1);
  xAxis_bottom->SetTitleOffset(1.20);
  xAxis_bottom->SetTitleSize(0.08);
  xAxis_bottom->SetLabelOffset(0.02);
  xAxis_bottom->SetLabelSize(0.08);
  xAxis_bottom->SetTickLength(0.055);
  
  TAxis* yAxis_bottom = dummyHistogram_bottom->GetYaxis();
  yAxis_bottom->SetTitle("#frac{Simulation}{Data}");
  yAxis_bottom->SetTitleOffset(0.70);
  yAxis_bottom->SetNdivisions(505);
  yAxis_bottom->CenterTitle();
  yAxis_bottom->SetTitleSize(0.08);
  yAxis_bottom->SetLabelSize(0.08);
  yAxis_bottom->SetTickLength(0.04); 

  dummyHistogram_bottom->Draw("axis");
  
  TGraph* graph2_div_ref = 0;
  if ( graph_with_Fit2 && graph_with_Fit2->graph_ ) {
    std::string graphName2_div_ref = std::string(graph_with_Fit2->graph_->GetName()).append("_div_").append(graph_with_Fit_ref->graph_->GetName());
    graph2_div_ref = compRatioGraph(graphName2_div_ref, graph_with_Fit2->graph_, graph_with_Fit_ref->graph_, false);
    graph2_div_ref->Draw("p");
  }

  TGraph* graph3_div_ref = 0;
  if ( graph_with_Fit3 && graph_with_Fit3->graph_ ) {
    std::string graphName3_div_ref = std::string(graph_with_Fit3->graph_->GetName()).append("_div_").append(graph_with_Fit_ref->graph_->GetName());
    graph3_div_ref = compRatioGraph(graphName3_div_ref, graph_with_Fit3->graph_, graph_with_Fit_ref->graph_, false);
    graph3_div_ref->Draw("p");
  }

  TGraph* graph4_div_ref = 0;
  if ( graph_with_Fit4 && graph_with_Fit4->graph_ ) {
    std::string graphName4_div_ref = std::string(graph_with_Fit4->graph_->GetName()).append("_div_").append(graph_with_Fit_ref->graph_->GetName());
    graph4_div_ref = compRatioGraph(graphName4_div_ref, graph_with_Fit4->graph_, graph_with_Fit_ref->graph_, false);
    graph4_div_ref->Draw("p");
  }

  if ( graph_with_Fit_ratio1 && graph_with_Fit_ratio1->graph_ ) {
    graph_with_Fit_ratio1->graph_->Draw("p");
    if ( graph_with_Fit_ratio1->fit_ ) {
      graph_with_Fit_ratio1->fit_->Draw("same");
    }
  }

  if ( graph_with_Fit_ratio2 && graph_with_Fit_ratio2->graph_ ) {
    graph_with_Fit_ratio2->graph_->Draw("p");
    if ( graph_with_Fit_ratio2->fit_ ) {
      graph_with_Fit_ratio2->fit_->Draw("same");
    }
  }

  if ( graph_with_Fit_ratio3 && graph_with_Fit_ratio3->graph_ ) {
    graph_with_Fit_ratio3->graph_->Draw("p");
    if ( graph_with_Fit_ratio3->fit_ ) {
      graph_with_Fit_ratio3->fit_->Draw("same");
    }
  }
     
  canvas->Update();
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  if ( useLogScale ) outputFileName_plot.append("_log");
  else outputFileName_plot.append("_linear");
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  canvas->Print(std::string(outputFileName_plot).append(".root").data());
  
  delete legend;
  delete graph2_div_ref;
  delete graph3_div_ref;
  delete graph4_div_ref;
  delete dummyHistogram_top;
  delete topPad;
  delete dummyHistogram_bottom;
  delete bottomPad;
  delete canvas;  

  delete graph_with_Fit_ref;
  delete graph_with_Fit2;
  delete graph_with_Fit3;
  delete graph_with_Fit4;
  delete graph_with_Fit_ratio1;
  delete graph_with_Fit_ratio2;
  delete graph_with_Fit_ratio3;
}

void showGraphs(double canvasSizeX, double canvasSizeY,
		TGraphAsymmErrors* graph_ref, const std::string& legendEntry_ref,
		TGraphAsymmErrors* graph2, const std::string& legendEntry2,
		TGraphAsymmErrors* graph3, const std::string& legendEntry3,
		TGraphAsymmErrors* graph4, const std::string& legendEntry4,
		TGraphAsymmErrors* graphRatio1, const std::string& legendEntryRatio1,
		TGraphAsymmErrors* graphRatio2, const std::string& legendEntryRatio2,
		TGraphAsymmErrors* graphRatio3, const std::string& legendEntryRatio3,
		double xMin, double xMax, unsigned numBinsX, const std::string& xAxisTitle, double xAxisOffset,
		bool useLogScale, double yMin, double yMax, double yMin_ratio, double yMax_ratio, const std::string& yAxisTitle, double yAxisOffset,
		double legendX0, double legendY0, 
		const std::string& outputFileName)
{
  showGraphs_with_Fits(canvasSizeX, canvasSizeY,
		       new graph_with_Fit(graph_ref, 0, legendEntry_ref),
		       new graph_with_Fit(graph2, 0, legendEntry2),
		       new graph_with_Fit(graph3, 0, legendEntry3),
		       new graph_with_Fit(graph4, 0, legendEntry4),
		       new graph_with_Fit(graphRatio1, 0, legendEntryRatio1),
		       new graph_with_Fit(graphRatio2, 0, legendEntryRatio2),
		       new graph_with_Fit(graphRatio3, 0, legendEntryRatio3),
		       xMin, xMax, numBinsX, xAxisTitle, xAxisOffset,
		       useLogScale, yMin, yMax, yMin_ratio, yMax_ratio, yAxisTitle, yAxisOffset,
		       legendX0, legendY0, 
		       outputFileName);
}
//-------------------------------------------------------------------------------

//-------------------------------------------------------------------------------
TGraphAsymmErrors* sumGraphs(const std::vector<const TGraphAsymmErrors*>& graphs)
{
  const TGraphAsymmErrors* graphRef = graphs[0];
  int numPoints = graphRef->GetN();
  TGraphAsymmErrors* graphSum = new TGraphAsymmErrors(numPoints);
  for ( int iPoint = 0; iPoint < numPoints; ++iPoint ) {
    double xRef, dummy;
    graphRef->GetPoint(iPoint, xRef, dummy);
    double xErrLow = graphRef->GetErrorXlow(iPoint);
    double xErrHigh = graphRef->GetErrorXhigh(iPoint);
    double y_sum = 0.;
    double yErrLow2_sum = 0.;
    double yErrHigh2_sum = 0.;
    for ( std::vector<const TGraphAsymmErrors*>::const_iterator graph = graphs.begin();
	  graph != graphs.end(); ++graph ) {
      double x, y;
      (*graph)->GetPoint(iPoint, x, y);
      y_sum += y;
      double yErrLow = (*graph)->GetErrorYlow(iPoint);
      double yErrHigh = (*graph)->GetErrorYhigh(iPoint);
      yErrLow2_sum += (yErrLow*yErrLow);
      yErrHigh2_sum += (yErrHigh*yErrHigh);
    }
    graphSum->SetPoint(iPoint, xRef, y_sum);
    graphSum->SetPointError(iPoint, xErrLow, xErrHigh, TMath::Sqrt(yErrLow2_sum), TMath::Sqrt(yErrHigh2_sum));
  }
  return graphSum;
}

TGraphAsymmErrors* sumGraphs(const TGraphAsymmErrors* graph1, const TGraphAsymmErrors* graph2, const TGraphAsymmErrors* graph3 = 0, const TGraphAsymmErrors* graph4 = 0)
{
  std::vector<const TGraphAsymmErrors*> graphs;
  graphs.push_back(graph1);
  graphs.push_back(graph2);
  if ( graph3 ) graphs.push_back(graph3);
  if ( graph4 ) graphs.push_back(graph4);
  return sumGraphs(graphs);
}

TH1* convertToHistogram(TGraphAsymmErrors* graph)
{
  int numBins = graph->GetN();
  TArrayF binning(numBins + 1);
  for ( int iBin = 0; iBin < numBins; ++iBin ) {
    double x, y;
    graph->GetPoint(iBin, x, y);

    double binEdgeLow  = x - graph->GetErrorXlow(iBin);
    double binEdgeHigh = x + graph->GetErrorXhigh(iBin);

    binning[iBin] = binEdgeLow;
    binning[iBin + 1] = binEdgeHigh;
  }

  std::string histogramName = std::string(graph->GetName()).append("_histogram");
  TH1F* histogram = new TH1F(histogramName.data(), histogramName.data(), numBins, binning.GetArray());
  for ( int iBin = 0; iBin < numBins; ++iBin ) {
    double x, y;
    graph->GetPoint(iBin, x, y);

    double yErrLow = graph->GetErrorYlow(iBin);
    double yErrHigh = graph->GetErrorYhigh(iBin);
    double yErr = TMath::Sqrt(0.5*(yErrLow*yErrLow + yErrHigh*yErrHigh));

    histogram->SetBinContent(iBin + 1, y);
    histogram->SetBinError(iBin + 1, yErr);
  }

  return histogram;
}

TGraphAsymmErrors* fixJetsResponse(TGraphAsymmErrors* graph)
{
  int numPoints = graph->GetN();
  TGraphAsymmErrors* graph_fixed = new TGraphAsymmErrors(numPoints);
  for ( int iPoint = 0; iPoint < numPoints; ++iPoint ) {
    double x, y;
    graph->GetPoint(iPoint, x, y);
    double xErrUp   = graph->GetErrorXhigh(iPoint);
    double xErrDown = graph->GetErrorXlow(iPoint);
    double yErrUp   = graph->GetErrorYhigh(iPoint);
    double yErrDown = graph->GetErrorYlow(iPoint);
    graph_fixed->SetPoint(iPoint, x, -y + 1.);
    graph_fixed->SetPointError(iPoint, xErrDown, xErrUp, yErrDown, yErrUp);
  }
  return graph_fixed;
}

TGraphAsymmErrors* fixResponse(TGraphAsymmErrors* graph)
{
  int numPoints = graph->GetN();
  TGraphAsymmErrors* graph_fixed = new TGraphAsymmErrors(numPoints);
  for ( int iPoint = 0; iPoint < numPoints; ++iPoint ) {
    double x, y;
    graph->GetPoint(iPoint, x, y);
    double xErrUp   = graph->GetErrorXhigh(iPoint);
    double xErrDown = graph->GetErrorXlow(iPoint);
    double yErrUp   = graph->GetErrorYhigh(iPoint);
    double yErrDown = graph->GetErrorYlow(iPoint);
    graph_fixed->SetPoint(iPoint, x, -y);
    graph_fixed->SetPointError(iPoint, xErrDown, xErrUp, yErrDown, yErrUp);
  }
  return graph_fixed;
}

void showResponse_in_Components(double canvasSizeX, double canvasSizeY,
				TGraphAsymmErrors* graph_unclustered_mc, TGraphAsymmErrors* graph_unclustered_data, 
				TGraphAsymmErrors* graph_unclCorr_mc, TGraphAsymmErrors* graph_unclCorr_data, 
				TGraphAsymmErrors* graph_jets_mc, TGraphAsymmErrors* graph_jets_data, 
				TGraphAsymmErrors* graph_jetCorr_mc, TGraphAsymmErrors* graph_jetCorr_data, 
				bool makeCumulative,
				double xMin, double xMax, unsigned numBinsX, const std::string& xAxisTitle, double xAxisOffset,
				bool useLogScale, double yMin, double yMax, double yMin_ratio, double yMax_ratio, const std::string& yAxisTitle, double yAxisOffset,
				double legendX0, double legendY0, 
				const std::string& outputFileName)
{
  //graph_unclustered_mc = fixResponse(graph_unclustered_mc);
  //graph_unclCorr_mc = fixResponse(graph_unclCorr_mc);
  //graph_jets_mc = fixJetsResponse(graph_jets_mc);
  //graph_jetCorr_mc = fixResponse(graph_jetCorr_mc);

  //graph_unclustered_data = fixResponse(graph_unclustered_data);
  //graph_unclCorr_data = fixResponse(graph_unclCorr_data);
  //graph_jets_data = fixJetsResponse(graph_jets_data);
  //graph_jetCorr_data = fixResponse(graph_jetCorr_data);

  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
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
  topPad->SetLogy(useLogScale);

  TPad* bottomPad = new TPad("bottomPad", "bottomPad", 0.00, 0.00, 1.00, 0.35);
  bottomPad->SetFillColor(10);
  bottomPad->SetTopMargin(0.02);
  bottomPad->SetLeftMargin(0.15);
  bottomPad->SetBottomMargin(0.24);
  bottomPad->SetRightMargin(0.05);
  bottomPad->SetLogy(false);

  canvas->cd();
  topPad->Draw();
  topPad->cd();

  int colors[7] = { 1, 2, 3, 4, 6, 7, 8 };
  int markerStyles[7] = { 20, 21, 22, 23, 25, 26, 32 };

  TLegend* legend = new TLegend(legendX0, legendY0, legendX0 + 0.43, legendY0 + 0.25, "", "brNDC"); 
  legend->SetBorderSize(0);
  legend->SetFillColor(0);

  TH1* dummyHistogram_top = new TH1D("dummyHistogram_top", "dummyHistogram_top", numBinsX, xMin, xMax);
  dummyHistogram_top->SetTitle("");
  dummyHistogram_top->SetStats(false);
  dummyHistogram_top->SetMinimum(yMin);
  dummyHistogram_top->SetMaximum(yMax);

  TAxis* xAxis_top = dummyHistogram_top->GetXaxis();
  xAxis_top->SetTitle(xAxisTitle.data());
  xAxis_top->SetTitleOffset(xAxisOffset);
  xAxis_top->SetLabelColor(10);
  xAxis_top->SetTitleColor(10);

  TAxis* yAxis_top = dummyHistogram_top->GetYaxis();
  yAxis_top->SetTitle(yAxisTitle.data());
  yAxis_top->SetTitleSize(0.045);
  yAxis_top->SetTitleOffset(yAxisOffset);

  dummyHistogram_top->Draw("axis");

  TGraphAsymmErrors* graph1_mc = graph_unclustered_mc;
  TGraphAsymmErrors* graph2_mc = 0;
  TGraphAsymmErrors* graph3_mc = 0;
  TGraphAsymmErrors* graph4_mc = 0;
  if ( makeCumulative ) {
    graph2_mc = sumGraphs(graph_jets_mc, graph_unclustered_mc);
    graph3_mc = sumGraphs(graph_jetCorr_mc, graph_unclustered_mc, graph_jets_mc);
    graph4_mc = sumGraphs(graph_unclCorr_mc, graph_unclustered_mc, graph_jets_mc, graph_jetCorr_mc);
  } else {
    graph2_mc = graph_jets_mc;
    graph3_mc = graph_jetCorr_mc;
    graph4_mc = graph_unclCorr_mc;
  }
  TH1* histogram1_mc = convertToHistogram(graph1_mc);
  TH1* histogram2_mc = convertToHistogram(graph2_mc);
  TH1* histogram3_mc = convertToHistogram(graph3_mc);
  TH1* histogram4_mc = convertToHistogram(graph4_mc);

  TGraphAsymmErrors* graph1_data = graph_unclustered_data;
  TGraphAsymmErrors* graph2_data = 0;
  TGraphAsymmErrors* graph3_data = 0;
  TGraphAsymmErrors* graph4_data = 0;
  if ( makeCumulative ) {
    graph2_data = sumGraphs(graph_jets_data, graph_unclustered_data);
    graph3_data = sumGraphs(graph_jetCorr_data, graph_unclustered_data, graph_jets_data);
    graph4_data = sumGraphs(graph_unclCorr_data, graph_unclustered_data, graph_jets_data, graph_jetCorr_data);
  } else {
    graph2_data = graph_jets_data;
    graph3_data = graph_jetCorr_data;
    graph4_data = graph_unclCorr_data;
  }

  histogram1_mc->SetLineColor(colors[0]);
  histogram1_mc->SetLineWidth(2);
  histogram1_mc->Draw("histsame");
  legend->AddEntry(histogram1_mc, "Unclustered MC", "l");
  graph1_data->SetLineColor(colors[0]);
  graph1_data->SetLineWidth(1);
  graph1_data->SetMarkerColor(colors[0]);
  graph1_data->SetMarkerStyle(markerStyles[0]);
  graph1_data->SetMarkerSize(1);
  graph1_data->Draw("p");
  legend->AddEntry(graph1_data, "Unclustered Data", "p");

  histogram2_mc->SetLineColor(colors[1]);
  histogram2_mc->SetLineWidth(2);
  histogram2_mc->Draw("histsame");
  legend->AddEntry(histogram2_mc, "Jets MC", "l");
  graph2_data->SetLineColor(colors[1]);
  graph2_data->SetLineWidth(1);
  graph2_data->SetMarkerColor(colors[1]);
  graph2_data->SetMarkerStyle(markerStyles[1]);
  graph2_data->SetMarkerSize(1);
  graph2_data->Draw("p");
  legend->AddEntry(graph2_data, "Jets Data", "p");

  histogram3_mc->SetLineColor(colors[2]);
  histogram3_mc->SetLineWidth(2);
  histogram3_mc->Draw("histsame");
  legend->AddEntry(histogram3_mc, "Jet Corr. MC", "l");
  graph3_data->SetLineColor(colors[2]);
  graph3_data->SetLineWidth(1);
  graph3_data->SetMarkerColor(colors[2]);
  graph3_data->SetMarkerStyle(markerStyles[2]);
  graph3_data->SetMarkerSize(1);
  graph3_data->Draw("p");
  legend->AddEntry(graph3_data, "Jet Corr. Data", "p");

  histogram4_mc->SetLineColor(colors[3]);
  histogram4_mc->SetLineWidth(2);
  histogram4_mc->Draw("histsame");
  legend->AddEntry(histogram4_mc, "Uncl. Corr MC", "l");
 
  legend->Draw();

  canvas->cd();
  bottomPad->Draw();
  bottomPad->cd();

  TH1* dummyHistogram_bottom = new TH1D("dummyHistogram_bottom", "dummyHistogram_bottom", numBinsX, xMin, xMax);
  dummyHistogram_bottom->SetTitle("");
  dummyHistogram_bottom->SetStats(false);
  dummyHistogram_bottom->SetMinimum(yMin_ratio);
  dummyHistogram_bottom->SetMaximum(yMax_ratio);

  TAxis* xAxis_bottom = dummyHistogram_bottom->GetXaxis();
  xAxis_bottom->SetTitle(xAxis_top->GetTitle());
  xAxis_bottom->SetLabelColor(1);
  xAxis_bottom->SetTitleColor(1);
  xAxis_bottom->SetTitleOffset(1.20);
  xAxis_bottom->SetTitleSize(0.08);
  xAxis_bottom->SetLabelOffset(0.02);
  xAxis_bottom->SetLabelSize(0.08);
  xAxis_bottom->SetTickLength(0.055);
  
  TAxis* yAxis_bottom = dummyHistogram_bottom->GetYaxis();
  yAxis_bottom->SetTitle("#frac{Simulation}{Data}");
  yAxis_bottom->SetTitleOffset(0.70);
  yAxis_bottom->SetNdivisions(505);
  yAxis_bottom->CenterTitle();
  yAxis_bottom->SetTitleSize(0.08);
  yAxis_bottom->SetLabelSize(0.08);
  yAxis_bottom->SetTickLength(0.04); 

  dummyHistogram_bottom->Draw("axis");
  
  std::string graphName1_mc_div_data = std::string(graph1_mc->GetName()).append("_div_").append(graph1_data->GetName());
  TGraph* graph1_mc_div_data = compRatioGraph(graphName1_mc_div_data, graph1_mc, graph1_data, false);
  graph1_mc_div_data->Draw("p");
  graph1_mc_div_data->SetLineColor(colors[0]);
  graph1_mc_div_data->SetLineWidth(1);
  graph1_mc_div_data->SetMarkerColor(colors[0]);
  graph1_mc_div_data->SetMarkerStyle(markerStyles[0]);
  graph1_mc_div_data->SetMarkerSize(1);

  std::string graphName2_mc_div_data = std::string(graph2_mc->GetName()).append("_div_").append(graph2_data->GetName());
  TGraph* graph2_mc_div_data = compRatioGraph(graphName2_mc_div_data, graph2_mc, graph2_data, false);
  graph2_mc_div_data->Draw("p");
  graph2_mc_div_data->SetLineColor(colors[1]);
  graph2_mc_div_data->SetLineWidth(1);
  graph2_mc_div_data->SetMarkerColor(colors[1]);
  graph2_mc_div_data->SetMarkerStyle(markerStyles[1]);
  graph2_mc_div_data->SetMarkerSize(1);

  std::string graphName3_mc_div_data = std::string(graph3_mc->GetName()).append("_div_").append(graph3_data->GetName());
  TGraph* graph3_mc_div_data = compRatioGraph(graphName3_mc_div_data, graph3_mc, graph3_data, false);
  graph3_mc_div_data->Draw("p");
  graph3_mc_div_data->SetLineColor(colors[2]);
  graph3_mc_div_data->SetLineWidth(1);
  graph3_mc_div_data->SetMarkerColor(colors[2]);
  graph3_mc_div_data->SetMarkerStyle(markerStyles[2]);
  graph3_mc_div_data->SetMarkerSize(1);

  TGraph* graph4_mc_div_data = 0;
  if ( makeCumulative ) {
    std::string graphName4_mc_div_data = std::string(graph4_mc->GetName()).append("_div_").append(graph4_data->GetName());
    graph4_mc_div_data = compRatioGraph(graphName4_mc_div_data, graph4_mc, graph4_data, false);
    graph4_mc_div_data->Draw("p");
    graph4_mc_div_data->SetLineColor(colors[3]);
    graph4_mc_div_data->SetLineWidth(1);
    graph4_mc_div_data->SetMarkerColor(colors[3]);
    graph4_mc_div_data->SetMarkerStyle(markerStyles[3]);
    graph4_mc_div_data->SetMarkerSize(1);
  }
     
  canvas->Update();
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  if ( useLogScale ) outputFileName_plot.append("_log");
  else outputFileName_plot.append("_linear");
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  canvas->Print(std::string(outputFileName_plot).append(".root").data());
  
  delete legend;
  if ( makeCumulative ) {
    delete graph1_data;
    delete graph2_data;
    delete graph3_data;
    delete graph4_data;
    delete graph1_mc;
    delete graph2_mc;
    delete graph3_mc;
    delete graph4_mc;
  }
  delete histogram1_mc;
  delete histogram2_mc;
  delete histogram3_mc;
  delete histogram4_mc;
  delete graph1_mc_div_data;
  delete graph2_mc_div_data;
  delete graph3_mc_div_data;
  delete graph4_mc_div_data;
  delete dummyHistogram_top;
  delete topPad;
  delete dummyHistogram_bottom;
  delete bottomPad;
  delete canvas;  
}
//-------------------------------------------------------------------------------

TH1* getHistogramProjY(TH2* histogram, int iBinX)
{
  std::string histogramNameProjY = Form("%sProjY_bin%i", histogram->GetName(), iBinX);
  TH1* histogramProjY = histogram->ProjectionY(histogramNameProjY.data(), iBinX, iBinX);
  return histogramProjY;
}

TGraphAsymmErrors* makeGraph_vs_qT(std::map<int, TGraphAsymmErrors*>& graphs_uParl_vs_eta, const std::vector<double>& qTbinning, int etaBin)
{
  int qTnumBins = qTbinning.size() - 1;

  TGraphAsymmErrors* graphs_uParl_vs_qT = new TGraphAsymmErrors(qTnumBins);

  for ( int qTbin = 0; qTbin < qTnumBins; ++qTbin ) {
    double qTmin = qTbinning[qTbin];
    double qTmax = qTbinning[qTbin + 1];

    double x = 0.5*(qTmin + qTmax);
    double xErr = 0.5*(qTmax - qTmin);

    TGraphAsymmErrors* graph_uParl_vs_eta_data = graphs_uParl_vs_eta[qTbin];
    if ( !graph_uParl_vs_eta_data ) return 0;
    assert(etaBin >= 0 && etaBin < graph_uParl_vs_eta_data->GetN());

    double dummy, y;
    graph_uParl_vs_eta_data->GetPoint(etaBin, dummy, y);
    double yErrUp   = graph_uParl_vs_eta_data->GetErrorYhigh(etaBin);
    double yErrDown = graph_uParl_vs_eta_data->GetErrorYlow(etaBin);
    graphs_uParl_vs_qT->SetPoint(qTbin, x, y);
    graphs_uParl_vs_qT->SetPointError(qTbin, xErr, xErr, yErrDown, yErrUp);
  }

  return graphs_uParl_vs_qT;
}

std::string dqmSubDirectory(const std::string& dqmDirectory, double numPileUpMin, double numPileUpMax, double qTmin, double qTmax)
{
  TString dqmSubDirectory_tstring = dqmDirectory.data();
  if ( !dqmSubDirectory_tstring.EndsWith("/") ) dqmSubDirectory_tstring.Append("/");
  dqmSubDirectory_tstring.Append(Form("numPileUp%1.1fto%1.1f_qT%1.1fto%1.1f", numPileUpMin, numPileUpMax, qTmin, qTmax));
  dqmSubDirectory_tstring.ReplaceAll(".", "_");
  return std::string(dqmSubDirectory_tstring.Data());
}

void copyGraphPoint(TGraphAsymmErrors* graph_source, int iPoint_source, 
		    TGraphAsymmErrors* graph_target, int iPoint_target, double x_target, double xErrLow_target, double xErrHigh_target)
{
  double x_source, y;
  graph_source->GetPoint(iPoint_source, x_source, y);

  double yErrLow = graph_source->GetErrorYlow(iPoint_source);
  double yErrHigh = graph_source->GetErrorYhigh(iPoint_source);

  graph_target->SetPoint(iPoint_target, x_target, y);
  graph_target->SetPointError(iPoint_target, xErrLow_target, xErrHigh_target, yErrLow, yErrHigh);
}

std::vector<std::string> getResidualCorrections(const TGraphAsymmErrors* graph, int etaNumBins, double* etaBinning)
{
  std::vector<std::string> residualCorrections;
  residualCorrections.push_back("{1         JetEta              1          JetPt               [0]     Correction     L2Relative}");

  for ( int etaBin = 0; etaBin < etaNumBins; ++etaBin ) {            
    double etaMin = etaBinning[etaBin];
    double etaMax = etaBinning[etaBin + 1];

    double x, y;
    graph->GetPoint(etaBin, x, y);
    residualCorrections.push_back(Form("%1.3f         %1.3f              3              3           3500        %1.5f        ", etaMin, etaMax, y));
  }

  return residualCorrections;
}

void saveResidualCorrections(const std::vector<std::string>& residualCorrections, const std::string& outputFileName)
{
  std::cout << "<saveResidualCorrectionsst>:" << std::endl;
  
  std::ofstream* outputFile = new std::ofstream(outputFileName.data(), std::ios::out);
  for ( std::vector<std::string>::const_iterator line = residualCorrections.begin();
	line != residualCorrections.end(); ++line ) {
    (*outputFile) << (*line) << std::endl;
  }
  delete outputFile;
}

void makeUnclusteredEnergyPlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  //std::string type_string = "caloTowers";
  std::string type_string = "pfCands";
  //std::string type_string = "pfJetsPlusCandsNotInJet";
  //std::string type_string = "tracks";
  //std::string type_string = "genParticles";
  
  TFile* inputFile = 0;
  std::map<std::string, std::string> dqmDirectories;
  int type = 0;
  double yMin_response = 0.;
  double type1JetPtThreshold = 0.;
  if ( type_string == "caloTowers" ) {
    std::string inputFilePath = "/data2/veelken/CMSSW_5_3_x/Ntuples/unclEnCalibration/plots/v9_07_2013May01/";
    std::string inputFileName = "UnclusteredEnergyAnalyzer_all_caloTowers.root";
    inputFile = new TFile(Form("%s/%s", inputFilePath.data(), inputFileName.data()));
    dqmDirectories["data"]                 = "Data_caloTowers_central";
    dqmDirectories["mc_central"]           = "ZplusJets_madgraph_caloTowers_central";
    dqmDirectories["mc_shiftUp"]           = "ZplusJets_madgraph_caloTowers_shiftUp";
    dqmDirectories["mc_shiftDown"]         = "ZplusJets_madgraph_caloTowers_shiftDown";
    type = type_Calo;
    yMin_response = 0.2;
  } else if ( type_string == "caloTowersNoHF" ) {
    std::string inputFilePath = "/data2/veelken/CMSSW_5_3_x/Ntuples/unclEnCalibration/plots/v9_07_2013May01/";
    std::string inputFileName = "UnclusteredEnergyAnalyzer_all_caloTowersNoHF.root";
    inputFile = new TFile(Form("%s/%s", inputFilePath.data(), inputFileName.data()));
    dqmDirectories["data"]                 = "Data_caloTowersNoHF_central";
    dqmDirectories["mc_central"]           = "ZplusJets_madgraph_caloTowersNoHF_central";
    dqmDirectories["mc_shiftUp"]           = "ZplusJets_madgraph_caloTowersNoHF_shiftUp";
    dqmDirectories["mc_shiftDown"]         = "ZplusJets_madgraph_caloTowersNoHF_shiftDown";
    type = type_Calo;
    yMin_response = 0.2;
  } else if ( type_string == "pfCands" ) {
    std::string inputFilePath = "/data2/veelken/CMSSW_5_3_x/Ntuples/unclEnCalibration/plots/v2_0_1/";
    std::string inputFileName = "UnclusteredEnergyAnalyzer_all_pfCands.root";
    inputFile = new TFile(Form("%s/%s", inputFilePath.data(), inputFileName.data()));
    dqmDirectories["data"]                 = "Data_pfCands_central";
    dqmDirectories["mc_central"]           = "ZplusJets_madgraph_pfCands_central";
    dqmDirectories["mc_shiftUp"]           = "ZplusJets_madgraph_pfCands_shiftUp";
    dqmDirectories["mc_shiftDown"]         = "ZplusJets_madgraph_pfCands_shiftDown";
    dqmDirectories["mc_wCalibration"]      = "ZplusJets_madgraph_pfCands_wCalibration";
    //dqmDirectories["mc_wCalibration"]      = dqmDirectories["mc_central"];
    type = type_PF;
    yMin_response = 0.5;
  } else if ( type_string == "pfJetsPlusCandsNotInJet" ) {
    std::string inputFilePath = "/data2/veelken/CMSSW_5_3_x/Ntuples/unclEnCalibration/plots/v2_0_1/";
    std::string inputFileName = "UnclusteredEnergyAnalyzer_all_pfJetsPlusCandsNotInJet.root";
    inputFile = new TFile(Form("%s/%s", inputFilePath.data(), inputFileName.data()));
    type1JetPtThreshold = 20.;
    dqmDirectories["data"]                 = TString(Form("Data_pfJetsPlusCandsNotInJet_central_jetPtThreshold%1.1f", type1JetPtThreshold)).ReplaceAll(".", "_").Data();
    dqmDirectories["mc_central"]           = TString(Form("ZplusJets_madgraph_pfJetsPlusCandsNotInJet_central_jetPtThreshold%1.1f", type1JetPtThreshold)).ReplaceAll(".", "_").Data();
    dqmDirectories["mc_shiftUp"]           = TString(Form("ZplusJets_madgraph_pfJetsPlusCandsNotInJet_shiftUp_jetPtThreshold%1.1f", type1JetPtThreshold)).ReplaceAll(".", "_").Data();
    dqmDirectories["mc_shiftDown"]         = TString(Form("ZplusJets_madgraph_pfJetsPlusCandsNotInJet_shiftDown_jetPtThreshold%1.1f", type1JetPtThreshold)).ReplaceAll(".", "_").Data();
    dqmDirectories["mc_wCalibration"]      = TString(Form("ZplusJets_madgraph_pfJetsPlusCandsNotInJet_wCalibration_jetPtThreshold%1.1f", type1JetPtThreshold)).ReplaceAll(".", "_").Data();
    type = type_PF;
    yMin_response = 0.5;
  } else if ( type_string == "tracks" ) {
    std::string inputFilePath = "/data2/veelken/CMSSW_5_3_x/Ntuples/unclEnCalibration/plots/v9_07_2013May01/";
    std::string inputFileName = "UnclusteredEnergyAnalyzer_all_tracks.root";
    inputFile = new TFile(Form("%s/%s", inputFilePath.data(), inputFileName.data()));
    dqmDirectories["data"]                 = "Data_tracks_central";
    dqmDirectories["mc_central"]           = "ZplusJets_madgraph_tracks_central";
    dqmDirectories["mc_shiftUp"]           = "";
    dqmDirectories["mc_shiftDown"]         = "";
    type = type_Track;
    yMin_response = 0.5;
  } else if ( type_string == "genParticles" ) {
    std::string inputFilePath = "/data2/veelken/CMSSW_5_3_x/Ntuples/unclEnCalibration/plots/v9_07_2013May01/";
    std::string inputFileName = "UnclusteredEnergyAnalyzer_all_pfCands.root";
    inputFile = new TFile(Form("%s/%s", inputFilePath.data(), inputFileName.data()));
    dqmDirectories["data"]                 = "";
    dqmDirectories["mc_central"]           = "ZplusJets_madgraph_genParticles_central";
    dqmDirectories["mc_shiftUp"]           = "";
    dqmDirectories["mc_shiftDown"]         = "";
    type = type_PF;
    yMin_response = 0.5;
  } else {
    std::cerr << "Invalid Configuration Parameter 'type' = " << type  << " !!" << std::endl;
    assert(0);
  }

  std::cout << "inputFile = " << inputFile->GetName() << std::endl;

  int mode = mode_uParl_div_qT;
  std::string histogramName_uParl_vs_qT;
  if ( mode == mode_uParl ) {
    histogramName_uParl_vs_qT = "uParl_vs_qT";
  } else if ( mode == mode_uParl_div_qT ) {
    histogramName_uParl_vs_qT = "uParlDivQt_vs_qT";
  } else assert(0);

  TH1* histogram_qT_data = getHistogram(inputFile, dqmDirectories["data"], "qT");
  TH1* histogram_qT_mc = getHistogram(inputFile, dqmDirectories["mc_central"], "qT");

  TH2* histogram_uParl_vs_qT_data            = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["data"], histogramName_uParl_vs_qT));
  TH2* histogram_uParl_vs_qT_mc_central      = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["mc_central"], histogramName_uParl_vs_qT));
  TH2* histogram_uParl_vs_qT_mc_shiftUp      = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["mc_shiftUp"], histogramName_uParl_vs_qT));
  TH2* histogram_uParl_vs_qT_mc_shiftDown    = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["mc_shiftDown"], histogramName_uParl_vs_qT));
  TH2* histogram_uParl_vs_qT_mc_wCalibration = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["mc_wCalibration"], histogramName_uParl_vs_qT));

  TGraphAsymmErrors* graph_uParlDivQt_vs_qT_data = makeGraph_uParl_div_qT( 
    "graphs_uParlDivQt_vs_qT_data", 
    "u_{#parallel}/q_{T} vs q_{T}", histogram_uParl_vs_qT_data, histogram_qT_data, mode, false);
  TGraphAsymmErrors* graph_uParlDivQt_vs_qT_mc_central = makeGraph_uParl_div_qT( 
    "graphs_uParlDivQt_vs_qT_mc_central", 
    "u_{#parallel}/q_{T} vs q_{T}", histogram_uParl_vs_qT_mc_central, histogram_qT_mc, mode, false);
  TGraphAsymmErrors* graph_uParlDivQt_vs_qT_mc_shiftUp = makeGraph_uParl_div_qT( 
    "graphs_uParlDivQt_vs_qT_mc_shiftUp", 
    "u_{#parallel}/q_{T} vs q_{T}", histogram_uParl_vs_qT_mc_shiftUp, histogram_qT_mc, mode, false);
  TGraphAsymmErrors* graph_uParlDivQt_vs_qT_mc_shiftDown = makeGraph_uParl_div_qT( 
    "graphs_uParlDivQt_vs_qT_mc_shiftDown", 
    "u_{#parallel}/q_{T} vs q_{T}", histogram_uParl_vs_qT_mc_shiftDown, histogram_qT_mc, mode, false);
  TGraphAsymmErrors* graph_uParlDivQt_vs_qT_mc_wCalibration = makeGraph_uParl_div_qT( 
    "graphs_uParlDivQt_vs_qT_mc_wCalibration", 
    "u_{#parallel}/q_{T} vs q_{T}", histogram_uParl_vs_qT_mc_wCalibration, histogram_qT_mc, mode, false);

  showGraphs(800, 900,
	     graph_uParlDivQt_vs_qT_data, "Data",	     
	     graph_uParlDivQt_vs_qT_mc_central, "Simulation",
	     //graph_uParlDivQt_vs_qT_mc_shiftUp, "Simulation +10%",
	     //graph_uParlDivQt_vs_qT_mc_shiftDown, "Simulation -10%",
	     graph_uParlDivQt_vs_qT_mc_wCalibration, "Simulation, calibrated",
	     //0, "",
	     0, "",
	     0, "",
	     0, "",
	     0, "",
	     0., 300., 10, "q_{T}", 1.2,
	     false, yMin_response, 1.4, 0.75, 1.25, "<-u_{#parallel}>/q_{T}", 1.2, 
	     0.21, 0.69, 
	     TString(Form(
               "plots/makeUnclusteredEnergyPlots_%s_uParlDivQt_vs_qT_jetPtThreshold%1.1f.png", 
	       type_string.data(), type1JetPtThreshold)).ReplaceAll(".", "_").Data());

  delete graph_uParlDivQt_vs_qT_data;
  delete graph_uParlDivQt_vs_qT_mc_central;
  delete graph_uParlDivQt_vs_qT_mc_shiftUp;
  delete graph_uParlDivQt_vs_qT_mc_shiftDown;

  std::string histogramName_unclustered_vs_qT;
  std::string histogramName_unclCorr_vs_qT;
  std::string histogramName_jets_vs_qT;
  std::string histogramName_jetCorr_vs_qT;
  if ( mode == mode_uParl ) {
    histogramName_unclustered_vs_qT = "uParl_unclustered_vs_qT";
    histogramName_unclCorr_vs_qT    = "uParl_unclCorr_vs_qT";
    histogramName_jets_vs_qT        = "uParl_jets_vs_qT";
    histogramName_jetCorr_vs_qT     = "uParl_jetCorr_vs_qT";
  } else if ( mode == mode_uParl_div_qT ) {
    histogramName_unclustered_vs_qT = "uParlDivQt_unclustered_vs_qT";
    histogramName_unclCorr_vs_qT    = "uParlDivQt_unclCorr_vs_qT";
    histogramName_jets_vs_qT        = "uParlDivQt_jets_vs_qT";
    histogramName_jetCorr_vs_qT     = "uParlDivQt_jetCorr_vs_qT";
  } else assert(0);

  TH2* histogram_uParl_unclustered_vs_qT_data = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["data"], histogramName_unclustered_vs_qT));
  TH2* histogram_uParl_unclCorr_vs_qT_data    = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["data"], histogramName_unclCorr_vs_qT));
  TH2* histogram_uParl_jets_vs_qT_data        = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["data"], histogramName_jets_vs_qT));
  TH2* histogram_uParl_jetCorr_vs_qT_data     = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["data"], histogramName_jetCorr_vs_qT));
  TH2* histogram_uParl_unclustered_vs_qT_mc   = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["mc_central"], histogramName_unclustered_vs_qT));
  TH2* histogram_uParl_unclCorr_vs_qT_mc      = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["mc_wCalibration"], histogramName_unclCorr_vs_qT));
  TH2* histogram_uParl_jets_vs_qT_mc          = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["mc_central"], histogramName_jets_vs_qT));
  TH2* histogram_uParl_jetCorr_vs_qT_mc       = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["mc_central"], histogramName_jetCorr_vs_qT));
   
  TGraphAsymmErrors* graph_uParlDivQt_unclustered_vs_qT_data = makeGraph_uParl_div_qT( 
    "graphs_uParlDivQt_unclustered_vs_qT_data", 
    "u_{#parallel}/q_{T} vs q_{T}", histogram_uParl_unclustered_vs_qT_data, histogram_qT_data, mode, false);
  TGraphAsymmErrors* graph_uParlDivQt_unclCorr_vs_qT_data = makeGraph_uParl_div_qT( 
    "graphs_uParlDivQt_unclCorr_vs_qT_data", 
    "u_{#parallel}/q_{T} vs q_{T}", histogram_uParl_unclCorr_vs_qT_data, histogram_qT_data, mode, false);
  TGraphAsymmErrors* graph_uParlDivQt_jets_vs_qT_data = makeGraph_uParl_div_qT( 
    "graphs_uParlDivQt_jets_vs_qT_data", 
    "u_{#parallel}/q_{T} vs q_{T}", histogram_uParl_jets_vs_qT_data, histogram_qT_data, mode, false);
  TGraphAsymmErrors* graph_uParlDivQt_jetCorr_vs_qT_data = makeGraph_uParl_div_qT( 
    "graphs_uParlDivQt_jetCorr_vs_qT_data", 
    "u_{#parallel}/q_{T} vs q_{T}", histogram_uParl_jetCorr_vs_qT_data, histogram_qT_data, mode, false);
  TGraphAsymmErrors* graph_uParlDivQt_unclustered_vs_qT_mc = makeGraph_uParl_div_qT( 
    "graphs_uParlDivQt_unclustered_vs_qT_mc", 
    "u_{#parallel}/q_{T} vs q_{T}", histogram_uParl_unclustered_vs_qT_mc, histogram_qT_mc, mode, false);
  TGraphAsymmErrors* graph_uParlDivQt_unclCorr_vs_qT_mc = makeGraph_uParl_div_qT( 
    "graphs_uParlDivQt_unclCorr_vs_qT_mc", 
    "u_{#parallel}/q_{T} vs q_{T}", histogram_uParl_unclCorr_vs_qT_mc, histogram_qT_mc, mode, false);
  TGraphAsymmErrors* graph_uParlDivQt_jets_vs_qT_mc = makeGraph_uParl_div_qT( 
    "graphs_uParlDivQt_jets_vs_qT_mc", 
    "u_{#parallel}/q_{T} vs q_{T}", histogram_uParl_jets_vs_qT_mc, histogram_qT_mc, mode, false);
  TGraphAsymmErrors* graph_uParlDivQt_jetCorr_vs_qT_mc = makeGraph_uParl_div_qT( 
    "graphs_uParlDivQt_jetCorr_vs_qT_mc", 
    "u_{#parallel}/q_{T} vs q_{T}", histogram_uParl_jetCorr_vs_qT_mc, histogram_qT_mc, mode, false);

  showResponse_in_Components(800, 900,
			     graph_uParlDivQt_unclustered_vs_qT_mc, graph_uParlDivQt_unclustered_vs_qT_data, 
			     graph_uParlDivQt_unclCorr_vs_qT_mc, graph_uParlDivQt_unclCorr_vs_qT_data, 
			     graph_uParlDivQt_jets_vs_qT_mc, graph_uParlDivQt_jets_vs_qT_data, 
			     graph_uParlDivQt_jetCorr_vs_qT_mc, graph_uParlDivQt_jetCorr_vs_qT_data, 
			     false,
			     0., 300., 10, "q_{T}", 1.2,
			     false, -0.3, 1.8, 0.5, 1.5, "<-u_{#parallel}>/q_{T}", 1.2, 
			     0.21, 0.69, 
			     TString(Form(
                               "plots/makeUnclusteredEnergyPlots_%s_uParlDivQt_vs_qT_in_Components_jetPtThreshold%1.1f.png", 
			       type_string.data(), type1JetPtThreshold)).ReplaceAll(".", "_").Data());
  showResponse_in_Components(800, 900,
			     graph_uParlDivQt_unclustered_vs_qT_mc, graph_uParlDivQt_unclustered_vs_qT_data, 
			     graph_uParlDivQt_unclCorr_vs_qT_mc, graph_uParlDivQt_unclCorr_vs_qT_data, 
			     graph_uParlDivQt_jets_vs_qT_mc, graph_uParlDivQt_jets_vs_qT_data, 
			     graph_uParlDivQt_jetCorr_vs_qT_mc, graph_uParlDivQt_jetCorr_vs_qT_data, 
			     true,
			     0., 300., 10, "q_{T}", 1.2,
			     false, -0.3, 1.8, 0.5, 1.5, "<-u_{#parallel}>/q_{T}", 1.2, 
			     0.21, 0.69, 
			     TString(Form(
                               "plots/makeUnclusteredEnergyPlots_%s_uParlDivQt_vs_qT_in_Components_cumulative_jetPtThreshold%1.1f.png", 
			       type_string.data(), type1JetPtThreshold)).ReplaceAll(".", "_").Data());
  
  //return; // CV: make control plots of calibrated response only
  if ( type_string == "pfJetsPlusCandsNotInJet" ) return; // CV: make control plots of calibrated response only

  const int numPileUpNumBins_data = 6;
  double numPileUpBinning_data[numPileUpNumBins_data + 1] = {
    0., 12.5, 15., 17.5, 20., 22.5, 1.e+3
  };
  const int numPileUpNumBins_mc = 6;
  double numPileUpBinning_mc[numPileUpNumBins_mc + 1] = {
    0., 12.5, 15., 17.5, 20., 22.5, 1.e+3
  };

  const int qTnumBins = 34;
  double qTbinning[qTnumBins + 1] = { 
    0., 2.5, 5., 7.5, 10., 12.5, 15., 17.5, 20., 22.5, 25., 27.5, 30., 35., 40., 45., 50., 
    60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 200., 220., 240., 260., 300.
  };

  const int qTnumBins_rebinned = 3;
  double qTbinning_rebinned[qTnumBins_rebinned + 1] = { 
    5., 40., 100., 300. 
  };

  std::map<int, std::map<int, TGraphAsymmErrors*> > graphs_uParl_vs_eta_data; // key = (qTbin, numPileUpBin)
  std::map<int, std::map<int, TGraphAsymmErrors*> > graphs_uParl_vs_eta_mc_central;
  std::map<int, std::map<int, TGraphAsymmErrors*> > graphs_uParl_vs_eta_mc_shiftUp;
  std::map<int, std::map<int, TGraphAsymmErrors*> > graphs_uParl_vs_eta_mc_shiftDown;

  std::map<int, std::map<int, TGraphAsymmErrors*> > graphs_uParl_vs_numPileUp_data; // key = (qTbin, etaBin)
  std::map<int, std::map<int, TGraphAsymmErrors*> > graphs_uParl_vs_numPileUp_mc_central;
  std::map<int, std::map<int, TGraphAsymmErrors*> > graphs_uParl_vs_numPileUp_mc_shiftUp;
  std::map<int, std::map<int, TGraphAsymmErrors*> > graphs_uParl_vs_numPileUp_mc_shiftDown;

  std::map<int, TGraphAsymmErrors*> graphs_offset_vs_qT_data; // key = (etaBin)
  std::map<int, TGraphAsymmErrors*> graphs_slope_vs_qT_data;
  std::map<int, TGraphAsymmErrors*> graphs_offset_vs_qT_mc_central;
  std::map<int, TGraphAsymmErrors*> graphs_slope_vs_qT_mc_central;
  std::map<int, TGraphAsymmErrors*> graphs_offset_vs_qT_mc_shiftUp;
  std::map<int, TGraphAsymmErrors*> graphs_slope_vs_qT_mc_shiftUp;
  std::map<int, TGraphAsymmErrors*> graphs_offset_vs_qT_mc_shiftDown;
  std::map<int, TGraphAsymmErrors*> graphs_slope_vs_qT_mc_shiftDown;

  const int etaNumBins = 32;
  double etaBinning[etaNumBins + 1] = { 
    -5.191, -3.489, -3.139, -2.964, -2.853, -2.5, -2.411, -2.322, -1.93, -1.479, -1.305, -1.131, -0.957, -0.783, -0.522, -0.261, 0.,
    0.261, 0.522, 0.783, 0.957, 1.131, 1.305, 1.479, 1.93, 2.322, 2.411, 2.5, 2.853, 2.964, 3.139, 3.489, 5.191
  };

  const int etaNumBinsJEC = 32;
  double caloResidualJEC[etaNumBinsJEC] = { // CV: values taken from GR_P_V42_AN3
    1.07944, 1.0677, 1.17153, 1.07494, 1.05095, 1.04558, 1.05412, 1.04263, 1.03258, 1.03162, 1.03287, 1.02467, 1.0234, 1.02081, 1.01314, 1.00617,
    1.01001, 1.01757, 1.02256, 1.0253, 1.02713, 1.03329, 1.0311, 1.02803, 1.04663, 1.06123, 1.05364, 1.05095, 1.07494, 1.17153, 1.0677, 1.07944 
  };
  double pfResidualJEC[etaNumBinsJEC] = { // CV: values taken from GR_P_V42_AN3
    1.11195, 1.12798, 1.18941, 1.12886, 1.0827, 1.07055, 1.07341, 1.05673, 1.03475, 1.01773, 1.01945, 1.02173, 1.02334, 1.0208, 1.01602, 1.01115,
    1.01411, 1.01732, 1.02277, 1.02486, 1.02556, 1.01977, 1.01731, 1.03052, 1.06081, 1.07791, 1.07608, 1.0827, 1.12886, 1.18941, 1.12798, 1.11195
  };

  TGraphAsymmErrors* graphResidualJEC_vs_eta = 0;
  if ( type == type_Calo || type == type_PF ) {    
    graphResidualJEC_vs_eta = new TGraphAsymmErrors(etaNumBins);
    for ( int etaBin = 0; etaBin < etaNumBinsJEC; ++etaBin ) {            
      double etaMin = etaBinning[etaBin];
      double etaMax = etaBinning[etaBin + 1];
      
      double x = 0.5*(etaMin + etaMax);
      double xErr = 0.5*(etaMax - etaMin);
      double y = 1.;
      if      ( type == type_Calo ) y = caloResidualJEC[etaBin];
      else if ( type == type_PF   ) y = pfResidualJEC[etaBin];
      else assert(0);
      graphResidualJEC_vs_eta->SetPoint(etaBin, x, y);
      graphResidualJEC_vs_eta->SetPointEXhigh(etaBin, xErr);
      graphResidualJEC_vs_eta->SetPointEXlow(etaBin, xErr);
    }
  }
    
  for ( int qTbin_rebinned = 0; qTbin_rebinned < qTnumBins_rebinned; ++qTbin_rebinned ) {
    double qTmin_rebinned = qTbinning_rebinned[qTbin_rebinned];
    double qTmax_rebinned = qTbinning_rebinned[qTbin_rebinned + 1];

    double qTav_data = getHistogramAvWithinRange(histogram_qT_data, qTmin_rebinned, qTmax_rebinned);
    double qTav_mc = getHistogramAvWithinRange(histogram_qT_mc, qTmin_rebinned, qTmax_rebinned);

    for ( int numPileUpBin = 0; numPileUpBin < numPileUpNumBins_data; ++numPileUpBin ) {
      double numPileUpMin = numPileUpBinning_data[numPileUpBin];
      double numPileUpMax = numPileUpBinning_data[numPileUpBin + 1];

      TH2* histogram_uParl_vs_eta_data = 0;
      for ( int qTbin = 0; qTbin < qTnumBins; ++qTbin ) {
	double qTmin = qTbinning[qTbin];
	double qTmax = qTbinning[qTbin + 1];

	if ( qTmin >= qTmin_rebinned && qTmax <= qTmax_rebinned ) {
	  TString dqmSubDirectory_data = dqmSubDirectory(dqmDirectories["data"], numPileUpMin, numPileUpMax, qTmin, qTmax);
	  TH2* histogram_uParl_vs_eta_data_i = dynamic_cast<TH2*>(getHistogram(inputFile, dqmSubDirectory_data.Data(), "uParl_vs_eta"));
	  if ( !histogram_uParl_vs_eta_data ) {
	    histogram_uParl_vs_eta_data = (TH2*)histogram_uParl_vs_eta_data_i->Clone(std::string(histogram_uParl_vs_eta_data_i->GetName()).append("_cloned").data());
	  } else {
	    histogram_uParl_vs_eta_data->Add(histogram_uParl_vs_eta_data_i);
	  }
	}
      }

      TGraphAsymmErrors* graph_uParl_vs_eta_data = makeGraph_mean_or_rms(
        "uParl_vs_eta_data",         
	"u_{#parallel} vs #eta", histogram_uParl_vs_eta_data, histogram_qT_data, qTmin_rebinned, qTmax_rebinned, mode, kMean, true);
      graphs_uParl_vs_eta_data[qTbin_rebinned][numPileUpBin] = graph_uParl_vs_eta_data;
    }
    for ( int numPileUpBin = 0; numPileUpBin < numPileUpNumBins_mc; ++numPileUpBin ) {
      double numPileUpMin = numPileUpBinning_mc[numPileUpBin];
      double numPileUpMax = numPileUpBinning_mc[numPileUpBin + 1];

      TH2* histogram_uParl_vs_eta_mc_central   = 0;
      TH2* histogram_uParl_vs_eta_mc_shiftUp   = 0;
      TH2* histogram_uParl_vs_eta_mc_shiftDown = 0;
      for ( int qTbin = 0; qTbin < qTnumBins; ++qTbin ) {
	double qTmin = qTbinning[qTbin];
	double qTmax = qTbinning[qTbin + 1];

	if ( qTmin >= qTmin_rebinned && qTmax <= qTmax_rebinned ) {
	  TString dqmSubDirectory_mc_central = dqmSubDirectory(dqmDirectories["mc_central"], numPileUpMin, numPileUpMax, qTmin, qTmax);
	  TH2* histogram_uParl_vs_eta_mc_central_i = dynamic_cast<TH2*>(getHistogram(inputFile, dqmSubDirectory_mc_central.Data(), "uParl_vs_eta"));
	  if ( !histogram_uParl_vs_eta_mc_central ) {
	    histogram_uParl_vs_eta_mc_central = (TH2*)histogram_uParl_vs_eta_mc_central_i->Clone(std::string(histogram_uParl_vs_eta_mc_central_i->GetName()).append("_cloned").data());
	  } else {
	    histogram_uParl_vs_eta_mc_central->Add(histogram_uParl_vs_eta_mc_central_i);
	  }
	  TString dqmSubDirectory_mc_shiftUp = dqmSubDirectory(dqmDirectories["mc_shiftUp"], numPileUpMin, numPileUpMax, qTmin, qTmax);
	  TH2* histogram_uParl_vs_eta_mc_shiftUp_i = dynamic_cast<TH2*>(getHistogram(inputFile, dqmSubDirectory_mc_shiftUp.Data(), "uParl_vs_eta"));
	  if ( !histogram_uParl_vs_eta_mc_shiftUp ) {
	    histogram_uParl_vs_eta_mc_shiftUp = (TH2*)histogram_uParl_vs_eta_mc_shiftUp_i->Clone(std::string(histogram_uParl_vs_eta_mc_shiftUp_i->GetName()).append("_cloned").data());
	  } else {
	    histogram_uParl_vs_eta_mc_shiftUp->Add(histogram_uParl_vs_eta_mc_shiftUp_i);
	  }
	  TString dqmSubDirectory_mc_shiftDown = dqmSubDirectory(dqmDirectories["mc_shiftDown"], numPileUpMin, numPileUpMax, qTmin, qTmax);
	  TH2* histogram_uParl_vs_eta_mc_shiftDown_i = dynamic_cast<TH2*>(getHistogram(inputFile, dqmSubDirectory_mc_shiftDown.Data(), "uParl_vs_eta"));
	  if ( !histogram_uParl_vs_eta_mc_shiftDown ) {
	    histogram_uParl_vs_eta_mc_shiftDown = (TH2*)histogram_uParl_vs_eta_mc_shiftDown_i->Clone(std::string(histogram_uParl_vs_eta_mc_shiftDown_i->GetName()).append("_cloned").data());
	  } else {
	    histogram_uParl_vs_eta_mc_shiftDown->Add(histogram_uParl_vs_eta_mc_shiftDown_i);
	  }
	}
      }

      TGraphAsymmErrors* graph_uParl_vs_eta_mc_central = makeGraph_mean_or_rms(
        "uParl_vs_eta_mc_central",         
	"u_{#parallel} vs #eta", histogram_uParl_vs_eta_mc_central, histogram_qT_mc, qTmin_rebinned, qTmax_rebinned, mode, kMean, true);
      graphs_uParl_vs_eta_mc_central[qTbin_rebinned][numPileUpBin] = graph_uParl_vs_eta_mc_central;

      TGraphAsymmErrors* graph_uParl_vs_eta_mc_shiftUp = makeGraph_mean_or_rms(
        "uParl_vs_eta_mc_shiftUp",         
	"u_{#parallel} vs #eta", histogram_uParl_vs_eta_mc_shiftUp, histogram_qT_mc, qTmin_rebinned, qTmax_rebinned, mode, kMean, true);
      graphs_uParl_vs_eta_mc_shiftUp[qTbin_rebinned][numPileUpBin] = graph_uParl_vs_eta_mc_shiftUp;

      TGraphAsymmErrors* graph_uParl_vs_eta_mc_shiftDown = makeGraph_mean_or_rms(
        "uParl_vs_eta_mc_shiftDown",         
	"u_{#parallel} vs #eta", histogram_uParl_vs_eta_mc_shiftDown, histogram_qT_mc, qTmin_rebinned, qTmax_rebinned, mode, kMean, true);
      graphs_uParl_vs_eta_mc_shiftDown[qTbin_rebinned][numPileUpBin] = graph_uParl_vs_eta_mc_shiftDown;
    }
    
    // WARNING: these control plots make sense only if similar numPileUpBinning is used for MC and Data!
    for ( int numPileUpBin = 0; numPileUpBin < TMath::Min(numPileUpNumBins_data, numPileUpNumBins_mc); ++numPileUpBin ) {
      double numPileUpMin = numPileUpBinning_data[numPileUpBin];
      double numPileUpMax = numPileUpBinning_data[numPileUpBin + 1];
      
      TString outputFileName = Form(
        "plots/makeUnclusteredEnergyPlots_%s_jetPtThreshold%1.1f_uParl_vs_eta_for_qT%1.1fto%1.1f_numPileUp%1.1fto%1.1f", 
	type_string.data(), type1JetPtThreshold, qTmin_rebinned, qTmax_rebinned, numPileUpMin, numPileUpMax);
      outputFileName.ReplaceAll(".", "_");
      outputFileName.Append(".png");
      showGraphs(800, 900,
		 graphs_uParl_vs_eta_data[qTbin_rebinned][numPileUpBin], "Data",
		 graphs_uParl_vs_eta_mc_central[qTbin_rebinned][numPileUpBin], "Simulation",
		 graphs_uParl_vs_eta_mc_shiftUp[qTbin_rebinned][numPileUpBin], "Simulation +10%",
		 graphs_uParl_vs_eta_mc_shiftDown[qTbin_rebinned][numPileUpBin], "Simulation -10%",
		 //0, "",
		 //0, "",
		 graphResidualJEC_vs_eta, "Residual JEC",
		 0, "",
		 0, "",
		 -5.191, +5.191, 10, "#eta", 1.2,
		 false, -0.05, 0.50, 0.5, 2., "<-u_{#parallel}>/q_{T}", 1.2, 
		 0.21, 0.69, 
		 outputFileName.Data());
    }

    TH1* histogram_numPileUp_data = 0;
    TH1* histogram_numPileUp_mc   = 0;

    for ( int numPileUpBin = 0; numPileUpBin < numPileUpNumBins_data; ++numPileUpBin ) {
      double numPileUpMin = numPileUpBinning_data[numPileUpBin];
      double numPileUpMax = numPileUpBinning_data[numPileUpBin + 1];
      
      for ( int qTbin = 0; qTbin < qTnumBins; ++qTbin ) {
	double qTmin = qTbinning[qTbin];
	double qTmax = qTbinning[qTbin + 1];

	if ( qTmin >= qTmin_rebinned && qTmax <= qTmax_rebinned ) {
	  TString dqmSubDirectory_data = dqmSubDirectory(dqmDirectories["data"], numPileUpMin, numPileUpMax, qTmin, qTmax);
	  TH1* histogram_numPileUp = getHistogram(inputFile, dqmSubDirectory_data.Data(), "numPileUp");
	  if ( !histogram_numPileUp_data ) {
	    std::string histogramName_numPileUp_data = Form("%s_dataSum_qT%1.0fto%1.0f", histogram_numPileUp->GetName(), qTmin_rebinned, qTmax_rebinned);
	    histogram_numPileUp_data = (TH1*)histogram_numPileUp->Clone(histogramName_numPileUp_data.data());
	  } else {
	    histogram_numPileUp_data->Add(histogram_numPileUp);
	  }
	}
      }
      double numPileUpAv = getHistogramAvWithinRange(histogram_numPileUp_data, numPileUpMin, numPileUpMax);
      
      for ( int etaBin = 0; etaBin < etaNumBins; ++etaBin ) {
	TGraphAsymmErrors* graph_uParl_vs_numPileUp_data = graphs_uParl_vs_numPileUp_data[qTbin_rebinned][etaBin];
	if ( !graph_uParl_vs_numPileUp_data ) {
	  graph_uParl_vs_numPileUp_data = new TGraphAsymmErrors(numPileUpNumBins_data);
	  graphs_uParl_vs_numPileUp_data[qTbin_rebinned][etaBin] = graph_uParl_vs_numPileUp_data;
	}
	copyGraphPoint(graphs_uParl_vs_eta_data[qTbin_rebinned][numPileUpBin], etaBin, 
		       graph_uParl_vs_numPileUp_data, numPileUpBin, numPileUpAv, numPileUpAv - numPileUpMin, numPileUpMax - numPileUpAv);
      }
    }
    for ( int numPileUpBin = 0; numPileUpBin < numPileUpNumBins_mc; ++numPileUpBin ) {
      double numPileUpMin = numPileUpBinning_mc[numPileUpBin];
      double numPileUpMax = numPileUpBinning_mc[numPileUpBin + 1];
      
      for ( int qTbin = 0; qTbin < qTnumBins; ++qTbin ) {
	double qTmin = qTbinning[qTbin];
	double qTmax = qTbinning[qTbin + 1];

	if ( qTmin >= qTmin_rebinned && qTmax <= qTmax_rebinned ) {
	  TString dqmSubDirectory_mc = dqmSubDirectory(dqmDirectories["mc_central"], numPileUpMin, numPileUpMax, qTmin, qTmax);
	  TH1* histogram_numPileUp = getHistogram(inputFile, dqmSubDirectory_mc.Data(), "numPileUp");
	  if ( !histogram_numPileUp_mc ) {
	    std::string histogramName_numPileUp_mc = Form("%s_mcSum_qT%1.0fto%1.0f", histogram_numPileUp->GetName(), qTmin_rebinned, qTmax_rebinned);
	    histogram_numPileUp_mc = (TH1*)histogram_numPileUp->Clone(histogramName_numPileUp_mc.data());
	  } else {
	    histogram_numPileUp_mc->Add(histogram_numPileUp);
	  }
	}
      }
      double numPileUpAv = getHistogramAvWithinRange(histogram_numPileUp_mc, numPileUpMin, numPileUpMax);
      
      for ( int etaBin = 0; etaBin < etaNumBins; ++etaBin ) {
	TGraphAsymmErrors* graph_uParl_vs_numPileUp_mc_central = graphs_uParl_vs_numPileUp_mc_central[qTbin_rebinned][etaBin];
	if ( !graph_uParl_vs_numPileUp_mc_central ) {
	  graph_uParl_vs_numPileUp_mc_central = new TGraphAsymmErrors(numPileUpNumBins_mc);
	  graphs_uParl_vs_numPileUp_mc_central[qTbin_rebinned][etaBin] = graph_uParl_vs_numPileUp_mc_central;
	}
	copyGraphPoint(graphs_uParl_vs_eta_mc_central[qTbin_rebinned][numPileUpBin], etaBin, 
		       graph_uParl_vs_numPileUp_mc_central, numPileUpBin, numPileUpAv, numPileUpAv - numPileUpMin, numPileUpMax - numPileUpAv);

	TGraphAsymmErrors* graph_uParl_vs_numPileUp_mc_shiftUp = graphs_uParl_vs_numPileUp_mc_shiftUp[qTbin_rebinned][etaBin];
	if ( !graph_uParl_vs_numPileUp_mc_shiftUp ) {
	  graph_uParl_vs_numPileUp_mc_shiftUp = new TGraphAsymmErrors(numPileUpNumBins_mc);
	  graphs_uParl_vs_numPileUp_mc_shiftUp[qTbin_rebinned][etaBin] = graph_uParl_vs_numPileUp_mc_shiftUp;
	}
	copyGraphPoint(graphs_uParl_vs_eta_mc_shiftUp[qTbin_rebinned][numPileUpBin], etaBin, 
		       graph_uParl_vs_numPileUp_mc_shiftUp, numPileUpBin, numPileUpAv, numPileUpAv - numPileUpMin, numPileUpMax - numPileUpAv);

	TGraphAsymmErrors* graph_uParl_vs_numPileUp_mc_shiftDown = graphs_uParl_vs_numPileUp_mc_shiftDown[qTbin_rebinned][etaBin];
	if ( !graph_uParl_vs_numPileUp_mc_shiftDown ) {
	  graph_uParl_vs_numPileUp_mc_shiftDown = new TGraphAsymmErrors(numPileUpNumBins_mc);
	  graphs_uParl_vs_numPileUp_mc_shiftDown[qTbin_rebinned][etaBin] = graph_uParl_vs_numPileUp_mc_shiftDown;
	}
	copyGraphPoint(graphs_uParl_vs_eta_mc_shiftDown[qTbin_rebinned][numPileUpBin], etaBin, 
		       graph_uParl_vs_numPileUp_mc_shiftDown, numPileUpBin, numPileUpAv, numPileUpAv - numPileUpMin, numPileUpMax - numPileUpAv);
      }
    }

    assert(histogram_numPileUp_data);
    double numPileUpAv_data = histogram_numPileUp_data->GetMean();
    assert(histogram_numPileUp_mc);
    double numPileUpAv_mc = histogram_numPileUp_mc->GetMean();

    double numPileUpMinFit =  5.;
    double numPileUpMaxFit = 40.;

    for ( int etaBin = 0; etaBin < etaNumBins; ++etaBin ) {
      TF1* fit_uParl_vs_numPileUp_data = new TF1("fit_uParl_vs_numPileUp_data", Form("[0] + [1]*(x - %1.1f)", numPileUpAv_data), numPileUpMinFit, numPileUpMaxFit);
      graphs_uParl_vs_numPileUp_data[qTbin_rebinned][etaBin]->Fit(fit_uParl_vs_numPileUp_data, "EN");
      double offset_data = fit_uParl_vs_numPileUp_data->GetParameter(0);
      double offsetErr_data = fit_uParl_vs_numPileUp_data->GetParError(0);
      TGraphAsymmErrors* graph_offset_vs_qT_data = graphs_offset_vs_qT_data[etaBin];
      if ( !graph_offset_vs_qT_data ) {
	graph_offset_vs_qT_data = new TGraphAsymmErrors(qTnumBins_rebinned);
	graphs_offset_vs_qT_data[etaBin] = graph_offset_vs_qT_data;
      }
      graph_offset_vs_qT_data->SetPoint(qTbin_rebinned, qTav_data, offset_data);
      graph_offset_vs_qT_data->SetPointError(qTbin_rebinned, qTav_data - qTmin_rebinned, qTmax_rebinned - qTav_data, offsetErr_data, offsetErr_data);
      double slope_data = fit_uParl_vs_numPileUp_data->GetParameter(1);
      double slopeErr_data = fit_uParl_vs_numPileUp_data->GetParError(1);
      TGraphAsymmErrors* graph_slope_vs_qT_data = graphs_slope_vs_qT_data[etaBin];
      if ( !graph_slope_vs_qT_data ) {
	graph_slope_vs_qT_data = new TGraphAsymmErrors(qTnumBins_rebinned);
	graphs_slope_vs_qT_data[etaBin] = graph_slope_vs_qT_data;
      }
      graph_slope_vs_qT_data->SetPoint(qTbin_rebinned, qTav_data, slope_data);
      graph_slope_vs_qT_data->SetPointError(qTbin_rebinned, qTav_data - qTmin_rebinned, qTmax_rebinned - qTav_data, slopeErr_data, slopeErr_data);

      TF1* fit_uParl_vs_numPileUp_mc_central = new TF1("fit_uParl_vs_numPileUp_mc_central", Form("[0] + [1]*(x - %1.1f)", numPileUpAv_mc), numPileUpMinFit, numPileUpMaxFit);
      graphs_uParl_vs_numPileUp_mc_central[qTbin_rebinned][etaBin]->Fit(fit_uParl_vs_numPileUp_mc_central, "EN");
      double offset_mc_central = fit_uParl_vs_numPileUp_mc_central->GetParameter(0);
      double offsetErr_mc_central = fit_uParl_vs_numPileUp_mc_central->GetParError(0);
      TGraphAsymmErrors* graph_offset_vs_qT_mc_central = graphs_offset_vs_qT_mc_central[etaBin];
      if ( !graph_offset_vs_qT_mc_central ) {
	graph_offset_vs_qT_mc_central = new TGraphAsymmErrors(qTnumBins_rebinned);
	graphs_offset_vs_qT_mc_central[etaBin] = graph_offset_vs_qT_mc_central;
      }
      graph_offset_vs_qT_mc_central->SetPoint(qTbin_rebinned, qTav_mc, offset_mc_central);
      graph_offset_vs_qT_mc_central->SetPointError(qTbin_rebinned, qTav_mc - qTmin_rebinned, qTmax_rebinned - qTav_mc, offsetErr_mc_central, offsetErr_mc_central);
      double slope_mc_central = fit_uParl_vs_numPileUp_mc_central->GetParameter(1);
      double slopeErr_mc_central = fit_uParl_vs_numPileUp_mc_central->GetParError(1);
      TGraphAsymmErrors* graph_slope_vs_qT_mc_central = graphs_slope_vs_qT_mc_central[etaBin];
      if ( !graph_slope_vs_qT_mc_central ) {
	graph_slope_vs_qT_mc_central = new TGraphAsymmErrors(qTnumBins_rebinned);
	graphs_slope_vs_qT_mc_central[etaBin] = graph_slope_vs_qT_mc_central;
      }
      graph_slope_vs_qT_mc_central->SetPoint(qTbin_rebinned, qTav_mc, slope_mc_central);
      graph_slope_vs_qT_mc_central->SetPointError(qTbin_rebinned, qTav_mc - qTmin_rebinned, qTmax_rebinned - qTav_mc, slopeErr_mc_central, slopeErr_mc_central);
      
      TF1* fit_uParl_vs_numPileUp_mc_shiftUp = new TF1("fit_uParl_vs_numPileUp_mc_shiftUp", Form("[0] + [1]*(x - %1.1f)", numPileUpAv_mc), numPileUpMinFit, numPileUpMaxFit);
      graphs_uParl_vs_numPileUp_mc_shiftUp[qTbin_rebinned][etaBin]->Fit(fit_uParl_vs_numPileUp_mc_shiftUp, "EN");
      double offset_mc_shiftUp = fit_uParl_vs_numPileUp_mc_shiftUp->GetParameter(0);
      double offsetErr_mc_shiftUp = fit_uParl_vs_numPileUp_mc_shiftUp->GetParError(0);
      TGraphAsymmErrors* graph_offset_vs_qT_mc_shiftUp = graphs_offset_vs_qT_mc_shiftUp[etaBin];
      if ( !graph_offset_vs_qT_mc_shiftUp ) {
	graph_offset_vs_qT_mc_shiftUp = new TGraphAsymmErrors(qTnumBins_rebinned);
	graphs_offset_vs_qT_mc_shiftUp[etaBin] = graph_offset_vs_qT_mc_shiftUp;
      }
      graph_offset_vs_qT_mc_shiftUp->SetPoint(qTbin_rebinned, qTav_mc, offset_mc_shiftUp);
      graph_offset_vs_qT_mc_shiftUp->SetPointError(qTbin_rebinned, qTav_mc - qTmin_rebinned, qTmax_rebinned - qTav_mc, offsetErr_mc_shiftUp, offsetErr_mc_shiftUp);
      double slope_mc_shiftUp = fit_uParl_vs_numPileUp_mc_shiftUp->GetParameter(1);
      double slopeErr_mc_shiftUp = fit_uParl_vs_numPileUp_mc_shiftUp->GetParError(1);
      TGraphAsymmErrors* graph_slope_vs_qT_mc_shiftUp = graphs_slope_vs_qT_mc_shiftUp[etaBin];
      if ( !graph_slope_vs_qT_mc_shiftUp ) {
	graph_slope_vs_qT_mc_shiftUp = new TGraphAsymmErrors(qTnumBins_rebinned);
	graphs_slope_vs_qT_mc_shiftUp[etaBin] = graph_slope_vs_qT_mc_shiftUp;
      }
      graph_slope_vs_qT_mc_shiftUp->SetPoint(qTbin_rebinned, qTav_mc, slope_mc_shiftUp);
      graph_slope_vs_qT_mc_shiftUp->SetPointError(qTbin_rebinned, qTav_mc - qTmin_rebinned, qTmax_rebinned - qTav_mc, slopeErr_mc_shiftUp, slopeErr_mc_shiftUp);

      TF1* fit_uParl_vs_numPileUp_mc_shiftDown = new TF1("fit_uParl_vs_numPileUp_mc_shiftDown", Form("[0] + [1]*(x - %1.1f)", numPileUpAv_mc), numPileUpMinFit, numPileUpMaxFit);
      graphs_uParl_vs_numPileUp_mc_shiftDown[qTbin_rebinned][etaBin]->Fit(fit_uParl_vs_numPileUp_mc_shiftDown, "EN");
      double offset_mc_shiftDown = fit_uParl_vs_numPileUp_mc_shiftDown->GetParameter(0);
      double offsetErr_mc_shiftDown = fit_uParl_vs_numPileUp_mc_shiftDown->GetParError(0);
      TGraphAsymmErrors* graph_offset_vs_qT_mc_shiftDown = graphs_offset_vs_qT_mc_shiftDown[etaBin];
      if ( !graph_offset_vs_qT_mc_shiftDown ) {
	graph_offset_vs_qT_mc_shiftDown = new TGraphAsymmErrors(qTnumBins_rebinned);
	graphs_offset_vs_qT_mc_shiftDown[etaBin] = graph_offset_vs_qT_mc_shiftDown;
      }
      graph_offset_vs_qT_mc_shiftDown->SetPoint(qTbin_rebinned, qTav_mc, offset_mc_shiftDown);
      graph_offset_vs_qT_mc_shiftDown->SetPointError(qTbin_rebinned, qTav_mc - qTmin_rebinned, qTmax_rebinned - qTav_mc, offsetErr_mc_shiftDown, offsetErr_mc_shiftDown);
      double slope_mc_shiftDown = fit_uParl_vs_numPileUp_mc_shiftDown->GetParameter(1);
      double slopeErr_mc_shiftDown = fit_uParl_vs_numPileUp_mc_shiftDown->GetParError(1);
      TGraphAsymmErrors* graph_slope_vs_qT_mc_shiftDown = graphs_slope_vs_qT_mc_shiftDown[etaBin];
      if ( !graph_slope_vs_qT_mc_shiftDown ) {
	graph_slope_vs_qT_mc_shiftDown = new TGraphAsymmErrors(qTnumBins_rebinned);
	graphs_slope_vs_qT_mc_shiftDown[etaBin] = graph_slope_vs_qT_mc_shiftDown;
      }
      graph_slope_vs_qT_mc_shiftDown->SetPoint(qTbin_rebinned, qTav_mc, slope_mc_shiftDown);
      graph_slope_vs_qT_mc_shiftDown->SetPointError(qTbin_rebinned, qTav_mc - qTmin_rebinned, qTmax_rebinned - qTav_mc, slopeErr_mc_shiftDown, slopeErr_mc_shiftDown);      

      double etaMin = etaBinning[etaBin];
      double etaMax = etaBinning[etaBin + 1];

      TString outputFileName = Form( 
        "plots/makeUnclusteredEnergyPlots_%s_jetPtThreshold%1.1f_uParl_vs_eta_for_qT%1.1fto%1.1f_eta%1.1fto%1.1f", 
	type_string.data(), type1JetPtThreshold, qTmin_rebinned, qTmax_rebinned, etaMin, etaMax);
      outputFileName.ReplaceAll(".", "_");
      outputFileName.Append(".png");
      showGraphs_with_Fits(800, 900,
			   new graph_with_Fit(graphs_uParl_vs_numPileUp_data[qTbin_rebinned][etaBin], fit_uParl_vs_numPileUp_data, "Data"),
			   new graph_with_Fit(graphs_uParl_vs_numPileUp_mc_central[qTbin_rebinned][etaBin], fit_uParl_vs_numPileUp_mc_central, "Simulation"),
			   new graph_with_Fit(graphs_uParl_vs_numPileUp_mc_shiftUp[qTbin_rebinned][etaBin], fit_uParl_vs_numPileUp_mc_shiftUp, "Simulation +10%"),
			   new graph_with_Fit(graphs_uParl_vs_numPileUp_mc_shiftDown[qTbin_rebinned][etaBin], fit_uParl_vs_numPileUp_mc_shiftDown, "Simulation -10%"),
			   //0,
			   //0,
			   0,
			   0,
			   0,
			   0., 50., 10, "N_{PU} (mean)", 1.2,
			   false, -0.05, 0.50, 0.5, 2., "<-u_{#parallel}>/q_{T}", 1.2, 
			   0.21, 0.69, 
			   outputFileName.Data());
    }
  }

  TGraphAsymmErrors* graph_offset_vs_eta_data         = new TGraphAsymmErrors(etaNumBins); 
  TGraphAsymmErrors* graph_slope_vs_eta_data          = new TGraphAsymmErrors(etaNumBins); 
  TGraphAsymmErrors* graph_offset_vs_eta_mc_central   = new TGraphAsymmErrors(etaNumBins); 
  TGraphAsymmErrors* graph_slope_vs_eta_mc_central    = new TGraphAsymmErrors(etaNumBins);
  TGraphAsymmErrors* graph_offset_vs_eta_mc_shiftUp   = new TGraphAsymmErrors(etaNumBins); 
  TGraphAsymmErrors* graph_slope_vs_eta_mc_shiftUp    = new TGraphAsymmErrors(etaNumBins);
  TGraphAsymmErrors* graph_offset_vs_eta_mc_shiftDown = new TGraphAsymmErrors(etaNumBins); 
  TGraphAsymmErrors* graph_slope_vs_eta_mc_shiftDown  = new TGraphAsymmErrors(etaNumBins); 

  double qTminFit = 10.;
  double qTmaxFit = 50.;

  for ( int etaBin = 0; etaBin < etaNumBins; ++etaBin ) {
    double etaMin = etaBinning[etaBin];
    double etaMax = etaBinning[etaBin + 1];
    double etaAv = 0.5*(etaMin + etaMax);
    
    TGraphAsymmErrors* graph_offset_vs_qT_data = graphs_offset_vs_qT_data[etaBin];
    TF1* fit_offset_vs_qT_data = new TF1("fit_offset_vs_qT_data", "[0]", qTminFit, qTmaxFit);
    graph_offset_vs_qT_data->Fit(fit_offset_vs_qT_data, "EN");
    double offset_data = fit_offset_vs_qT_data->GetParameter(0);
    double offsetErr_data = fit_offset_vs_qT_data->GetParError(0);
    graph_offset_vs_eta_data->SetPoint(etaBin, etaAv, offset_data);
    graph_offset_vs_eta_data->SetPointError(etaBin, etaAv - etaMin, etaMax - etaAv, offsetErr_data, offsetErr_data);
    TGraphAsymmErrors* graph_slope_vs_qT_data = graphs_slope_vs_qT_data[etaBin];
    TF1* fit_slope_vs_qT_data = new TF1("fit_slope_vs_qT_data", "[0]", qTminFit, qTmaxFit);
    graph_slope_vs_qT_data->Fit(fit_slope_vs_qT_data, "EN");
    double slope_data = fit_slope_vs_qT_data->GetParameter(0);
    double slopeErr_data = fit_slope_vs_qT_data->GetParError(0);
    graph_slope_vs_eta_data->SetPoint(etaBin, etaAv, slope_data);
    graph_slope_vs_eta_data->SetPointError(etaBin, etaAv - etaMin, etaMax - etaAv, slopeErr_data, slopeErr_data);

    TGraphAsymmErrors* graph_offset_vs_qT_mc_central = graphs_offset_vs_qT_mc_central[etaBin];
    TF1* fit_offset_vs_qT_mc_central = new TF1("fit_offset_vs_qT_mc_central", "[0]", qTminFit, qTmaxFit);
    graph_offset_vs_qT_mc_central->Fit(fit_offset_vs_qT_mc_central, "EN");
    double offset_mc_central = fit_offset_vs_qT_mc_central->GetParameter(0);
    double offsetErr_mc_central = fit_offset_vs_qT_mc_central->GetParError(0);
    graph_offset_vs_eta_mc_central->SetPoint(etaBin, etaAv, offset_mc_central);
    graph_offset_vs_eta_mc_central->SetPointError(etaBin, etaAv - etaMin, etaMax - etaAv, offsetErr_mc_central, offsetErr_mc_central);
    TGraphAsymmErrors* graph_slope_vs_qT_mc_central = graphs_slope_vs_qT_mc_central[etaBin];
    TF1* fit_slope_vs_qT_mc_central = new TF1("fit_slope_vs_qT_mc_central", "[0]", qTminFit, qTmaxFit);
    graph_slope_vs_qT_mc_central->Fit(fit_slope_vs_qT_mc_central, "EN");
    double slope_mc_central = fit_slope_vs_qT_mc_central->GetParameter(0);
    double slopeErr_mc_central = fit_slope_vs_qT_mc_central->GetParError(0);
    graph_slope_vs_eta_mc_central->SetPoint(etaBin, etaAv, slope_mc_central);
    graph_slope_vs_eta_mc_central->SetPointError(etaBin, etaAv - etaMin, etaMax - etaAv, slopeErr_mc_central, slopeErr_mc_central);

    TGraphAsymmErrors* graph_offset_vs_qT_mc_shiftUp = graphs_offset_vs_qT_mc_shiftUp[etaBin];
    TF1* fit_offset_vs_qT_mc_shiftUp = new TF1("fit_offset_vs_qT_mc_shiftUp", "[0]", qTminFit, qTmaxFit);
    graph_offset_vs_qT_mc_shiftUp->Fit(fit_offset_vs_qT_mc_shiftUp, "EN");
    double offset_mc_shiftUp = fit_offset_vs_qT_mc_shiftUp->GetParameter(0);
    double offsetErr_mc_shiftUp = fit_offset_vs_qT_mc_shiftUp->GetParError(0);
    graph_offset_vs_eta_mc_shiftUp->SetPoint(etaBin, etaAv, offset_mc_shiftUp);
    graph_offset_vs_eta_mc_shiftUp->SetPointError(etaBin, etaAv - etaMin, etaMax - etaAv, offsetErr_mc_shiftUp, offsetErr_mc_shiftUp);
    TGraphAsymmErrors* graph_slope_vs_qT_mc_shiftUp = graphs_slope_vs_qT_mc_shiftUp[etaBin];
    TF1* fit_slope_vs_qT_mc_shiftUp = new TF1("fit_slope_vs_qT_mc_shiftUp", "[0]", qTminFit, qTmaxFit);
    graph_slope_vs_qT_mc_shiftUp->Fit(fit_slope_vs_qT_mc_shiftUp, "EN");
    double slope_mc_shiftUp = fit_slope_vs_qT_mc_shiftUp->GetParameter(0);
    double slopeErr_mc_shiftUp = fit_slope_vs_qT_mc_shiftUp->GetParError(0);
    graph_slope_vs_eta_mc_shiftUp->SetPoint(etaBin, etaAv, slope_mc_shiftUp);
    graph_slope_vs_eta_mc_shiftUp->SetPointError(etaBin, etaAv - etaMin, etaMax - etaAv, slopeErr_mc_shiftUp, slopeErr_mc_shiftUp);

    TGraphAsymmErrors* graph_offset_vs_qT_mc_shiftDown = graphs_offset_vs_qT_mc_shiftDown[etaBin];
    TF1* fit_offset_vs_qT_mc_shiftDown = new TF1("fit_offset_vs_qT_mc_shiftDown", "[0]", qTminFit, qTmaxFit);
    graph_offset_vs_qT_mc_shiftDown->Fit(fit_offset_vs_qT_mc_shiftDown, "EN");
    double offset_mc_shiftDown = fit_offset_vs_qT_mc_shiftDown->GetParameter(0);
    double offsetErr_mc_shiftDown = fit_offset_vs_qT_mc_shiftDown->GetParError(0);
    graph_offset_vs_eta_mc_shiftDown->SetPoint(etaBin, etaAv, offset_mc_shiftDown);
    graph_offset_vs_eta_mc_shiftDown->SetPointError(etaBin, etaAv - etaMin, etaMax - etaAv, offsetErr_mc_shiftDown, offsetErr_mc_shiftDown);
    TGraphAsymmErrors* graph_slope_vs_qT_mc_shiftDown = graphs_slope_vs_qT_mc_shiftDown[etaBin];
    TF1* fit_slope_vs_qT_mc_shiftDown = new TF1("fit_slope_vs_qT_mc_shiftDown", "[0]", qTminFit, qTmaxFit);
    graph_slope_vs_qT_mc_shiftDown->Fit(fit_slope_vs_qT_mc_shiftDown, "EN");
    double slope_mc_shiftDown = fit_slope_vs_qT_mc_shiftDown->GetParameter(0);
    double slopeErr_mc_shiftDown = fit_slope_vs_qT_mc_shiftDown->GetParError(0);
    graph_slope_vs_eta_mc_shiftDown->SetPoint(etaBin, etaAv, slope_mc_shiftDown);
    graph_slope_vs_eta_mc_shiftDown->SetPointError(etaBin, etaAv - etaMin, etaMax - etaAv, slopeErr_mc_shiftDown, slopeErr_mc_shiftDown);

    TString outputFileName_offset = Form(
      "plots/makeUnclusteredEnergyPlots_%s_jetPtThreshold%1.1f_offset_vs_qT_eta%1.1fto%1.1f", 
      type_string.data(), type1JetPtThreshold, etaMin, etaMax);
    outputFileName_offset.ReplaceAll(".", "_");
    outputFileName_offset.Append(".png");
    showGraphs_with_Fits(800, 900,
			 new graph_with_Fit(graph_offset_vs_qT_data, fit_offset_vs_qT_data, "Data"),
			 new graph_with_Fit(graph_offset_vs_qT_mc_central, fit_offset_vs_qT_mc_central, "Simulation"),
			 new graph_with_Fit(graph_offset_vs_qT_mc_shiftUp, fit_offset_vs_qT_mc_shiftUp, "Simulation +10%"),
			 new graph_with_Fit(graph_offset_vs_qT_mc_shiftDown, fit_offset_vs_qT_mc_shiftDown, "Simulation -10%"),
			 //0, 
			 //0,  
			 0, 
			 0,  
			 0,  
			 0., 300., 30, "q_{T} / GeV", 1.2,
			 false, -0.05, 0.50, 0.5, 2., "offset", 1.2, 
			 0.21, 0.69, 
			 outputFileName_offset.Data());
    TString outputFileName_slope = Form(
      "plots/makeUnclusteredEnergyPlots_%s_jetPtThreshold%1.1f_slope_vs_qT_eta%1.1fto%1.1f", 
      type_string.data(), type1JetPtThreshold, etaMin, etaMax);
    outputFileName_slope.ReplaceAll(".", "_");
    outputFileName_slope.Append(".png");
    showGraphs_with_Fits(800, 900,
			 new graph_with_Fit(graph_slope_vs_qT_data, fit_slope_vs_qT_data, "Data"),
			 new graph_with_Fit(graph_slope_vs_qT_mc_central, fit_slope_vs_qT_mc_central, "Simulation"),
			 new graph_with_Fit(graph_slope_vs_qT_mc_shiftUp, fit_slope_vs_qT_mc_shiftUp, "Simulation +10%"),
			 new graph_with_Fit(graph_slope_vs_qT_mc_shiftDown, fit_slope_vs_qT_mc_shiftDown, "Simulation -10%"),
			 //0, 
			 //0,  
			 0, 
			 0,  
			 0,  
			 0., 300., 30, "q_{T} / GeV", 1.2,
			 false, -0.0015, +0.0035, 0.5, 2., "slope", 1.2, 
			 0.21, 0.69, 
			 outputFileName_slope.Data());
  }
  
  showGraphs(800, 900,
	     graph_offset_vs_eta_data, "Data", 
	     graph_offset_vs_eta_mc_central, "Simulation", 
	     graph_offset_vs_eta_mc_shiftUp, "Simulation +10%",
	     graph_offset_vs_eta_mc_shiftDown, "Simulation -10%",
	     //0, "",
	     //0, "", 
	     graphResidualJEC_vs_eta, "Residual JEC",
	     0, "",
	     0, "", 
	     -5.191, +5.191, 10, "#eta", 1.2,
	     false, -0.05, 0.50, 0.5, 2., "offset", 1.2, 
	     0.21, 0.69, 
	     TString(Form("plots/makeUnclusteredEnergyPlots_%s_jetPtThreshold%1.1f_offset_vs_eta", 
			  type_string.data(), type1JetPtThreshold)).ReplaceAll(".", "_").Data());

  TGraphAsymmErrors* graph_offset_vs_eta_mc_central_div_data   = compRatioGraph("offset_vs_eta_mc_central_div_data",   graph_offset_vs_eta_mc_central,   graph_offset_vs_eta_data, true);
  TGraphAsymmErrors* graph_offset_vs_eta_mc_shiftUp_div_data   = compRatioGraph("offset_vs_eta_mc_shiftUp_div_data",   graph_offset_vs_eta_mc_shiftUp,   graph_offset_vs_eta_data, true);
  TGraphAsymmErrors* graph_offset_vs_eta_mc_shiftDown_div_data = compRatioGraph("offset_vs_eta_mc_shiftDown_div_data", graph_offset_vs_eta_mc_shiftDown, graph_offset_vs_eta_data, true);  
  showGraphs(800, 900,
	     graph_offset_vs_eta_mc_central_div_data, "Simulation/Data", 
	     graph_offset_vs_eta_mc_shiftUp_div_data, "Simulation +10%/Data",
	     graph_offset_vs_eta_mc_shiftDown_div_data, "Simulation -10%/Data",
	     graphResidualJEC_vs_eta, "Residual JEC",
	     //0, "", 
	     //0, "", 
	     0, "",
	     0, "", 
	     0, "", 
	     -5.191, +5.191, 10, "#eta", 1.2,
	     false, 0.50, 2., 0.5, 2., "offset", 1.2, 
	     0.21, 0.69, 
	     TString(Form(
               "plots/makeUnclusteredEnergyPlots_%s_jetPtThreshold%1.1f_offset_vs_eta_mc_div_data", 
	       type_string.data(), type1JetPtThreshold)).ReplaceAll(".", "_").Data());
  delete graph_offset_vs_eta_mc_central_div_data;
  delete graph_offset_vs_eta_mc_shiftUp_div_data;
  delete graph_offset_vs_eta_mc_shiftDown_div_data;

  showGraphs(800, 900,
	     graph_slope_vs_eta_data, "Data", 
	     graph_slope_vs_eta_mc_central, "Simulation", 
	     graph_slope_vs_eta_mc_shiftUp, "Simulation +10%",
	     graph_slope_vs_eta_mc_shiftDown, "Simulation -10%",
	     //0, "", 
	     //0, "", 
	     graphResidualJEC_vs_eta, "Residual JEC",
	     0, "", 
	     0, "", 
	     -5.191, +5.191, 10, "#eta", 1.2,
	     false, -0.0015, +0.0035, 0.5, 2., "slope", 1.2, 
	     0.21, 0.69, 
	     TString(Form("plots/makeUnclusteredEnergyPlots_%s_jetPtThreshold%1.1f_slope_vs_eta", 
			  type_string.data(), type1JetPtThreshold)).ReplaceAll(".", "_").Data());

  TGraphAsymmErrors* graph_slope_vs_eta_mc_central_div_data   = compRatioGraph("slope_vs_eta_mc_central_div_data",   graph_slope_vs_eta_mc_central,   graph_slope_vs_eta_data, true);
  TGraphAsymmErrors* graph_slope_vs_eta_mc_shiftUp_div_data   = compRatioGraph("slope_vs_eta_mc_shiftUp_div_data",   graph_slope_vs_eta_mc_shiftUp,   graph_slope_vs_eta_data, true);
  TGraphAsymmErrors* graph_slope_vs_eta_mc_shiftDown_div_data = compRatioGraph("slope_vs_eta_mc_shiftDown_div_data", graph_slope_vs_eta_mc_shiftDown, graph_slope_vs_eta_data, true);  
  showGraphs(800, 900,
	     graph_slope_vs_eta_mc_central_div_data, "Simulation/Data", 
	     graph_slope_vs_eta_mc_shiftUp_div_data, "Simulation +10%/Data",
	     graph_slope_vs_eta_mc_shiftDown_div_data, "Simulation -10%/Data",
	     graphResidualJEC_vs_eta, "Residual JEC",
	     //0, "", 
	     //0, "", 
	     0, "", 
	     0, "", 
	     0, "", 
	     -5.191, +5.191, 10, "#eta", 1.2,
	     false, 0.50, 2., 0.5, 2., "slope", 1.2, 
	     0.21, 0.69, 
	     TString(Form("plots/makeUnclusteredEnergyPlots_%s_jetPtThreshold%1.1f_slope_vs_eta_mc_div_data", 
			  type_string.data(), type1JetPtThreshold)).ReplaceAll(".", "_").Data());
  delete graph_slope_vs_eta_mc_central_div_data;
  delete graph_slope_vs_eta_mc_shiftUp_div_data;
  delete graph_slope_vs_eta_mc_shiftDown_div_data;

  std::vector<std::string> residualCorrections_offset_data = getResidualCorrections(graph_offset_vs_eta_data, etaNumBins, etaBinning);
  saveResidualCorrections(residualCorrections_offset_data, TString(Form(
    "unclEnResidualCorr_Data_runs190456to208686_%s_jetPtThreshold%1.1f_offset", 
    type_string.data(), type1JetPtThreshold)).ReplaceAll(".", "_").Append(".txt").Data());
  std::vector<std::string> residualCorrections_slope_data = getResidualCorrections(graph_slope_vs_eta_data, etaNumBins, etaBinning);
  saveResidualCorrections(residualCorrections_slope_data, TString(Form(
    "unclEnResidualCorr_Data_runs190456to208686_%s_jetPtThreshold%1.1f_slope", 
    type_string.data(), type1JetPtThreshold)).ReplaceAll(".", "_").Append(".txt").Data());

  std::vector<std::string> residualCorrections_offset_mc = getResidualCorrections(graph_offset_vs_eta_mc_central, etaNumBins, etaBinning);
  saveResidualCorrections(residualCorrections_offset_mc, TString(Form(
    "unclEnResidualCorr_ZplusJets_madgraph_%s_jetPtThreshold%1.1f_offset", 
    type_string.data(), type1JetPtThreshold)).ReplaceAll(".", "_").Append(".txt").Data());
  std::vector<std::string> residualCorrections_slope_mc = getResidualCorrections(graph_slope_vs_eta_mc_central, etaNumBins, etaBinning);
  saveResidualCorrections(residualCorrections_slope_mc, TString(Form(
    "unclEnResidualCorr_ZplusJets_madgraph_%s_jetPtThreshold%1.1f_slope", 
    type_string.data(), type1JetPtThreshold)).ReplaceAll(".", "_").Append(".txt").Data());

  delete inputFile;
}
