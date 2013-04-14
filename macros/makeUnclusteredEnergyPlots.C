
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

enum { kUndefined, kCalo, kPF };

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
    unsigned numBins = histogram->GetNbinsX();
    unsigned numBins_withinRange = 0;
    for ( unsigned iBin = 1; iBin <= numBins; ++iBin ) {
      double binCenter = histogram->GetBinCenter(iBin);
      if ( binCenter >= xMin && binCenter <= xMax ) ++numBins_withinRange;
    }

    //std::cout << "histogram = " << histogram->GetName() << ":" 
    //          << " numBins(" << xMin << ".." << "xMax) = " << numBins_withinRange << ", integral = " << histogram->Integral() << std::endl;
    
    unsigned numBins_rebinned = numBins_withinRange;

    for ( int combineNumBins = 5; combineNumBins >= 2; --combineNumBins ) {
      if ( numBins_withinRange >= (combineNumBins*numBinsMin_rebinned) && (numBins % combineNumBins) == 0 ) {
        numBins_rebinned /= combineNumBins;
        numBins_withinRange /= combineNumBins;
      }
    }

    std::string histogramName_rebinned = std::string(histogram->GetName()).append("_rebinned");
    histogram_rebinned = new TH1D(histogramName_rebinned.data(), histogram->GetTitle(), numBins_rebinned, xMin, xMax);
    if ( !histogram_rebinned->GetSumw2N() ) histogram_rebinned->Sumw2();

    TAxis* xAxis = histogram_rebinned->GetXaxis();
      
    unsigned iBin = 1;
    for ( unsigned iBin_rebinned = 1; iBin_rebinned <= numBins_rebinned; ++iBin_rebinned ) {
      double binContent_rebinned = 0.;
      double binError2_rebinned = 0.;

      double xMin_rebinnedBin = xAxis->GetBinLowEdge(iBin_rebinned);
      double xMax_rebinnedBin = xAxis->GetBinUpEdge(iBin_rebinned);

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

    //std::cout << "histogram(rebinned) = " << histogram_rebinned->GetName() << ":" 
    //          << " numBins = " << histogram_rebinned->GetNbinsX() << ", integral = " << histogram_rebinned->Integral() << std::endl;
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

  TLegend* legend = new TLegend(legendX0, legendY0, legendX0 + 0.32, legendY0 + 0.22, "", "brNDC"); 
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

TGraphAsymmErrors* makeGraph_mean_or_rms(const std::string& name, const std::string& title, 
					 const TH3* histogram_uParl_vs_eta_vs_qT, const TH1* histogram_qT, int qTbin, int mode, bool divideByBinWidth,
					 std::map<int, std::map<int, TH1*> >& histograms_uParl_proj) 
{
  //std::cout << "<makeGraph_mean_or_rms>:" << std::endl;
  //std::cout << " name = " << name << std::endl;

  if ( !(histogram_uParl_vs_eta_vs_qT && histogram_qT) ) return 0;

  const TAxis* qTaxis = histogram_uParl_vs_eta_vs_qT->GetXaxis();
  int qTnumBins = qTaxis->GetNbins();
  assert(qTbin >= 1 && qTbin <= qTnumBins);

  double qTmin = qTaxis->GetBinLowEdge(qTbin);
  double qTmax = qTaxis->GetBinUpEdge(qTbin);

  int binLowIndex = const_cast<TH1*>(histogram_qT)->FindBin(qTmin);
  int binUpIndex  = const_cast<TH1*>(histogram_qT)->FindBin(qTmax);
  histogram_qT->GetXaxis()->SetRange(binLowIndex, binUpIndex);

  double qTav = histogram_qT->GetMean();
  //std::cout << "qT: min = " << qTmin << ", max = " << qTmax << ", av = " << qTav << std::endl;
  
  const TAxis* etaAxis = histogram_uParl_vs_eta_vs_qT->GetYaxis();
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

    TString histogramName_uParl_proj = Form("%s_pz_%i_%i", histogram_uParl_vs_eta_vs_qT->GetName(), qTbin, etaBin);
    TH1D* histogram_uParl_proj = histogram_uParl_vs_eta_vs_qT->ProjectionZ(histogramName_uParl_proj.Data(), qTbin, qTbin, etaBin, etaBin, "e");
    histograms_uParl_proj[qTbin][etaBin] = histogram_uParl_proj;
    // CV: skip (qT, eta) bins with limited event statistics
    //if ( !(histogram_uParl_proj->GetEntries() >= 100) ) continue;

    double y = 0.;
    double yErr = 0.;
    if ( mode == kMean ) {
      y = -histogram_uParl_proj->GetMean()/qTav;
      yErr = histogram_uParl_proj->GetMeanError()/qTav;
      //std::cout << " -uParl/qT = " << y << " +/- " << yErr << " (#entries = " << histogram_uParl_proj->GetEntries() << ")" << std::endl;
    } else if ( mode == kRMS ) {
      y = histogram_uParl_proj->GetRMS()/qTav;
      yErr = histogram_uParl_proj->GetRMSError()/qTav;
      //std::cout << " sigma(-uParl/qT) = " << y << " +/- " << yErr << " (#entries = " << histogram_uParl_proj->GetEntries() << ")" << std::endl;
    } else assert (0);
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

//-------------------------------------------------------------------------------
TGraphAsymmErrors* compRatioGraph(const std::string& ratioGraphName, const TGraph* numerator, const TGraph* denominator)
{
  //std::cout << "<compRatioGraph>:" << std::endl;
  //std::cout << " ratioGraphName = " << ratioGraphName << std::endl;

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
    if ( TMath::Abs(x_denominator - x_numerator) > 1.e-3 ) {
      std::cerr << "Incompatible binning of numerator and denominator histograms !!" << std::endl;
      std::cout << "point #" << iPoint << ": x(numerator) = " << x_numerator << ", x(denominator) = " << x_denominator << std::endl;
      assert(0);
    }
    assert(x_denominator == x_numerator);
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

void showGraphs(double canvasSizeX, double canvasSizeY,
		TGraphAsymmErrors* graph_ref, const std::string& legendEntry_ref,
		TGraphAsymmErrors* graph2, const std::string& legendEntry2,
		TGraphAsymmErrors* graph3, const std::string& legendEntry3,
		TGraphAsymmErrors* graph4, const std::string& legendEntry4,
		TGraphAsymmErrors* graph5, const std::string& legendEntry5,
		TGraphAsymmErrors* graph6, const std::string& legendEntry6,
		TF1* fitFunction, const std::string& legendEntryFit,
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

  int colors[6] = { 1, 2, 3, 4, 6, 7 };
  int markerStyles[6] = { 22, 32, 20, 24, 21, 25 };

  TLegend* legend = new TLegend(legendX0, legendY0, legendX0 + 0.33, legendY0 + 0.25, "", "brNDC"); 
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

  graph_ref->SetLineColor(colors[0]);
  graph_ref->SetLineWidth(1);
  graph_ref->SetMarkerColor(colors[0]);
  graph_ref->SetMarkerStyle(markerStyles[0]);
  graph_ref->SetMarkerSize(1);
  graph_ref->Draw("p");
  legend->AddEntry(graph_ref, legendEntry_ref.data(), "p");

  if ( graph2 ) {
    graph2->SetLineColor(colors[1]);
    graph2->SetLineWidth(1);
    graph2->SetMarkerColor(colors[1]);
    graph2->SetMarkerStyle(markerStyles[1]);
    graph2->SetMarkerSize(1);
    graph2->Draw("p");
    legend->AddEntry(graph2, legendEntry2.data(), "p");
  }
  
  if ( graph3 ) {
    graph3->SetLineColor(colors[2]);
    graph3->SetLineWidth(1);
    graph3->SetMarkerColor(colors[2]);
    graph3->SetMarkerStyle(markerStyles[2]);
    graph3->SetMarkerSize(1);
    graph3->Draw("p");
    legend->AddEntry(graph3, legendEntry3.data(), "p");
  }

  if ( graph4 ) {
    graph4->SetLineColor(colors[3]);
    graph4->SetLineWidth(1);
    graph4->SetMarkerColor(colors[3]);
    graph4->SetMarkerStyle(markerStyles[3]);
    graph4->SetMarkerSize(1);
    graph4->Draw("p");
    legend->AddEntry(graph4, legendEntry4.data(), "p");
  }

  if ( graph5 ) {
    graph5->SetLineColor(colors[4]);
    graph5->SetLineWidth(1);
    graph5->SetMarkerColor(colors[4]);
    graph5->SetMarkerStyle(markerStyles[4]);
    graph5->SetMarkerSize(1);
    graph5->Draw("p");
    legend->AddEntry(graph5, legendEntry5.data(), "p");
  }

  if ( graph6 ) {
    graph6->SetLineColor(colors[5]);
    graph6->SetLineWidth(1);
    graph6->SetMarkerColor(colors[5]);
    graph6->SetMarkerStyle(markerStyles[5]);
    graph6->SetMarkerSize(1);
    //graph6->Draw("p");
    legend->AddEntry(graph6, legendEntry6.data(), "p");
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
  if ( graph2 ) {
    std::string graphName2_div_ref = std::string(graph2->GetName()).append("_div_").append(graph_ref->GetName());
    graph2_div_ref = compRatioGraph(graphName2_div_ref, graph2, graph_ref);
    graph2_div_ref->Draw("p");
  }

  TGraph* graph3_div_ref = 0;
  if ( graph3 ) {
    std::string graphName3_div_ref = std::string(graph3->GetName()).append("_div_").append(graph_ref->GetName());
    graph3_div_ref = compRatioGraph(graphName3_div_ref, graph3, graph_ref);
    graph3_div_ref->Draw("p");
  }

  TGraph* graph4_div_ref = 0;
  if ( graph4 ) {
    std::string graphName4_div_ref = std::string(graph4->GetName()).append("_div_").append(graph_ref->GetName());
    graph4_div_ref = compRatioGraph(graphName4_div_ref, graph4, graph_ref);
    graph4_div_ref->Draw("p");
  }

  TGraph* graph5_div_ref = 0;
  if ( graph5 ) {
    std::string graphName5_div_ref = std::string(graph5->GetName()).append("_div_").append(graph_ref->GetName());
    graph5_div_ref = compRatioGraph(graphName5_div_ref, graph5, graph_ref);
    graph5_div_ref->Draw("p");
  }
 
  //TGraph* graph6_div_ref = 0;
  //if ( graph6 ) {
  // std::string graphName6_div_ref = std::string(graph6->GetName()).append("_div_").append(graph_ref->GetName());
  //  graph6_div_ref = compRatioGraph(graphName6_div_ref, graph6, graph_ref);
  //  graph6_div_ref->Draw("p");
  //}
  if ( graph6 ) {
    graph6->Draw("p");
  }

  if ( fitFunction ) {
    fitFunction->SetLineColor(8);
    fitFunction->SetLineWidth(1);
    fitFunction->Draw("same");
    legend->AddEntry(fitFunction, legendEntryFit.data(), "l");
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
  delete graph2_div_ref;
  delete graph3_div_ref;
  delete graph4_div_ref;
  delete graph5_div_ref;
  //delete graph6_div_ref;
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

void makeUnclusteredEnergyPlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  //std::string type_string = "caloTowers";
  std::string type_string = "pfCands";
  //std::string type_string = "tracks";
  
  TFile* inputFile = 0;
  std::map<std::string, std::string> dqmDirectories;
  int type = 0;
  if ( type_string == "caloTowers" ) {
    std::string inputFilePath = "/afs/cern.ch/user/v/veelken/scratch0/CMSSW_5_3_3_patch2/src/TauAnalysis/RecoTools/test";
    std::string inputFileName = "UnclusteredEnergyAnalyzer_all_caloTowers.root";
    inputFile = new TFile(Form("%s/%s", inputFilePath.data(), inputFileName.data()));
    dqmDirectories["data"]         = "Data_runs203894to208686_caloTowers_central";
    dqmDirectories["mc_central"]   = "ZplusJets_madgraph_caloTowers_central";
    dqmDirectories["mc_shiftUp"]   = "ZplusJets_madgraph_caloTowers_shiftUp";
    dqmDirectories["mc_shiftDown"] = "ZplusJets_madgraph_caloTowers_shiftDown";
    type = kCalo;
  } else if ( type_string == "caloTowersNoHF" ) {
    std::string inputFilePath = "/afs/cern.ch/user/v/veelken/scratch0/CMSSW_5_3_3_patch2/src/TauAnalysis/RecoTools/test";
    std::string inputFileName = "UnclusteredEnergyAnalyzer_all_caloTowersNoHF.root";
    inputFile = new TFile(Form("%s/%s", inputFilePath.data(), inputFileName.data()));
    dqmDirectories["data"]         = "Data_runs203894to208686_caloTowersNoHF_central";
    dqmDirectories["mc_central"]   = "ZplusJets_madgraph_caloTowersNoHF_central";
    dqmDirectories["mc_shiftUp"]   = "ZplusJets_madgraph_caloTowersNoHF_shiftUp";
    dqmDirectories["mc_shiftDown"] = "ZplusJets_madgraph_caloTowersNoHF_shiftDown";
    type = kCalo;
  } else if ( type_string == "pfCands" ) {
    std::string inputFilePath = "/afs/cern.ch/user/v/veelken/scratch0/CMSSW_5_3_3_patch2/src/TauAnalysis/RecoTools/test";
    std::string inputFileName = "UnclusteredEnergyAnalyzer_all_pfCands.root";
    inputFile = new TFile(Form("%s/%s", inputFilePath.data(), inputFileName.data()));
    dqmDirectories["data"]         = "Data_runs203894to208686_pfCands_central";
    dqmDirectories["mc_central"]   = "ZplusJets_madgraph_pfCands_central";
    dqmDirectories["mc_shiftUp"]   = "ZplusJets_madgraph_pfCands_shiftUp";
    dqmDirectories["mc_shiftDown"] = "ZplusJets_madgraph_pfCands_shiftDown";
    type = kPF;
  } else if ( type_string == "tracks" ) {
    std::string inputFilePath = "/afs/cern.ch/user/v/veelken/scratch0/CMSSW_5_3_3_patch2/src/TauAnalysis/RecoTools/test";
    std::string inputFileName = "UnclusteredEnergyAnalyzer_all_tracks.root";
    inputFile = new TFile(Form("%s/%s", inputFilePath.data(), inputFileName.data()));
    dqmDirectories["data"]         = "Data_runs203894to208686_tracks_central";
    dqmDirectories["mc_central"]   = "ZplusJets_madgraph_tracks_central";
    dqmDirectories["mc_shiftUp"]   = "";
    dqmDirectories["mc_shiftDown"] = "";
    type = kPF;
  } else if ( type_string == "genParticles" ) {
    std::string inputFilePath = "/afs/cern.ch/user/v/veelken/scratch0/CMSSW_5_3_3_patch2/src/TauAnalysis/RecoTools/test";
    std::string inputFileName = "UnclusteredEnergyAnalyzer_all_pfCands.root";
    inputFile = new TFile(Form("%s/%s", inputFilePath.data(), inputFileName.data()));
    dqmDirectories["data"]         = "";
    dqmDirectories["mc_central"]   = "ZplusJets_madgraph_genParticles_central";
    dqmDirectories["mc_shiftUp"]   = "";
    dqmDirectories["mc_shiftDown"] = "";
    type = kPF;
  } else {
    std::cerr << "Invalid Configuration Parameter 'type' = " << type  << " !!" << std::endl;
    assert(0);
  }

/*
  TH2* histogram_uParl_vs_qT_absZvtxLt5_data   = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["data"], "uParl_vs_qT_absZvtxLt5"));
  TH2* histogram_uParl_vs_qT_absZvtx5to15_data = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["data"], "uParl_vs_qT_absZvtx5to15"));
  TH2* histogram_uParl_vs_qT_absZvtxGt15_data  = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["data"], "uParl_vs_qT_absZvtxGt15"));
  TH2* histogram_uParl_vs_qT_data              = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["data"], "uParl_vs_qT"));

  TH2* histogram_uParl_vs_qT_absZvtxLt5_mc     = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["mc_central"],      "uParl_vs_qT_absZvtxLt5"));
  TH2* histogram_uParl_vs_qT_absZvtx5to15_mc   = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["mc_central"],      "uParl_vs_qT_absZvtx5to15"));
  TH2* histogram_uParl_vs_qT_absZvtxGt15_mc    = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["mc_central"],      "uParl_vs_qT_absZvtxGt15"));
  TH2* histogram_uParl_vs_qT_mc                = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["mc_central"],      "uParl_vs_qT"));
  
  for ( int iBinX = 1; iBinX <= 8; ++iBinX ) {
    TH1* histogram_uParl_absZvtxLt5_data   = getHistogramProjY(histogram_uParl_vs_qT_absZvtxLt5_data, iBinX);
    TH1* histogram_uParl_absZvtx5to15_data = getHistogramProjY(histogram_uParl_vs_qT_absZvtx5to15_data, iBinX);
    TH1* histogram_uParl_absZvtxGt15_data  = getHistogramProjY(histogram_uParl_vs_qT_absZvtxGt15_data, iBinX);
    TH1* histogram_uParl_data              = getHistogramProjY(histogram_uParl_vs_qT_data, iBinX);
    
    TH1* histogram_uParl_absZvtxLt5_mc     = getHistogramProjY(histogram_uParl_vs_qT_absZvtxLt5_mc, iBinX);
    TH1* histogram_uParl_absZvtx5to15_mc   = getHistogramProjY(histogram_uParl_vs_qT_absZvtx5to15_mc, iBinX);
    TH1* histogram_uParl_absZvtxGt15_mc    = getHistogramProjY(histogram_uParl_vs_qT_absZvtxGt15_mc, iBinX);
    TH1* histogram_uParl_mc                = getHistogramProjY(histogram_uParl_vs_qT_mc, iBinX);

    TAxis* xAxis = histogram_uParl_vs_qT_absZvtxLt5_data->GetXaxis();
    double qTmin = xAxis->GetBinLowEdge(iBinX);
    double qTmax = xAxis->GetBinUpEdge(iBinX);

    TString outputFileName = Form("plots/makeUnclusteredEnergyPlots_uParl_for_qT%1.1fto%1.1f", qTmin, qTmax);
    outputFileName.ReplaceAll(".", "_");
    outputFileName.Append(".png");
    showDistribution(800, 900,
		     //histogram_uParl_absZvtxLt5_data, "Data: |z_{Vtx}| < 5cm",
		     //histogram_uParl_absZvtx5to15_data, "Data: 5 < |z_{Vtx}| < 15cm",
		     //histogram_uParl_absZvtxGt15_data, "Data: |z_{Vtx}| > 15cm",
		     histogram_uParl_data, "Data",
		     //histogram_uParl_absZvtxLt5_mc, "MC: |z_{Vtx}| < 5cm",
		     //histogram_uParl_absZvtx5to15_mc, "MC: 5 < |z_{Vtx}| < 15cm",
		     //histogram_uParl_absZvtxGt15_mc, "MC: |z_{Vtx}| > 15cm",
		     histogram_uParl_mc, "MC",
		     0, "",
		     0, "",
		     0, "",
		     0, "",
		     -100.0, +75.0, 70, "u_{#parallel} / GeV", 1.2, 
		     true, 1.e-6, 1.e+1, "a.u.", 1.2,
		     0.64, 0.72, 
		     outputFileName.Data());

    delete histogram_uParl_absZvtxLt5_data;
    delete histogram_uParl_absZvtx5to15_data;
    delete histogram_uParl_absZvtxGt15_data;
    delete histogram_uParl_data;
    delete histogram_uParl_absZvtxLt5_mc;
    delete histogram_uParl_absZvtx5to15_mc;
    delete histogram_uParl_absZvtxGt15_mc;
    delete histogram_uParl_mc;
  }
  
  TH2* histogram_uParlDivQt_vs_qT_absZvtxLt5_data   = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["data"], "uParlDivQt_vs_qT_absZvtxLt5"));
  TH2* histogram_uParlDivQt_vs_qT_absZvtx5to15_data = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["data"], "uParlDivQt_vs_qT_absZvtx5to15"));
  TH2* histogram_uParlDivQt_vs_qT_absZvtxGt15_data  = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["data"], "uParlDivQt_vs_qT_absZvtxGt15"));
  TH2* histogram_uParlDivQt_vs_qT_data              = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["data"], "uParlDivQt_vs_qT"));

  TH2* histogram_uParlDivQt_vs_qT_absZvtxLt5_mc     = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["mc_central"],      "uParlDivQt_vs_qT_absZvtxLt5"));
  TH2* histogram_uParlDivQt_vs_qT_absZvtx5to15_mc   = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["mc_central"],      "uParlDivQt_vs_qT_absZvtx5to15"));
  TH2* histogram_uParlDivQt_vs_qT_absZvtxGt15_mc    = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["mc_central"],      "uParlDivQt_vs_qT_absZvtxGt15"));
  TH2* histogram_uParlDivQt_vs_qT_mc                = dynamic_cast<TH2*>(getHistogram(inputFile, dqmDirectories["mc_central"],      "uParlDivQt_vs_qT"));
  
  for ( int iBinX = 1; iBinX <= 8; ++iBinX ) {
    TH1* histogram_uParlDivQt_absZvtxLt5_data   = getHistogramProjY(histogram_uParl_vs_qT_absZvtxLt5_data, iBinX);
    TH1* histogram_uParlDivQt_absZvtx5to15_data = getHistogramProjY(histogram_uParl_vs_qT_absZvtx5to15_data, iBinX);
    TH1* histogram_uParlDivQt_absZvtxGt15_data  = getHistogramProjY(histogram_uParl_vs_qT_absZvtxGt15_data, iBinX);
    TH1* histogram_uParlDivQt_data              = getHistogramProjY(histogram_uParl_vs_qT_data, iBinX);
    
    TH1* histogram_uParlDivQt_absZvtxLt5_mc     = getHistogramProjY(histogram_uParl_vs_qT_absZvtxLt5_mc, iBinX);
    TH1* histogram_uParlDivQt_absZvtx5to15_mc   = getHistogramProjY(histogram_uParl_vs_qT_absZvtx5to15_mc, iBinX);
    TH1* histogram_uParlDivQt_absZvtxGt15_mc    = getHistogramProjY(histogram_uParl_vs_qT_absZvtxGt15_mc, iBinX);
    TH1* histogram_uParlDivQt_mc                = getHistogramProjY(histogram_uParl_vs_qT_mc, iBinX);

    TAxis* xAxis = histogram_uParlDivQt_vs_qT_absZvtxLt5_data->GetXaxis();
    double qTmin = xAxis->GetBinLowEdge(iBinX);
    double qTmax = xAxis->GetBinUpEdge(iBinX);

    TString outputFileName = Form("plots/makeUnclusteredEnergyPlots_uParlDivQt_for_qT%1.1fto%1.1f", qTmin, qTmax);
    outputFileName.ReplaceAll(".", "_");
    outputFileName.Append(".png");
    showDistribution(800, 900,
		     //histogram_uParlDivQt_absZvtxLt5_data, "Data: |z_{Vtx}| < 5cm",
		     //histogram_uParlDivQt_absZvtx5to15_data, "Data: 5 < |z_{Vtx}| < 15cm",
		     //histogram_uParlDivQt_absZvtxGt15_data, "Data: |z_{Vtx}| > 15cm",
		     histogram_uParlDivQt_data, "Data",
		     //histogram_uParlDivQt_absZvtxLt5_mc, "MC: |z_{Vtx}| < 5cm",
		     //histogram_uParlDivQt_absZvtx5to15_mc, "MC: 5 < |z_{Vtx}| < 15cm",
		     //histogram_uParlDivQt_absZvtxGt15_mc, "MC: |z_{Vtx}| > 15cm",
		     histogram_uParlDivQt_mc, "MC",
		     0, "",
		     0, "",
		     0, "",
		     0, "",
		     -100.0, +75.0, 70, "u_{#parallel} / GeV", 1.2, 
		     true, 1.e-6, 1.e+1, "a.u.", 1.2,
		     0.61, 0.72, 
		     outputFileName.Data());

    delete histogram_uParlDivQt_absZvtxLt5_data;
    delete histogram_uParlDivQt_absZvtx5to15_data;
    delete histogram_uParlDivQt_absZvtxGt15_data;
    delete histogram_uParlDivQt_data;
    delete histogram_uParlDivQt_absZvtxLt5_mc;
    delete histogram_uParlDivQt_absZvtx5to15_mc;
    delete histogram_uParlDivQt_absZvtxGt15_mc;
    delete histogram_uParlDivQt_mc;
  }
 */  
  TH3* histogram_uParl_vs_eta_vs_qT_data         = dynamic_cast<TH3*>(getHistogram(inputFile, dqmDirectories["data"], "uParl_vs_eta_vs_qT"));
  TH1* histogram_qT_data                         = getHistogram(inputFile, dqmDirectories["data"], "qT");
  TH3* histogram_uParl_vs_eta_vs_qT_mc_central   = dynamic_cast<TH3*>(getHistogram(inputFile, dqmDirectories["mc_central"], "uParl_vs_eta_vs_qT"));
  TH1* histogram_qT_mc_central                   = getHistogram(inputFile, dqmDirectories["mc_central"], "qT");
  TH3* histogram_uParl_vs_eta_vs_qT_mc_shiftUp   = dynamic_cast<TH3*>(getHistogram(inputFile, dqmDirectories["mc_shiftUp"], "uParl_vs_eta_vs_qT"));
  TH1* histogram_qT_mc_shiftUp                   = getHistogram(inputFile, dqmDirectories["mc_shiftUp"], "qT");
  TH3* histogram_uParl_vs_eta_vs_qT_mc_shiftDown = dynamic_cast<TH3*>(getHistogram(inputFile, dqmDirectories["mc_shiftDown"], "uParl_vs_eta_vs_qT"));
  TH1* histogram_qT_mc_shiftDown                 = getHistogram(inputFile, dqmDirectories["mc_shiftDown"], "qT");
    
  std::map<int, TGraphAsymmErrors*> graphs_uParl_vs_eta_data; // key = qTbin
  std::map<int, TGraphAsymmErrors*> graphs_uParl_vs_eta_mc_central;
  std::map<int, TGraphAsymmErrors*> graphs_uParl_vs_eta_mc_shiftUp;
  std::map<int, TGraphAsymmErrors*> graphs_uParl_vs_eta_mc_shiftDown;

  std::map<int, std::map<int, TH1*> > histograms_uParl_proj_data; // key = (qTbin, etaBin)
  std::map<int, std::map<int, TH1*> > histograms_uParl_proj_mc_central;
  std::map<int, std::map<int, TH1*> > histograms_uParl_proj_mc_shiftUp;
  std::map<int, std::map<int, TH1*> > histograms_uParl_proj_mc_shiftDown;

  const int numEtaBins = 32;
  double caloResidualJEC[numEtaBins] = { // CV: values taken from GR_P_V42_AN3
    1.07944, 1.0677, 1.17153, 1.07494, 1.05095, 1.04558, 1.05412, 1.04263, 1.03258, 1.03162, 1.03287, 1.02467, 1.0234, 1.02081, 1.01314, 1.00617,
    1.01001, 1.01757, 1.02256, 1.0253, 1.02713, 1.03329, 1.0311, 1.02803, 1.04663, 1.06123, 1.05364, 1.05095, 1.07494, 1.17153, 1.0677, 1.07944 
  };
  double pfResidualJEC[numEtaBins] = { // CV: values taken from GR_P_V42_AN3
    1.11195, 1.12798, 1.18941, 1.12886, 1.0827, 1.07055, 1.07341, 1.05673, 1.03475, 1.01773, 1.01945, 1.02173, 1.02334, 1.0208, 1.01602, 1.01115,
    1.01411, 1.01732, 1.02277, 1.02486, 1.02556, 1.01977, 1.01731, 1.03052, 1.06081, 1.07791, 1.07608, 1.0827, 1.12886, 1.18941, 1.12798, 1.11195
  };
  
  const int qTnumBins = 34;
  const int etaNumBins = histogram_uParl_vs_eta_vs_qT_data->GetNbinsY();

  TGraphAsymmErrors* graphResidualJEC_vs_eta = new TGraphAsymmErrors(etaNumBins);
  for ( int etaBin = 1; etaBin <= etaNumBins; ++etaBin ) {
    
    TAxis* etaAxis = histogram_uParl_vs_eta_vs_qT_data->GetYaxis();
    double etaMin = etaAxis->GetBinLowEdge(etaBin);
    double etaMax = etaAxis->GetBinUpEdge(etaBin);

    double x = 0.5*(etaMin + etaMax);
    double xErr = 0.5*(etaMax - etaMin);
    double y = 1.;
    if      ( type == kCalo ) y = caloResidualJEC[etaBin - 1];
    else if ( type == kPF   ) y = pfResidualJEC[etaBin - 1];
    graphResidualJEC_vs_eta->SetPoint(etaBin - 1, x, y);
    graphResidualJEC_vs_eta->SetPointEXhigh(etaBin - 1, xErr);
    graphResidualJEC_vs_eta->SetPointEXlow(etaBin - 1, xErr);
  }

  std::vector<double> qTbinning;
  for ( int qTbin = 1; qTbin <= qTnumBins; ++qTbin ) {
    TGraphAsymmErrors* graph_uParl_vs_eta_data = makeGraph_mean_or_rms(
      "uParl_vs_eta_data",         
      "u_{#parallel} vs #eta", histogram_uParl_vs_eta_vs_qT_data, histogram_qT_data, qTbin, kMean, true,
      histograms_uParl_proj_data);
    TGraphAsymmErrors* graph_uParl_vs_eta_mc_central = makeGraph_mean_or_rms(
      "uParl_vs_eta_mc_central",   
      "u_{#parallel} vs #eta", histogram_uParl_vs_eta_vs_qT_mc_central, histogram_qT_mc_central, qTbin, kMean, true,
      histograms_uParl_proj_mc_central);
    TGraphAsymmErrors* graph_uParl_vs_eta_mc_shiftUp = makeGraph_mean_or_rms(
      "uParl_vs_eta_mc_shiftUp",   
      "u_{#parallel} vs #eta", histogram_uParl_vs_eta_vs_qT_mc_shiftUp, histogram_qT_mc_shiftUp, qTbin, kMean, true,
      histograms_uParl_proj_mc_shiftUp);
    TGraphAsymmErrors* graph_uParl_vs_eta_mc_shiftDown = makeGraph_mean_or_rms(
      "uParl_vs_eta_mc_shiftDown", 
      "u_{#parallel} vs #eta", histogram_uParl_vs_eta_vs_qT_mc_shiftDown, histogram_qT_mc_shiftDown, qTbin, kMean, true,
      histograms_uParl_proj_mc_shiftDown);
    
    TAxis* qTaxis = histogram_uParl_vs_eta_vs_qT_data->GetXaxis();
    double qTmin = qTaxis->GetBinLowEdge(qTbin);
    double qTmax = qTaxis->GetBinUpEdge(qTbin);
    qTbinning.push_back(qTmin);
    if ( qTbin == qTnumBins ) qTbinning.push_back(qTmax);

    TString outputFileName = Form("plots/makeUnclusteredEnergyPlots_%s_uParl_vs_eta_for_qT%1.1fto%1.1f", type_string.data(), qTmin, qTmax);
    outputFileName.ReplaceAll(".", "_");
    outputFileName.Append(".png");
    showGraphs(800, 900,
	       graph_uParl_vs_eta_data, "Data",
	       graph_uParl_vs_eta_mc_central, "Simulation",
	       //graph_uParl_vs_eta_mc_shiftUp, "Simulation +10%",
	       //graph_uParl_vs_eta_mc_shiftDown, "Simulation -10%",
	       0, "",
	       0, "",
	       0, "",
	       graphResidualJEC_vs_eta, "Residual JEC",
	       0, "",
	       -5.191, +5.191, 10, "#eta", 1.2,
	       false, -0.05, 0.50, 0.5, 2., "<-u_{#parallel}>/q_{T}", 1.2, 
	       0.61, 0.69, 
	       outputFileName.Data());
    
    //delete graph_uParl_vs_eta_data;
    //delete graph_uParl_vs_eta_mc_central;
    //delete graph_uParl_vs_eta_mc_shiftUp;
    //delete graph_uParl_vs_eta_mc_shiftDown;
    graphs_uParl_vs_eta_data[qTbin - 1]         = graph_uParl_vs_eta_data;
    graphs_uParl_vs_eta_mc_central[qTbin - 1]   = graph_uParl_vs_eta_mc_central;
    graphs_uParl_vs_eta_mc_shiftUp[qTbin - 1]   = graph_uParl_vs_eta_mc_shiftUp;
    graphs_uParl_vs_eta_mc_shiftDown[qTbin - 1] = graph_uParl_vs_eta_mc_shiftDown;
  }
  assert(qTbinning.size() == (qTnumBins + 1));

  delete graphResidualJEC_vs_eta;

  std::vector<std::string> residualCorrections;
  residualCorrections.push_back("{1         JetEta              1          JetPt               [0]     Correction     L2Relative}");

  for ( int etaBin = 1; etaBin <= etaNumBins; ++etaBin ) {
    TGraphAsymmErrors* graph_uParl_vs_qT_data         = makeGraph_vs_qT(graphs_uParl_vs_eta_data, qTbinning, etaBin - 1);
    TGraphAsymmErrors* graph_uParl_vs_qT_mc_central   = makeGraph_vs_qT(graphs_uParl_vs_eta_mc_central, qTbinning, etaBin - 1);
    TGraphAsymmErrors* graph_uParl_vs_qT_mc_shiftUp   = makeGraph_vs_qT(graphs_uParl_vs_eta_mc_shiftUp, qTbinning, etaBin - 1);
    TGraphAsymmErrors* graph_uParl_vs_qT_mc_shiftDown = makeGraph_vs_qT(graphs_uParl_vs_eta_mc_shiftDown, qTbinning, etaBin - 1);
    
    TAxis* etaAxis = histogram_uParl_vs_eta_vs_qT_data->GetYaxis();
    double etaMin = etaAxis->GetBinLowEdge(etaBin);
    double etaMax = etaAxis->GetBinUpEdge(etaBin);

    TGraphAsymmErrors* graphResidualJEC_vs_qT = new TGraphAsymmErrors(etaNumBins);
    for ( int qTbin = 1; qTbin <= qTnumBins; ++qTbin ) {
    
      TAxis* qTaxis = histogram_uParl_vs_eta_vs_qT_data->GetXaxis();
      double qTmin = qTaxis->GetBinLowEdge(qTbin);
      double qTmax = qTaxis->GetBinUpEdge(qTbin);

      double x = 0.5*(qTmin + qTmax);
      double xErr = 0.5*(qTmax - qTmin);
      double y = 1.;
      if      ( type == kCalo ) y = caloResidualJEC[etaBin - 1];
      else if ( type == kPF   ) y = pfResidualJEC[etaBin - 1];
      graphResidualJEC_vs_qT->SetPoint(qTbin - 1, x, y);
      graphResidualJEC_vs_qT->SetPointEXhigh(qTbin - 1, xErr);
      graphResidualJEC_vs_qT->SetPointEXlow(qTbin - 1, xErr);
    }
    
    TF1* fitResidualCorr = new TF1("residualCorr","[0]", 5., 150.);
    TGraphAsymmErrors* graph_uParl_vs_qT_mc_div_data = compRatioGraph("uParl_vs_qT_mc_div_data", graph_uParl_vs_qT_mc_central, graph_uParl_vs_qT_data);    
    graph_uParl_vs_qT_mc_div_data->Fit(fitResidualCorr, "E");
    std::cout << "eta = " << etaMin << ".." << etaMax << ": residual Corr. = " << fitResidualCorr->GetParameter(0) << " +/- " << fitResidualCorr->GetParError(0)
	      <<" (JEC value = " << graphResidualJEC_vs_qT->GetY()[etaBin - 1] << ")" << std::endl;
    residualCorrections.push_back(Form("%1.3f         %1.3f              3              3           3500        %1.5f        ", etaMin, etaMax, fitResidualCorr->GetParameter(0)));

    TString outputFileName = Form("plots/makeUnclusteredEnergyPlots_%s_uParl_vs_qT_for_eta%1.3fto%1.3f", type_string.data(), etaMin, etaMax);
    outputFileName.ReplaceAll(".", "_");
    outputFileName.Append(".png");
    showGraphs(800, 900,
	       graph_uParl_vs_qT_data, "Data",
	       graph_uParl_vs_qT_mc_central, "Simulation",
	       //graph_uParl_vs_qT_mc_shiftUp, "Simulation +10%",
	       //graph_uParl_vs_qT_mc_shiftDown, "Simulation -10%",
	       0, "",
	       0, "",
	       0, "",
	       graphResidualJEC_vs_qT, "Residual JEC",
	       fitResidualCorr, "Fit",
	       0., 300., 10, "q_{T} / GeV", 1.2,
	       false, -0.05, 0.50, 0.5, 2., "<-u_{#parallel}>/q_{T}", 1.2, 
	       0.61, 0.69, 
	       outputFileName.Data());
    
    delete graphResidualJEC_vs_qT;
    delete fitResidualCorr;
  }
/*
  for ( int qTbin = 1; qTbin <= qTnumBins; ++qTbin ) {
    for ( int etaBin = 1; etaBin <= etaNumBins; ++etaBin ) {
      TAxis* qTaxis = histogram_uParl_vs_eta_vs_qT_data->GetXaxis();
      double qTmin = qTaxis->GetBinLowEdge(qTbin);
      double qTmax = qTaxis->GetBinUpEdge(qTbin);

      TAxis* etaAxis = histogram_uParl_vs_eta_vs_qT_data->GetYaxis();
      double etaMin = etaAxis->GetBinLowEdge(etaBin);
      double etaMax = etaAxis->GetBinUpEdge(etaBin);

      TString outputFileName = Form("plots/makeUnclusteredEnergyPlots_uParl_for_qT%1.1fto%1.1f_and_eta%1.3fto%1.3f", qTmin, qTmax, etaMin, etaMax);
      outputFileName.ReplaceAll(".", "_");
      outputFileName.Append(".png");
      showDistribution(800, 900,
		       histograms_uParl_proj_data[qTbin][etaBin], "Data",
		       histograms_uParl_proj_mc_central[qTbin][etaBin], "Simulation",
		       //histograms_uParl_proj_mc_shiftUp[qTbin][etaBin], "Simulation +10%",
		       //histograms_uParl_proj_mc_shiftDown[qTbin][etaBin], "Simulation -10%",
		       0, "",
		       0, "",
		       0, "",
		       0, "",
		       -(0.5*(qTmax + qTmin) + 75.0), +75.0, TMath::Nint((0.5*(qTmax + qTmin) + 2.*75.0)/2.5), "u_{#parallel} / GeV", 1.2,
		       true, 1.e-6, 1.e+2, -0.50, +0.50, "a.u.", 1.2,
		       0.61, 0.72, 
		       outputFileName.Data());
    }
  }
 */
  showDistribution(800, 900,
		   histogram_qT_data, "Data",
		   histogram_qT_mc_central, "Simulation",
		   0, "",
		   0, "",
		   0, "",
		   0, "",
		   0.0, 50.0, 100, "q_{T} / GeV", 1.2, 
		   true, 1.e-6, 1.e+1, -0.10, +0.10, "a.u.", 1.2,
		   0.61, 0.72, 
		   Form("plots/makeUnclusteredEnergyPlots_%s_qT.png", type_string.data()));

  TGraphAsymmErrors* graph_qTmean_vs_qTbinCenter_data = new TGraphAsymmErrors(qTnumBins);
  TGraphAsymmErrors* graph_qTmean_vs_qTbinCenter_mc   = new TGraphAsymmErrors(qTnumBins);
  for ( int qTbin = 1; qTbin <= qTnumBins; ++qTbin ) {

    TAxis* qTaxis = histogram_uParl_vs_eta_vs_qT_data->GetXaxis();
    double qTmin = qTaxis->GetBinLowEdge(qTbin);
    double qTmax = qTaxis->GetBinUpEdge(qTbin);
    double qTbinCenter = 0.5*(qTmax + qTmin);
    double qTbinWidth = qTmax - qTmin;

    int binLowIndex = const_cast<TH1*>(histogram_qT_data)->FindBin(qTmin);
    int binUpIndex  = const_cast<TH1*>(histogram_qT_data)->FindBin(qTmax);
    histogram_qT_data->GetXaxis()->SetRange(binLowIndex, binUpIndex);
    histogram_qT_mc_central->GetXaxis()->SetRange(binLowIndex, binUpIndex);

    double qTav_data = histogram_qT_data->GetMean();
    double qTavErr_data = histogram_qT_data->GetMeanError();
    double qTav_mc = histogram_qT_mc_central->GetMean();
    double qTavErr_mc = histogram_qT_mc_central->GetMeanError();
    //std::cout << "qT: min = " << qTmin << ", max = " << qTmax << "," 
    //          << " av(data) = " << qTav_data << " +/- " << qTavErr_data << "," 
    //	        << " av(mc) = " << qTav_mc << " +/- " << qTavErr_mc << std::endl;

    graph_qTmean_vs_qTbinCenter_data->SetPoint(qTbin - 1, qTbinCenter, qTav_data/qTbinCenter);
    graph_qTmean_vs_qTbinCenter_data->SetPointError(qTbin - 1, 0.5*qTbinWidth, 0.5*qTbinWidth, qTavErr_data/qTbinCenter, qTavErr_data/qTbinCenter);
    graph_qTmean_vs_qTbinCenter_mc->SetPoint(qTbin - 1, qTbinCenter, qTav_mc/qTbinCenter);
    graph_qTmean_vs_qTbinCenter_mc->SetPointError(qTbin - 1, 0.5*qTbinWidth, 0.5*qTbinWidth, qTavErr_mc/qTbinCenter, qTavErr_mc/qTbinCenter);
  }

  // reset x-axis range selection 
  histogram_qT_data->GetXaxis()->SetRange(1., 0.);
  histogram_qT_mc_central->GetXaxis()->SetRange(1., 0.);

  TString outputFileName = Form("plots/makeUnclusteredEnergyPlots_%s_qTbinCenterCorrections", type_string.data());
  outputFileName.ReplaceAll(".", "_");
  outputFileName.Append(".png");
  showGraphs(800, 900,
	     graph_qTmean_vs_qTbinCenter_data, "Data",
	     graph_qTmean_vs_qTbinCenter_mc, "Simulation",
	     0, "",
	     0, "",
	     0, "",
	     0, "",
	     0, "",
	     0., 300., 10, "q_{T} / GeV", 1.2,
	     false, 0.8, 1.7, 0.975, 1.025, "<q_{T}>/q_{T}^{bin-center}", 1.2, 
	     0.61, 0.69, 
	     outputFileName.Data());

  TH1* histogram_zVtx_data = getHistogram(inputFile, dqmDirectories["data"], "zVtx");
  TH1* histogram_zVtx_mc_central = getHistogram(inputFile, dqmDirectories["mc_central"],      "zVtx");
  showDistribution(800, 900,
		   histogram_zVtx_data, "Data",
		   histogram_zVtx_mc_central, "Simulation",
		   0, "",
		   0, "",
		   0, "",
		   0, "",
		   -25.0, +25.0, 500, "z_{Vtx} / cm", 1.2, 
		   true, 1.e-6, 1.e+1, -0.10, +0.10, "a.u.", 1.2,
		   0.61, 0.72, 
		   Form("plots/makeUnclusteredEnergyPlots_%s_zVtx.png", type_string.data()));

  TFile* outputFile = new TFile("zVtxReweight_runs203894to208686_vs_Summer12mc.root", "RECREATE");
  histogram_zVtx_data->Scale(1./histogram_zVtx_data->Integral());
  histogram_zVtx_mc_central->Scale(1./histogram_zVtx_mc_central->Integral());
  TH1* histogram_zVtx_reweight = compRatioHistogram("zVtxReweight", histogram_zVtx_data, histogram_zVtx_mc_central, false);
  histogram_zVtx_reweight->Write();
  delete outputFile;

  std::ofstream* outputFile_txt = new std::ofstream(Form("unclEnResidualCorr_2012RunABCDruns190456to208686_%s.txt", type_string.data()), std::ios::out);
  for ( std::vector<std::string>::const_iterator line = residualCorrections.begin();
	line != residualCorrections.end(); ++line ) {
    (*outputFile_txt) << (*line) << std::endl;
  }
  delete outputFile_txt;

  delete inputFile;
}
