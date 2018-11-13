
#ifndef CERNfunctions_H
#define CERNfunctions_H

#if !defined (__CINT__) || defined (__MAKECINT__)
#include "Rtypes.h"
#endif

extern "C" Float_t denlan_(Float_t *x);

#endif



#include <TH1D.h>
#include <TStyle.h>
#include <TString.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TRandom3.h>
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "Math/Integrator.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/Functor.h"
#include "Math/GaussIntegrator.h"
#include "TSpline.h"


//C, C++
#include <string>
#include <sstream>
#include <iostream>
#include <cctype>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

double GlobalParamSum = -999;
TSpline3 GlobalSpline;



struct FitedParams{
  string FileName;
  Double_t MaxPosition;
  double Sigma;
  double SigmaError;
  double Angle;

  double ASU;
  double MP;
  double MP2;
  double MPError;
  double MP2Error;

  bool SigmaFix;
  double Median;

  double Width;
  double WidthError;
  double Width2;
  double Width2Error;

  double Area;
  double Area2;


  double ErrCentre;
  double ErrCentreError;
  double ErrWidth;
  double ErrWidthError;

  double AngleCoef;
  double AngleCoefError;


  double Pedestal;


};

struct InputData{
  double ADC;
  double RealADC;
  int number;
  int bx;
  string name;
  int sca;
  int channel;
};




int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

bool CompareFunction (InputData i, InputData j)
{
  int fun1 =  i.channel*100 + i.sca;
  int fun2 =  j.channel*100 + j.sca;


  return (fun1 < fun2);

}



Double_t landaufun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //par[4]=Width 2 (scale) parameter of Landau density
   //par[5]=Most Probable 2 (MP, location) parameter of Landau density
   //par[6]=Total area 2 (integral -inf to inf, normalization constant)
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t mpc2;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      mpc = par[1] - mpshift * par[0];
      mpc2 = par[5] - mpshift * par[4];

      sum = par[2]*TMath::Landau(x[0],mpc,par[0]) / par[0] + par[6]*TMath::Landau(x[0],mpc2,par[4]) / par[4];


      return ( sum * invsq2pi );
}





Double_t langaufun(Double_t *x, Double_t *par) {
  // Control constants
  /*
  Double_t np = 60;      // number of convolution steps
  Double_t sc =   3.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step ;
  Double_t i;
  xupp = 600;
  xlow = 0;
  // Range of convolution integral


  step = (xupp - x[0]) / np;

  for(i=1; i <= np; i++) {
     xx = x[0] + i*step;
     fland =landaufun(&xx,par);
     sum += fland/(xx-xlow)*par[7]*step;

  }






   sum += landaufun(&x[0],par)*(1 - par[7]);

  return sum*((TMath::Erf ((x[0] - par[8])/par[9])+1)/2);
  */
  double paramsum = 0;
  double answer = 0;
  //cout << "Calculate Sum" << endl;
  for (int kk = 0; kk < 10; kk++){
    paramsum += par[kk];
    //cout << kk << " " << par[kk] << endl;
  }

  if (paramsum == GlobalParamSum){
       answer = (GlobalSpline.Eval(x[0])*(TMath::Erf ((x[0] - par[8])/par[9])+1)/2);
       //cout << "Paras Stay" << endl;
  }

  if (paramsum != GlobalParamSum){
  int xlow = 0;
  int xupp = 400;
  double digitisation = 1;
  double step = 1/digitisation;

   TH2F *LandauHist = new TH2F("landauHist", "LandauHist", digitisation*(xupp-xlow) , xlow , xupp, digitisation*(xupp-xlow) , xlow , xupp);


   for (int i = 0; i < digitisation*(xupp-xlow) ; i ++ ){
   double xx = (i-0.5) * step;
   double landau = landaufun(&xx,par);
    LandauHist->SetBinContent(i, i, landau*(1 - par[7]));
    for (int j = 0; j < i; j ++){
       LandauHist->SetBinContent(i, j, landau*par[7]/i);
    }
  }

  TH1D * LandauHist1D = LandauHist->ProjectionY();


    TH2F *LanGaus = new TH2F("LanGauss", "Langauss", digitisation*(xupp-xlow) , xlow , xupp, digitisation*(xupp-xlow) , xlow , xupp);

    for (int i = 0; i < digitisation*(xupp-xlow) ; i ++ ){

      double xx = (i-0.5) * step;
      double landau = LandauHist1D->GetBinContent(i) ;
      double sigma = sqrt(par[3]*par[3]*xx + 9);
      double nsigm = 5;
      int jmin = 0;
      int jmax = digitisation*(xupp-xlow);
      double jstep = step;

      if (jmin*jstep < (xx - sigma*nsigm)) {
        jmin = (xx - sigma*nsigm)/jstep;
      }

      if (jmax*jstep > (xx + sigma*nsigm)) {
        jmax = (xx + sigma*nsigm)/jstep;
      }

      for (int j = jmin; j < jmax; j ++){
        double xy = (j-0.5)*jstep;
        Float_t vrem = landau*TMath::Gaus(xy,xx,sigma,kTRUE)*jstep;
        LanGaus->SetBinContent(i, j, vrem);
      }
    }

     TH1D *LanGaus1D = LanGaus->ProjectionY();







   TSpline3 s =  TSpline3(LanGaus1D);
   delete LanGaus1D;
   delete LanGaus;
   delete LandauHist1D;
   delete LandauHist;
   GlobalSpline = s;

   answer = (s.Eval(x[0])*(TMath::Erf ((x[0] - par[8])/par[9])+1)/2);

   GlobalParamSum = paramsum;

   //cout << "Paras Changed " << GlobalParamSum << answer << endl;


  //cout << "here" << " ";

}
return answer;

}

Double_t langaufunAng(Double_t *x, Double_t *par) {

  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //par[4]=Width 2 (scale) parameter of Landau density
  //par[5]=Most Probable 2 (MP, location) parameter of Landau density
  //par[6]=Total area 2 (integral -inf to inf, normalization constant)
  //In the Landau distribution (represented by the CERNLIB approximation),
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.


      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;
      Double_t sigma = sqrt(par[3]*par[3]*x[0] + 9);


      // MP shift correction
      mpc = par[1] - mpshift * par[0];


      // Range of convolution integral

      xlow = x[0] - sc * sigma;
      xupp = x[0] + sc * sigma;

      xlow = xlow * (xlow > 0);

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum

      for(i=1.0; i<=np/2; i++) {

         xx = xlow + (i-.5) * step;
         fland = langaufunAng(&xx,par);
         sum += fland * TMath::Gaus(xx,x[0],sqrt(par[3]*par[3]*xx + 9));

         xx = xupp - (i-.5) * step;
         fland = langaufunAng(&xx,par);
         sum += fland * TMath::Gaus(xx,x[0],sqrt(par[3]*par[3]*xx + 9));
      }
      //return (par[2] * step * sum * invsq2pi / par[3]);
      return sum*step/sigma*(TMath::Erf ((x[0] - par[8])/par[9])+1)/2;
}




TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //par[4]=Width 2 (scale) parameter of Landau density
  //par[5]=Most Probable 2 (MP, location) parameter of Landau density
  //par[6]=Total area 2 (integral -inf to inf, normalization constant)
  //In the Landau distribution (represented by the CERNLIB approximation),
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

   // Variables for langaufit call:
   //   his             histogram to fit
   //   fitrange[2]     lo and hi boundaries of fit range
   //   startvalues[4]  reasonable start values for the fit
   //   parlimitslo[4]  lower parameter limits
   //   parlimitshi[4]  upper parameter limits
   //   fitparams[4]    returns the final fit parameters
   //   fiterrors[4]    returns the final fit errors
   //   ChiSqr          returns the chi square
   //   NDF             returns ndf

   Int_t i;
   Char_t FunName[100];

   sprintf(FunName,"Fitfcn_%s",his->GetName());

   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   if (ffitold) delete ffitold;

   TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],10);
   cout << "Eval" << endl;
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MP","Area","GSigma","Width2","MP2", "Area2", "AngleCoef", "Errentre", "ErrWidth");

   for (i=0; i<10; i++) {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }

   his->Fit(FunName,"RBL0");   // fit within specified range, use ParLimits, do not plot
   cout << "Fitted" << endl;

   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<10; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf

   return (ffit);              // return fit function

}


Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {

   // Seaches for the location (x value) at the maximum of the
   // Landau-Gaussian convolute and its full width at half-maximum.
   //
   // The search is probably not very efficient, but it's a first try.

   Double_t p,x,fy,fxr,fxl;
   Double_t step;
   Double_t l,lold;
   Int_t i = 0;
   Int_t MAXCALLS = 10000;


   // Search for maximum

   p = params[1] - 0.1 * params[0];
   step = 0.05 * params[0];
   lold = -2.0;
   l    = -1.0;


   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = langaufun(&x,params);

      if (l < lold)
         step = -step/10;

      p += step;
   }

   if (i == MAXCALLS)
      return (-1);

   maxx = x;

   fy = l/2;


   // Search for right x location of fy

   p = maxx + params[0];
   step = params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;


   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);

      if (l > lold)
         step = -step/10;

      p += step;
   }

   if (i == MAXCALLS)
      return (-2);

   fxr = x;


   // Search for left x location of fy

   p = maxx - 0.5 * params[0];
   step = -params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;

   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);

      if (l > lold)
         step = -step/10;

      p += step;
   }

   if (i == MAXCALLS)
      return (-3);


   fxl = x;

   FWHM = fxr - fxl;
   return (0);
}

FitedParams FitMacro(vector<InputData> InputVector, int event, TFile *MyFile, Int_t max, Int_t min, int chip, string basic, double sigma, double width ){
   FitedParams output;

   string angle = "angle";
   string ASU = "ASU";
   string DIF = "dif";
   output.Angle = 0;

   for (unsigned int i=0; i < basic.size() - 6; i++){
     if(basic[i] == ASU[0] && basic[i+1] == ASU[1] && basic[i+2] == ASU[2]){
           string temp;
           temp += basic[i+3];
           output.ASU = atof(temp.c_str());
     }
     if(basic[i] == DIF[0] && basic[i+1] == DIF[1] && basic[i+2] == DIF[2]){
           string temp;
           temp += basic[i+4];
           temp += basic[i+5];
           temp += basic[i+6];
           temp += basic[i+7];
           temp += basic[i+8];

           output.ASU = atof(temp.c_str());
     }

      if(basic[i] == angle[0] && basic[i+1] == angle[1] && basic[i+2] == angle[2]){
        string temp;
        temp += basic[i+5];
        if (isdigit(basic[i+6])) temp += basic[i+6];
        output.Angle = atof(temp.c_str());
      }
     //else output.Angle = 0;

   }

   int realchip = output.ASU*16 - 6;

   vector<double> MedianSerch;

for (int i = 0; i < event; i++){
  if(InputVector[i].ADC >= 30 && InputVector[i].ADC < 80) MedianSerch.push_back(InputVector[i].ADC);
}

sort(MedianSerch.begin(), MedianSerch.end());

if(MedianSerch.size() > 0){
output.Median = MedianSerch[MedianSerch.size()/2];
}

   TH1F *h1 = new TH1F(basic.c_str(), basic.c_str(), (max-min) , min , max );

vector<double> MedianPedestal;

   for (int i = 0; i < event; i++){
    if (InputVector[i].number == realchip && ((output.ASU != 8) || ( InputVector[i].bx > 350))){
      MedianPedestal.push_back( InputVector[i].RealADC - InputVector[i].ADC);

    }
   if (InputVector[i].number == realchip && InputVector[i].ADC >= min && InputVector[i].ADC < max && ((output.ASU != 8) || ( InputVector[i].bx > 350))){
     h1 ->Fill(InputVector[i].ADC);

   }
 }

 sort (MedianPedestal.begin(), MedianPedestal.end());


if(MedianPedestal.size() > 0){
output.Pedestal = MedianPedestal[MedianPedestal.size()/2];
}



     Double_t fr[2];
     Double_t sv[10], pllo[10], plhi[10], fp[10], fpe[10];
     fr[0]=10;
     fr[1]=190/cos(output.Angle/180*3.14159);
     cout <<"Angle = " << output.Angle << endl;
     cout <<"ASU = " << output.ASU << endl;
     //Fit parameters:
     //par[0]=Width (scale) parameter of Landau density
     //par[1]=Most Probable (MP, location) parameter of Landau density
     //par[2]=Total area (integral -inf to inf, normalization constant)
     //par[3]=Width (sigma) of convoluted Gaussian function
     //
     //par[4]=Width 2 (scale) parameter of Landau density
     //par[5]=Most Probable 2 (MP, location) parameter of Landau density
     //par[6]=Total area 2 (integral -inf to inf, normalization constant)
     //In the Landau distribution (represented by the CERNLIB approximation),
     //the maximum is located at x=-0.22278298 with the location parameter=0.
     //This shift is corrected within this function, so that the actual
     //maximum is identical to the MP parameter.

     // Variables for langaufit call:
     //   his             histogram to fit
     //   fitrange[2]     lo and hi boundaries of fit range
     //   startvalues[4]  reasonable start values for the fit
     //   parlimitslo[4]  lower parameter limits
     //   parlimitshi[4]  upper parameter limits
     //   fitparams[4]    returns the final fit parameters
     //   fiterrors[4]    returns the final fit errors
     //   ChiSqr          returns the chi square
     //   NDF             returns ndf


     if (sigma == 0 && width == 0){
     pllo[0]=1; pllo[1]=1/cos(output.Angle/180*3.14159); pllo[2]=h1->GetEntries()/1000; pllo[3]=0.1;
     plhi[0]=20; plhi[1]=300/cos(output.Angle/180*3.14159); plhi[2]=h1->GetEntries()*2000; plhi[3]=10;
     sv[0]=6*(1/cos(output.Angle/180*3.14159)); sv[1]=64/cos((output.Angle+2)/180*3.14159) - output.ASU/1.5; sv[2]=h1->GetEntries()*2.1; sv[3]=0.8;
     output.SigmaFix = 0;
     }

     else if (sigma != 0 && width != 0){
     pllo[0]=width - 0.001; pllo[1]=1; pllo[2]=h1->GetEntries()/1000; pllo[3]=sigma + 0.001;
     plhi[0]=width + 0.001; plhi[1]=300; plhi[2]=h1->GetEntries()*2000; plhi[3]=sigma + 0.001;
     sv[0]=width; sv[1]=60/cos(output.Angle/180*3.14159); sv[2]=h1->GetEntries()/2; sv[3]=sigma;
      output.SigmaFix = 1;
     }

     else if (sigma != 0 && width == 0){
     pllo[0]=1 ; pllo[1]=1; pllo[2]=h1->GetEntries()/1000; pllo[3]=sigma + 0.001;
     plhi[0]=200 ; plhi[1]=300; plhi[2]=h1->GetEntries()*2000; plhi[3]=sigma + 0.001;
     sv[0]=10; sv[1]=60/cos(output.Angle/180*3.14159); sv[2]=h1->GetEntries()/2; sv[3]=sigma;
      output.SigmaFix = 1;
     }

     pllo[4] = 2; pllo[5] = 2.1*40/cos(output.Angle/180*3.14159); pllo[6] = h1->GetEntries()*0.001;
     plhi[4] = 40; plhi[5] = 2.1*70/cos(output.Angle/180*3.14159); plhi[6] = h1->GetEntries()*100;
     sv[4] =6*(2/cos(output.Angle/180*3.14159)); sv[5] = 2.1*64/cos(output.Angle/180*3.14159); sv[6] = sv[2]/6;

     sv[7] = 0.1/cos(output.Angle/180*3.14159);
     pllo[7] = 0;
     plhi[7] = 0.5;

     sv[8] = 35;
     pllo[8] = 20;
     plhi[8] = 80;

     if (output.Angle == 1){
        sv[8] = 38;
        if(output.ASU == 4) sv[8] = 49;
        if(output.ASU == 8) sv[8] = 42;
        pllo[8] = sv[8] - 1;
        plhi[8] = sv[8] + 1;

     }


     sv[9] = 10;
     pllo[9] = 1;
     plhi[9] = 20;


     /*
     sv[7] = 0.11;
     pllo[7] = 0.00;
     plhi[7] = 0.5;

     sv[8] = 140/cos(output.Angle/180*3.14159);
     pllo[8] =  pllo[1]*2;
     plhi[8] =  plhi[1]*2;
   */


     Double_t chisqr;
     Int_t    ndf;
     TF1 *fitsnr = langaufit(h1,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);


     Double_t SNRPeak, SNRFWHM;
     langaupro(fp,SNRPeak,SNRFWHM);

     printf("Fitting done\nPlotting results...\n");

     h1->SetBinErrorOption(TH1::kPoisson);
     cout << "GlobalVarSet" << endl;
     //GlobalParamSum = -9999;



     // Global style settings
     gStyle->SetOptStat(1111);
     gStyle->SetOptFit(111);
     gStyle->SetLabelSize(0.03,"x");
     gStyle->SetLabelSize(0.03,"y");



     output.MaxPosition = fitsnr->GetMaximumX(-100, 600);
     output.FileName = basic;
     output.Sigma = fp[3];
     output.SigmaError = fpe[3];

     output.MP = fp[1];
     output.MPError = fpe[1];
     output.MP2 = fp[5];
     output.MP2Error = fpe[5];

     output.Width = fp[0];
     output.WidthError = fpe[0];
     output.Width2 = fp[4];
     output.Width2Error = fpe[4];

     output.Area = fp[2];
     output.Area2 = fp[6];

     output.AngleCoef = fp[7];
     output.AngleCoefError = fpe[7];

     output.ErrCentre = fp[8];
     output.ErrCentreError = fpe[8];

     output.ErrWidth = fp[9];
     output.ErrWidthError = fpe[9];





     h1->GetListOfFunctions()->Add(fitsnr);




     h1->Draw();


   if ( MyFile->IsOpen() ){
     h1->Write();

   }

   return output;
}


int main(int argc, char *argv[]){


vector<InputData> InputVector;

InputData InputVectorElement;


    ifstream input("MIP_all.dat"); //put your program together with thsi file in the same folder.

    if(input.is_open()){
      string just;
        getline(input, just);
        getline(input, just);
        while(!input.eof()){

              string str;
              int data;
              getline(input, str); //read number

              int counter = 0;

              string str1 = " ";
              string str2 = " ";
              string str3 = "\"";
              string str4 = ".";
              string str5 = "-";
              double swap[100];
              vector<string> buffer;
              int i = 0;

              while ( i < str.size() ){

                if((str[i] == str3[0] || str[i] == str2[0]) && i < str.size() -1 ){
                  i++;
                  string temporary;
                  while (str[i] != str1[0] && i < str.size() ){
                    if (str[i] != str3[0]) temporary += str[i];
                    i++;
                  }
                  buffer.push_back(temporary);
                 }


                 if (str[i] != str1[0] && i < str.size()) i++;

                 }


                 double num = 0;
                 if (buffer.size() > 20){
                 string temp = buffer[0];
                 cout << buffer[0] << endl;
                 InputVectorElement.name = temp;


                 temp = buffer[6];
                 num = atof(temp.c_str());

                 InputVectorElement.bx = num;


                 temp = buffer[7];
                 num = atof(temp.c_str());

                 InputVectorElement.RealADC = num;


                 temp = buffer[1];
                 num = atof(temp.c_str());

                 InputVectorElement.number = num;

                 temp = buffer[4];
                 num = atof(temp.c_str());

                 InputVectorElement.sca = num;



                 temp = buffer[2];
                 num = atof(temp.c_str());

                 InputVectorElement.channel = num;


                 temp = buffer[buffer.size()-1];
                 num = atof(temp.c_str());
                 InputVectorElement.ADC = num;


                 InputVector.push_back(InputVectorElement);

               }

                 }


       }
  /*
   for (int i = 0; i < number.size(); i++){

   }
*/



TFile *MyFile = new TFile("MIP.root","NEW");

vector<FitedParams>  FitParVector;

int NumOfFiles = 0;

while (InputVector.size()>10){

GlobalParamSum = -999.0;
int event = 0;
while (InputVector[event].name == InputVector[event+1].name && (event < InputVector.size()-1)){
event++;

}
event ++;



sort (InputVector.begin(), InputVector.begin() + event, CompareFunction);
cout << event << endl;

/*
event = 0;
while (InputVector[event].channel == InputVector[event+1].channel && (event < InputVector.size()-1)){
event++;

}
event ++;
*/
stringstream ss1;
ss1 << InputVector[0].name << "Ch" << InputVector[0].channel;
string str1 = ss1.str();

string basic = str1;
Int_t max = 600;
Int_t min = 0;


cout << "max " << max  << endl;
cout << "min " << min << endl;


int chip = 122;



FitedParams result = FitMacro(InputVector, event, MyFile, max,  min, chip, basic, 0, 0);
cout << "1" << endl;
FitParVector.push_back(result);



InputVector.erase(InputVector.begin(), InputVector.begin() + event);

NumOfFiles ++;

}


Double_t Angle[1000];
Double_t Sigma[1000];
Double_t SigmaError[1000];
Double_t SigmaCos[1000];
Double_t ASU[1000];
Double_t ASUError[1000];
Double_t mp[1000];
Double_t mperror[1000];
Double_t mp2[1000];
Double_t mp2error[1000];

Double_t width[1000];
Double_t widtherror[1000];
Double_t width2[1000];
Double_t width2error[1000];

Double_t area[1000];
Double_t area2[1000];
Double_t areas[1000];

Double_t errorcentre[1000];
Double_t errorcentreerror[1000];
Double_t errorwidth[1000];
Double_t errorwidtherror[1000];

Double_t anglecoef[1000];
Double_t anglecoeferror[1000];

Double_t avg_sigma = 0;
Double_t avg_width = 0;

Double_t mpFix[1000];
Double_t mperrorFix[1000];
Double_t ASUFix[1000];
Double_t ASUErrorFix[1000];
Double_t Median[1000];
Double_t Pedestal[1000];




int CountFree = 0;
int CountFix = 0;

ofstream out;
out.open ("output.txt");


for (int i = 0; i < FitParVector.size(); i++){
  if (FitParVector[i].SigmaFix == 0 &&  FitParVector[i].MP > 10 && FitParVector[i].Area > 20){// && FitParVector[i].ASU != 8){
  cout << "Angle = " << FitParVector[i].Angle <<" MaxX = " << FitParVector[i].MaxPosition << " Sigma  = " << FitParVector[i].Sigma ; //<< " SigmaError = " << FitParVector[i].SigmaError ;
  cout << endl;
  Angle[CountFree] = FitParVector[i].Angle;

  ASU[CountFree] = FitParVector[i].ASU;
  ASUError[CountFree] = 0;
  Sigma[CountFree] = FitParVector[i].Sigma;
  SigmaError[CountFree] = FitParVector[i].SigmaError;


  mp[CountFree] = FitParVector[i].MP;
  mperror[CountFree] = FitParVector[i].MPError;
  mp2[CountFree] = FitParVector[i].MP2;
  mp2error[CountFree] = FitParVector[i].MP2Error;

  width[CountFree] = FitParVector[i].Width;
  widtherror[CountFree] = FitParVector[i].WidthError;
  width2[CountFree] = FitParVector[i].Width2;
  width2error[CountFree] = FitParVector[i].Width2Error;

  Median[CountFree] = FitParVector[i].Median;
  Pedestal[CountFree] = FitParVector[i].Pedestal;


  area[CountFree] = FitParVector[i].Area;
  area2[CountFree] = FitParVector[i].Area2;
  areas[CountFree] = area2[CountFree] / area[CountFree];

  errorcentre[CountFree] = FitParVector[i].ErrCentre;
  errorcentreerror[CountFree] = FitParVector[i].ErrCentreError;
  errorwidth[CountFree] = FitParVector[i].ErrWidth;
  errorwidtherror[CountFree] = FitParVector[i].ErrWidthError;

  anglecoef[CountFree] = FitParVector[i].AngleCoef;
  anglecoeferror[CountFree] = FitParVector[i].AngleCoefError;


  SigmaCos[CountFree] = Sigma[CountFree]*cos(Angle[CountFree]/360*2*3.14159);

  avg_sigma += FitParVector[i].Sigma;
  avg_width += FitParVector[i].Width;

  out << mp[CountFree] << " " << mperror[CountFree] << " " << width[CountFree]<< " " << widtherror[CountFree] << endl;
  out << mp2[CountFree] << " " << mp2error[CountFree] << " " << width2[CountFree]<< " " << width2error[CountFree] << endl;
  out << anglecoef[CountFree] << " " << anglecoeferror[CountFree] << " " << errorcentre[CountFree] << " " << errorcentreerror[CountFree] << " " << errorwidth[CountFree] << " " << errorwidtherror[CountFree] << endl;
  out << areas[CountFree] << " " << Sigma[CountFree] << " " << SigmaError[CountFree] << " " << Angle[CountFree] << " " << ASU[CountFree] << endl;
  CountFree ++;
  }
  if (FitParVector[i].SigmaFix == 1){

    mpFix[CountFix] = FitParVector[i].MP;
    mperrorFix[CountFix] = FitParVector[i].MPError;
    ASUErrorFix[CountFix] = 0;
    ASUFix[CountFix] = FitParVector[i].ASU;
    CountFix ++;
  }
}
out.close();


  TGraph* sigma = new TGraph(CountFree,Angle,Sigma);

  sigma->SetLineWidth(0);
  sigma->SetMarkerColor(4);
  sigma->SetMarkerSize(0.5);
  sigma->SetMarkerStyle(21);

  sigma->SetTitle("Sigma");
  sigma->SetName("Sigma");

  sigma->GetXaxis()->SetTitle("Angle (Degree)");
  sigma->GetYaxis()->SetTitle("Sigma");

  sigma->Draw("ap");



  TGraph* sigmaCos = new TGraph(CountFree,Angle,SigmaCos);


  sigmaCos->SetLineWidth(0);
  sigmaCos->SetMarkerColor(4);
  sigmaCos->SetMarkerSize(0.5);
  sigmaCos->SetMarkerStyle(21);

  sigmaCos->SetTitle("Sigma*Cos(angle)");
  sigmaCos->SetName("Sigma*Cos(angle)");

  sigmaCos->GetXaxis()->SetTitle("Angle (Degree)");
  sigmaCos->GetYaxis()->SetTitle("Sigma*Cos(angle)");

  sigmaCos->Draw("ap");

  TGraph* MedianASU = new TGraph(CountFree,ASU,Median);


  MedianASU->SetLineWidth(0);
  MedianASU->SetMarkerColor(4);
  MedianASU->SetMarkerSize(0.5);
  MedianASU->SetMarkerStyle(21);

  MedianASU->SetTitle("MedianASU");
  MedianASU->SetName("MedianASU");

  MedianASU->GetXaxis()->SetTitle("ASU");
  MedianASU->GetYaxis()->SetTitle("Median");

  MedianASU->Draw("ap");


  TGraph* PedestalAverage = new TGraph(CountFree,ASU,Pedestal);


  PedestalAverage->SetLineWidth(0);
  PedestalAverage->SetMarkerColor(4);
  PedestalAverage->SetMarkerSize(0.5);
  PedestalAverage->SetMarkerStyle(21);

  PedestalAverage->SetTitle("Pedestal");
  PedestalAverage->SetName("Pedestal");

  PedestalAverage->GetXaxis()->SetTitle("ASU");
  PedestalAverage->GetYaxis()->SetTitle("Pedestal");

  PedestalAverage->Draw("ap");

  TGraph* Areas = new TGraph(CountFree,ASU,areas);


  Areas->SetLineWidth(0);
  Areas->SetMarkerColor(4);
  Areas->SetMarkerSize(0.5);
  Areas->SetMarkerStyle(21);

  Areas->SetTitle("Area2/Area1");
  Areas->SetName("Areas");

  Areas->GetXaxis()->SetTitle("ASU");
  Areas->GetYaxis()->SetTitle("Area2/Area1");

  Areas->Draw("ap");



  TGraphErrors* MPerr = new TGraphErrors(CountFree,ASU,mp,ASUError,mperror);

  MPerr->SetLineWidth(0);
  MPerr->SetMarkerColor(4);
  MPerr->SetMarkerSize(0.5);
  MPerr->SetMarkerStyle(21);

  MPerr->SetTitle("MP ASU");
  MPerr->SetName("MP ASU");

  MPerr->GetXaxis()->SetTitle("ASU");
  MPerr->GetYaxis()->SetTitle("MP");

  TGraphErrors* Centre = new TGraphErrors(CountFree,ASU,errorcentre,ASUError,errorcentreerror);

  Centre->SetLineWidth(0);
  Centre->SetMarkerColor(4);
  Centre->SetMarkerSize(0.5);
  Centre->SetMarkerStyle(21);

  Centre->SetTitle("ErfCentre");
  Centre->SetName("ErfCentre");

  Centre->GetXaxis()->SetTitle("ASU");
  Centre->GetYaxis()->SetTitle("ErfCentre");

  TGraphErrors* ErrWidth = new TGraphErrors(CountFree,ASU,errorwidth,ASUError,errorwidtherror);

  ErrWidth->SetLineWidth(0);
  ErrWidth->SetMarkerColor(4);
  ErrWidth->SetMarkerSize(0.5);
  ErrWidth->SetMarkerStyle(21);

  ErrWidth->SetTitle("ErfWidth");
  ErrWidth->SetName("ErfWidth");

  ErrWidth->GetXaxis()->SetTitle("ASU");
  ErrWidth->GetYaxis()->SetTitle("ErfWidth");

  TGraphErrors* AngleCoef = new TGraphErrors(CountFree,ASU,anglecoef,ASUError,anglecoeferror);

  AngleCoef->SetLineWidth(0);
  AngleCoef->SetMarkerColor(4);
  AngleCoef->SetMarkerSize(0.5);
  AngleCoef->SetMarkerStyle(21);

  AngleCoef->SetTitle("AngleCoef");
  AngleCoef->SetName("AngleCoef");

  AngleCoef->GetXaxis()->SetTitle("ASU");
  AngleCoef->GetYaxis()->SetTitle("AngleCoef");

  TGraphErrors* Width = new TGraphErrors(CountFree,ASU,width,ASUError,widtherror);

  Width->SetLineWidth(0);
  Width->SetMarkerColor(4);
  Width->SetMarkerSize(0.5);
  Width->SetMarkerStyle(21);

  Width->SetTitle("Width");
  Width->SetName("Width");

  Width->GetXaxis()->SetTitle("ASU");
  Width->GetYaxis()->SetTitle("Width");

  TGraphErrors* Width2 = new TGraphErrors(CountFree,ASU,width2,ASUError,width2error);

  Width2->SetLineWidth(0);
  Width2->SetMarkerColor(4);
  Width2->SetMarkerSize(0.5);
  Width2->SetMarkerStyle(21);

  Width2->SetTitle("Width2");
  Width2->SetName("Width2");

  Width2->GetXaxis()->SetTitle("ASU");
  Width2->GetYaxis()->SetTitle("Width2");

  TGraphErrors* MP2err = new TGraphErrors(CountFree,ASU,mp2,ASUError,mp2error);

  MP2err->SetLineWidth(0);
  MP2err->SetMarkerColor(4);
  MP2err->SetMarkerSize(0.5);
  MP2err->SetMarkerStyle(21);

  MP2err->SetTitle("MP2 ASU");
  MP2err->SetName("MP2 ASU");

  MP2err->GetXaxis()->SetTitle("ASU");
  MP2err->GetYaxis()->SetTitle("MP2");




  TGraphErrors* MPerrFix = new TGraphErrors(CountFree,ASUFix,mpFix,ASUErrorFix,mperrorFix);

  MPerrFix->SetLineWidth(0);
  MPerrFix->SetMarkerColor(4);
  MPerrFix->SetMarkerSize(0.5);
  MPerrFix->SetMarkerStyle(21);

  MPerrFix->SetTitle("MP ASU Sigma Fix");
  MPerrFix->SetName("MP ASU Sigma Fix");

  MPerrFix->GetXaxis()->SetTitle("ASU");
  MPerrFix->GetYaxis()->SetTitle("MP");



  TGraphErrors* SigmaASU = new TGraphErrors(CountFree,ASU,Sigma,ASUError,SigmaError);

  SigmaASU->SetLineWidth(0);
  SigmaASU->SetMarkerColor(4);
  SigmaASU->SetMarkerSize(0.5);
  SigmaASU->SetMarkerStyle(21);

  SigmaASU->SetTitle("Sigma ASU");
  SigmaASU->SetName("Sigma ASU");

  SigmaASU->GetXaxis()->SetTitle("ASU");
  SigmaASU->GetYaxis()->SetTitle("Sigma");

  cout << "write" << endl;


  PedestalAverage->Write();

  Areas->Write();

  MedianASU->Write();
  SigmaASU->Write();

  MPerr->Write();
  MP2err->Write();

  Width->Write();
  Width2->Write();
  AngleCoef->Write();
  ErrWidth->Write();
  Centre->Write();

  MPerrFix->Write();

  sigmaCos->Write();
  sigma->Write();

  cout << "Average Sigma = "<< avg_sigma/NumOfFiles << endl;
  cout << "Average Width = "<< avg_width/NumOfFiles << endl;




  //TF1 *convolution = new TF1("Convolution",Conv,fr[0],fr[1],7);
  //convolution->SetParameters(parm);

  //convolution->Draw();
  //convolution->Write();


  MyFile->Close();
  return 0;



}
