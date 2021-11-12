#include "greenfunction.h"


//__________________________________________________________________________________________
//##########################################################################################
//  Class constructor
//    Create empty GreensFunctions
//##########################################################################################
GreensFunctions::GreensFunctions()
{

}

GreensFunctions::GreensFunctions(string backgroundAttractorFile, string greensfunctionsfile, int backgroundPoints, int greensFunctionsPoints, int greensFunctionsChuncks, double cInfinity, double etaOverS, double tauHydro)
{
  //  Assinging input parameters
  background_attractor_file = backgroundAttractorFile;
  greens_functions_file = greensfunctionsfile;

  background_points = backgroundPoints;
  greens_functions_points = greensFunctionsPoints;
  greens_functions_chuncks = greensFunctionsChuncks;
  c_infinity = cInfinity;
  eta_over_s = etaOverS;
  tau_hydro = tauHydro;

  //  Constructing nuEff for background attractor
  double Nc = 3.0;
  double Nf = 3.0;
  double nuG = 2.0*(Nc*Nc - 1.0);
  double nuQ = 2.0*Nc*Nf;
  nuEff = nuG + 7.0/4.0*nuQ;

  SetupBackgroundAttractor();
  SetupGreensFunctions();
//  cout << "Greens Function setup complete" << endl;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Class deconstructor
//##########################################################################################
GreensFunctions::~GreensFunctions()
{

}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Implicit Copy
//##########################################################################################
GreensFunctions::GreensFunctions(const GreensFunctions &original)
{
  CopyGreensFunctions(original);
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  GreensFunctions Copy Function
//##########################################################################################
void GreensFunctions::CopyGreensFunctions(const GreensFunctions &e)
{
  background_attractor_file = e.background_attractor_file;
  greens_functions_file = e.greens_functions_file;

  background_points = e.background_points;
  greens_functions_points = e.greens_functions_points;
  greens_functions_chuncks = e.greens_functions_chuncks;
  c_infinity = e.c_infinity;
  eta_over_s = e.eta_over_s;
  tau_hydro = e.tau_hydro;
//cout << "c_infinity " << c_infinity << endl;
  nuEff = e.nuEff;

  // EVALUATE BETWEEN dXdTMin,dXdTMax AND wTMin,wTMax //
  dXdTMin = e.dXdTMin;
  dXdTMax = e.dXdTMax;
  wTMin = e.wTMin;
  wTMax = e.wTMax;

  // GSL INTERPOLATION OBJECTS //
  FsswTAcc = e.FsswTAcc;
  GsswTAcc = e.GsswTAcc;
  FssdXdTAcc = e.FssdXdTAcc;
  GssdXdTAcc = e.GssdXdTAcc;
  FssInt = e.FssInt;
  GssInt = e.GssInt;

  wTValues = e.wTValues;
  dXdTValues = e.dXdTValues;

  FssValues = e.FssValues;
  GssValues = e.GssValues;

  // GSL INTERPOLATION OBJECTS //
  EAcc = e.EAcc;
  EInt = e.EInt;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Overide operator=
//##########################################################################################
GreensFunctions& GreensFunctions::operator= (const GreensFunctions& original)
{
	CopyGreensFunctions(original);
	return *this;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Read in data for background attractor and setup relevant variables
//##########################################################################################
void GreensFunctions::SetupBackgroundAttractor()
{
  // SET DATA //
  vector<double> wTildeValues;
  vector<double> EValues;

  ifstream InStream;
  InStream.open(background_attractor_file);

  int i = 0;

  while(!InStream.eof())
  {
      double wT;
      double EVal;

      InStream >> wT;
      InStream >> EVal;

      wTildeValues.push_back(wT);
      EValues.push_back(EVal);
//      cout << "wT " << wT << " EVal " << EVal << endl;

  }
//  cout << "read in background attractor file" << endl;

  // SETUP SPLINE //
  int NumberOfOpenMPThreads = omp_get_max_threads();
  EAcc = new gsl_interp_accel*[NumberOfOpenMPThreads];

  #pragma omp parallel for
  for(int i = 0; i < NumberOfOpenMPThreads; i++)
  {
      EAcc[i] = gsl_interp_accel_alloc();
  }

  double wTildeValuesSetUp[wTildeValues.size()];
  double EValuesSetUp[wTildeValues.size()];
  for (int i = 0; i < wTildeValues.size(); i++)
  {
    wTildeValuesSetUp[i] = wTildeValues[i];
    EValuesSetUp[i] = EValues[i];
  }

  EInt = gsl_spline_alloc(gsl_interp_cspline, wTildeValues.size());
  gsl_spline_init(EInt, wTildeValuesSetUp, EValuesSetUp,  wTildeValues.size());
//  cout << "set up interpolators for background attractor" << endl;

  // SET BOUNDARIES //
  wTMin = wTildeValues[0];
  wTMax = wTildeValues[wTildeValues.size() - 1];
  InStream.close();
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Read in data for greens functions and setup relevant variables
//##########################################################################################
void GreensFunctions::SetupGreensFunctions()
{
//  cout << "Setting up Greens Function" << endl;

  // ALLOCATE //
  wTValues = new double[greens_functions_chuncks];
  dXdTValues = new double[greens_functions_points];

  FssValues = new double[greens_functions_chuncks*greens_functions_points];
  GssValues = new double[greens_functions_chuncks*greens_functions_points];

  // SETUP GSL INTERPOLATION //
  int NumberOfOpenMPThreads = omp_get_max_threads();

  FsswTAcc = new gsl_interp_accel*[NumberOfOpenMPThreads];
  GsswTAcc = new gsl_interp_accel*[NumberOfOpenMPThreads];

  FssdXdTAcc = new gsl_interp_accel*[NumberOfOpenMPThreads];
  GssdXdTAcc = new gsl_interp_accel*[NumberOfOpenMPThreads];

  #pragma omp parallel for
  for (int i = 0; i < NumberOfOpenMPThreads; i++)
  {

      FsswTAcc[i] = gsl_interp_accel_alloc();
      GsswTAcc[i] = gsl_interp_accel_alloc();

      FssdXdTAcc[i] = gsl_interp_accel_alloc();
      GssdXdTAcc[i] = gsl_interp_accel_alloc();
  }

  FssInt = gsl_spline2d_alloc(gsl_interp2d_bilinear, greens_functions_chuncks, greens_functions_points);
  GssInt = gsl_spline2d_alloc(gsl_interp2d_bilinear, greens_functions_chuncks, greens_functions_points);

  //  SetValues Function
  double FssVal[greens_functions_chuncks*greens_functions_points];
  double GssVal[greens_functions_chuncks*greens_functions_points];

//  cout << "Reading in Greens Function" << endl;

  ifstream InStream;
  InStream.open(greens_functions_file);

  int wCounter = 0;
  int xCounter = 0;

  // READ INPUT FILE //
  while (InStream.good())
  {

      double wT;
      double dXdT;
      double Fss;
      double Gss;
      double Fsv;
      double Gsv;
      double Gstd;
      double Gstr;

      InStream >> wT;
      InStream >> dXdT;
      InStream >> Fss;
      InStream >> Gss;
      InStream >> Fsv;
      InStream >> Gsv;
      InStream >> Gstd;
      InStream >> Gstr;

//      cout << wT << " " << dXdT << " " << Fss << " " << Gss << endl;
      // WRITE EACH POSITION-STEP ONCE INTO dXdTValues //
      dXdTValues[xCounter] = dXdT;

      // WRITE EACH TIME-STEP ONCE INTO wTValues //
      if ((wCounter%greens_functions_points) == 0)
      {
          wTValues[wCounter/greens_functions_points] = wT;
      }

      // WRITE VALUES OF GREENS FUNCTIONS INTO ARRAYS //
      FssVal[wCounter] = Fss;
      GssVal[wCounter] = Gss;


      wCounter++;

      // RESET xCounter WHEN NEW BLOCK IS REACHED IN INPUT FILE //
      if (xCounter < greens_functions_points)
      {
          xCounter++;
      }

      if (xCounter == greens_functions_points)
      {
          xCounter = 0;
      }

  }

//  cout << "setting grid values Greens Function" << endl;

  // SET GRID VALUES FOR INTERPOLATION //
  int FIndex = 0;
  for (int wTIndex = 0; wTIndex < greens_functions_chuncks; wTIndex++)
  {

      for (int xIndex = 0; xIndex < greens_functions_points; xIndex++)
      {

          gsl_spline2d_set(FssInt, FssValues, wTIndex, xIndex, FssVal[FIndex]);
          gsl_spline2d_set(GssInt, GssValues, wTIndex, xIndex, GssVal[FIndex]);

          FIndex++;

      }

  }

  // SetupInterpolators function
  // SET BOUNDARIES //
  wTMin = wTValues[0];
  wTMax = wTValues[greens_functions_chuncks - 1];
  dXdTMin = dXdTValues[0];
  dXdTMax = dXdTValues[greens_functions_points - 1];

  // INITIALIZE INTERPOLATOR //
  gsl_spline2d_init(FssInt, wTValues, dXdTValues, FssValues, greens_functions_chuncks, greens_functions_points);
  gsl_spline2d_init(GssInt, wTValues, dXdTValues, GssValues, greens_functions_chuncks, greens_functions_points);

  // CLEAN-UP //
  delete[] wTValues;
  delete[] dXdTValues;

  delete[] FssValues;
  delete[] GssValues;

}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//
//##########################################################################################
double GreensFunctions::EVALUATE_GSL_INTERPOLATOR_2D(gsl_spline2d* Interpolator, double xValue, double yValue, gsl_interp_accel* xAccelerator, gsl_interp_accel* yAccelerator, double xMinValue, double xMaxValue, double yMinValue, double yMaxValue)
{
  if (yValue > yMaxValue)
  {
    return 0.0;
  }
  else if (xValue < xMinValue || xValue > xMaxValue)
  {
    if (xValue < xMinValue)
    {
      return gsl_spline2d_eval(Interpolator, xMinValue, yValue, xAccelerator, yAccelerator);
    }
    else
    {
      return gsl_spline2d_eval(Interpolator, xMaxValue, yValue, xAccelerator, yAccelerator);
    }
  }
  else
  {
    return gsl_spline2d_eval(Interpolator, xValue, yValue, xAccelerator, yAccelerator);
  }



/*  if (xValue < xMinValue || xValue > xMaxValue || yValue < yMinValue || yValue > yMaxValue)
  {
    if (yValue >= yMinValue && yValue <= yMaxValue)
    {
      if (xValue < xMinValue)
      {
//        cout << "test 2" << endl;

        return gsl_spline2d_eval(Interpolator,xMinValue,yValue,xAccelerator,yAccelerator);
      }
      else
      {
        return 0.0;
      }
    }
    else if (yValue > yMaxValue)
    {
      if (xValue < xMinValue)
      {
//        cout << "test 3" << endl;

        return gsl_spline2d_eval(Interpolator, xMinValue, yMaxValue, xAccelerator, yAccelerator);
      }
      else
      {
//        cout << "test 4" << endl;

        return gsl_spline2d_eval(Interpolator, xValue, yMaxValue, xAccelerator, yAccelerator);
      }
    }
    else
    {
//      cout << "test 5" << endl;

//      cerr << "#WARNING " << xValue << " " << xMinValue  << " " << xMaxValue << " " << yValue << " " << yMinValue  << " " << yMaxValue  << endl;
      return gsl_spline2d_eval(Interpolator, xValue, yMinValue, xAccelerator, yAccelerator);
    }
  }
  else
  {
//    cout << "test 6" << endl;

    return gsl_spline2d_eval(Interpolator, xValue, yValue, xAccelerator, yAccelerator);
  }
*/
}
//__________________________________________________________________________________________


//__________________________________________________________________________________________
//##########################################################################################
//
//##########################################################################################
double GreensFunctions::FssScalingCurve(double wT, double dXdT)
{
  int tID = omp_get_thread_num();

  if (dXdT < dXdTMin)
  {
    if (wT < wTMin)
    {
      return gsl_spline2d_eval(FssInt, wTMin, dXdTMin, FsswTAcc[tID], FssdXdTAcc[tID]);
    }
    else if (wT > wTMax)
    {
      return gsl_spline2d_eval(FssInt, wTMax, dXdTMin, FsswTAcc[tID], FssdXdTAcc[tID]);
    }
    else
    {
      return gsl_spline2d_eval(FssInt, wT, dXdTMin, FsswTAcc[tID], FssdXdTAcc[tID]);
    }
  }
  else
  {
    EVALUATE_GSL_INTERPOLATOR_2D(FssInt, wT, dXdT, FsswTAcc[tID], FssdXdTAcc[tID], wTMin, wTMax, dXdTMin, dXdTMax);
  }

/*
  int tID = omp_get_thread_num();
  return EVALUATE_GSL_INTERPOLATOR_2D(FssInt, wT, dXdT, FsswTAcc[tID], FssdXdTAcc[tID], wTMin, wTMax, dXdTMin, dXdTMax);
*/
} // FssScalingCurve

//__________________________________________________________________________________________
//##########################################################################################
//
//##########################################################################################
double GreensFunctions::Fss(double wT, double dXdT)
{
  return FssScalingCurve(wT, dXdT);
} // Fss


//__________________________________________________________________________________________
//##########################################################################################
//
//##########################################################################################
double GreensFunctions::GssScalingCurve(double wT, double dXdT)
{
  int tID = omp_get_thread_num();

  if (dXdT < dXdTMin)
  {
    if (wT < wTMin)
    {
      return gsl_spline2d_eval(GssInt, wTMin, dXdTMin, GsswTAcc[tID], GssdXdTAcc[tID]);
    }
    else if (wT > wTMax)
    {
      return gsl_spline2d_eval(GssInt, wTMax, dXdTMin, GsswTAcc[tID], GssdXdTAcc[tID]);
    }
    else
    {
      return gsl_spline2d_eval(GssInt, wT, dXdTMin, GsswTAcc[tID], GssdXdTAcc[tID]);
    }
  }
  else
  {
    EVALUATE_GSL_INTERPOLATOR_2D(GssInt, wT, dXdT, GsswTAcc[tID], GssdXdTAcc[tID], wTMin, wTMax, dXdTMin, dXdTMax);
  }



/*
  int tID = omp_get_thread_num();

  return EVALUATE_GSL_INTERPOLATOR_2D(GssInt, wT, dXdT, GsswTAcc[tID], GssdXdTAcc[tID], wTMin, wTMax, dXdTMin, dXdTMax);
*/
} // GssScalingCurve

//__________________________________________________________________________________________
//##########################################################################################
//
//##########################################################################################
double GreensFunctions::Gss(double wT, double dXdT)
{
//  cout << "Gss " << wT << " " << dXdT << endl;
  return GssScalingCurve(wT, dXdT);
} // Gss




//__________________________________________________________________________________________
//##########################################################################################
//
//##########################################################################################
double GreensFunctions::E(double wT)
{

//  cout << "In e wt " << wT << endl;

  if (wT < wTMin)
  {
    return 1.0/c_infinity*pow(wT, 4.0/9.0);
  }
  else if (wT > wTMax)
  {
    return 1.0 - 2.0/(3.0*M_PI*wT);
  }
  else
  {
    int tID = omp_get_thread_num();
    return gsl_spline_eval(EInt, wT, EAcc[tID]);
  }

} // E

//__________________________________________________________________________________________
//##########################################################################################
//
//##########################################################################################
void GreensFunctions::GetValues(double eTau0, double Tau, double etaOverS, double &e, double &wTilde)
{

//cout << "c_infinity " << c_infinity << endl;
  // DETERMINE (e(tau) tau^{4/3})_{infty} //
  double eTau43Infty = pow(4.0*M_PI*etaOverS, 4.0/9.0)*pow(M_PI*M_PI*nuEff/30.0, 1.0/9.0)*c_infinity*pow(eTau0, 8.0/9.0);
//  cout << "etaOverS " << etaOverS << " nuEff " << nuEff << " c_infinity " << c_infinity << " eTau0 " << eTau0 << endl;
//  cout << "eTau43Infty " << eTau43Infty << endl;

  //////////////////////////////////////////////////////////
  // DETERMINE TEMPERATURE SELF-CONSISTENTLY ACCORDING TO //
  // e(T)tau^{4/3} = E(wTilde) (e(tau) tau^{4/3})_{infty} //
  //          wTilde= (T tau)/(4pi eta/s)                 //
  //////////////////////////////////////////////////////////

  double TLow = 0.0;
  double THigh = pow(eTau43Infty/((M_PI*M_PI/30.0)*nuEff*pow(1.0, 4.0)*pow(Tau, 4.0/3.0)), 1.0/4.0);
//  cout << "THigh " << THigh << endl;

  double TMid = (THigh + TLow)/2.0;
  double wTildeMid = (TMid*Tau)/(4.0*M_PI*etaOverS);
//  cout << "wTildeMid " << wTildeMid << endl;


  while (THigh - TLow > 1E-6*TMid)
  {
//    cout << "THigh-TLow " << THigh-TLow << endl;
//    cout << "1E-6*TMid " << 1E-6*TMid << endl;
//    cout << "Output " << TLow << " " << THigh << " " << TMid << " " << wTildeMid << endl;
    if (E(wTildeMid)/pow(TMid, 4) > (M_PI*M_PI/30.0)*nuEff*pow(1.0, 4.0)*pow(Tau, 4.0/3.0)/eTau43Infty)
    {
      TLow = TMid;
    }
    else
    {
      THigh = TMid;
    }

    TMid = (THigh + TLow)/2.0;
    wTildeMid = (TMid*Tau)/(4.0*M_PI*etaOverS);
  }

  // CHECK THAT eEq(T) == E(wTilde) (e(tau) tau^{4/3})_{infty} IS SOLVED //
  //std::cerr << "wT=" << wTilde << " " << "eEq=" << (M_PI*M_PI/30.0)*nuEff*std::pow(TMid,4.0) << " " << "e=" << eTau43Infty*E((TMid*Tau)/(4.0*M_PI*etaOverS))/std::pow(Tau,4.0/3.0) << std::endl;

  // SET FINAL VALUE OF wTilde //
  wTilde = (TMid*Tau)/(4.0*M_PI*etaOverS);
//  cout << "wTilde " << wTilde << endl;

  // SET FINAL VALUES OF T,e IN GeV //
  e = (M_PI*M_PI/30.0)*nuEff*pow(TMid, 4.0);
//  cout << "e " << e << endl;
//  cout << "at end of get values" << endl;


} // GetValues
