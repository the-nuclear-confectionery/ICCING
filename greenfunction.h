#ifndef GreensFunctions_H
#define GreensFunctions_H
//__________________________________________________________________________________________
//##########################################################################################
//  C++ Libraries
//##########################################################################################
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include <omp.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  ICCING Header Files
//##########################################################################################

//__________________________________________________________________________________________

using namespace std;

class GreensFunctions
{
private:
  //__________________________________________________________________________________________
  //##########################################################################################
  //  GreensFunctions Input Parameters
  //##########################################################################################
  string background_attractor_file;
  string greens_functions_file;

  int background_points;
  int greens_functions_points;
  int greens_functions_chuncks;
  double c_infinity;
  double eta_over_s;
  double tau_hydro;



  // EVALUATE BETWEEN dXdTMin,dXdTMax AND wTMin,wTMax //
  double dXdTMin, dXdTMax, wTMin, wTMax;

  // GSL INTERPOLATION OBJECTS //
  gsl_interp_accel **FsswTAcc, **GsswTAcc;
  gsl_interp_accel **FssdXdTAcc, **GssdXdTAcc;
  gsl_interp_accel **FsvwTAcc;
  gsl_interp_accel **FsvdXdTAcc;
  gsl_spline2d *FsvInt;
  gsl_spline2d *FssInt, *GssInt;

  double *wTValues,*dXdTValues;

  double *FssValues;
  double *GssValues;
  double *FsvValues;
  // GSL INTERPOLATION OBJECTS //
  gsl_interp_accel **EAcc;
  gsl_spline *EInt;

//    double wTMin; double wTMax;
   double CInfty;

   double nuEff;
  //__________________________________________________________________________________________

  //__________________________________________________________________________________________
  //##########################################################################################
  //  Internal Functions
  //##########################################################################################
  //  Copy function for GreensFunctions class, called by operator= and implicit copy functions
  void CopyGreensFunctions(const GreensFunctions &e);

  double EVALUATE_GSL_INTERPOLATOR_2D(gsl_spline2d* Interpolator, double xValue, double yValue, gsl_interp_accel* xAccelerator, gsl_interp_accel* yAccelerator, double xMinValue, double xMaxValue, double yMinValue, double yMaxValue);

  void SetupBackgroundAttractor();
  void SetupGreensFunctions();

  double FssScalingCurve(double wT, double dXdT);




  double GssScalingCurve(double wT, double dXdT);



  double E(double wT);
  //__________________________________________________________________________________________

public:

  //__________________________________________________________________________________________
  //##########################################################################################
  //  Basic Class Functions
  //##########################################################################################
  GreensFunctions(string backgroundAttractorFile, string greensFunctionsFile, int backgroundPoints, int greensFunctionsPoints, int greensFunctionsChuncks, double cInfinity, double etaOverS, double tauHydro);  // Class Constructor
  ~GreensFunctions(); //  Class Destructor

  GreensFunctions(const GreensFunctions &original); //  Implicit copy function, newIOObject(oldIOObject)
  GreensFunctions& operator=(const GreensFunctions& original);  //  Defines what happens when you use = operator on class
  //__________________________________________________________________________________________

  //__________________________________________________________________________________________
  //##########################################################################################
  //  GreensFunctions Specific Functions
  //##########################################################################################
  void GetValues(double eTau0, double Tau, double etaOverS, double &e, double &wTilde);

  double Fss(double wT, double dXdT);
  double Gss(double wT, double dXdT);
  //__________________________________________________________________________________________
};
#endif
