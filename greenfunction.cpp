#include "greenfunction.h"

namespace GreensFunctions{

    ///////////////////////////////////////////////
    // EVALUTE GREENS FUNCTIONS BY INTERPOLATION //
    ///////////////////////////////////////////////
    // EVALUATE BETWEEN dXdTMin,dXdTMax AND wTMin,wTMax //
    double dXdTMin,dXdTMax,wTMin,wTMax;

    // GSL INTERPOLATION OBJECTS //
    gsl_interp_accel **FsswTAcc,**GsswTAcc;
    gsl_interp_accel **FssdXdTAcc,**GssdXdTAcc;

    gsl_spline2d *FssInt,*GssInt;


    double FssScalingCurve(double wT,double dXdT){
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(FssInt,wT,dXdT,FsswTAcc[tID],FssdXdTAcc[tID],wTMin,wTMax,dXdTMin,dXdTMax);
    } // FssScalingCurve

    double Fss(double wT,double dXdT){
        return FssScalingCurve(wT,dXdT);
    } // Fss


    double GssScalingCurve(double wT,double dXdT){
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(GssInt,wT,dXdT,GsswTAcc[tID],GssdXdTAcc[tID],wTMin,wTMax,dXdTMin,dXdTMax);
    } // GssScalingCurve

    double Gss(double wT,double dXdT){
        return GssScalingCurve(wT,dXdT);
    } // Gss


    /////////////////////////////////////////////////////////
    //    INTERPOLATE COORDINATE SPACE GRRENS FUNCTIONS    //
    // AS FUNCTIONS OF wTilde and (\Delta x)/(\Delta \tau) //
    /////////////////////////////////////////////////////////
//    double *wTValues,*dXdTValues;

//    double *FssValues;
//    double *GssValues;


/*    void Setup(int NumberOfTimes,int NumberOfPoints){

        // ALLOCATE //
        wTValues=new double[NumberOfTimes];
        dXdTValues=new double[NumberOfPoints];

        FssValues=new double[NumberOfTimes*NumberOfPoints];
        GssValues=new double[NumberOfTimes*NumberOfPoints];

        // SETUP GSL INTERPOLATION //
        int NumberOfOpenMPThreads=omp_get_max_threads();

        FsswTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GsswTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];

        FssdXdTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GssdXdTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];

        #pragma omp parallel for
        for(int i=0;i<NumberOfOpenMPThreads;i++){

            FsswTAcc[i] = gsl_interp_accel_alloc ();
            GsswTAcc[i] = gsl_interp_accel_alloc ();

            FssdXdTAcc[i] = gsl_interp_accel_alloc ();
            GssdXdTAcc[i] = gsl_interp_accel_alloc ();
            }

        FssInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfTimes,NumberOfPoints);
        GssInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfTimes,NumberOfPoints);

    } // Setup
*/
    // READ INPUT FILE //
    // INPUT FILE MUST HAVE FOLLOWING STRUCTURE: 1:wTilde 2:|x-x_0|/|tau-tau_0| 3:|tau-tau_0|^2*Fss 4:|tau-tau_0|^2*Gss //
    // AVOID EMPTY LINES IN INPUT FILE! //
/*    void SetValues(std::string fname,int NumberOfTimes,int NumberOfPoints){

        double FssVal[NumberOfTimes*NumberOfPoints];
        double GssVal[NumberOfTimes*NumberOfPoints];

        std::ifstream InStream;
        InStream.open(fname);

        int wCounter=0; int xCounter=0;

        // READ INPUT FILE //
        while(InStream.good()){

            double wT; double dXdT; double Fss; double Gss;

            InStream >> wT; InStream >> dXdT; InStream >> Fss; InStream >> Gss;


            // WRITE EACH POSITION-STEP ONCE INTO dXdTValues //
            dXdTValues[xCounter]=dXdT;

            // WRITE EACH TIME-STEP ONCE INTO wTValues //
            if((wCounter%NumberOfPoints)==0){
                wTValues[wCounter/NumberOfPoints]=wT;
            }

            // WRITE VALUES OF GREENS FUNCTIONS INTO ARRAYS //
            FssVal[wCounter]=Fss; GssVal[wCounter]=Gss;


            wCounter++;

            // RESET xCounter WHEN NEW BLOCK IS REACHED IN INPUT FILE //
            if(xCounter<NumberOfPoints){
                xCounter++;
            }

            if(xCounter==NumberOfPoints){
                xCounter=0;
            }

        }

        // CHECK WHETHER TIMES wTilde AND POSITIONS dXdT ARE READ CORRECTLY //
        /* int k=0;
        while (k < NumberOfTimes) {
           std::cout << wTValues[k] << " ";
           printf("\n");
           k++;
        }

        int l=0;
        while (l < NumberOfPoints) {
           std::cout << dXdTValues[l] << " ";
           printf("\n");
           l++;
        } */

/*
        // SET GRID VALUES FOR INTERPOLATION //
        int FIndex=0;
        for(int wTIndex=0;wTIndex<NumberOfTimes;wTIndex++){

            for(int xIndex=0;xIndex<NumberOfPoints;xIndex++){

                gsl_spline2d_set(FssInt,FssValues,wTIndex,xIndex,FssVal[FIndex]);
                gsl_spline2d_set(GssInt,GssValues,wTIndex,xIndex,GssVal[FIndex]);

                FIndex++;

            }

        }

    } // SetValues
*/
    /////////////////////////
    // SETUP INTERPOLATORS //
    /////////////////////////

/*    void SetupInterpolators(int NumberOfTimes,int NumberOfPoints){

        // SET BOUNDARIES //
        wTMin=wTValues[0]; wTMax=wTValues[NumberOfTimes-1];
        dXdTMin=dXdTValues[0]; dXdTMax=dXdTValues[NumberOfPoints-1];

        // INITIALIZE INTERPOLATOR //
        gsl_spline2d_init(FssInt,wTValues,dXdTValues,FssValues,NumberOfTimes,NumberOfPoints);
        gsl_spline2d_init(GssInt,wTValues,dXdTValues,GssValues,NumberOfTimes,NumberOfPoints);

        // CLEAN-UP //
        delete[] wTValues;
        delete[] dXdTValues;

        delete[] FssValues;
        delete[] GssValues;

    } // SetupInterpolators
*/

    // CREATE OUTPUT //
/*    void Output(std::string fname,int NwT,int NdXdT){


        std::ofstream Outstream;
        Outstream.open(fname.c_str());
        Outstream << "# 1:wTilde=tau T(tau)/(eta/s) 2:|x-x_0|/|tau-tau0| 3:|tau-tau0|^2Fss 4:|tau-tau0|^2Gss" << std::endl;

        for(int tIndex=0;tIndex<NwT;tIndex++){

            double wTi=wTMin + (tIndex+0.5)*(wTMax-wTMin)/(NwT);

            for(int xIndex=0;xIndex<NdXdT;xIndex++){

                double dXdTi=dXdTMin + (xIndex+0.5)*(dXdTMax-dXdTMin)/(NdXdT);

                Outstream << wTi << " " << dXdTi << " " << Fss(wTi,dXdTi) << " " << Gss(wTi,dXdTi) << std::endl;

            }

            Outstream << std::endl;
            Outstream << std::endl;

        }

        Outstream.close();
    } // Output
*/

// GSL INTERPOLATION OBJECTS //
    gsl_interp_accel **EAcc;
    gsl_spline *EInt;

//    double wTMin; double wTMax;
     double CInfty;

    // ENERGY ATTRACTOR CURVE //
    double E(double wT){

       double Nc=3.0;
      double Nf=3.0;
      double nuG=2.0*(Nc*Nc-1.0);
      double nuQ=2.0*Nc*Nf;
      double nuEff=nuG+7.0/4.0*nuQ;

        if(wT<wTMin){
            return 1.0/CInfty*std::pow(wT,4.0/9.0);
        }
        else if(wT>wTMax){
            return 1.0-2.0/(3.0*M_PI*wT);
        }
        else{
            int tID=omp_get_thread_num();
            return gsl_spline_eval(EInt,wT,EAcc[tID]);
        }

    } // E

    void GetValues(double eTau0,double Tau,double etaOverS,double &e,double &wTilde){

       double Nc=3.0;
      double Nf=3.0;
      double nuG=2.0*(Nc*Nc-1.0);
      double nuQ=2.0*Nc*Nf;
      double nuEff=nuG+7.0/4.0*nuQ;
        // DETERMINE (e(tau) tau^{4/3})_{infty} //
        double eTau43Infty=std::pow(4.0*M_PI*etaOverS,4.0/9.0)*std::pow(M_PI*M_PI*nuEff/30.0,1.0/9.0)*CInfty*std::pow(eTau0,8.0/9.0);

        //////////////////////////////////////////////////////////
        // DETERMINE TEMPERATURE SELF-CONSISTENTLY ACCORDING TO //
        // e(T)tau^{4/3} = E(wTilde) (e(tau) tau^{4/3})_{infty} //
        //          wTilde= (T tau)/(4pi eta/s)                 //
        //////////////////////////////////////////////////////////

        double TLow=0.0; double THigh=std::pow(eTau43Infty/((M_PI*M_PI/30.0)*nuEff*std::pow(1.0,4.0)*std::pow(Tau,4.0/3.0)),1.0/4.0);

        double TMid=(THigh+TLow)/2.0;
        double wTildeMid=(TMid*Tau)/(4.0*M_PI*etaOverS);


        while(THigh-TLow>1E-6*TMid){

            if(E(wTildeMid)/std::pow(TMid,4)>(M_PI*M_PI/30.0)*nuEff*std::pow(1.0,4.0)*std::pow(Tau,4.0/3.0)/eTau43Infty){
                TLow=TMid;
            }
            else{
                THigh=TMid;
            }

            TMid=(THigh+TLow)/2.0;
            wTildeMid=(TMid*Tau)/(4.0*M_PI*etaOverS);

        }

        // CHECK THAT eEq(T) == E(wTilde) (e(tau) tau^{4/3})_{infty} IS SOLVED //
        //std::cerr << "wT=" << wTilde << " " << "eEq=" << (M_PI*M_PI/30.0)*nuEff*std::pow(TMid,4.0) << " " << "e=" << eTau43Infty*E((TMid*Tau)/(4.0*M_PI*etaOverS))/std::pow(Tau,4.0/3.0) << std::endl;

        // SET FINAL VALUE OF wTilde //
        wTilde=(TMid*Tau)/(4.0*M_PI*etaOverS);

        // SET FINAL VALUES OF T,e IN GeV //
        e=(M_PI*M_PI/30.0)*nuEff*std::pow(TMid,4.0);



    } // GetValues

} // GreensFunctions
