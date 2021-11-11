#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>

namespace GreensFunctions{


    // EVALUATION MACRO FOR GSL INTERPOLATION FUNCTIONS //
    #define EVALUATE_GSL_INTERPOLATOR_2D(Interpolator,xValue,yValue,xAccelerator,yAccelerator,xMinValue,xMaxValue,yMinValue,yMaxValue) \
    if((yValue)>(yMaxValue)){ \
        return 0.0; \
    } \
    else if((xValue)<(xMinValue) || (xValue)>(xMaxValue)){ \
        if((xValue)<(xMinValue)){ \
            return gsl_spline2d_eval(Interpolator,xMinValue,yValue,xAccelerator,yAccelerator); \
        } \
        else{ \
            return gsl_spline2d_eval(Interpolator,xMaxValue,yValue,xAccelerator,yAccelerator); \
        } \
    } \
    else{ \
	    return gsl_spline2d_eval(Interpolator,xValue,yValue,xAccelerator,yAccelerator); \
    } \
    
    ///////////////////////////////////////////////
    // EVALUTE GREENS FUNCTIONS BY INTERPOLATION //
    ///////////////////////////////////////////////

    // EVALUATE BETWEEN dXdTMin,dXdTMax AND wTMin,wTMax //
    double dXdTMin,dXdTMax,wTMin,wTMax;
    
    // GSL INTERPOLATION OBJECTS //
    gsl_interp_accel **FsswTAcc,**GsswTAcc,**FsvwTAcc,**GsvwTAcc,**GstdwTAcc,**GstrwTAcc;
    gsl_interp_accel **FssdXdTAcc,**GssdXdTAcc,**FsvdXdTAcc,**GsvdXdTAcc,**GstddXdTAcc,**GstrdXdTAcc;

    gsl_spline2d *FssInt,*GssInt,*FsvInt,*GsvInt,*GstdInt,*GstrInt;
    
    double FssScalingCurve(double wT,double dXdT){
        int tID=omp_get_thread_num();
        
        if((dXdT)<(dXdTMin)){
            if((wT)<(wTMin)){
                return gsl_spline2d_eval(FssInt,wTMin,dXdTMin,FsswTAcc[tID],FssdXdTAcc[tID]);
            }
            else if((wT)>(wTMax)){
                return gsl_spline2d_eval(FssInt,wTMax,dXdTMin,FsswTAcc[tID],FssdXdTAcc[tID]);
            }
            else{
                return gsl_spline2d_eval(FssInt,wT,dXdTMin,FsswTAcc[tID],FssdXdTAcc[tID]);
            }
        }
        else{
            EVALUATE_GSL_INTERPOLATOR_2D(FssInt,wT,dXdT,FsswTAcc[tID],FssdXdTAcc[tID],wTMin,wTMax,dXdTMin,dXdTMax);
        }

    } // FssScalingCurve
    
    double Fss(double wT,double dXdT){
        return FssScalingCurve(wT,dXdT);
    } // Fss
    

    double GssScalingCurve(double wT,double dXdT){
        int tID=omp_get_thread_num();

        if((dXdT)<(dXdTMin)){
            if((wT)<(wTMin)){
                return gsl_spline2d_eval(GssInt,wTMin,dXdTMin,GsswTAcc[tID],GssdXdTAcc[tID]);;
            }
            else if((wT)>(wTMax)){
                return gsl_spline2d_eval(GssInt,wTMax,dXdTMin,GsswTAcc[tID],GssdXdTAcc[tID]);
            }
            else{
                return gsl_spline2d_eval(GssInt,wT,dXdTMin,GsswTAcc[tID],GssdXdTAcc[tID]);
            }
        }
        else{
            EVALUATE_GSL_INTERPOLATOR_2D(GssInt,wT,dXdT,GsswTAcc[tID],GssdXdTAcc[tID],wTMin,wTMax,dXdTMin,dXdTMax);
        }
        

    } // GssScalingCurve

    double Gss(double wT,double dXdT){
        return GssScalingCurve(wT,dXdT);
    } // Gss


    double FsvScalingCurve(double wT,double dXdT){
        int tID=omp_get_thread_num();

        if((dXdT)<(dXdTMin)){
            return 0.0;
        }
        else{
            EVALUATE_GSL_INTERPOLATOR_2D(FsvInt,wT,dXdT,FsvwTAcc[tID],FsvdXdTAcc[tID],wTMin,wTMax,dXdTMin,dXdTMax);
        }

    } // FsvScalingCurve
    
    double Fsv(double wT,double dXdT){
        return FsvScalingCurve(wT,dXdT); 
    } // Fsv


    double GsvScalingCurve(double wT,double dXdT){
        int tID=omp_get_thread_num();

        if((dXdT)<(dXdTMin)){
            return 0.0;
        }
        else{
            EVALUATE_GSL_INTERPOLATOR_2D(GsvInt,wT,dXdT,GsvwTAcc[tID],GsvdXdTAcc[tID],wTMin,wTMax,dXdTMin,dXdTMax);
        }

    } // GsvScalingCurve

    double Gsv(double wT,double dXdT){
        return GsvScalingCurve(wT,dXdT);
    } // Gsv


    double GstdScalingCurve(double wT,double dXdT){
        int tID=omp_get_thread_num();

        if((dXdT)<(dXdTMin)){
            if((wT)<(wTMin)){
                return gsl_spline2d_eval(GstdInt,wTMin,dXdTMin,GstdwTAcc[tID],GstddXdTAcc[tID]);;
            }
            else if((wT)>(wTMax)){
                return gsl_spline2d_eval(GstdInt,wTMax,dXdTMin,GstdwTAcc[tID],GstddXdTAcc[tID]);
            }
            else{
                return gsl_spline2d_eval(GstdInt,wT,dXdTMin,GstdwTAcc[tID],GstddXdTAcc[tID]);
            }
        }
        else{
            EVALUATE_GSL_INTERPOLATOR_2D(GstdInt,wT,dXdT,GstdwTAcc[tID],GstddXdTAcc[tID],wTMin,wTMax,dXdTMin,dXdTMax);
        }

    } // GstdScalingCurve

    double Gstd(double wT,double dXdT){
        return GstdScalingCurve(wT,dXdT);
    } // Gstd


    double GstrScalingCurve(double wT,double dXdT){
        int tID=omp_get_thread_num();

        if((dXdT)<(dXdTMin)){
            return 0.0;
        }
        else{
            EVALUATE_GSL_INTERPOLATOR_2D(GstrInt,wT,dXdT,GstrwTAcc[tID],GstrdXdTAcc[tID],wTMin,wTMax,dXdTMin,dXdTMax);
        }

    } // GstrScalingCurve

    double Gstr(double wT,double dXdT){
        return GstrScalingCurve(wT,dXdT);
    } // Gstr
    
    
    /////////////////////////////////////////////////////////
    //    INTERPOLATE COORDINATE SPACE GRRENS FUNCTIONS    //
    // AS FUNCTIONS OF wTilde and (\Delta x)/(\Delta \tau) //
    /////////////////////////////////////////////////////////
    
    double *wTValues,*dXdTValues;
    
    double *FssValues;
    double *GssValues;
    double *FsvValues;
    double *GsvValues;
    double *GstdValues;
    double *GstrValues;
    
    void Setup(int NumberOfTimes,int NumberOfPoints){
        
        // ALLOCATE //
        wTValues=new double[NumberOfTimes];
        dXdTValues=new double[NumberOfPoints];
        
        FssValues=new double[NumberOfTimes*NumberOfPoints];
        GssValues=new double[NumberOfTimes*NumberOfPoints];
        FsvValues=new double[NumberOfTimes*NumberOfPoints];
        GsvValues=new double[NumberOfTimes*NumberOfPoints];
        GstdValues=new double[NumberOfTimes*NumberOfPoints];
        GstrValues=new double[NumberOfTimes*NumberOfPoints];
        
        // SETUP GSL INTERPOLATION //
        int NumberOfOpenMPThreads=omp_get_max_threads();
        
        FsswTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GsswTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        FsvwTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GsvwTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GstdwTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GstrwTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        
        FssdXdTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GssdXdTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        FsvdXdTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GsvdXdTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GstddXdTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GstrdXdTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        
        #pragma omp parallel for
        for(int i=0;i<NumberOfOpenMPThreads;i++){
            
            FsswTAcc[i] = gsl_interp_accel_alloc ();
            GsswTAcc[i] = gsl_interp_accel_alloc ();
            FsvwTAcc[i] = gsl_interp_accel_alloc ();
            GsvwTAcc[i] = gsl_interp_accel_alloc ();
            GstdwTAcc[i] = gsl_interp_accel_alloc ();
            GstrwTAcc[i] = gsl_interp_accel_alloc ();

            FssdXdTAcc[i] = gsl_interp_accel_alloc ();
            GssdXdTAcc[i] = gsl_interp_accel_alloc ();
            FsvdXdTAcc[i] = gsl_interp_accel_alloc ();
            GsvdXdTAcc[i] = gsl_interp_accel_alloc ();
            GstddXdTAcc[i] = gsl_interp_accel_alloc ();
            GstrdXdTAcc[i] = gsl_interp_accel_alloc ();

            }
        
        FssInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfTimes,NumberOfPoints);
        GssInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfTimes,NumberOfPoints);
        FsvInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfTimes,NumberOfPoints);
        GsvInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfTimes,NumberOfPoints);
        GstdInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfTimes,NumberOfPoints);
        GstrInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfTimes,NumberOfPoints);
        
    } // Setup
    
    // READ INPUT FILE //
    // INPUT FILE MUST HAVE FOLLOWING STRUCTURE: 1:wTilde 2:|x-x_0|/|tau-tau_0| 3:|tau-tau_0|^2*Fss 4:|tau-tau_0|^2*Gss 5:|tau-tau_0|^2*Fsv 6:|tau-tau_0|^2*Gsv 7:|tau-tau_0|^2*Gstd 8:|tau-tau_0|^2*Gstr //
    // AVOID EMPTY LINES IN INPUT FILE! //
    void SetValues(std::string fname,int NumberOfTimes,int NumberOfPoints){
        
        double FssVal[NumberOfTimes*NumberOfPoints];
        double GssVal[NumberOfTimes*NumberOfPoints];
        double FsvVal[NumberOfTimes*NumberOfPoints];
        double GsvVal[NumberOfTimes*NumberOfPoints];
        double GstdVal[NumberOfTimes*NumberOfPoints];
        double GstrVal[NumberOfTimes*NumberOfPoints];

        std::ifstream InStream;
        InStream.open(fname);
        
        int wCounter=0; int xCounter=0;
        
        // READ INPUT FILE //
        while(InStream.good()){
            
            double wT; double dXdT; double Fss; double Gss; double Fsv; double Gsv; double Gstd; double Gstr;
                
            InStream >> wT; InStream >> dXdT; InStream >> Fss; InStream >> Gss; InStream >> Fsv; InStream >> Gsv; InStream >> Gstd; InStream >> Gstr;


            // WRITE EACH POSITION-STEP ONCE INTO dXdTValues //
            dXdTValues[xCounter]=dXdT;

            // WRITE EACH TIME-STEP ONCE INTO wTValues //
            if((wCounter%NumberOfPoints)==0){
                wTValues[wCounter/NumberOfPoints]=wT;
            }

            // WRITE VALUES OF GREENS FUNCTIONS INTO ARRAYS //
            FssVal[wCounter]=Fss; GssVal[wCounter]=Gss; FsvVal[wCounter]=Fsv; GsvVal[wCounter]=Gsv; GstdVal[wCounter]=Gstd; GstrVal[wCounter]=Gstr;


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


        // SET GRID VALUES FOR INTERPOLATION //
        int FIndex=0;
        for(int wTIndex=0;wTIndex<NumberOfTimes;wTIndex++){

            for(int xIndex=0;xIndex<NumberOfPoints;xIndex++){

                gsl_spline2d_set(FssInt,FssValues,wTIndex,xIndex,FssVal[FIndex]);
                gsl_spline2d_set(GssInt,GssValues,wTIndex,xIndex,GssVal[FIndex]);
                gsl_spline2d_set(FsvInt,FsvValues,wTIndex,xIndex,FsvVal[FIndex]);
                gsl_spline2d_set(GsvInt,GsvValues,wTIndex,xIndex,GsvVal[FIndex]);
                gsl_spline2d_set(GstdInt,GstdValues,wTIndex,xIndex,GstdVal[FIndex]);
                gsl_spline2d_set(GstrInt,GstrValues,wTIndex,xIndex,GstrVal[FIndex]);
                
                FIndex++;

            }

        }
        
    } // SetValues
    
    /////////////////////////
    // SETUP INTERPOLATORS //
    /////////////////////////
    
    void SetupInterpolators(int NumberOfTimes,int NumberOfPoints){

        // SET BOUNDARIES //
        wTMin=wTValues[0]; wTMax=wTValues[NumberOfTimes-1];
        dXdTMin=dXdTValues[0]; dXdTMax=dXdTValues[NumberOfPoints-1];
        
        // INITIALIZE INTERPOLATOR //
        gsl_spline2d_init(FssInt,wTValues,dXdTValues,FssValues,NumberOfTimes,NumberOfPoints);
        gsl_spline2d_init(GssInt,wTValues,dXdTValues,GssValues,NumberOfTimes,NumberOfPoints);
        gsl_spline2d_init(FsvInt,wTValues,dXdTValues,FsvValues,NumberOfTimes,NumberOfPoints);
        gsl_spline2d_init(GsvInt,wTValues,dXdTValues,GsvValues,NumberOfTimes,NumberOfPoints);
        gsl_spline2d_init(GstdInt,wTValues,dXdTValues,GstdValues,NumberOfTimes,NumberOfPoints);
        gsl_spline2d_init(GstrInt,wTValues,dXdTValues,GstrValues,NumberOfTimes,NumberOfPoints);
        
        // CLEAN-UP //
        delete[] wTValues;
        delete[] dXdTValues;
        
        delete[] FssValues;
        delete[] GssValues;
        delete[] FsvValues;
        delete[] GsvValues;
        delete[] GstdValues;
        delete[] GstrValues;
        
    } // SetupInterpolators
    
    
    // CREATE OUTPUT //
    void Output(std::string fname,int NwT,int NdXdT){
        
        
        std::ofstream Outstream;
        Outstream.open(fname.c_str());
        Outstream << "# 1:wTilde=tau T(tau)/(eta/s) 2:|x-x_0|/|tau-tau0| 3:|tau-tau0|^2Fss 4:|tau-tau0|^2Gss 5:|tau-tau_0|^2Fsv 6:|tau-tau_0|^2Gsv 7:|tau-tau_0|^2Gstd 8:|tau-tau_0|^2Gstr" << std::endl;
        
        for(int tIndex=0;tIndex<NwT;tIndex++){
            
            double wTi=wTMin + (tIndex+0.5)*(wTMax-wTMin)/(NwT);
            
            for(int xIndex=0;xIndex<NdXdT;xIndex++){
                
                double dXdTi=dXdTMin + (xIndex+0.5)*(dXdTMax-dXdTMin)/(NdXdT);
                
                Outstream << wTi << " " << dXdTi << " " << Fss(wTi,dXdTi) << " " << Gss(wTi,dXdTi) << " " << Fsv(wTi,dXdTi) << " " << Gsv(wTi,dXdTi) << " " << Gstd(wTi,dXdTi) << " " << Gstr(wTi,dXdTi) << std::endl;
                
            }
            
            Outstream << std::endl;
            Outstream << std::endl;
            
        }
        
        Outstream.close();
    } // Output
    
    
} // GreensFunctions