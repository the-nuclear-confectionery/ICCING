#include "mask.h"

//__________________________________________________________________________________________
//##########################################################################################
//  Class constructor
//    Create empty Mask
//##########################################################################################
Mask::Mask(string type, double rad, double gridstep, double tau)
{
  mask_type = type;
  dist_rad = rad;
  grid_step = gridstep;
  tau_ = tau;

//  cout << mask_type << " " << dist_rad << " " << tau_ << " " << endl;
  if (mask_type == "Uniform")
  {
    UniformMask();
  }
  else if (mask_type == "Gaussian")
  {
    GaussianMask();
  }
  else if (mask_type == "Kernel")
  {
    KernelMask();
  }
  else if (mask_type == "Greens")
  {
    GreensFunctionMask();
  }

//  for (int i=0; i <dist_mask.size(); i++)
//  {
//    for (int j=0; j <dist_mask.size(); j++)
//    {cout << dist_mask[i][j] << " ";}

//    cout << endl;
//  }
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Class deconstructor
//##########################################################################################
Mask::~Mask()
{

}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Implicit Copy
//##########################################################################################
Mask::Mask(const Mask &original)
{
  CopyMask(original);
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Mask Copy Function
//##########################################################################################
void Mask::CopyMask(const Mask &e)
{
//  corr = e.corr;
  mask_type = e.mask_type;
  dist_rad = e.dist_rad;
  grid_step = e.grid_step;
  tau_ = e.tau_;

  dist_mask = e.dist_mask;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Overide operator=
//##########################################################################################
Mask& Mask::operator= (const Mask& original)
{
	CopyMask(original);
	return *this;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Uniform Mask
//##########################################################################################
void Mask::UniformMask()
{
  //******************************************************************************************
  //  Initialze Gluon Distribution for sampling
  // `  Using this simple mask, makes repeat calculations quick and easy since there is no
  //    need to calculate the circle multiple times.
  //******************************************************************************************
  //  Set size of gluon_dist grid used to sample energy from initial_energy
  dist_mask.resize(2*dist_rad + 1, vector<double>(2*dist_rad + 1, 0.));

  //  Initialize gluon distribution
  int ox = dist_rad;  //  x-value of gluon_dist center
  int oy = dist_rad;  //  y-value of gluon_dist center
  double point;
  double normalization = 0;

  //  Loop through only points in radius of gluon and set to 1
  for (int i = -dist_rad; i <= dist_rad; i++) //  This goes -radius to radius in x
  {
    // This calculates the hight of the gluon_dist at a given x-value
    int height = round(sqrt(pow(dist_rad, 2) - pow(i, 2)));
    for (int j = -height; j <= height; j++) //  This loops over the points in circle at given x
    {
      dist_mask[i + ox][j + oy] = 1;  //  Set points in circle to 1 for calculations
    }
  }

}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Uniform Mask
//##########################################################################################
void Mask::GaussianMask()
{
  //******************************************************************************************
  //  Initialze Gluon Distribution for sampling
  // `  Using this simple mask, makes repeat calculations quick and easy since there is no
  //    need to calculate the circle multiple times.
  //******************************************************************************************
  //  Set size of gluon_dist grid used to sample energy from initial_energy
  dist_mask.resize(2*dist_rad + 1, vector<double>(2*dist_rad + 1, 0.));

  //  Initialize gluon distribution
  int ox = dist_rad;  //  x-value of gluon_dist center
  int oy = dist_rad;  //  y-value of gluon_dist center
  double point;
  double normalization = 0;

  //  Loop through only points in radius of quark to get normalization factor
  for (int i = -dist_rad; i <= dist_rad; i++)  //  This goes -radius to radius in x
  {
    // This calculates the hight of the quark_dist at a given x-value
    int height = round(sqrt(pow(dist_rad, 2) - pow(i, 2)));
    for (int j = -height; j <= height; j++) //  This loops over the points in circle at given x
    {
      point = sqrt(pow((i), 2) + pow((j), 2));  //  Get distance of point from center of circle
      normalization += exp(-((pow(point, 2))/(2*pow(dist_rad, 2))));  //  Add value at poinnt to a normalization factor
    }
  }


  //  Loop through only points in radius of gluon and set to 1
  for (int i = -dist_rad; i <= dist_rad; i++) //  This goes -radius to radius in x
  {
    // This calculates the hight of the gluon_dist at a given x-value
    int height = round(sqrt(pow(dist_rad, 2) - pow(i, 2)));
    for (int j = -height; j <= height; j++) //  This loops over the points in circle at given x
    {
      point = sqrt(pow(i, 2) + pow(j, 2));  //  Get distance of point from center of circle
      //  Calculate value of gaussian at point in circle
      dist_mask[i + ox][j + oy] = 1/(normalization*pow(grid_step, 2)*tau_)*exp(-((pow(point, 2))/(2*pow(dist_rad, 2))));
    }
  }
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  KernelMask Mask
//##########################################################################################
void Mask::KernelMask()
{
  //******************************************************************************************
  //  Initialze Gluon Distribution for sampling
  // `  Using this simple mask, makes repeat calculations quick and easy since there is no
  //    need to calculate the circle multiple times.
  //******************************************************************************************
  //  Set size of gluon_dist grid used to sample energy from initial_energy
  dist_mask.resize(2*dist_rad + 1, vector<double>(2*dist_rad + 1, 0.));

  //  Initialize gluon distribution
  int ox = dist_rad;  //  x-value of gluon_dist center
  int oy = dist_rad;  //  y-value of gluon_dist center
  double point;
  double normalization = 0;
//  cout << "test 1" << endl;

    for (int i = -dist_rad; i <= dist_rad; i++)  //  This goes -radius to radius in x
    {
      for (int j = -dist_rad; j <= dist_rad; j++)  //  This goes -radius to radius in x
      {
        double q = sqrt(pow(i, 2) + pow(j, 2))/(dist_rad*0.5);
        double norm = 10.0/(7.0*M_PI*pow(dist_rad*0.5, 2)*pow(grid_step, 2)*tau_);
//        cout << "test 2 q: "  << q << endl;

        if (q >= 2.0)
        { dist_mask[i + ox][j + oy] = 0.0; }
        else if (q >= 1.0)
        {
          dist_mask[i + ox][j + oy] = 0.25*norm*pow(2.0 - q, 3);
        }
        else
        {
          dist_mask[i + ox][j + oy] = norm*(1 - 1.5*pow(q, 2) + 0.75*pow(q, 3));
        }
      }
    }
//    cout << "test 3" << endl;

}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  GreensFunctionMask Mask
//##########################################################################################
void Mask::GreensFunctionMask()
{
  //******************************************************************************************
  //  Initialze Gluon Distribution for sampling
  // `  Using this simple mask, makes repeat calculations quick and easy since there is no
  //    need to calculate the circle multiple times.
  //******************************************************************************************
  //  Set size of gluon_dist grid used to sample energy from initial_energy
  dist_mask.resize(2*dist_rad + 1, vector<double>(2*dist_rad + 1, 0.));

  //  Initialize gluon distribution
  int ox = dist_rad;  //  x-value of gluon_dist center
  int oy = dist_rad;  //  y-value of gluon_dist center
  double point;
  double normalization = 0;

  //  Loop through only points in radius of quark to get normalization factor
  for (int i = -dist_rad; i <= dist_rad; i++)  //  This goes -radius to radius in x
  {
    // This calculates the hight of the quark_dist at a given x-value
    int height = round(sqrt(pow(dist_rad, 2) - pow(i, 2)));
    for (int j = -height; j <= height; j++) //  This loops over the points in circle at given x
    {
      point = sqrt(pow(i, 2) + pow(j, 2));  //  Get distance of point from center of circle
      normalization += exp(-((pow(point, 2))/(2*pow(dist_rad, 2))));  //  Add value at poinnt to a normalization factor
    }
  }

  //  Loop through only points in radius of greens distribution and set to 1
  for (int i = -dist_rad; i <= dist_rad; i++) //  This goes -radius to radius in x
  {
    // This calculates the hight of the greens_dist at a given x-value
    int height = round(sqrt(pow(dist_rad, 2) - pow(i, 2)));
    for (int j = -height; j <= height; j++) //  This loops over the points in circle at given x
    {
      point = sqrt(pow(i, 2) + pow(j, 2));  //  Get distance of point from center of circle

      //  Calculate value of gaussian at point in circle
      dist_mask[i + ox][j + oy] = 1/(normalization*pow(grid_step, 2)*tau_)*exp(-pow(point, 2)/(2*pow(dist_rad, 2)));
    }

  }

}
//__________________________________________________________________________________________

