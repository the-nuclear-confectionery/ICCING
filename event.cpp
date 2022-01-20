#include "event.h"

//__________________________________________________________________________________________
//##########################################################################################
//  Class constructor
//    Create empty Event
//##########################################################################################
Event::Event()
{

}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Class deconstructor
//##########################################################################################
Event::~Event()
{
  CleanEvent();
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Implicit Copy
//##########################################################################################
Event::Event(const Event &original)
{
	CopyEvent(original);
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Event Copy Function
//##########################################################################################
void Event::CopyEvent(const Event &e)
{
  kappa_ = e.kappa_;
  gluon_rad = e.gluon_rad;
  quark_rad = e.quark_rad;
  lambda_ = e.lambda_;
  grid_max = e.grid_max;
  grid_step = e.grid_step;
  tau_0 = e.tau_0;
  e_thresh = e.e_thresh;
  grid_points = e.grid_points;
  test_ = e.test_;
  get_grid_point = e.get_grid_point;
  output_dir = e.output_dir;

  initial_energy = e.initial_energy;
  t_a = e.t_a;
  t_b = e.t_b;
  density = e.density;
  gluon_dist = e.gluon_dist;
  quark_dist = e.quark_dist;
  valued_points = e.valued_points;
  eccentricities = e.eccentricities;
  ecc = e.ecc;
  number_gluon = e.number_gluon;
  number_up = e.number_up;
  number_down = e.number_down;
  number_strange = e.number_strange;
  number_charm = e.number_charm;
  x_center = e.x_center;
  y_center = e.y_center;

  total_initial_energy = e.total_initial_energy;
  total_initial_entropy = e.total_initial_entropy;
  total_energy = e.total_energy;
  total_entropy = e.total_entropy;
  out_sample = e.out_sample;
  seed = e.seed;

  up_chop = e.up_chop;
  down_chop = e.down_chop;
  strange_chop = e.strange_chop;
  charm_chop = e.charm_chop;

  evolution = e.evolution;
  w_tilde = e.w_tilde;
  tau_hydro = e.tau_hydro;
  eta_over_s = e.eta_over_s;
  greens_rad = e.greens_rad;
  greens_dist = e.greens_dist;
  initial_energy_backup = e.initial_energy_backup;
  final_energy_backup = e.final_energy_backup;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Overide operator=
//##########################################################################################
Event& Event::operator= (const Event& original)
{
	CopyEvent(original);
	return *this;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Select energy of gluon
//##########################################################################################
Sample Event::GetGlue()
{
  double q_s, e_tot;
  Sample samp;
  int total_points = 0;

  //  Get bounds of gluon using center point as defined by SampleEnergy
  //    Makes sure calculations are only done on points in initial_energy
  vector<int> gluon_bounds = GetIntegrationBounds(gluon_dist.size(), gluon_rad, x_center, y_center);

  //  Loop over gluon_dist using gluon_bounds
  for (int i = gluon_bounds[0]; i < gluon_bounds[2]; i++)
  {
    for (int j = gluon_bounds[1]; j < gluon_bounds[3]; j++)
    {
      //  Reminder: gluon_dist is a circular mask of 1's for ease of calculation

      //  Sum up q_s of gluon region
      // think about shifting calculation of qs to here
      q_s += kappa_*sqrt(t_b[x_center - gluon_rad + i][y_center - gluon_rad + j])*gluon_dist[i][j];

      //  Sum up total energy from gluon region
      e_tot += initial_energy[x_center - gluon_rad + i][y_center - gluon_rad + j]*gluon_dist[i][j];

      if(gluon_dist[i][j] == 1)
      {
        total_points++; //  Calculate total for normalization of q_s
      }
    }
  }

  //  Set normalized q_s and e_tot for output
  samp.q_s = q_s/total_points;
  samp.e_tot = pow(grid_step,2)*tau_0*e_tot;

  return samp;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Select energy of gluon
//##########################################################################################
double Event::GetOriginalEnergy()
{
  double e_tot;
  int total_points = 0;

  //  Get bounds of gluon using center point as defined by SampleEnergy
  //    Makes sure calculations are only done on points in initial_energy
  vector<int> gluon_bounds = GetIntegrationBounds(gluon_dist.size(), gluon_rad, x_center, y_center);

  //  Loop over gluon_dist using gluon_bounds
  for (int i = gluon_bounds[0]; i < gluon_bounds[2]; i++)
  {
    for (int j = gluon_bounds[1]; j < gluon_bounds[3]; j++)
    {
      //  Reminder: gluon_dist is a circular mask of 1's for ease of calculation

      //  Sum up total energy from gluon region
      e_tot += initial_energy_backup[x_center - gluon_rad + i][y_center - gluon_rad + j]*gluon_dist[i][j];

      if(gluon_dist[i][j] == 1)
      {
        total_points++; //  Calculate total for normalization of q_s
      }
    }
  }

  //  Set normalized q_s and e_tot for output
  e_tot = pow(grid_step,2)*tau_0*e_tot/total_points;

  return e_tot;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Sample Initial Energy for ICCING algorithm
//##########################################################################################
Sample Event::SampleEnergy()
{
  int point;
  bool got_point = false;

  //  Initialize distribution for selecting points from initial_energy
  get_grid_point = uniform_int_distribution<int>(0, valued_points.size()-1);

  //  Loop until a point gives e_tot > e_thresh
  while (!got_point)
  {
    //  Get random valued point for gluon center
    point = get_grid_point(get_random_number);

    //  Store x and y of gluon center in class variables for use in other functions
    x_center = valued_points[point][0];
    y_center = valued_points[point][1];

    //  Get the total energy and q_s using center picked above
    out_sample = GetGlue();

    //  Test if e_tot of gluon is less than e_thresh
    //    If true, then just copy energy over to density[0]
    //    Else, set got_point = true and pass out_sample to Splitter
    if (out_sample.e_tot < e_thresh)
    {
      //  Copy over energy to density[0], not enough to run though algorithm
      UpdateEnergy(1.);

      //  If total energy is less than e_thresh then no more quarks can be made
      //  so copy all energy left to density[0]
      if (total_initial_energy < e_thresh)
      {
        for (int i = 0; i < initial_energy.size(); i++)
        {
          for (int j = 0; j < initial_energy.size(); j++)
          {

            if (test_ == "GreensFunction")
            {
              // skip this since all of the energy is already copied over
              initial_energy[i][j] = 0;
              continue;
            }
            //  Copy all energy left to density[0] and set all initial_energy = 0
            density[0][i][j] += initial_energy[i][j];
            initial_energy[i][j] = 0;
          }
        }

        total_initial_energy = 0; //  Set total_energy = 0
        out_sample.q_s = -100;  //  Set flag to end event (q_s = -100, unphysical value)
        got_point = true; //  Set got_point to exit loop
      }

    }
    else
    {
      got_point = true; //  Set got_point to exit loop
    }
  }

  return out_sample;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Propogates Results of Splitter
//##########################################################################################
bool Event::UpdateDensity(Quarks quark_density)
{
  //******************************************************************************************
  //  If quark charge is an gluon, copy that energy over to output
  //******************************************************************************************
  if (quark_density.GetCharge()[0] == 0.)
  {
    number_gluon++;
    UpdateEnergy(quark_density.GetEnergyFraction());
  }
  //******************************************************************************************
  //  If quark charge is not a gluon, move quark to output
  //******************************************************************************************
  else {
    int temp_x;
    int temp_y;
    double energy;
    double original_energy;

    //******************************************************************************************
    //  Calculate centers of Quark and Anti-Quark
    //******************************************************************************************
    int quark_x = x_center + round((1 - quark_density.GetAlpha())*quark_density.GetPosition()[0]);
    int quark_y =  y_center + round((1 - quark_density.GetAlpha())*quark_density.GetPosition()[1]);
    int antiquark_x = x_center - round(quark_density.GetAlpha()*quark_density.GetPosition()[0]);
    int antiquark_y = y_center - round(quark_density.GetAlpha()*quark_density.GetPosition()[1]);

    vector<int> quark_bounds;
    vector<int> antiquark_bounds;
    //******************************************************************************************
    //  Test if Quark is in bounds
    //******************************************************************************************

    if (test_ == "GreensFunction")
    {
//      cout << "test 1" << endl;
      quark_bounds = GetIntegrationBounds(greens_dist.size(), greens_rad, quark_x, quark_y);
      if (abs(quark_bounds[0] - quark_bounds[2]) < greens_dist.size() || abs(quark_bounds[1] - quark_bounds[3]) < greens_dist.size())
      { return false; }
    }
    else
    {
      quark_bounds = GetIntegrationBounds(quark_dist.size(), quark_rad, quark_x, quark_y);
      if (abs(quark_bounds[0] - quark_bounds[2]) < quark_dist.size() || abs(quark_bounds[1] - quark_bounds[3]) < quark_dist.size())
      { return false; }
    }

    //******************************************************************************************
    //  Test if Anti-Quark is in bounds
    //******************************************************************************************

    if (test_ == "GreensFunction")
    {
//      cout << "test 2" << endl;
      antiquark_bounds = GetIntegrationBounds(greens_dist.size(), greens_rad, antiquark_x, antiquark_y);
      if (abs(antiquark_bounds[0] - antiquark_bounds[2]) < greens_dist.size() || abs(antiquark_bounds[1] - antiquark_bounds[3]) < greens_dist.size())
      { return false; }
    }
    else
    {
      antiquark_bounds = GetIntegrationBounds(quark_dist.size() , quark_rad, antiquark_x, antiquark_y);
      if (antiquark_bounds[0] - antiquark_bounds[2] < quark_dist.size() || antiquark_bounds[1] - antiquark_bounds[3] < quark_dist.size())
      { return false; }
    }


//    original_energy = GetOriginalEnergy();

//    cout << total_initial_energy << " original_energy " << original_energy << " up_chop " << up_chop << endl;
    if (quark_density.GetCharge()[0] == 0.0023)
    {
  //    if (original_energy < up_chop) {  return true;  }
      number_up++;
    }
    else if (quark_density.GetCharge()[0] == 0.0048)
    {
//      if (original_energy < down_chop) {  return true;  }
      number_down++;
    }
    else if (quark_density.GetCharge()[0] == 0.095)
    {
//      if (original_energy < strange_chop) {  return true;  }
      number_strange++;
    }
    else if (quark_density.GetCharge()[0] == 1.29)
    {
//      if (original_energy < charm_chop) {  return true;  }
      number_charm++;
    }

    if (test_ == "hotspots")
    {
      ofstream output;
      output.open(output_dir + "hotspottracking.dat", ios::app);  //  Append event to end of file
      output << quark_density.GetCharge()[0] << " ";

      int total_points_gluon = 0;
      double gluon_energy = 0;
    }


    //******************************************************************************************
    //  Update Total energies and initial_energy
    //******************************************************************************************
    vector<int> gluon_bounds = GetIntegrationBounds(gluon_dist.size(), gluon_rad, x_center, y_center);

    for (int i = gluon_bounds[0]; i < gluon_bounds[2]; i++)
    {
      for (int j = gluon_bounds[1]; j < gluon_bounds[3]; j++)
      {
        temp_x = x_center - gluon_rad + i;
        temp_y = y_center - gluon_rad + j;
        energy = gluon_dist[i][j]*quark_density.GetEnergyFraction()*initial_energy[temp_x][temp_y];

        total_initial_energy -= energy;
        total_energy += energy;
        initial_energy[temp_x][temp_y] -= energy;

        if (test_ == "hotspots")
        {
          gluon_energy += energy;
          if(gluon_dist[i][j] == 1) { total_points_gluon++; }
        }
        // This removes the gluon from the final state which was chosen to split and is now being redistributed
        if (test_ == "GreensFunction")
        {
          double evolved_energy = 0, wtilde = 0;
          evolution.GetValues(tau_0*energy, tau_hydro, eta_over_s, evolved_energy, wtilde);

          density[0][temp_x][temp_y] -= evolved_energy;
        }
      }
    }

    if (test_ == "hotspots")
    {
      output << GetOriginalEnergy() << " " << GetOriginalEnergy()/total_points_gluon << " ";
      output << out_sample.e_tot << " " << out_sample.e_tot/total_points_gluon << " ";
      output << gluon_energy << " " << gluon_energy/total_points_gluon << endl;
      output.close();
    }
    //******************************************************************************************
    //  Update Output Densities
    //******************************************************************************************
    for (int i = quark_bounds[0]; i < quark_bounds[2]; i++)
    {
      for (int j = quark_bounds[1]; j < quark_bounds[3]; j++)
      {

        if (test_ == "GreensFunction")
        {
//          cout << "test 3" << endl;

          //  Deposit Quark Energy and Charges
          temp_x = quark_x - greens_rad + i;
          temp_y = quark_y - greens_rad + j;
          double quark_distance = sqrt(pow((quark_x - temp_x)*grid_step, 2) + pow((quark_y - temp_y)*grid_step, 2));
/*          cout << "test 4 " << greens_dist.size() << endl;
          cout << quark_bounds[0] << " " << quark_bounds[1] << " " << quark_bounds[2] << " " << quark_bounds[3] << endl;
          cout << i << " " << j << endl;
          cout << greens_dist[i][j] << endl;
          cout << final_energy_backup[quark_x][quark_y] << endl;
          cout << initial_energy_backup[quark_x][quark_y] << endl;
          cout << w_tilde[quark_x][quark_y] << endl;
          cout << quark_distance << endl;
          cout << (tau_hydro - tau_0) << endl;
          cout << evolution.Gss(w_tilde[quark_x][quark_y], quark_distance/(tau_hydro - tau_0)) << endl;
  */        //  Energy = alpha*(E_glueon/E_tot)*E_tot*quark_dist
          density[0][temp_x][temp_y] += quark_density.GetAlpha()*(quark_density.GetEnergyFraction()*out_sample.e_tot)*greens_dist[i][j]
                                        *(final_energy_backup[quark_x][quark_y]/initial_energy_backup[quark_x][quark_y])
                                        *evolution.Gss(w_tilde[quark_x][quark_y], quark_distance/(tau_hydro - tau_0));
//          cout << "test 5" << endl;
//  Baryon = baron_number*quark_dist
          density[1][temp_x][temp_y] += quark_density.GetCharge()[1]*greens_dist[i][j]
                                        *(tau_0/tau_hydro)*evolution.Fss(w_tilde[quark_x][quark_y], quark_distance/(tau_hydro - tau_0));
          //  Strangeness = strangeness*quark_dist
//          cout << "test 6" << endl;
          density[2][temp_x][temp_y] += quark_density.GetCharge()[2]*greens_dist[i][j]
                                        *(tau_0/tau_hydro)*evolution.Fss(w_tilde[quark_x][quark_y], quark_distance/(tau_hydro - tau_0));;
          //  EM_charge = em_charge*quark_dist
//          cout << "test 7" << endl;
          density[3][temp_x][temp_y] += quark_density.GetCharge()[3]*greens_dist[i][j]
                                        *(tau_0/tau_hydro)*evolution.Fss(w_tilde[quark_x][quark_y], quark_distance/(tau_hydro - tau_0));

          //  Deposit Anti-Quark Energy and Charges
          temp_x = antiquark_x - greens_rad + i;
          temp_y = antiquark_y - greens_rad + j;
          double antiquark_distance = sqrt(pow((antiquark_x - temp_x)*grid_step, 2) + pow((antiquark_y - temp_y)*grid_step, 2));

          //  Energy = alpha*(E_glueon/E_tot)*E_tot*quark_dist
          density[0][temp_x][temp_y] += (1 - quark_density.GetAlpha())*(quark_density.GetEnergyFraction()*out_sample.e_tot)*greens_dist[i][j]
                                        *(final_energy_backup[antiquark_x][antiquark_y]/initial_energy_backup[antiquark_x][antiquark_y])
                                        *evolution.Gss(w_tilde[antiquark_x][antiquark_y], antiquark_distance/(tau_hydro - tau_0));
          //  Baryon = baron_number*quark_dist
          density[1][temp_x][temp_y] -= quark_density.GetCharge()[1]*greens_dist[i][j]
                                        *(tau_0/tau_hydro)*evolution.Fss(w_tilde[antiquark_x][antiquark_y], antiquark_distance/(tau_hydro - tau_0));
          //  Strangeness = strangeness*quark_dist
          density[2][temp_x][temp_y] -= quark_density.GetCharge()[2]*greens_dist[i][j]
                                        *(tau_0/tau_hydro)*evolution.Fss(w_tilde[antiquark_x][antiquark_y], antiquark_distance/(tau_hydro - tau_0));
          //  EM_charge = em_charge*quark_dist
          density[3][temp_x][temp_y] -= quark_density.GetCharge()[3]*greens_dist[i][j]
                                        *(tau_0/tau_hydro)*evolution.Fss(w_tilde[antiquark_x][antiquark_y], antiquark_distance/(tau_hydro - tau_0));
        }
        else
        {
          //  Deposit Quark Energy and Charges
          temp_x = quark_x - quark_rad + i;
          temp_y = quark_y - quark_rad + j;

          //  Energy = alpha*(E_glueon/E_tot)*E_tot*quark_dist
          density[0][temp_x][temp_y] += quark_density.GetAlpha()*(quark_density.GetEnergyFraction()*out_sample.e_tot)*quark_dist[i][j];
          //  Baryon = baron_number*quark_dist
          density[1][temp_x][temp_y] += quark_density.GetCharge()[1]*quark_dist[i][j];
          //  Strangeness = strangeness*quark_dist
          density[2][temp_x][temp_y] += quark_density.GetCharge()[2]*quark_dist[i][j];
          //  EM_charge = em_charge*quark_dist
          density[3][temp_x][temp_y] += quark_density.GetCharge()[3]*quark_dist[i][j];

          //  Deposit Anti-Quark Energy and Charges
          temp_x = antiquark_x - quark_rad + i;
          temp_y = antiquark_y - quark_rad + j;
          //  Energy = alpha*(E_glueon/E_tot)*E_tot*quark_dist
          density[0][temp_x][temp_y] += (1 - quark_density.GetAlpha())*(quark_density.GetEnergyFraction()*out_sample.e_tot)*quark_dist[i][j];
          //  Baryon = baron_number*quark_dist
          density[1][temp_x][temp_y] -= quark_density.GetCharge()[1]*quark_dist[i][j];
          //  Strangeness = strangeness*quark_dist
          density[2][temp_x][temp_y] -= quark_density.GetCharge()[2]*quark_dist[i][j];
          //  EM_charge = em_charge*quark_dist
          density[3][temp_x][temp_y] -= quark_density.GetCharge()[3]*quark_dist[i][j];
        }

      }
    }
  }

  return true;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Calculate Eccentricities from density grids
//##########################################################################################
void Event::CalculateEccentricities()
{
  eccentricities = ecc.CalculateEccentricities(grid_max, grid_step, density);
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Subtracts energy from initial_energy and adds it to density[0]
//##########################################################################################
void Event::UpdateEnergy(double ratio)
{
  //  Get bounds of gluon using center point as defined by SampleEnergy
  //    Makes sure calculations are only done on points in initial_energy
  vector<int> gluon_bounds = GetIntegrationBounds(gluon_dist.size(), gluon_rad, x_center, y_center);

  for (int i = gluon_bounds[0]; i < gluon_bounds[2]; i++)
  {
    for (int j = gluon_bounds[1]; j < gluon_bounds[3]; j++)
    {
      if (test_ == "GreensFunction")
      {
        // this must be updated to only subtract energy from the initial condition since the energy is already present in the output
        //  Change energy totals to reflect change in energy grids
        total_initial_energy -= gluon_dist[i][j]*ratio*initial_energy[x_center - gluon_rad + i][y_center - gluon_rad + j];
        total_energy += gluon_dist[i][j]*ratio*initial_energy[x_center - gluon_rad + i][y_center - gluon_rad + j];

        //  Subtract energy proportional to ratio from initial_energy and add it to density[0]
        initial_energy[x_center - gluon_rad + i][y_center - gluon_rad + j] -= gluon_dist[i][j]*ratio*initial_energy[x_center - gluon_rad + i][y_center - gluon_rad + j];

        continue;
      }

      //  Change energy totals to reflect change in energy grids
      total_initial_energy -= gluon_dist[i][j]*ratio*initial_energy[x_center - gluon_rad + i][y_center - gluon_rad + j];
      total_energy += gluon_dist[i][j]*ratio*initial_energy[x_center - gluon_rad + i][y_center - gluon_rad + j];

      //  Subtract energy proportional to ratio from initial_energy and add it to density[0]
      density[0][x_center - gluon_rad + i][y_center - gluon_rad + j] += gluon_dist[i][j]*ratio*initial_energy[x_center - gluon_rad + i][y_center - gluon_rad + j];
      initial_energy[x_center - gluon_rad + i][y_center - gluon_rad + j] -= gluon_dist[i][j]*ratio*initial_energy[x_center - gluon_rad + i][y_center - gluon_rad + j];
    }
  }
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Gets intigration bounds for density grid manipulations
//##########################################################################################
vector<int> Event::GetIntegrationBounds(int size, double raduis, double xcenter, double ycenter)
{
  vector<int> bounds;

  //  Set default values for integration bounds
  bounds.push_back(0);
  bounds.push_back(0);
  bounds.push_back(size);
  bounds.push_back(size);

  //  Update lower bounds of integration to avoid accessing elements out of bounds of density grids
  if (xcenter - raduis < 0)
  { bounds[0] = -(xcenter - raduis);  }
  if (ycenter - gluon_rad < 0)
  { bounds[1] = -(ycenter - raduis);  }

  //  Update upper bounds of integration to avoid accessing elements out of bounds of density grids
  if (xcenter + raduis > t_b.size())
  { bounds[2] = grid_points - (xcenter + raduis);  }
  if (ycenter + gluon_rad > t_b.size())
  { bounds[3] = grid_points - (xcenter + raduis);  }

  return bounds;
}
//__________________________________________________________________________________________

//__________________________________________________________________________________________
//##########################################################################################
//  Checks Event totals and returns true when initial_total is below a threshold
//##########################################################################################
bool Event::IsEventDone()
{

//  if (test_ == "GreensFunction" && number_strange == 1)
//  { return true;  }

  if (test_ == "SChop")
  {
    for (int i = 0; i < valued_points.size(); i++)
    {
      density[0][valued_points[i][0]][valued_points[i][1]] = initial_energy[valued_points[i][0]][valued_points[i][1]];
      initial_energy[valued_points[i][0]][valued_points[i][1]] = 0;
      valued_points.erase(valued_points.begin() + i);
    }
    return true;
  }

  //  Check all points in valued points and remove ones that are now 0
  for (int i = 0; i < valued_points.size(); i++)
  {
    if (initial_energy[valued_points[i][0]][valued_points[i][1]] == 0)
    { valued_points.erase(valued_points.begin() + i); }
  }

  //  If there are no valued points left, exit algorithm loop and write event to file
  if (valued_points.size() == 0)
  { return true;  }
  else
  { return false; }
}

//__________________________________________________________________________________________
//##########################################################################################
//  Clears Event variables as a cautionary measure
//##########################################################################################
void Event::CleanEvent()
{
  initial_energy.clear();
  t_a.clear();
  t_b.clear();
  density.clear();
  gluon_dist.clear();
  quark_dist.clear();
  eccentricities.clear();
  ecc.CleanEccentricity();

  total_initial_energy = 0;
  total_energy = 0;
  seed = 0;

  initial_energy_backup.clear();
  final_energy_backup.clear();
}
//__________________________________________________________________________________________
