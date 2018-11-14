#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include "GlobalVariables.C"
#include <fun4all/Fun4AllServer.h>
#include <g4detectors/PHG4MapsCellReco.h>
#include <g4detectors/PHG4MapsSubsystem.h>
#include <g4detectors/PHG4SiliconTrackerCellReco.h>
#include <g4detectors/PHG4SiliconTrackerDefs.h>
#include <g4detectors/PHG4SiliconTrackerSubsystem.h>
#include <g4detectors/PHG4TPCSpaceChargeDistortion.h>
#include <g4eval/SvtxEvaluator.h>
#include <g4hough/PHG4GenFitTrackProjection.h>
#include <g4hough/PHG4KalmanPatRec.h>
#include <g4hough/PHG4SiliconTrackerDigitizer.h>
#include <g4hough/PHG4SvtxClusterizer.h>
#include <g4hough/PHG4SvtxDeadArea.h>
#include <g4hough/PHG4SvtxDigitizer.h>
#include <g4hough/PHG4SvtxThresholds.h>
#include <g4hough/PHG4TPCClusterizer.h>
#include <g4hough/PHG4TrackKalmanFitter.h>
#include <g4hough/PHG4TruthPatRec.h>
#include <g4main/PHG4Reco.h>
#include <g4tpc/PHG4TPCElectronDrift.h>
#include <g4tpc/PHG4TPCPadPlane.h>
#include <g4tpc/PHG4TPCPadPlaneReadout.h>
#include <g4tpc/PHG4TPCSubsystem.h>
R__LOAD_LIBRARY(libg4tpc.so)
R__LOAD_LIBRARY(libg4hough.so)
R__LOAD_LIBRARY(libg4eval.so)
#endif

#include <vector>
// define INTTLADDER8, INTTLADDER6, INTTLADDER4_ZP or INTTLADDER4_PP, INTTLADDER0 to get 8, 6, 4 or 0 layers
// one and only one of these has to be defined, because #elseif does not seem to work properly in the interpreter
#define INTTLADDER4_PP

// Dead map options for INTT
enum enu_INTTDeadMapType
{
  kINTTNoDeadMap = 0,       // All channel in INTT is alive
  kINTT4PercentDeadMap = 4, // 4% of dead/masked area (2% sensor + 2% chip) as a typical FVTX Run14 production run.
  kINTT8PercentDeadMap = 8  // 8% dead/masked area (6% sensor + 2% chip) as threshold of operational
};

// Choose INTT deadmap here
enu_INTTDeadMapType INTTDeadMapOption = kINTTNoDeadMap;

// ONLY if backward compatibility with hits files already generated with 8 inner TPC layers is needed, you can set this to "true"
bool tpc_layers_40  = false;

// if true, refit tracks with primary vertex included in track fit  - good for analysis of prompt tracks only
// Adds second node to node tree, keeps original track node undisturbed
// Adds second evaluator to process refitted tracks and outputs separate ntuples
bool use_primary_vertex = false;

const int n_maps_layer = 3;  // must be 0-3, setting it to zero removes MVTX completely, n < 3 gives the first n layers

// Configure the INTT layers
// offsetphi is in deg, every other layer is offset by one half of the phi spacing between ladders
#ifdef INTTLADDER8
int n_intt_layer = 8;  
// default layer configuration
int laddertype[8] = {PHG4SiliconTrackerDefs::SEGMENTATION_Z, 
		     PHG4SiliconTrackerDefs::SEGMENTATION_Z, 
		     PHG4SiliconTrackerDefs::SEGMENTATION_PHI,
		     PHG4SiliconTrackerDefs::SEGMENTATION_PHI,
		     PHG4SiliconTrackerDefs::SEGMENTATION_PHI,
		     PHG4SiliconTrackerDefs::SEGMENTATION_PHI,
		     PHG4SiliconTrackerDefs::SEGMENTATION_PHI,
		     PHG4SiliconTrackerDefs::SEGMENTATION_PHI};  // default
int nladder[8] = {17,  17, 15, 15, 18, 18, 21, 21};  // default
double sensor_radius[8] = {6.876, 7.462, 8.987, 9.545, 10.835, 11.361, 12.676, 13.179};  // radius of center of sensor for layer default
double offsetphi[8] = {0.0, 0.5 * 360.0 / nladder[1] , 0.0, 0.5 * 360.0 / nladder[3], 0.0, 0.5 * 360.0 / nladder[5], 0.0, 0.5 * 360.0 / nladder[7]};
#endif
#ifdef INTTLADDER6
int n_intt_layer = 6;
int laddertype[6] = {PHG4SiliconTrackerDefs::SEGMENTATION_Z, 
		       PHG4SiliconTrackerDefs::SEGMENTATION_Z, 
		       PHG4SiliconTrackerDefs::SEGMENTATION_PHI,
		       PHG4SiliconTrackerDefs::SEGMENTATION_PHI,
		       PHG4SiliconTrackerDefs::SEGMENTATION_PHI,
		       PHG4SiliconTrackerDefs::SEGMENTATION_PHI };
int nladder[6] = {17,  17, 15, 15, 18, 18}; 
double sensor_radius[6] = {6.876, 7.462, 8.987, 9.545, 10.835, 11.361};  // radius of center of sensor for layer default
double offsetphi[6] = {0.0, 0.5 * 360.0 / nladder[1] , 0.0, 0.5 * 360.0 / nladder[3], 0.0, 0.5 * 360.0 / nladder[5]};
#endif
#ifdef INTTLADDER4_ZP
int n_intt_layer = 4;
int laddertype[4] = {PHG4SiliconTrackerDefs::SEGMENTATION_Z, 
		       PHG4SiliconTrackerDefs::SEGMENTATION_Z, 
		       PHG4SiliconTrackerDefs::SEGMENTATION_PHI,
		       PHG4SiliconTrackerDefs::SEGMENTATION_PHI};
int nladder[4] = {17,  17, 18, 18}; 
double sensor_radius[6] = {6.876, 7.462, 10.835, 11.361};  // radius of center of sensor for layer default
double offsetphi[6] = {0.0, 0.5 * 360.0 / nladder[1] , 0.0, 0.5 * 360.0 / nladder[3]};
#endif
#ifdef INTTLADDER4_PP
int n_intt_layer = 4;
int laddertype[4] = {PHG4SiliconTrackerDefs::SEGMENTATION_PHI, 
		       PHG4SiliconTrackerDefs::SEGMENTATION_PHI, 
		       PHG4SiliconTrackerDefs::SEGMENTATION_PHI,
		       PHG4SiliconTrackerDefs::SEGMENTATION_PHI};
int nladder[4] = {15,  15, 18, 18}; 
double sensor_radius[6] = { 8.987, 9.545, 10.835, 11.361};  // radius of center of sensor for layer default
double offsetphi[6] = {0.0, 0.5 * 360.0 / nladder[1] , 0.0, 0.5 * 360.0 / nladder[3]};
#endif
#ifdef INTTLADDER0
int n_intt_layer = 0;
#endif

int n_tpc_layer_inner = 16;
int tpc_layer_rphi_count_inner = 1152;
int n_tpc_layer_mid = 16;
int n_tpc_layer_outer = 16;
int n_gas_layer = n_tpc_layer_inner + n_tpc_layer_mid + n_tpc_layer_outer;

int Max_si_layer;

void TrackingInit(int verbosity = 0)
{
  Max_si_layer = n_maps_layer + n_intt_layer + n_gas_layer;
}

double Tracking(PHG4Reco* g4Reco, double radius,
            const int absorberactive = 0,
            int verbosity = 0)
{
  // create the three tracker subsystems 

  PHG4CylinderSubsystem* cyl;

  if (n_maps_layer > 0)
    {
      bool maps_overlapcheck = false;  // set to true if you want to check for overlaps
      
      // MAPS inner barrel layers
      //======================================================
      
      double maps_layer_radius[3] = {24.61, 32.59, 39.88}; // mm - numbers from Walt 6 Aug 2018
      
      // D. McGlinchey 6Aug2018 - type no longer is used, included here because I was too lazy to remove it from the code
      int stave_type[3] = {0, 0, 0};
      int staves_in_layer[3] = {12, 16, 20};       // Number of staves per layer in sPHENIX MVTX
      double phi_tilt[3] = {0.300, 0.305, 0.300}; // radians - numbers from Walt 6 Aug 2018
      
      for (int ilayer = 0; ilayer < n_maps_layer; ilayer++)
	{
//if (verbosity)
//	    cout << "Create Maps layer " << ilayer << " with radius " << maps_layer_radius[ilayer] << " mm, stave type " << stave_type[ilayer]
//		 << " pixel size 30 x 30 microns "
//		 << " active pixel thickness 0.0018 microns" << endl;
//
//	  PHG4MapsSubsystem* lyr = new PHG4MapsSubsystem("MAPS", ilayer, stave_type[ilayer]);
//	  lyr->Verbosity(verbosity);
//
//	  lyr->set_double_param("layer_nominal_radius", maps_layer_radius[ilayer]);  // thickness in cm
//	  lyr->set_int_param("N_staves", staves_in_layer[ilayer]);       // uses fixed number of staves regardless of radius, if set. Otherwise, calculates optimum number of staves
//
//	  // The cell size is used only during pixilization of sensor hits, but it is convemient to set it now because the geometry object needs it
//	  lyr->set_double_param("pixel_x", 0.0030);          // pitch in cm
//	  lyr->set_double_param("pixel_z", 0.0030);          // length in cm
//	  lyr->set_double_param("pixel_thickness", 0.0018);  // thickness in cm
//	  lyr->set_double_param("phitilt", phi_tilt[ilayer]);
//
//	  lyr->set_int_param("active", 1);
//	  lyr->OverlapCheck(maps_overlapcheck);
//
//	  lyr->set_string_param("stave_geometry_file", string(getenv("CALIBRATIONROOT")) + string("/Tracking/geometry/mvtx_stave_v01.gdml"));
//
//	  g4Reco->registerSubsystem(lyr);

        cyl = new PHG4CylinderSubsystem("MVTX",ilayer);
        cyl->set_double_param("radius",  maps_layer_radius[ilayer] / 10);
        cyl->set_int_param("lengthviarapidity", 0);
        cyl->set_double_param("length", 30);
        cyl->set_string_param("material", "G4_AIR");
        cyl->set_double_param("thickness", 50e-4);
        cyl->SetActive();
        cyl->SuperDetector("MVTX");
        cyl->Verbosity(0);
        g4Reco->registerSubsystem(cyl);

	  radius = maps_layer_radius[ilayer];
	}
    }
  
  if (n_intt_layer > 0)
    {
      //-------------------
      // INTT ladders
      //-------------------
      
      bool intt_overlapcheck = false;  // set to true if you want to check for overlaps
      
      // instantiate the Silicon tracker subsystem and register it
      // We make one instance of PHG4TrackerSubsystem for all four layers of tracker
      // dimensions are in mm, angles are in radians
      
      // PHG4SiliconTrackerSubsystem creates the detetor layer using PHG4SiliconTrackerDetector
      // and instantiates the appropriate PHG4SteppingAction
      const double intt_radius_max = 140.;  // including stagger radius (mm)
      
      // The length of vpair is used to determine the number of layers
      std::vector<std::pair<int, int>> vpair;  // (sphxlayer, inttlayer)
      for (int i = 0; i < n_intt_layer; i++)
	{
	  // We want the sPHENIX layer numbers for the INTT to be from n_maps_layer to n_maps_layer+n_intt_layer - 1
	  vpair.push_back(std::make_pair(n_maps_layer + i, i));  // sphxlayer=n_maps_layer+i corresponding to inttlayer=i
	  if (verbosity) cout << "Create strip tracker layer " << vpair[i].second << " as  sphenix layer  " << vpair[i].first << endl;
	}
      
//      PHG4SiliconTrackerSubsystem* sitrack = new PHG4SiliconTrackerSubsystem("SILICON_TRACKER", vpair);
//      sitrack->Verbosity(verbosity);
//      sitrack->SetActive(1);
//      sitrack->OverlapCheck(intt_overlapcheck);
//      g4Reco->registerSubsystem(sitrack);
      
      // Set the laddertype and ladder spacing configuration
      cout << "INTT has " << n_intt_layer << " layers with layer setup:" << endl;
      for(int i=0;i<n_intt_layer;i++)
	{
//	  cout << " INTT layer " << i << " laddertype " << laddertype[i] << " nladders " << nladder[i]
//	       << " sensor radius " << sensor_radius[i] << " offsetphi " << offsetphi[i] << endl;
//	  sitrack->set_int_param(i, "laddertype", laddertype[i]);
//	  sitrack->set_int_param(i, "nladder", nladder[i]);
//	  sitrack->set_double_param(i,"sensor_radius", sensor_radius[i]);  // expecting cm
//	  sitrack->set_double_param(i,"offsetphi",offsetphi[i]);  // expecting degrees


        cyl = new PHG4CylinderSubsystem("SILICON_TRACKER",n_maps_layer + i);
        cyl->set_double_param("radius",  sensor_radius[i]);
        cyl->set_int_param("lengthviarapidity", 0);
        cyl->set_double_param("length", 60);
        cyl->set_string_param("material", "G4_AIR"); // define material here
        cyl->set_double_param("thickness", 300e-4); // define thickness here, e.g. 300um
        cyl->SuperDetector("SILICON_TRACKER");
        cyl->SetActive();
        cyl->Verbosity(0);
        g4Reco->registerSubsystem(cyl);


	}
      
      // outer radius marker (translation back to cm)
      radius = intt_radius_max * 0.1;
    }

  // The TPC - always present!
  //================================

  double inner_cage_radius = 20.;
  double inner_readout_radius = 30.;


  radius = inner_cage_radius;

  double cage_length = 211.0;  // From TPC group, gives eta = 1.1 at 78 cm
  double n_rad_length_cage = 1.13e-02;
  double cage_thickness = 28.6 * n_rad_length_cage;  // Kapton X_0 = 28.6 cm  // mocks up Kapton + carbon fiber structure


  double inner_readout_radius = radius;
//  if (inner_readout_radius < radius) inner_readout_radius = radius;
//
  string tpcgas = "sPHENIX_TPC_Gas";  //  Ne(90%) CF4(10%) - defined in g4main/PHG4Reco.cc
//  radius = inner_readout_radius;

  double outer_radius = 78.;

  int n_tpc_layer_inner = 16;
  double tpc_layer_thick_inner = 1.25; // EIC- recover default inner radius of TPC vol.
  int tpc_layer_rphi_count_inner = 1152;

  int n_tpc_layer_mid = 16;
  double tpc_layer_thick_mid = 1.25;
  int tpc_layer_rphi_count_mid = 1536;

  int n_tpc_layer_outer = 16;
  double tpc_layer_thick_outer = 1.125; // outer later reach from 60-78 cm (instead of 80 cm), that leads to radial thickness of 1.125 cm
  int tpc_layer_rphi_count_outer = 2304;

  int n_gas_layer = n_tpc_layer_inner + n_tpc_layer_mid + n_tpc_layer_outer;

  double inner_cage_radius = 20.;
  double inner_readout_radius = 30.;

  // inner field cage
  cyl = new PHG4CylinderSubsystem("SVTXSUPPORT", n_maps_layer + n_intt_layer);
  cyl->set_double_param("radius", radius);
  cyl->set_int_param("lengthviarapidity", 0);
  cyl->set_double_param("length", cage_length);
  cyl->set_string_param("material", "G4_KAPTON");
  cyl->set_double_param("thickness", cage_thickness);
  cyl->SuperDetector("SVTXSUPPORT");
  cyl->Verbosity(0);
  g4Reco->registerSubsystem(cyl);

  radius += cage_thickness;
  
  // Active layers of the TPC from 30-40 cm (inner layers)

  for (int ilayer = n_maps_layer + n_intt_layer; ilayer < (n_maps_layer + n_intt_layer + n_tpc_layer_inner); ++ilayer)
  {
    if (verbosity)
      cout << "Create TPC gas layer " << ilayer << " with inner radius " << radius << " cm "
           << " thickness " << tpc_layer_thick_inner - 0.01 << " length " << cage_length << endl;

    cyl = new PHG4CylinderSubsystem("SVTX", ilayer);
    cyl->set_double_param("radius", radius);
    cyl->set_int_param("lengthviarapidity", 0);
    cyl->set_double_param("length", cage_length);
    cyl->set_string_param("material", tpcgas.c_str());
    cyl->set_double_param("thickness", tpc_layer_thick_inner - 0.01);
    cyl->SetActive();
    cyl->SuperDetector("SVTX");
    g4Reco->registerSubsystem(cyl);

    radius += tpc_layer_thick_inner;
  }

  // Active layers of the TPC from 40-60 cm (mid layers)

  for (int ilayer = n_maps_layer + n_intt_layer + n_tpc_layer_inner; ilayer < (n_maps_layer + n_intt_layer + n_tpc_layer_inner + n_tpc_layer_mid); ++ilayer)
  {
    if (verbosity)
      cout << "Create TPC gas layer " << ilayer << " with inner radius " << radius << " cm "
           << " thickness " << tpc_layer_thick_mid - 0.01 << " length " << cage_length << endl;

    cyl = new PHG4CylinderSubsystem("SVTX", ilayer);
    cyl->set_double_param("radius", radius);
    cyl->set_int_param("lengthviarapidity", 0);
    cyl->set_double_param("length", cage_length);
    cyl->set_string_param("material", tpcgas.c_str());
    cyl->set_double_param("thickness", tpc_layer_thick_mid - 0.01);
    cyl->SetActive();
    cyl->SuperDetector("SVTX");
    g4Reco->registerSubsystem(cyl);

    radius += tpc_layer_thick_mid;
  }

  // Active layers of the TPC from 60-80 cm (outer layers)

  for (int ilayer = n_maps_layer + n_intt_layer + n_tpc_layer_inner + n_tpc_layer_mid; ilayer < (n_maps_layer + n_intt_layer + n_tpc_layer_inner + n_tpc_layer_mid + n_tpc_layer_outer); ++ilayer)
  {
    if (verbosity)
      cout << "Create TPC gas layer " << ilayer << " with inner radius " << radius << " cm "
           << " thickness " << tpc_layer_thick_outer - 0.01 << " length " << cage_length << endl;

    cyl = new PHG4CylinderSubsystem("SVTX", ilayer);
    cyl->set_double_param("radius", radius);
    cyl->set_int_param("lengthviarapidity", 0);
    cyl->set_double_param("length", cage_length);
    cyl->set_string_param("material", tpcgas.c_str());
    cyl->set_double_param("thickness", tpc_layer_thick_outer - 0.01);
    cyl->SetActive();
    cyl->SuperDetector("SVTX");
    g4Reco->registerSubsystem(cyl);

    radius += tpc_layer_thick_outer;
  }

  // outer field cage
  cyl = new PHG4CylinderSubsystem("SVTXSUPPORT", n_maps_layer + n_intt_layer + n_gas_layer);
  cyl->set_double_param("radius", radius);
  cyl->set_int_param("lengthviarapidity", 0);
  cyl->set_double_param("length", cage_length);
  cyl->set_string_param("material", "G4_KAPTON");
  cyl->set_double_param("thickness", cage_thickness);  // Kapton X_0 = 28.6 cm
  cyl->SuperDetector("SVTXSUPPORT");
  g4Reco->registerSubsystem(cyl);

  radius += cage_thickness;
  radius += no_overlapp;
  
  return radius; 
}

void Tracking_Cells(int verbosity = 0)
{
  return;
}

void Tracking_Reco(int verbosity = 0)
{
  //---------------
  // Load libraries
  //---------------

  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4hough.so");

  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer* se = Fun4AllServer::instance();


  //---------------------
  // Kalman Filter
  //---------------------

  PHG4TrackFastSim* kalman = new PHG4TrackFastSim("PHG4TrackFastSim");
  kalman->Verbosity(0);

  kalman->set_use_vertex_in_fitting(false);
//  kalman->set_vertex_xy_resolution(50E-4);
//  kalman->set_vertex_z_resolution(50E-4);

  std::string phg4hits_names[] = {"G4HIT_MVTX", "G4HIT_SILICON_TRACKER"};
  const PHG4TrackFastSim::DETECTOR_TYPE dettypes[] = {PHG4TrackFastSim::Cylinder,PHG4TrackFastSim::Cylinder};
  float rad[] = {5.0,  86}; // setting for cluster radial resolution in um
  float phi[] = {5.0,  24}; // setting for cluster azimuthal resolution in um
  float lon[] = {5.0,  4600}; // setting for cluster z (longitudinal) resolution in um
  for(int i=0; i<2; ++i) { // from um to cm
    rad[i] *= 1e-4;
    phi[i] *= 1e-4;
    lon[i] *= 1e-4;
  }
  float eff[] = {1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0}; // efficiency
  float noi[] = {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}; // Noise
  kalman->set_phg4hits_names(phg4hits_names, dettypes, rad, phi, lon, eff, noi, 2);

  //std::string phg4hits_names[] = {"G4HIT_EGEM_0","G4HIT_EGEM_1","G4HIT_EGEM_2","G4HIT_EGEM_3","G4HIT_FGEM_0","G4HIT_FGEM_1","G4HIT_FGEM_2","G4HIT_FGEM_3","G4HIT_FGEM_4"};
  //kalman->set_phg4hits_names(phg4hits_names, 9);
  kalman->set_sub_top_node_name("SVTX");
  kalman->set_trackmap_out_name("SvtxTrackMap");

  //std::string state_names[] = {"FEMC","FHCAL"};
  //kalman->set_state_names(state_names, 2);

  kalman->set_fit_alg_name("KalmanFitterRefTrack");//
  kalman->set_primary_assumption_pid(13);
  kalman->set_do_evt_display(false);

  se->registerSubsystem(kalman);


  return;
}

void Tracking_Eval(std::string outputfile, int verbosity = 0)
{
  //---------------
  // Load libraries
  //---------------

  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4detectors.so");
  gSystem->Load("libg4hough.so");
  gSystem->Load("libg4eval.so");

  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer* se = Fun4AllServer::instance();

  //----------------
  // Tracking evaluation
  //----------------

//  SvtxEvaluator* eval;
//  eval = new SvtxEvaluator("SVTXEVALUATOR", outputfile.c_str(), "SvtxTrackMap", n_maps_layer, n_intt_layer, n_gas_layer);
//  eval->do_cluster_eval(false);
//  eval->do_g4hit_eval(false);
//  eval->do_hit_eval(false);  // enable to see the hits that includes the chamber physics...
//  eval->do_gpoint_eval(false);
//  eval->scan_for_embedded(false);  // take all tracks if false - take only embedded tracks if true
//  eval->Verbosity(0);
//  se->registerSubsystem(eval);


  if (gSystem->Load("libFastTrackingEval.so") !=0)
  {
    cout <<"-------------------------------------------------------------------------------------------"<<endl;
    cout <<"Tracking_Eval - please build the FastTrackingEval module for analyzing the fast tracking results."
        <<" The module is located at analysis repository : https://github.com/sPHENIX-Collaboration/analysis/tree/master/Tracking/FastTrackingEval"
        <<endl;
    cout <<"-------------------------------------------------------------------------------------------"<<endl;
    exit(1);
  }

  Fun4AllServer *se = Fun4AllServer::instance();

  FastTrackingEval *fast_sim_eval = new FastTrackingEval("FastTrackingEval");
  fast_sim_eval->set_filename( outputfile.c_str() );
  se->registerSubsystem( fast_sim_eval );


  return;
}
