// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#if defined USE_GEAR

// EUTelescope includes:
#include "EUTelAnalysisCMSPixel.h"
#include "EUTelAnalysisCMSPixelInstance.h"
#include "EUTELESCOPE.h"
//#include "EUTelSparseDataImpl.h"
#include "EUTelBaseSparsePixel.h"
#include "EUTelGenericSparsePixel.h"
#include "EUTelTrackerDataInterfacerImpl.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelPStream.h" // process streams redi::ipstream

// AIDA histogram package (on top of ROOT):
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IAxis.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/IProfile2D.h>
#include <AIDA/ITree.h>
#endif

// GBL:
#include "include/GblTrajectory.h"
#include "include/MilleBinary.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/ProcessorMgr.h"
#include "marlin/Global.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// LCIO includes <.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>

// system includes <>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <memory>
#include <string.h>
#include <map>
#include <cstdlib>
#include <limits>
#include <iterator>

// ROOT includes ".h"
#include <TMath.h>
#include <TVectorD.h>
#include <TF1.h>
#include <TMatrixD.h>
#include <TVector3.h>
#include <TRotation.h>
#include "TH1D.h"

// Include the CMSPixelDecoder class:
#include "CMSPixelDecoder.h"

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;


int n_clusters_dut = 0;
int n_clusters_ref = 0;
int n_matched_clusters_dut = 0;
int n_matched_clusters_ref = 0;
int naldut = 0;
int ngbl = 0;
bool ldb = 0;
gbl::MilleBinary * mille; // for producing MillePede-II binary file

double turn = 0;
double tilt = 0;
double DUTz = 40;
double DUTalignx = 0;
double DUTaligny = 0;
double DUTrot = 0;


EUTelAnalysisCMSPixel::EUTelAnalysisCMSPixel() : Processor("EUTelAnalysisCMSPixel"), _siPlanesParameters(), _siPlanesLayerLayout(), _inputCollectionTelescope(""), _inputCollectionDUT(""), _inputCollectionREF(""), _inputTrackCollection(""), _isFirstEvent(0), _eBeam(0), _nEvt(0), _nTelPlanes(0), time_event0(0), time_event1(0), time_reference(0), fTLU(0), gTLU(0), _DUT_chip(0), _DUT_gain(""), _DUT_conversion(0), _DUT_calibration_type(""), dut_calibration(), _DUTalignx(0), _DUTaligny(0), _DUTz(0), _DUTrot(0), _DUTtilt(0), _DUTturn(0), _REF_chip(0), _REF_gain(""), _REF_calibration_type(""), ref_calibration(), _REFalignx(0), _REFaligny(0), _REFz(0), _REFrot(0), _cutx(0.15), _cuty(0.1), _skew_db(""), _have_skew_db(false), skew_par0(0), skew_par1(0), _CMS_gain_path(""), _adc_correction(0), _gearfile(""), _alignmentrun(""), _planeSort(), _planeID(), _planePosition(), _planeThickness(), _planeX0(), _planeResolution(), _skip_dut(0), _skip_ref(0), _skip_tel(0), dut_event_buffer(), ref_event_buffer(), tel_event_buffer(), ClustDUT(), ClustREF(), m_millefilename(""), _skip_db("") , _have_skip_db(false), _uptimestart(), _uptimeend() {

  // modify processor description
  _description = "Analysis for CMS PSI46 Pixel Detectors as DUT in AIDA telescopes ";

  // processor parameters
  registerInputCollection( LCIO::TRACKERHIT,
                           "InputCollectionTelescope" ,
                           "Name of the input TrackerHit collection of the telescope",
                           _inputCollectionTelescope,
                           std::string("alignedhits") );
  registerInputCollection( LCIO::TRACKERHIT,
                           "InputCollectionDUT" ,
                           "Name of the input TrackerHit collection of the DUT",
                           _inputCollectionDUT,
                           std::string("duthits") );
  registerInputCollection( LCIO::TRACKERHIT,
                           "InputCollectionREF" ,
                           "Name of the input TrackerHit collection of the REF",
                           _inputCollectionREF,
                           std::string("refhits") );
  registerInputCollection( LCIO::TRACKERDATA,
                           "InputPixelsDUT" ,
                           "Name of the input TrackerData (pixels) collection of the DUT",
                           _inputPixelsDUT,
                           std::string("dutpixels"));
  registerInputCollection( LCIO::TRACKERDATA,
                           "InputPixelsREF" ,
                           "Name of the input TrackerData (pixels) collection of the REF",
                           _inputPixelsREF,
                           std::string("refpixels"));
  registerOutputCollection(LCIO::TRACK,"InputTrackCollectionTelescope",
                           "Name of the input Track collection of the telescope",
                           _inputTrackCollection, std::string( "fittracks"));

  registerProcessorParameter( "Ebeam",
                              "Beam energy [GeV]",
			      _eBeam, static_cast < double >( 4.0));

  registerProcessorParameter("CMS_gain_path",
                             "Path to the CMS gain calibration files",
                             _CMS_gain_path, std::string("/dev/null"));


  registerProcessorParameter("DUT_chip",
                             "Serial number of the DUT used in this run",
                             _DUT_chip, static_cast< int > (0) );
  registerProcessorParameter("DUT_gain",
                             "CMS DUT gain file to be used",
                             _DUT_gain, std::string("/dev/null"));
  registerProcessorParameter("DUT_conversion",
                             "CMS DUT VCAL->electrons conversion factor. Set to 0 if default is to be used.",
                             _DUT_conversion, static_cast<double>(0.));
  registerProcessorParameter("DUT_calibration",
                             "Choose DUT calibration type (psi_tanh, desy_tanh, desy_weibull)",
                             _DUT_calibration_type, std::string("desy_weibull"));

  registerProcessorParameter("REF_chip",
                             "Serial number of the REF used in this run",
                             _REF_chip, static_cast< int > (0) );
  registerProcessorParameter("REF_gain",
                             "CMS REF gain file to be used",
                             _REF_gain, std::string("/dev/null"));
  registerProcessorParameter("REF_calibration",
                             "Choose REF calibration type (psi_tanh, desy_tanh, desy_weibull)",
                             _REF_calibration_type, std::string("desy_weibull"));

  // Alignment stuff:

  registerProcessorParameter( "DUT_align_x",
                              "DUT alignment in x",
			      _DUTalignx, static_cast < double >(0));
  registerProcessorParameter( "DUT_align_y",
                              "DUT alignment in y",
			      _DUTaligny, static_cast < double >(0));
  registerProcessorParameter( "DUT_rot",
                              "DUT rotation in xy",
			      _DUTrot, static_cast < double >(0));
  registerProcessorParameter( "DUT_tilt",
                              "DUT tilt angle",
			      _DUTtilt, static_cast < double >(0));
  registerProcessorParameter( "DUT_turn",
                              "DUT turn angle",
			      _DUTturn, static_cast < double >(0));
  registerProcessorParameter( "DUT_pos_z",
                              "DUT position z in mm after the telescope plane",
			      _DUTz, static_cast < double >(0));

  registerProcessorParameter( "REF_align_x",
                              "REF alignment in x",
			      _REFalignx, static_cast < double >(0));
  registerProcessorParameter( "REF_align_y",
                              "REF alignment in y",
			      _REFaligny, static_cast < double >(0));
  registerProcessorParameter( "REF_rot",
                              "REF rotation in xy",
			      _REFrot, static_cast < double >(0));
  registerProcessorParameter( "REF_pos_z",
                              "REF position z in mm after the last telescope plane",
			      _REFz, static_cast < double >(0));

  registerProcessorParameter( "matching_cut_x",
                              "cut for matching in x coordinate in mm",
			      _cutx, static_cast < double >(0.15));
  registerProcessorParameter( "matching_cut_y",
                              "cut for matching in y coordinate in mm",
			      _cuty, static_cast < double >(0.10));

  registerProcessorParameter( "slope_cut_x",
                              "cut for track slopes in x coordinate in rad",
			      _slope_x, static_cast < double >(0.002));
  registerProcessorParameter( "slope_cut_y",
                              "cut for track slopes in y coordinate in rad",
			      _slope_y, static_cast < double >(0.002));

  registerProcessorParameter( "skew_database",
                              "database file for skew corrections",
			      _skew_db, std::string("none"));

  registerProcessorParameter( "adc_correction",
                              "ADC value c orrection for all but first pixel",
			      _adc_correction, static_cast < double >(2.6));

  registerOptionalParameter( "skip_dut",
			     "Skip N events from DUT data stream at beginning",
			     _skip_dut, static_cast < int >(0) );
  registerOptionalParameter( "skip_ref",
			     "Skip N events from REL data stream at beginning",
			     _skip_ref, static_cast < int >(0) );
  registerOptionalParameter( "skip_tel",
			     "Skip N events from TEL data stream at beginning",
			     _skip_tel, static_cast < int >(0) );

  registerOptionalParameter( "skip_database",
			     "Database for skipping of faulty time regions",
			     _skip_db, std::string("none") );

  // Stuff only needed for the printout of the updated runlist line:
  registerOptionalParameter( "gearfile",
			     "Again, the gear file of the used telescope configuration",
			     _gearfile, std::string("none") );
  registerOptionalParameter( "alignmentrun",
			     "Run number of the telescope alignment run",
			     _alignmentrun, std::string("none") );
}


void EUTelAnalysisCMSPixel::init() {

  srand( time(NULL) ); //Randomize seed initialization

  // usually a good idea to
  printParameters();

  _nEvt = 0;

  _isFirstEvent = true;

  streamlog_out(MESSAGE0) << "Beam energy " << _eBeam << " GeV" <<  std::endl;

  // Read geometry information from GEAR

  streamlog_out( MESSAGE0 ) << "Reading telescope geometry description from GEAR " << std::endl;

  _siPlanesParameters  = const_cast<gear::SiPlanesParameters* >( &(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*>(  &(_siPlanesParameters->getSiPlanesLayerLayout() ));


  // Take all layers defined in GEAR geometry
  _nTelPlanes = _siPlanesLayerLayout->getNLayers();

  _planeSort = new int[_nTelPlanes];
  _planePosition = new double[_nTelPlanes]; // z pos
  _planeID         = new int[_nTelPlanes];
  _planeThickness  = new double[_nTelPlanes];
  _planeX0         = new double[_nTelPlanes];
  _planeResolution = new double[_nTelPlanes];

  for(int i = 0; i < _nTelPlanes; i++) {
    
    _planeID[i]=_siPlanesLayerLayout->getID(i);
    _planePosition[i]=_siPlanesLayerLayout->getLayerPositionZ(i);
    _planeThickness[i]=_siPlanesLayerLayout->getLayerThickness(i);
    _planeX0[i]=_siPlanesLayerLayout->getLayerRadLength(i);
    _planeResolution[i] = _siPlanesLayerLayout->getSensitiveResolution(i);
  }

  streamlog_out( MESSAGE2 ) <<  "Telescope configuration with " << _nTelPlanes << " planes" << std::endl;

  for( int ipl = 0; ipl < _nTelPlanes; ipl++) {
    std::stringstream ss;
    ss << "  ID = " << _planeID[ipl]
       << "  at Z [mm] = " << _planePosition[ipl]
       << " dZ [um] = " << _planeThickness[ipl]*1000.;

    ss << "  Res [um] = " << _planeResolution[ipl]*1000.;

    streamlog_out( MESSAGE2 ) <<  ss.str() << std::endl;
    
  }

  // Check for abnormalities in DUT position:
  if(_DUTz <= 0.001) { 
    streamlog_out(ERROR) << "Your DUT seems to be at the same z position as the 3rd telescope plane." << std:: endl << "Please adjust DUT_pos_z." << std::endl;
    throw InvalidGeometryException("GEAR manager is not initialised");
  }

  // Book histograms:
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  bookHistos();
#endif

  // Read database with skew corrections:
  std::ifstream in;
  streamlog_out(MESSAGE2) << "Attempting to open skew correction database from " << endl;
  streamlog_out(MESSAGE2) << _skew_db << endl;
  in.open(_skew_db.c_str(), std::ifstream::in);

  bool have_lower_tilt = false;
  bool have_higher_tilt = false;
  double lower_skew_par0 = 0, lower_skew_par1 = 0, lower_tilt = 0;
  double higher_skew_par0 = 0, higher_skew_par1 = 0, higher_tilt = 0;

  if (in.is_open()) {
    string line;
    // Skip first line:
    getline(in,line);

    while(std::getline(in,line)){
      // Skip reading comments:
      if (line.empty() || line[0] == '#') continue;

      istringstream s(line);
      int i = 0;
      while (s) {
	string str;
	if(!getline( s, str, ',' )) break;
	if(i == 0) {  // Read tilt angle
	  if(atof(str.c_str()) < _DUTtilt) { have_lower_tilt = true; lower_tilt = atof(str.c_str()); }
	  if(atof(str.c_str()) > _DUTtilt) { have_higher_tilt = true; higher_tilt = atof(str.c_str()); }
	}
	if(have_higher_tilt) {
	  if(i == 1) { higher_skew_par0 = atof(str.c_str()); }
	  if(i == 2) { higher_skew_par1 = atof(str.c_str()); break; }
	}
	else if(have_lower_tilt) {
	  if(i == 1) { lower_skew_par0 = atof(str.c_str()); }
	  if(i == 2) { lower_skew_par1 = atof(str.c_str()); break; }
	}
	i++;
      }
      if(have_higher_tilt) { break; }
    }
    in.close();

    if(have_higher_tilt && have_lower_tilt) {
      streamlog_out(MESSAGE2) << "Skew correction parameters for tilt " << lower_tilt << " < " << _DUTtilt 
			      << ": " << lower_skew_par0 << " " << lower_skew_par1 << endl;
      streamlog_out(MESSAGE2) << "Skew correction parameters for tilt " << higher_tilt << " > " << _DUTtilt
			      << ": " << higher_skew_par0 << " " << higher_skew_par1 << endl;

      double slope1 = (higher_skew_par1-lower_skew_par1)/(higher_tilt-lower_tilt);
      skew_par1 = lower_skew_par1 + (_DUTtilt-lower_tilt)*slope1;

      double slope0 = (higher_skew_par0-lower_skew_par0)/(higher_tilt-lower_tilt);
      skew_par0 = lower_skew_par0 + (_DUTtilt-lower_tilt)*slope0;

      streamlog_out(MESSAGE2) << "Interpolated skew correction parameters for tilt " << _DUTtilt 
			      << ": " << skew_par0 << " " << skew_par1 << endl;
      _have_skew_db = true;
    }
  }
  else {
    streamlog_out(WARNING) << "Could not open skew correction database, no correction applied." << endl;    
  }

}//init

//------------------------------------------------------------------------------
void EUTelAnalysisCMSPixel::processRunHeader( LCRunHeader* runHeader) {

  std::auto_ptr<EUTelRunHeaderImpl> eutelHeader( new EUTelRunHeaderImpl( runHeader ) );
  eutelHeader->addProcessor( type() );

  // Decode and print out Run Header information - just a check

  _nRun = runHeader->getRunNumber();

  streamlog_out( MESSAGE2 )  << "Processing run header"
                             << ", run nr " << runHeader->getRunNumber() << std::endl;

  // Prepare the mille binary file:
  std::stringstream name;
  name << "run" << _nRun << "-mille.bin";
  m_millefilename = name.str();
  mille = new gbl::MilleBinary( m_millefilename.c_str() );


  // Scale cuty with the expected resolution (depending on tilt angle)
  // should always be around 10sigma
  if(_DUT_chip == 506) {
    if     (_DUTtilt < 14) { _cuty = 0.40; }
    else if(_DUTtilt < 19) { _cuty = 0.25; }
    else if(_DUTtilt < 35) { _cuty = 0.10; }
    else { _cuty = 0.1 + 0.00625*(_DUTtilt-35); }
  }
  else {
    if     (_DUTtilt < 7) { _cuty = 0.30; }
    else if(_DUTtilt < 15) { _cuty = 0.15; }
    else if(_DUTtilt < 25) { _cuty = 0.05; }
    else { _cuty = 0.10; }
  }
  streamlog_out(MESSAGE0) << "Set matching cut Y to " << _cuty << "mm" << std::endl;

  const std::string detectorName = runHeader->getDetectorName();
  const std::string detectorDescription = runHeader->getDescription();

  streamlog_out( MESSAGE0 ) << detectorName << " : " << detectorDescription << std::endl;

  // DUT calibration:
  if(!InitializeCalibration(_CMS_gain_path + "/chip" + ZeroPadNumber(_DUT_chip,3) + "/" + _DUT_gain,
			    _DUT_chip, _DUT_calibration_type, dut_calibration))
    throw marlin::ParseException("Error in processing DUT calibration file.");

  // REF calibration:
  if(!InitializeCalibration(_CMS_gain_path + "/chip" + ZeroPadNumber(_REF_chip,3) + "/" + _REF_gain,
			    _REF_chip, _REF_calibration_type, ref_calibration))
    throw marlin::ParseException("Error in processing REF calibration file.");

  // Get conversion factors to kiloelectrons:
  _DUT_conversion = GetConversionFactor(dut_calibration,_DUT_conversion);
  _REF_conversion = GetConversionFactor(ref_calibration);

  streamlog_out(MESSAGE0) << "DUT Conversion: " << _DUT_conversion << "ke" << endl;
  streamlog_out(MESSAGE0) << "REF Conversion: " << _REF_conversion << "ke" << endl;

  if(_DUT_chip == 506 && _nRun <= 15268) {
    // Deal with backwards-incompatible change in DAQ software:
    // before run 15268 the PCBTYPE flag has not been set and
    // col/row are not inverted:
    streamlog_out(WARNING) << "Old data run without PCBTYPE set by EUDAQ producer. Assuming rotated chip..." << endl;
  }

  
  // Read database from uptime recognition:

  if(_skip_db.at(_skip_db.size()-1) != '/'){
    _skip_db.append("/");
  }
  char buffer[100];
  sprintf(buffer, "run%06d-uptime.dat", _nRun);
  _skip_db.append(buffer);

  std::ifstream in2;
  streamlog_out(MESSAGE2) << "Attempting to open uptime database from " << endl;
  streamlog_out(MESSAGE2) << _skip_db << endl;

  in2.open(_skip_db.c_str(), std::ifstream::in);

  double upt, downt;
  if (in2.is_open()) {
    while(in2 >> upt >> downt){
      _uptimestart.push_back(upt);
      _uptimeend.push_back(downt);
    }
    _have_skip_db = true;

    streamlog_out(DEBUG2) << "uptime\tdowntime" << endl;

    for(int i = 0; i < _uptimestart.size(); i++){
      
      streamlog_out(DEBUG2) << _uptimestart.at(i) << "\t" << _uptimeend.at(i) << endl;
      
    }

    
  }else{
    streamlog_out(MESSAGE2) << "Could not open uptime database, no events skipped." << endl;
  }
  

} // processRunHeader

//----------------------------------------------------------------------------
void EUTelAnalysisCMSPixel::processEvent( LCEvent * event ) {

  EUTelEventImpl * euEvent = static_cast<EUTelEventImpl*>( event );

  if( euEvent->getEventType() == kEORE ) {
    streamlog_out( DEBUG5 ) <<  "EORE found: nothing else to do." << std::endl;
    return;
  }

  int runNumber = event->getRunNumber();

  if( _isFirstEvent ) {

    // Initialize the timing calculations:
    time_event0 = event->getTimeStamp();
    time_event1 = time_event0;
    time_reference = time_event0;

    // apply all GEAR/alignment offsets to get corrected X,Y,Z
    // position of the sensor center

    for( int iplane = 0; iplane < _siPlanesLayerLayout->getNLayers(); iplane++ ) {

      std::map< unsigned int , double > _planeCenter;
      std::map< unsigned int , double > _planeNormal;

      // start filling the map with Gear values:
      _planeCenter[ 0 ] =  _siPlanesLayerLayout->getLayerPositionX(iplane); // X
      _planeCenter[ 1 ] =  _siPlanesLayerLayout->getLayerPositionY(iplane); // Y
      _planeCenter[ 2 ] =  _siPlanesLayerLayout->getLayerPositionZ(iplane); // Z
      _planeNormal[ 0 ] =  0.; // X
      _planeNormal[ 1 ] =  0.; // Y
      _planeNormal[ 2 ] =  1.; // Z

      TVector3  _normalTVec( _planeNormal[0], _planeNormal[1], _planeNormal[2]); 

      // do initial rotation from GEAR
      try{
 	double gRotation[3] = { 0., 0., 0.}; // not rotated
                         
	gRotation[0] = _siPlanesLayerLayout->getLayerRotationXY(iplane); // Euler alpha
	gRotation[1] = _siPlanesLayerLayout->getLayerRotationZX(iplane); // Euler alpha
	gRotation[2] = _siPlanesLayerLayout->getLayerRotationZY(iplane); // Euler alpha
                          
	// input angles are in degree, translate into radian
	gRotation[0] =  gRotation[0]*3.1415926/180.;
	gRotation[1] =  gRotation[1]*3.1415926/180.;
	gRotation[2] =  gRotation[2]*3.1415926/180.;

	TRotation r;
	r.RotateX( gRotation[2] );
	r.RotateY( gRotation[1] );
	r.RotateZ( gRotation[0] );
	_normalTVec.Transform( r );
                                  
	_planeNormal[0] = _normalTVec[0];
	_planeNormal[1] = _normalTVec[1];
	_planeNormal[2] = _normalTVec[2];
      }
      catch(...) {
	streamlog_out(DEBUG5) << "no sensor rotation is given in the GEAR steering file, assume NONE" << std::endl;
      }
    }
    if( isFirstEvent() ) _isFirstEvent = false;

  } // FirstEvent
  
 
  // ################## TIMING CALCULATIONS ###########################
  static int64_t time_prev_tlu = event->getTimeStamp(); //init

  // Initialize variables:
  int64_t time_now_tlu = event->getTimeStamp();
  if( _nEvt == 1 ) time_event1 = time_now_tlu;

  // phase from TLU time between events:

  // TLU time stamp: count 384 MHz clock (384 = 48*8)
  // CMS timestamp: count 40 MHz clock, take as reference (HP pulser)

  // 31.10.2013, 80 MHz, 78*12.6 ns
  gTLU = 0.384678355; // [GHz]

  // fTLU = Frequency Trigger Logic Unit; Clock freq of the TLU box in Hz
  // gTLU = gigahertz Frequency Trigger Logic Unit; Freq in GHz
  //fTLU = 375.056 * 1.024 * 1E6; // [Hz] = 384.05734 MHz, from Umlauftakt, see below
  fTLU = gTLU * 1E9;

  double eventTime =  (time_now_tlu-time_event0)/fTLU;
  
  if(_have_skip_db){
    
    bool is_in_uptime = false;
    
    for(unsigned int i = 0; i < _uptimestart.size(); i++){
      if(eventTime > _uptimestart.at(i) && eventTime < _uptimeend.at(i) ){
	is_in_uptime = true;
	break;
      }
      
    }

    if(!is_in_uptime) return;
    
  }


  if( _nEvt > 0 ) {
    // Fill the dt plots (scanning for Umlauftakt etc.)
    FillDeltaTPlots(time_now_tlu,time_prev_tlu,time_event0);

    // Detection of long gaps between two telescope events:
    if( ( time_now_tlu - time_prev_tlu ) / fTLU > 10 ) {
      streamlog_out(WARNING5) << "long gap event "
			      << std::setw( 6) << std::setiosflags(std::ios::right)
			      << event->getEventNumber()
			      << " in run " << std::setw( 6) << std::setiosflags(std::ios::right)
			      << runNumber
			      << ", time " << std::setw(10) << std::setprecision(9)
			      << ( time_now_tlu - time_event0 ) /fTLU << " s"
			      << ", dt = " << int( (time_now_tlu - time_prev_tlu)/384E0 ) << " us"
			      << std::endl;
    }// Gap detection
  }//dt plots

  // Some more timing histograms:
  if(time_now_tlu > time_event0 ) {

    // Time of arrival TLU triggers:
    tlutrigvstusHisto->fill(int((time_now_tlu)/gTLU/1E6 )); // [ms]

  } // > time_event0
   
  // Informational message about event number and timing:
  if( _nEvt < 20 || _nEvt % 1000 == 0 ) {

    streamlog_out( MESSAGE2 ) << "Processing event "
			      << std::setw( 7) << std::setiosflags(std::ios::right)
			      << event->getEventNumber()
			      << " in run " << std::setw( 4) << std::setiosflags(std::ios::right)
			      << runNumber
			      << ", time " << " = " << std::fixed << std::setprecision(1)
			      << ( time_now_tlu - time_event0 )/fTLU << " s"
			      << ", dt = " << int( (time_now_tlu - time_prev_tlu)/gTLU/1E3 ) << " us"
			      << ", rate "
			      << (_nEvt+1) / ( time_now_tlu - time_event0 + 0.1 ) * fTLU
			      << " Hz"
			      << std::endl;
  }
  
  // Renew the timestamps for the next event:
  time_prev_tlu = time_now_tlu;
  
  // Increase event count
  _nEvt++;
     
  //----------------------------------------------------------------------------
  // check input collection (aligned hits):

#if 0
  {
    LCCollection* hits = event->getCollection("hit");
    LCCollection* prealignedHits = event->getCollection("prealignedhits");
    LCCollection* alignedHits = event->getCollection("alignedhits");

    std::cout << "Raw hits: " << std::endl;
    for(unsigned int i = 0; i < hits->getNumberOfElements(); ++i)
    {
      TrackerHit * meshit = dynamic_cast<TrackerHit*>( hits->getElementAt(i) );
      std::cout << i << ": " << meshit->getPosition()[0] << ", " << meshit->getPosition()[1] << ", " << meshit->getPosition()[2] << std::endl;
    }

    std::cout << "Prealigned hits: " << std::endl;
    for(unsigned int i = 0; i < prealignedHits->getNumberOfElements(); ++i)
    {
      TrackerHit * meshit = dynamic_cast<TrackerHit*>( prealignedHits->getElementAt(i) );
      std::cout << i << ": " << meshit->getPosition()[0] << ", " << meshit->getPosition()[1] << ", " << meshit->getPosition()[2] << std::endl;
    }

    std::cout << "Aligned hits: " << std::endl;
    for(unsigned int i = 0; i < alignedHits->getNumberOfElements(); ++i)
    {
      TrackerHit * meshit = dynamic_cast<TrackerHit*>( alignedHits->getElementAt(i) );
      std::cout << i << ": " << meshit->getPosition()[0] << ", " << meshit->getPosition()[1] << ", " << meshit->getPosition()[2] << std::endl;
    }
  }
#endif

  LCCollection *collection, *dutpx, *refpx;
  try {
    collection = event->getCollection( _inputCollectionTelescope );
    dutpx = event->getCollection( _inputPixelsDUT );
    refpx = event->getCollection( _inputPixelsREF );
  }
  catch( lcio::DataNotAvailableException& e) {
    streamlog_out( DEBUG1 ) << "Not able to get collections "
			    << _inputCollectionTelescope << " "
			    << _inputPixelsDUT << " "
			    << _inputPixelsREF
			    << "\nfrom event " << event->getEventNumber()
			    << " in run " << runNumber  << std::endl;

    return;
  }


  // DUT chip is hanging (0,0) in upper right corner looking downstream:
  bool hanging = 1;

  // Horizontal DUT alignment for FPIX test beam:
  bool FPIX = 0;

  // Rotatet chip on desytb-pcb-rot:
  bool rot90 = false;
  if(_DUT_chip == 506) rot90 = true;

  // ETH and KIT test beam Mar 2013 analog ROCs
  // Horizontal DUT alignment for ETH test beam:
  bool ETHh = 0;

  EUTelGenericSparsePixel* pixel = new EUTelGenericSparsePixel;

  // Read the DUT event:
  std::vector<CMSPixel::pixel> * dutPixels = new std::vector<CMSPixel::pixel>;
  LCCollectionVec * dutCollectionVec  = dynamic_cast < LCCollectionVec * > (dutpx);
  CellIDDecoder<TrackerDataImpl> cellDecoder( dutCollectionVec );

  // Only one DUT plane!
  if(dutCollectionVec->size() != 1) { throw StopProcessingException(this); }

  TrackerDataImpl * dutData = dynamic_cast< TrackerDataImpl * > ( dutCollectionVec->getElementAt(0) );
  unsigned int sensorID = static_cast<unsigned int > ( cellDecoder( dutData )["sensorID"] );
  streamlog_out ( DEBUG5 ) << "evt " << event->getEventNumber() << " SensorID " << sensorID << endl;

  auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> > pixelDUTData( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>( dutData ));
  streamlog_out ( DEBUG5 ) << "Processing data on detector " << sensorID << ", " << pixelDUTData->size() << " pixels " << endl;

  //This for-loop loads all the hits of the given event and detector plane and stores them
  for(size_t i = 0; i < pixelDUTData->size(); ++i ) {
    // Load the information of the hit pixel into genericPixel
    pixelDUTData->getSparsePixelAt(i, pixel);
      
    CMSPixel::pixel px;
    px.roc = 8;

    if(rot90) {
      // Deal with backwards-incompatible change in DAQ software:
      // before run 15268 the PCBTYPE flag has not been set and
      // col/row are not inverted:
      if(_nRun > 15268) {
	px.col = pixel->getYCoord();
	px.row = pixel->getXCoord();
      } else {
	px.col = pixel->getXCoord();
	px.row = pixel->getYCoord();
      }
    }
    else {
      px.col = pixel->getXCoord();
      px.row = pixel->getYCoord();
    }
    px.raw = pixel->getSignal();

    dutpxcolHisto->fill(px.col);
    dutpxrowHisto->fill(px.row);
    dutpxadcHisto->fill(px.raw);
    //and push this pixel back
    dutPixels->push_back(px);
  }
  //streamlog_out( WARNING ) << "Evt " << event->getEventNumber() << ": " << dutPixels->size() << " on DUT";

  // Read the REF event:
  std::vector<CMSPixel::pixel> * refPixels = new std::vector<CMSPixel::pixel>;
  LCCollectionVec * refCollectionVec  = dynamic_cast < LCCollectionVec * > (refpx);
  CellIDDecoder<TrackerDataImpl> cellDecoderRef( refCollectionVec );

  // Only one REF plane!
  if(refCollectionVec->size() != 1) { throw StopProcessingException(this); }

  TrackerDataImpl * refData = dynamic_cast< TrackerDataImpl * > ( refCollectionVec->getElementAt(0) );
  unsigned int sensorIDref = static_cast<unsigned int > ( cellDecoder( refData )["sensorID"] );
  streamlog_out ( DEBUG5 ) << "evt " << event->getEventNumber() << " SensorID " << sensorIDref << endl;

  auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> > pixelREFData( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>( refData ));
  streamlog_out ( DEBUG5 ) << "Processing data on detector " << sensorIDref << ", " << pixelREFData->size() << " pixels " << endl;

  //This for-loop loads all the hits of the given event and detector plane and stores them
  for(size_t i = 0; i < pixelREFData->size(); ++i ) {
    // Load the information of the hit pixel into genericPixel
    pixelREFData->getSparsePixelAt(i, pixel);
      
    CMSPixel::pixel px;
    px.roc = 8;
    px.col = pixel->getXCoord();
    px.row = pixel->getYCoord();
    px.raw = pixel->getSignal();

    //and push this pixel back
    refPixels->push_back(px);
  }
  //streamlog_out( WARNING ) << "Evt " << event->getEventNumber() << ": " << refPixels->size() << " on REF";


  //----------------------------------------------------------------------------
  // Copy hits to local table
  // Assign hits to sensor planes

  streamlog_out( DEBUG4 )  << "Total of " << collection->getNumberOfElements() << " tracker hits in input collection " << std::endl;

  nAllHitHisto->fill(collection->getNumberOfElements());

  //----------------------------------------------------------------------------

  std::vector<hit> * hits = new std::vector<hit>;

  // Extract all hits from the LCIO collection and calculate the
  // position uncertainty if necessary:
  for( int ihit = 0; ihit < collection->getNumberOfElements(); ihit++ ) {

    hit newhit;
    TrackerHit * meshit = dynamic_cast<TrackerHit*>( collection->getElementAt(ihit) );
    const double * pos = meshit->getPosition();
    const EVENT::FloatVec cov = meshit->getCovMatrix();

    // Write the position:
    newhit.x = pos[0];
    newhit.y = pos[1];
    newhit.z = pos[2];

    // Find Plane ID to which the hit belongs by minimizing its
    // distance in Z
    //FIXME to be fixed with more general geometry description!
    double distMin = 99; // [mm]

    bool foundplane = false;
    for( int ipl = 0; ipl < _nTelPlanes; ipl++ ) {
      if( abs(newhit.z - _planePosition[ipl]) < distMin ) {
	newhit.plane = ipl;
	foundplane = true;
	distMin = abs(newhit.z - _planePosition[ipl]);
      }
    }

    if(!foundplane) {
      streamlog_out(ERROR5) << "Could not associate hit with a telescope plane. Skipping event." << std::endl;
      return;
    }

    // ...and uncertainty:
    if( cov.at(0) > 0. ) newhit.ex = sqrt(cov.at(0));
    else newhit.ex = _planeResolution[newhit.plane];

    if( cov.at(2) > 0. ) newhit.ey = sqrt(cov.at(2));
    else newhit.ey = _planeResolution[newhit.plane];

    streamlog_out(DEBUG3) << "Hit " << ihit << ": " << newhit << std::endl;

    // Add it to the vector of telescope hits:
    hits->push_back(newhit);

  } // loop over hits

  streamlog_out(DEBUG4) << "Event " << event->getEventNumber()
			<< " contains " << hits->size() << " hits" << std::endl;



  // Do the event shifting if necessary:
  int nskip = max(_skip_dut,max(_skip_ref,_skip_tel));
  if(nskip > 0) {
    streamlog_out ( DEBUG4 ) << "Evt " << event->getEventNumber() << " read: DUT " << dutPixels->size() << "px, REF " << refPixels->size() << "px, TEL " << hits->size() << "px" << endl;

    // Do the skipping magic:
    if(event->getEventNumber() < nskip) {
      if(event->getEventNumber() < _skip_dut) dut_event_buffer.push_back(*dutPixels);
      if(event->getEventNumber() < _skip_ref) ref_event_buffer.push_back(*refPixels);
      if(event->getEventNumber() < _skip_tel) tel_event_buffer.push_back(*hits);
      throw SkipEventException(this);
    }
    streamlog_out ( DEBUG4 ) << "Evt " << event->getEventNumber() << " buffers: DUT " << dut_event_buffer.size() 
			     << "ev, REF " << ref_event_buffer.size() 
			     << "ev, TEL" << tel_event_buffer.size() << "ev" << endl;

    // Retrieve info back from buffer:

    // Continue pushing and reading the shifted data streams:
    if(!dut_event_buffer.empty()) {
      dut_event_buffer.push_back(*dutPixels);
      dutPixels->clear();
      std::vector<CMSPixel::pixel> temp  = dut_event_buffer.front();
      for(size_t px = 0; px < temp.size(); px++) { dutPixels->push_back(temp.at(px)); }
      dut_event_buffer.erase(dut_event_buffer.begin());
    }

    if(!ref_event_buffer.empty()) {
      ref_event_buffer.push_back(*refPixels);
      refPixels->clear();
      std::vector<CMSPixel::pixel> temp2  = ref_event_buffer.front();
      for(size_t px = 0; px < temp2.size(); px++) { refPixels->push_back(temp2.at(px)); }
      ref_event_buffer.erase(ref_event_buffer.begin());
    }

    if(!tel_event_buffer.empty()) {
      tel_event_buffer.push_back(*hits);
      hits->clear();
      std::vector<hit> temp3  = tel_event_buffer.front();
      for(size_t px = 0; px < temp3.size(); px++) { hits->push_back(temp3.at(px)); }
      tel_event_buffer.erase(tel_event_buffer.begin());
    }

    streamlog_out ( DEBUG4 ) << "Evt " << event->getEventNumber() << " buffers post-read: DUT " << dut_event_buffer.size() 
			     << "ev, REF " << ref_event_buffer.size()
			     << "ev, TEL" << tel_event_buffer.size() << "ev" << endl;

    streamlog_out ( DEBUG4 ) << "Evt " << event->getEventNumber() << " fetched: DUT " << dutPixels->size() << "px, REF " << refPixels->size() << "px, TEL " << hits->size() << "px" << endl;
  }

  // Calibrate the pixel hits with the initialized calibration data:
  if(!CalibratePixels(dutPixels,dut_calibration,_DUT_conversion))
    throw StopProcessingException(this);
  for(size_t i = 0; i < dutPixels->size(); i++) { dutpxqHisto->fill(dutPixels->at(i).vcal); }
  ClustDUT = GetClusters(dutPixels,_DUT_chip);

  if(!CalibratePixels(refPixels,ref_calibration,_REF_conversion))
    throw StopProcessingException(this);
  ClustREF = GetClusters(refPixels,_REF_chip);


  streamlog_out(DEBUG4) << std::setw(6) << std::setiosflags(std::ios::right) << event->getEventNumber() << ". DUT pixels " << dutPixels->size();

  if( dutPixels->size() > 0 ) {
    streamlog_out(DEBUG3) << ", clusters at";
    for( std::vector<cluster>::iterator c = ClustDUT.begin(); c != ClustDUT.end(); c++ ){
      streamlog_out(DEBUG3) << "  (" << c->col << ", " << c->row << ")";
    }
  }
  streamlog_out(DEBUG4) << std::endl;


  // Fill some informational plots about the REF and DUT cluster distributions:
  FillClusterStatisticsPlots( ClustDUT, dutPixels->size(), ClustREF, refPixels->size() );

  // Fill the telescope plane correlation plots:
  TelescopeCorrelationPlots(hits);


  // Compare first pixels of first and second cluster:
  if(ClustDUT.size() > 1) {
    cmspxq2c1stHisto->fill(ClustDUT.at(0).vpix.at(0).vcal);
    cmspxq2c2ndHisto->fill(ClustDUT.at(1).vpix.at(0).vcal);
  }


  //----------------------------------------------------------------------------
  // DUT alignment: relative to _planePosition[2]
  // REF alignment: relative to _planePosition[5]

  bool lDUT = 1; // DUT present
  DUTalignx = _DUTalignx;
  DUTaligny = _DUTaligny;

  DUTz = _DUTz + _planePosition[2]; // default

  DUTrot = _DUTrot;
  tilt = _DUTtilt;
  turn = _DUTturn;

  double cutx = _cutx;
  double cuty = _cuty;

  double slope_x = _slope_x;
  double slope_y = _slope_y;

  double wt = atan(1.0) / 45.0; // pi/180 deg
  double DUTX0 = 0.07; // init, for analog ROC with thick carrier socket

  if( runNumber >= 4685 ) { // since Aug 2012: digital ROCs on thin PCB
    DUTX0 = 0.01;
  }
  DUTX0 = DUTX0 / cos( tilt * wt ) / cos( turn*wt );

  double REFz = _REFz + _planePosition[5]; // timing reference plane
  double REFalignx =  _REFalignx;
  double REFaligny = _REFaligny;
  double REFrot = _REFrot;

  double pitchcol = 0.150; // [mm]
  double pitchrow = 0.100; // [mm]
  double pitchcolt = 0.150 * cos( turn * wt ); // [mm] in telescope system
  double pitchrowt = 0.100 * cos( tilt * wt ); // [mm] in telescope system

  if( FPIX || ETHh) {
    // FPIX rot at 90 deg:
    // col tilted
    // row turned
    pitchcolt = 0.150 * cos( tilt * wt ); // [mm] in telescope system
    pitchrowt = 0.100 * cos( turn * wt ); // [mm] in telescope system
  }

  // normal vector on DUT surface:
  // N = ( 0, 0, -1 ) on DUT
  // transform into tele system:
  // tilt alpha around x
  // turn omega around y

  double Nx =-sin( turn*wt )*cos( tilt*wt );
  double Ny = sin( tilt*wt );
  double Nz =-cos( turn*wt )*cos( tilt*wt );

  double co = cos( turn*wt );
  double so = sin( turn*wt );
  double ca = cos( tilt*wt );
  double sa = sin( tilt*wt );
  double cf = cos( DUTrot );
  double sf = sin( DUTrot );

  double dNxdo =-co*ca;
  double dNxda = so*sa;
  double dNydo = 0;
  double dNyda = ca;
  double dNzdo = so*ca;
  double dNzda = co*sa;

  double norm = fabs(cos( turn*wt )) * fabs(cos( tilt*wt )); // length of Nz. Landau 21.6
  //double path = sqrt( 1 + tan(turn*wt)*tan(turn*wt) + tan(tilt*wt)*tan(tilt*wt) );
  //norm = 1/path; // test. Landau 22.4


  // -----------------------------------------------------------------
  // Downstream Telescope Triplets ("driplets")
  int ndriref = 0;
  int ndrilk = 0;

  // Generate new triplet set for the Telescope Downstream Arm:
  std::vector<triplet> * downstream_triplets = new std::vector<triplet>;
  downstream_triplets = FindTriplets(hits, 3, 4, 5);

  // Iterate over all found downstream triplets to fill histograms and match them to the REF and DUT:
  for( std::vector<triplet>::iterator drip = downstream_triplets->begin(); drip != downstream_triplets->end(); drip++ ){

    // Fill some histograms for downstream triplets:
    dridxHisto->fill( (*drip).getdx(4)*1E3 ); //sig = 7 um at 5 GeV
    dridxvsx->fill( (*drip).base().x, (*drip).getdx(4)*1E3 ); // check for rot
    dridxvsy->fill( (*drip).base().y, (*drip).getdx(4)*1E3 );
    dridxvstx->fill( (*drip).slope().x*1E3, (*drip).getdx(4)*1E3 ); // check for z shift
    dridxvsty->fill( (*drip).slope().y*1E3, (*drip).getdx(4)*1E3 );

    dridyHisto->fill( (*drip).getdy(4)*1E3 );
    dridyvsx->fill( (*drip).base().x, (*drip).getdy(4)*1E3 );
    dridyvsy->fill( (*drip).base().y, (*drip).getdy(4)*1E3 );
    dridyvstx->fill( (*drip).slope().x*1E3, (*drip).getdy(4)*1E3 );
    dridyvsty->fill( (*drip).slope().y*1E3, (*drip).getdy(4)*1E3 );

    drixHisto->fill( -(*drip).gethit(4).x ); // -x = x_DP = out
    driyHisto->fill( -(*drip).gethit(4).y ); // -y = y_DP =  up
    drixyHisto->fill( -(*drip).gethit(4).x, -(*drip).gethit(4).y );
    dritxHisto->fill( (*drip).slope().x*1E3 );
    drityHisto->fill( (*drip).slope().y*1E3 );


    // ################ HEAVY TRANSFORMATIONS DOWN HERE ##############################
    // CHANGE TO TRUE 3D ROTATION

    // extrapolate to REF timing plane: [mm]

    // telescope coordinates (right handed, look along beam):
    //   x_tele = inward (towards DESY)
    //   y_tele = down
    //   z = beam
    // REF coordinates (mounted upside down in old black socket carrier):
    //   col = -x = outward
    //   row =  y = down
    // DP coordinates (right handed, look at beam):
    //   x_DP = outward (away from DESY)
    //   y_DP = up
    //   z = beam
    //   x_DP = -x_tele
    //   y_DP = -y_tele, together = 180 deg rot

    // Get impact point on the REF plane:
    double xB = (*drip).getx_at(REFz);
    double yB = (*drip).gety_at(REFz);

    // transform into REF system:

    double xBm = -xB + REFalignx; // invert x, shift to mid
    double yBm =  yB + REFaligny; // shift to mid
    double xBr = xBm - DUTrot*yBm; // rotate CMS pixel
    double yBr = yBm + DUTrot*xBm;
    //double xBt = xBr + 0.075; // shift by half bin?
    double xBt = xBr;
    double yBt = yBr; // no tilt for REF

    double xmod = fmod( 20 + xBt, 0.3 ) * 1E3; // [0,300] um, 2 pixel wide
    double ymod = fmod( 20 + yBt, 0.2 ) * 1E3; // [0,200] um

    drixrefHisto->fill( -xB ); // -xB = x_DP = out
    driyrefHisto->fill( -yB ); // -yB = y_DP = up
    drixyrefHisto->fill( -xB, -yB );

    // Correlate Downstream Telescope Triplet with the REF:
    if( refPixels->size() > 0 ) {

      ndriref++;

      for( std::vector<cluster>::iterator c = ClustREF.begin(); c != ClustREF.end(); c++ ){

	refxxHisto->fill( c->col, xB ); // anti-correlation
	refyyHisto->fill( c->row, yB ); // correlation

	double refx = ( c->col - 26 ) * 0.15; // -3.9..3.9 mm
	double refy = ( c->row - 40 ) * 0.10; // -4..4 mm

	refsxaHisto->fill( refx + xB );
	refdxaHisto->fill( refx - xB );
	refsyaHisto->fill( refy + yB );
	refdyaHisto->fill( refy - yB );

	double refsx = refx + REFrot*refy - REFalignx + xB; // residual x
	double refdy = refy - REFrot*refx - REFaligny - yB; // residual y

	if( _REF_chip >= 200 ) { // dig = REF
	  refsx = refx + REFrot*refy - REFalignx - xB; // residual x
	  refdy = refy - REFrot*refx - REFaligny + yB; // residual y
	}

	refsxHisto->fill( refsx );
	refdyHisto->fill( refdy );

	if( abs( refdy ) < 1 ) refsxcHisto->fill( refsx );
	if( abs( refsx ) < 1 ) refdycHisto->fill( refdy );

	if( abs( refdy ) < 1 && abs( refsx ) < 1 ) {

	  drixlkHisto->fill( -xB ); // -xB = x_DP = out
	  driylkHisto->fill( -yB ); // -yB = y_DP = up
	  drixylkHisto->fill( -xB, -yB );
	  refpixvsxmym->fill( xmod, ymod ); // occupancy map
	  refqHisto->fill( c->charge );
	  refqvsxmym->fill( xmod, ymod, c->charge ); // cluster charge profile
	  refdyvsx->fill( refx, refdy*1E3 ); // dy vs x: rot?
	  refdyvsy->fill( refy, refdy*1E3 ); // dy vs y: tilt?
	  refdyvsty->fill( (*drip).slope().y*1E3, refdy*1E3 ); // dy vs tety: z shift?
	}

	// Checking for a match:
	if( abs( refdy ) < 0.3 && abs( refsx ) < 0.3 ) {
	  n_matched_clusters_ref++;
	  //FIXME link to REF instead of DUT?
	  (*drip).linked_dut = true;
	  reflkcolHisto->fill( c->col );
	  reflkrowHisto->fill( c->row );
	  ndrilk++;
	}
      }//loop over CMS clusters
    }// REF Correlation


     // Correlate Downstream Telescope Triplet with the  DUT:
    if( dutPixels->size() > 0 ) {

      //extrapolate to DUT plane: [mm]

      // Get downstream triplet impact point on DUT plane:
      double xs = (*drip).getx_at(DUTz);
      double ys = (*drip).gety_at(DUTz);

      for( std::vector<cluster>::iterator c = ClustDUT.begin(); c != ClustDUT.end(); c++ ){

	double bacx = ( c->col - 26 ) * pitchcolt; // -3.9..3.9 mm
	double bacy = ( c->row - 40 ) * pitchrowt; // -4..4 mm

	  if (ETHh){
	    //bacx = ( c->row - 40 ) * pitchrowt; // -4..4 mm
	    bacx = ( 40 - c->row ) * pitchrowt; // -4..4 mm
	    bacy = ( c->col - 26 ) * pitchcolt; // -3.9..3.9 mm
	  }

	bacsxaHisto->fill( bacx + xs );
	bacdyaHisto->fill( bacy - ys );

	double bacsx = bacx + DUTrot*bacy - DUTalignx + xs; // residual x
	double bacdy = bacy - DUTrot*bacx - DUTaligny - ys; // residual y

	if( _DUT_chip >= 100 ) { // xdb and dig mounted upright
	  bacsx = bacx + DUTrot*bacy - DUTalignx - xs; // residual x
	  bacdy = bacy - DUTrot*bacx - DUTaligny + ys; // residual y
	}

	if( abs( bacdy ) < 1 ) bacsxcHisto->fill( bacsx*1E3 );
	if( abs( bacsx ) < 1 ) bacdycHisto->fill( bacdy*1E3 );

	if( c->charge * norm > 18 ) {
	  if( abs( bacdy ) < 1 ) bacsxcqHisto->fill( bacsx*1E3 );
	  if( abs( bacsx ) < 1 ) bacdycqHisto->fill( bacdy*1E3 ); // 18 um @ 5 GeV
	}

      }//loop over BAC clusters

    }//have telescope driplet

  }//loop over driplets

  ndriHisto->fill( downstream_triplets->size() );
  if( dutPixels->size() > 0 ) ndrirefHisto->fill( ndriref );
  lkBvst->fill( (time_now_tlu-time_event0)/fTLU, ndrilk );


  //##########################  0-1-2 upstream triplets #######################
  // here again a telescope hit triplet is formed, from planes 0 and 2
  // and the correlation to plane 1. the found triplets are
  // extrapolated and matched to the DUT plane and exhaustive
  // histograms are written.
  // Furthermore the hits of a triplet matched to the DUT are
  // collected and fed into GBL for a track fit.

  int ntrilk = 0;

  // Generate new triplet set for the Telescope Downstream Arm:
  std::vector<triplet> * upstream_triplets = new std::vector<triplet>;
  upstream_triplets = FindTriplets(hits, 0, 1, 2);

  // Iterate over all found upstream triplets to fill histograms and match them to the REF and DUT:
  for( std::vector<triplet>::iterator trip = upstream_triplets->begin(); trip != upstream_triplets->end(); trip++ ) {

    // Fill some histograms for the upstream triplets:
    tridxHisto->fill( (*trip).getdx(1)*1E3 );
    tridyHisto->fill( (*trip).getdy(1)*1E3 );
    tridx1Histo->fill( (*trip).getdx(1)*1E0 );
    tridy1Histo->fill( (*trip).getdy(1)*1E0 );
    tridxvsx->fill( (*trip).base().x, (*trip).getdx(1)*1E3 ); // check for rot
    tridxvsy->fill( (*trip).base().y, (*trip).getdx(1)*1E3 );
    tridxvstx->fill( (*trip).slope().x*1E3, (*trip).getdx(1)*1E3 ); // check for z shift
    tridxvsty->fill( (*trip).slope().y*1E3, (*trip).getdx(1)*1E3 );
    tridyvsx->fill( (*trip).base().x, (*trip).getdy(1)*1E3 );
    tridyvsy->fill( (*trip).base().y, (*trip).getdy(1)*1E3 );
    tridyvstx->fill( (*trip).slope().x*1E3, (*trip).getdy(1)*1E3 );
    tridyvsty->fill( (*trip).slope().y*1E3, (*trip).getdy(1)*1E3 );
    trixHisto->fill( -(*trip).gethit(1).x );
    triyHisto->fill( -(*trip).gethit(1).y );
    trixyHisto->fill( -(*trip).gethit(1).x, -(*trip).gethit(1).y );
    tritxHisto->fill( (*trip).slope().x*1E3 );
    trityHisto->fill( (*trip).slope().y*1E3 );

    // scale hit to CMS pixel:
    // sizeX="21.2"  sizeY="10.6"
    // npixelX="1152"  npixelY="576" 

    streamlog_out(DEBUG2) << " x " << int( ( (*trip).gethit(2).x + 10.6 ) / 0.15 )
			  << ", y " << int( ( (*trip).gethit(2).y +  5.3 ) / 0.10 ) << std::endl;

    // ######################## HEAVY TRANSFORMATIONS DOWN HERE ##############################
    //extrapolate to CMS: [mm]

    // telescope coordinates (right handed, look along beam):
    //   x_tele = inward (towards DESY)
    //   y_tele = down
    //   z = beam
    // REF coordinates (mounted upside down in old black socket carrier):
    //   col = -x = outward
    //   row =  y = down
    // DP coordinates (right handed, look at beam):
    //   x_DP = outward (away from DESY)
    //   y_DP = up
    //   z = beam
    //   x_DP = -x_tele
    //   y_DP = -y_tele, together = 180 deg rot around z

    double xA = (*trip).getx_at(DUTz); // triplet impact point on CMS
    double yA = (*trip).gety_at(DUTz);
    double avx = (*trip).base().x;
    double avy = (*trip).base().y;
    double avz = (*trip).base().z;
    double tx = (*trip).slope().x;
    double ty = (*trip).slope().y;
    double zA = DUTz - avz; // z CMS from mid of triplet

    trixdutHisto->fill( -xA ); // -xA = x_DP = out
    triydutHisto->fill( -yA ); // -yA = y_DP = up
    trixydutHisto->fill( -xA, -yA );

    bool isolatedTrip = true;

    double ddAMin = -1.0;
    for( std::vector<triplet>::iterator tripIsoCheck = upstream_triplets->begin(); tripIsoCheck != upstream_triplets->end(); tripIsoCheck++ ) {
      if(trip != tripIsoCheck){
	double xAIsoCheck = (*tripIsoCheck).getx_at(DUTz);
	double xyIsoCheck = (*tripIsoCheck).gety_at(DUTz);
	double ddA = sqrt( fabs(xAIsoCheck - xA)*fabs(xAIsoCheck - xA) + fabs(xyIsoCheck - yA)*fabs(xyIsoCheck - yA) );
	if(ddAMin < 0 || ddA < ddAMin)
	  ddAMin = ddA;
      }
    }
    

    triddaMindutHisto->fill(ddAMin);
    if(ddAMin < 0.3)
      isolatedTrip = false;

    // intersect inclined track with tilted DUT plane:

    double zc = (Nz*zA - Ny*avy - Nx*avx) / (Nx*tx + Ny*ty + Nz); // from avz
    double yc = avy + ty * zc;
    double xc = avx + tx * zc;

    double dxcdz = tx*Nz / (Nx*tx + Ny*ty + Nz);
    double dycdz = ty*Nz / (Nx*tx + Ny*ty + Nz);

    double dzc = zc + avz - DUTz; // from DUT z0 [-8,8] mm

    double dzcda = ( dNzda*zA - dNyda*avy - dNxda*avx ) / (Nx*tx + Ny*ty + Nz) - zc*( dNxda*tx + dNyda*ty + dNzda ) / (Nx*tx + Ny*ty + Nz);
    double dzcdo = ( dNzdo*zA - dNydo*avy - dNxdo*avx ) / (Nx*tx + Ny*ty + Nz) - zc*( dNxdo*tx + dNydo*ty + dNzdo ) / (Nx*tx + Ny*ty + Nz);

    double ddzcdz = Nz / (Nx*tx + Ny*ty + Nz) - 1; // ddzc/dDUTz

    dzcvsxy->fill( xc, yc, dzc ); // tilted plane

    // transform into DUT system: (passive).
    // large rotations don't commute: careful with order

    double x1 = co*xc - so*dzc; // turn
    double y1 = yc;
    double z1 = so*xc + co*dzc;

    double x2 = x1;
    double y2 = ca*y1 + sa*z1; // tilt
    double z2 =-sa*y1 + ca*z1; // should be zero (in DUT plane). is zero

    double x3 = cf*x2 + sf*y2; // rot
    double y3 =-sf*x2 + cf*y2;
    double z3 = z2; // should be zero (in DUT plane). is zero

    z3vsxy->fill( xc, yc, z3 ); // should be zero. is zero

    double upsign = -1; // analog psi46: upside down
    if( _DUT_chip >= 100 )
      upsign = 1; // xdb and dig mounted upright

    double x4 = upsign*x3 + DUTalignx; // shift to mid
    if(ETHh) x4 = upsign*x3 - DUTalignx; // shift to mid

    double y4 =-upsign*y3 + DUTaligny; // invert y, shift to mid

    double xAt = x4;
    double yAt = y4;


    // ######################################################################################

    // reduce to 2x2 pixel region:

    double xmod = fmod( 9.075 + xAt, 0.3 ) * 1E3; // [0,300] um, 2 pixel wide
    double ymod = fmod( 9.050 + yAt, 0.2 ) * 1E3; // [0,200] um
    double ymd3 = fmod( 9.050 + yAt, 0.3 ) * 1E3; // [0,300] um
    double ymd6 = fmod( 9.050 + yAt, 0.6 ) * 1E3; // [0,600] um
    if(hanging) {
      xmod = fmod( 9.075 - xAt, 0.3 ) * 1E3; // [0,300] um, 2 pixel wide
      ymod = fmod( 9.050 - yAt, 0.2 ) * 1E3; // [0,200] um
      ymd3 = fmod( 9.050 - yAt, 0.3 ) * 1E3; // [0,300] um
      ymd6 = fmod( 9.050 - yAt, 0.6 ) * 1E3; // [0,600] um
    }
    if( FPIX || ETHh || rot90) { // x = col = yt, y = row = xt
      xmod = fmod( 9.075 - yAt, 0.3 ) * 1E3; // [0,300] um, 2 pixel wide
      ymod = fmod( 9.050 + xAt, 0.2 ) * 1E3; // [0,200] um
      ymd3 = fmod( 9.050 - xAt, 0.3 ) * 1E3; // [0,300] um
      ymd6 = fmod( 9.050 - xAt, 0.6 ) * 1E3; // [0,600] um
    }

    // 0 deg:
    bool ldot = 1; // bias dot, from cmsqvsxmym
    // Not rotated, cut in xmod (and ymod at small tilt)
    if(!rot90) {
      if( xmod < 105 ) ldot = 0; // dot at x = 125
      if( xmod > 195 ) ldot = 0; // and at x = 175
      if( tilt < 6 ) {
	if( ymod <  55 ) ldot = 0; // dot at y =  75
	if( ymod > 195 ) ldot = 0; // dot at y = 175
	if( ymod >  95 && ymod < 155 ) ldot = 0; // band between dots
      }
    }
    // Not rotated, cut in ymod (and xmod at small tilt)
    else {
      if( ymod > 195 ) ldot = 0; // dot at y = 175
      if( ymod <  55 ) ldot = 0; // dot at y =  75
      if( ymod >  95 && ymod < 155 ) ldot = 0; // band between dots
      if( tilt < 6 ) {
	if( xmod < 105 ) ldot = 0; // dot at y =  75
	if( xmod > 195 ) ldot = 0; // dot at y = 175
      }
    }

    bool lcore = 1; // pixel core, 2x2 pixel region
    if( xmod <  20 ) lcore = 0; // outer edge, see cmsncolvsxm
    if( xmod > 280 ) lcore = 0; // outer edge
    if( ymod <  20 ) lcore = 0; // outer edge, see cmsnrowvsym
    if( ymod > 180 ) lcore = 0; // outer edge
    if( xmod > 130 && xmod < 170 ) lcore = 0; // inner edge
    if( ymod >  80 && ymod < 120 ) lcore = 0; // inner edge


    // Extrapolate Upstream triplet to Downstream planes 3,4,5: Resolution studies
    for( std::vector<hit>::iterator lhit = hits->begin(); lhit != hits->end(); lhit++ ) {

      if( (*lhit).plane <= 2 ) continue; // want 3,4, or 5

      // Fill residuals of triplet and hit in the selected plane:
      if( (*lhit).plane == 3 ) {
	tridx3Histo->fill( (*trip).getdx((*lhit))*1E0 );
	tridy3Histo->fill( (*trip).getdy((*lhit))*1E0 ); // 65 um at 4.7 GeV with CMS
	tridx3bHisto->fill( (*trip).getdx((*lhit))*1E3 ); // finer binning
	tridy3bHisto->fill( (*trip).getdy((*lhit))*1E3 ); // 
      }
      else if( (*lhit).plane == 4 ) {
	tridx4Histo->fill( (*trip).getdx((*lhit))*1E0 );
	tridy4Histo->fill( (*trip).getdy((*lhit))*1E0 ); //174 um at 4.7 GeV
	tridx4bHisto->fill( (*trip).getdx((*lhit))*1E3 ); // finer binning
	tridy4bHisto->fill( (*trip).getdy((*lhit))*1E3 ); // 
      }
      else if( (*lhit).plane == 5 ) {
	tridx5Histo->fill( (*trip).getdx((*lhit))*1E0 );
	tridy5Histo->fill( (*trip).getdy((*lhit))*1E0 ); //273 um at 4.7 GeV
	tridx5bHisto->fill( (*trip).getdx((*lhit))*1E3 ); // finer binning
	tridy5bHisto->fill( (*trip).getdy((*lhit))*1E3 ); // 
      }
    }// Resolution studies


     //----------------------------------------------------------------------
     // correlate with CMS DUT:
    if( dutPixels->size() > 0 ) {

      cmstimingcut->fill(eventTime);

      // CMS pixel clusters:
      bool trackHasLostSeedPixel = false;
      int nLinkedClusters = 0;
      for( std::vector<cluster>::iterator c = ClustDUT.begin(); c != ClustDUT.end(); c++ ){

	if(ETHh || FPIX || rot90) {
	  // for ETH/KIT orientation
	  cmsxxHisto->fill( c->col, yA ); // anti-correlation x-x
	  cmsyyHisto->fill( c->row, xA ); // correlation y-y
	}
	else {
	  cmsxxHisto->fill( c->col, xA ); // anti-correlation x-x
	  cmsyyHisto->fill( c->row, yA ); // correlation y-y
	}

	// pix in clus:
	int colmin = 99;
	int colmax = -1;
	int rowmin = 99;
	int rowmax = -1;

	double qcol[52];
	for( int icol = 0; icol < 52; ++icol ) qcol[icol] = 0;

	double qrow[80];
	for( int irow = 0; irow < 80; ++irow ) qrow[irow] = 0;

	for( std::vector<CMSPixel::pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); px++ ){

	  if(px->vcal < 0) continue;
	  
	  qcol[px->col] += px->vcal; // project cluster onto cols
	  qrow[px->row] += px->vcal; // project cluster onto rows
	  if( px->col < colmin ) { colmin = px->col; }
	  if( px->col > colmax ) { colmax = px->col; }
	  if( px->row < rowmin ) { rowmin = px->row; }
	  if( px->row > rowmax ) { rowmax = px->row; }

	}//pix
 
	bool fiducial = 1;
	if( rowmin ==  0 ) fiducial = 0;
	if( rowmax == 79 ) fiducial = 0;
	if( colmin ==  0 ) fiducial = 0;
	if( colmax == 51 ) fiducial = 0;

	bool lowClusterCharge = false;
	if(c->charge < 8) { lowClusterCharge = true; }

	
	int ncol = colmax - colmin + 1;
	int nrow = rowmax - rowmin + 1;

	double a1 = 0; //1st
	double a2 = 0; //2nd
	int i1 = 99;
	int i2 = 99;

	// eta-algo in rows:

	for( int irow = 0; irow < 80; ++irow ) {
	  if( qrow[irow] > a1 ) {
	    a2 = a1;
	    a1 = qrow[irow];
	    i2 = i1;
	    i1 = irow;
	  }
	  else if( qrow[irow] > a2 ) {
	    a2 = qrow[irow];
	    i2 = irow;
	  }
	}
	double a12 = a1 + a2;
	double eta = 0;
	if( a12 > 1 ) eta = ( a1 - a2 ) / a12;
	if( i1 > i2 ) eta = -eta;

      
	if(ncol > 1) {
	  // First pixel, second pixel?
	  cmspxq1stHisto->fill(c->vpix.at(0).vcal);
	  cmspxq2ndHisto->fill(c->vpix.at(1).vcal);
	}
	if(ncol > 2) {
	  cmspxq3rdHisto->fill(c->vpix.at(2).vcal);
	}

	// skew calculation (3rd moment):
	double skw = 0;
	double skw_unnorm = 0;
	double tst = 0;
	if(rot90) {
	  // Sum third powers of x-x_cog:
	  for(int col = colmin; col <= colmax; col++) { skw += pow((col - c->col),3)*qcol[col]; tst += qcol[col];}
	  skw_unnorm = 8*skw/pow(ncol,3);
	  // Normalize to total charge and cluster length/2 ^3:
	  skw /= (c->charge*pow(ncol,3));
	  skw *= 8;
	}
	else {
	  // Sum third powers of x-x_cog:
	  for(int row = rowmin; row <= rowmax; row++) { skw += pow((row - c->row),3)*qrow[row]; }
	  skw_unnorm = 8*skw/pow(nrow,3);
	  // Normalize to total charge and cluster length/2 ^3:
	  skw /= (c->charge*pow(nrow,3));
	  skw *= 8;
	}
	cmsskwHisto->fill(skw);

	// Detect huge skew values:
	if(fabs(skw) > 0.2) {
	  cmsskw1colcogHisto->fill(c->col);
	  cmsskw1rowcogHisto->fill(c->row);
	  cmsskw1qHisto->fill(c->charge);
	  cmsskw1ncolHisto->fill(ncol);
	  cmsskw1nrowHisto->fill(nrow);
	}
	else {
	  cmsskw0qHisto->fill(c->charge);
	  cmsskw0ncolHisto->fill(ncol);
	  cmsskw0nrowHisto->fill(nrow);
	}

	// lq: Cut on cluster charge, checking whether lies inside the Landau peak
	bool lq = 0;
	double Q0 = c->charge * norm; // cluster charge normalized to vertical incidence
	//if( Q0 > 18 &&  Q0 < 35 ) lq = 1;
	if( Q0 > 17 &&  Q0 < 30 ) lq = 1;

	if(lq) cmslq0Histo->fill(Q0);

	// DUT - triplet:
	// Move from chip coordinates (col, row) to physical
	// coordinates: cmsx, cmsy
	double cmsx = ( c->col - 26.0 ) * pitchcol; // -3.9..3.9 mm
	double cmsy = ( c->row - 40.0 ) * pitchrow; // -4..4 mm

	if(rot90) {
	  cmsx = ( c->row - 40.0 ) * pitchrow; // -4..4 mm
	  cmsy = ( 26.0 - c->col ) * pitchcol; // -3.9..3.9 mm
	}
	else if(hanging) {
	  cmsx = ( 26.0 - c->col ) * pitchcol; // -3.9..3.9 mm
	  cmsy = ( 40.0 - c->row ) * pitchrow; // -4..4 mm
	}
	else if( FPIX ) {
	  // FPIX rot at 90 deg:
	  // col = y
	  // row = x
	  cmsx = ( c->row - 40 ) * pitchrowt; // old style
	  cmsy = ( c->col - 26 ) * pitchcolt; // -3.9..3.9 mm
	}
	else if (ETHh) {
	  cmsx = ( 40 - c->row ) * pitchrowt; // KIT orientation
	  cmsy = ( c->col - 26 ) * pitchcolt; // -3.9..3.9 mm
	}

	cmssxaHisto->fill( cmsx + x3 ); // rot, tilt and turn but no shift
	cmsdxaHisto->fill( cmsx - x3 );

	cmssyaHisto->fill( cmsy + y3 );
	cmsdyaHisto->fill( cmsy - y3 );

	double dx4 = cmsx - x4;
	double dy4 = cmsy - y4;

	cmsdx4Histo->fill( dx4 );
	cmsdy4Histo->fill( dy4 );

	double dy5 = dy4;
	// Skew correction, limit to multi-pixel clusters:
	if(_have_skew_db) {
	  if((rot90 && ncol > 1) || (!rot90 && nrow > 1)) {
	    dy5 += (skew_par0 + skew_par1*skw)*1E-3;
	    cmsskwcorrHisto->fill((skew_par0 + skew_par1*skw));
	  }
	}

	cmsdx5Histo->fill( dx4 );
	cmsdy5Histo->fill( dy5 );

	double cmsdx = dx4;
	double cmsdy = dy5;
	// Comparison: new skew correction:
	double cmsdy0 = dy4;

	// Skew only for fiducial clusters, avoid edge effects:
	if(fiducial && 
	   fabs( cmsdx ) < cutx && 
	   isolatedTrip &&
	   fabs( ty-0.000 ) < slope_y &&
	   fabs( tx-0.000 ) < slope_x) {
	  cmsskwfctHisto->fill(skw);
	  cmsskwufctHisto->fill(skw_unnorm);
	  cmsskwuvsq->fill(Q0,skw_unnorm);
	  cmsskwvsq->fill(Q0,skw);
	  if(rot90) cmsskwvsqn->fill(Q0/ncol,skw);
	  else cmsskwvsqn->fill(Q0/nrow,skw);
	  cmsqvsskwu->fill(skw_unnorm,Q0);
	  cmsqvsskw->fill(skw,Q0);
	}

	// Even/odd effect:
	int iodd = static_cast<int>(floor( fmod( c->col, 2 ) ));
	if(iodd) { cmsdyoddHisto->fill(cmsdy*1E3); }
	else { cmsdyevenHisto->fill(cmsdy*1E3); }

	double oddcorr = 1.0E-3;
	if(_DUT_chip >= 500) { oddcorr = 1.0E-3; }
	else if( _DUT_chip >= 200 ) { // even/odd col effect for dig
	  oddcorr = 1.5E-3;
	}
	oddcorr = 0;

	if(!rot90) {
	  if( iodd ) {// odd
	    cmsdy  = cmsdy  - oddcorr;
	    cmsdy0 = cmsdy0 - oddcorr;
	  }
	  else {
	    cmsdy  = cmsdy  + oddcorr;
	    cmsdy0 = cmsdy0 + oddcorr;
	  }
	}
	// For rotated: other way around:
	else {
	  if( iodd ) {// odd
	    cmsdy  = cmsdy  + oddcorr;
	    cmsdy0 = cmsdy0 + oddcorr;
	  }
	  else {
	    cmsdy  = cmsdy  - oddcorr;
	    cmsdy0 = cmsdy0 - oddcorr;
	  }
	}

	bool seedPixelLost = false;     //Check if the seedPixel was
					//probably lost
	if(fiducial && fabs( cmsdx ) < cutx && fabs( ty-0.000 ) < 0.002 &&  fabs( tx-0.000 ) < 0.002 ) { //Same
													  //requirements 
													  //as for cmsdyfct
	  if(c->size == 1 && fabs( cmsdy ) > 0.04){
	    seedPixelLost = true;
	    trackHasLostSeedPixel = true;
	  }

	}

	cmsdxHisto->fill( cmsdx*1E3 );
	cmsdyHisto->fill( cmsdy*1E3 );
	cmsdy0Histo->fill( cmsdy0*1E3 );
	
	if( fiducial ) {

	  cmsdxfHisto->fill( cmsdx*1E3 );
	  cmsdyfHisto->fill( cmsdy*1E3 );

	  // Check for difference in 2col clusters:
	  if(ncol == 2) {
	    // Cluster in one double column?
	    if(colmin%2 == 0) { 
	      cmsdyfdc1Histo->fill(cmsdy*1E3);
	      cmsdyfdc1qHisto->fill(Q0);
	    }
	    // Cluster spread over two double columns
	    else {  
	      cmsdyfdc2qHisto->fill(Q0); 
	    }
	  }

	  if( fabs( cmsdy ) < cuty ) cmsdxfcHisto->fill( cmsdx*1E3 );
	  if( fabs( cmsdx ) < cutx ) cmsdyfcHisto->fill( cmsdy*1E3 );
	  
	  
	}//CMS fiducial

	 // accumulate cuts for y:
	if( fiducial && fabs( cmsdx ) < cutx && isolatedTrip) {

	  if(      nrow == 1 )
	    cmsdyfc1Histo->fill( cmsdy*1E3 ); // 3972: 7.7
	  else if( nrow == 2 )
	    cmsdyfc2Histo->fill( cmsdy*1E3 ); // 3972: 9.5
	  else
	    cmsdyfc3Histo->fill( cmsdy*1E3 ); // 3872: 68

	  if(      ncol == 1 )
	    cmsdyfc1cHisto->fill( cmsdy*1E3 ); // 3972: 7.7
	  else if( ncol == 2 )
	    cmsdyfc2cHisto->fill( cmsdy*1E3 ); // 3972: 9.5
	  else
	    cmsdyfc3cHisto->fill( cmsdy*1E3 ); // 3872: 68

	  if(      Q0 < 18 )
	    cmsdyq0Histo->fill( cmsdy*1E3 );

	  else if( Q0 < 35 ){
	    cmsdyq1Histo->fill( cmsdy*1E3 );
	    if( eta < 0 ) {
	      cmsdyeta0Histo->fill( cmsdy*1E3 );
	    }
	    else {
	      cmsdyeta1Histo->fill( cmsdy*1E3 );
	    }
	  }
	  else
	    cmsdyq2Histo->fill( cmsdy*1E3 );

	  if( fabs( ty-0.000 ) < slope_y &&  fabs( tx-0.000 ) < slope_x ) {
	    cmsdyfctHisto->fill( cmsdy*1E3 );
	    if(lowClusterCharge)
	      cmsdyfctLowChargeHisto->fill( cmsdy*1E3 );
	    if(!lowClusterCharge)
	      cmsdyfctHighChargeHisto->fill( cmsdy*1E3 );
	    if(c->size == 1)
	      cmsdyfctOnePixelHisto->fill( cmsdy*1E3 );
	    if(c->size == 1 && !lowClusterCharge)
	      cmsdyfctOnePixelHighChargeHisto->fill( cmsdy*1E3 );
	    if(c->size == 1 && lowClusterCharge)
	      cmsdyfctOnePixelLowChargeHisto->fill( cmsdy*1E3 );
	    if( nrow <= 2 ) cmsdyfcntHisto->fill( cmsdy*1E3 );
	  }

	  if( fabs( ty-0.000 ) < slope_y &&
	      fabs( tx-0.000 ) < slope_x &&
	      Q0 > 16 ) {

	    cmsdyfctqHisto->fill( cmsdy*1E3 ); // 7.8 um @ 4 GeV, 19 deg
	    if( nrow <= 2 ) cmsdyfcntqHisto->fill( cmsdy*1E3 );

	    if( Q0 < 35 ) {
	      cmsdyfctq1Histo->fill( cmsdy*1E3 ); // 7.8 um @ 4 GeV, 19 deg, more Gaussian
	      if( nrow <= 2 ) cmsdyfcntq1Histo->fill( cmsdy*1E3 );
	      if( c->col < 26 )
		cmsdyfctq1lHisto->fill( cmsdy*1E3 ); // xdb
	      else
		cmsdyfctq1rHisto->fill( cmsdy*1E3 ); // xdb
	    }
	    if( Q0 < 30 ) {
	      cmsdyfctq2Histo->fill( cmsdy*1E3 ); // inserted 26.12.2012
	    }
	    if( Q0 < 25 ) {
	      cmsdy0fctq3Histo->fill( cmsdy0*1E3 ); // for comparison, without skew correction
	      cmsdyfctq3Histo->fill( cmsdy*1E3 ); // was fctq2. 7.4 um @ 4 GeV, 19 deg
	      if( ldot ) 
		cmsdyfctqdotHisto->fill( cmsdy*1E3 ); // 8.1 um in run 5234
	      else
		cmsdyfctq3dHisto->fill( cmsdy*1E3 ); // 7.2 um in run 5234
	    }
	    bool inner_cut = false;
	    if(rot90) {
	      if(Q0 < 25 && Q0 > 19) inner_cut = true;
	    }
	    else if(Q0 < 23) inner_cut = true;

	    if(inner_cut) {
	      cmsskwfctq4Histo->fill(skw);
	      cmsdy0vsskwfctq4->fill( skw, cmsdy0*1E3 ); //skew vs uncorrected residual
	      cmsskwvsdy0fctq4->fill( cmsdy0*1E3, skw ); //skew vs uncorrected residual
	      cmsdyvsskwfctq4->fill( skw, cmsdy*1E3 ); //skew vs uncorrected residual

	      cmsqfctq4Histo->fill(Q0);
	      cmsdy0fctq4Histo->fill( cmsdy0*1E3 ); // for comparison, without skew correction
	      cmsdyfctq4Histo->fill( cmsdy*1E3 ); // was fctq2. 7.4 um
						  // @ 4 GeV, 19 deg

	      // Check for difference in 2col clusters:
	      if(ncol == 2) {
		// Cluster in one double column?
		if(colmin%2 == 0) { 
		  cmsdyfctq4dc1Histo->fill(cmsdy*1E3);
		  cmsdyfctq4dc1qHisto->fill(Q0);
		}
		// Cluster spread over two double columns
		else {  
		  cmsdyfctq4dc2Histo->fill(cmsdy*1E3);
		  cmsdyfctq4dc2qHisto->fill(Q0);
		}
	      }

	      if( !ldot ) {
		cmsdy0fctq4dHisto->fill( cmsdy0*1E3 ); // for comparison, without skew correction
		cmsdyfctq4dHisto->fill( cmsdy*1E3 ); // 7.2 um in run 5234
	      }
	    }
	  }

	} // CMS fiducial for y

	// accumulate cuts for x:

	if( fiducial && fabs( cmsdy ) < cuty && isolatedTrip) {

	  if( fabs( ty-0.000 ) < slope_y &&  fabs( tx-0.000 ) < slope_x ){
	    cmsdxfctHisto->fill( cmsdx*1E3 );
	    if(lowClusterCharge)
	      cmsdxfctLowChargeHisto->fill( cmsdx*1E3 );
	  }

	  if( fabs( ty-0.000 ) < slope_y &&
	      fabs( tx-0.000 ) < slope_y &&
	      Q0 > 18 ) {

	    cmsdxfctqHisto->fill( cmsdx*1E3 );

	    if( Q0 < 40 ) {
	      cmsdxfctq1Histo->fill( cmsdx*1E3 );
	    }

	    if( Q0 < 35 ) {
	      cmsdxfctq2Histo->fill( cmsdx*1E3 );
	    }

	    if( Q0 < 30 ) {
	      cmsdxfctq3Histo->fill( cmsdx*1E3 );
	    }

	  }

	} // CMS fiducial for x

	// Match CMS cluster and Upstream telescope triplet:
	if( fabs( cmsdx ) < cutx && fabs( cmsdy ) < cuty  && isolatedTrip) {
	  
	  nLinkedClusters++;
	  n_matched_clusters_dut++;
	  correvt100Histo->fill(event->getEventNumber());
	  correvt300Histo->fill(event->getEventNumber());
	  correvt1000Histo->fill(event->getEventNumber());
	  correvt4000Histo->fill(event->getEventNumber());

	  cmscolHisto->fill( c->col );
	  cmsrowHisto->fill( c->row );

	  // Fill statistics for first and last pixel in the cluster:
	  // Used to monitor VIColOr effect
	  cmscolfirstHisto->fill( c->vpix.front().col );
	  cmscollastHisto->fill( c->vpix.back().col );
	  // Statistics for first DC, charge of left and right columns:
	  if(c->vpix.front().col%2 == 0) { cmscol1stevenqHisto->fill( c->vpix.front().vcal ); }
	  else { cmscol1stoddqHisto->fill( c->vpix.front().vcal ); }

	  cmsqHisto->fill( c->charge );
	  cmsq0Histo->fill( Q0 );
	  

	  trixlkHisto->fill( -xA ); // -xA = x_DP = out
	  triylkHisto->fill( -yA ); // -yA = y_DP = up
	  trixylkHisto->fill( -xA, -yA );

	  if( lq ) {
	    cmsdxvsx->fill( xAt, cmsdx*1E3 ); // dx vs x: turn?
	    cmsdxvsy->fill( yAt, cmsdx*1E3 ); // dx vs y: rot?
	    cmsdyvsx->fill( xAt, cmsdy*1E3 ); // dy vs x: rot?
	    cmsdyvsy->fill( yAt, cmsdy*1E3 ); // dy vs y: tilt?
	    cmsdxvstx->fill( tx*1E3, cmsdx*1E3 ); // dx vs tetx: z shift?
	    cmsdyvsty->fill( ty*1E3, cmsdy*1E3 ); // dy vs tety: z shift?

	    //FIXME not correct for ETHh:
	    // should be : cmsdyvsxHisto->fill( yAt, cmsdy*1E3 );
	    // Check above plots, too
	    cmsdyvsxHisto->fill( xAt, cmsdy*1E3 );
	  }

	  if( fiducial ) {

	    cmsnpxHisto->fill( c->size );
	    cmsncolHisto->fill( ncol );
	    cmsnrowHisto->fill( nrow );
	    if(lowClusterCharge) {
	      cmsnpxLowChargeHisto->fill( c->size );
	      cmsncolLowChargeHisto->fill( ncol );
	      cmsnrowLowChargeHisto->fill( nrow );
	    }
	    cmsnrowvst1->fill( (time_now_tlu-time_event0)/fTLU, nrow ); // cluster rows vs time

	    if( lq ) cmsnrowqHisto->fill( nrow ); // no 3-rows anymore

	    if( nrow == 2 ) {
	      // For cluster only occupying one column:
	      if( ncol == 1) {
		// ETA only in odd column numbers (readout token going upwards):
		// (Check col of first pixel, since we only have one column...)
		if(c->vpix.at(0).col%2 != 0) { cmsetaoddHisto->fill( eta ); }
		// ETA only in even column numbers (readout token going downwards):
		else { cmsetaevenHisto->fill( eta ); }
	      }

	      cmsetaHisto->fill( eta );
	    }
	    
	    float qseed = 0;
	    int oddcol = int(c->col) % 2;

	    for( std::vector<CMSPixel::pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); px++ ){
	      cmspxqHisto->fill( px->vcal);
	      if( c->size == 2)
		cmspxqcl2Histo->fill( px->vcal);
	      if( nrow == 2)
		cmspxqrow2Histo->fill( px->vcal);
	      if(px->vcal > qseed)
		qseed =  px->vcal;
	      if( oddcol )
		cmspxqoddHisto->fill( px->vcal );
	      else
		cmspxqeveHisto->fill( px->vcal );
	      //if( c->size == 1 ) cmspxq1Histo->fill( px->vcal );
	      //if( c->size == 2 ) cmspxq2Histo->fill( px->vcal ); // flat
	      //if( c->size >= 3 ) cmspxq3Histo->fill( px->vcal ); // hump at low q

	      cmspxqvsq->fill( Q0, px->vcal );
	      cmspxqvsqv->fill( Q0, px->vcal );
	      cmspxqvsxm->fill( xmod, px->vcal );
	      cmspxqvsym->fill( ymod, px->vcal );
	    }
	    cmsqseedfHisto->fill(qseed);

	    cmsqfHisto->fill( c->charge );
	    cmsq0fHisto->fill( Q0 );
	    cmsq0fHistoRoot->Fill( Q0 );

	    if(seedPixelLost)
	      cmsqfOnePixeldyCutHisto->fill( c->charge );
	    else
	      cmsqfNotOnePixeldyCutHisto->fill( c->charge );
		

	    if( c->size == 1){
	      cmsqfcl1Histo->fill( c->charge );
	    }
	    if( c->size == 2){
	      cmsqfcl2Histo->fill( c->charge );
	    }

	    if( nrow == 1){
	      cmsqfrow1Histo->fill( c->charge );
	    }
	    if( nrow == 2){
	      cmsqfrow2Histo->fill( c->charge );
	    }

	    if( ldot ) { // sensor bias dot
	      cmsnpx0Histo->fill( c->size );
	      cmsqf0Histo->fill( c->charge );
	    }
	    else {
	      cmsnpx1Histo->fill( c->size );
	      cmsqf1Histo->fill( c->charge );
	      if( lcore ) { // pixel core
		cmsnpx2Histo->fill( c->size );
		cmsqf2Histo->fill( c->charge );
	      }
	      else {
		cmsnpx3Histo->fill( c->size );
		cmsqf3Histo->fill( c->charge );
	      }//core
	    }//dot

	    cmsxyHitMap->fill( xAt, yAt );
	    if(lowClusterCharge) {
	      cmsxyHitMapLowCharge->fill( xAt, yAt );
	    }
	    
	    cmsqvsx->fill( xAt, c->charge ); // cluster charge profile
	    cmsqvsy->fill( yAt, c->charge ); // cluster charge profile
	    cmsqvsxm->fill( xmod, c->charge ); //q within pixel
	    cmsqvsym->fill( ymod, c->charge ); //q within pixel
	    if(c->seed_id > -1) cmspxqsvsym->fill( ymod, c->vpix.at(c->seed_id).vcal ); //q within pixel
	    cmsqvsxmym->fill( xmod, ymod, c->charge ); // cluster charge profile
	    if(c->seed_id > -1) cmspxqvsxmym->fill( xmod, ymod, c->vpix.at(c->seed_id).vcal ); // cluster charge profile

	    if(ldot) { cmsqvsxmymdot->fill( xmod, ymod, c->charge ); } // cluster charge profile

	    // KIT: added for efficiency analysis
	    double dotsize=10;
	    double cutsize=5;

	    if( xmod >= 50-cutsize  && xmod <= 50+cutsize )  cmsqvsxm50->fill( ymod, c->charge );
	    if( xmod >= 100-dotsize && xmod <= 100+dotsize ) cmsqvsxm100->fill( ymod, c->charge );
	    if( xmod >= 150-cutsize && xmod <= 150+cutsize ) cmsqvsxm150->fill( ymod, c->charge );
	    if( xmod >= 200-dotsize && xmod <= 200+dotsize ) cmsqvsxm200->fill( ymod, c->charge );
	    if( xmod >= 250-cutsize && xmod <= 250+cutsize ) cmsqvsxm250->fill( ymod, c->charge );

	    if( ymod >= 25-cutsize  && ymod <= 25+cutsize )  cmsqvsym25->fill( xmod, c->charge );
	    if( ymod >= 50-dotsize  && ymod <= 50+dotsize )  cmsqvsym50->fill( xmod, c->charge );
	    if( ymod >= 75-cutsize  && ymod <= 75+cutsize )  cmsqvsym75->fill( xmod, c->charge );
	    if( ymod >= 100-cutsize && ymod <= 100+cutsize ) cmsqvsym100->fill( xmod, c->charge );
	    if( ymod >= 125-cutsize && ymod <= 125+cutsize ) cmsqvsym125->fill( xmod, c->charge );
	    if( ymod >= 150-dotsize && ymod <= 150+dotsize ) cmsqvsym150->fill( xmod, c->charge );
	    if( ymod >= 175-cutsize && ymod <= 175+cutsize ) cmsqvsym175->fill( xmod, c->charge );
	    // KIT end

	    cmsqvst1->fill( (time_now_tlu-time_event0)/fTLU, c->charge ); // cluster charge vs time
	    cmsqvst2->fill( (time_now_tlu-time_event0)/fTLU, c->charge ); // cluster charge vs time
	    cmsqvst3->fill( (time_now_tlu-time_event0)/fTLU, c->charge ); // cluster charge vs time
	    cmsqvst4->fill( (time_now_tlu-time_event0)/fTLU, c->charge ); // cluster charge vs time

	    cmsrmsxvsq->fill( Q0, fabs(cmsdx)*1E3 ); //resolution vs charge
	    cmsrmsyvsq->fill( Q0, fabs(cmsdy)*1E3 ); //resolution vs charge

	    double pMoyal = exp( -exp( -( Q0 - 28 ) / 3.3 ) ); // fitMoyal
	    cmspMoyalvsq->fill( Q0, pMoyal );
	    cmspMoyalHisto->fill( pMoyal );
	    cmsrmsyvsp->fill( pMoyal, fabs(cmsdy)*1E3 ); // resolution vs charge

	    if(lowClusterCharge)
	      cmspixvsxmymLowCharge->fill( xmod, ymod );
	    if(c->size == 1)
	      cmspix1vsxmym->fill( xmod, ymod );

	    if( lq ) {

	      cmsdyvsxm->fill( xmod, cmsdy*1E3 );
	      cmsdy0vsxm->fill( xmod, cmsdy0*1E3 );
	      cmsdyvsym->fill( ymod, cmsdy*1E3 );
	      cmsdy0vsym->fill( ymod, cmsdy0*1E3 );

	      cmsdxvsxm->fill( xmod, cmsdx*1E3 );
	      cmsdxvsym->fill( ymod, cmsdx*1E3 );

	      cmspixvsxmym->fill( xmod, ymod ); // occupancy map

	      // KIT: added for efficiency analysis
	      dotsize=10;
	      cutsize=5;

	      if( xmod >= 50-cutsize  && xmod <= 50+cutsize )  cmspixvsxm50->fill( ymod );
	      if( xmod >= 100-dotsize && xmod <= 100+dotsize ) cmspixvsxm100->fill( ymod );
	      if( xmod >= 150-cutsize && xmod <= 150+cutsize ) cmspixvsxm150->fill( ymod );
	      if( xmod >= 200-dotsize && xmod <= 200+dotsize ) cmspixvsxm200->fill( ymod );
	      if( xmod >= 250-cutsize && xmod <= 250+cutsize ) cmspixvsxm250->fill( ymod );

	      if( ymod >= 25-cutsize  && ymod <= 25+cutsize )  cmspixvsym25->fill( xmod );
	      if( ymod >= 50-dotsize  && ymod <= 50+dotsize )  cmspixvsym50->fill( xmod );
	      if( ymod >= 75-cutsize  && ymod <= 75+cutsize )  cmspixvsym75->fill( xmod );
	      if( ymod >= 100-cutsize && ymod <= 100+cutsize ) cmspixvsym100->fill( xmod );
	      if( ymod >= 125-cutsize && ymod <= 125+cutsize ) cmspixvsym125->fill( xmod );
	      if( ymod >= 150-dotsize && ymod <= 150+dotsize ) cmspixvsym150->fill( xmod );
	      if( ymod >= 175-cutsize && ymod <= 175+cutsize ) cmspixvsym175->fill( xmod );
	      // KIT end

	      cmsrmsxvsx->fill( xAt, fabs(cmsdx)*1E3 ); //resolution across cols
	      cmsrmsyvsx->fill( xAt, fabs(cmsdy)*1E3 ); //resolution across cols
	      cmsrmsxvsy->fill( yAt, fabs(cmsdx)*1E3 ); //resolution across rows
	      cmsrmsyvsy->fill( yAt, fabs(cmsdy)*1E3 ); //resolution across rows
	      cmsrmsxvsxm->fill( xmod, fabs(cmsdx)*1E3 ); //resolution within pixel
	      cmsrmsyvsxm->fill( xmod, fabs(cmsdy)*1E3 ); //resolution within pixel
	      cmsrmsy0vsxm->fill( xmod, fabs(cmsdy0)*1E3 ); //resolution within pixel,noskewcorr
	      cmsncolvsxm->fill( xmod, ncol );
	      cmsnrowvsxm->fill( xmod, nrow );

	      // Resolution vs. position within pixel
	      cmsrmsxvsxmym->fill( xmod, ymod, fabs(cmsdx)*1E3 );
	      cmsrmsyvsxmym->fill( xmod, ymod, fabs(cmsdy)*1E3 );
	      cmsrmsxyvsxmym->fill( xmod, ymod, fabs(sqrt(cmsdx*cmsdx+cmsdy*cmsdy))*1E3 );
	      double xmdist = fmod( xmod, 150. ) - 75.;
	      double ymdist = fmod( ymod, 100. ) - 50.;
	      double xmymdist = sqrt(xmdist*xmdist+ymdist*ymdist);
	      cmsrmsxymposvsxmym->fill( xmod, ymod, (sqrt(cmsdx*cmsdx+cmsdy*cmsdy)*1E3-xmymdist) );

	      if( !ldot ) {
		cmsrmsxvsym->fill( ymod, fabs(cmsdx)*1E3 ); //resolution within pixel
		cmsrmsyvsym->fill( ymod, fabs(cmsdy)*1E3 ); //resolution within pixel
		cmsrmsyvsym3->fill( ymd3, fabs(cmsdy)*1E3 ); //resolution within pixel
		cmsrmsyvsym6->fill( ymd6, fabs(cmsdy)*1E3 ); //resolution within pixel
	      }
	      cmsrmsyvst->fill( (time_now_tlu-time_event0)/fTLU, fabs(cmsdy)*1E3 ); //resolution vs time

	      if( nrow <= 2 ) {
		cmsdyvseta->fill( eta, cmsdy*1E3 );
		cmsrmsyvseta->fill( eta, fabs(cmsdy)*1E3 );
	      }

	      cmsnpxvsxmym->fill( xmod, ymod, c->size ); // cluster
							 // size map
	      if(c->size ==  1) cmsnpx1vsxmym->fill( xmod, ymod, 1); // cluster
	      if(c->size ==  2) cmsnpx2vsxmym->fill( xmod, ymod, 1); // cluster size
	      if(c->size ==  3) cmsnpx3vsxmym->fill( xmod, ymod, 1); // cluster size
	      if(c->size ==  4) cmsnpx4vsxmym->fill( xmod, ymod, 1); // cluster size map
	      cmsncolvsym->fill( ymod, ncol ); // within pixel
	      cmsnrowvsym->fill( ymod, nrow ); // within pixel

	      // KIT: added for efficiency analysis
	      double dotsize=10;
	      double cutsize=5;

	      if( xmod >= 50-cutsize  && xmod <= 50+cutsize )  cmsnpxvsxm50->fill( ymod, c->size );
	      if( xmod >= 100-dotsize && xmod <= 100+dotsize ) cmsnpxvsxm100->fill( ymod, c->size );
	      if( xmod >= 150-cutsize && xmod <= 150+cutsize ) cmsnpxvsxm150->fill( ymod, c->size );
	      if( xmod >= 200-dotsize && xmod <= 200+dotsize ) cmsnpxvsxm200->fill( ymod, c->size );
	      if( xmod >= 250-cutsize && xmod <= 250+cutsize ) cmsnpxvsxm250->fill( ymod, c->size );

	      if( ymod >= 25-cutsize  && ymod <= 25+cutsize )  cmsnpxvsym25->fill( xmod, c->size );
	      if( ymod >= 50-dotsize  && ymod <= 50+dotsize )  cmsnpxvsym50->fill( xmod, c->size );
	      if( ymod >= 75-cutsize  && ymod <= 75+cutsize )  cmsnpxvsym75->fill( xmod, c->size );
	      if( ymod >= 100-cutsize && ymod <= 100+cutsize ) cmsnpxvsym100->fill( xmod, c->size );
	      if( ymod >= 125-cutsize && ymod <= 125+cutsize ) cmsnpxvsym125->fill( xmod, c->size );
	      if( ymod >= 150-dotsize && ymod <= 150+dotsize ) cmsnpxvsym150->fill( xmod, c->size );
	      if( ymod >= 175-cutsize && ymod <= 175+cutsize ) cmsnpxvsym175->fill( xmod, c->size );
	      // KIT end

	      cmsetavsym->fill( ymod, eta ); // eta within pixel
	      if( nrow == 2 ) cmsetavsym2->fill( ymod, eta ); // eta within pixel
	      cmsetavsym3->fill( ymd3, eta ); // eta within pixel

	      if(      nrow == 1 ) 
		cmsym1Histo->fill( ymod ); //where are 1-rows?
	      else if( nrow == 2 ) 
		cmsym2Histo->fill( ymod ); //where are 2-rows?

	      if(ncol > 1) { 
		cmsskwvsym->fill( ymod, skw ); //skew within pixel
		cmsskwvsxm->fill( xmod, skw ); //skew within pixel
		cmsdy0vsskw->fill( skw, cmsdy0*1E3 ); //skew vs uncorrected residual
		cmsskwvsdy0->fill( cmsdy0*1E3, skw ); //skew vs uncorrected residual
		cmsdyvsskw->fill( skw, cmsdy*1E3 ); //skew vs corrected residual
		cmsskwfcqHisto->fill(skw);
	      }

	    } // q Landau peak

	    if(ncol > 1) { 
	      cmsskwfcHisto->fill(skw);
	      if(ncol == 3) cmsskw3pxHisto->fill(skw);
	      if(ncol == 4) cmsskw4pxHisto->fill(skw);
	    }

	  } // fiducial

	} // CMS - triplet match

	//Plot 2 Cluster distance for Tracks with 2 Clusters
	if(ClustDUT.size() == 2 && c == (ClustDUT.begin()+1) ){
	  double cmsxFirstClus = ( ClustDUT.begin()->col - 26 ) * pitchcol; // -3.9..3.9 mm
	  double cmsyFirstClus = ( ClustDUT.begin()->row - 40 ) * pitchrow; // -4..4 mm
	  double twoClusterDistance = sqrt( fabs(cmsxFirstClus - cmsx)*fabs(cmsxFirstClus - cmsx) + fabs(cmsyFirstClus - cmsy)*fabs(cmsyFirstClus - cmsy) );
	  twoClusterDistanceHisto->fill(twoClusterDistance);
	  twoClusterXDistanceHisto->fill(fabs(cmsxFirstClus - cmsx));
	  twoClusterYDistanceHisto->fill(fabs(cmsyFirstClus - cmsy));
	  if(trackHasLostSeedPixel){
	    twoClusterDistanceLostSeedHisto->fill(twoClusterDistance);
	    twoClusterXDistanceLostSeedHisto->fill(fabs(cmsxFirstClus - cmsx));
	    twoClusterYDistanceLostSeedHisto->fill(fabs(cmsyFirstClus - cmsy));
	  }
	  if(nLinkedClusters >= 1){
	    twoClusterDistanceLinkedTrackHisto->fill(twoClusterDistance);
	    twoClusterXDistanceLinkedTrackHisto->fill(fabs(cmsxFirstClus - cmsx));
	    twoClusterYDistanceLinkedTrackHisto->fill(fabs(cmsyFirstClus - cmsy));
	    if(trackHasLostSeedPixel){
	      twoClusterDistanceLostSeedLinkedTrackHisto->fill(twoClusterDistance);
	      twoClusterXDistanceLostSeedLinkedTrackHisto->fill(fabs(cmsxFirstClus - cmsx));
	      twoClusterYDistanceLostSeedLinkedTrackHisto->fill(fabs(cmsyFirstClus - cmsy));
	    }
	  }
	}

	if( fabs( cmsdx ) < 0.5 && fabs( cmsdy ) < 0.5 ){ // link to CMS
	  (*trip).linked_dut = true;
	  (*trip).cmsdx = cmsdx;
	  (*trip).cmsdy = cmsdy;
	  ntrilk++;

	}

	// DUT alignment using GBL track fitting:
	if( fiducial && fabs( cmsdx ) < cutx && fabs( cmsdy ) < cuty ) {
	  //if( fiducial && abs( cmsdx ) < 0.2 && abs( cmsdy ) < 0.125) {

	  //FIXME move into function addGBLmeasurement(int planeID) with new geometry?

	  // GBL point vector for the trajectory
	  std::vector<gbl::GblPoint> traj_points;

	  TMatrixD proL2m(2,2);
	  proL2m.UnitMatrix();

	  TVectorD meas(2);

	  TVectorD measPrec(2); // precision = 1/resolution^2
	  measPrec[0] = 1.0 / telescope_resx / telescope_resx;
	  measPrec[1] = 1.0 / telescope_resy / telescope_resy;

	  // Adding scatter:
	  TVectorD scat(2);
	  scat.Zero(); //mean is zero

	  double X0Si = 65e-3 / 94; // Si + Kapton
	  double tetSi = 0.0136 * sqrt(X0Si) / _eBeam * ( 1 + 0.038*std::log(X0Si) );

	  TVectorD wscatSi(2);
	  wscatSi[0] = 1.0 / ( tetSi * tetSi ); //weight
	  wscatSi[1] = 1.0 / ( tetSi * tetSi );

	  // Adding the Upstream Telescope Layer Measurements:
	  double zprev = _planePosition[0];
	  double step = 0;

	  for( int ipl = 0; ipl < 3; ++ipl ){

	    step = _planePosition[ipl] - zprev;
	    zprev = _planePosition[ipl];

	    gbl::GblPoint * point = new gbl::GblPoint( JacobianPointToPoint( step ) );

	    meas[0] = (*trip).getdx(ipl);
	    meas[1] = (*trip).getdy(ipl);

	    point->addMeasurement( proL2m, meas, measPrec );
	    point->addScatterer( scat, wscatSi );
	    traj_points.push_back(*point);

	    delete point;

	  } // loop over planes

	  // Adding the DUT Measurement:
	  step = DUTz - zprev;

	  gbl::GblPoint * point = new gbl::GblPoint( JacobianPointToPoint( step ) );

	  double tetDUT = 0.0136 * sqrt(DUTX0) / _eBeam * ( 1 + 0.038*std::log(DUTX0) );

	  TVectorD wscatDUT(2);
	  wscatDUT[0] = 1.0 / ( tetDUT * tetDUT ); //weight
	  wscatDUT[1] = 1.0 / ( tetDUT * tetDUT );

	  point->addScatterer( scat, wscatDUT );


	  naldut++;
	  meas[0] = (*trip).cmsdx;
	  meas[1] = (*trip).cmsdy;

	  double resx = 48E-3; // [mm] CMS col resolution
	  if( abs( turn ) > 11 ) resx = 22E-3;
	  double resy = 12E-3; // [mm] CMS row resolution at 20 deg tilt
	  resy = 22E-3; // same weight as x

	  TVectorD measWeight(2);
	  measWeight[0] = 1.0 / resx / resx; // weight = 1/resolution^2
	  measWeight[1] = 1.0 / resy / resy;

	  // rot*tilt*turn:

	  proL2m[0][0] = cf*co+sf*sa*so;
	  proL2m[0][1] = sf*ca;
	  proL2m[1][0] =-sf*co+cf*sa*so;
	  proL2m[1][1] = cf*ca;

	  // inverse:

	  point->addMeasurement( proL2m, meas, measWeight );

	  // DUT alignment: dx, dy, drot, dtilt, dturn, (dz)

	  int nAlignmentParametersForDUT = 6;

	  // Adding the DUT alignment derivatives for the MillepedeII alignment:
	  TMatrixD DUTDer( 2, nAlignmentParametersForDUT ); // alignment derivatives
	  // residuals in DUT plane
	  // transformed track impact point coordinates
	  // residx = col - x4;
	  // residy = row - y4;
	  // alignment = linearized correction around initial values
	  // may need to iterate

	  DUTDer[0][0] =-1.0; // dresidx/ddeltax
	  DUTDer[1][0] = 0.0;

	  DUTDer[0][1] = 0.0;
	  DUTDer[1][1] =-1.0; // dresidy/ddeltay

	  DUTDer[0][2] =-upsign*y2; // dresidx/drot, linearized
	  DUTDer[1][2] =-upsign*x2; // dresidy/drot

	  DUTDer[0][3] =-upsign*( cf*(-so*dzcda)+sf*(-sa*y1 + ca*z1 + sa*co*dzcda)); // dresidx/dtilt a
	  DUTDer[1][3] = upsign*(-sf*(-so*dzcda)+cf*(-sa*y1 + ca*z1 + sa*co*dzcda)); // dresidy/dtilt a

	  derxtiltHisto->fill( DUTDer[0][3] );
	  derytiltHisto->fill( DUTDer[1][3] );

	  DUTDer[0][4] =-upsign*( cf*(-so*xc - co*dzc - so*dzcdo) + sf*sa*(co*xc-so*dzc+co*dzcdo)); // dresidx/dturn o
	  DUTDer[1][4] = upsign*(-sf*(-so*xc - co*dzc - so*dzcdo) + cf*sa*(co*xc-so*dzc+co*dzcdo)); // dresidy/dturn o

	  derxturnHisto->fill( DUTDer[0][4] );
	  deryturnHisto->fill( DUTDer[1][4] );

	  // dz:

	  DUTDer[0][5] =-upsign*( cf*(co*dxcdz-so*ddzcdz) + sf*(ca*dycdz+sa*(so*dxcdz+co*ddzcdz))); // dresidx/dz
	  DUTDer[1][5] = upsign*(-sf*(co*dxcdz-so*ddzcdz) + cf*(ca*dycdz+sa*(so*dxcdz+co*ddzcdz))); // dresidy/dz

	  // global labels for Pede:

	  std::vector<int> DUTLabels( nAlignmentParametersForDUT );
	  DUTLabels[0] = 1; // dx
	  DUTLabels[1] = 2; // dy
	  DUTLabels[2] = 3; // drot
	  DUTLabels[3] = 4; // dtilt
	  DUTLabels[4] = 5; // dturn
	  DUTLabels[5] = 6; // dz

	  point->addGlobals( DUTLabels, DUTDer ); // for MillePede alignment
	  traj_points.push_back(*point);
	  delete point;

	  double Chi2;
	  int Ndf;
	  double lostWeight;

	  gbl::GblTrajectory traj(traj_points, false); // curvature = false
	  traj.fit( Chi2, Ndf, lostWeight );

	  streamlog_out(DEBUG2) << "Triplet GBL Fit: Chi2/Ndf=" << Chi2/Ndf 
				<< " lostWeight=" << lostWeight << std::endl;

	  traj.milleOut( *mille ); // write to mille.bin

	} // linked triplet

      } // loop over CMS clusters

      nLinkedTripClus->fill(nLinkedClusters);
      nTripClus->fill(ClustDUT.size());
      if(trackHasLostSeedPixel)
	nTripClusLostSeed->fill(ClustDUT.size());
      int numberOfTripPixels = 0;
      for( std::vector<cluster>::iterator c = ClustDUT.begin(); c != ClustDUT.end(); c++ ){
	numberOfTripPixels += c->size;
      }
      nTripPixels->fill(numberOfTripPixels);
      if(trackHasLostSeedPixel)
	nTripPixelsLostSeed->fill(numberOfTripPixels);

      if(nLinkedClusters >= 1){
	nTripClusLinkedTrack->fill(ClustDUT.size());
	if(trackHasLostSeedPixel)
	  nTripClusLostSeedLinkedTrack->fill(ClustDUT.size());
	int numberOfTripPixels = 0;
	for( std::vector<cluster>::iterator c = ClustDUT.begin(); c != ClustDUT.end(); c++ ){
	  numberOfTripPixels += c->size;
	}
	nTripPixelsLinkedTrack->fill(numberOfTripPixels);
	if(trackHasLostSeedPixel)
	  nTripPixelsLostSeedLinkedTrack->fill(numberOfTripPixels);
      }

    }// have some CMS hit

  }// iterate over upstream triplets

    
  ntriHisto->fill( upstream_triplets->size() );
  lkAvst->fill( (time_now_tlu-time_event0)/fTLU, ntrilk );



  //----------------------------------------------------------------------------
  // telescope triplet efficiency w.r.t. CMS pixel:
  if( dutPixels->size() > 0 ){

    // CMS pixel clusters:

    for( std::vector<cluster>::iterator c = ClustDUT.begin(); c != ClustDUT.end(); c++ ){

      double cmsx = ( c->col - 26 ) * pitchcol; // -3.9..3.9 mm
      double cmsy = ( c->row - 40 ) * pitchrowt; // -4..4 mm

      if(rot90) {
	cmsx = ( c->row - 40 ) * pitchrow; // old style
	cmsy = ( 26 - c->col ) * pitchcol; // -3.9..3.9 mm
      }
      else if( FPIX) {
	// FPIX rot at 90 deg:
	// col = y
	// row = x
	cmsx = ( c->row - 40 ) * pitchrowt; // old style
	cmsy = ( c->col - 26 ) * pitchcolt; // -3.9..3.9 mm
      }
      else if( ETHh ) {
	cmsx = ( 40 - c->row ) * pitchrowt; // KIT orientation
	cmsy = ( c->col - 26 ) * pitchcolt; // -3.9..3.9 mm
      }

      double cmsxe = cmsx + DUTrot*cmsy - DUTalignx; // aligned w.r.t. telescope
      double cmsye = cmsy - DUTrot*cmsx - DUTaligny; // 

      cmsxeHisto->fill( cmsxe );
      cmsyeHisto->fill( cmsye );

      // triplets:

      int  nm = 0;
      bool im = 0;

      for( std::vector<triplet>::iterator trip = upstream_triplets->begin(); trip != upstream_triplets->end(); trip++ ) {

	double sx = cmsxe + (*trip).getx_at(DUTz);
	double dy = cmsye - (*trip).gety_at(DUTz);

	if( _DUT_chip >= 100 ) { // xdb and dig mounted upright
	  sx = cmsxe - (*trip).getx_at(DUTz);
	  dy = cmsye + (*trip).gety_at(DUTz);
	}

	cmsdxeHisto->fill( sx*1E3 );
	cmsdyeHisto->fill( dy*1E3 );

	if( fabs(sx) < 0.300 && fabs(dy) < 0.200 ) {
	  nm++;
	  im = 1;
	}

      }//triplets

      cmsnmHisto->fill( nm );
      trieffvsxy->fill( cmsxe, cmsye, im ); //efficiency profile

    }//clusters
  }//have CMS data



   //----------------------------------------------------------------------------
   // six: triplets A and driplets B
   // matching and GBL fit
   // kinks: triplets A vs driplets B
   // scattering point = DUT:


   // Match the Telescope Upstream and Downstream Arm triplets to get tracks:
  std::vector<track> * telescope_tracks = new std::vector<track>;
  telescope_tracks = MatchTriplets(upstream_triplets,downstream_triplets,DUTz);





  for( std::vector<track>::iterator tr = telescope_tracks->begin(); tr != telescope_tracks->end(); tr++ ){

    triplet trip = (*tr).get_upstream();
    triplet drip = (*tr).get_downstream();


    // Track kinks as difference in triplet slopes:
    //      double kx = (*drip).slope().x - (*trip).slope().x; //kink
    //      double ky = (*drip).slope().y - (*trip).slope().y;
    double kx = (*tr).kink_x();
    double ky = (*tr).kink_y();

    // Track impact position at DUT from Downstream:
    double xB = drip.getx_at(DUTz);
    double yB = drip.gety_at(DUTz);

    // Track impact position at DUT from Upstream:
    double xA = trip.getx_at(DUTz);
    double yA = trip.gety_at(DUTz);

    double dx = xB - xA; // driplet - triplet
    double dy = yB - yA;

    // Track impact position at REF from Downstream:
    double xR = drip.getx_at(REFz);
    double yR = drip.gety_at(REFz);



    double probchi = 0;

    double avx = trip.base().x;
    double avy = trip.base().y;
    double avz = trip.base().z;
    double tx = trip.slope().x;
    double ty = trip.slope().y;
    double zA = DUTz - avz; // z CMS from mid of triplet

    // intersect inclined track with tilted DUT plane:

    double zc = (Nz*zA - Ny*avy - Nx*avx) / (Nx*tx + Ny*ty + Nz); // from avz
    double yc = avy + ty * zc;
    double xc = avx + tx * zc;
    double dzc = zc + avz - DUTz; // from DUT z0 [-8,8] mm

    // transform into DUT system: (passive).
    // large rotations don't commute: careful with order

    double x1 = co*xc - so*dzc; // turn
    double y1 = yc;
    double z1 = so*xc + co*dzc;

    double x2 = x1;
    double y2 = ca*y1 + sa*z1; // tilt
    // double z2 =-sa*y1 + ca*z1; // should be zero (in DUT plane). is zero

    double x3 = cf*x2 + sf*y2; // rot
    double y3 =-sf*x2 + cf*y2;

    double upsign = -1; // analog psi46: upside down
    if( _DUT_chip >= 100 )
      upsign = 1; // xdb and dig mounted upright

    double x4 = upsign*x3 + DUTalignx; // shift to mid
    if(ETHh) x4 = upsign*x3 - DUTalignx;

    double y4 =-upsign*y3 + DUTaligny; // invert y, shift to mid

    // new style:
    double xAt = x4;
    double yAt = y4;

    double xmod = fmod( 9.075 + xAt, 0.3 ) * 1E3; // [0,300] um, 2 pixel wide
    double ymod = fmod( 9.050 + yAt, 0.2 ) * 1E3; // [0,200] um

    if(hanging) {
      xmod = fmod( 9.075 - xAt, 0.3 ) * 1E3; // [0,300] um, 2 pixel wide
      ymod = fmod( 9.050 - yAt, 0.2 ) * 1E3; // [0,200] um
    }
    if (FPIX || ETHh){
      xmod = fmod( 9.075 + yAt, 0.3 ) * 1E3; // [0,300] um, 2 pixel wide
      ymod = fmod( 9.030 + xAt, 0.2 ) * 1E3; // [0,200] um
    }

    // Get Kink angle

    double fidxmax =  3.8;
    double fidxmin = -3.8;
    double fidymax =  3.8;
    double fidymin = -3.8;

    if( xAt > fidxmin && xAt < fidxmax && yAt > fidymin && yAt < fidymax ) {
      double kinkangle = (kx*kx+ky*ky)*1E6;
      kinkvsxmym->fill( xmod, ymod, kinkangle ); // kink angle
      kinkvsxm->fill(xmod, kinkangle);
      kinkvsym->fill(ymod, kinkangle);
      kinksqdutl->fill(kinkangle);
      kinksqduth->fill(kinkangle);
      kinkxdut->fill(kx*1E3);
      kinkydut->fill(ky*1E3);
    }


    // GBL with triplet A as seed:

    std::vector<gbl::GblPoint> traj_points;
    //GblTrajectory traj( false ); // curvature = false

    // build up trajectory:

    std::vector<double> sPoint;

    // plane 0:

    double s = 0;

    TMatrixD proL2m(2,2);
    proL2m.UnitMatrix();

    TVectorD meas(2);

    double res = 3.5E-3; // [mm] Anemone telescope intrinsic resolution
    //res = 4.5E-3; // EUDET

    TVectorD measPrec(2); // precision = 1/resolution^2
    measPrec[0] = 1.0 / res / res;
    measPrec[1] = 1.0 / res / res;

    // scatter:
    TVectorD scat(2);
    scat.Zero(); //mean is zero

    double p = _eBeam; // beam momentum
    double X0Si = 65e-3 / 94; // Si + Kapton
    double tetSi = 0.0136 * sqrt(X0Si) / p * ( 1 + 0.038*std::log(X0Si) );

    TVectorD wscatSi(2);
    wscatSi[0] = 1.0 / ( tetSi * tetSi ); //weight
    wscatSi[1] = 1.0 / ( tetSi * tetSi );

    TMatrixD alDer( 2, 3 ); // alignment derivatives
    alDer[0][0] = 1.0; // dx/dx
    alDer[1][0] = 0.0; // dy/dx
    alDer[0][1] = 0.0; // dx/dy
    alDer[1][1] = 1.0; // dy/dy

    std::vector<unsigned int> ilab; // 0-5 = telescope, 6 = DUT, 7 = REF

    // plane 0-5:

    double rx[6];
    double ry[6];
    double zprev = _planePosition[0];

    for( int ipl = 0; ipl < 6; ++ipl ){

      // Get the corresponding hit from the track:
      hit trackhit = (*tr).gethit(ipl);

      double step = _planePosition[ipl] - zprev;
      zprev = _planePosition[ipl];

      gbl::GblPoint * point = new gbl::GblPoint( JacobianPointToPoint( step ) );
      s += step;

      double dz = trackhit.z - trip.base().z;
      double xs = trip.base().x + trip.slope().x * dz; // Ax at plane
      double ys = trip.base().y + trip.slope().y * dz; // Ay at plane

      rx[ipl] = trackhit.x - xs;
      ry[ipl] = trackhit.y - ys;

      if( ipl == 0 ) {
	sixx0Histo->fill( -trackhit.x );
	sixy0Histo->fill( -trackhit.y );
      }
      if( ipl == 1 ) {
	sixx1Histo->fill( -trackhit.x );
	sixy1Histo->fill( -trackhit.y );
      }
      if( ipl == 2 ) {
	sixx2Histo->fill( -trackhit.x );
	sixy2Histo->fill( -trackhit.y );
      }
      if( ipl == 3 ) {
	sixx3Histo->fill( -trackhit.x );
	sixy3Histo->fill( -trackhit.y );
      }
      if( ipl == 4 ) {
	sixx4Histo->fill( -trackhit.x );
	sixy4Histo->fill( -trackhit.y );
      }
      if( ipl == 5 ) {
	sixx5Histo->fill( -trackhit.x );
	sixy5Histo->fill( -trackhit.y );
      }

      meas[0] = rx[ipl];
      meas[1] = ry[ipl];

      point->addMeasurement( proL2m, meas, measPrec );

      point->addScatterer( scat, wscatSi );

      std::vector<int> globalLabels(3);
      globalLabels[0] = 10 + ipl; // dx
      globalLabels[1] = 20 + ipl; // dy
      globalLabels[2] = 40 + ipl; // rot
      alDer[0][2] = -ys; // dx/rot
      alDer[1][2] =  xs; // dy/rot

      //DP 2013 point->addGlobals( globalLabels, alDer ); // for MillePede alignment
      traj_points.push_back(*point);
      //unsigned int iLabel = traj.addPoint(*point);
      //ilab[ipl] = iLabel;
      sPoint.push_back( s );

      delete point;

      if( lDUT && ipl == 2 ) { // insert DUT

	step = DUTz - zprev;
	zprev = DUTz;

	point = new gbl::GblPoint( JacobianPointToPoint( step ) );
	s += step;

	double tetDUT = 0.0136 * sqrt(DUTX0) / p * ( 1 + 0.038*std::log(DUTX0) );

	TVectorD wscatDUT(2);
	wscatDUT[0] = 1.0 / ( tetDUT * tetDUT ); //weight
	wscatDUT[1] = 1.0 / ( tetDUT * tetDUT );

	point->addScatterer( scat, wscatDUT );

	// DUT measurement:

	if( trip.linked_dut ) {

	  meas[0] = trip.cmsdx;
	  meas[1] = trip.cmsdy;

	  double resx = 48E-3; // [mm] CMS col resolution
	  if( fabs( turn ) > 11 ) resx = 22E-3;
	  double resy =  8E-3; // [mm] CMS row resolution at 20 deg tilt

	  TVectorD measWeight(2);
	  measWeight[0] = 1.0 / resx / resx; // weight = 1/resolution^2
	  measWeight[1] = 1.0 / resy / resy;

	  point->addMeasurement( proL2m, meas, measWeight );

	}//lkA linked hit in DUT

	traj_points.push_back(*point);
	//iLabel = traj.addPoint(*point);

	//ilab[6] = iLabel;

	delete point;

      }//DUT present

    } // loop over planes

    // REF:

    // monitor what we put into GBL:

    selxHisto->fill( -xA ); // triplet at DUT
    selyHisto->fill( -yA );
    selaxHisto->fill( trip.slope().x*1E3 );
    selayHisto->fill( trip.slope().y*1E3 );
    seldxHisto->fill( dx*1E3 ); // triplet-driplet match
    seldyHisto->fill( dy*1E3 );
    selkxHisto->fill( kx*1E3 ); // triplet-driplet kink
    selkyHisto->fill( ky*1E3 );

    seldx1Histo->fill( rx[1]*1E3 ); // triplet interpol
    seldy1Histo->fill( ry[1]*1E3 );
    seldx3Histo->fill( rx[3]*1E3 ); // triplet extrapol
    seldy3Histo->fill( ry[3]*1E3 );
    seldx4Histo->fill( rx[4]*1E3 );
    seldy4Histo->fill( ry[4]*1E3 );
    seldx5Histo->fill( rx[5]*1E3 );
    seldy5Histo->fill( ry[5]*1E3 );
    if( trip.linked_dut ) {
      seldx6Histo->fill( trip.cmsdx*1E3 );
      seldy6Histo->fill( trip.cmsdy*1E3 );
    }

    double Chi2;
    int Ndf;
    double lostWeight;

    gbl::GblTrajectory traj(traj_points, false ); // curvature = false
    traj.fit( Chi2, Ndf, lostWeight );
    traj.getLabels(ilab);

    ngbl++;


    gblndfHisto->fill( Ndf );
    if( Ndf == 8 ) 
      gblchi2aHisto->fill( Chi2 );
    else
      gblchi2bHisto->fill( Chi2 );
    probchi = TMath::Prob( Chi2, Ndf );
    gblprbHisto->fill( probchi );

    // bad fits:

    if( probchi < 0.01 ) {

      badxHisto->fill( -xA ); // triplet at DUT
      badyHisto->fill( -yA );
      badaxHisto->fill( trip.slope().x*1E3 );
      badayHisto->fill( trip.slope().y*1E3 );
      baddxHisto->fill( dx*1E3 ); // triplet-driplet match
      baddyHisto->fill( dy*1E3 );
      badkxHisto->fill( kx*1E3 ); // triplet-driplet kink
      badkyHisto->fill( ky*1E3 );

      baddx1Histo->fill( rx[1]*1E3 ); // triplet interpol
      baddy1Histo->fill( ry[1]*1E3 );
      baddx3Histo->fill( rx[3]*1E3 ); // triplet extrapol
      baddy3Histo->fill( ry[3]*1E3 );
      baddx4Histo->fill( rx[4]*1E3 );
      baddy4Histo->fill( ry[4]*1E3 );
      baddx5Histo->fill( rx[5]*1E3 );
      baddy5Histo->fill( ry[5]*1E3 );
      if( trip.linked_dut ) {
	baddx6Histo->fill( trip.cmsdx*1E3 );
	baddy6Histo->fill( trip.cmsdy*1E3 );
      }

    }// bad fit

    else {

      goodx1Histo->fill( rx[1]*1E3 ); // triplet interpol
      goody1Histo->fill( ry[1]*1E3 );
      if( trip.linked_dut ) {
	goodx6Histo->fill( trip.cmsdx*1E3 );
	goody6Histo->fill( trip.cmsdy*1E3 );
      }

    } // OK fit

    // look at fit:

    TVectorD aCorrection(5);
    TMatrixDSym aCovariance(5);

    double ax[8];
    // double ay[8];
    unsigned int k = 0;

    // at plane 0:

    int ipos = ilab[0];
    traj.getResults( ipos, aCorrection, aCovariance );

    unsigned int ndim = 2;
    TVectorD aResiduals(ndim);
    TVectorD aMeasErrors(ndim);
    TVectorD aResErrors(ndim);
    TVectorD aDownWeights(ndim);

    traj.getMeasResults( static_cast<unsigned int>(ipos), ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );

    TVectorD aKinks(ndim);
    TVectorD aKinkErrors(ndim);
    TVectorD kResErrors(ndim);
    TVectorD kDownWeights(ndim);
    traj.getScatResults( static_cast<unsigned int>(ipos), ndim, aKinks, aKinkErrors, kResErrors, kDownWeights );

    //track = q/p, x', y', x, y
    //        0,   1,  2,  3, 4

    gblax0Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
    gbldx0Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
    gblrx0Histo->fill( ( rx[0] - aCorrection[3] ) * 1E3 ); // residual x [um]
    gblpx0Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
    gblqx0Histo->fill( aKinks[0]*1E3 ); // kink
    ax[k] = aCorrection[1]; // angle correction at plane, for kinks
    // ay[k] = aCorrection[2]; // angle correction at plane, for kinks
    k++;

    ipos = ilab[1];
    traj.getResults( ipos, aCorrection, aCovariance );
    traj.getMeasResults(static_cast<unsigned int>(ipos), ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );
    traj.getScatResults(static_cast<unsigned int>(ipos), ndim, aKinks, aKinkErrors, kResErrors, kDownWeights );
    gblax1Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
    gbldx1Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
    gblrx1Histo->fill( ( rx[1] - aCorrection[3] ) * 1E3 ); // residual x [um]
    gblpx1Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
    gblqx1Histo->fill( aKinks[0]*1E3 ); // kink
    gblsx1Histo->fill( aKinks[0]/aKinkErrors[0] ); // x kink pull
    gbltx1Histo->fill( aKinks[0]/kResErrors[0] ); // x kink pull
    ax[k] = aCorrection[1]; // angle correction at plane, for kinks
    //ay[k] = aCorrection[2]; // angle correction at plane, for kinks
    k++;

    ipos = ilab[2];
    traj.getResults( ipos, aCorrection, aCovariance );
    traj.getMeasResults(static_cast<unsigned int>(ipos), ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );
    traj.getScatResults(static_cast<unsigned int>(ipos), ndim, aKinks, aKinkErrors, kResErrors, kDownWeights );
    gblax2Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
    gbldx2Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
    gblrx2Histo->fill( ( rx[2] - aCorrection[3] ) * 1E3 ); // residual x [um]
    gblpx2Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
    gblqx2Histo->fill( aKinks[0]*1E3 ); // kink
    ax[k] = aCorrection[1]; // angle correction at plane, for kinks
    //ay[k] = aCorrection[2]; // angle correction at plane, for kinks
    k++;

    if( lDUT ) {
      ipos = ilab[6]; // 6 = DUT
      traj.getResults( ipos, aCorrection, aCovariance );
      traj.getScatResults(static_cast<unsigned int>(ipos), ndim, aKinks, aKinkErrors, kResErrors, kDownWeights );
      gblax6Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
      gbldx6Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
      gbldy6Histo->fill( aCorrection[4]*1E3 ); // shift y [um]
      gblqx6Histo->fill( aKinks[0]*1E3 ); // x kink
      gblsx6Histo->fill( aKinks[0]/aKinkErrors[0] ); // x kink pull
      gbltx6Histo->fill( aKinks[0]/kResErrors[0] ); // x kink pull

      if( trip.linked_dut ) {
	traj.getMeasResults(static_cast<unsigned int>(ipos), ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );
	gblrx6Histo->fill( ( trip.cmsdx - aCorrection[3] ) * 1E3 ); // residual x [um]
	gblry6Histo->fill( ( trip.cmsdy - aCorrection[4] ) * 1E3 ); // residual y [um]
	gblpx6Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
	gblpy6Histo->fill( aResiduals[1] / aResErrors[1] ); // pull
      }
      ax[k] = aCorrection[1]; // angle correction at plane, for kinks
      //ay[k] = aCorrection[2]; // angle correction at plane, for kinks
      k++;
    }//DUT

    ipos = ilab[3];
    traj.getResults( ipos, aCorrection, aCovariance );
    traj.getMeasResults(static_cast<unsigned int>(ipos), ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );
    traj.getScatResults(static_cast<unsigned int>(ipos), ndim, aKinks, aKinkErrors, kResErrors, kDownWeights );
    gblax3Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
    gbldx3Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
    gblrx3Histo->fill( ( rx[3] - aCorrection[3] ) * 1E3 ); // residual x [um]
    gblpx3Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
    gblqx3Histo->fill( aKinks[0]*1E3 ); // kink
    ax[k] = aCorrection[1]; // angle correction at plane, for kinks
    //ay[k] = aCorrection[2]; // angle correction at plane, for kinks
    k++;

    ipos = ilab[4];
    traj.getResults( ipos, aCorrection, aCovariance );
    traj.getMeasResults(static_cast<unsigned int>(ipos), ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );
    traj.getScatResults(static_cast<unsigned int>(ipos), ndim, aKinks, aKinkErrors, kResErrors, kDownWeights );
    gblax4Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
    gbldx4Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
    gblrx4Histo->fill( ( rx[4] - aCorrection[3] ) * 1E3 ); // residual x [um]
    gblpx4Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
    gblqx4Histo->fill( aKinks[0]*1E3 ); // kink
    ax[k] = aCorrection[1]; // angle correction at plane, for kinks
    //ay[k] = aCorrection[2]; // angle correction at plane, for kinks
    k++;

    ipos = ilab[5];
    traj.getResults( ipos, aCorrection, aCovariance );
    traj.getMeasResults(static_cast<unsigned int>(ipos), ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );
    traj.getScatResults(static_cast<unsigned int>(ipos), ndim, aKinks, aKinkErrors, kResErrors, kDownWeights );
    gblax5Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
    gbldx5Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
    gblrx5Histo->fill( ( rx[5] - aCorrection[3] ) * 1E3 ); // residual x [um]
    gblpx5Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
    gblqx5Histo->fill( aKinks[0]*1E3 ); // kink
    ax[k] = aCorrection[1]; // angle correction at plane, for kinks
    //ay[k] = aCorrection[2]; // angle correction at plane, for kinks
    k++;

    // kinks: 1,2 = tele, 3 = DUT, 4,5 = tele

    gblkx1Histo->fill( (ax[1] - ax[0])*1E3 ); // kink at 1 [mrad]
    gblkx2Histo->fill( (ax[2] - ax[1])*1E3 ); // kink at 2 [mrad]
    gblkx3Histo->fill( (ax[3] - ax[2])*1E3 ); // kink at 3 [mrad]
    gblkx4Histo->fill( (ax[4] - ax[3])*1E3 ); // kink at 4 [mrad]
    gblkx5Histo->fill( (ax[5] - ax[4])*1E3 ); // kink at 5 [mrad]
    gblkx6Histo->fill( (ax[6] - ax[5])*1E3 ); // kink at 6 [mrad]



    //------------------------------------------------------------------------
    // CMS DUT efficiency:

    if( fabs(dx) < 0.1 && fabs(dy) < 0.1 && drip.linked_dut ) { // six with link

      sixxylkHisto->fill( -xA, -yA ); // six-tracks with REF link at CMS DUT

      bool nm = 0;
      if( trip.linked_dut ) nm = 1;

      double fidxmax =  3.8;
      double fidxmin = -3.8;
      double fidymax =  3.8;
      double fidymin = -3.8;

      if( abs(turn) > 2 ) {
	fidxmin = -3.0 - ( yAt + 4 ) / 10;
	fidxmax =  3.8 - ( yAt + 4 ) / 10;
      }

      effxyHisto->fill( xAt, yAt );
      effvsxy->fill( xAt, yAt, nm ); // CMS DUT efficiency profile

      if( yAt > fidymin && yAt < fidymax ) {
	effvsx->fill( xAt, nm ); // CMS DUT efficiency profile
	if( probchi > 0.01 ) effvsxg->fill( xAt, nm );
      }

      if( xAt > fidxmin && xAt < fidxmax ) {
	effvsy->fill( yAt, nm ); // CMS DUT efficiency profile
      }

      if( xAt > fidxmin && xAt < fidxmax && yAt > fidymin && yAt < fidymax ) {
	eff300->fill( (time_now_tlu-time_event0)/fTLU, nm );
	eff600->fill( (time_now_tlu-time_event0)/fTLU, nm );
	eff1200->fill( (time_now_tlu-time_event0)/fTLU, nm );
	eff1800->fill( (time_now_tlu-time_event0)/fTLU, nm );
	eff3600->fill( (time_now_tlu-time_event0)/fTLU, nm );

	effvsxmym->fill( xmod, ymod, nm ); // CMS DUT efficiency profile
	effvsndri->fill(downstream_triplets->size(), nm); // CMS DUT efficiency vs n drownstream tracks	    
	
	// KIT: added for efficiency analysis
	double dotsize=10;
	double cutsize=5;

	if( xmod >= 50-cutsize  && xmod <= 50+cutsize )  effvsxm50->fill( ymod, nm );
	if( xmod >= 100-dotsize && xmod <= 100+dotsize ) effvsxm100->fill( ymod, nm );
	if( xmod >= 150-cutsize && xmod <= 150+cutsize ) effvsxm150->fill( ymod, nm );
	if( xmod >= 200-dotsize && xmod <= 200+dotsize ) effvsxm200->fill( ymod, nm );
	if( xmod >= 250-cutsize && xmod <= 250+cutsize ) effvsxm250->fill( ymod, nm );

	if( ymod >= 25-cutsize  && ymod <= 25+cutsize )  effvsym25->fill( xmod, nm );
	if( ymod >= 50-dotsize  && ymod <= 50+dotsize )  effvsym50->fill( xmod, nm );
	if( ymod >= 75-cutsize  && ymod <= 75+cutsize )  effvsym75->fill( xmod, nm );
	if( ymod >= 100-cutsize && ymod <= 100+cutsize ) effvsym100->fill( xmod, nm );
	if( ymod >= 125-cutsize && ymod <= 125+cutsize ) effvsym125->fill( xmod, nm );
	if( ymod >= 150-dotsize && ymod <= 150+dotsize ) effvsym150->fill( xmod, nm );
	if( ymod >= 175-cutsize && ymod <= 175+cutsize ) effvsym175->fill( xmod, nm );
	// KIT end

	//if( prevdutrefddt == 0 && dutrefddt == 0 )
	//effd600->fill( (time_now_tlu-time_event0)/fTLU, nm );
	if( downstream_triplets->size() < 3 ) effn600->fill( (time_now_tlu-time_event0)/fTLU, nm );
	if( downstream_triplets->size() < 2 ) effm600->fill( (time_now_tlu-time_event0)/fTLU, nm );

      }

    }// triplet-driplet match

     //------------------------------------------------------------------------
     // eff(REF) with DUT as timing plane:

    if( fabs(dx) < 0.1 && fabs(dy) < 0.1 && trip.linked_dut ) {

      bool nm = 0;
      if( drip.linked_dut ) nm = 1;

      rffvsxy->fill( -xR, -yR, nm ); // CMS REF efficiency profile

      if( fabs( yR ) < 3 ) {
	rffvsx->fill( -xR, nm ); // CMS REF efficiency profile
      }

    }// triplet-driplet match

     //------------------------------------------------------------------------
     // intersect point in z:
    hit intersect = (*tr).intersect();
    double zx = intersect.x;
    double zy = intersect.y;

    if( fabs(dy) < 0.1 ) { // no cut on dx
      if( fabs( kx ) > 0.003 ) { // 
	sixzx3Histo->fill( zx - _planePosition[2] );
      }
      if( fabs( kx ) > 0.002 ) { // 
	sixzx2Histo->fill( zx - _planePosition[2] );
      }
    }

    if( fabs(dx) < 0.1 ) { // no cut on dy
      if( fabs( ky ) > 0.003 ) { // 
	sixzy3Histo->fill( zy - _planePosition[2] );
      }
      if( fabs( ky ) > 0.002 ) { // 
	sixzy2Histo->fill( zy - _planePosition[2] );
      }
    }

    //------------------------------------------------------------------------
    // z intersect:

    if( fabs(dx) < 0.2 && fabs(dy) < 0.2 ) { // looser cut allows more z range

      // measure scattering angle x after cuts in y:
      // cut on ky creates bias in kx

      if( fabs( kx ) > 0.001 ) {
	sixzx1Histo->fill( zx - _planePosition[2] );
	if( fabs( zx - DUTz ) < 30 ) {
	  sixkyzxHisto->fill( ky*1E3 );
	  sixkxzxHisto->fill( kx*1E3 ); // plot with gap, fittp0g.C("sixkxzx")
	}
      }

      if( fabs( ky ) > 0.001 ) {
	sixzy1Histo->fill( zy - _planePosition[2] );
	if( fabs( zy - DUTz ) < 30 ) {
	  sixkxzyHisto->fill( kx*1E3 );
	  sixkyzyHisto->fill( ky*1E3 ); // plot with gap
	}
      }

    }//match

     //------------------------------------------------------------------------

  } // Loop over found tracks

  nsixHisto->fill( telescope_tracks->size() );

  //prevdutrefddt = dutrefddt;
  // Clear memory
  delete telescope_tracks;
  delete downstream_triplets;
  delete upstream_triplets;
  delete hits;
}


//------------------------------------------------------------------------------
void EUTelAnalysisCMSPixel::check( LCEvent * /* evt */  ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


//------------------------------------------------------------------------------
void EUTelAnalysisCMSPixel::end(){

  // Print the summary:
  // ...but don't print it twice for CMSPixelDecoder: muting it:
  CMSPixel::Log::ReportingLevel() = CMSPixel::Log::FromString("QUIET");

  streamlog_out(MESSAGE5)
    << "---------------------------------------------------------------------------------------------------------" << std::endl
    << std::endl
    << "Processed events:    "
    << std::setw(10) << std::setiosflags(std::ios::right)
    << _nEvt << std::resetiosflags(std::ios::right) << std::endl
    << "DUT clusters " << n_clusters_dut << std::endl
    << "DUT matches  " << n_matched_clusters_dut << " (" << ((double)n_matched_clusters_dut/n_clusters_dut)*100 << "%)" << std::endl << std::endl
    << "REF clusters " << n_clusters_ref << std::endl
    << "REF matches  " << n_matched_clusters_ref << " (" << ((double)n_matched_clusters_ref/n_clusters_ref)*100 << "%)" << std::endl
    << std::endl;

  // pre-align:
  streamlog_out(MESSAGE5) << "pre-align:" << std::endl;

  streamlog_out(MESSAGE5) << "DUT x   " << cmsdxaHisto->title() << std::endl;
  streamlog_out(MESSAGE5) << "Entries " << cmsdxaHisto->entries() << std::endl;
  streamlog_out(MESSAGE5) << "Maximum " << cmsdxaHisto->maxBinHeight() << std::endl;

  int nmax = 0;
  double DUTx = 0;
  for( int ib = 0; ib < cmsdxaHisto->axis().bins(); ++ib ) {
    if( cmsdxaHisto->binHeight(ib) > nmax ) {
      nmax = static_cast<int>(cmsdxaHisto->binHeight(ib));
      DUTx = cmsdxaHisto->binMean(ib);
    }
  }
  streamlog_out(MESSAGE5) << "DUTx    " << DUTx << std::endl;

  streamlog_out(MESSAGE5) << "DUT y   " << cmssyaHisto->title() << std::endl;
  streamlog_out(MESSAGE5) << "Entries " << cmssyaHisto->entries() << std::endl;
  streamlog_out(MESSAGE5) << "Maximum " << cmssyaHisto->maxBinHeight() << std::endl;

  nmax = 0;
  double DUTy = 0;
  for( int ib = 0; ib < cmssyaHisto->axis().bins(); ++ib ) {
    if( cmssyaHisto->binHeight(ib) > nmax ) {
      nmax = static_cast<int>(cmssyaHisto->binHeight(ib));
      DUTy = cmssyaHisto->binMean(ib);
    }
  }
  streamlog_out(MESSAGE5) << "DUTy    " << DUTy << std::endl;

  streamlog_out(MESSAGE5) << "REF x   " << refdxaHisto->title() << std::endl;
  streamlog_out(MESSAGE5) << "Entries " << refdxaHisto->entries() << std::endl;
  streamlog_out(MESSAGE5) << "Maximum " << refdxaHisto->maxBinHeight() << std::endl;

  nmax = 0;
  double REFx = 0;
  for( int ib = 0; ib < refdxaHisto->axis().bins(); ++ib ) {
    if( refdxaHisto->binHeight(ib) > nmax ) {
      nmax = static_cast<int>(refdxaHisto->binHeight(ib));
      REFx = refdxaHisto->binMean(ib);
    }
  }
  streamlog_out(MESSAGE5) << "REFx    " << REFx << std::endl;

  streamlog_out(MESSAGE5) << "REF y   " << refsyaHisto->title() << std::endl;
  streamlog_out(MESSAGE5) << "Entries " << refsyaHisto->entries() << std::endl;
  streamlog_out(MESSAGE5) << "Maximum " << refsyaHisto->maxBinHeight() << std::endl;

  nmax = 0;
  double REFy = 0;
  for( int ib = 0; ib < refsyaHisto->axis().bins(); ++ib ) {
    if( refsyaHisto->binHeight(ib) > nmax ) {
      nmax = static_cast<int>(refsyaHisto->binHeight(ib));
      REFy = refsyaHisto->binMean(ib);
    }
  }
  streamlog_out(MESSAGE5) << "REFy    " << REFy << std::endl;

  char cCurrentPath[150];
  int nRunNew;
  
  if (!getcwd(cCurrentPath, sizeof(cCurrentPath))){
    printf("No working directory found.");
    nRunNew = _nRun;
  }else{
    streamlog_out(MESSAGE5) << "Working in directory: " << cCurrentPath << endl;
    string path(cCurrentPath);
    istringstream(path.substr(path.find_last_of("/run")+1,6)) >> nRunNew;
    if(nRunNew!=_nRun)streamlog_out(MESSAGE5) << "Overwriting runnumber " << _nRun << " with " << nRunNew << endl;
  }


  streamlog_out(MESSAGE5) 
    << std::endl
    << "runlistPreAlign: "
    << nRunNew
    << "," << _alignmentrun
    << "," << _gearfile
    << "," <<_eBeam
    << "," << _DUT_chip
    << "," << _DUT_gain
    << "," << _DUT_calibration_type
    << "," << _DUT_conversion
    << "," << _REF_chip
    << "," << _REF_gain
    << "," << _REF_calibration_type
    << "," << DUTx
    << "," << DUTy
    << "," << DUTz - _planePosition[2]
    << "," << tilt
    << "," << turn
    << "," << DUTrot
    << "," << REFx
    << "," << REFy
    << "," << _REFz
    << "," << _REFrot
    << "," << _skip_dut
    << "," << _skip_ref
    << "," << _skip_tel
    << std::endl;
  streamlog_out(MESSAGE5) << std::endl;


  ofstream prealignrunfile;
  prealignrunfile.open("prelines-for-runlist.txt",ios::app);
  prealignrunfile << nRunNew << "," << _alignmentrun << "," << _gearfile << "," <<_eBeam << "," << _DUT_chip << "," << _DUT_gain  << "," << _DUT_calibration_type << "," << _DUT_conversion << "," << _REF_chip << "," << _REF_gain << "," << _REF_calibration_type << "," << DUTx << "," << DUTy << "," << DUTz - _planePosition[2] << "," << tilt << "," << turn << "," << DUTrot << "," << REFx << "," << REFy << "," << _REFz << "," << _REFrot << "," << _skip_dut << "," << _skip_ref << "," << _skip_tel << endl;
  prealignrunfile.close();

  // Clean memory:

  delete [] _planeSort;
  delete [] _planeID;
  delete [] _planeThickness;
  delete [] _planeX0;
  delete [] _planeResolution;

  // Pede:

  // close the output file:
  delete mille;

  streamlog_out( MESSAGE4 ) << std::endl << "Generating the steering file for the pede program..." << std::endl;

  streamlog_out( MESSAGE2 ) << "have " << ngbl << " GBL tracks" << std::endl;
  streamlog_out( MESSAGE2 ) << "have " << naldut << " DUT links" << std::endl;

  bool ldut = naldut > 99;

  ofstream steerFile;
  std::stringstream steername;
  steername << "run" << _nRun << "-steerPede.txt";
  steerFile.open( steername.str().c_str() );

  if( steerFile.is_open() && ldut ) {

    steerFile << "! generated by EUTelTestFitter" << std::endl;
    steerFile << "Cfiles" << std::endl;
    steerFile << m_millefilename << std::endl;
    steerFile << std::endl;

    //FIXME maybe reduce dimensions if tilt || turn == 0?
    // set the second value to -1.0 to exclude
    steerFile << "Parameter" << std::endl;
    steerFile << 1 << "  0.0  0.0" << std::endl; // dx
    steerFile << 2 << "  0.0  0.0" << std::endl; // dy
    steerFile << 3 << "  0.0  0.0" << std::endl; // rot
    steerFile << 4 << "  0.0  0.0" << std::endl; // dtilt
    steerFile << 5 << "  0.0  0.0" << std::endl; // dturn
    steerFile << 6 << "  0.0  0.0" << std::endl; // dz
    steerFile << std::endl;
    steerFile << "! chiscut 5.0 2.5" << std::endl;
    steerFile << "outlierdownweighting 4" << std::endl;
    steerFile << std::endl;
    steerFile << "method inversion 10  0.1" << std::endl;
    steerFile << "threads 10 1" << std::endl;
    steerFile << std::endl;
    steerFile << "! histprint" << std::endl;
    steerFile << std::endl;
    steerFile << "end" << std::endl;

    steerFile.close();
    streamlog_out( MESSAGE2 ) << "Pede steer file written." << std::endl;

    // before starting pede, let's check if it is in the path
    bool isPedeInPath = true;
    // create a new process
    redi::ipstream which("which pede");
    // wait for the process to finish
    which.close();

    // get the status
    // if it is 255 then the program wasn't found in the path
    isPedeInPath = !( which.rdbuf()->status() == 255 );

    if( !isPedeInPath ) {
      streamlog_out( ERROR ) << "Pede cannot be executed because not found in the path" << std::endl;
    }
    else {
      std::string command = "pede  " + steername.str();

      streamlog_out( MESSAGE2 ) << std::endl;
      streamlog_out( MESSAGE2 ) << "Starting pede..." << std::endl;
      streamlog_out( MESSAGE2 ) << command.c_str() << std::endl;

      // run pede and create a streambuf that reads its stdout and stderr
      redi::ipstream pede( command.c_str(), redi::pstreams::pstdout|redi::pstreams::pstderr ); 
      std::string output;
      while ( getline( pede, output ) ) {
	streamlog_out( MESSAGE2 ) << output << std::endl;
      }

      // wait for the pede execution to finish
      pede.close();
      // check the exit value of pede
      if( pede.rdbuf()->status() == 0 )
	streamlog_out( MESSAGE2 ) << "Pede successfully finished" << std::endl;

      // reading back the millepede.res file:
      std::string millepedeResFileName = "millepede.res";
      streamlog_out( MESSAGE2 ) << "Reading back the " << millepedeResFileName << std::endl;

      // open the millepede ASCII output file
      ifstream millepede( millepedeResFileName.c_str() );

      if( millepede.bad() || !millepede.is_open() ) {
	streamlog_out( ERROR4 )
	  << "Error opening the " << millepedeResFileName << std::endl
	  << "The alignment file cannot be saved" << std::endl;
      }
      else {
	std::vector<double > tokens;
	std::stringstream tokenizer;
	std::string line;
	double buffer;

	// get the first line and throw it away since it is a comment!
	getline( millepede, line );
	streamlog_out(DEBUG2) << "line: " <<  line  << std::endl;

	unsigned int numpars = 6; // DUT dx, dy, drot, dtilt, dturn, dz
	std::map< unsigned int, double > alpar; // map = associative array

	for( unsigned int ipar = 0; ipar < numpars; ++ipar ) {

	  if( millepede.eof() ) break;
	  getline( millepede, line );
	  if( line.empty() ) continue;

	  tokens.clear();
	  tokenizer.clear();
	  tokenizer.str( line );

	  while( tokenizer >> buffer ) tokens.push_back( buffer );
	  int lpar = static_cast< int >( tokens[0] + 0.5 ); // par label
	  bool isFixed = ( tokens.size() == 3 );

	  if( isFixed ) {
	    streamlog_out(MESSAGE5) << "Parameter " << lpar
				    << std::resetiosflags(std::ios::floatfield) << std::setprecision(9) 
				    << " is at " << tokens[1]
				    << " (fixed)"  << std::endl;
	  }
	  else {
	    streamlog_out(MESSAGE5) << "Parameter " << lpar
				    << std::resetiosflags(std::ios::floatfield) << std::setprecision(9) 
				    << " is at " << tokens[1]
				    << " +/- " << tokens[4] << std::endl;
	  }

	  alpar[lpar] = tokens[1];

	}//loop param

	if( ldut ){
	  streamlog_out(MESSAGE5) << std::endl << "DUT alignment corrections:" << std::endl;
	  streamlog_out(MESSAGE5) << std::resetiosflags(std::ios::floatfield) << std::setprecision(9) << "   dx    " << alpar[1]*1E3 << " um" << std::endl;
	  streamlog_out(MESSAGE5) << "   dy        = " << alpar[2]*1E3 << " um" << std::endl;
	  streamlog_out(MESSAGE5) << "   drot      = " << alpar[3]*1E3 << " mrad" << std::endl;
	  streamlog_out(MESSAGE5) << "   dtilt     = " << alpar[4]*180/3.141592654 << " deg" << std::endl;
	  streamlog_out(MESSAGE5) << "   dturn     = " << alpar[5]*180/3.141592654 << " deg" << std::endl;
	  streamlog_out(MESSAGE5) << "   dz        = " << alpar[6]*1E3 << " um" << std::endl << std::endl;
	  streamlog_out(MESSAGE5) << "DUT alignment:" << std::endl;
	  streamlog_out(MESSAGE5) << "   DUTalignx = " << DUTalignx-alpar[1] << ";" << std::endl;
	  streamlog_out(MESSAGE5) << "   DUTaligny = " << DUTaligny-alpar[2] << ";" << std::endl;
	  streamlog_out(MESSAGE5) << "   DUTrot    = " << DUTrot-alpar[3] << ";" << std::endl;
	  streamlog_out(MESSAGE5) << "   DUTtilt   = " << tilt-alpar[4]*180/3.141592654 << ";" << std::endl;
	  streamlog_out(MESSAGE5) << "   DUTturn   = " << turn-alpar[5]*180/3.141592654 << ";" << std::endl;
	  streamlog_out(MESSAGE5) << "   DUTz      = " << DUTz-alpar[6] - _planePosition[2] << " + _planePosition[2];" << std::endl;

	  // If the tilt correction is small (< 1mrad) we should be stable, correct the Landau peak:
	  // tilt angle correction below 0.2deg:
	  if(fabs(alpar[4]) < 0.003490659) {
	    streamlog_out(MESSAGE5) << std::endl << "Tilt angle converged." << std::endl;
	    // Update the DUT_conversion factor to get the landau peak
	    // to 22k electrons:
	    streamlog_out(MESSAGE5) << std::endl << "DUT charge calibration corrections:" << std::endl;
	    Double_t landau_peak = landau_gauss_peak(cmsq0fHistoRoot);
	    streamlog_out(MESSAGE5) << "   current landau peak = " << landau_peak << std::endl;
	    streamlog_out(MESSAGE5) << "   old conversion      = " << _DUT_conversion << std::endl;
	    _DUT_conversion = _DUT_conversion*22/landau_peak;
	    streamlog_out(MESSAGE5) << "   new conversion      = " << _DUT_conversion << std::endl;
	  }
	  else {
	    streamlog_out(MESSAGE5) << "Tilt angle not converged, no charge calibration correction applied." << std::endl;
	  }

	  streamlog_out(MESSAGE5)
	    << std::endl
	    << "for runlist.csv:" << std::endl
	    << "runlistFullAlign: "
	    << nRunNew
	    << "," << _alignmentrun
	    << "," << _gearfile
	    << std::setprecision(2) << "," << _eBeam
	    << "," << _DUT_chip
	    << "," << _DUT_gain
	    << "," << _DUT_calibration_type
	     << std::setprecision(3) << "," << _DUT_conversion
	    << "," << _REF_chip
	    << "," << _REF_gain
	    << "," << _REF_calibration_type
	    << std::setprecision(9) << "," << DUTalignx-alpar[1]
	    << "," << DUTaligny-alpar[2]
	    << "," << DUTz-alpar[6] - _planePosition[2]
	    << "," << tilt-alpar[4]*180/3.141592654
	    << "," << turn-alpar[5]*180/3.141592654
	    << "," << DUTrot-alpar[3]
	    << "," << _REFalignx
	    << "," << _REFaligny
	    << "," << _REFz
	    << "," << _REFrot
	    << "," << _skip_dut 
	    << "," << _skip_ref 
	    << "," << _skip_tel 
	    << std::endl;


	  
	  ofstream runfile;
	  runfile.open("lines-for-runlist.txt",ios::app);
	  runfile << nRunNew << "," << _alignmentrun << "," << _gearfile << "," << _eBeam << "," << _DUT_chip << "," << _DUT_gain << "," << _DUT_calibration_type << "," << _DUT_conversion << "," << _REF_chip << "," << _REF_gain << "," << _REF_calibration_type << "," << DUTalignx-alpar[1] << "," << DUTaligny-alpar[2] << "," << DUTz-alpar[6] - _planePosition[2] << "," << tilt-alpar[4]*180/3.141592654 << "," << turn-alpar[5]*180/3.141592654 << "," << DUTrot-alpar[3] << "," << _REFalignx << "," << _REFaligny << "," << _REFz << "," << _REFrot << "," << _skip_dut << "," << _skip_ref << "," << _skip_tel << endl;
	  runfile.close();

	} // ldut

      }//millepede OK

      millepede.close();

    }//PedeInPath

  }// pede steer file open

  else {
    streamlog_out( ERROR2 ) << "Could not open MillepedeII-steering file." << std::endl;
    //throw marlin::ParseException("Could not open MillepedeII-steering file.");
  }

  delete [] _planePosition;

} // end end


double GetCorrelationFactor(AIDA::IHistogram2D* histo)
{
  const int binsX = histo->xAxis().bins();
  const int binsY = histo->yAxis().bins();

  double sumw = 0;
  double sumwx = 0;
  double sumwxx= 0;
  double sumwy = 0;
  double sumwyy = 0;
  double sumwxy = 0;

  for(int y = 0; y < binsY; ++y)
    {
      const double yc = (histo->yAxis().binUpperEdge(y) + histo->yAxis().binLowerEdge(y))/2;
      for(int x = 0; x < binsX; ++x)
	{
	  const double xc = (histo->xAxis().binUpperEdge(x) + histo->xAxis().binLowerEdge(x))/2;
	  const double w = histo->binHeight(x, y);

	  sumw += w;
	  sumwx += xc*w;
	  sumwxx += xc*xc*w;
	  sumwy += yc*w;
	  sumwyy += yc*yc*w;
	  sumwxy += xc*yc*w;
	}
    }

  const double rmsX = sqrt(fabs(sumwxx/sumw - sumwx/sumw * sumwx/sumw));
  const double rmsY = sqrt(fabs(sumwyy/sumw - sumwy/sumw * sumwy/sumw));

  const double cov = sumwxy/sumw - sumwx/sumw*sumwy/sumw;
  return cov/rmsX/rmsY;
}

std::string EUTelAnalysisCMSPixel::ZeroPadNumber(int num, int len)
{
  std::ostringstream ss;
  ss << std::setw( len ) << std::setfill( '0' ) << num;
  return ss.str();
}

uint8_t EUTelAnalysisCMSPixel::GetChipTypeFromID(int chip_id) {

  // This function implements the DESY chip ID numbering scheme:
  if(chip_id < 100) return CMSPixel::ROC_PSI46V2;
  else if(chip_id < 200) return CMSPixel::ROC_PSI46XDB;
  else if(chip_id < 300) return CMSPixel::ROC_PSI46DIG;
  else return CMSPixel::ROC_PSI46DIGV2;
}


TMatrixD EUTelAnalysisCMSPixel::JacobianPointToPoint( double ds ) {
  /* for GBL:
     Jacobian for straight line track
     track = q/p, x', y', x, y
     0,   1,  2,  3, 4
  */
  TMatrixD jac( 5, 5 );
  jac.UnitMatrix();
  jac[3][1] = ds; // x = x0 + xp * ds
  jac[4][2] = ds; // y = y0 + yp * ds
  return jac;
}


void EUTelAnalysisCMSPixel::TelescopeCorrelationPlots(std::vector<hit> * telescopehits) {
  for( std::vector<hit>::iterator ihit = telescopehits->begin(); ihit != telescopehits->end(); ihit++ ){

    int ipl = (*ihit).plane;

    for( std::vector<hit>::iterator jhit = telescopehits->begin(); jhit != telescopehits->end(); jhit++ ){

      int jpl = (*jhit).plane;
      double dx = (*jhit).x - (*ihit).x;
      double dy = (*jhit).y - (*ihit).y;

      if( ipl == 0 ) {
	if( jpl == 1 ) {
	  dx01Histo->fill( dx );
	  dy01Histo->fill( dy );
	  if( abs(dy) < 1 ) du01Histo->fill( dx );
	}
	if( jpl == 2 ) {
	  dx02Histo->fill( dx );
	}
	if( jpl == 3 ) {
	  dx03Histo->fill( dx );
	}
	if( jpl == 4 ) {
	  dx04Histo->fill( dx );
	}
	if( jpl == 5 ) {
	  dx05Histo->fill( dx );
	}

      }//ipl==0

      if( ipl == 1 ) {
	if( jpl == 2 ) {
	  dx12Histo->fill( dx );
	  dy12Histo->fill( dy );
	  if( abs(dy) < 1 ) du12Histo->fill( dx );
	}
      }

      if( ipl == 2 ) {
	if( jpl == 3 ) {
	  dx23Histo->fill( dx );
	  dy23Histo->fill( dy );
	  if( abs(dy) < 1 ) du23Histo->fill( dx );
	}
      }

      if( ipl == 3 ) {
	if( jpl == 4 ) {
	  dx34Histo->fill( dx );
	  dy34Histo->fill( dy );
	  if( abs(dy) < 1 ) du34Histo->fill( dx );
	}
      }

      if( ipl == 4 ) {
	if( jpl == 5 ) {
	  dx45Histo->fill( dx );
	  dy45Histo->fill( dy );
	  if( abs(dy) < 1 ) du45Histo->fill( dx );
	}
      }
    }
  }
}


std::vector<EUTelAnalysisCMSPixel::triplet> * EUTelAnalysisCMSPixel::FindTriplets(std::vector<hit> * hits, unsigned int plane0, unsigned int plane1, unsigned int plane2) {

  std::vector<triplet> * triplets = new std::vector<triplet>;

  // Cut on the triplet track angle: +- 10 mrad
  double triplet_angle_cut = 0.010;
  // Cut on the triplet residual on the middle plane:
  double triplet_residual_cut = 0.1; // [mm]
  
  for( std::vector<hit>::iterator ihit = hits->begin(); ihit != hits->end(); ihit++ ){
    if( (*ihit).plane != plane0 ) continue; // First plane

    for( std::vector<hit>::iterator jhit = hits->begin(); jhit != hits->end(); jhit++ ){
      if( (*jhit).plane != plane2 ) continue; // Last plane

      for( std::vector<hit>::iterator khit = hits->begin(); khit != hits->end(); khit++ ){
	if( (*khit).plane != plane1 ) continue; // Middle plane

	// Create new preliminary triplet from the three hits:
	triplet new_triplet((*ihit),(*khit),(*jhit));

	// Setting cuts on the triplet track angle:
	if( fabs(new_triplet.getdy()) > triplet_angle_cut * new_triplet.getdz()) continue;
	if( fabs(new_triplet.getdx()) > triplet_angle_cut * new_triplet.getdz()) continue;

	// Setting cuts on the triplet residual on the middle plane
	if( fabs(new_triplet.getdx(plane1)) > triplet_residual_cut) continue;
	if( fabs(new_triplet.getdy(plane1)) > triplet_residual_cut) continue;

	// The triplet is accepted, push it back:
	triplets->push_back(new_triplet);
	streamlog_out(DEBUG2) << new_triplet;

      }//loop over hits
    }//loop over hits
  }// loop over hits

  return triplets;
}


std::vector<EUTelAnalysisCMSPixel::track> * EUTelAnalysisCMSPixel::MatchTriplets(std::vector<triplet> * up, std::vector<triplet> * down, double z_match) {
  
  std::vector<track> * tracks = new std::vector<track>;

  // Cut on the matching of two triplets [mm]
  double intersect_residual_cut = 0.1;

  for( std::vector<triplet>::iterator trip = up->begin(); trip != up->end(); trip++ ){

    for( std::vector<triplet>::iterator drip = down->begin(); drip != down->end(); drip++ ){

      // Build a track candidate from one upstream and one downstream triplet:
      track newtrack((*trip),(*drip));

      // Track kinks as difference in triplet slopes:
      double kx = newtrack.kink_x();
      double ky = newtrack.kink_y();

      // Track impact position at Matching Point from Downstream:
      double xB = (*drip).getx_at(z_match);
      double yB = (*drip).gety_at(z_match);

      // Track impact position at Matching Point from Upstream:
      double xA = (*trip).getx_at(z_match);
      double yA = (*trip).gety_at(z_match);
      
      double dx = xB - xA; // driplet - triplet
      double dy = yB - yA;

      sixkxHisto->fill( kx*1E3 );
      sixkyHisto->fill( ky*1E3 );
      sixdxHisto->fill( dx );
      sixdyHisto->fill( dy );
      
      
      if( fabs(dy) < 0.4 ) sixdxcHisto->fill( dx*1E3 ); // sig = 17 um at 5 GeV
      if( fabs(dx) < 0.4 ) sixdycHisto->fill( dy*1E3 );

      // match driplet and triplet:
      if( fabs(dx) > intersect_residual_cut) continue;
      if( fabs(dy) > intersect_residual_cut) continue;
	
      sixkxcHisto->fill( kx*1E3 );
      sixkycHisto->fill( ky*1E3 );
      sixxHisto->fill( -xA ); // -xA = x_DP = out
      sixyHisto->fill( -yA ); // -yA = y_DP = up
      sixxyHisto->fill( -xA, -yA ); // DP: x_out, y_up

      // Fill kink map histogram:
      if( fabs( kx ) > 0.002 || fabs( ky ) > 0.002 ) sixxycHisto->fill( -xA, -yA );

      kinkvsxy->fill( -xA, -yA, (kx*kx + ky*ky)*1E6 ); //<kink^2> [mrad^2]

      // Add the track to the vector, the triplets are matched:
      tracks->push_back(newtrack);

    } // Downstream
  } // Upstream

  streamlog_out(DEBUG2) << "Found " << tracks->size() << " tracks from matched triplets." << std::endl;
  return tracks;
}



EUTelAnalysisCMSPixel::track::track(triplet up, triplet down) : upstream(), downstream() {
  upstream = up;
  downstream = down;
}

double EUTelAnalysisCMSPixel::track::kink_x() {
  return (downstream.slope().x - upstream.slope().x);
}

double EUTelAnalysisCMSPixel::track::kink_y() {
  return (downstream.slope().y - upstream.slope().y);
}

EUTelAnalysisCMSPixel::hit EUTelAnalysisCMSPixel::track::intersect() {
  hit inter;
  // Re-check what this actually is...
  // and simplifie using triplet class members...
  inter.x = ( upstream.base().x - upstream.slope().x * upstream.base().z - downstream.base().x + downstream.slope().x * downstream.base().z ) / kink_x();
  inter.y = ( upstream.base().y - upstream.slope().y * upstream.base().z - downstream.base().y + downstream.slope().y * downstream.base().z ) / kink_y();
  return inter;
}

EUTelAnalysisCMSPixel::triplet EUTelAnalysisCMSPixel::track::get_upstream() {
  return upstream;
}

EUTelAnalysisCMSPixel::triplet EUTelAnalysisCMSPixel::track::get_downstream() {
  return downstream;
}

EUTelAnalysisCMSPixel::hit EUTelAnalysisCMSPixel::track::gethit(int plane) {
  if(plane < 3) return upstream.gethit(plane);
  else return downstream.gethit(plane);
}



EUTelAnalysisCMSPixel::triplet::triplet() : linked_dut(false), cmsdx(0), cmsdy(0), hits() {
  // Empty default constructor
}

EUTelAnalysisCMSPixel::triplet::triplet(hit hit0, hit hit1, hit hit2) : linked_dut(false), cmsdx(0), cmsdy(0), hits() {
  triplet();
  filltriplet(hit0, hit1, hit2);
}

EUTelAnalysisCMSPixel::hit EUTelAnalysisCMSPixel::triplet::getpoint_at(double z) {
  hit impact;
  impact.z = z - base().z;
  impact.x = base().x + slope().x * impact.z;
  impact.y = base().y + slope().y * impact.z;
  return impact;
}

double EUTelAnalysisCMSPixel::triplet::getx_at(double z) {
  return this->base().x + this->slope().x * (z - this->base().z);
}

double EUTelAnalysisCMSPixel::triplet::getdx() {
  return this->hits.rbegin()->second.x - this->hits.begin()->second.x;
}

double EUTelAnalysisCMSPixel::triplet::getdx(int ipl) {
  return this->hits[ipl].x - this->base().x - this->slope().x * (this->hits[ipl].z - this->base().z);
}

double EUTelAnalysisCMSPixel::triplet::getdx(hit point) {
  return point.x - this->base().x - this->slope().x * (point.z - this->base().z);
}

double EUTelAnalysisCMSPixel::triplet::gety_at(double z) {
  return this->base().y + this->slope().y * (z - this->base().z);
}

double EUTelAnalysisCMSPixel::triplet::getdy() {
  return this->hits.rbegin()->second.y - this->hits.begin()->second.y;
}

double EUTelAnalysisCMSPixel::triplet::getdy(int ipl) {
  return this->hits[ipl].y - this->base().y - this->slope().y * (this->hits[ipl].z - this->base().z);
}

double EUTelAnalysisCMSPixel::triplet::getdy(hit point) {
  return point.y - this->base().y - this->slope().y * (point.z - this->base().z);
}

double EUTelAnalysisCMSPixel::triplet::getdz() {
  return this->hits.rbegin()->second.z - this->hits.begin()->second.z;
}

EUTelAnalysisCMSPixel::hit EUTelAnalysisCMSPixel::triplet::gethit(int plane) {
  return this->hits[plane];
}

EUTelAnalysisCMSPixel::hit EUTelAnalysisCMSPixel::triplet::base() {
  hit center;
  center.x = 0.5*( this->hits.begin()->second.x + this->hits.rbegin()->second.x );
  center.y = 0.5*( this->hits.begin()->second.y + this->hits.rbegin()->second.y );
  center.z = 0.5*( this->hits.begin()->second.z + this->hits.rbegin()->second.z );
  return center;
}

EUTelAnalysisCMSPixel::hit EUTelAnalysisCMSPixel::triplet::slope() {
  hit sl;
  double dz = (this->hits.rbegin()->second.z - this->hits.begin()->second.z);
  sl.x = (this->hits.rbegin()->second.x - this->hits.begin()->second.x) / dz;
  sl.y = (this->hits.rbegin()->second.y - this->hits.begin()->second.y) / dz;
  return sl;
}


void EUTelAnalysisCMSPixel::FillDeltaTPlots(double time_now_tlu, double time_before, double time_event0) {

  //Filling event timing histograms for the telescope:
  t1Histo->fill( ( time_now_tlu - time_event0 ) / fTLU ); //event time
  t10Histo->fill( ( time_now_tlu - time_event0 ) / fTLU ); //event time
  t100Histo->fill( ( time_now_tlu - time_event0 ) / fTLU ); //event time
  t300Histo->fill( ( time_now_tlu - time_event0 ) / fTLU ); //event time
  t600Histo->fill( ( time_now_tlu - time_event0 ) / fTLU ); //event time
  t1000Histo->fill( ( time_now_tlu - time_event0 ) / fTLU ); //event time
  t1800Histo->fill( ( time_now_tlu - time_event0 ) / fTLU ); //event time
  t3600Histo->fill( ( time_now_tlu - time_event0 ) / fTLU ); //event time
  
  // dt plots:
  dtHisto->fill( ( time_now_tlu - time_before ) / 384E0 ); //us
  dtmsHisto->fill( ( time_now_tlu - time_before ) / 384E3 ); //ms
  logdtHisto->fill( std::log( ( time_now_tlu - time_before ) / 384E3 ) / std::log(10.0) );
  
  //FIXME fill this only if cluster exist in the DUT!
  logdtcmsHisto->fill( std::log( ( time_now_tlu - time_before ) / gTLU/1E6 ) / std::log(10.0) );
  // scan for Umlauftakt:

  for( double tau = 372; tau < 378; tau += 0.01 ) {
    double ganz;
    double frac = std::modf( (time_now_tlu - time_before) / tau, &ganz );
    dtfvstau->fill( tau, fabs( frac - 0.5 ) ); // peak at 375.055 for run 5513
    frac = std::modf( (time_now_tlu - time_event0) / tau, &ganz ); // absolute phase
    tfvstau->fill( tau, fabs( frac - 0.5 ) );
  }

  double tau = 375.06106; // Umlauftakt [TLU ticks] for run 5359

  double ganz;
  double frac = std::modf( (time_now_tlu - time_before) / tau, &ganz );
  if( frac > 0.5 ) frac = frac - 1; // complement for peak around zero
  dtfHisto->fill( frac ); // phase [-0.5,0.5]
  dtfvsdt->fill( ( time_now_tlu - time_before ) / gTLU/1E3, frac ); // new_tau = old*(1+slope)
  dtfvst->fill( ( time_now_tlu - time_event0 ) / fTLU, frac ); // phase [-0.5,0.5]

  // absolute phase:
  frac = std::modf( (time_now_tlu - time_event0) / tau, &ganz ); // ganz = number of turns
  if( frac > 0.5 ) frac = frac - 1; // complement
  tfHisto->fill( frac ); // phase [-0.5,0.5]
  tfvst->fill( ganz, frac ); // phase [-0.5,0.5]
  tfvst1->fill( ganz, frac ); // phase [-0.5,0.5]
  tfvst10->fill( ganz, frac ); // phase [-0.5,0.5]
  tfvst100->fill( ganz, frac ); // phase [-0.5,0.5]
  tfvst300->fill( ganz, frac ); // phase [-0.5,0.5]
  
}

void EUTelAnalysisCMSPixel::FillClusterStatisticsPlots(std::vector<cluster> dutclusters, int dutpix, std::vector<cluster> refclusters, int refpix) {

  // DUT cluster statistics:
  dutnclusHisto->fill( dutclusters.size() );
  if( dutpix > 0 ) {
    for( std::vector<cluster>::iterator c = dutclusters.begin(); c != dutclusters.end(); c++ ){
      n_clusters_dut++;
      dutcolHisto->fill( c->col );

      // Fill statistics for first and last pixel in the cluster:
      // Used to monitor VIColOr effect, only look at long clusters:
      if(c->size > 2) {
	int colmin = 99;
	int rowmin = 99;
	double qmin = 0;
	int colmax = -1;
	int rowmax = -1;
	double qmax = 0;

	for( std::vector<CMSPixel::pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); px++ ){
	  if( px->col < colmin ) { colmin = px->col; qmin = px->vcal; }
	  if( px->col > colmax ) { colmax = px->col; qmax = px->vcal; }
	  if( px->row < rowmin ) { rowmin = px->row; }
	  if( px->row > rowmax ) { rowmax = px->row; }
	}//pix

	dutcolfirstHisto->fill( colmin );
	dutcollastHisto->fill( colmax );
	dutrowfirstHisto->fill( rowmin );
	dutrowlastHisto->fill( rowmax );
	// Statistics for first DC, charge of left and right columns:
	if(colmin%2 == 0) { dutcol1stevenqHisto->fill( qmin ); }
	else { dutcol1stoddqHisto->fill( qmin ); }
      }
      
      dutrowHisto->fill( c->row );
      dutnpxHisto->fill( c->size );
      dutadcHisto->fill( c->charge );
    }//DUT clusters
  }//pix


  // REF cluster statistics:
  refnclusHisto->fill( refclusters.size() );
  if( refpix > 0 ) {
    for( std::vector<cluster>::iterator c = refclusters.begin(); c != refclusters.end(); c++ ){
      n_clusters_ref++;
      refcolHisto->fill( c->col );
      refrowHisto->fill( c->row );
      refnpxHisto->fill( c->size );
      refadcHisto->fill( c->charge );
    }//REF clusters
  }//pixREF

}


std::vector<EUTelAnalysisCMSPixel::cluster> EUTelAnalysisCMSPixel::GetClusters(std::vector<CMSPixel::pixel> * pixels, int chip) {

  // FIXME: BEWARE! this currently handles only single ROCs!

  std::vector<cluster> clusters;
  if(pixels->empty()) return clusters;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Tsunami correction
  double eps = 0; // for analog ROCs

  if(chip >= 200) {
    //eps = 0.03; // tsunami correction for psi46digV1
    eps = 0.6; // [ke] additive
  }
  if(chip >= 400  ) {
    //eps = 0.06; // tsunami correction for psi46digV2.1
    eps = 1.1; // [ke] additive
  }

  // Apply tsunami correction for all pixels beside first:
  std::vector<CMSPixel::pixel>::iterator px = pixels->begin();
  std::advance(px, 1);
  for(; px != pixels->end(); px++ ) {
    //px->vcal -= eps; // additive
  }

  // Alternative by Daniel: correct by 3% of precursor:
  /*
  eps = 0.03;
  double previous_charge = 0;
  for(std::vector<CMSPixel::pixel>::iterator px = pixels->begin(); px != pixels->end(); px++) {
    double temp_charge = px->vcal;
    px->vcal -= eps*previous_charge;
    previous_charge = temp_charge;
  }
  */

  // clustering: 1 = no gap, 2 = gap of 1px
  int fCluCut = 1;
  //int fCluCut = 2;

  int* gone = new int[pixels->size()];
  for(size_t i = 0; i < pixels->size(); i++) gone[i] = 0;

  uint seed = 0;

  while( seed < pixels->size() ) {

    // Start a new cluster
    cluster c;
    c.vpix.push_back(pixels->at(seed));
    gone[seed]=1;
    c.sumA = 0;
    c.charge = 0;
    c.size = 0;
    c.col = 0;
    c.row = 0;
    //    c.xy[0] = 0;
    //c.xy[1] = 0;

    // Let it grow as much as possible:
    int growing;
    do{
      growing = 0;
      for(size_t i = 0; i < pixels->size(); i++ ){
        if( !gone[i] ){//unused pixel
	  for(std::vector<CMSPixel::pixel>::iterator p = c.vpix.begin(); p != c.vpix.end(); p++ ) {//vpix in cluster so far
	    int dr = p->row - pixels->at(i).row;
            int dc = p->col - pixels->at(i).col;
            if( dr >= -fCluCut && dr <= fCluCut &&
	        dc >= -fCluCut && dc <= fCluCut) {
              c.vpix.push_back(pixels->at(i));
	      gone[i] = 1;
              growing = 1;
              break;
            }
          }//loop over vpix
        }//not gone
      }//loop over all pix
    }while(growing);

    // sort pixels along col*80 + row:
    //c.vpix.sort();

    // Added all we could. determine position and append it to the list of clusters:
    for( std::vector<CMSPixel::pixel>::iterator p=c.vpix.begin();  p!=c.vpix.end();  p++){
      double Qpix = p->vcal; // calibrated [keV]
      if( Qpix < 0 ) Qpix = 1; // DP 1.7.2012
      c.charge += Qpix;
      c.sumA += p->raw; //Raw Chip pulse height value
      c.col += (*p).col*Qpix;
      c.row += (*p).row*Qpix;
    }

    c.size = c.vpix.size();

    streamlog_out(DEBUG4) << "(cluster with " << c.vpix.size() << " pixels)" << std::endl;

    if( !c.charge == 0 ) {
      c.col = c.col/c.charge;
      c.row = c.row/c.charge;
    }
    else {
      c.col = (*c.vpix.begin()).col;
      c.row = (*c.vpix.begin()).row;
      streamlog_out(DEBUG3) << "GetHits: cluster with zero charge" << std::endl;
    }

    // Mark the pixel containing the COG as seed pixel:
    int cog_col = static_cast<int>(round(c.col));
    int cog_row = static_cast<int>(round(c.row));
    for(size_t px = 0; px < c.vpix.size(); px++){
      if(c.vpix.at(px).col == cog_col && c.vpix.at(px).row == cog_row) {
	c.seed_id = px;
      }
    }

    clusters.push_back(c);//add cluster to vector

    //look for a new seed = unused pixel:
    while((++seed < pixels->size()) && (gone[seed]));
  }//while over seeds

  // Fill the pixel charge / readout order plots:
  if(chip == _DUT_chip) {
    // Select single 2px clusters:
    if(clusters.size() == 1 && clusters.front().size == 2) {
      if(pixels->size() == 2) {
	cms2pxq0Histo->fill(pixels->at(0).vcal);
	cms2pxq1Histo->fill(pixels->at(1).vcal);
      }
    }

    // Select single 4px clusters:
    if(clusters.size() == 1 && clusters.front().size == 4) {
      if(pixels->size() == 4) {
	cms4pxq0Histo->fill(pixels->at(0).vcal);
	cms4pxq1Histo->fill(pixels->at(1).vcal);
	cms4pxq2Histo->fill(pixels->at(2).vcal);
	cms4pxq3Histo->fill(pixels->at(3).vcal);
      }
    }
  }

  // nothing left to do, return clusters
  delete gone;
  return clusters;
}

//------------------------------------------------------------------------------
double EUTelAnalysisCMSPixel::GetConversionFactor(EUTelAnalysisCMSPixel::calibration cal, double conversion) {

  double keV = 0.45; // keV / large Vcal DAC
  // Check if we have an external VCal->Electrons conversion factor:
  if(conversion > 0.1) { keV = conversion; }
  // Old chip dependent calibration:
  else {
    if( cal.chip_id ==  10 ) keV = 0.36; // Landau peak at 24 ke, May 2012
    if( cal.chip_id ==  22 ) keV = 0.50; // Landau peak at 24 ke
    if( cal.chip_id >= 100 ) keV = 0.35; // xdb 15.7.2012 (7*0.05 = 0.35)
    if( cal.chip_id >= 110 && cal.chip_id <= 120) keV = 0.45; // AB -- sync to b2h, 5.10.2013
    //if( cal.chip_id == 202) keV = 0.307;
    if( cal.chip_id == 202) keV = 0.288; // VS -- data vs pixelav
    // comparison, Landau peak
    // Feb 2014
    if( cal.chip_id == 203) keV = 0.324;

    if( cal.chip_id == 400 ) keV = 0.288; // 12593 to get q0f peak at 22 ke
    if( cal.chip_id == 401 ) keV = 0.275; // 12603 to get q0f peak at 22 ke
    if( cal.chip_id == 402 ) keV = 0.295; // 12556 to get q0f peak at 22 ke
    if( cal.chip_id == 403 ) keV = 0.280; // copy  to get q0f peak at 22 ke
    if( cal.chip_id == 404 ) keV = 0.262; // 12525 to get q0f peak at 22 ke
    if( cal.chip_id == 405 ) keV = 0.301; // 12544 to get q0f peak at 22 ke

    //if( cal.chip_id == 500 ) keV = 0.312; // 14470 to get q0f peak at 22 ke
    if( cal.chip_id == 500 ) keV = 0.305; // 14393 to get q0f peak at 22 ke no eps in Q

    //if( cal.chip_id == 504 && _nRun < 20253) keV = 0.258; // 19045 to get q0f peak at 22 ke
    //else if( cal.chip_id == 504 ) keV = 0.2501; // 19045 to get q0f peak at 22 ke

    if( cal.chip_id == 502) keV = 0.225;
    

    if( cal.chip_id == 504) keV = 0.25721;
    if( cal.chip_id == 504 && _nRun >= 20785) keV = 0.246;
    if( cal.chip_id == 504 && _nRun >= 20811) keV = 0.24958;

    if( cal.chip_id == 504 && _nRun >= 20386 && _nRun <= 20398) keV = 0.245;

    if( cal.chip_id == 506 ) keV = 0.295; // 19582 chiller off, tilt 28

    if( cal.chip_id == 506 && _nRun >= 14697 && _nRun <= 14706) keV = 0.2897;

    if( cal.chip_id == 506 && _nRun >= 19714 ) keV = 0.282;
    if( cal.chip_id == 506 && _nRun >= 19714 && _nRun <= 19735) { keV = 0.275; }
    if( cal.chip_id == 506 && _nRun >= 19845 && _nRun <= 19857) { keV = 0.287; }
    if( cal.chip_id == 506 && _nRun >= 19866 && _nRun <= 19958) { keV = 0.2779; }
    if( cal.chip_id == 506 && _nRun >= 20165 && _nRun <= 20179) { keV = 0.2736; }
    if( cal.chip_id == 506 && _nRun >= 20749 && _nRun <= 20772) { keV = 0.2785; }
    
    if( cal.chip_id == 603 ) { keV = 0.289; }
  }
  return keV;
}

bool EUTelAnalysisCMSPixel::CalibratePixels(std::vector<CMSPixel::pixel> * pixels, EUTelAnalysisCMSPixel::calibration cal, double keV) {

  if(cal.chip_id == _DUT_chip) {
    // Fill the readout order plots:
    if(pixels->size() > 0) dutpx0adcHisto->fill(pixels->at(0).raw);
    if(pixels->size() > 1) dutpx1adcHisto->fill(pixels->at(1).raw);
    if(pixels->size() > 2) dutpx2adcHisto->fill(pixels->at(2).raw);
    if(pixels->size() > 3) dutpx3adcHisto->fill(pixels->at(3).raw);
    if(pixels->size() > 4) dutpx4adcHisto->fill(pixels->at(4).raw);
    if(pixels->size() > 5) dutpx5adcHisto->fill(pixels->at(5).raw);
    if(pixels->size() > 6) dutpx6adcHisto->fill(pixels->at(6).raw);
    if(pixels->size() > 7) dutpx7adcHisto->fill(pixels->at(7).raw);
    // 2px Clusters only
    if(pixels->size() == 2) dutpx8adcHisto->fill(pixels->at(0).raw);
    if(pixels->size() == 2) dutpx9adcHisto->fill(pixels->at(1).raw);
  }

  for( std::vector<CMSPixel::pixel>::iterator pix = pixels->begin(); pix != pixels->end(); pix++) {

    // try to average out the pixel 1 ADC response:
    if(pix == (pixels->begin()+1)) pix->raw -= (rand()%2);
    
    // Correct ADC shift:
    if(pix != pixels->begin()) pix->raw -= _adc_correction;

    size_t col = (*pix).col;
    size_t row = (*pix).row;

    int ff = 1;
    if(cal.chip_id >= 110 && cal.chip_id <= 113)
      ff = 4; // AB -- sync to b2h with tbm=true... 5.10.2013

    double Ared = (*pix).raw / ff - cal.fitParameter[3][col][row]; // sub vert off, see gaintanh2ps.C
    double ma9 = cal.fitParameter[0][col][row];
    if( Ared >  ma9-1 ) ma9 = 1.000001 * Ared;
    if( Ared < -ma9+1 ) ma9 = Ared - 1;
    
    // PSI Tanh Calibration (psi46expert vanilla):
    if(strcmp(cal.type.c_str(),"psi_tanh") == 0) {
      (*pix).vcal = (TMath::ATanH(((*pix).raw - cal.fitParameter[3][col][row])/cal.fitParameter[2][col][row])
		     + cal.fitParameter[1][col][row])/cal.fitParameter[0][col][row] *65E-3; // 65e/VcalDAC
    }
    // Weibull Calibration:
    else if(strcmp(cal.type.c_str(),"desy_weibull") == 0) {
      (*pix).vcal = (std::pow( - std::log( 1.0 - Ared / ma9 ), 1.0/cal.fitParameter[4][col][row]) 
		     * cal.fitParameter[1][col][row] + cal.fitParameter[2][col][row] ) * keV;
    }
    // Decorrelated Weibull Calibration:
    else if(strcmp(cal.type.c_str(),"desy_weibull_decorr") == 0) {

      // phcal2ps decorrelated: a = p4 + p3*exp(-t^p2), t = p0 + x/p1

      Ared = (*pix).raw - cal.fitParameter[4][col][row]; // p4 is asymptotic maximum
      if(Ared > 0) { Ared = -0.1; } // avoid overflow

      ma9 = cal.fitParameter[3][col][row]; // ma9 must be negative,
      if(ma9 > 0) ma9 = -ma9;             // if not, change sign

      (*pix).vcal = cal.fitParameter[1][col][row] * ( std::pow( -std::log( Ared / ma9 ), 1./cal.fitParameter[2][col][row] ) - cal.fitParameter[0][col][row] ) * keV;
      // q = ( (-ln((A-p4)/p3))^1/p2 - p0 )*p1
    }
    // Decorrelated Weibull with inverted sign:
    else if(strcmp(cal.type.c_str(),"desy_weibull_decorr_pos") == 0){

      // phroc2ps decorrelated: a = p4 - p3*exp(-t^p2), t = p0 + x/p1

      Ared = (*pix).raw - cal.fitParameter[4][col][row]; // p4 is asymptotic maximum
      if(Ared > 0) { Ared = -0.1; } // avoid overflow

      ma9 = cal.fitParameter[3][col][row]; // ma9 must be positive,
      if(ma9 < 0) ma9 = -ma9;             // if not, change sign

      (*pix).vcal = cal.fitParameter[1][col][row] * ( std::pow( -std::log( -Ared / ma9 ), 1./cal.fitParameter[2][col][row] ) - cal.fitParameter[0][col][row] ) * keV;
      // q = ( (-ln(-(A-p4)/p3))^1/p2 - p0 )*p1
    }
    // DESY TanH Calibration:
    else {
      (*pix).vcal = (TMath::ATanH( Ared / ma9 ) 
		      * cal.fitParameter[1][col][row] + cal.fitParameter[2][col][row]) * keV; // [ke]
    }
  }

  if(cal.chip_id == _DUT_chip) {
    // Fill the readout order plots:
    if(pixels->size() > 0) dutpx0cadcHisto->fill(pixels->at(0).raw);
    if(pixels->size() > 1) dutpx1cadcHisto->fill(pixels->at(1).raw);
    if(pixels->size() > 2) dutpx2cadcHisto->fill(pixels->at(2).raw);
    if(pixels->size() > 3) dutpx3cadcHisto->fill(pixels->at(3).raw);
    if(pixels->size() > 4) dutpx4cadcHisto->fill(pixels->at(4).raw);
    if(pixels->size() > 5) dutpx5cadcHisto->fill(pixels->at(5).raw);
    if(pixels->size() > 6) dutpx6cadcHisto->fill(pixels->at(6).raw);
    if(pixels->size() > 7) dutpx7cadcHisto->fill(pixels->at(7).raw);
    // 2px Clusters only
    if(pixels->size() == 2) dutpx8cadcHisto->fill(pixels->at(0).raw);
    if(pixels->size() == 2) dutpx9cadcHisto->fill(pixels->at(1).raw);
  }

  return true;
}

bool EUTelAnalysisCMSPixel::InitializeCalibration(std::string gainfilename, int chip_id, std::string CalibrationType, EUTelAnalysisCMSPixel::calibration & cal) {

  //FIXME initialize cal with 0?

  streamlog_out(MESSAGE4) << "CMSPixel gainfile: " << gainfilename 
			  << " of type " << CalibrationType << std::endl;
  std::ifstream gainFile(gainfilename.c_str());

  if(!gainFile) {
    streamlog_out(ERROR0) << "CMSPixel Pulse height calibration file not found!" << std::endl;
    return false;
  }

  char string[500];
  int icol;
  int irow;
  int a,b;
  double am;
  double ho;
  double ga;
  double vo;
  double ex = -1000; // non-Weibull flag

  int i = 0;

  // PSI Tanh Calibration (psi46expert vanilla):
  if(strcmp(CalibrationType.c_str(),"psi_tanh") == 0){
    gainFile.getline(string,500);
    gainFile.getline(string,500);
    gainFile.getline(string,500);
      
    for (int icol = 0; icol < 52; icol++) {
	for (int irow = 0; irow < 80; irow++) {
	    gainFile >> cal.fitParameter[0][icol][irow] >> cal.fitParameter[1][icol][irow]
		     >> cal.fitParameter[2][icol][irow] >> cal.fitParameter[3][icol][irow]
		     >> string >> a >> b;
	    i++;
	  }	  
      }
  }
  else if(strcmp(CalibrationType.c_str(),"desy_weibull") == 0) {
    while( gainFile >> string ) {
      gainFile >> icol;
      gainFile >> irow;

      gainFile >> cal.fitParameter[2][icol][irow];//horz offset
      gainFile >> cal.fitParameter[1][icol][irow];//width
      gainFile >> cal.fitParameter[4][icol][irow];//exponent
      gainFile >> cal.fitParameter[0][icol][irow];//gain
      gainFile >> cal.fitParameter[3][icol][irow];//vert offset

      i++;
    }
  }
  else if(strcmp(CalibrationType.c_str(),"desy_weibull_decorr") == 0 || strcmp(CalibrationType.c_str(),"desy_weibull_decorr_pos") == 0) {
    while( gainFile >> string ) {
      gainFile >> icol;
      gainFile >> irow;

      gainFile >> cal.fitParameter[0][icol][irow];
      gainFile >> cal.fitParameter[1][icol][irow];
      gainFile >> cal.fitParameter[2][icol][irow];
      gainFile >> cal.fitParameter[3][icol][irow];
      gainFile >> cal.fitParameter[4][icol][irow];
      gainFile >> cal.fitParameter[5][icol][irow];

      i++;
    }
  }
  else if(strcmp(CalibrationType.c_str(),"desy_tanh") == 0) {
    while( gainFile >> string ) {
      gainFile >> icol;
      gainFile >> irow;
      
      gainFile >> am;//Amax
      gainFile >> ho;//horz offset
      gainFile >> ga;//gain [ADC/large Vcal]
      gainFile >> vo;//vert offset
      
      cal.fitParameter[0][icol][irow] = am; // amax[icol][irow] = am*aa; // gain for Weib
      cal.fitParameter[1][icol][irow] = ga; // Gain[icol][irow] = ga; // width for Weib
      cal.fitParameter[2][icol][irow] = ho; // horz[icol][irow] = ho;
      cal.fitParameter[3][icol][irow] = vo; // vert[icol][irow] = vo;
      cal.fitParameter[4][icol][irow] = ex; // expo[icol][irow] = ex;

      i++;
    }
  }
  else { 
    streamlog_out(ERROR) << "Selected calibration method " << CalibrationType << " not supported!" << std::endl;
    return false;
  }

  streamlog_out(MESSAGE2) << "Read calibration parameters for " << i-1 << " pixels, type " 
			  << CalibrationType << std::endl;

  // Store calibration type and chip id:
  cal.type = CalibrationType;
  cal.chip_id = chip_id;

  return true;
}

Double_t fitLandauGauss( Double_t *x, Double_t *par ) {

  static int nn=0;
  nn++;
  static double xbin = 1;
  static double b1 = 0;
  if( nn == 1 ) { b1 = x[0]; }
  if( nn == 2 ) { xbin = x[0] - b1; } // bin width needed for normalization

  // Landau:
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location

  // MP shift correction:
  double mpc = par[0] - mpshift * par[1]; //most probable value (peak pos)

  //Fit parameters:
  //par[0] = Most Probable (MP, location) parameter of Landau density
  //par[1] = Width (scale) parameter of Landau density
  //par[2] = Total area (integral -inf to inf, normalization constant)
  //par[3] = Gaussian smearing

  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Range of convolution integral
  double xlow = x[0] - sc * par[3];
  double xupp = x[0] + sc * par[3];

  double step = (xupp-xlow) / np;

  // Convolution integral of Landau and Gaussian by sum

  double sum = 0;
  double xx;
  double fland;

  for( int i = 1; i <= np/2; i++ ) {

    xx = xlow + ( i - 0.5 ) * step;
    fland = TMath::Landau( xx, mpc, par[1] ) / par[1];
    sum += fland * TMath::Gaus( x[0], xx, par[3] );

    xx = xupp - ( i - 0.5 ) * step;
    fland = TMath::Landau( xx, mpc, par[1] ) / par[1];
    sum += fland * TMath::Gaus( x[0], xx, par[3] );
  }

  return( par[2] * invsq2pi * xbin * step * sum / par[3] );
}

Double_t EUTelAnalysisCMSPixel::landau_gauss_peak(TH1* h) {

  double aa = h->GetEntries();//normalization

  // find peak:
  int ipk = h->GetMaximumBin();
  double xpk = h->GetBinCenter(ipk);
  double sm = xpk / 9; // sigma
  double ns = sm; // noise

  // fit range:
  int ib0 = h->FindBin(18);
  int ib9 = h->FindBin(40);
  double x0 = h->GetBinLowEdge(ib0);
  double x9 = h->GetBinLowEdge(ib9) + h->GetBinWidth(ib9);

  // create a TF1 with the range from x0 to x9 and 4 parameters
  TF1 *fitFcn = new TF1( "fitFcn", fitLandauGauss, x0, x9, 4 );

  // set start values:
  fitFcn->SetParameter( 0, xpk ); // peak position, defined above
  fitFcn->SetParameter( 1, sm ); // width
  fitFcn->SetParameter( 2, aa ); // area
  fitFcn->SetParameter( 3, ns ); // noise

  h->Fit("fitFcn", "R Q", "ep" );// R = range from fitFcn
  TF1 *fit = h->GetFunction("fitFcn");
  return fit->GetParameter(0);
}

#endif
