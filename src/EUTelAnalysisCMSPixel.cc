
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
#include "EUTELESCOPE.h"
#include "EUTelSparseDataImpl.h"
#include "EUTelBaseSparsePixel.h"
#include "EUTelSimpleSparsePixel.h"
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

// ROOT includes ".h"
#include <TMath.h>
#include <TVectorD.h>
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


int dutnsync = 0;
int refnsync = 0;
int telnsync = 0;
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


EUTelAnalysisCMSPixel::EUTelAnalysisCMSPixel() : Processor("EUTelAnalysisCMSPixel"), _siPlanesParameters(), _siPlanesLayerLayout(), _inputCollectionTelescope(""), _inputCollectionDUT(""), _inputTrackCollection(""), _isFirstEvent(0), _nSkipTelescope(0), _gTLU(0), _nEvtBeg(0), _nRefBeg(0), _eBeam(0), _nRun(0), _nEvt(0), _leff_val(0), _nTelPlanes(0), time_event0(0), time_event1(0), time_reference(0), fTLU(0), gTLU(0), _TEL_run(0), _DUT_run(0), _DUT_data(""), _nSkipDUT(0), _DUT_chip(0), _DUT_gain(""), _DUT_address(""), _DUT_calibration_type(""), _DUT_decoding_flags(0), dut_calibration(), _DUTalignx(0), _DUTaligny(0), _DUTz(0), _DUTrot(0), _DUTtilt(0), _DUTturn(0), _REF_run(0), _REF_data(""), _nSkipRef(0), _REF_chip(0), _REF_gain(""), _REF_address(""), _REF_calibration_type(""), _REF_decoding_flags(0), ref_calibration(), _REFalignx(0), _REFaligny(0), _REFz(0), _REFrot(0), _CMS_gain_path(""), _CMS_data_path(""), _DATE_run(""), _gearfile(""), _DUT_board(""), _REF_board(""), _TimingRun(false), _planeSort(), _planeID(), _planePosition(), _planeThickness(), _planeX0(), _planeResolution(), _dutDecoder(), _refDecoder(), ClustDUT(), ClustREF(), time_now_dut(0), time_now_ref(0), n_uncorrelated_triggers_dut(0), n_uncorrelated_triggers_ref(0), n_uncorrelated_triggers_tlu(0) {

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
  registerOutputCollection(LCIO::TRACK,"InputTrackCollectionTelescope",
                           "Name of the input Track collection of the telescope",
                           _inputTrackCollection, std::string( "fittracks"));

  registerProcessorParameter( "Ebeam",
                              "Beam energy [GeV]",
			      _eBeam, static_cast < double >( 4.0));

  registerOptionalParameter( "nSkipEventDUT",
			     "Skip n CMS DUT events at begin of file",
			     _nSkipDUT, static_cast < int >(0) );
  registerOptionalParameter( "nSkipEventREF",
			     "Skip n CMS REF events at begin of file",
			     _nSkipRef, static_cast < int >(0) );
  registerOptionalParameter( "nSkipEventTel",
			     "Skip n Telescope events at begin of file",
			     _nSkipTelescope, static_cast < int >(0) );
  registerOptionalParameter( "SkipScanRangeDUT",
                             "Range for DUT skip scanner",
			     _nSkipScanRangeDUT, static_cast < int >(0) );

  registerProcessorParameter( "gTLU",
                              "TLU clock in GHz",
			      _gTLU, static_cast < double >(0.3846776));

  registerProcessorParameter("CMS_data_path",
                             "Path to the native CMS data files",
                             _CMS_data_path, std::string("/dev/null"));
  registerProcessorParameter("CMS_gain_path",
                             "Path to the CMS gain calibration files",
                             _CMS_gain_path, std::string("/dev/null"));


  registerProcessorParameter("DUT_run",
                             "Runnumber of the DUT associated with this telescope run",
                             _DUT_run, static_cast< int > (0) );
  registerProcessorParameter("REF_run",
                             "Runnumber of the REF associated with this telescope run",
                             _REF_run, static_cast< int > (0) );

  registerProcessorParameter("DUT_chip",
                             "Serial number of the DUT used in this run",
                             _DUT_chip, static_cast< int > (0) );
  registerProcessorParameter("DUT_gain",
                             "CMS DUT gain file to be used",
                             _DUT_gain, std::string("/dev/null"));
  registerProcessorParameter("DUT_address",
                             "CMS DUT address levels file to be used",
                             _DUT_address, std::string("/dev/null"));
  registerProcessorParameter("DUT_calibration",
                             "Choose DUT calibration type (psi_tanh, desy_tanh, desy_weibull)",
                             _DUT_calibration_type, std::string("desy_weibull"));
  registerProcessorParameter("DUT_decoding_flags",
                             "Specify flags DUT for data decoding (give summarized int value)",
                             _DUT_decoding_flags, static_cast<int>(0));

  registerProcessorParameter("REF_chip",
                             "Serial number of the REF used in this run",
                             _REF_chip, static_cast< int > (0) );
  registerProcessorParameter("REF_gain",
                             "CMS REF gain file to be used",
                             _REF_gain, std::string("/dev/null"));
  registerProcessorParameter("REF_address",
                             "CMS REF address levels file to be used",
                             _REF_address, std::string("/dev/null"));
  registerProcessorParameter("REF_calibration",
                             "Choose REF calibration type (psi_tanh, desy_tanh, desy_weibull)",
                             _REF_calibration_type, std::string("desy_weibull"));
  registerProcessorParameter("REF_decoding_flags",
                             "Specify flags REF for data decoding (give summarized int value)",
                             _REF_decoding_flags, static_cast<int>(0));

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

  // Stuff only needed for the printout of the updated runlist line:
  registerOptionalParameter( "DUT_board",
			     "PSI testboard ID of the DUT chip",
			     _DUT_board, std::string("0") );
  registerOptionalParameter( "REF_board",
			     "PSI testboard ID of the REF chip",
			     _REF_board, std::string("0") );
  registerOptionalParameter( "date",
			     "Date when the test beam was taken",
			     _DATE_run, std::string("1970-01-01") );
  registerOptionalParameter( "gearfile",
			     "Again, the gear fiel of the used telescope configuration",
			     _gearfile, std::string("none") );
  
  registerOptionalParameter( "TimingRun",
			     "Switch on and off timing run mode. This will enable or disable trigger syncronization and can be used to determine TLU clock frequencies",
			     _TimingRun, static_cast<bool>(false) );
  
}


void EUTelAnalysisCMSPixel::init() {

  // usually a good idea to
  printParameters();

  _nRun = 0;
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

  if(_nSkipScanRangeDUT != 0)
    _dutSkipScanner.reset(new CMSSkipScanner);

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  bookHistos();
#endif
  mille = new gbl::MilleBinary( "mille.bin" );

}//init

//------------------------------------------------------------------------------
void EUTelAnalysisCMSPixel::processRunHeader( LCRunHeader* runHeader) {

  std::auto_ptr<EUTelRunHeaderImpl> eutelHeader( new EUTelRunHeaderImpl( runHeader ) );
  eutelHeader->addProcessor( type() );

  _nRun++;

  // Decode and print out Run Header information - just a check

  streamlog_out( MESSAGE2 )  << "Processing run header " << _nRun
                             << ", run nr " << runHeader->getRunNumber() << std::endl;
  _TEL_run = runHeader->getRunNumber();

  const std::string detectorName = runHeader->getDetectorName();
  const std::string detectorDescription = runHeader->getDescription();

  streamlog_out( MESSAGE0 ) << detectorName << " : " << detectorDescription << std::endl;


  // ################### INIT CMS PIXEL DATA STREAMS #######################

  // Setting the CMSPixelDecoder Logging mechanism to report only errors:
  CMSPixel::Log::ReportingLevel() = CMSPixel::Log::FromString("ERROR");

  std::string dut_file = _CMS_data_path + "/bt05r" + ZeroPadNumber(_DUT_run,6) + "/mtb.bin";
  std::string ref_file = _CMS_data_path + "/bt05r" + ZeroPadNumber(_REF_run,6) + "/mtb.bin";

  std::string dut_addresslevels = _CMS_gain_path + "/chip" + ZeroPadNumber(_DUT_chip,3) + "/" + _DUT_address;
  std::string ref_addresslevels = _CMS_gain_path + "/chip" + ZeroPadNumber(_REF_chip,3) + "/" + _REF_address;

  streamlog_out(MESSAGE3) << "DUT chiptype detected: " << static_cast<int>(GetChipTypeFromID(_DUT_chip)) 
			  << ", loading address levels from " << dut_addresslevels
			  << std::endl;
  streamlog_out(MESSAGE3) << "REF chiptype detected: " << static_cast<int>(GetChipTypeFromID(_REF_chip))
			  << ", loading address levels from " << ref_addresslevels
			  << std::endl;

  _dutDecoder = new CMSPixel::CMSPixelFileDecoderPSI_ATB(dut_file.c_str(),1,_DUT_decoding_flags,
							 GetChipTypeFromID(_DUT_chip),dut_addresslevels.c_str());
  _refDecoder = new CMSPixel::CMSPixelFileDecoderPSI_ATB(ref_file.c_str(),1,_REF_decoding_flags,
							 GetChipTypeFromID(_REF_chip),ref_addresslevels.c_str());

  // DUT calibration:
  if(!InitializeCalibration(_CMS_gain_path + "/chip" + ZeroPadNumber(_DUT_chip,3) + "/" + _DUT_gain,
			    _DUT_chip, _DUT_calibration_type, dut_calibration))
    throw marlin::ParseException("Error in processing DUT calibration file.");

  // REF calibration:
  if(!InitializeCalibration(_CMS_gain_path + "/chip" + ZeroPadNumber(_REF_chip,3) + "/" + _REF_gain,
			    _REF_chip, _REF_calibration_type, ref_calibration))
    throw marlin::ParseException("Error in processing REF calibration file.");

  n_uncorrelated_triggers_tlu = 0;
} // processRunHeader

//----------------------------------------------------------------------------
void EUTelAnalysisCMSPixel::processEvent( LCEvent * event ) {

  EUTelEventImpl * euEvent = static_cast<EUTelEventImpl*>( event );
  std::vector<CMSPixel::pixel> * dutPixels = new std::vector<CMSPixel::pixel>;
  std::vector<CMSPixel::pixel> * refPixels = new std::vector<CMSPixel::pixel>;

  if( euEvent->getEventType() == kEORE ) {
    streamlog_out( DEBUG5 ) <<  "EORE found: nothing else to do." << std::endl;
    return;
  }

  int runNumber = event->getRunNumber();

  //----------------------------------------------------------------------------
  // skip first Telescope events
  if(_nEvt < _nSkipTelescope) {
    if(_nEvt == 0) streamlog_out(MESSAGE5) << "Skip first " << _nSkipTelescope << " Telescope events" << std::endl;
    _nEvt++;
    throw SkipEventException(this);
  }
  //----------------------------------------------------------------------------

  if( _isFirstEvent ) {

    // Initialize the timing calculations:
    time_event0 = event->getTimeStamp();
    time_event1 = time_event0;
    time_reference = time_event0;

    // Temporarily mute CMS Pixel decoder, we don't want to hear anything for skipped events:
    CMSPixel::Log::ReportingLevel() = CMSPixel::Log::FromString("QUIET");
    
    // skip first DUT events
    streamlog_out(MESSAGE5) << "Skip first " << _nSkipDUT << " DUT events" << std::endl;
    for( int iskip = 0; iskip < _nSkipDUT; ++iskip ){
      if(_dutDecoder->get_event(dutPixels, evt_time_dut) <= DEC_ERROR_NO_MORE_DATA) break;
    }

    // skip first REF events
    streamlog_out(MESSAGE5) << "Skip first " << _nSkipRef << " REF events" << std::endl;
    for( int iskip = 0; iskip < _nSkipRef; ++iskip ) {
      if(_refDecoder->get_event(refPixels, evt_time_ref) <= DEC_ERROR_NO_MORE_DATA) break;
    }

    // Restore error reporting of the decoders:
    CMSPixel::Log::ReportingLevel() = CMSPixel::Log::FromString("ERROR");

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
  // dutddt: TLU/CMS = 1.00007
  // => TLU/1.00007 = CMS
  // dutddtns: <TLU/1.00007344 - CMS> = 0

  gTLU = _gTLU; // [GHz]

  if( runNumber >  3000 ) gTLU = 0.384666; // [GHz] 2012 TB 21 Apr/May 2012
  if( runNumber >  3672 ) gTLU = 0.384677; // [GHz] 2012 TB 21 xdb and dig, new TLU box
  if( runNumber >  5000 ) gTLU = 0.3846775; // [GHz] Sep 2012 TB 21
  if( runNumber >  5100 ) gTLU = 0.384678; // [GHz] Dec 2012 TB 21
  if( runNumber >  6120 ) gTLU = 0.3846776; // [GHz] Feb 2013 TB 21
  if( runNumber >= 7709 ) gTLU = 0.3859129; // [GHz] 22.3.2013 ETH 25ns
  if( runNumber >= 7965 ) gTLU = 0.3846776; // [GHz] 26.3.2013 KIT 24.9 ns
  if( runNumber >= 9600 ) gTLU = 0.3846745; // [GHz] 10.6.2013 new TLU  AA002 TLUv0.2c
  if( runNumber >= 9649 ) gTLU = 0.3846776; // [GHz] 11.6.2013 old TLU

  if( runNumber >= 11111 && runNumber <= 11132 )
    gTLU = 0.3840282; // [GHz] internal clock

  if( runNumber >= 11251 && runNumber <= 11264 )
    gTLU = 0.3840282; // [GHz] internal clock

  if( runNumber == 11347 ) // 80 MHz, 78*12.4 ns
    gTLU = 0.3840282/0.98549; // [GHz]

  if( runNumber >= 11351 && runNumber <= 11381 ) // 4.10.2013, 80 MHz, 78*12.4 ns
    gTLU = 0.3840282*1.0017*(1-70/10e6); // [GHz]

  // 7./8.10/2013: chip 202i, 203i

  if( runNumber >= 11415 && runNumber <= 11476 )  // 8.10.2013, 80 MHz, 78*12.6 ns
    gTLU = 0.384678355; // [GHz]

  if( runNumber >= 11495 && runNumber <= 11600 )  // 10.10.2013, 78 MHz, 76*12.9 ns
    gTLU = 0.394800638; // [GHz]

  if( runNumber >= 12016 && runNumber <= 12090 )  // 28.10.2013, 78 MHz, 76*12.9 ns
    gTLU = 0.394800638; // [GHz] FPIX

  if( runNumber >= 12096 && runNumber <= 12211 )  // 31.10.2013, 80 MHz, 78*12.6 ns
    gTLU = 0.384678355; // [GHz]

  // Timing run: disable skipping!
  // from dutddtvsdt: gtlu_new = gtlu_old * ( 1 + slope )

  double refTLU = gTLU;
  if( runNumber >= 8515 ) {
    refTLU = 0.3840282; // ref uses internal clock!
  }

  //DP double fTLU = gTLU * 1E9; // [Hz] w.r.t. 40 MHz
  // fTLU = Frequency Trigger Logic Unit; Clock freq of the TLU box in Hz
  // gTLU = gigahertz Frequency Trigger Logic Unit; Freq in GHz
  //fTLU = 375.056 * 1.024 * 1E6; // [Hz] = 384.05734 MHz, from Umlauftakt, see below
  fTLU = gTLU * 1E9;

  double eventTime =  (time_now_tlu-time_event0)/fTLU;

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
  
  //----------------------------------------------------------------------------
  // CMS pixel:

  // Horizontal DUT alignment for FPIX test beam:
  bool FPIX = 0;
  if( runNumber >= 5657 && runNumber < 6140 ) FPIX = 1;

  // ETH and KIT test beam Mar 2013 analog ROCs
  // Horizontal DUT alignment for ETH test beam:
  bool ETHh = 0;
  if( runNumber >= 7405 && runNumber <= 8500 ) ETHh = 1; // horizontal
  if( runNumber >= 9695 && runNumber <= 9933 ) ETHh = 1;


  //-----------------------------------------------------------------------
  // SKIPPING EVENTS
  //-----------------------------------------------------------------------
  /*
  if( n_uncorrelated_triggers_dut > 2 ) { 
    streamlog_out(MESSAGE5) << "Event " << event->getEventNumber() 
			    << ": Resync DUT" << std::endl;
    _dutDecoder->get_event( dutPixels, evt_time_dut);
    dutnsync++;
    n_uncorrelated_triggers_dut = 0; // reset
  }

  if( n_uncorrelated_triggers_ref > 2 ) {
    streamlog_out(MESSAGE5) << "Event " << event->getEventNumber() 
			    << ": Resync REF" << std::endl;
    _refDecoder->get_event(refPixels, evt_time_ref);
    refnsync++;
    n_uncorrelated_triggers_ref = 0; // reset
  }
  // Daniel : add else ?
  if( n_uncorrelated_triggers_tlu > 2 ) {
    streamlog_out(MESSAGE5) << "Event " << event->getEventNumber() 
			    << ": Resync Telescope" << std::endl;
    telnsync++;
    n_uncorrelated_triggers_tlu = 0;
    throw SkipEventException(this);
  }
  */
  //---------------------------------------------------------------------
  // Read useful events
  //---------------------------------------------------------------------

  // Read the DUT event:

  if( n_uncorrelated_triggers_dut < 3 ) { // DP 20.1.2014
    if(_dutDecoder->get_event(dutPixels, evt_time_dut) <= DEC_ERROR_NO_MORE_DATA)
      throw StopProcessingException(this);
    time_now_dut = evt_time_dut.timestamp;
  }
  else {
    streamlog_out(MESSAGE5) << "Event " << event->getEventNumber() 
			    << ": Resync DUT" << std::endl;
    dutnsync++;
    n_uncorrelated_triggers_dut = 0; // reset
  }
  // Calibrate the pixel hits with the initialized calibration data:
  if(!CalibratePixels(dutPixels,dut_calibration))
    throw StopProcessingException(this);

  ClustDUT = GetClusters(dutPixels);

  // Read the REF event:
  if( n_uncorrelated_triggers_ref < 3 ) { // DP 20.1.2014
    if(_refDecoder->get_event(refPixels, evt_time_ref) <= DEC_ERROR_NO_MORE_DATA)
      throw StopProcessingException(this);
    time_now_ref = evt_time_ref.timestamp;
  }
  else {
    streamlog_out(MESSAGE5) << "Event " << event->getEventNumber() 
			    << ": Resync REF" << std::endl;
    refnsync++;
    n_uncorrelated_triggers_ref = 0; // reset
  }
  // Calibrate the pixel hits with the initialized calibration data:
  if(!CalibratePixels(refPixels,ref_calibration))
    throw StopProcessingException(this);
  ClustREF = GetClusters(refPixels);

  // Stop processing if we can't find any correlation:
  if(!_TimingRun &&  (refnsync > 1000 || dutnsync > 1000)) {
    streamlog_out(ERROR) << "Too many re-syncs, probably either starting event does not match or relative clock is incorrect (gTLU)." << std::endl;
    throw StopProcessingException( this ) ;
  }

  static long dutTimeStamp0 = time_now_dut; // init
  static long refTimeStamp0 = time_now_ref; // init
  static long time_prev_dut = time_now_dut; // init
  static long time_old_dut = time_now_dut; // init
  static long time_prev_ref = time_now_ref; // init
  static long time_old_ref = time_now_ref; // init

  static long dutTimeReference = time_now_dut;
  static long refTimeReference = time_now_ref;


  static long prevgapTimeStamp = time_now_tlu; // init

  streamlog_out(DEBUG4) << std::setw(6) << std::setiosflags(std::ios::right) << event->getEventNumber() << ". DUT pixels " << dutPixels->size();

  if( dutPixels->size() > 0 ) {
    streamlog_out(DEBUG3) << ", clusters at";
    for( std::vector<cluster>::iterator c = ClustDUT.begin(); c != ClustDUT.end(); c++ ){
      streamlog_out(DEBUG3) << "  (" << c->col << ", " << c->row << ")";
    }
  }
  streamlog_out(DEBUG4) << std::endl;

  // FIXME Not understood: for digital chips seems to be necessary to check the
  // time gap before the actual event...?!

  long dutdt = time_now_dut - time_prev_dut; // [clocks]

  if( GetChipTypeFromID(_DUT_chip) >= CMSPixel::ROC_PSI46DIG )
    dutdt = time_prev_dut - time_old_dut;

  long refdt = time_now_ref - time_prev_ref; // [clocks]

  if( GetChipTypeFromID(_REF_chip) >= CMSPixel::ROC_PSI46DIG )
    refdt = time_prev_ref - time_old_ref;

  long dutrefddt = (dutdt/refTLU - refdt/gTLU) / 0.040; // [clocks]
  static long prevdutrefddt = dutrefddt;

  // ! Important:
  // Calculating the timing difference between TLU and DUT/REF for the
  // given event:

  double reference_tlu = ( time_now_tlu - time_reference ) / gTLU/1E3;

  double dutddtns = ( time_now_tlu - time_prev_tlu ) / gTLU - dutdt / 0.040;
  /*double dutddtnsX =
    ( time_now_tlu - time_prev_tlu ) / gTLU -
    ( time_now_dut - time_prev_dut ) / 0.040;
  */

  double reference_dut = (time_now_dut - dutTimeReference) / 40;

  // digital chips were different until 16.9.2013

  if( _DUT_chip > 100 && runNumber < 11098 )
    reference_dut = (time_prev_dut - dutTimeReference) / 40;

  // ref chip:

  double refddtns = ( time_now_tlu - time_prev_tlu )/refTLU - refdt / 0.040; // REF

  double reference_ref = (time_now_ref - refTimeReference) / 40;

  if( _REF_chip > 100 && runNumber < 11098 )
    reference_ref = (time_prev_ref - refTimeReference) / 40;
  
  if( time_now_dut > dutTimeStamp0 && time_now_tlu > time_event0 ) {

    cmsdtHisto->fill( dutdt / 40E0 ); // [us]
    dutrefddtHisto->fill( dutrefddt ); // [clocks]

    // Time of arrival TLU triggers:
    tlutrigvstusHisto->fill(int((time_now_tlu)/gTLU/1E6 )); // [ms]

    // time between events, TLU - DUT:
    dutddtnsHisto->fill( dutddtns ); // [ns]
    refddtnsHisto->fill( refddtns ); // [ns]
    dutddtusHisto->fill( dutddtns/1E3 ); // [us]. peak at 25 us in run 5348
    dutddtmsHisto->fill( dutddtns/1E6 ); // [ms]
    dutddtvst->fill( ( time_now_tlu - time_event0 ) / fTLU, dutddtns );
    dutddtvsdt->fill( ( time_now_tlu - time_prev_tlu ) / gTLU/1E3, dutddtns );

    if( abs( dutddtns ) < 99 ) {
      ddtvst->fill( ( time_now_tlu - time_event0 ) / fTLU, dutddtns );
      ddtvstms->fill( ( time_now_tlu - time_event0 ) / fTLU, dutddtns ); // no trend
      ddtvsdt->fill( ( time_now_tlu - time_prev_tlu ) / gTLU/1E3, dutddtns );
    }
    else {
      long dtgap = time_now_tlu - prevgapTimeStamp;
      gapdtHisto->fill( dtgap / fTLU ); // [s]
      prevgapTimeStamp = time_now_tlu;
    }

    if( _nSkipDUT == 0 ) {
      sysrtHisto->fill( ( time_now_tlu - time_event0 ) / fTLU /
			( time_now_dut - dutTimeStamp0 ) * 40E6 ); // mean = 1.00007
      sysrdtHisto->fill( ( time_now_tlu - time_prev_tlu) / fTLU /
			 ( time_now_dut - time_prev_dut ) * 40E6 );
    }
    else {
      sysrtHisto->fill( ( time_now_tlu - time_event0 ) / fTLU /
			( time_prev_dut - dutTimeStamp0 ) * 40E6 );
      sysrdtHisto->fill( ( time_now_tlu - time_prev_tlu) / fTLU /
			 ( time_prev_dut - time_old_dut ) * 40E6 );
    }
  } // > time_event0
   

  // Renew the timing reference point when having good correlation:
  if( abs(dutddtns) < 10 && abs(refddtns) < 10 && abs(dutrefddt) < 10 ) {
    time_reference = time_now_tlu;

    if( _DUT_chip >= 100 ) dutTimeReference = time_prev_dut;
    else  dutTimeReference = time_now_dut;

    if( _REF_chip >= 100 ) refTimeReference = time_prev_ref;
    else refTimeReference = time_now_ref;

    streamlog_out(DEBUG5) << "Refreshed reference timing." << std::endl;
  }

  // Only check for uncorrelated triggers when not in timing mode
  if( !_TimingRun && !_dutSkipScanner.get() ) {

    // Important step for correlation:
    // Check the timing differences for DUT, REF and TEL
    /*DP
    // Case: DUT timing is off:  abs added DP 15.11.2013
    if( ( dutddtns  > 99 && abs( reference_tlu - reference_dut ) > 10 ) ||
	( dutrefddt > 99 && abs( reference_ref - reference_dut ) > 10 ) ) n_uncorrelated_triggers_dut++;

    // Case: REF timing is off:  abs added DP 15.11.2013
    if( ( refddtns  > 99 && abs( reference_tlu - reference_ref ) > 10 ) ||
	( dutrefddt > 99 && abs( reference_dut - reference_ref ) > 10 ) ) n_uncorrelated_triggers_ref++;

    // Case: TLU timing is off: abs added DP 15.11.2013
    if( ( dutddtns > 99 && abs( reference_dut - reference_tlu ) > 10 ) ||
	( refddtns > 99 && abs( reference_ref - reference_tlu ) > 10 ) ) n_uncorrelated_triggers_tlu++;
    */
    //DP version:

    if( abs(dutddtns) > 99 ) n_uncorrelated_triggers_dut++;
    if( abs(refddtns) > 99 ) n_uncorrelated_triggers_ref++;
    if( abs(dutddtns) > 99 && abs(refddtns) > 99 ) n_uncorrelated_triggers_tlu++;

    // Reset counter if timing is okay with respect to the others:
    if( abs(dutddtns) < 99 && abs(dutrefddt) < 99 ) n_uncorrelated_triggers_dut = 0;
    if( abs(refddtns) < 99 && abs(dutrefddt) < 99 ) n_uncorrelated_triggers_ref = 0;
    if( abs(dutddtns) < 99 && abs(refddtns)  < 99 ) n_uncorrelated_triggers_tlu = 0;

    // Reset all counters for long gaps between events, clocks
    // deviate too much:
    if( ( time_now_tlu - time_prev_tlu ) / fTLU > 0.010 ) { // [s] long gap
      n_uncorrelated_triggers_dut = 0;
      n_uncorrelated_triggers_ref = 0;
      n_uncorrelated_triggers_tlu = 0;
      streamlog_out(DEBUG5) << "Long Gap detected, resetting trigger count in event " << _nEvt << "." << std::endl;
    }

    //FIXME this is a rare occasion: mixture of long gap and DUT missed trigger,
    // Need to think about a solution...
    if( runNumber == 11059 && event->getEventNumber() == 8152 ) {
      n_uncorrelated_triggers_tlu = 3;
      n_uncorrelated_triggers_ref = 3;
    }

    if(runNumber == 11198 && (event->getEventNumber() == 17882 || event->getEventNumber() == 17884)) {
      n_uncorrelated_triggers_tlu = 3;
      n_uncorrelated_triggers_ref = 3;
    }

    if(runNumber == 11196 && (event->getEventNumber() == 6084 || event->getEventNumber() == 6086 || event->getEventNumber() == 6088)) {
      n_uncorrelated_triggers_tlu = 3;
      n_uncorrelated_triggers_ref = 3;
    }

    if(runNumber == 11196 && (event->getEventNumber() == 154730 )) {
      n_uncorrelated_triggers_tlu = 3;
    }

  } // timing check wanted

  // Informational message about event number and timing:

  if( streamlog_level(DEBUG5) && 
      ( n_uncorrelated_triggers_tlu != 0 
	|| n_uncorrelated_triggers_dut != 0 
	|| n_uncorrelated_triggers_ref != 0 
	|| _nEvt < 100000
	|| _nEvt % 1000 == 0 ) ) {

    streamlog_out( DEBUG5 ) << "Processing event "
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

    streamlog_out( DEBUG5 ) << "       DUT event "
			    << std::setw( 7) << std::setiosflags(std::ios::right)
			    << event->getEventNumber()
			    << " in run " << std::setw( 4) << std::setiosflags(std::ios::right)
			    << runNumber
			    << ", time " << " = " << std::fixed << std::setprecision(1)
			    << ( time_now_dut - dutTimeStamp0 )/fTLU*10 << " s"
			    << ", dt = " << int(dutdt/40) << " us"
			    << ", ddt "
			    << int((time_now_tlu - time_prev_tlu)/gTLU/1E3 - dutdt/40)
			    << " or "  
			    << int( (time_now_tlu - time_prev_tlu)/gTLU/1E3 -
				    (time_now_dut - time_prev_dut)/40 )
      //<< " us (" << dutddtns << ")"
			    << std::endl;

    streamlog_out( DEBUG5 ) << "       REF event "
			    << std::setw( 7) << std::setiosflags(std::ios::right)
			    << event->getEventNumber()
			    << " in run " << std::setw( 4) << std::setiosflags(std::ios::right)
			    << runNumber
			    << ", time " << " = " << std::fixed << std::setprecision(1)
			    << ( time_now_ref - refTimeStamp0 )/fTLU*10 << " s"
			    << ", dt = " << int(refdt/40) << " us"
			    << ", ddt "
			    << int((time_now_tlu - time_prev_tlu)/gTLU/1E3 - refdt/40)
			    << " us"
			    << " (" << refddtns << ")"
			    << std::endl;

    streamlog_out( DEBUG5 ) << "   Reference Time tlu "
			    << std::setw( 7) << std::setiosflags(std::ios::right)
			    << reference_tlu
			    << " dut " << std::setw( 4) << std::setiosflags(std::ios::right)
			    << reference_dut
			    << " ref " << std::setw( 4) << std::setiosflags(std::ios::right)
			    << reference_ref  
			    << ", scalers: "
			    << " " << n_uncorrelated_triggers_tlu 
			    << "|" << n_uncorrelated_triggers_dut
			    << "|" << n_uncorrelated_triggers_ref
			    << " dutrefddt " << dutrefddt
			    << std::endl;
     
  }
  else if( _nEvt < 20 || _nEvt % 1000 == 0 ) {

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
  
  // Fill some informational plots about the REF and DUT cluster distributions:

  FillClusterStatisticsPlots( ClustDUT, dutPixels->size(), ClustREF, refPixels->size() );
  if( _dutSkipScanner.get() ) {
    _dutSkipScanner->newEvent(event->getEventNumber());
    for( std::vector<cluster>::iterator c = ClustDUT.begin(); c != ClustDUT.end(); c++ )
      _dutSkipScanner->cmsHit(c->col, c->row);
  }

  // Renew the timestamps for the next event:
  time_old_dut = time_prev_dut;
  time_old_ref = time_prev_ref;
  time_prev_dut = time_now_dut;
  time_prev_ref = time_now_ref;
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

  LCCollection* collection;
  try {
    collection = event->getCollection( _inputCollectionTelescope );
  }
  catch( lcio::DataNotAvailableException& e) {
    streamlog_out( DEBUG5 ) << "Not able to get collection "
			    << _inputCollectionTelescope
			    << "\nfrom event " << event->getEventNumber()
			    << " in run " << runNumber  << std::endl;

    return;
  }

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


  // Fill the telescope plane correlation plots:
  TelescopeCorrelationPlots(hits);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Tsunami correction

  double eps = 0; // for analog ROCs

  if( _DUT_chip >= 200  )
    //eps = 0.03; // tsunami correction for psi46digV1
    eps = 0.6; // [ke] additive

  if( _DUT_chip >= 400  )
    //eps = 0.06; // tsunami correction for psi46digV2.1
    eps = 1.2; // [ke] additive

  if( dutPixels->size() > 0 ) {

    for( std::vector<cluster>::iterator c = ClustDUT.begin(); c != ClustDUT.end(); c++ ){

      if( c->size == 1 ) continue;

      // look if some other pixel was readout earlier in the same DC:

      for( std::vector<CMSPixel::pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); px++ ) {

	int icol = px->col;
	int dcp = icol / 2; // 0..25 double column number

	bool hasPrevious = 0;
	double qprv = 0;

	for( std::vector<CMSPixel::pixel>::iterator qx = c->vpix.begin(); qx != c->vpix.end(); qx++ ) {

	  if( px == qx ) continue; // the same pixel

	  if( qx->col > px->col ) continue; // qx later than px

	  int dcq = qx->col / 2;

	  if( dcp != dcq ) continue; // want same double column

	  if( qx->col < px->col ) {
	    hasPrevious = 1;
	    qprv = qx->vcal;
	    continue;
	  }

	  // px and qx are in same column.

	  if( icol%2 == 0 ) { // ascending readout
	    if( qx->row < px->row ) {
	      hasPrevious = 1;
	      qprv = qx->vcal;
	    }
	  }
	  else
	    if( qx->row > px->row ) {
	      hasPrevious = 1;
	      qprv = qx->vcal;
	    }

	} // qx

	// p is 2nd px in DC readout: apply tsunami correction:

	if( hasPrevious )
	  // px->vcal -= eps * qprv; // proportional, overwrite!
	  px->vcal -= eps; // subtract constant, overwrite!

      } // pix

      // sum up corrected pixel charges:

      double q = 0;
      for( std::vector<CMSPixel::pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); px++ )
	q += px->vcal;

      c->charge = q; // overwritten !

    } // DUT clusters

  } // have DUT cust

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

  if( FPIX || ETHh ) {
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

  double norm = cos( turn*wt ) * cos( tilt*wt ); // length of Nz. Landau 21.6
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

  if(_dutSkipScanner.get())
    for( std::vector<triplet>::iterator trip = upstream_triplets->begin(); trip != upstream_triplets->end(); trip++ )
      _dutSkipScanner->telescopeTriplet(trip->getdx(1)*1E3, trip->getdy(1)*1E3);

  // Iterate over all found downstream triplets to fill histograms and match them to the REF and DUT:
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
    if( FPIX || ETHh ) { // x = col = yt, y = row = xt
      xmod = fmod( 9.075 + yAt, 0.3 ) * 1E3; // [0,300] um, 2 pixel wide
      ymod = fmod( 9.030 + xAt, 0.2 ) * 1E3; // [0,200] um
      ymd3 = fmod( 9.030 + xAt, 0.3 ) * 1E3; // [0,300] um
      ymd6 = fmod( 9.030 + xAt, 0.6 ) * 1E3; // [0,600] um
    }

    // 0 deg:
    bool ldot = 1; // bias dot, from cmsqvsxmym
    if( xmod < 105 ) ldot = 0; // dot at x = 125
    if( xmod > 195 ) ldot = 0; // and at x = 175
    if( tilt < 6 ) {
      if( ymod <  55 ) ldot = 0; // dot at y =  75
      if( ymod > 195 ) ldot = 0; // dot at y = 175
      if( ymod >  95 && ymod < 155 ) ldot = 0; // band between dots
    }

    bool lcore = 1; // pixel core, 2x2 pixel region
    if( xmod <  20 ) lcore = 0; // outer edge, see cmsncolvsxm
    if( xmod > 280 ) lcore = 0; // outer edge
    if( ymod <  20 ) lcore = 0; // outer edge, see cmsnrowvsym
    if( ymod > 180 ) lcore = 0; // outer edge
    if( xmod > 130 && xmod < 170 ) lcore = 0; // inner edge
    if( ymod >  80 && ymod < 120 ) lcore = 0; // inner edge


    // Extrapolate Upstream triplet to Downstream planes 3,4,5: Resolution studies
    for( std::vector<hit>::iterator lhit = hits->begin(); lhit != hits->end(); lhit++ ){

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

      bool leff = 1;

      if( runNumber == 11269 && eventTime >   45 && eventTime <  160 ) leff = 0;
      if( runNumber == 11269 && eventTime >  240                     ) leff = 0;
      
      if( runNumber == 11271 && eventTime <   10                     ) leff = 0;
      if( runNumber == 11271 && eventTime >   95 && eventTime <  140 ) leff = 0;
      if( runNumber == 11271 && eventTime >  195                     ) leff = 0;
      
      if( runNumber == 11272 && eventTime <   60                     ) leff = 0;
      if( runNumber == 11272 && eventTime >  160 && eventTime <  280 ) leff = 0;
      if( runNumber == 11272 && eventTime >  340                     ) leff = 0;
      
      if( runNumber == 11278 && eventTime <   10                     ) leff = 0;
      if( runNumber == 11278 && eventTime >   40 && eventTime <   90 ) leff = 0;
      if( runNumber == 11278 && eventTime >  140                     ) leff = 0;
      
      if( runNumber == 11279 && eventTime >   60 && eventTime <   90 ) leff = 0;
      if( runNumber == 11279 && eventTime >  120                     ) leff = 0;
      
      if( runNumber == 11280 && eventTime <   20                     ) leff = 0;
      if( runNumber == 11280 && eventTime >  200 && eventTime <  260 ) leff = 0;
      if( runNumber == 11280 && eventTime >  420                     ) leff = 0;
      
      if( runNumber == 11281 && eventTime <   60                     ) leff = 0;
      if( runNumber == 11281 && eventTime >  100 && eventTime <  150 ) leff = 0;
      if( runNumber == 11281 && eventTime >  420                     ) leff = 0;
      
      if( runNumber == 11282 && eventTime >   30 && eventTime <  100 ) leff = 0;
      if( runNumber == 11282 && eventTime >  160 && eventTime <  280 ) leff = 0;
      if( runNumber == 11282 && eventTime >  360 && eventTime <  420 ) leff = 0;
      
      if( runNumber == 11283 && eventTime <   30                     ) leff = 0;
      if( runNumber == 11283 && eventTime >  110 && eventTime <  220 ) leff = 0;
      if( runNumber == 11283 && eventTime >  260 && eventTime <  300 ) leff = 0;
      if( runNumber == 11283 && eventTime >  510                     ) leff = 0;
      
      if( runNumber == 11284 && eventTime <   20                     ) leff = 0;
      if( runNumber == 11284 && eventTime >  150 && eventTime <  200 ) leff = 0;
      if( runNumber == 11284 && eventTime >  340 && eventTime <  410 ) leff = 0;
      
      if( runNumber == 11285 && eventTime <   20                     ) leff = 0;
      if( runNumber == 11285 && eventTime >  190 && eventTime <  230 ) leff = 0;
      if( runNumber == 11285 && eventTime >  340 && eventTime <  390 ) leff = 0;
      if( runNumber == 11285 && eventTime >  420 && eventTime <  450 ) leff = 0;
      
      if( runNumber == 11286 && eventTime <   50                     ) leff = 0;
      if( runNumber == 11286 && eventTime >  260 && eventTime <  280 ) leff = 0;
      if( runNumber == 11286 && eventTime >  490 && eventTime <  520 ) leff = 0;
      
      if( runNumber == 11287 && eventTime <   30                     ) leff = 0;
      if( runNumber == 11287 && eventTime >  170 && eventTime <  190 ) leff = 0;
      if( runNumber == 11287 && eventTime >  270 && eventTime <  300 ) leff = 0;
      if( runNumber == 11287 && eventTime >  450 && eventTime <  530 ) leff = 0;
      
      if( runNumber == 11288 && eventTime <   10                     ) leff = 0;
      if( runNumber == 11288 && eventTime >   90 && eventTime <  110 ) leff = 0;
      if( runNumber == 11288 && eventTime >  200 && eventTime <  230 ) leff = 0;
      if( runNumber == 11288 && eventTime >  390 && eventTime <  410 ) leff = 0;
      if( runNumber == 11288 && eventTime >  500 && eventTime <  530 ) leff = 0;
      
      if( runNumber == 11289 && eventTime >  110 && eventTime <  200 ) leff = 0;
      if( runNumber == 11289 && eventTime >  400 && eventTime <  420 ) leff = 0;
      if( runNumber == 11289 && eventTime >  540 && eventTime <  560 ) leff = 0;
      
      if(leff)
	cmstimingcut->fill(eventTime);

      // CMS pixel clusters:

      for( std::vector<cluster>::iterator c = ClustDUT.begin(); c != ClustDUT.end(); c++ ){

	if(ETHh || FPIX ) {
	  // for ETH/KIT orientation
	  cmsxxHisto->fill( c->row, xA ); // anti-correlation x-x
	  cmsyyHisto->fill( c->col, yA ); // correlation y-y
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

	  qcol[px->col] += fabs(px->vcal); // project cluster onto cols
	  qrow[px->row] += fabs(px->vcal); // project cluster onto rows
	  if( px->col < colmin ) colmin = px->col;
	  if( px->col > colmax ) colmax = px->col;
	  if( px->row < rowmin ) rowmin = px->row;
	  if( px->row > rowmax ) rowmax = px->row;

	}//pix

	bool fiducial = 1;
	if( rowmin ==  0 ) fiducial = 0;
	if( rowmax == 79 ) fiducial = 0;
	if( colmin ==  0 ) fiducial = 0;
	if( colmax == 51 ) fiducial = 0;

	bool lowClusterCharge = false;
	if(c->charge < 8)
	  lowClusterCharge = true;

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

	// head-tail in cols:

	if( ncol > 2 ) {
	  double colmid = colmin + 0.5*(ncol-1);
	  double qht  = qcol[colmin] + qcol[colmax];
	  double asy = ( qcol[colmax] - qcol[colmin] ) / qht; // -1..1
	  double range = 0.5; // range 0..1. 0.5 is best at 46 deg turn
	  c->col = colmid + range * 0.5*(ncol-1) * asy;
	}

	// lq: Cut on cluster charge, checking whether lies inside the Landau peak
	bool lq = 0;
	double Q0 = c->charge * norm; // cluster charge normalized to vertical incidence
	if( Q0 > 18 &&  Q0 < 35 ) lq = 1;

	// DUT - triplet:

	double cmsx = ( c->col - 26 ) * pitchcol; // -3.9..3.9 mm
	double cmsy = ( c->row - 40 ) * pitchrow; // -4..4 mm

	if( FPIX ) {
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

	double cmsdx = dx4;
	double cmsdy = dy4;

	if( FPIX ) {
	  // FPIX: resid = sums
	  cmsdx = cmsx + DUTrot*cmsy - DUTalignx + xA; // residual x
	  cmsdy = cmsy - DUTrot*cmsx - DUTaligny + yA; // residual y
	}
	
	//if( ETHh ) {
	//  cmsdx = cmsx + DUTrot*cmsy - DUTalignx - xA; // residual x
	//  cmsdy = cmsy - DUTrot*cmsx - DUTaligny - yA; // residual y
	//}

	if( _DUT_chip >= 200 ) { // even/odd col effect for dig
	  int iodd = static_cast<int>(floor( fmod( c->col, 2 ) ));
	  if( iodd ) // odd
	    cmsdy = cmsdy + 1.5E-3;
	  else
	    cmsdy = cmsdy - 1.5E-3;
	}

	if( leff ){
	  cmsdxHisto->fill( cmsdx*1E3 );
	  cmsdyHisto->fill( cmsdy*1E3 );
	}

	if( leff && fiducial ) {

	  cmsdxfHisto->fill( cmsdx*1E3 );
	  cmsdyfHisto->fill( cmsdy*1E3 );

	  if( abs( cmsdy ) < 0.10 ) cmsdxfcHisto->fill( cmsdx*1E3 );
	  if( abs( cmsdx ) < 0.15 ) cmsdyfcHisto->fill( cmsdy*1E3 );
	  

	}//CMS fiducial

	 // accumulate cuts for y:

	if( leff && fiducial && abs( cmsdx ) < 0.15 ) {

	  if(      nrow == 1 )
	    cmsdyfc1Histo->fill( cmsdy*1E3 ); // 3972: 7.7
	  else if( nrow == 2 )
	    cmsdyfc2Histo->fill( cmsdy*1E3 ); // 3972: 9.5
	  else
	    cmsdyfc3Histo->fill( cmsdy*1E3 ); // 3872: 68

	  if(      Q0 < 18 )
	    cmsdyq0Histo->fill( cmsdy*1E3 );

	  else if( Q0 < 40 ){
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

	  if( abs( ty-0.000 ) < 0.002 &&  abs( tx-0.000 ) < 0.002 ) {
	    cmsdyfctHisto->fill( cmsdy*1E3 );
	    if(lowClusterCharge)
	      cmsdyfctLowChargeHisto->fill( cmsdy*1E3 );
	    if( nrow <= 2 ) cmsdyfcntHisto->fill( cmsdy*1E3 );
	  }

	  if( abs( ty-0.000 ) < 0.002 &&
	      abs( tx-0.000 ) < 0.002 &&
	      Q0 > 18 ) {

	    cmsdyfctqHisto->fill( cmsdy*1E3 ); // 7.8 um @ 4 GeV, 19 deg
	    if( nrow <= 2 ) cmsdyfcntqHisto->fill( cmsdy*1E3 );

	    if( Q0 < 40 ) {
	      cmsdyfctq1Histo->fill( cmsdy*1E3 ); // 7.8 um @ 4 GeV, 19 deg, more Gaussian
	      if( nrow <= 2 ) cmsdyfcntq1Histo->fill( cmsdy*1E3 );
	      if( c->col < 26 )
		cmsdyfctq1lHisto->fill( cmsdy*1E3 ); // xdb
	      else
		cmsdyfctq1rHisto->fill( cmsdy*1E3 ); // xdb
	    }
	    if( Q0 < 35 ) {
	      cmsdyfctq2Histo->fill( cmsdy*1E3 ); // inserted 26.12.2012
	    }
	    if( Q0 < 30 ) {
	      cmsdyfctq3Histo->fill( cmsdy*1E3 ); // was fctq2. 7.4 um @ 4 GeV, 19 deg
	      if( ldot ) 
		cmsdyfctqdotHisto->fill( cmsdy*1E3 ); // 8.1 um in run 5234
	      else
		cmsdyfctq3dHisto->fill( cmsdy*1E3 ); // 7.2 um in run 5234
	    }
	  }

	} // CMS fiducial for y

	// accumulate cuts for x:

	if( leff &&  fiducial && abs( cmsdy ) < 0.100 ) {

	  if( abs( ty-0.000 ) < 0.002 &&  abs( tx-0.000 ) < 0.002 ){
	    cmsdxfctHisto->fill( cmsdx*1E3 );
	    if(lowClusterCharge)
	      cmsdxfctLowChargeHisto->fill( cmsdx*1E3 );
	  }

	  if( abs( ty-0.000 ) < 0.002 &&
	      abs( tx-0.000 ) < 0.002 &&
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
	if( leff &&  abs( cmsdx ) < 0.3 && abs( cmsdy ) < 0.3 ) {

	  n_matched_clusters_dut++;
	  correvt100Histo->fill(event->getEventNumber());
	  correvt300Histo->fill(event->getEventNumber());
	  correvt1000Histo->fill(event->getEventNumber());
	  correvt4000Histo->fill(event->getEventNumber());

	  cmscolHisto->fill( c->col );
	  cmsrowHisto->fill( c->row );

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
	    for( std::vector<CMSPixel::pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); px++ ){
	      cmspxqHisto->fill( px->vcal);
	      if( c->size == 2)
		cmspxqcl2Histo->fill( px->vcal);
	      if( nrow == 2)
		cmspxqrow2Histo->fill( px->vcal);
	      if(px->vcal > qseed)
		qseed =  px->vcal;
	    }
	    cmsqseedfHisto->fill(qseed);

	    cmsqfHisto->fill( c->charge );
	    cmsq0fHisto->fill( Q0 );

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
	    cmsqvsxmym->fill( xmod, ymod, c->charge ); // cluster charge profile
	    
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

	    cmsqvsddt->fill( dutddtns, c->charge ); // cluster charge vs phase
	    cmsqvst1->fill( (time_now_tlu-time_event0)/fTLU, c->charge ); // cluster charge vs time
	    cmsqvst2->fill( (time_now_tlu-time_event0)/fTLU, c->charge ); // cluster charge vs time
	    cmsqvst3->fill( (time_now_tlu-time_event0)/fTLU, c->charge ); // cluster charge vs time
	    cmsqvst4->fill( (time_now_tlu-time_event0)/fTLU, c->charge ); // cluster charge vs time

	    cmsrmsxvsq->fill( Q0, abs(cmsdx)*1E3 ); //resolution vs charge
	    cmsrmsyvsq->fill( Q0, abs(cmsdy)*1E3 ); //resolution vs charge

	    double pMoyal = exp( -exp( -( Q0 - 28 ) / 3.3 ) ); // fitMoyal
	    cmspMoyalvsq->fill( Q0, pMoyal );
	    cmspMoyalHisto->fill( pMoyal );
	    cmsrmsyvsp->fill( pMoyal, abs(cmsdy)*1E3 ); // resolution vs charge

	    if(lowClusterCharge)
	      cmspixvsxmymLowCharge->fill( xmod, ymod );

	    if( lq ) {

	      cmsdyvsxm->fill( xmod, cmsdy*1E3 );
	      cmsdyvsym->fill( ymod, cmsdy*1E3 );

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

	      cmsrmsxvsx->fill( xAt, abs(cmsdx)*1E3 ); //resolution across cols
	      cmsrmsyvsx->fill( xAt, abs(cmsdy)*1E3 ); //resolution across cols
	      cmsrmsxvsy->fill( yAt, abs(cmsdx)*1E3 ); //resolution across rows
	      cmsrmsyvsy->fill( yAt, abs(cmsdy)*1E3 ); //resolution across rows
	      cmsrmsxvsxm->fill( xmod, abs(cmsdx)*1E3 ); //resolution within pixel
	      cmsrmsyvsxm->fill( xmod, abs(cmsdy)*1E3 ); //resolution within pixel
	      cmsncolvsxm->fill( xmod, ncol );
	      cmsnrowvsxm->fill( xmod, nrow );
	      if( !ldot ) {
		cmsrmsxvsym->fill( ymod, abs(cmsdx)*1E3 ); //resolution within pixel
		cmsrmsyvsym->fill( ymod, abs(cmsdy)*1E3 ); //resolution within pixel
		cmsrmsyvsym3->fill( ymd3, abs(cmsdy)*1E3 ); //resolution within pixel
		cmsrmsyvsym6->fill( ymd6, abs(cmsdy)*1E3 ); //resolution within pixel
	      }
	      cmsrmsyvst->fill( (time_now_tlu-time_event0)/fTLU, abs(cmsdy)*1E3 ); //resolution vs time
	      cmsrmsyvsddt->fill( dutddtns, abs(cmsdy)*1E3 ); // flat

	      if( nrow <= 2 ) {
		cmsdyvseta->fill( eta, cmsdy*1E3 );
		cmsrmsyvseta->fill( eta, abs(cmsdy)*1E3 );
	      }

	      cmsnpxvsxmym->fill( xmod, ymod, c->size ); // cluster size map
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

	    } // q Landau peak

	  } // fiducial

	} // CMS - triplet match

	if( leff && abs( cmsdx ) < 0.5 && abs( cmsdy ) < 0.5 ){ // link to CMS
	  (*trip).linked_dut = true;
	  (*trip).cmsdx = cmsdx;
	  (*trip).cmsdy = cmsdy;
	  ntrilk++;

	}

	// DUT alignment using GBL track fitting:
	if( fiducial && abs( cmsdx ) < 0.300 && abs( cmsdy ) < 0.300 ) {
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

      if( FPIX) {
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

	if( abs(sx) < 0.300 && abs(dy) < 0.200 ) {
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

    if (FPIX || ETHh){
      xmod = fmod( 9.075 + yAt, 0.3 ) * 1E3; // [0,300] um, 2 pixel wide
      ymod = fmod( 9.030 + xAt, 0.2 ) * 1E3; // [0,200] um
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
	  if( abs( turn ) > 11 ) resx = 22E-3;
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

    if( abs(dx) < 0.1 && abs(dy) < 0.1 && drip.linked_dut ) { // six with link

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

      bool leff = 1;

      if( leff ) effxyHisto->fill( xAt, yAt );
      if( leff ) effvsxy->fill( xAt, yAt, nm ); // CMS DUT efficiency profile

      if( leff && yAt > fidymin && yAt < fidymax ) {
	effvsx->fill( xAt, nm ); // CMS DUT efficiency profile
	if( probchi > 0.01 ) effvsxg->fill( xAt, nm );
      }

      if( leff && xAt > fidxmin && xAt < fidxmax ) {
	effvsy->fill( yAt, nm ); // CMS DUT efficiency profile
      }

      if( xAt > fidxmin && xAt < fidxmax && yAt > fidymin && yAt < fidymax ) {
	eff300->fill( (time_now_tlu-time_event0)/fTLU, nm );
	eff600->fill( (time_now_tlu-time_event0)/fTLU, nm );
	eff1200->fill( (time_now_tlu-time_event0)/fTLU, nm );
	eff1800->fill( (time_now_tlu-time_event0)/fTLU, nm );
	eff3600->fill( (time_now_tlu-time_event0)/fTLU, nm );

	if( leff ) { 
        effvsddt->fill( dutddtns, nm );
	    effvsxmym->fill( xmod, ymod, nm ); // CMS DUT efficiency profile
	    
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
	  }
	if( prevdutrefddt == 0 && dutrefddt == 0 )
	  effd600->fill( (time_now_tlu-time_event0)/fTLU, nm );
	if( downstream_triplets->size() < 3 ) effn600->fill( (time_now_tlu-time_event0)/fTLU, nm );
	if( downstream_triplets->size() < 2 ) effm600->fill( (time_now_tlu-time_event0)/fTLU, nm );

      }

    }// triplet-driplet match

     //------------------------------------------------------------------------
     // eff(REF) with DUT as timing plane:

    if( abs(dx) < 0.1 && abs(dy) < 0.1 && trip.linked_dut ) {

      bool nm = 0;
      if( drip.linked_dut ) nm = 1;

      rffvsxy->fill( -xR, -yR, nm ); // CMS REF efficiency profile

      if( abs( yR ) < 3 ) {
	rffvsx->fill( -xR, nm ); // CMS REF efficiency profile
      }

    }// triplet-driplet match

     //------------------------------------------------------------------------
     // intersect point in z:
    hit intersect = (*tr).intersect();
    double zx = intersect.x;
    double zy = intersect.y;

    if( abs(dy) < 0.1 ) { // no cut on dx
      if( abs( kx ) > 0.003 ) { // 
	sixzx3Histo->fill( zx - _planePosition[2] );
      }
      if( abs( kx ) > 0.002 ) { // 
	sixzx2Histo->fill( zx - _planePosition[2] );
      }
    }

    if( abs(dx) < 0.1 ) { // no cut on dy
      if( abs( ky ) > 0.003 ) { // 
	sixzy3Histo->fill( zy - _planePosition[2] );
      }
      if( abs( ky ) > 0.002 ) { // 
	sixzy2Histo->fill( zy - _planePosition[2] );
      }
    }

    //------------------------------------------------------------------------
    // z intersect:

    if( abs(dx) < 0.2 && abs(dy) < 0.2 ) { // looser cut allows more z range

      // measure scattering angle x after cuts in y:
      // cut on ky creates bias in kx

      if( abs( kx ) > 0.001 ) {
	sixzx1Histo->fill( zx - _planePosition[2] );
	if( abs( zx - DUTz ) < 30 ) {
	  sixkyzxHisto->fill( ky*1E3 );
	  sixkxzxHisto->fill( kx*1E3 ); // plot with gap, fittp0g.C("sixkxzx")
	}
      }

      if( abs( ky ) > 0.001 ) {
	sixzy1Histo->fill( zy - _planePosition[2] );
	if( abs( zy - DUTz ) < 30 ) {
	  sixkxzyHisto->fill( kx*1E3 );
	  sixkyzyHisto->fill( ky*1E3 ); // plot with gap
	}
      }

    }//match

     //------------------------------------------------------------------------

  } // Loop over found tracks

  nsixHisto->fill( telescope_tracks->size() );

  prevdutrefddt = dutrefddt;

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
    << "---------------------------------------------------------------------------------------------------------" 
    << std::endl << "DUT Data Decoding Statistics:" << std::endl
    << "---------------------------------------------------------------------------------------------------------"
    << std::endl;
  streamlog_out(MESSAGE5) << _dutDecoder->statistics.get() << std::endl;
  delete _dutDecoder;

  streamlog_out(MESSAGE5)
    << "---------------------------------------------------------------------------------------------------------" 
    << std::endl << "REF Data Decoding Statistics:" << std::endl
    << "---------------------------------------------------------------------------------------------------------"
    << std::endl;
  streamlog_out(MESSAGE5) << _refDecoder->statistics.get() << std::endl;
  delete _refDecoder;


  streamlog_out(MESSAGE5)
    << "---------------------------------------------------------------------------------------------------------" << std::endl
    << std::endl
    << "Processed events:    "
    << std::setw(10) << std::setiosflags(std::ios::right)
    << _nEvt << std::resetiosflags(std::ios::right) << std::endl
    << "TEL re-syncs " << telnsync << std::endl << std::endl
    << "DUT re-syncs " << dutnsync << std::endl
    << "DUT clusters " << n_clusters_dut << std::endl
    << "DUT matches  " << n_matched_clusters_dut << " (" << ((double)n_matched_clusters_dut/n_clusters_dut)*100 << "%)" << std::endl << std::endl
    << "REF re-syncs " << refnsync << std::endl
    << "REF clusters " << n_clusters_ref << std::endl
    << "REF matches  " << n_matched_clusters_ref << " (" << ((double)n_matched_clusters_ref/n_clusters_ref)*100 << "%)" << std::endl
    << std::endl;

  // pre-align:
  streamlog_out(MESSAGE5) << "pre-align:" << std::endl;

  streamlog_out(MESSAGE5) << "DUT x   " << cmsdxaHisto->title() << std::endl;
  streamlog_out(MESSAGE5) << "Entries " << cmsdxaHisto->entries() << std::endl;
  streamlog_out(MESSAGE5) << "Maximum " << cmsdxaHisto->maxBinHeight() << std::endl;

  int nb = cmsdxaHisto->axis().bins();
  int nmax = 0;
  double DUTx = 0;
  for( int ib = 0; ib < nb; ++ib ) {
    if( cmsdxaHisto->binHeight(ib) > nmax ) {
      nmax = static_cast<int>(cmsdxaHisto->binHeight(ib));
      DUTx = cmsdxaHisto->binMean(ib);
    }
  }
  streamlog_out(MESSAGE5) << "DUTx    " << DUTx << std::endl;

  streamlog_out(MESSAGE5) << "DUT y   " << cmssyaHisto->title() << std::endl;
  streamlog_out(MESSAGE5) << "Entries " << cmssyaHisto->entries() << std::endl;
  streamlog_out(MESSAGE5) << "Maximum " << cmssyaHisto->maxBinHeight() << std::endl;
  nb = cmssyaHisto->axis().bins();
  nmax = 0;
  double DUTy = 0;
  for( int ib = 0; ib < nb; ++ib ) {
    if( cmssyaHisto->binHeight(ib) > nmax ) {
      nmax = static_cast<int>(cmssyaHisto->binHeight(ib));
      DUTy = cmssyaHisto->binMean(ib);
    }
  }
  streamlog_out(MESSAGE5) << "DUTy    " << DUTy << std::endl;

  streamlog_out(MESSAGE5) << "REF x   " << refdxaHisto->title() << std::endl;
  streamlog_out(MESSAGE5) << "Entries " << refdxaHisto->entries() << std::endl;
  streamlog_out(MESSAGE5) << "Maximum " << refdxaHisto->maxBinHeight() << std::endl;

  nb = refdxaHisto->axis().bins();
  nmax = 0;
  double REFx = 0;
  for( int ib = 0; ib < nb; ++ib ) {
    if( refdxaHisto->binHeight(ib) > nmax ) {
      nmax = static_cast<int>(refdxaHisto->binHeight(ib));
      REFx = refdxaHisto->binMean(ib);
    }
  }
  streamlog_out(MESSAGE5) << "REFx    " << REFx << std::endl;

  streamlog_out(MESSAGE5) << "REF y   " << refsyaHisto->title() << std::endl;
  streamlog_out(MESSAGE5) << "Entries " << refsyaHisto->entries() << std::endl;
  streamlog_out(MESSAGE5) << "Maximum " << refsyaHisto->maxBinHeight() << std::endl;
  nb = refsyaHisto->axis().bins();
  nmax = 0;
  double REFy = 0;
  for( int ib = 0; ib < nb; ++ib ) {
    if( refsyaHisto->binHeight(ib) > nmax ) {
      nmax = static_cast<int>(refsyaHisto->binHeight(ib));
      REFy = refsyaHisto->binMean(ib);
    }
  }
  streamlog_out(MESSAGE5) << "REFy    " << REFy << std::endl;

  streamlog_out(MESSAGE5) 
    << std::endl
    << "runlistPreAlign: "
    << _TEL_run
    << "," << _DATE_run
    << "," << _gearfile
    << "," <<_eBeam
    << "," << _DUT_run
    << "," << _DUT_chip
    << "," << _DUT_decoding_flags
    << "," << _DUT_board
    << "," << _DUT_gain
    << "," << _DUT_calibration_type
    << "," << _REF_run
    << "," << _REF_chip
    << "," << _REF_decoding_flags
    << "," << _REF_board
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
    << "," << _nSkipDUT
    << "," << _nSkipRef
    << "," << _nSkipTelescope
    << std::endl;
  streamlog_out(MESSAGE5) << std::endl;


  ofstream prealignrunfile;
  prealignrunfile.open("prelines-for-runlist.txt",ios::app);
  prealignrunfile << _TEL_run << "," << _DATE_run << "," << _gearfile << "," <<_eBeam << "," << _DUT_run << "," << _DUT_chip << "," << _DUT_decoding_flags << "," << _DUT_board << "," << _DUT_gain  << "," << _DUT_calibration_type << "," << _REF_run << "," << _REF_chip << "," << _REF_decoding_flags << "," << _REF_board << "," << _REF_gain << "," << _REF_calibration_type << "," << DUTx << "," << DUTy << "," << DUTz - _planePosition[2] << "," << tilt << "," << turn << "," << DUTrot << "," << REFx << "," << REFy << "," << _REFz << "," << _REFrot << "," << _nSkipDUT << "," << _nSkipRef << "," << _nSkipTelescope << endl;
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
  steerFile.open( "steerPede.txt" );

  if( steerFile.is_open() && ldut ) {

    steerFile << "! generated by EUTelTestFitter" << std::endl;
    steerFile << "Cfiles" << std::endl;
    steerFile << "mille.bin" << std::endl;
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
      std::string command = "pede  steerPede.txt";

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
	  streamlog_out(MESSAGE5) << std::resetiosflags(std::ios::floatfield) << std::setprecision(9) << "dx    " << alpar[1]*1E3 << " um" << std::endl;
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
	  streamlog_out(MESSAGE5)
	    << std::endl
	    << "for runlist.csv:" << std::endl
	    << "runlistFullAlign: "
	    << _TEL_run
	    << "," << _DATE_run
	    << "," << _gearfile
	    << "," << _eBeam
	    << "," << _DUT_run
	    << "," << _DUT_chip
	    << "," << _DUT_decoding_flags
	    << "," << _DUT_board
	    << "," << _DUT_gain
	    << "," << _DUT_calibration_type
	    << "," << _REF_run
	    << "," << _REF_chip
	    << "," << _REF_decoding_flags
	    << "," << _REF_board
	    << "," << _REF_gain
	    << "," << _REF_calibration_type
	    << "," << DUTalignx-alpar[1]
	    << "," << DUTaligny-alpar[2]
	    << "," << DUTz-alpar[6] - _planePosition[2]
	    << "," << tilt-alpar[4]*180/3.141592654
	    << "," << turn-alpar[5]*180/3.141592654
	    << "," << DUTrot-alpar[3]
	    << "," << _REFalignx
	    << "," << _REFaligny
	    << "," << _REFz
	    << "," << _REFrot
	    << "," << _nSkipDUT
	    << "," << _nSkipRef
	    << "," << _nSkipTelescope
	    << std::endl;

	  ofstream runfile;
	  runfile.open("lines-for-runlist.txt",ios::app);
	  runfile << _TEL_run << "," << _DATE_run << "," << _gearfile << "," << _eBeam << "," << _DUT_run << "," << _DUT_chip << "," << _DUT_decoding_flags << "," << _DUT_board << "," << _DUT_gain << "," << _DUT_calibration_type << "," << _REF_run << "," << _REF_chip << "," << _REF_decoding_flags << "," << _REF_board << "," << _REF_gain << "," << _REF_calibration_type << "," << DUTalignx-alpar[1] << "," << DUTaligny-alpar[2] << "," << DUTz-alpar[6] - _planePosition[2] << "," << tilt-alpar[4]*180/3.141592654 << "," << turn-alpar[5]*180/3.141592654 << "," << DUTrot-alpar[3] << "," << _REFalignx << "," << _REFaligny << "," << _REFz << "," << _REFrot << "," << _nSkipDUT << "," << _nSkipRef << "," << _nSkipTelescope << endl;
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


//------------------------------------------------------------------------------
void EUTelAnalysisCMSPixel::bookHistos()
{

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)


  // event time:
  AIDAProcessor::tree(this)->mkdir("Timing");

  t1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/t1", 125, 0, 1 );
  t1Histo->setTitle( "telescope event time;telescope event time [s];events/s" );

  t10Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/t10", 250, 0, 10 );
  t10Histo->setTitle( "telescope event time;telescope event time [s];events/s" );

  t100Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/t100", 100, 0, 100 );
  t100Histo->setTitle( "telescope event time;telescope event time [s];events/s" );

  t300Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/t300", 350, 0, 350 );
  t300Histo->setTitle( "telescope event time;telescope event time [s];events/s" );

  t600Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/t600", 700, 0, 700 );
  t600Histo->setTitle( "telescope event time;telescope event time [s];events/s" );

  t1000Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/t1000", 100, 0, 1000 );
  t1000Histo->setTitle( "telescope event time;telescope event time [s];events/10s" );

  t1800Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/t1800", 200, 0, 2000 );
  t1800Histo->setTitle( "telescope event time;telescope event time [s];events/10s" );

  t3600Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/t3600", 400, 0, 4000 );
  t3600Histo->setTitle( "telescope event time;telescope event time [s];events/10s" );

  dtHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/dt", 500, 0, 5000 );
  dtHisto->setTitle( "telescope time between events;telescope time between events [#mus];events" );

  dtmsHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/dtms", 500, 0, 100 );
  dtmsHisto->setTitle( "telescope time between events;telescope time between events [ms];events" );

  logdtHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/logdt", 500, -1, 4 );
  logdtHisto->setTitle( "telescope time between events;telescope time between events log_{10}(#Deltat [ms]);events" );

  logdtcmsHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/logdtcms", 500, -1, 4 );
  logdtcmsHisto->setTitle( "telescope time between events with DUT;telescope time between events log_{10}(#Deltat [ms]);events" );

  dtfvstau = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/dtfvstau", 600, 372, 378, 0, 1 );
  dtfvstau->setTitle( "dt/tau fractional part;tau [ticks];dt/tau [fractional part]" );

  tfvstau = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/tfvstau", 600, 372, 378, 0, 1 );
  tfvstau->setTitle( "t/tau fractional part;tau [ticks];t/tau [fractional part]" );

  dtfHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/dtf", 500, -0.5, 0.5 );
  dtfHisto->setTitle( "phase of telescope time between events;phase of time between events [turns];telescope events" );

  dtfvst = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "dtfvst", 300, 0, 300, -1, 1 );
  dtfvst->setTitle( "dt phase vs t;run time [s];<TLU phase> [turns];" );

  dtfvsdt = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/dtfvsdt", 100, 0, 10000, -1, 1 );
  dtfvsdt->setTitle( "dt phase vs dt;delta time [us];<TLU phase> [turns];" );

  tfHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/tf", 100, -0.5, 0.5 );
  tfHisto->setTitle( "phase of telescope event time;phase of event time [turns];telescope events" );

  tfvst = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/tfvst", 300, 0, 300E6, -1, 1 );
  tfvst->setTitle( "t phase vs t;turns;<TLU phase> [turns];" );

  tfvst1 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/tfvst1", 100, 0, 1E6, -1, 1 );
  tfvst1->setTitle( "t phase vs t;turns;<TLU phase> [turns];" );

  tfvst10 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/tfvst10", 1000, 0, 10E6, -1, 1 );
  tfvst10->setTitle( "t phase vs t;turns;<TLU phase> [turns];" );

  tfvst100 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/tfvst100", 1000, 0, 100E6, -1, 1 );
  tfvst100->setTitle( "t phase vs t;turns;<TLU phase> [turns];" );

  tfvst300 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/tfvst300", 1000, 0, 300E6, -1, 1 );
  tfvst300->setTitle( "t phase vs t;turns;<TLU phase> [turns];" );



  // DUT clusters:

  dutnclusHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "dutnclus", 11, -0.5, 11.5 );
  dutnclusHisto->setTitle( "DUT clusters/event;DUT clusters/event;events" );

  dutcolHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "dutcol", 52, -0.5, 51.5 );
  dutcolHisto->setTitle( "DUT column;DUT cluster col;DUT clusters" );

  dutrowHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "dutrow", 80, -0.5, 79.5 );
  dutrowHisto->setTitle( "DUT row;DUT cluster row;DUT clusters" );

  dutnpxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "dutnpx", 11, -0.5, 10.5 );
  dutnpxHisto->setTitle( "DUT cluster size;DUT pixel per cluster;DUT clusters" );

  dutadcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "dutadc", 100, 0, 100 );
  dutadcHisto->setTitle( "DUT cluster charge;DUT cluster charge [ke];DUT clusters" );


  refnclusHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refnclus", 11, -0.5, 11.5 );
  refnclusHisto->setTitle( "REF clusters/event;REF clusters/event;events" );

  refcolHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refcol", 52, -0.5, 51.5 );
  refcolHisto->setTitle( "REF column;REF cluster col;REF clusters" );

  refrowHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refrow", 80, -0.5, 79.5 );
  refrowHisto->setTitle( "REF row;REF cluster row;REF clusters" );

  refnpxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refnpx", 11, -0.5, 10.5 );
  refnpxHisto->setTitle( "REF cluster size;REF pixel per cluster;REF clusters" );

  refadcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refadc", 100, 0, 100 );
  refadcHisto->setTitle( "REF cluster charge;REF cluster charge [ke];REF clusters" );

  // telescope hits per plane:
  AIDAProcessor::tree(this)->mkdir("Telescope");

  nAllHitHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/nallhit", 201, -0.5, 200.5 );
  nAllHitHisto->setTitle( "Telescope hits/event;telescope hits;events" );

  // telescope dx:

  dx01Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx01", 100, -1, 1 );
  dx01Histo->setTitle( "x1-x0;x_{1}-x_{0} [mm];hit pairs" );

  dy01Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dy01", 100, -1, 1 );
  dy01Histo->setTitle( "y1-y0;y_{1}-y_{0} [mm];hit pairs" );

  du01Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/du01", 100, -1, 1 );
  du01Histo->setTitle( "x1-x0, |dy| < 1;x_{1}-x_{0} [mm];hit pairs" );

  dx02Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx02", 100, -1, 1 );
  dx02Histo->setTitle( "x2-x0;x_{2}-x_{0} [mm];hit pairs" );

  dx03Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx03", 100, -1, 1 );
  dx03Histo->setTitle( "x3-x0;x_{3}-x_{0} [mm];hit pairs" );

  dx04Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx04", 100, -1, 1 );
  dx04Histo->setTitle( "x4-x0;x_{4}-x_{0} [mm];hit pairs" );

  dx05Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx05", 100, -1, 1 );
  dx05Histo->setTitle( "x5-x0;x_{5}-x_{0} [mm];hit pairs" );

  dx12Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx12", 100, -1, 1 );
  dx12Histo->setTitle( "x2-x1;x_{2}-x_{1} [mm];hit pairs" );

  dy12Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dy12", 100, -1, 1 );
  dy12Histo->setTitle( "y2-y1;y_{2}-y_{1} [mm];hit pairs" );

  du12Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/du12", 100, -1, 1 );
  du12Histo->setTitle( "x2-x1, |dy| < 1;x_{2}-x_{1} [mm];hit pairs" );

  dx23Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx23", 100, -1, 1 );
  dx23Histo->setTitle( "x3-x2;x_{3}-x_{2} [mm];hit pairs" );

  dy23Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dy23", 100, -1, 1 );
  dy23Histo->setTitle( "y3-y2;y_{3}-y_{2} [mm];hit pairs" );

  du23Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/du23", 100, -1, 1 );
  du23Histo->setTitle( "x3-x2, |dy| < 1;x_{3}-x_{2} [mm];hit pairs" );

  dx34Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx34", 100, -1, 1 );
  dx34Histo->setTitle( "x4-x3;x_{4}-x_{3} [mm];hit pairs" );

  dy34Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dy34", 100, -1, 1 );
  dy34Histo->setTitle( "y4-y3;y_{4}-y_{3} [mm];hit pairs" );

  du34Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/du34", 100, -1, 1 );
  du34Histo->setTitle( "x4-x3, |dy| < 1;x_{4}-x_{3} [mm];hit pairs" );

  dx45Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx45", 100, -1, 1 );
  dx45Histo->setTitle( "x5-x4;x_{5}-x_{4} [mm];hit pairs" );

  dy45Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dy45", 100, -1, 1 );
  dy45Histo->setTitle( "y5-y4;y_{5}-y_{4} [mm];hit pairs" );

  du45Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/du45", 100, -1, 1 );
  du45Histo->setTitle( "x5-x4, |dy| < 1;x_{5}-x_{4} [mm];hit pairs" );

  // triplets:
  AIDAProcessor::tree(this)->mkdir("Upstream");

  dzcvsxy = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "Upstream/dzcvsxy", 120, -12, 12, 60, -6, 6, -999, 999 );
  dzcvsxy->setTitle( "DUT plane;telescope track x_{DUT} [mm];telescope track y_{DUT} [mm];z_{DUT} [mm]" );

  z3vsxy = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "Upstream/z3vsxy", 120, -12, 12, 60, -6, 6, -999, 999 );
  z3vsxy->setTitle( "DUT plane;telescope track x_{DUT} [mm];telescope track y_{DUT} [mm];z_{DUT} [mm]" );

  tridxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx", 100, -100, 100 );
  tridxHisto->setTitle( "triplet dx;x_{1}-x_{m} [#mum];telescope triplets" );

  tridyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy", 100, -100, 100 );
  tridyHisto->setTitle( "triplet dy;y_{1}-y_{m} [#mum];telescope triplets" );

  tridx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx1", 100, -1, 1 );
  tridx1Histo->setTitle( "triplet dx;x_{1}-x_{t} [mm];telescope triplets" );

  tridy1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy1", 100, -1, 1 );
  tridy1Histo->setTitle( "triplet dy;y_{1}-y_{t} [mm];telescope triplets" );

  tridxvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridxvsx", 100, -10, 10, -100, 100 );
  tridxvsx->setTitle( "triplet x resid vs x;x [mm];triplet <#Deltax> [#mum]" );

  tridxvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridxvsy", 50, -5, 5, -100, 100 );
  tridxvsy->setTitle( "triplet x resid vs y;y [mm];triplet <#Deltax> [#mum]" );

  tridxvstx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridxvstx", 80, -2, 2, -100, 100 );
  tridxvstx->setTitle( "triplet x resid vs tx;t_{x} [mrad];triplet <#Deltax> [#mum]" );

  tridxvsty = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridxvsty", 80, -2, 2, -100, 100 );
  tridxvsty->setTitle( "triplet x resid vs ty;t_{y} [mrad];triplet <#Deltax> [#mum]" );

  tridyvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridyvsx", 100, -10, 10, -100, 100 );
  tridyvsx->setTitle( "triplet y resid vs x;x [mm];triplet <#Deltay> [#mum]" );

  tridyvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridyvsy", 50, -5, 5, -100, 100 );
  tridyvsy->setTitle( "triplet y resid vs y;y [mm];triplet <#Deltay> [#mum]" );

  tridyvstx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridyvstx", 80, -2, 2, -100, 100 );
  tridyvstx->setTitle( "triplet y resid vs tx;t_{x} [mrad];triplet <#Deltay> [#mum]" );

  tridyvsty = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridyvsty", 80, -2, 2, -100, 100 );
  tridyvsty->setTitle( "triplet y resid vs ty;t_{y} [mrad];triplet <#Deltay> [#mum]" );

  tridx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx3", 100, -1, 1 );
  tridx3Histo->setTitle( "triplet dx;x_{3}-x_{t} [mm];telescope triplets" );

  tridy3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy3", 100, -1, 1 );
  tridy3Histo->setTitle( "triplet dy;y_{3}-y_{t} [mm];telescope triplets" );

  tridx3bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx3b", 100, -250, 250 );
  tridx3bHisto->setTitle( "triplet dx;x_{3}-x_{t} [um];telescope triplets" );

  tridy3bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy3b", 100, -250, 250 );
  tridy3bHisto->setTitle( "triplet dy;y_{3}-y_{t} [um];telescope triplets" );

  tridx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx4", 100, -2, 2 );
  tridx4Histo->setTitle( "triplet dx;x_{4}-x_{t} [mm];telescope triplets" );

  tridy4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy4", 100, -2, 2 );
  tridy4Histo->setTitle( "triplet dy;y_{4}-y_{t} [mm];telescope triplets" );

  tridx4bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx4b", 100, -400, 400 );
  tridx4bHisto->setTitle( "triplet dx;x_{4}-x_{t} [um];telescope triplets" );

  tridy4bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy4b", 100, -400, 400 );
  tridy4bHisto->setTitle( "triplet dy;y_{4}-y_{t} [um];telescope triplets" );

  tridx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx5", 100, -3, 3 );
  tridx5Histo->setTitle( "triplet dx;x_{5}-x_{t} [mm];telescope triplets" );

  tridy5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy5", 100, -3, 3 );
  tridy5Histo->setTitle( "triplet dy;y_{5}-y_{t} [mm];telescope triplets" );

  tridx5bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx5b", 100, -1000, 1000 );
  tridx5bHisto->setTitle( "triplet dx;x_{5}-x_{t} [um];telescope triplets" );

  tridy5bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy5b", 100, -1000, 1000 );
  tridy5bHisto->setTitle( "triplet dy;y_{5}-y_{t} [um];telescope triplets" );

  trixHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/trix", 240, -12, 12 );
  trixHisto->setTitle( "triplet x1;x1_{out} [mm];telescope triplets" );

  triyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/triy", 120, -6, 6 );
  triyHisto->setTitle( "triplet y1;y1_{up} [mm];telescope triplets" );

  trixyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "Upstream/trixy", 240, -12, 12, 120, -6, 6 );
  trixyHisto->setTitle( "triplet y1 vs x1;x1_{out} [mm];y1_{up} [mm];telescope triplets" );

  tritxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tritx", 100, -10, 10 );
  tritxHisto->setTitle( "triplet slope x;#theta_{x} [mrad];telescope triplets" );

  trityHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/trity", 100, -10, 10 );
  trityHisto->setTitle( "triplet slope y;#theta_{y} [mrad];telescope triplets" );

  trixdutHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/trixdut", 240, -12, 12 );
  trixdutHisto->setTitle( "triplet at DUT;triplet x_{out} at DUT [mm];telescope triplets" );

  triydutHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/triydut", 120, -6, 6 );
  triydutHisto->setTitle( "triplet at DUT;triplet y_{up} at DUT [mm];telescope triplets" );

  trixydutHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "Upstream/trixydut", 240, -12, 12, 120, -6, 6 );
  trixydutHisto->setTitle( "triplet at DUT;triplet x_{out} at DUT [mm];triplet y_{up} at DUT [mm];telescope triplets" );

  // DUT pixel vs triplets:

  cmstimingcut = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmstimingcut", 140, 0, 700 );
  cmsnpxHisto->setTitle( "check if timing cut was applied;time[s]" );
  

  //FIXME cmsxx and cmsyy need swapped dimensions for ETHh and FPIX:
  cmsxxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "cmsxx", 52, -0.5, 51.5, 110, -11, 11 );
  cmsxxHisto->setTitle( "x correlation;DUT cluster col;telescope triplet x [mm];clusters" );

  cmsyyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "cmsyy", 80, -0.5, 79.5, 55, -5.5, 5.5 );
  cmsyyHisto->setTitle( "y correlation;DUT cluster row;telescope triplet y [mm];clusters" );

  cmspxqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmspxq", 100, 0, 25 );
  cmspxqHisto->setTitle( "DUT pixel charge linked;DUT pixel charge [ke];DUT linked pixels" );

  cmspxqcl2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmspxqcl2", 100, 0, 25 );
  cmspxqcl2Histo->setTitle( "DUT pixel charge linked 2 pixel cluster;DUT pixel charge [ke];DUT linked pixels" );

  cmspxqrow2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmspxqrow2", 100, 0, 25 );
  cmspxqrow2Histo->setTitle( "DUT pixel charge linked 2 row cluster;DUT pixel charge [ke];DUT linked pixels" );

  cmssxaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmssxa", 440, -11, 11 );
  cmssxaHisto->setTitle( "Pixel + Telescope x;cluster + triplet #Sigmax [mm];clusters" );

  cmsdyaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdya", 220, -5.5, 5.5 );
  cmsdyaHisto->setTitle( "Pixel - telescope y;cluster - triplet #Deltay [mm];clusters" );

  cmsdxaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdxa", 440, -11, 11 );
  cmsdxaHisto->setTitle( "Pixel - Telescope x;cluster - triplet #Deltax [mm];clusters" );

  cmssyaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmssya", 220, -5.5, 5.5 );
  cmssyaHisto->setTitle( "Pixel + telescope y;cluster + triplet #Sigmay [mm];clusters" );

  cmsdx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdx4", 440, -11, 11 );
  cmsdx4Histo->setTitle( "Pixel - Telescope x;cluster + triplet #Sigmax [mm];clusters" );

  cmsdy4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdy4", 220, -5.5, 5.5 );
  cmsdy4Histo->setTitle( "Pixel - telescope y;cluster - triplet #Deltay [mm];clusters" );

  cmsdxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdx", 200, -500, 500 );
  cmsdxHisto->setTitle( "Pixel + Telescope x;cluster + triplet #Sigmax [#mum];clusters" );

  cmsdyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdy", 500, -500, 500 );
  cmsdyHisto->setTitle( "Pixel - telescope y;cluster - triplet #Deltay [#mum];clusters" );

  cmsdxfHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdxf", 200, -500, 500 );
  cmsdxfHisto->setTitle( "fiducial Pixel - telescope x;fiducial cluster - triplet #Deltax [#mum];fiducial clusters" );

  cmsdyfHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyf", 500, -500, 500 );
  cmsdyfHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdxfcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdxfc", 200, -500, 500 );
  cmsdxfcHisto->setTitle( "fiducial Pixel - telescope x;fiducial cluster - triplet #Deltax [#mum];fiducial clusters" );

  cmsdyfcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfc", 500, -500, 500 );
  cmsdyfcHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfc1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfc1", 500, -500, 500 );
  cmsdyfc1Histo->setTitle( "Pixel - telescope y;1-row cluster - triplet #Deltay [#mum];fiducial 1-row clusters" );

  cmsdyfc2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfc2", 500, -500, 500 );
  cmsdyfc2Histo->setTitle( "Pixel - telescope y;2-row cluster - triplet #Deltay [#mum];fiducial 2-row clusters" );

  cmsdyfc3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfc3", 500, -500, 500 );
  cmsdyfc3Histo->setTitle( "Pixel - telescope y;3-row cluster - triplet #Deltay [#mum];fiducial 3-row clusters" );

  cmsdyq0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyq0", 500, -500, 500 );
  cmsdyq0Histo->setTitle( "Pixel - telescope y, Q < 18;cluster - triplet #Deltay [#mum];low dE/dx clusters" );

  cmsdyq1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyq1", 500, -500, 500 );
  cmsdyq1Histo->setTitle( "Pixel - telescope y, 18 < Q < 40;cluster - triplet #Deltay [#mum];med dE/dx clusters" );

  cmsdyeta0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyeta0", 500, -500, 500 );
  cmsdyeta0Histo->setTitle( "Pixel - telescope y, neg eta;cluster - triplet #Deltay [#mum];neg eta clusters" );

  cmsdyeta1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyeta1", 500, -500, 500 );
  cmsdyeta1Histo->setTitle( "Pixel - telescope y, pos eta;cluster - triplet #Deltay [#mum];pos eta clusters" );

  cmsdyq2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyq2", 500, -500, 500 );
  cmsdyq2Histo->setTitle( "Pixel - telescope y, Q > 40;cluster - triplet #Deltay [#mum];high dE/dx clusters" );

  cmsdxfctHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdxfct", 500, -500, 500 );
  cmsdxfctHisto->setTitle( "fiducial Pixel - telescope x;fiducial cluster - triplet #Deltax [#mum];fiducial clusters" );

  cmsdxfctLowChargeHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdxfctLowCharge", 500, -500, 500 );
  cmsdxfctLowChargeHisto->setTitle( "fiducial Pixel - telescope x;fiducial cluster - triplet #Deltax [#mum];fiducial clusters" );

  cmsdyfctHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfct", 500, -500, 500 );
  cmsdyfctHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfctLowChargeHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctLowCharge", 500, -500, 500 );
  cmsdyfctLowChargeHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfcntHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfcnt", 500, -500, 500 );
  cmsdyfcntHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial 1-,2-row clusters" );

  cmsdxfctqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdxfctq", 500, -500, 500 );
  cmsdxfctqHisto->setTitle( "fiducial Pixel - telescope x;fiducial cluster - triplet #Deltax [#mum];fiducial clusters" );

  cmsdyfctqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctq", 500, -500, 500 );
  cmsdyfctqHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfcntqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfcntq", 500, -500, 500 );
  cmsdyfcntqHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial 1-,2-row clusters" );

  cmsdxfctq1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdxfctq1", 500, -500, 500 );
  cmsdxfctq1Histo->setTitle( "fiducial Pixel - telescope x;fiducial cluster - triplet #Deltax [#mum];fiducial clusters" );

  cmsdyfctq1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctq1", 500, -500, 500 );
  cmsdyfctq1Histo->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfcntq1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfcntq1", 500, -500, 500 );
  cmsdyfcntq1Histo->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial 1-,2-row clusters" );

  cmsdyfctq1lHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctq1l", 500, -500, 500 );
  cmsdyfctq1lHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfctq1rHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctq1r", 500, -500, 500 );
  cmsdyfctq1rHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdxfctq2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdxfctq2", 500, -500, 500 );
  cmsdxfctq2Histo->setTitle( "fiducial Pixel - telescope x;fiducial cluster - triplet #Deltax [#mum];fiducial clusters" );

  cmsdyfctq2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctq2", 500, -500, 500 );
  cmsdyfctq2Histo->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdxfctq3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdxfctq3", 500, -500, 500 );
  cmsdxfctq3Histo->setTitle( "fiducial Pixel - telescope x;fiducial cluster - triplet #Deltax [#mum];fiducial clusters" );

  cmsdyfctq3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctq3", 500, -500, 500 );
  cmsdyfctq3Histo->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfctqdotHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctqdot", 500, -500, 500 );
  cmsdyfctqdotHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfctq3dHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctq3d", 500, -500, 500 );
  cmsdyfctq3dHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );


  cmscolHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmscol", 52, -0.5, 51.5 );
  cmscolHisto->setTitle( "DUT linked column;DUT linked cluster col;DUT linked clusters" );

  cmsrowHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsrow", 80, -0.5, 79.5 );
  cmsrowHisto->setTitle( "DUT linked row;DUT linked cluster row;DUT linked clusters" );

  cmsqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsq", 100, 0, 100 );
  cmsqHisto->setTitle( "DUT cluster charge linked;DUT cluster charge [ke];DUT linked clusters" );

  cmsq0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsq0", 100, 0, 100 );
  cmsq0Histo->setTitle( "DUT cluster charge linked;normal DUT cluster charge [ke];DUT linked clusters" );

  trixlkHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "trixlk", 240, -12, 12 );
  trixlkHisto->setTitle( "DUT - triplet link;triplet x_{out} [mm];linked clusters" );

  triylkHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "triylk", 120, -6, 6 );
  triylkHisto->setTitle( "DUT - triplet link;triplet y_{up} [mm];linked clusters" );

  trixylkHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "trixylk", 240, -12, 12, 120, -6, 6 );
  trixylkHisto->setTitle( "DUT - triplet link;triplet x_{out} [mm];triplet y_{up} [mm];linked clusters" );


  cmsdxvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdxvsx", 52, -26*0.15, 26*0.15, -150, 150 ); // 16.2.2013
  cmsdxvsx->setTitle( "DUT x resid vs x;x [mm];<DUT cluster - telescope triplet #Deltax> [#mum]" );

  cmsdyvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdyvsx", 52, -26*0.15, 26*0.15, -150, 150 );
  cmsdyvsx->setTitle( "DUT y resid vs x;x [mm];<DUT cluster - telescope triplet #Deltay> [#mum]" );

  cmsdxvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdxvsy", 80, -4, 4, -150, 150 );
  cmsdxvsy->setTitle( "DUT x resid vs y;y [mm];<DUT cluster - telescope triplet #Deltax> [#mum]" );

  cmsdyvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdyvsy", 80, -4, 4, -150, 150 );
  cmsdyvsy->setTitle( "DUT y resid vs y;y [mm];<DUT cluster - telescope triplet #Deltay> [#mum]" );

  cmsdxvstx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdxvstx", 80, -2, 2, -150, 150 );
  cmsdxvstx->setTitle( "DUT x resid vs tet x;#theta_{x} [mrad];<DUT cluster - telescope triplet #Deltax> [#mum]" );

  cmsdyvsty = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdyvsty", 80, -2, 2, -150, 150 );
  cmsdyvsty->setTitle( "DUT y resid vs tet y;#theta_{y} [mrad];<DUT cluster - telescope triplet #Deltay> [#mum]" );

  cmsdyvsxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "cmsdyvsxh", 52, -26*0.15, 26*0.15, 50, -250, 250 );
  cmsdyvsxHisto->setTitle( "DUT resolution;x [mm];#Deltay [#mum];clusters" );

  cmsnpxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsnpx", 11, -0.5, 10.5 );
  cmsnpxHisto->setTitle( "linked CMS cluster size;CMS pixel per linked cluster;linked CMS clusters" );

  cmsnpxLowChargeHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsnpxLowCharge", 11, -0.5, 10.5 );
  cmsnpxLowChargeHisto->setTitle( "linked low charge CMS cluster size;CMS pixel per linked cluster;linked CMS clusters" );

  cmsnpx0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsnpx0", 11, -0.5, 10.5 );
  cmsnpx0Histo->setTitle( "linked CMS cluster size, bias dot;CMS pixel per linked cluster;linked CMS clusters in bias dot" );

  cmsnpx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsnpx1", 11, -0.5, 10.5 );
  cmsnpx1Histo->setTitle( "linked CMS cluster size, no bias dot;CMS pixel per linked cluster;linked CMS clusters no bias dot" );

  cmsnpx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsnpx2", 11, -0.5, 10.5 );
  cmsnpx2Histo->setTitle( "linked CMS cluster size, no bias dot, core;CMS pixel per linked cluster;linked CMS clusters no bias dot, core" );

  cmsnpx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsnpx3", 11, -0.5, 10.5 );
  cmsnpx3Histo->setTitle( "linked CMS cluster size, no bias dot, edge;CMS pixel per linked cluster;linked CMS clusters no bias dot, edge" );

  cmsncolHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsncol", 11, -0.5, 10.5 );
  cmsncolHisto->setTitle( "DUT cluster col size;columns per cluster;DUT clusters" );

  cmsncolLowChargeHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsncolLowCharge", 11, -0.5, 10.5 );
  cmsncolLowChargeHisto->setTitle( "linked low charge CMS col size;CMS pixel per linked cluster;linked CMS clusters" );

  cmsnrowHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsnrow", 11, -0.5, 10.5 );
  cmsnrowHisto->setTitle( "DUT cluster row size;rows per cluster;DUT clusters" );

  cmsnrowLowChargeHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsnrowLowCharge", 11, -0.5, 10.5 );
  cmsnrowLowChargeHisto->setTitle( "linked low charge CMS row size;CMS pixel per linked cluster;linked CMS clusters" );

  cmsnrowqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsnrowq", 11, -0.5, 10.5 );
  cmsnrowqHisto->setTitle( "DUT Landau peak cluster row size;rows per cluster in Landau peak ;DUT Landau peak clusters" );

  cmsnrowvst1 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsnrowvst1", 70, 0, 350, 0, 5.5 );
  cmsnrowvst1->setTitle( "DUT cluster rows vs time;run time [s];<cluster rows>" );

  cmsetaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmseta", 100, -1, 1 );
  cmsetaHisto->setTitle( "DUT 2-row eta;DUT cluster 2-row eta;DUT 2-row clusters" );

  cmsetaevenHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsetaeven", 100, -1, 1 );
  cmsetaevenHisto->setTitle( "DUT 2-row eta (even columns);DUT cluster 2-row eta;DUT 2-row clusters" );

  cmsetaoddHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsetaodd", 100, -1, 1 );
  cmsetaoddHisto->setTitle( "DUT 2-row eta (odd columns);DUT cluster 2-row eta;DUT 2-row clusters" );

  cmsqfHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsqf", 100, 0, 100 );
  cmsqfHisto->setTitle( "DUT cluster charge linked fiducial;DUT cluster charge [ke];DUT linked fiducial clusters" );

  cmsqfcl1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsqfcl1", 100, 0, 100 );
  cmsqfcl1Histo->setTitle( "DUT 1 pixel cluster charge linked fiducial;DUT cluster charge [ke];DUT linked fiducial clusters" );

  cmsqfcl2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsqfcl2", 100, 0, 100 );
  cmsqfcl2Histo->setTitle( "DUT 2 pixel cluster charge linked fiducial;DUT cluster charge [ke];DUT linked fiducial clusters" );

  cmsqfrow1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsqfrow1", 100, 0, 100 );
  cmsqfrow1Histo->setTitle( "DUT 1 row cluster charge linked fiducial;DUT cluster charge [ke];DUT linked fiducial clusters" );

  cmsqfrow2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsqfrow2", 100, 0, 100 );
  cmsqfrow2Histo->setTitle( "DUT 2 row cluster charge linked fiducial;DUT cluster charge [ke];DUT linked fiducial clusters" );

  cmsq0fHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsq0f", 100, 0, 100 );
  cmsq0fHisto->setTitle( "DUT cluster charge linked fiducial;normal DUT cluster charge [ke];DUT linked fiducial clusters" );

  cmsqf0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsqf0", 100, 0, 100 );
  cmsqf0Histo->setTitle( "DUT cluster charge linked fiducial bias dot;DUT cluster charge [ke];DUT linked fiducial clusters bias dot" );

  cmsqf1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsqf1", 100, 0, 100 );
  cmsqf1Histo->setTitle( "DUT cluster charge linked fiducial no dot;DUT cluster charge [ke];DUT linked fiducial clusters no dot" );

  cmsqf2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsqf2", 100, 0, 100 );
  cmsqf2Histo->setTitle( "DUT cluster charge linked fiducial no dot core;DUT cluster charge [ke];DUT linked fiducial clusters no dot core" );

  cmsqf3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsqf3", 100, 0, 100 );
  cmsqf3Histo->setTitle( "DUT cluster charge linked fiducial no dot edge;DUT cluster charge [ke];DUT linked fiducial clusters no dot edge" );
  
  cmsqseedfHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsqseedf", 100, 0, 100 );
  cmsqseedfHisto->setTitle( "DUT cluster seed charge linked fiducial;DUT seed charge [ke];DUT linked fiducial clusters" );

  cmsdyvsxm = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdyvsxm", 60, 0, 300, -150, 150 );
  cmsdyvsxm->setTitle( "DUT y resid vs xmod;telescope x mod 300 [#mum];<DUT cluster - telescope triplet #Deltay> [#mum]" );

  cmsdyvsym = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdyvsym", 40, 0, 200, -150, 150 );
  cmsdyvsym->setTitle( "DUT y resid vs ymod;telescope y mod 200 [#mum];<DUT cluster - telescope triplet #Deltay> [#mum]" );

  cmsdxvsxm = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdxvsxm", 60, 0, 300, -150, 150 );
  cmsdxvsxm->setTitle( "DUT x resid vs xmod;telescope x mod 300 [#mum];<DUT cluster - telescope triplet #Deltax> [#mum]" );

  cmsdxvsym = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdxvsym", 40, 0, 200, -150, 150 );
  cmsdxvsym->setTitle( "DUT x resid vs ymod;telescope y mod 200 [#mum];<DUT cluster - telescope triplet #Deltax> [#mum]" );

  cmspixvsxmym = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "cmspixvsxmym", 60, 0, 300, 40, 0, 200 );
  cmspixvsxmym->setTitle( "DUT pixel occupancy;x mod 300 #mum;y mod 200 #mum;clusters" );

  cmspixvsxmymLowCharge = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "cmspixvsxmymLowCharge", 60, 0, 300, 40, 0, 200 );
  cmspixvsxmymLowCharge->setTitle( "DUT pixel occupancy;x mod 300 #mum;y mod 200 #mum;clusters" );

  // KIT: added for efficiency analysis - occupancy

  cmspixvsxm50 = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmspixvsxm50", 100, 0, 200 );
  cmspixvsxm50->setTitle( "DUT pixel occupancy vs ymod (xm50);telescope y_{DUT} mod at xm50 [#mum];clusters" );

  cmspixvsxm100 = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmspixvsxm100", 100, 0, 200 );
  cmspixvsxm100->setTitle( "DUT pixel occupancy vs ymod (xm100);telescope y_{DUT} mod at xm100 [#mum];clusters" );

  cmspixvsxm150 = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmspixvsxm150", 100, 0, 200 );
  cmspixvsxm150->setTitle( "DUT pixel occupancy vs ymod (xm150);telescope y_{DUT} mod at xm150 [#mum];clusters" );

  cmspixvsxm200 = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmspixvsxm200", 100, 0, 200 );
  cmspixvsxm200->setTitle( "DUT pixel occupancy vs ymod (xm200);telescope y_{DUT} mod at xm200 [#mum];clusters" );

  cmspixvsxm250 = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmspixvsxm250", 100, 0, 200 );
  cmspixvsxm250->setTitle( "DUT pixel occupancy vs ymod (xm250);telescope y_{DUT} mod at xm250 [#mum];clusters" );

  // ------

  cmspixvsym25 = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmspixvsym25", 150, 0, 300 );
  cmspixvsym25->setTitle( "DUT pixel occupancy vs xmod (ym25);telescope x_{DUT} mod at ym25 [#mum];clusters" );

  cmspixvsym50 = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmspixvsym50", 150, 0, 300 );
  cmspixvsym50->setTitle( "DUT pixel occupancy vs xmod (ym50);telescope x_{DUT} mod at ym50 [#mum];clusters" );

  cmspixvsym75 = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmspixvsym75", 150, 0, 300 );
  cmspixvsym75->setTitle( "DUT pixel occupancy vs xmod (ym75);telescope x_{DUT} mod at ym75 [#mum];clusters" );

  cmspixvsym100 = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmspixvsym100", 150, 0, 300 );
  cmspixvsym100->setTitle( "DUT pixel occupancy vs xmod (ym100);telescope x_{DUT} mod at ym100 [#mum];clusters" );

  cmspixvsym125 = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmspixvsym125", 150, 0, 300 );
  cmspixvsym125->setTitle( "DUT pixel occupancy vs xmod (ym125);telescope x_{DUT} mod at ym125 [#mum];clusters" );

  cmspixvsym150 = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmspixvsym150", 150, 0, 300 );
  cmspixvsym150->setTitle( "DUT pixel occupancy vs xmod (ym150);telescope x_{DUT} mod at ym150 [#mum];clusters" );

  cmspixvsym175 = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmspixvsym175", 150, 0, 300 );
  cmspixvsym175->setTitle( "DUT pixel occupancy vs xmod (ym175);telescope x_{DUT} mod at ym175 [#mum];clusters" );

  // KIT end

  cmsxyHitMap = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "cmsxyHitMap", 170, -8.5, 8.5, 170, -8.5, 8.5 );
  cmsxyHitMap->setTitle( "HitMap for chip;triplet x_{DUT} [mm];triplet y_{DUT} [mm];tracks" );

  cmsxyHitMapLowCharge = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "cmsxyHitMapLowCharge", 170, -8.5, 8.5, 170, -8.5, 8.5 );
  cmsxyHitMapLowCharge->setTitle( "HitMap for chip with low cluster charges;triplet x_{DUT} [mm];triplet y_{DUT} [mm];tracks" );  

  cmsqvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvsx", 380, -3.8, 3.8, 0, 100 );
  cmsqvsx->setTitle( "DUT q vs x;telescope x at DUT [mm];<q> [ke]" );

  cmsqvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvsy", 380, -3.8, 3.8, 0, 100 );
  cmsqvsy->setTitle( "DUT q vs y;telescope y at DUT [mm];<q> [ke]" );

  cmsqvsxm = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvsxm", 60, 0, 300, 0, 100 );
  cmsqvsxm->setTitle( "DUT q vs xmod;telescope x_{DUT} mod 300 [#mum];<q> [ke]" );

  cmsqvsym = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvsym", 40, 0, 200, 0, 100 );
  cmsqvsym->setTitle( "DUT q vs ymod;telescope y_{DUT} mod 200 [#mum];<q> [ke]" );

  cmsqvsxmym = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "cmsqvsxmym", 60, 0, 300, 40, 0, 200, 0, 250 );
  cmsqvsxmym->setTitle( "DUT cluster charge map;x_{track} mod 300 #mum;y_{track} mod 200 #mum;<cluster charge> [ke]" );

  // KIT: added for efficiency analysis - charge

  cmsqvsxm50 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvsxm50", 100, 0, 200, 0, 100 );
  cmsqvsxm50->setTitle( "DUT q vs ymod (xm50);telescope y_{DUT} mod at xm50 [#mum];<q> [ke]" );

  cmsqvsxm100 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvsxm100", 100, 0, 200, 0, 100 );
  cmsqvsxm100->setTitle( "DUT q vs ymod (xm100);telescope y_{DUT} mod at xm100 [#mum];<q> [ke]" );

  cmsqvsxm150 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvsxm150", 100, 0, 200, 0, 100 );
  cmsqvsxm150->setTitle( "DUT q vs ymod (xm150);telescope y_{DUT} mod at xm150 [#mum];<q> [ke]" );

  cmsqvsxm200 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvsxm200", 100, 0, 200, 0, 100 );
  cmsqvsxm200->setTitle( "DUT q vs ymod (xm200);telescope y_{DUT} mod at xm200 [#mum];<q> [ke]" );

  cmsqvsxm250 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvsxm250", 100, 0, 200, 0, 100 );
  cmsqvsxm250->setTitle( "DUT q vs ymod (xm250);telescope y_{DUT} mod at xm250 [#mum];<q> [ke]" );

  // ------

  cmsqvsym25 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvsym25", 150, 0, 300, 0, 100 );
  cmsqvsym25->setTitle( "DUT q vs xmod (ym25);telescope x_{DUT} mod at ym25 [#mum];<q> [ke]" );

  cmsqvsym50 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvsym50", 150, 0, 300, 0, 100 );
  cmsqvsym50->setTitle( "DUT q vs xmod (ym50);telescope x_{DUT} mod at ym50 [#mum];<q> [ke]" );

  cmsqvsym75 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvsym75", 150, 0, 300, 0, 100 );
  cmsqvsym75->setTitle( "DUT q vs xmod (ym75);telescope x_{DUT} mod at ym75 [#mum];<q> [ke]" );

  cmsqvsym100 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvsym100", 150, 0, 300, 0, 100 );
  cmsqvsym100->setTitle( "DUT q vs xmod (ym100);telescope x_{DUT} mod at ym100 [#mum];<q> [ke]" );

  cmsqvsym125 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvsym125", 150, 0, 300, 0, 100 );
  cmsqvsym125->setTitle( "DUT q vs xmod (ym125);telescope x_{DUT} mod at ym125 [#mum];<q> [ke]" );

  cmsqvsym150 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvsym150", 150, 0, 300, 0, 100 );
  cmsqvsym150->setTitle( "DUT q vs xmod (ym150);telescope x_{DUT} mod at ym150 [#mum];<q> [ke]" );

  cmsqvsym175 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvsym175", 150, 0, 300, 0, 100 );
  cmsqvsym175->setTitle( "DUT q vs xmod (ym175);telescope x_{DUT} mod at ym175 [#mum];<q> [ke]" );

  // KIT end

  cmsqvsddt = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvsddt", 40, -40, 40, 0, 250 );
  cmsqvsddt->setTitle( "DUT cluster charge vs phase;TLU-DUT #delta#Deltat [ns];<cluster charge> [ke]" );

  cmsqvst1 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvst1", 70, 0, 350, 0, 250 );
  cmsqvst1->setTitle( "DUT cluster charge vs time;run time [s];<cluster charge> [ke]" );

  cmsqvst2 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvst2", 140, 0, 700, 0, 250 );
  cmsqvst2->setTitle( "DUT cluster charge vs time;run time [s];<cluster charge> [ke]" );

  cmsqvst3 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvst3", 200, 0, 2000, 0, 250 );
  cmsqvst3->setTitle( "DUT cluster charge vs time;run time [s];<cluster charge> [ke]" );

  cmsqvst4 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvst4", 400, 0, 4000, 0, 250 );
  cmsqvst4->setTitle( "DUT cluster charge vs time;run time [s];<cluster charge> [ke]" );

  cmsrmsxvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsxvsx", 152, -3.8, 3.8, 0, 100 );
  cmsrmsxvsx->setTitle( "DUT x resolution vs x;telescope x [mm];RMS(#Deltax) [#mum]" );

  cmsrmsyvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsyvsx", 152, -3.8, 3.8, 0, 100 );
  cmsrmsyvsx->setTitle( "DUT y resolution vs x;telescope x [mm];RMS(#Deltay) [#mum]" );

  cmsrmsxvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsxvsy", 381, -3.81, 3.81, 0, 100 ); // shift by 10 um
  cmsrmsxvsy->setTitle( "DUT x resolution vs y;telescope y [mm];RMS(#Deltax) [#mum]" );

  cmsrmsyvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsyvsy", 381, -3.81, 3.81, 0, 100 ); // shift by 10 um
  cmsrmsyvsy->setTitle( "DUT y resolution vs y;telescope y [mm];RMS(#Deltay) [#mum]" );

  cmsrmsxvsxm = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsxvsxm", 60, 0, 300, 0, 100 );
  cmsrmsxvsxm->setTitle( "DUT x resolution;telescope x mod 300 [#mum];RMS(#Deltax) [#mum]" );

  cmsrmsyvsxm = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsyvsxm", 60, 0, 300, 0, 100 );
  cmsrmsyvsxm->setTitle( "DUT y resolution vs x mod 300;telescope x mod 300 [#mum];RMS(#Deltay) [#mum]" );

  cmsncolvsxm = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsncolvsxm", 60, 0, 300, 0, 10 );
  cmsncolvsxm->setTitle( "DUT cluster cols vs x mod 300;telescope x mod 300 [#mum];<cols/cluster>" );

  cmsnrowvsxm = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsnrowvsxm", 60, 0, 300, 0, 10 );
  cmsnrowvsxm->setTitle( "DUT cluster rows vs x mod 300;telescope x mod 300 [#mum];<rows/cluster>" );

  cmsrmsxvsym = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsxvsym", 40, 0, 200, 0, 100 );
  cmsrmsxvsym->setTitle( "DUT x resolution vs y mod 200;telescope y mod 200 [#mum];RMS(#Deltax) [#mum]" );

  cmsrmsyvsym = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsyvsym", 40, 0, 200, 0, 100 );
  cmsrmsyvsym->setTitle( "DUT y resolution vs y mod 200;telescope y mod 200 [#mum];RMS(#Deltay) [#mum]" );

  cmsrmsyvsym3 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsyvsym3", 60, 0, 300, 0, 100 );
  cmsrmsyvsym3->setTitle( "DUT y resolution vs y mod 300;telescope y mod 300 [#mum];RMS(#Deltay) [#mum]" );

  cmsrmsyvsym6 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsyvsym6", 120, 0, 600, 0, 100 );
  cmsrmsyvsym6->setTitle( "DUT y resolution vs y mod 600;telescope y mod 600 [#mum];RMS(#Deltay) [#mum]" );

  cmsrmsyvst = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsyvst", 60, 0, 3600, 0, 100 );
  cmsrmsyvst->setTitle( "DUT y resolution;time [s];RMS(#Deltay) [#mum]" );

  cmsrmsyvsddt = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsyvsddt", 50, -25, 25, 0, 100 );
  cmsrmsyvsddt->setTitle( "DUT y resolution vs ddt;TLU-DUT #delta#Deltat [ns];RMS(#Deltay) [#mum]" );

  cmsrmsxvsq = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsxvsq", 150, 0, 150, 0, 100 );
  cmsrmsxvsq->setTitle( "DUT x resolution vs q;cluster charge [ke];RMS(#Deltax) [#mum]" );

  cmsrmsyvsq = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsyvsq", 150, 0, 150, 0, 100 );
  cmsrmsyvsq->setTitle( "DUT y resolution vs q;cluster charge [ke];RMS(#Deltay) [#mum]" );

  cmsdyvseta = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdyvseta", 100, -1, 1, -200, 200 );
  cmsdyvseta->setTitle( "DUT dy vs eta;charge sharing eta;<#Deltay> [#mum]" );

  cmsrmsyvseta = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsyvseta", 100, -1, 1, 0, 100 );
  cmsrmsyvseta->setTitle( "DUT y resolution vs eta;charge sharing eta;RMS(#Deltay) [#mum]" );

  cmspMoyalvsq = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmspMoyalvsq", 100, 0, 100, 0, 1 );
  cmspMoyalvsq->setTitle( "cluster Moyal prob vs q;cluster charge [ke];<Moyal prob>" );

  cmspMoyalHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmspMoyal", 100, 0, 1 );
  cmspMoyalHisto->setTitle( "cluster Moyal prob;Moyal prob;DUT linked clusters" );

  cmsrmsyvsp = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsyvsp", 100, 0, 1, 0, 100 );
  cmsrmsyvsp->setTitle( "DUT y resolution vs prob;cluster charge [Moyal prob];RMS(#Deltay) [#mum]" );

  cmsnpxvsxmym = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "cmsnpxvsxmym", 60, 0, 300, 40, 0, 200, 0, 4.5 );
  cmsnpxvsxmym->setTitle( "DUT cluster size vs xm-ym;telescope x mod 300 [#mum];telescope y mod 200 [#mum];<pixels/cluster>" );

  cmsncolvsym = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsncolvsym", 40, 0, 200, 0, 4.5 );
  cmsncolvsym->setTitle( "DUT cols vs xmod;telescope x mod 300 [#mum];<cols/cluster>" );

  cmsnrowvsym = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsnrowvsym", 40, 0, 200, 0, 4.5 );
  cmsnrowvsym->setTitle( "DUT rows vs ymod;telescope y mod 200 [#mum];<rows/cluster>" );
  // KIT: added for efficiency analysis - cluster size

  cmsnpxvsxm50 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsnpxvsxm50", 100, 0, 200, 0, 100 );
  cmsnpxvsxm50->setTitle( "DUT cluster size vs ymod (xm50);telescope y_{DUT} mod at xm50 [#mum];pixels/cluster" );

  cmsnpxvsxm100 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsnpxvsxm100", 100, 0, 200, 0, 100 );
  cmsnpxvsxm100->setTitle( "DUT cluster size vs ymod (xm100);telescope y_{DUT} mod at xm100 [#mum];pixels/cluster" );

  cmsnpxvsxm150 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsnpxvsxm150", 100, 0, 200, 0, 100 );
  cmsnpxvsxm150->setTitle( "DUT cluster size vs ymod (xm150);telescope y_{DUT} mod at xm150 [#mum];pixels/cluster" );

  cmsnpxvsxm200 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsnpxvsxm200", 100, 0, 200, 0, 100 );
  cmsnpxvsxm200->setTitle( "DUT cluster size vs ymod (xm200);telescope y_{DUT} mod at xm200 [#mum];pixels/cluster" );

  cmsnpxvsxm250 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsnpxvsxm250", 100, 0, 200, 0, 100 );
  cmsnpxvsxm250->setTitle( "DUT cluster size vs ymod (xm250);telescope y_{DUT} mod at xm250 [#mum];pixels/cluster" );

  // ------

  cmsnpxvsym25 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsnpxvsym25", 150, 0, 300, 0, 100 );
  cmsnpxvsym25->setTitle( "DUT cluster size vs xmod (ym25);telescope x_{DUT} mod at ym25 [#mum];pixels/cluster" );

  cmsnpxvsym50 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsnpxvsym50", 150, 0, 300, 0, 100 );
  cmsnpxvsym50->setTitle( "DUT cluster size vs xmod (ym50);telescope x_{DUT} mod at ym50 [#mum];pixels/cluster" );

  cmsnpxvsym75 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsnpxvsym75", 150, 0, 300, 0, 100 );
  cmsnpxvsym75->setTitle( "DUT cluster size vs xmod (ym75);telescope x_{DUT} mod at ym75 [#mum];pixels/cluster" );

  cmsnpxvsym100 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsnpxvsym100", 150, 0, 300, 0, 100 );
  cmsnpxvsym100->setTitle( "DUT cluster size vs xmod (ym100);telescope x_{DUT} mod at ym100 [#mum];pixels/cluster" );

  cmsnpxvsym125 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsnpxvsym125", 150, 0, 300, 0, 100 );
  cmsnpxvsym125->setTitle( "DUT cluster size vs xmod (ym125);telescope x_{DUT} mod at ym125 [#mum];pixels/cluster" );

  cmsnpxvsym150 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsnpxvsym150", 150, 0, 300, 0, 100 );
  cmsnpxvsym150->setTitle( "DUT cluster size vs xmod (ym150);telescope x_{DUT} mod at ym150 [#mum];pixels/cluster" );

  cmsnpxvsym175 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsnpxvsym175", 150, 0, 300, 0, 100 );
  cmsnpxvsym175->setTitle( "DUT cluster size vs xmod (ym175);telescope x_{DUT} mod at ym175 [#mum];pixels/cluster" );

  // KIT end

  cmsetavsym = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsetavsym", 40, 0, 200, -1.5, 1.5 );
  cmsetavsym->setTitle( "DUT eta vs ymod;telescope y mod 200 [#mum];<eta>" );

  cmsetavsym2 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsetavsym2", 40, 0, 200, -1.5, 1.5 );
  cmsetavsym2->setTitle( "DUT eta 2-row vs ymod;telescope y mod 200 [#mum];<eta> 2-row" );

  cmsetavsym3 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsetavsym3", 60, 0, 300, -1.5, 1.5 );
  cmsetavsym3->setTitle( "DUT eta vs ymod;telescope y mod 200 [#mum];<eta>" );

  cmsym1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsym1", 40, 0, 200 );
  cmsym1Histo->setTitle( "ymod at DUT 1-row clus;telescope y mod 200 [#mum];1-row clusters" );

  cmsym2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsym2", 40, 0, 200 );
  cmsym2Histo->setTitle( "ymod at DUT 2-row clus;telescope y mod 200 [#mum];2-row clusters" );


  effxyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "effxy", 60, -4.5, 4.5, 90, -4.5, 4.5 );
  effxyHisto->setTitle( "six with REF link;triplet x_{DUT} [mm];triplet y_{DUT} [mm];tracks" );

  effvsxy = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "effvsxy", 60, -4.5, 4.5, 90, -4.5, 4.5, -1, 2 );
  effvsxy->setTitle( "DUT efficiency;telescope track x_{DUT} [mm];telescope track y_{DUT} [mm];CMS DUT efficiency" );

  effvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effvsx", 60, -4.5, 4.5, -1, 2 );
  effvsx->setTitle( "DUT efficiency;telescope track x_{DUT} [mm];CMS DUT efficiency" );

  effvsxg = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effvsxg", 60, -4.5, 4.5, -1, 2 );
  effvsxg->setTitle( "DUT efficiency;telescope track x_{DUT} [mm];CMS DUT efficiency" );

  effvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effvsy", 90, -4.5, 4.5, -1, 2 );
  effvsy->setTitle( "DUT efficiency;telescope track y_{DUT} [mm];<CMS DUT efficiency>" );

  eff300 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "eff300",  70, 0, 350, -1, 2 );
  eff300->setTitle( "DUT efficiency;time [s];<CMS DUT efficiency>" );

  eff600 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "eff600",  70, 0, 700, -1, 2 );
  eff600->setTitle( "DUT efficiency;time [s];<CMS DUT efficiency>" );

  effd600 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effd600",  70, 0, 700, -1, 2 );
  effd600->setTitle( "DUT efficiency DUT-REF ddt=0;time [s];<CMS DUT efficiency>" );

  effn600 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effn600",  70, 0, 700, -1, 2 );
  effn600->setTitle( "DUT efficiency ndri < 3;time [s];<CMS DUT efficiency>" );

  effm600 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effm600",  70, 0, 700, -1, 2 );
  effm600->setTitle( "DUT efficiency ndri = 1;time [s];<CMS DUT efficiency>" );

  eff1200 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "eff1200", 130, 0, 1300, -1, 2 );
  eff1200->setTitle( "DUT efficiency;time [s];<CMS DUT efficiency>" );

  eff1800 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "eff1800", 200, 0, 2000, -1, 2 );
  eff1800->setTitle( "DUT efficiency;time [s];<CMS DUT efficiency>" );

  eff3600 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "eff3600", 400, 0, 4000, -1, 2 );
  eff3600->setTitle( "DUT efficiency;time [s];<CMS DUT efficiency>" );

  effvsddt = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effvsddt", 40, -40, 40, -1, 2 );
  effvsddt->setTitle( "DUT efficiency;TLU-DUT #delta#Deltat [ns];<CMS DUT efficiency>" );

  effvsxmym = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "effvsxmym", 60, 0, 300, 40, 0, 200, -1, 2 );
  effvsxmym->setTitle( "DUT efficiency;telescope x' mod 300 [um];telescope y' mod 200 [um];efficiency" );

  // KIT: added for efficiency analysis - efficiency

  effvsxm50 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effvsxm50", 100, 0, 200, 0, 100 );
  effvsxm50->setTitle( "DUT efficiency vs ymod (xm50);telescope y_{DUT} mod at xm50 [#mum];efficiency" );

  effvsxm100 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effvsxm100", 100, 0, 200, 0, 100 );
  effvsxm100->setTitle( "DUT efficiency vs ymod (xm100);telescope y_{DUT} mod at xm100 [#mum];efficiency" );

  effvsxm150 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effvsxm150", 100, 0, 200, 0, 100 );
  effvsxm150->setTitle( "DUT efficiency vs ymod (xm150);telescope y_{DUT} mod at xm150 [#mum];efficiency" );

  effvsxm200 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effvsxm200", 100, 0, 200, 0, 100 );
  effvsxm200->setTitle( "DUT efficiency vs ymod (xm200);telescope y_{DUT} mod at xm200 [#mum];efficiency" );

  effvsxm250 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effvsxm250", 100, 0, 200, 0, 100 );
  effvsxm250->setTitle( "DUT efficiency vs ymod (xm250);telescope y_{DUT} mod at xm250 [#mum];efficiency" );

  // ------

  effvsym25 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effvsym25", 150, 0, 300, 0, 100 );
  effvsym25->setTitle( "DUT efficiency vs xmod (ym25);telescope x_{DUT} mod at ym25 [#mum];efficiency" );

  effvsym50 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effvsym50", 150, 0, 300, 0, 100 );
  effvsym50->setTitle( "DUT efficiency vs xmod (ym50);telescope x_{DUT} mod at ym50 [#mum];efficiency" );

  effvsym75 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effvsym75", 150, 0, 300, 0, 100 );
  effvsym75->setTitle( "DUT efficiency vs xmod (ym75);telescope x_{DUT} mod at ym75 [#mum];efficiency" );

  effvsym100 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effvsym100", 150, 0, 300, 0, 100 );
  effvsym100->setTitle( "DUT efficiency vs xmod (ym100);telescope x_{DUT} mod at ym100 [#mum];efficiency" );

  effvsym125 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effvsym125", 150, 0, 300, 0, 100 );
  effvsym125->setTitle( "DUT efficiency vs xmod (ym125);telescope x_{DUT} mod at ym125 [#mum];efficiency" );

  effvsym150 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effvsym150", 150, 0, 300, 0, 100 );
  effvsym150->setTitle( "DUT efficiency vs xmod (ym150);telescope x_{DUT} mod at ym150 [#mum];efficiency" );

  effvsym175 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effvsym175", 150, 0, 300, 0, 100 );
  effvsym175->setTitle( "DUT efficiency vs xmod (ym175);telescope x_{DUT} mod at ym175 [#mum];efficiency" );

  // KIT end

  rffvsxy = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "rffvsxy", 120, -12, 12, 60, -6, 6, -1, 2 );
  rffvsxy->setTitle( "REF efficiency;telescope x [mm];telescope y [mm];DUT REF efficiency" );

  rffvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "rffvsx", 100, -2, 8, -1, 2 );
  rffvsx->setTitle( "REF efficiency;telescope x [mm];DUT REF efficiency" );


  lkAvst = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "lkAvst", 360, 0, 3600, -0.5, 2.5 );
  lkAvst->setTitle( "DUT links/event;time [s];<tri-DUT links/event>" );

  ntriHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "ntri", 31, -0.5, 30.5 );
  ntriHisto->setTitle( "telescope triplets;0-1-2 triplets;events" );

  // triplet efficiency w.r.t. DUT:

  cmsxeHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsxe", 240, -12, 12 );
  cmsxeHisto->setTitle( "DUT hits;aligned x [mm];clusters" );

  cmsyeHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsye", 240, -6, 6 );
  cmsyeHisto->setTitle( "DUT hits;aligned y [mm];clusters" );

  cmsdxeHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdxe", 200, -500, 500 );
  cmsdxeHisto->setTitle( "Pixel + Triplet x;cluster + triplet #Sigmax [#mum];clusters" );

  cmsdyeHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdye", 500, -500, 500 );
  cmsdyeHisto->setTitle( "Pixel - triplet y;cluster - triplet #Deltay [#mum];clusters" );

  cmsnmHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsnm", 11, -0.5, 10.5 );
  cmsnmHisto->setTitle( "triplets linked to clusters;linked triplets;clusters" );

  trieffvsxy = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "trieffvsxy", 120, -12, 12, 60, -6, 6, -1, 2 );
  trieffvsxy->setTitle( "triplet efficiency;DUT cluster x [mm];DUT cluster y [mm];efficiency" );

  // driplets:
  AIDAProcessor::tree(this)->mkdir("Downstream");

  dridxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/dridx", 100, -100, 100 );
  dridxHisto->setTitle( "driplet dx;x_{4}-x_{m} [#mum];driplets" );

  dridyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/dridy", 100, -100, 100 );
  dridyHisto->setTitle( "driplet dy;y_{4}-y_{m} [#mum];driplets" );

  dridxvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridxvsx", 100, -10, 10, -100, 100 );
  dridxvsx->setTitle( "driplet x resid vs x;x [mm];driplet <#Deltax> [#mum]" );

  dridxvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridxvsy", 50, -5, 5, -100, 100 );
  dridxvsy->setTitle( "driplet x resid vs y;y [mm];driplet <#Deltax> [#mum]" );

  dridxvstx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridxvstx", 80, -2, 2, -100, 100 );
  dridxvstx->setTitle( "driplet x resid vs tx;t_{x} [mrad];driplet <#Deltax> [#mum]" );

  dridxvsty = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridxvsty", 80, -2, 2, -100, 100 );
  dridxvsty->setTitle( "driplet x resid vs ty;t_{y} [mrad];driplet <#Deltax> [#mum]" );

  dridyvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridyvsx", 100, -10, 10, -100, 100 );
  dridyvsx->setTitle( "driplet y resid vs x;x [mm];driplet <#Deltay> [#mum]" );

  dridyvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridyvsy", 50, -5, 5, -100, 100 );
  dridyvsy->setTitle( "driplet y resid vs y;y [mm];driplet <#Deltay> [#mum]" );

  dridyvstx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridyvstx", 80, -2, 2, -100, 100 );
  dridyvstx->setTitle( "driplet y resid vs tx;t_{x} [mrad];driplet <#Deltay> [#mum]" );

  dridyvsty = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridyvsty", 80, -2, 2, -100, 100 );
  dridyvsty->setTitle( "driplet y resid vs ty;t_{y} [mrad];driplet <#Deltay> [#mum]" );

  drixHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/drix", 240, -12, 12 );
  drixHisto->setTitle( "driplet x4;x4_{out} [mm];driplets" );

  driyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/driy", 120, -6, 6 );
  driyHisto->setTitle( "driplet y4;y4_{up} [mm];driplets" );

  drixyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "Downstream/drixy", 240, -12, 12, 120, -6, 6 );
  drixyHisto->setTitle( "driplet y4 vs x4;x4_{out} [mm];y4_{up} [mm];driplets" );

  dritxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/dritx", 100, -10, 10 );
  dritxHisto->setTitle( "driplet slope x;#theta_{x} [mrad];driplets" );

  drityHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/drity", 100, -10, 10 );
  drityHisto->setTitle( "driplet slope y;#theta_{y} [mrad];driplets" );


  drixrefHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "drixref", 240, -12, 12 );
  drixrefHisto->setTitle( "driplet at REF;driplet x_{out} at REF [mm];driplets" );

  driyrefHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "driyref", 120, -6, 6 );
  driyrefHisto->setTitle( "driplet at REF;driplet y_{up} at REF [mm];driplets" );

  drixyrefHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "drixyref", 240, -12, 12, 120, -6, 6 );
  drixyrefHisto->setTitle( "driplet at REF;driplet x_{out} at REF [mm];driplet y_{up} at REF [mm];driplets" );


  drixlkHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "drixlk", 240, -12, 12 );
  drixlkHisto->setTitle( "driplet with REF link;x_{REF} [mm];linked driplets" );

  driylkHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "driylk", 120, -6, 6 );
  driylkHisto->setTitle( "driplet with REF link;y_{REF} [mm];linked driplets" );

  drixylkHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "drixylk", 240, -12, 12, 120, -6, 6 );
  drixylkHisto->setTitle( "driplet with REF link;driplet x_{out} at REF [mm];driplet y_{up} at REF [mm];linked driplets" );

  refpixvsxmym = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "refpixvsxmym", 60, 0, 300, 40, 0, 200 );
  refpixvsxmym->setTitle( "REF pixel occupancy;x mod 300 #mum;y mod 200 #mum;clusters" );

  refqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refq", 100, 0, 100 );
  refqHisto->setTitle( "REF cluster charge linked;REF cluster charge [ke];REF linked clusters" );

  refqvsxmym = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "refqvsxmym", 60, 0, 300, 40, 0, 200, 0, 100 );
  refqvsxmym->setTitle( "REF cluster charge map;x mod 300 #mum;y mod 200 #mum;<cluster charge> [ke]" );


  refxxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "refxx", 52, -0.5, 51.5, 110, -11, 11 );
  refxxHisto->setTitle( "x correlation;REF pixel cluster col;driplet x [mm];REF clusters" );

  refyyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "refyy", 80, -0.5, 79.5, 55, -5.5, 5.5 );
  refyyHisto->setTitle( "y correlation;REF pixel cluster row;driplet y [mm];REF clusters" );

  refsxaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refsxa", 200, -10, 10 );
  refsxaHisto->setTitle( "REF + driplet x;cluster + driplet #Sigmax [mm];REF clusters" );

  refdxaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refdxa", 200, -10, 10 );
  refdxaHisto->setTitle( "REF - driplet x;cluster - driplet #Deltax [mm];REF clusters" );

  refsyaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refsya", 200, -10, 10 );
  refsyaHisto->setTitle( "REF + driplet y;cluster + driplet #Sigmay [mm];REF clusters" );

  refdyaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refdya", 200, -10, 10 );
  refdyaHisto->setTitle( "REF - driplet y;cluster - driplet #Deltay [mm];REF clusters" );

  refsxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refsx", 100, -5, 5 );
  refsxHisto->setTitle( "REF + driplet x;cluster + driplet #Sigmax [mm];REF clusters" );

  refdyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refdy", 100, -5, 5 );
  refdyHisto->setTitle( "REF - driplet y;cluster - driplet #Deltay [mm];REF clusters" );

  refsxcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refsxc", 100, -1, 1 );
  refsxcHisto->setTitle( "REF + driplet x;cluster + driplet #Sigmax [mm];REF clusters" );

  refdycHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refdyc", 100, -1, 1 );
  refdycHisto->setTitle( "REF - driplet y;cluster - driplet #Deltay [mm];REF clusters" );

  refdyvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "refdyvsx", 52, -26*0.15, 26*0.15, -150, 150 );
  refdyvsx->setTitle( "REF y resid vs x;x [mm];<REF cluster - driplet #Deltay> [#mum]" );

  refdyvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "refdyvsy", 80, -4, 4, -150, 150 );
  refdyvsy->setTitle( "REF y resid vs y;y [mm];<REF cluster - driplet #Deltay> [#mum]" );

  refdyvsty = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "refdyvsty", 80, -2, 2, -150, 150 );
  refdyvsty->setTitle( "REF y resid vs tet y;#theta_{y} [mrad];<REF cluster - driplet #Deltay> [#mum]" );

  reflkcolHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "reflkcol", 52, -0.5, 51.5 );
  reflkcolHisto->setTitle( "REF linked column;REF linked cluster col;REF linked clusters" );

  reflkrowHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "reflkrow", 80, -0.5, 79.5 );
  reflkrowHisto->setTitle( "REF linked row;REF linked cluster row;REF linked clusters" );


  bacsxaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacsxa", 100, -10, 10 );
  bacsxaHisto->setTitle( "DUT + driplet x;DUT cluster + driplet #Sigmax [mm];DUT clusters" );

  bacdyaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacdya", 100, -10, 10 );
  bacdyaHisto->setTitle( "DUT - driplet y;DUT cluster - driplet #Deltay [mm];DUT clusters" );

  bacsxcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacsxc", 200, -500, 500 );
  bacsxcHisto->setTitle( "DUT + driplet x;DUT cluster + driplet #Sigmax [mm];DUT clusters" );

  bacdycHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacdyc", 200, -500, 500 );
  bacdycHisto->setTitle( "DUT - driplet y;DUT cluster - driplet #Deltay [mm];DUT clusters" );

  bacsxcqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacsxcq", 200, -500, 500 );
  bacsxcqHisto->setTitle( "DUT + driplet x;DUT cluster + driplet #Sigmax [mm];DUT clusters" );

  bacdycqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacdycq", 100, -500, 500 );
  bacdycqHisto->setTitle( "DUT - driplet y;DUT cluster - driplet #Deltay [mm];DUT clusters" );


  ndriHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/ndri", 31, -0.5, 30.5 );
  ndriHisto->setTitle( "telescope driplets;3-4-5 driplets;events" );

  ndrirefHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "ndriref", 31, -0.5, 30.5 );
  ndrirefHisto->setTitle( "telescope driplets and DUT hit in event;3-4-5 driplets and DUT;events" );

  lkBvst = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "lkBvst", 360, 0, 3600, -0.5, 2.5 );
  lkBvst->setTitle( "REF links/event;time [s];<dri-REF links/event>" );

  //driplets-triplets
  // Tracks
  AIDAProcessor::tree(this)->mkdir("Tracks");

  nsixHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/nsix", 21, -0.5, 20.5 );
  nsixHisto->setTitle( "telescope six-plane-tracks;six-plane-tracks;events" );

  sixkxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkx", 100, -10, 10 );
  sixkxHisto->setTitle( "kink x;kink x [mrad];triplet-driplet pairs" );

  sixkyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixky", 100, -10, 10 );
  sixkyHisto->setTitle( "kink y;kink y [mrad];triplet-driplet pairs" );

  sixdxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixdx", 100, -1, 1 );
  sixdxHisto->setTitle( "six match x;match x [mm];triplet-driplet pairs" );

  sixdyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixdy", 100, -1, 1 );
  sixdyHisto->setTitle( "six match y;match y [mm];triplet-driplet pairs" );

  sixdxcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixdxc", 100, -250, 250 );
  sixdxcHisto->setTitle( "six match x;track #Deltax[#mum];triplet-driplet pairs" );

  sixdycHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixdyc", 100, -250, 250 );
  sixdycHisto->setTitle( "six match y;track #Deltay[#mum];triplet-driplet pairs" );

  sixkxcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkxc", 100, -10, 10 );
  sixkxcHisto->setTitle( "kink x, x-y matched;kink x [mrad];triplet-driplet pairs" );

  sixkycHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkyc", 100, -10, 10 );
  sixkycHisto->setTitle( "kink y, x-y matched;kink y [mrad];triplet-driplet pairs" );

  sixxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx", 240, -12, 12 );
  sixxHisto->setTitle( "six x at DUT;six x_{out} at DUT [mm];six-plane tracks" );

  sixyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy", 120, -6, 6 );
  sixyHisto->setTitle( "six y at DUT;six y_{up} at DUT [mm];six-plane tracks" );

  sixxyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "Tracks/sixxy", 240, -12, 12, 120, -6, 6 );
  sixxyHisto->setTitle( "six at DUT;six x_{out} at DUT [mm];six y_{up} at DUT [mm];six-plane tracks" );

  sixxycHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "Tracks/sixxyc", 240, -12, 12, 120, -6, 6 );
  sixxycHisto->setTitle( "six large kink;six x_{out} at DUT [mm];six y_{up} at DUT [mm];large kink tracks" );

  sixxylkHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "Tracks/sixxylk", 240, -12, 12, 120, -6, 6 );
  sixxylkHisto->setTitle( "six with REF link at DUT;six x_{out} at DUT [mm];six y_{up} at DUT [mm];six-plane tracks with REF link" );

  kinkvsxy = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "Tracks/kinkvsxy", 120, -12, 12, 60, -6, 6, 0, 100 );
  kinkvsxy->setTitle( "kink;six x_{out} at DUT [mm];six y_{up} at DUT [mm];<kink^{2}> [mrad^{2}]" );

  sixx0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx0", 240, -12, 12 );
  sixx0Histo->setTitle( "six x at 0;six x_{out} at 0 [mm];six-plane tracks" );

  sixy0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy0", 120, -6, 6 );
  sixy0Histo->setTitle( "six y at 0;six y_{up} at 0 [mm];six-plane tracks" );

  sixx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx1", 240, -12, 12 );
  sixx1Histo->setTitle( "six x at 1;six x_{out} at 1 [mm];six-plane tracks" );

  sixy1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy1", 120, -6, 6 );
  sixy1Histo->setTitle( "six y at 1;six y_{up} at 1 [mm];six-plane tracks" );

  sixx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx2", 240, -12, 12 );
  sixx2Histo->setTitle( "six x at 2;six x_{out} at 2 [mm];six-plane tracks" );

  sixy2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy2", 120, -6, 6 );
  sixy2Histo->setTitle( "six y at 2;six y_{up} at 2 [mm];six-plane tracks" );

  sixx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx3", 240, -12, 12 );
  sixx3Histo->setTitle( "six x at 3;six x_{out} at 3 [mm];six-plane tracks" );

  sixy3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy3", 120, -6, 6 );
  sixy3Histo->setTitle( "six y at 3;six y_{up} at 3 [mm];six-plane tracks" );

  sixx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx4", 240, -12, 12 );
  sixx4Histo->setTitle( "six x at 4;six x_{out} at 4 [mm];six-plane tracks" );

  sixy4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy4", 120, -6, 6 );
  sixy4Histo->setTitle( "six y at 4;six y_{up} at 4 [mm];six-plane tracks" );

  sixx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx5", 240, -12, 12 );
  sixx5Histo->setTitle( "six x at 5;six x_{out} at 5 [mm];six-plane tracks" );

  sixy5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy5", 120, -6, 6 );
  sixy5Histo->setTitle( "six y at 5;six y_{up} at 5 [mm];six-plane tracks" );

  // GBL:
  AIDAProcessor::tree(this)->mkdir("GBL");

  derxtiltHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/derxtilt", 100, -0.1, 0.1 );
  derxtiltHisto->setTitle( "ddx/dtilt;ddx/dtilt [mm/rad];align hits" );

  derytiltHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/derytilt", 100, -10, 10 );
  derytiltHisto->setTitle( "ddy/dtilt;ddy/dtilt [mm/rad];align hits" );

  derxturnHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/derxturn", 100, -10, 10 );
  derxturnHisto->setTitle( "ddx/dturn;ddx/dturn [mm/rad];align hits" );

  deryturnHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/deryturn", 100, -1, 1 );
  deryturnHisto->setTitle( "ddy/dturn;ddy/dturn [mm/rad];align hits" );

  selxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/selx", 240, -12, 12 );
  selxHisto->setTitle( "x at DUT, sel GBL;six x_{out} at DUT [mm];selected tracks" );

  selyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/sely", 120, -6, 6 );
  selyHisto->setTitle( "y at DUT, sel GBL;six y_{up} at DUT [mm];selected tracks" );

  selaxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/selax", 100, -5, 5 );
  selaxHisto->setTitle( "track angle x, sel GBL;x angle [mrad];tracks" );

  selayHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/selay", 100, -5, 5 );
  selayHisto->setTitle( "track angle y, sel GBL;y angle [mrad];tracks" );

  seldxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx", 100, -150, 150 );
  seldxHisto->setTitle( "track match x, sel GBL;#Deltax [#mum];tracks" );

  seldyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy", 100, -150, 150 );
  seldyHisto->setTitle( "track match y, sel GBL;#Deltay [#mum];tracks" );

  selkxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/selkx", 100, -10, 10 );
  selkxHisto->setTitle( "kink x, sel GBL;kink x [mrad];tracks" );

  selkyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/selky", 100, -10, 10 );
  selkyHisto->setTitle( "kink y, sel GBL;kink y [mrad];tracks" );

  seldx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx1", 100, -100, 100 );
  seldx1Histo->setTitle( "triplet resid x at 1, sel GBL;#Deltax [#mum];tracks" );

  seldy1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy1", 100, -100, 100 );
  seldy1Histo->setTitle( "triplet resid y at 1, sel GBL;#Deltay [#mum];tracks" );

  seldx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx3", 100, -1000, 1000 );
  seldx3Histo->setTitle( "triplet resid x at 3, sel GBL;#Deltax [#mum];tracks" );

  seldy3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy3", 100, -1000, 1000 );
  seldy3Histo->setTitle( "triplet resid y at 3, sel GBL;#Deltay [#mum];tracks" );

  seldx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx4", 100, -1500, 1500 );
  seldx4Histo->setTitle( "triplet resid x at 4, sel GBL;#Deltax [#mum];tracks" );

  seldy4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy4", 100, -1500, 1500 );
  seldy4Histo->setTitle( "triplet resid y at 4, sel GBL;#Deltay [#mum];tracks" );

  seldx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx5", 100, -3000, 3000 );
  seldx5Histo->setTitle( "triplet resid x at 5, sel GBL;#Deltax [#mum];tracks" );

  seldy5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy5", 100, -3000, 3000 );
  seldy5Histo->setTitle( "triplet resid y at 5, sel GBL;#Deltay [#mum];tracks" );

  seldx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx6", 100, -500, 500 );
  seldx6Histo->setTitle( "triplet resid x at DUT, sel GBL;#Deltax [#mum];tracks" );

  seldy6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy6", 100, -500, 500 );
  seldy6Histo->setTitle( "triplet resid y at DUT, sel GBL;#Deltay [#mum];tracks" );

  gblndfHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblndf", 16, -0.5, 15.5 );
  gblndfHisto->setTitle( "GBL fit NDF;GBL NDF;tracks" );

  gblchi2aHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblchi2a", 100, 0, 100 );
  gblchi2aHisto->setTitle( "GBL fit chi2, DoF 8;GBL chi2;tracks" );

  gblchi2bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblchi2b", 100, 0, 100 );
  gblchi2bHisto->setTitle( "GBL fit chi2;GBL chi2;tracks" );

  gblprbHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblprb", 100, 0, 1 );
  gblprbHisto->setTitle( "GBL fit probability;GBL fit probability;tracks" );

  // bad fits:

  badxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/badx", 240, -12, 12 );
  badxHisto->setTitle( "x at DUT, bad GBL;six x_{out} at DUT [mm];bad tracks" );

  badyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/bady", 120, -6, 6 );
  badyHisto->setTitle( "y at DUT, bad GBL;six y_{up} at DUT [mm];bad tracks" );

  badaxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/badax", 100, -5, 5 );
  badaxHisto->setTitle( "track angle x, bad GBL;x angle [mrad];tracks" );

  badayHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baday", 100, -5, 5 );
  badayHisto->setTitle( "track angle y, bad GBL;y angle [mrad];tracks" );

  baddxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx", 100, -150, 150 );
  baddxHisto->setTitle( "track match x, bad GBL;#Deltax [#mum];tracks" );

  baddyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy", 100, -150, 150 );
  baddyHisto->setTitle( "track match y, bad GBL;#Deltay [#mum];tracks" );

  badkxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/badkx", 100, -10, 10 );
  badkxHisto->setTitle( "kink x, bad GBL;kink x [mrad];tracks" );

  badkyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/badky", 100, -10, 10 );
  badkyHisto->setTitle( "kink y, bad GBL;kink y [mrad];tracks" );

  baddx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx1", 100, -100, 100 );
  baddx1Histo->setTitle( "triplet resid x at 1, bad GBL;#Deltax [#mum];tracks" );

  baddy1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy1", 100, -100, 100 );
  baddy1Histo->setTitle( "triplet resid y at 1, bad GBL;#Deltay [#mum];tracks" );

  baddx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx3", 100, -1000, 1000 );
  baddx3Histo->setTitle( "triplet resid x at 3, bad GBL;#Deltax [#mum];tracks" );

  baddy3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy3", 100, -1000, 1000 );
  baddy3Histo->setTitle( "triplet resid y at 3, bad GBL;#Deltay [#mum];tracks" );

  baddx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx4", 100, -1500, 1500 );
  baddx4Histo->setTitle( "triplet resid x at 4, bad GBL;#Deltax [#mum];tracks" );

  baddy4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy4", 100, -1500, 1500 );
  baddy4Histo->setTitle( "triplet resid y at 4, bad GBL;#Deltay [#mum];tracks" );

  baddx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx5", 100, -3000, 3000 );
  baddx5Histo->setTitle( "triplet resid x at 5, bad GBL;#Deltax [#mum];tracks" );

  baddy5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy5", 100, -3000, 3000 );
  baddy5Histo->setTitle( "triplet resid y at 5, bad GBL;#Deltay [#mum];tracks" );

  baddx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx6", 100, -250, 250 );
  baddx6Histo->setTitle( "triplet resid x at DUT, bad GBL;#Deltax [#mum];tracks" );

  baddy6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy6", 100, -250, 250 );
  baddy6Histo->setTitle( "triplet resid y at DUT, bad GBL;#Deltay [#mum];tracks" );

  // good fits:

  goodx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/goodx1", 100, -100, 100 );
  goodx1Histo->setTitle( "triplet resid x at 1, good GBL;#Deltax [#mum];tracks" );

  goody1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/goody1", 100, -100, 100 );
  goody1Histo->setTitle( "triplet resid y at 1, good GBL;#Deltay [#mum];tracks" );

  goodx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/goodx6", 100, -250, 250 );
  goodx6Histo->setTitle( "triplet resid x at 6, good GBL;#Deltax [#mum];tracks" );

  goody6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/goody6", 100, -250, 250 );
  goody6Histo->setTitle( "triplet resid y at 6, good GBL;#Deltay [#mum];tracks" );

  // look at fit:

  gblax0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax0", 100, -1, 1 );
  gblax0Histo->setTitle( "GBL x angle at plane 0;x angle at plane 0 [mrad];tracks" );

  gbldx0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx0", 100, -10, 10 );
  gbldx0Histo->setTitle( "GBL x shift at plane 0;x shift at plane 0 [#mum];tracks" );

  gblrx0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx0", 100, -25, 25 );
  gblrx0Histo->setTitle( "GBL x resid at plane 0;x resid at plane 0 [#mum];tracks" );

  gblpx0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx0", 100, -10, 10 );
  gblpx0Histo->setTitle( "GBL x pull at plane 0;x pull at plane 0 [#sigma];tracks" );

  gblqx0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx0", 100, -1, 1 );
  gblqx0Histo->setTitle( "GBL x kink at plane 0;x kink at plane 0 [mrad];tracks" );


  gblax1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax1", 100, -1, 1 );
  gblax1Histo->setTitle( "GBL x angle at plane 1;x angle at plane 1 [mrad];tracks" );

  gbldx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx1", 100, -100, 100 );
  gbldx1Histo->setTitle( "GBL x shift at plane 1;x shift at plane 1 [#mum];tracks" );

  gblrx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx1", 100, -25, 25 );
  gblrx1Histo->setTitle( "GBL x resid at plane 1;x resid at plane 1 [#mum];tracks" );

  gblpx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx1", 100, -10, 10 );
  gblpx1Histo->setTitle( "GBL x pull at plane 1;x pull at plane 1 [#sigma];tracks" );

  gblqx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx1", 100, -1, 1 );
  gblqx1Histo->setTitle( "GBL x kink at plane 1;x kink at plane 1 [mrad];tracks" );

  gblsx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblsx1", 100, -10, 10 );
  gblsx1Histo->setTitle( "GBL x kink at plane 1/error;x kink at plane 1/error;tracks" );

  gbltx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbltx1", 100, -10, 10 );
  gbltx1Histo->setTitle( "GBL x kink pull at plane 1;x kink pull at plane 1;tracks" );


  gblax2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax2", 100, -1, 1 );
  gblax2Histo->setTitle( "GBL x angle at plane 2;x angle at plane 2 [mrad];tracks" );

  gbldx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx2", 100, -20, 20 );
  gbldx2Histo->setTitle( "GBL x shift at plane 2;x shift at plane 2 [#mum];tracks" );

  gblrx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx2", 100, -25, 25 );
  gblrx2Histo->setTitle( "GBL x resid at plane 2;x resid at plane 2 [#mum];tracks" );

  gblpx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx2", 100, -10, 10 );
  gblpx2Histo->setTitle( "GBL x pull at plane 2;x pull at plane 2 [#sigma];tracks" );

  gblqx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx2", 100, -1, 1 );
  gblqx2Histo->setTitle( "GBL x kink at plane 2;x kink at plane 2 [mrad];tracks" );


  gblax3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax3", 100, -1, 1 );
  gblax3Histo->setTitle( "GBL x angle at plane 3;x angle at plane 3 [mrad];tracks" );

  gbldx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx3", 100, -250, 250 );
  gbldx3Histo->setTitle( "GBL x shift at plane 3;x shift at plane 3 [#mum];tracks" );

  gblrx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx3", 100, -25, 25 );
  gblrx3Histo->setTitle( "GBL x resid at plane 3;x resid at plane 3 [#mum];tracks" );

  gblpx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx3", 100, -10, 10 );
  gblpx3Histo->setTitle( "GBL x pull at plane 3;x pull at plane 3 [#sigma];tracks" );

  gblqx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx3", 100, -1, 1 );
  gblqx3Histo->setTitle( "GBL x kink at plane 3;x kink at plane 3 [mrad];tracks" );


  gblax4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax4", 100, -1, 1 );
  gblax4Histo->setTitle( "GBL x angle at plane 4;x angle at plane 4 [mrad];tracks" );

  gbldx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx4", 100, -500, 500 );
  gbldx4Histo->setTitle( "GBL x shift at plane 4;x shift at plane 4 [#mum];tracks" );

  gblrx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx4", 100, -25, 25 );
  gblrx4Histo->setTitle( "GBL x resid at plane 4;x resid at plane 4 [#mum];tracks" );

  gblpx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx4", 100, -10, 10 );
  gblpx4Histo->setTitle( "GBL x pull at plane 4;x pull at plane 4 [#sigma];tracks" );

  gblqx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx4", 100, -1, 1 );
  gblqx4Histo->setTitle( "GBL x kink at plane 4;x kink at plane 4 [mrad];tracks" );


  gblax5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax5", 100, -1, 1 );
  gblax5Histo->setTitle( "GBL x angle at plane 5;x angle at plane 5 [mrad];tracks" );

  gbldx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx5", 100, -1000, 1000 );
  gbldx5Histo->setTitle( "GBL x shift at plane 5;x shift at plane 5 [#mum];tracks" );

  gblrx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx5", 100, -25, 25 );
  gblrx5Histo->setTitle( "GBL x resid at plane 5;x resid at plane 5 [#mum];tracks" );

  gblpx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx5", 100, -10, 10 );
  gblpx5Histo->setTitle( "GBL x pull at plane 5;x pull at plane 5 [#sigma];tracks" );

  gblqx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx5", 100, -1, 1 );
  gblqx5Histo->setTitle( "GBL x kink at plane 5;x kink at plane 5 [mrad];tracks" );


  gblax6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax6", 100, -1, 1 );
  gblax6Histo->setTitle( "GBL x angle at DUT;x angle at DUT [mrad];tracks" );

  gbldx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx6", 100, -1000, 1000 );
  gbldx6Histo->setTitle( "GBL x shift at DUT;x shift at DUT [#mum];tracks" );

  gbldy6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldy6", 100, -1000, 1000 );
  gbldy6Histo->setTitle( "GBL y shift at DUT;y shift at DUT [#mum];tracks" );

  gblrx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx6", 100, -250, 250 );
  gblrx6Histo->setTitle( "GBL x resid at DUT;x resid at DUT [#mum];tracks" );

  gblry6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblry6", 100, -100, 100 );
  gblry6Histo->setTitle( "GBL y resid at DUT;y resid at DUT [#mum];tracks" );

  gblpx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx6", 100, -10, 10 );
  gblpx6Histo->setTitle( "GBL x pull at DUT;x pull at DUT [#sigma];tracks" );

  gblpy6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy6", 100, -10, 10 );
  gblpy6Histo->setTitle( "GBL y pull at DUT;y pull at DUT [#sigma];tracks" );

  gblqx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx6", 100, -10, 10 );
  gblqx6Histo->setTitle( "GBL x kink at DUT;x kink at DUT [mrad];tracks" );

  gblsx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblsx6", 100, -10, 10 );
  gblsx6Histo->setTitle( "GBL x kink at DUT/error;x kink at DUT/error;tracks" );

  gbltx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbltx6", 100, -10, 10 );
  gbltx6Histo->setTitle( "GBL x kink pull at DUT;x kink pull at DUT;tracks" );


  gblkx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx1", 100, -1, 1 );
  gblkx1Histo->setTitle( "GBL kink angle at plane 1;plane 1 kink [mrad];tracks" );

  gblkx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx2", 100, -1, 1 );
  gblkx2Histo->setTitle( "GBL kink angle at plane 2;plane 2 kink [mrad];tracks" );

  gblkx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx3", 100, -10, 10 );
  gblkx3Histo->setTitle( "GBL kink angle at plane 3;plane 3 kink [mrad];tracks" );

  gblkx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx4", 100, -1, 1 );
  gblkx4Histo->setTitle( "GBL kink angle at plane 4;plane 4 kink [mrad];tracks" );

  gblkx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx5", 100, -1, 1 );
  gblkx5Histo->setTitle( "GBL kink angle at plane 5;plane 5 kink [mrad];tracks" );

  gblkx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx6", 100, -1, 1 );
  gblkx6Histo->setTitle( "GBL kink angle at plane 6;plane 6 kink [mrad];tracks" );

  // z intersect:

  sixzx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzx3", 100, -50, 250 );
  sixzx3Histo->setTitle( "intersect z-x, kink > 3 mrad;intersect z(x) - z_{2} [mm];tracks" );

  sixzy3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzy3", 100, -50, 250 );
  sixzy3Histo->setTitle( "intersect z-y, kink > 3 mrad;intersect z(y) - z_{2} [mm];tracks" );

  sixzx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzx2", 100, -50, 250 );
  sixzx2Histo->setTitle( "intersect z-x, kink > 2 mrad;intersect z(x) - z_{2} [mm];tracks" );

  sixzy2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzy2", 100, -50, 250 );
  sixzy2Histo->setTitle( "intersect z-y, kink > 2 mrad;intersect z(y) - z_{2} [mm];tracks" );

  sixzx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzx1", 100, -50, 250 );
  sixzx1Histo->setTitle( "intersect z-x, kink > 1 mrad;intersect z(x) - z_{2} [mm];tracks" );

  sixzy1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzy1", 100, -50, 250 );
  sixzy1Histo->setTitle( "intersect z-y, kink > 1 mrad;intersect z(y) - z_{2} [mm];tracks" );

  sixkxzyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkxzy", 200, -20, 20 );
  sixkxzyHisto->setTitle( "kink x at DUT;kink x [mrad];tracks" );

  sixkyzxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkyzx", 200, -20, 20 );
  sixkyzxHisto->setTitle( "kink y at DUT;kink y [mrad];tracks" );

  sixkxzxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkxzx", 200, -20, 20 );
  sixkxzxHisto->setTitle( "kink x at DUT;kink x [mrad];tracks" );

  sixkyzyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkyzy", 200, -20, 20 );
  sixkyzyHisto->setTitle( "kink y at DUT;kink y [mrad];tracks" );

  // system times:
  tlutrigvstusHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/tlutrigvstus", 10000, 0, 10000 );
  tlutrigvstusHisto->setTitle( "TLU triggers;TLU time t [us];triggers" );

  cmsdtHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/cmsdt", 500, 0, 5000 );
  cmsdtHisto->setTitle( "DUT time between events;DUT time between events [us];events" );

  dutrefddtHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/dutrefddt", 21, -10.5, 10.5 );
  dutrefddtHisto->setTitle( "DUT - REF #delta#Deltat;DUT - REF #delta#Deltat[clocks];events" );

  sysrtHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/sysrt", 200, 0.9998, 1.0002 );
  sysrtHisto->setTitle( "TLU time / DUT time;event time ratio;events" );

  sysrdtHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/sysrdt", 200, 0.999, 1.001 );
  sysrdtHisto->setTitle( "TLU time / DUT time;time between events ratio;events" );

  dutddtnsHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/dutddtns", 100, -100, 100 );
  dutddtnsHisto->setTitle( "TLU - DUT time;TLU - DUT #delta#Deltat [ns];events" );

  refddtnsHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/refddtns", 100, -100, 100 );
  refddtnsHisto->setTitle( "TLU - REF time;TLU - REF #delta#Deltat [ns];events" );

  dutddtusHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/dutddtus", 100, -100, 100 );
  dutddtusHisto->setTitle( "TLU - DUT time;TLU - DUT #delta#Deltat [us];events" );

  dutddtmsHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/dutddtms", 200, -1, 1 );
  dutddtmsHisto->setTitle( "TLU - DUT time;TLU - DUT #delta#Deltat [ms];events" );

  dutddtvst = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/dutddtvst", 100, 0, 1000, -1E9, 1E9 );
  dutddtvst->setTitle( "TLU - DUT time;time [s];<#delta#Deltat> [ns]" );

  dutddtvsdt = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/dutddtvsdt", 100, 0, 10000, -1E9, 1E9 );
  dutddtvsdt->setTitle( "TLU - DUT time;time between events [us];<#delta#Deltat> [ns]" );

  ddtvst = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/ddtvst", 3000, 0, 300, -99, 99 );
  ddtvst->setTitle( "TLU - DUT time;time [s];<#delta#Deltat> [ns]" );

  ddtvstms = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/ddtvstms", 4000, 8, 12, -99, 99 ); // ms bins
  ddtvstms->setTitle( "TLU - DUT time;time [s];<#delta#Deltat> [ns]" );

  ddtvsdt = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/ddtvsdt", 200, 0, 2000, -99, 99 );
  ddtvsdt->setTitle( "TLU - DUT time;time between events [us];<#delta#Deltat> [ns]" );

  gapdtHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/gapdt", 100, 0, 1 );
  gapdtHisto->setTitle( "time between long clock gaps;time between long clock gaps [s];long clock gaps" );



  // Correlated events vs. time:
  correvt100Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "correvt100", 100, 0, 100000 );
  correvt100Histo->setTitle( "Correlated events (with matched DUT cluster);matched events/1000 events;events" );
  
  correvt300Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "correvt300", 300, 0, 300000 );
  correvt300Histo->setTitle( "Correlated events (with matched DUT cluster);matched events/1000 events;events" );
  
  correvt1000Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "correvt1000", 1000, 0, 1000000 );
  correvt1000Histo->setTitle( "Correlated events (with matched DUT cluster);matched events/1000 events;events" );
  
  correvt4000Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "correvt4000", 4000, 0, 4000000 );
  correvt4000Histo->setTitle( "Correlated events (with matched DUT cluster);matched events/1000 events;events" );

  // Skip Scanner
  if(_dutSkipScanner.get())
    _dutSkipScanner->bookHistos(this);

  streamlog_out( DEBUG2 ) << "Histogram booking completed \n\n" << std::endl;

#else

  streamlog_out( MESSAGE4 ) << "No histogram produced because Marlin doesn't use AIDA" << std::endl;

#endif

  return;
}

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

CMSSkipScanner::CMSSkipScanner(unsigned int range):
  _range(range), _eventNumber(0) {}

void CMSSkipScanner::bookHistos(marlin::Processor* processor)
{
  cmsxxHisto.resize(_range);
  cmsyyHisto.resize(_range);

  AIDAProcessor::tree(processor)->mkdir("SkipScan");

  for(unsigned int i = 0; i < _range; ++i)
    {
      std::stringstream skipNumS;
      skipNumS << i;
      std::string skipNum = skipNumS.str();

      // TODO: Find out how to put these into a subdirectory
      cmsxxHisto[i] = AIDAProcessor::histogramFactory(processor)->
	createHistogram2D( "SkipScan/cmsxx_skip" + skipNum, 52, -0.5, 51.5, 110, -11, 11 );
      cmsxxHisto[i]->setTitle( "x correlation (skip " + skipNum + ");CMS cluster col;telescope triplet x [mm];clusters" );

      cmsyyHisto[i] = AIDAProcessor::histogramFactory(processor)->
	createHistogram2D( "SkipScan/cmsyy_skip" + skipNum, 80, -0.5, 79.5, 55, -5.5, 5.5 );
      cmsyyHisto[i]->setTitle( "y correlation (skip " + skipNum + ");CMS cluster row;telescope triplet y [mm];clusters" );

      /*		cmssxHisto = AIDAProcessor::histogramFactory(fitter)->
			createHistogram1D( "skipScan/cmssx", 200, -500, 500 );
			cmssxHisto->setTitle( "Pixel + Telescope x;cluster + triplet #Sigmax [#mum];clusters" );

			cmsdyHisto = AIDAProcessor::histogramFactory(fitter)->
			createHistogram1D( "skipScan/cmsdy", 200, -500, 500 );
			cmsdyHisto->setTitle( "Pixel - telescope y;cluster - triplet #Deltay [#mum];clusters" );*/
    }
}

void CMSSkipScanner::newEvent(unsigned int eventNumber)
{
  _eventNumber = eventNumber;

  _triplets.push_front(std::vector<TelescopeTriplet>());
  while(_triplets.size() > _range) _triplets.pop_back();
  _cmsHits.clear();
}

void CMSSkipScanner::telescopeTriplet(double xs, double ys)
{
  TelescopeTriplet triplet;
  triplet.xs = xs;
  triplet.ys = ys;
  _triplets[0].push_back(triplet);

  // Correlate CMS hits with the new triplet
  for(unsigned int i = 0; i < _cmsHits.size(); ++i)
    {
      cmsxxHisto[0]->fill(_cmsHits[i].first, xs);
      cmsyyHisto[0]->fill(_cmsHits[i].second, ys);
    }
}

void CMSSkipScanner::cmsHit(double cmsCol, double cmsRow)
{
  // Some x-check...
  for(unsigned int i = 0; i < _cmsHits.size(); ++i)
    if(_cmsHits[i].first == cmsCol && _cmsHits[i].second == cmsRow)
      throw std::runtime_error("CMSSkipScanner: Double CMS hit!");
  std::pair<double, double> hit;
  hit.first = cmsCol;
  hit.second = cmsRow;
  _cmsHits.push_back(hit);

  // Correlate this hit with the triplets from the last _range telescope events
  for(unsigned i = 0; i < _triplets.size(); ++i)
    {
      for(unsigned int j = 0; j < _triplets[i].size(); ++j)
	{
	  cmsxxHisto[i]->fill(cmsCol, _triplets[i][j].xs);
	  cmsyyHisto[i]->fill(cmsRow, _triplets[i][j].ys);
	}
    }
}

int CMSSkipScanner::scan() const
{
  double maxXX = fabs(GetCorrelationFactor(cmsxxHisto[0])); unsigned int maxXSkip = 0;
  double maxYY = fabs(GetCorrelationFactor(cmsyyHisto[0])); unsigned int maxYSkip = 0;

  for(unsigned int i = 1; i < _triplets.size(); ++i)
    {
      double corrXX = fabs(GetCorrelationFactor(cmsxxHisto[i]));
      double corrYY = fabs(GetCorrelationFactor(cmsyyHisto[i]));

      if(corrXX > maxXX) { maxXX = corrXX; maxXSkip = i; }
      if(corrYY > maxYY) { maxYY = corrYY; maxYSkip = i; }
    }

  if(maxXX > 1.0e-1)
    std::cout << "Threshold reached in XX: " << maxXX << ", skip " << maxXSkip << std::endl;
  if(maxYY > 1.0e-1)
    std::cout << "Threshold reached in YY: " << maxYY << ", skip " << maxYSkip << std::endl;

  // No correlations found
  if(maxXX < 1.0e-1 && maxYY < 1.0e-1)
    return -1;

  if(maxXX < 1.0e-1 && maxYY > 1.0e-1)
    return maxYSkip;
  if(maxXX > 1.0e-1 && maxYY < 1.0e-1)
    return maxXSkip;

  if(maxXSkip != maxYSkip)
    {
      std::cout << "WARNING: Max correlations differ in X and Y!!" << std::endl;
      return -1;
    }

  return maxXSkip;
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
	if( abs(new_triplet.getdy()) > triplet_angle_cut * new_triplet.getdz()) continue;
	if( abs(new_triplet.getdx()) > triplet_angle_cut * new_triplet.getdz()) continue;

	// Setting cuts on the triplet residual on the middle plane
	if( abs(new_triplet.getdx(plane1)) > triplet_residual_cut) continue;
	if( abs(new_triplet.getdy(plane1)) > triplet_residual_cut) continue;

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
      
      
      if( abs(dy) < 0.4 ) sixdxcHisto->fill( dx*1E3 ); // sig = 17 um at 5 GeV
      if( abs(dx) < 0.4 ) sixdycHisto->fill( dy*1E3 );

      // match driplet and triplet:
      if( abs(dx) > intersect_residual_cut) continue;
      if( abs(dy) > intersect_residual_cut) continue;
	
      sixkxcHisto->fill( kx*1E3 );
      sixkycHisto->fill( ky*1E3 );
      sixxHisto->fill( -xA ); // -xA = x_DP = out
      sixyHisto->fill( -yA ); // -yA = y_DP = up
      sixxyHisto->fill( -xA, -yA ); // DP: x_out, y_up

      // Fill kink map histogram:
      if( abs( kx ) > 0.002 || abs( ky ) > 0.002 ) sixxycHisto->fill( -xA, -yA );

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
    dtfvstau->fill( tau, abs( frac - 0.5 ) ); // peak at 375.055 for run 5513
    frac = std::modf( (time_now_tlu - time_event0) / tau, &ganz ); // absolute phase
    tfvstau->fill( tau, abs( frac - 0.5 ) );
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


std::vector<EUTelAnalysisCMSPixel::cluster> EUTelAnalysisCMSPixel::GetClusters(std::vector<CMSPixel::pixel> * pixels) {

  // FIXME: BEWARE! this currently handles only single ROCs!

  int fCluCut = 1; // clustering: 1 = no gap (15.7.2012)

  std::vector<cluster> clusters;
  if(pixels->empty()) return clusters;

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
          for( unsigned int p = 0; p < c.vpix.size(); p++ ){//vpix in cluster so far
            int dr = c.vpix.at(p).row - pixels->at(i).row;
            int dc = c.vpix.at(p).col - pixels->at(i).col;
            if( (   dr>=-fCluCut) && (dr<=fCluCut) 
		&& (dc>=-fCluCut) && (dc<=fCluCut) ) {
              c.vpix.push_back(pixels->at(i));
	      gone[i] = 1;
              growing = 1;
              break;
            }
          }//loop over vpix
        }//not gone
      }//loop over all pix
    }while(growing);

    // Added all we could. determine position and append it to the list of clusters:
    for( std::vector<CMSPixel::pixel>::iterator p=c.vpix.begin();  p!=c.vpix.end();  p++){
      double Qpix = p->vcal; // calibrated [keV]
      if( Qpix < 0 ) Qpix = 1; // DP 1.7.2012
      c.charge += Qpix;
      c.sumA += p->raw; //Raw Chip pulse height value
      c.col += (*p).col*Qpix;
      c.row += (*p).row*Qpix;
      //c.xy[0] += (*p).xy[0]*Qpix;
      //c.xy[1] += (*p).xy[1]*Qpix;
    }

    c.size = c.vpix.size();

    streamlog_out(DEBUG4) << "(cluster with " << c.vpix.size() << " pixels)" << std::endl;

    if( !c.charge == 0 ) {
      c.col = c.col/c.charge;
      c.row = c.row/c.charge;
      //c.xy[0] /= c.charge;
      //c.xy[1] /= c.charge;
    }
    else {
      c.col = (*c.vpix.begin()).col;
      c.row = (*c.vpix.begin()).row;
      //c.xy[0] = (*c.vpix.begin()).xy[0];
      //c.xy[1] = (*c.vpix.begin()).xy[1];
      streamlog_out(DEBUG3) << "GetHits: cluster with zero charge" << std::endl;
    }

    clusters.push_back(c);//add cluster to vector

    //look for a new seed = unused pixel:
    while((++seed < pixels->size()) && (gone[seed]));
  }//while over seeds

  // nothing left to do, return clusters
  delete gone;
  return clusters;
}

//------------------------------------------------------------------------------
bool EUTelAnalysisCMSPixel::CalibratePixels(std::vector<CMSPixel::pixel> * pixels, EUTelAnalysisCMSPixel::calibration cal) {

  for( std::vector<CMSPixel::pixel>::iterator pix = pixels->begin(); pix != pixels->end(); pix++) {
    
    size_t col = (*pix).col;
    size_t row = (*pix).row;

    int ff = 1;
    if(cal.chip_id >= 110 && cal.chip_id <= 113)
      ff = 4; // AB -- sync to b2h with tbm=true... 5.10.2013

    double Ared = (*pix).raw / ff - cal.fitParameter[3][col][row]; // sub vert off, see gaintanh2ps.C
    double ma9 = cal.fitParameter[0][col][row];
    if( Ared >  ma9-1 ) ma9 = 1.000001 * Ared;
    if( Ared < -ma9+1 ) ma9 = Ared - 1;

    double keV = 0.45; // keV / large Vcal DAC
    if( cal.chip_id ==  10 ) keV = 0.36; // Landau peak at 24 ke, May 2012
    if( cal.chip_id ==  22 ) keV = 0.50; // Landau peak at 24 ke
    if( cal.chip_id >= 100 ) keV = 0.35; // xdb 15.7.2012 (7*0.05 = 0.35)
    if( cal.chip_id >= 110 && cal.chip_id <= 120) keV = 0.45; // AB -- sync to b2h, 5.10.2013
    //if( cal.chip_id == 202) keV = 0.307;
    if( cal.chip_id == 202) keV = 0.288; // VS -- data vs pixelav
					 // comparison, Landau peak
 					 // Feb 2014
    if( cal.chip_id == 203) keV = 0.324;
    
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

      ma9 = cal.fitParameter[3][col][row]; // p3 is negative

      (*pix).vcal = cal.fitParameter[1][col][row] * ( std::pow( -std::log( Ared / ma9 ), 1/cal.fitParameter[2][col][row] ) - cal.fitParameter[0][col][row] ) * keV;
      // q = ( (-ln((A-p4)/p3))^1/p2 - p0 )*p1

    }
    // DESY TanH Calibration:
    else {
      (*pix).vcal = (TMath::ATanH( Ared / ma9 ) 
		      * cal.fitParameter[1][col][row] + cal.fitParameter[2][col][row]) * keV; // [ke]
    }
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
  else if(strcmp(CalibrationType.c_str(),"desy_weibull_decorr") == 0) {
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

#endif
