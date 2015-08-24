// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// EUTelescope includes:
#include "EUTelAnalysisCMSPixel.h"
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

// ROOT includes ".h"
#include <TMath.h>
#include <TVectorD.h>
#include <TF1.h>
#include <TMatrixD.h>
#include <TVector3.h>
#include <TRotation.h>
#include "TH1D.h"

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;


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

  dutcolcorrHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "dutcolcorr", 52, -0.5, 51.5 );
  dutcolcorrHisto->setTitle( "DUT column Tsunami corrected;DUT cluster col;DUT clusters" );

  dutrowcorrHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "dutrowcorr", 80, -0.5, 79.5 );
  dutrowcorrHisto->setTitle( "DUT row Tsunami corrected;DUT cluster row;DUT clusters" );

  dutnpxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "dutnpx", 21, -0.5, 20.5 );
  dutnpxHisto->setTitle( "DUT cluster size;DUT pixel per cluster;DUT clusters" );

  dutadcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "dutadc", 150, 0, 150 );
  dutadcHisto->setTitle( "DUT cluster charge;DUT cluster charge [ke];DUT clusters" );

  dutpxcolHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "dutpxcol", 52, -0.5, 51.5 );
  dutpxcolHisto->setTitle( "DUT column;DUT pixel col;DUT pixels" );

  dutpxrowHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "dutpxrow", 80, -0.5, 79.5 );
  dutpxrowHisto->setTitle( "DUT row;DUT pixel row;DUT pixels" );

  dutpxadcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "dutpxadc", 120, -100, 500 );
  dutpxadcHisto->setTitle( "DUT pixel ADC;DUT pixel PH [ADC];DUT pixels" );

  dutpxqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "dutpxq", 100, 0, 25 );
  dutpxqHisto->setTitle( "DUT pixel ADC;DUT pixel PH [ke];DUT pixels" );


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

  triddaMindutHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "triddaMindut", 120, 0, 10 );
  triddaMindutHisto->setTitle( "minimal triplet distance at DUT;triplet distance at DUT [mm];telescope triplets" );

  // DUT pixel vs triplets:

  cmstimingcut = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmstimingcut", 140, 0, 700 );
  cmstimingcut->setTitle( "check if timing cut was applied;time[s]" );
  
  twoClusterDistanceHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "twoClusterDistance", 200, 0, 5 );
  twoClusterDistanceHisto->setTitle( "Two Cluster Distance;cluster distance [mm];clusters" );

  twoClusterXDistanceHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "twoClusterXDistance", 200, 0, 5 );
  twoClusterXDistanceHisto->setTitle( "Two Cluster X Distance;cluster distance [mm];clusters" );

  twoClusterYDistanceHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "twoClusterYDistance", 200, 0, 5 );
  twoClusterYDistanceHisto->setTitle( "Two Cluster Y Distance;cluster distance [mm];clusters" );

  twoClusterDistanceLostSeedHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "twoClusterDistanceLostSeed", 200, 0, 5 );
  twoClusterDistanceLostSeedHisto->setTitle( "Two Cluster Distance Lost Seed;cluster distance [mm];clusters" );

  twoClusterXDistanceLostSeedHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "twoClusterXDistanceLostSeed", 200, 0, 5 );
  twoClusterXDistanceLostSeedHisto->setTitle( "Two Cluster X Distance Lost Seed;cluster distance [mm];clusters" );

  twoClusterYDistanceLostSeedHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "twoClusterYDistanceLostSeed", 200, 0, 5 );
  twoClusterYDistanceLostSeedHisto->setTitle( "Two Cluster Y Distance Lost Seed;cluster distance [mm];clusters" );

  twoClusterDistanceLinkedTrackHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "twoClusterDistanceLinkedTrack", 200, 0, 5 );
  twoClusterDistanceLinkedTrackHisto->setTitle( "Two Cluster Distance for linked tracks;cluster distance [mm];clusters" );

  twoClusterXDistanceLinkedTrackHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "twoClusterXDistanceLinkedTrack", 200, 0, 5 );
  twoClusterXDistanceLinkedTrackHisto->setTitle( "Two Cluster X  Distance for linked tracks;cluster distance [mm];clusters" );

  twoClusterYDistanceLinkedTrackHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "twoClusterYDistanceLinkedTrack", 200, 0, 5 );
  twoClusterYDistanceLinkedTrackHisto->setTitle( "Two Cluster Y Distance for linked tracks;cluster distance [mm];clusters" );

  twoClusterDistanceLostSeedLinkedTrackHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "twoClusterDistanceLostSeedLinkedTrack", 200, 0, 5 );
  twoClusterDistanceLostSeedLinkedTrackHisto->setTitle( "Two Cluster Distance Lost Seed for linked tracks;cluster distance [mm];clusters" );

  twoClusterXDistanceLostSeedLinkedTrackHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "twoClusterXDistanceLostSeedLinkedTrack", 200, 0, 5 );
  twoClusterXDistanceLostSeedLinkedTrackHisto->setTitle( "Two Cluster X Distance Lost Seed for linked tracks;cluster distance [mm];clusters" );

  twoClusterYDistanceLostSeedLinkedTrackHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "twoClusterYDistanceLostSeedLinkedTrack", 200, 0, 5 );
  twoClusterYDistanceLostSeedLinkedTrackHisto->setTitle( "Two Cluster Y Distance Lost Seed for linked tracks;cluster distance [mm];clusters" );

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

  cmsskwHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsskw", 80, -0.2, 0.2 );
  cmsskwHisto->setTitle( "DUT cluster skew (tilt);DUT cluster skew;clusters" );

  cmsskw3pxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsskw3px", 80, -0.2, 0.2 );
  cmsskw3pxHisto->setTitle( "DUT cluster skew (tilt) 3px;DUT cluster skew;clusters" );

  cmsskw4pxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsskw4px", 80, -0.2, 0.2 );
  cmsskw4pxHisto->setTitle( "DUT cluster skew (tilt) 4px;DUT cluster skew;clusters" );

  cmsskwfcqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsskwfcq", 80, -0.2, 0.2 );
  cmsskwfcqHisto->setTitle( "DUT linked cluster skew (tilt);DUT linked cluster skew;clusters" );

  cmsskwfcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsskwfc", 80, -0.2, 0.2 );
  cmsskwfcHisto->setTitle( "DUT linked cluster skew (tilt);DUT linked cluster skew;clusters" );

  cmsskwcorrHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsskwcorr", 2000, -1000, +1000 );
  cmsskwcorrHisto->setTitle( "DUT cluster skew correction;DUT cluster skew correction [#mum];clusters" );

  cmsskw1colcogHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsskw1colcog", 52, -0.5, 51.5 );
  cmsskw1colcogHisto->setTitle( "DUT cluster column COG w/ large skew;DUT cluster col cog;clusters" );

  cmsskw1rowcogHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsskw1rowcog", 80, -0.5, 79.5 );
  cmsskw1rowcogHisto->setTitle( "DUT cluster row COG w/ large skew;DUT cluster row cog;clusters" );

  cmsskw1qHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsskw1q", 100, 0, 100 );
  cmsskw1qHisto->setTitle( "DUT cluster charge linked w/ large skew;DUT cluster charge [ke];DUT linked clusters" );

  cmsskw1ncolHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsskw1ncol", 21, -0.5, 20.5 );
  cmsskw1ncolHisto->setTitle( "DUT cluster col size /w large skew;columns per cluster;DUT clusters" );

  cmsskw1nrowHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsskw1nrow", 11, -0.5, 10.5 );
  cmsskw1nrowHisto->setTitle( "DUT cluster row size w/ large skew;rows per cluster;DUT clusters" );

  cmsskw0qHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsskw0q", 100, 0, 100 );
  cmsskw0qHisto->setTitle( "DUT cluster charge linked w/ small skew;DUT cluster charge [ke];DUT linked clusters" );

  cmsskw0ncolHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsskw0ncol", 21, -0.5, 20.5 );
  cmsskw0ncolHisto->setTitle( "DUT cluster col size /w small skew;columns per cluster;DUT clusters" );

  cmsskw0nrowHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsskw0nrow", 11, -0.5, 10.5 );
  cmsskw0nrowHisto->setTitle( "DUT cluster row size w/ small skew;rows per cluster;DUT clusters" );

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

  cmsdx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdx5", 440, -11, 11 );
  cmsdx5Histo->setTitle( "Pixel - Telescope x, skw corr;cluster + triplet #Sigmax [mm];clusters" );

  cmsdy5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdy5", 220, -5.5, 5.5 );
  cmsdy5Histo->setTitle( "Pixel - telescope y, skw corr;cluster - triplet #Deltay [mm];clusters" );

  cmsdxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdx", 200, -500, 500 );
  cmsdxHisto->setTitle( "Pixel + Telescope x;cluster + triplet #Sigmax [#mum];clusters" );

  cmsdyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdy", 500, -500, 500 );
  cmsdyHisto->setTitle( "Pixel - telescope y;cluster - triplet #Deltay [#mum];clusters" );

  cmsdy0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdy0", 500, -500, 500 );
  cmsdy0Histo->setTitle( "Pixel - telescope y, no skew corr;cluster - triplet #Deltay [#mum];clusters" );

  cmsdxfHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdxf", 200, -500, 500 );
  cmsdxfHisto->setTitle( "fiducial Pixel - telescope x;fiducial cluster - triplet #Deltax [#mum];fiducial clusters" );

  cmsdyfHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyf", 500, -500, 500 );
  cmsdyfHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfdc1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfdc1", 500, -500, 500 );
  cmsdyfdc1Histo->setTitle( "fiducial Pixel - telescope y, clusters contained in 1 DC;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfdc2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfdc2", 500, -500, 500 );
  cmsdyfdc2Histo->setTitle( "fiducial Pixel - telescope y, clusters spread over 2 DC;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

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

  cmsdyfc1cHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfc1c", 500, -500, 500 );
  cmsdyfc1cHisto->setTitle( "Pixel - telescope y;1-col cluster - triplet #Deltay [#mum];fiducial 1-col clusters" );

  cmsdyfc2cHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfc2c", 500, -500, 500 );
  cmsdyfc2cHisto->setTitle( "Pixel - telescope y;2-col cluster - triplet #Deltay [#mum];fiducial 2-col clusters" );

  cmsdyfc3cHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfc3c", 500, -500, 500 );
  cmsdyfc3cHisto->setTitle( "Pixel - telescope y;3-col cluster - triplet #Deltay [#mum];fiducial 3-col clusters" );

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
  cmsdyfctLowChargeHisto->setTitle( "fiducial Pixel - telescope y;fiducial low charge clusters - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfctHighChargeHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctHighCharge", 500, -500, 500 );
  cmsdyfctHighChargeHisto->setTitle( "fiducial Pixel - telescope y;fiducial low charge clusters - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfctOnePixelHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctOnePixel", 500, -500, 500 );
  cmsdyfctOnePixelHisto->setTitle( "fiducial Pixel - telescope y;fiducial one pixel clusters - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfctOnePixelLowChargeHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctOnePixelLowCharge", 500, -500, 500 );
  cmsdyfctOnePixelLowChargeHisto->setTitle( "fiducial Pixel - telescope y;fiducial one pixel, low charge cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfctOnePixelHighChargeHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctOnePixelHighCharge", 500, -500, 500 );
  cmsdyfctOnePixelHighChargeHisto->setTitle( "fiducial Pixel - telescope y;fiducial one pixel high charge cluster - triplet #Deltay [#mum];fiducial clusters" );

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

  cmsdy0fctq3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdy0fctq3", 500, -500, 500 );
  cmsdy0fctq3Histo->setTitle( "fiducial Pixel - telescope y, no skw corr;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfctqdotHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctqdot", 500, -500, 500 );
  cmsdyfctqdotHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfctq3dHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctq3d", 500, -500, 500 );
  cmsdyfctq3dHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );


  cmsdyfctq4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctq4", 500, -500, 500 );
  cmsdyfctq4Histo->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdy0fctq4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdy0fctq4", 500, -500, 500 );
  cmsdy0fctq4Histo->setTitle( "fiducial Pixel - telescope y, no skw corr;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfctq4dc1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctq4dc1", 500, -500, 500 );
  cmsdyfctq4dc1Histo->setTitle( "fiducial Pixel - telescope y, clusters contained in 1 DC;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfctq4dc2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfdc2", 500, -500, 500 );
  cmsdyfctq4dc2Histo->setTitle( "fiducial Pixel - telescope y, clusters spread over 2 DC;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfctq4dHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctq4d", 500, -500, 500 );
  cmsdyfctq4dHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdy0fctq4dHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdy0fctq4d", 500, -500, 500 );
  cmsdy0fctq4dHisto->setTitle( "fiducial Pixel - telescope y, no skw;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );


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

  cmslq0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmslq0", 100, 0, 100 );
  cmslq0Histo->setTitle( "DUT cluster charge linked within Q0 cuts;normal DUT cluster charge [ke];DUT linked clusters" );

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
    createHistogram1D( "cmsncol", 21, -0.5, 20.5 );
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
    createHistogram1D( "cmsqf", 400, 0, 400 );
  cmsqfHisto->setTitle( "DUT cluster charge linked fiducial;DUT cluster charge [ke];DUT linked fiducial clusters" );

  cmsqfOnePixeldyCutHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsqfOnePixeldyCut", 100, 0, 100 );
  cmsqfOnePixeldyCutHisto->setTitle( "DUT cluster charge linked fiducial;DUT cluster charge [ke];DUT linked fiducial clusters" );
  
  cmsqfNotOnePixeldyCutHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsqfNotOnePixeldyCut", 100, 0, 100 );
  cmsqfNotOnePixeldyCutHisto->setTitle( "DUT cluster charge linked fiducial;DUT cluster charge [ke];DUT linked fiducial clusters" );
  

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

  cmsq0fHistoRoot = new TH1D("cmsq0f_hist","DUT cluster charge linked fiducial;normal DUT cluster charge [ke];DUT linked fiducial clusters", 100, 0, 100 );

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

  cmsdy0vsxm = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdy0vsxm", 60, 0, 300, -150, 150 );
  cmsdy0vsxm->setTitle( "DUT y resid vs xmod, no skew corr;telescope x mod 300 [#mum];<DUT cluster - telescope triplet #Deltay> [#mum]" );

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

  cmspix1vsxmym = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "cmspix1vsxmym", 60, 0, 300, 40, 0, 200 );
  cmspix1vsxmym->setTitle( "DUT pixel occupancy for one pixel clusters;x mod 300 #mum;y mod 200 #mum;clusters" );

  cmspixvsxmymLowCharge = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "cmspixvsxmymLowCharge", 60, 0, 300, 40, 0, 200 );
  cmspixvsxmymLowCharge->setTitle( "DUT pixel occupancy for low charge clusters;x mod 300 #mum;y mod 200 #mum;clusters" );

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

  cmsqvsxmymdot = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "cmsqvsxmymdot", 60, 0, 300, 40, 0, 200, 0, 250 );
  cmsqvsxmymdot->setTitle( "DUT cluster charge map - Dot;x_{track} mod 300 #mum;y_{track} mod 200 #mum;<cluster charge> [ke]" );

  cmsskwvsym = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsskwvsym", 40, 0, 200, -0.2, 0.2 );
  cmsskwvsym->setTitle( "DUT skew vs ymod;telescope y_{DUT} mod 200 [#mum];cluster column skew" );

  cmsskwvsxm = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsskwvsxm", 60, 0, 300, -0.2, 0.2 );
  cmsskwvsxm->setTitle( "DUT skew vs xmod;telescope x_{DUT} mod 300 [#mum];cluster column skew" );

  cmsdyvsskw = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdyvsskw", 80, -0.2, 0.2, -200, 200);
  cmsdyvsskw->setTitle( "DUT skew vs ymod;cluster column skew;<#Deltay> [#mum]" );

  cmsdy0vsskw = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdy0vsskw", 80, -0.2, 0.2, -200, 200);
  cmsdy0vsskw->setTitle( "DUT skew vs ymod, no skw corr;cluster column skew;<#Deltay> [#mum]" );

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

  nTripClus = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "nTripClus", 11, -0.5, 10.5 );
  nTripClus->setTitle( "Clusters per triplet; number of clusters; events" );

  nTripClusLostSeed = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "nTripClusLostSeed", 11, -0.5, 10.5 );
  nTripClusLostSeed->setTitle( "Clusters per triplet for tracks with lost seed Cluster; number of clusters; events" );

  nTripPixels = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "nTripPixels", 11, -0.5, 10.5 );
  nTripPixels->setTitle( "Pixels per triplet; number of Pixels; events" );

  nTripPixelsLostSeed = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "nTripPixelsLostSeed",11, -0.5, 10.5 );
  nTripPixelsLostSeed->setTitle( "Pixels per triplet for tracks with lost seed Cluster; number of Pixels; events" );

  nLinkedTripClus = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "nLinkedTripClus", 11, -0.5, 10.5 );
  nLinkedTripClus->setTitle( "Linked Clusters per triplet; number of clusters; events" );

  nTripClusLinkedTrack = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "nTripClusLinkedTrack", 11, -0.5, 10.5 );
  nTripClusLinkedTrack->setTitle( "Clusters per triplet for linked Triplets; number of clusters; events" );

  nTripClusLostSeedLinkedTrack = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "nTripClusLostSeedLinkedTrack", 11, -0.5, 10.5 );
  nTripClusLostSeedLinkedTrack->setTitle( "Clusters per triplet for tracks with lost seed Cluster and linked Triplets; number of clusters; events" );

  nTripPixelsLinkedTrack = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "nTripPixelsLinkedTrack", 11, -0.5, 10.5 );
  nTripPixelsLinkedTrack->setTitle( "Pixels per triplet for linked Triplets; number of Pixels; events" );

  nTripPixelsLostSeedLinkedTrack = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "nTripPixelsLostSeedLinkedTrack",11, -0.5, 10.5 );
  nTripPixelsLostSeedLinkedTrack->setTitle( "Pixels per triplet for tracks with lost seed Cluster and linked Triplets; number of Pixels; events" );

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

  refddtnsHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/refddtns", 100, -100, 100 );
  refddtnsHisto->setTitle( "TLU - REF time;TLU - REF #delta#Deltat [ns];events" );

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

  streamlog_out( DEBUG2 ) << "Histogram booking completed \n\n" << std::endl;

#else

  streamlog_out( MESSAGE4 ) << "No histogram produced because Marlin doesn't use AIDA" << std::endl;

#endif

  return;
}

