// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTelAnalysisCMSPixel_h
#define EUTelAnalysisCMSPixel_h 1

#include <memory>
#include "CMSPixelDecoder.h"

#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include "lcio.h"
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>
#include <IMPL/TrackImpl.h>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IAxis.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/IProfile2D.h>
#include <AIDA/ITree.h>
#endif

// ROOT includes
#include <TMatrixD.h>
#include "TH1D.h"

// system includes <>
#include <string>
#include <vector>
#include <map>
#include <deque>


// Intrinsic nominal telescope resolution (MIMOSA26)
#define telescope_resx 3.5E-3
#define telescope_resy 3.5E-3

// Instrinsic nominal DUT resolution (CMS Pixel)
#define dut_resx 3.5E-3
#define dut_resy 3.5E-3

namespace eutelescope {

  class EUTelAnalysisCMSPixel : public marlin::Processor {

  public:

    //! Returns a new instance of EUTelAnalysisCMSPixel
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelAnalysisCMSPixel
     */
    virtual Processor*  newProcessor() { return new EUTelAnalysisCMSPixel; }
    EUTelAnalysisCMSPixel(const EUTelAnalysisCMSPixel&); 
    void operator=(EUTelAnalysisCMSPixel const&); 


    //! Default constructor
    EUTelAnalysisCMSPixel();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution.
     *
     */
    virtual void init();

    //! Called for every run.
    /*!
     * @param run the LCRunHeader of the current run
     */
    virtual void processRunHeader( LCRunHeader* run );

    //! Called every event
    /*! This is called for each event in the file.
     *
     *  @param evt the current LCEvent event
     */
    virtual void processEvent( LCEvent * evt );

    //! Check event method
    /*! This method is called by the Marlin execution framework as
     *  soon as the processEvent is over. It can be used to fill check
     *  plots. For the time being there is nothing to check and do in
     *  this slot.
     *
     *  @param evt The LCEvent event as passed by the ProcessMgr
     */
    virtual void check( LCEvent * evt );

    //! Book histograms
    /*! This method is used to books all required
     *  histograms. Histogram pointers are stored into
     *  _aidaHistoMap so that they can be recalled and filled
     * from anywhere in the code.
     */
    void bookHistos();

    //! Called after data processing for clean up.
    /*! Used to release memory allocated in init() step
     */
    virtual void end();

    //! CMS Pixel Cluster class:
    class cluster {
    public:
      cluster() : vpix(), seed_id(-1), size(0), sumA(0), charge(0), col(0), row(0), layer(0) {};
      std::vector <CMSPixel::pixel> vpix;
      int seed_id;
      int size;
      int sumA;
      float charge;
      float col,row;
      int layer;
    };

    //! Definition for calibration data format:
    class calibration {
    public:
      calibration() : fitParameter(), type(""), chip_id(0) {};
      double fitParameter[6][52][80];
      std::string type;
      int chip_id;
    };

    class hit {
    public:
      // Coordinates and their position uncertainty
      double x;
      double ex;
      double y;
      double ey;
      double z;
      double ez;
      // Plane to which the hit belongs
      unsigned int plane;
      // Overloading ostream operator for printing hits:
      friend std::ostream& operator << (std::ostream& out, const hit& point) // output
      {
	out << "(" << point.plane << ", " 
	    << point.x << " +/- " << point.ex << " | " 
	    << point.y << " +/- " << point.ey << " | " 
	    << point.z << " +/- " << point.ez << ")";
	return out;
      }
    };

    class triplet {
    public:
      triplet();
      triplet(hit hit0, hit hit1, hit hit2);

      // Keep track of linking status to DUT and REF:
      bool linked_dut;
      //bool linked_ref;

      hit getpoint_at(double z);


      // Returns x coordinate of the triplet at given z:
      double getx_at(double z);

      // Return dx = (x_end - x_start) for the full triplet:
      double getdx();

      // Returns dx = (x_measure - x_triplet) in the given plane ipl:
      double getdx(int ipl);

      // Returns dx for a given point:
      double getdx(hit point);


      // Returns y coordinate of the triplet at given z:
      double gety_at(double z);

      // Return dy = (y_end - y_start) for the full triplet:
      double getdy();

      // Returns dy = (y_measure - y_triplet) in the given plane ipl:
      double getdy(int ipl);

      // Returns dy for a given point:
      double getdy(hit point);


      // Return dz = (z_end - z_start) for the full triplet:
      double getdz();

      // Returning the hit for the given plane ID
      hit gethit(int plane);

      //! Returning the center point of the triplet:
      hit base();

      //! Returning the slope of the triplet (x,y):
      hit slope();

      friend std::ostream& operator << (std::ostream& out, triplet trip)
      {
	out << "Triplet: " << std::endl;
	for( std::map<unsigned int,hit>::iterator itr = trip.hits.begin(); itr != trip.hits.end(); itr++) {
	  out << "    " << itr->second << std::endl;
	}
	return out;
      };

      //FIXME stupid DUT data I don't know where to store else.
      // Contains lots of coord tranformation stuff and is likely
      // to be removed when we have proper coord trafo...
      double cmsdx;
      double cmsdy;
    private:
      void filltriplet(hit hit0, hit hit1, hit hit2) {
	hits.insert( std::pair<unsigned int,hit>(hit0.plane,hit0));
	hits.insert( std::pair<unsigned int,hit>(hit1.plane,hit1));
	hits.insert( std::pair<unsigned int,hit>(hit2.plane,hit2));
      };
      //! The hits belonging to the triplet:
      /* Use map since it's already ordered according to plane IDs.
       * We rely on begin() and rbegin() to deliver pointers to the first and last plane of the triplet.
       */
      std::map<unsigned int,hit> hits;   
    };

    class track {
    public:
      //! Default Track constructor. To be called with two triplets.
      track(triplet up, triplet down);

      //! Return the track kink angle in x
      double kink_x();

      //! Return the track kink angle in y
      double kink_y();

      //! Return the intersection point of the triplets
      hit intersect();

      //! Return the track upstream triplet
      triplet get_upstream();

      //! Return the track downstream triplet
      triplet get_downstream();

      //! Return the track hit in a given plane
      hit gethit(int plane);

    private:
      //! Members to store the up- and downstream triplets
      triplet upstream;
      triplet downstream;
    };

  private:

    //! Add padded zeros to given integer number to fill length
    std::string ZeroPadNumber(int num, int len);

    //! Helper function implementing the DESY chip ID scheme,
    //! returning the PSI46 chip type
    uint8_t GetChipTypeFromID(int chip_id);

    //! Timing plot being filled:
    void FillDeltaTPlots(double time_now_tlu, double time_before, double time_event0);

    //! Cluster statistics for DUT and REF detector:
    void FillClusterStatisticsPlots(std::vector<cluster> dutclusters, int dutpix, std::vector<cluster> refclusters, int refpix);

    // Calculate Point-To-Point Jacobian Transport Matrix for distance "ds"
    TMatrixD JacobianPointToPoint( double ds );

    //! Fill the telescope plane correlation plots:
    void TelescopeCorrelationPlots(std::vector<hit> * telescopehits);

    //! Find hit triplets from three telescope planes
    /*! This runs over all hits in the planes of the telescope and
     * tries to match triplets by comparing with the middle planes.
     * Two cut criteria can be set:
     * @param triplet_residual_cut Cut on the hit residual in the middle plane with respect to the triplet defined by first and last plane
     * @param triplet_angle_cut Cut on the triplet track angle
     *
     * @return a vector of found triplets among the given set of hits.
     */
    std::vector<triplet> * FindTriplets(std::vector<hit> * hits, unsigned int plane0, unsigned int plane1, unsigned int plane2);

    //! Match the upstream and downstream triplets to tracks
    std::vector<track> * MatchTriplets(std::vector<triplet> * up, std::vector<triplet> * down, double z_match);

    //! Do local clustering on the decoded pixel data:
    std::vector<cluster> GetClusters(std::vector<CMSPixel::pixel> * pixels, int chip);

    //! Calibrate the raw ROC pulse height into kev / VCal units:
    bool CalibratePixels(std::vector<CMSPixel::pixel> * pixels, calibration cal, double keV);

    //! Get VCAL->Electron conversion factor - either external or from
    //! chip list
    double GetConversionFactor(EUTelAnalysisCMSPixel::calibration cal, double conversion = 0);

    //! Fit returning the landau peak position, input is a cluster
    //! charge histogram such as cmsq0f:
    Double_t landau_gauss_peak(TH1* h);

    // Fit function: convoluted landau and gauss
    //Double_t fitLandauGauss( Double_t *x, Double_t *par );

    //! Initialize the calibration data, read in the files:
    bool InitializeCalibration(std::string gainfilename, int chip_id, std::string CalibrationType, calibration & cal);
  
  protected:
    //! Silicon planes parameters as described in GEAR
    /*! This structure actually contains the following:
     *  @li A reference to the telescope geoemtry and layout
     *  @li An integer number saying if the telescope is w/ or w/o DUT
     *  @li An integer number saying the number of planes in the
     *  telescope.
     *
     *  This object is provided by GEAR during the init() phase and
     *  stored here for local use.
     */
    gear::SiPlanesParameters * _siPlanesParameters;

    //! Silicon plane layer layout
    /*! This is the real geoemetry description. For each layer
     *  composing the telescope the relevant information are
     *  available.
     *
     *  This object is taken from the _siPlanesParameters during the
     *  init() phase and stored for local use
     */
    gear::SiPlanesLayerLayout * _siPlanesLayerLayout;

    std::string _inputCollectionTelescope;

    std::string _inputCollectionDUT;
    std::string _inputCollectionREF;

    std::string _inputPixelsDUT;
    std::string _inputPixelsREF;

    std::string _inputTrackCollection;

    //Analysis parameters
    bool _isFirstEvent;
    double _eBeam;
    int _nRun;
    int _nEvt;
    int _nTelPlanes;
    
    //! Timestamp of event number 0
    int64_t time_event0;

    //! Timestamp of event number 1
    int64_t time_event1;

    //! Reference timing, renewed whenever perfect match is found for correlation:
    int64_t time_reference;

    //! Timing: TLU frequency:
    double fTLU; // [Hz]
    double gTLU; // [GHz]

    //! DUT specific variables and parameters
    int _DUT_chip;
    std::string _DUT_gain;
    double _DUT_conversion;
    std::string _DUT_calibration_type;

    calibration dut_calibration;
    double _DUTalignx; //from cmsdxa
    double _DUTaligny; //from cmssya
    double _DUTz; // + _planePosition[2] // from -cmsdyvsty
    double _DUTrot; //compromise
    double _DUTtilt; // from cmsdyvsy 
    double _DUTturn; // from cmsdxvsx  flat

    //! Timing REF specific variables and parameters
    int _REF_chip;
    std::string _REF_gain;
    std::string _REF_calibration_type;
    double _REF_conversion;

    calibration ref_calibration;
    double _REFalignx; // [mm] refsxa
    double _REFaligny; // [mm] refdya
    double _REFz;
    double _REFrot; //from refdyvsx

    // Cuts for matching:
    double _cutx;
    double _cuty;
    double _slope_x;
    double _slope_y;

    std::string _skew_db;
    bool _have_skew_db;
    double skew_par0, skew_par1;
    std::string _CMS_gain_path;

    double _adc_correction;
    // Other variables only to print out the full runlist.csv line at
    // the end:
    std::string _gearfile;
    std::string _alignmentrun;

    // Partly outdated GEAR readings:
    int * _planeSort;
    int * _planeID;
    double * _planePosition;
    double * _planeThickness;
    double * _planeX0;
    double * _planeResolution;

    // Event buffer for event shifting in DUT:
    int _skip_dut, _skip_ref, _skip_tel;
    std::vector<std::vector<CMSPixel::pixel> > dut_event_buffer;
    std::vector<std::vector<CMSPixel::pixel> > ref_event_buffer;
    std::vector<std::vector<hit> > tel_event_buffer;

    // Reading in CMS clusters:
    CMSPixel::timing evt_time_dut;
    CMSPixel::timing evt_time_ref;
    std::vector <cluster> ClustDUT;
    std::vector <cluster> ClustREF;

    // Millepede binary file name:
    std::string m_millefilename;

    // definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

    AIDA::IHistogram1D *t1Histo, *t10Histo, *t100Histo, *t300Histo, *t600Histo, *t1000Histo, *t1800Histo, *t3600Histo;
    AIDA::IHistogram1D *dtHisto, *dtmsHisto, *logdtHisto, *logdtcmsHisto;
    AIDA::IProfile1D *dtfvstau, *tfvstau;
    AIDA::IHistogram1D * dtfHisto;
    AIDA::IProfile1D * dtfvst, *dtfvsdt;
    AIDA::IHistogram1D * tfHisto;
    AIDA::IProfile1D * tfvst, * tfvst1, * tfvst10, * tfvst100, * tfvst300;

    AIDA::IHistogram1D * nAllHitHisto;

    AIDA::IHistogram1D * dutnclusHisto, * dutcolHisto, * dutrowHisto, * dutnpxHisto, * dutadcHisto,* dutcolcorrHisto, * dutrowcorrHisto;
    AIDA::IHistogram1D * dutcolfirstHisto, * dutcollastHisto,* dutrowfirstHisto, * dutrowlastHisto, * dutcol1stevenqHisto, *dutcol1stoddqHisto;
    AIDA::IHistogram1D * cmscolfirstHisto, * cmscollastHisto, * cmscol1stevenqHisto, *cmscol1stoddqHisto;
    AIDA::IHistogram1D * dutpxcolHisto, * dutpxrowHisto, * dutpxadcHisto, * dutpxqHisto;
    AIDA::IHistogram1D * dutpx0adcHisto,* dutpx1adcHisto,* dutpx2adcHisto,* dutpx3adcHisto,* dutpx4adcHisto,* dutpx5adcHisto,* dutpx6adcHisto,* dutpx7adcHisto,* dutpx8adcHisto,* dutpx9adcHisto;
    AIDA::IHistogram1D * dutpx0cadcHisto,* dutpx1cadcHisto,* dutpx2cadcHisto,* dutpx3cadcHisto,* dutpx4cadcHisto,* dutpx5cadcHisto,* dutpx6cadcHisto,* dutpx7cadcHisto,* dutpx8cadcHisto,* dutpx9cadcHisto;
    AIDA::IHistogram1D * refnclusHisto, * refcolHisto, * refrowHisto, * refnpxHisto, * refadcHisto;

    AIDA::IHistogram1D * tlutrigvstusHisto;
    AIDA::IHistogram1D * cmsdtHisto;
    AIDA::IHistogram1D * dutrefddtHisto;
    AIDA::IHistogram1D * sysrtHisto;
    AIDA::IHistogram1D * sysrdtHisto;
    AIDA::IHistogram1D * dutddtnsHisto;
    AIDA::IHistogram1D * refddtnsHisto;
    AIDA::IHistogram1D * dutddtusHisto;
    AIDA::IHistogram1D * dutddtmsHisto;
    AIDA::IProfile1D * dutddtvst;
    AIDA::IProfile1D * dutddtvsdt;
    AIDA::IProfile1D * ddtvst;
    AIDA::IProfile1D * ddtvstms;
    AIDA::IProfile1D * ddtvsdt;
    AIDA::IHistogram1D * gapdtHisto;


    // Correlation plots for telescope planes
    AIDA::IHistogram1D * dx01Histo, * dy01Histo, * du01Histo, * dx02Histo, * dx03Histo, * dx04Histo, * dx05Histo, * dx12Histo, * dy12Histo, * du12Histo, * dx23Histo, * dy23Histo, * du23Histo, * dx34Histo, * dy34Histo, * du34Histo, * dx45Histo, * dy45Histo, * du45Histo;


    // triplets 0-1-2:

    AIDA::IHistogram1D * da02Histo;
    AIDA::IHistogram1D * db02Histo;

    AIDA::IProfile2D * dzcvsxy;
    AIDA::IProfile2D * z3vsxy;

    AIDA::IHistogram1D * tridxHisto;
    AIDA::IHistogram1D * tridyHisto;

    AIDA::IProfile1D * tridxvsx;
    AIDA::IProfile1D * tridxvsy;
    AIDA::IProfile1D * tridxvstx;
    AIDA::IProfile1D * tridxvsty;
    AIDA::IProfile1D * tridyvsx;
    AIDA::IProfile1D * tridyvsy;
    AIDA::IProfile1D * tridyvstx;
    AIDA::IProfile1D * tridyvsty;

    AIDA::IHistogram1D * tridx1Histo, * tridy1Histo, * tridx3Histo, * tridy3Histo, * tridx3bHisto, * tridy3bHisto, * tridx4Histo, * tridy4Histo, * tridx4bHisto, * tridy4bHisto, * tridx5Histo, * tridy5Histo, * tridx5bHisto, * tridy5bHisto, * trixHisto, * triyHisto, * tritxHisto, * trityHisto;
    AIDA::IHistogram2D * trixyHisto;

    AIDA::IHistogram1D * trixdutHisto; // at DUT
    AIDA::IHistogram1D * triydutHisto;
    AIDA::IHistogram2D * trixydutHisto;

    AIDA::IHistogram1D * triddaMindutHisto;

    AIDA::IHistogram2D * cmsxxHisto;
    AIDA::IHistogram2D * cmsyyHisto;

    AIDA::IHistogram1D * cmspxqHisto,* cmspxq1stHisto,* cmspxq2ndHisto,* cmspxq3rdHisto,* cmspxq2c1stHisto,* cmspxq2c1stcHisto,* cmspxq2c2ndHisto, * cmspxqhp1Histo,* cmspxqhp1cHisto,* cmspxqhp0Histo, * cmspxqcl2Histo, * cmspxqrow2Histo,* cmspxqeveHisto, * cmspxqoddHisto;
    AIDA::IHistogram1D * cms2pxq0Histo,* cms2pxq1Histo, * cms4pxq0Histo, * cms4pxq1Histo,* cms4pxq2Histo, * cms4pxq3Histo;
    AIDA::IProfile1D * cmspxqvsq;
    AIDA::IProfile1D * cmspxqvsqv;
    AIDA::IProfile1D * cmspxqvsxm;
    AIDA::IProfile1D * cmspxqvsym;
    AIDA::IHistogram1D *cmsskwHisto, *cmsskw3pxHisto, *cmsskw4pxHisto, *cmsskwfcqHisto, *cmsskwfctq4Histo, *cmsskwfctHisto, *cmsskwufctHisto, *cmsskwfcHisto, *cmsskwcorrHisto, *cmsskw1colcogHisto, *cmsskw1rowcogHisto, *cmsskw1qHisto, *cmsskw1ncolHisto, *cmsskw1nrowHisto, *cmsskw0qHisto, *cmsskw0ncolHisto, *cmsskw0nrowHisto, * cmssxaHisto, * cmsdyaHisto, * cmsdxaHisto, * cmssyaHisto, * cmsdx4Histo, * cmsdy4Histo, * cmsdx5Histo, *cmsdy5Histo, * cmsdxHisto, * cmsdyHisto, * cmsdyevenHisto, * cmsdyoddHisto, * cmsdy0Histo, * cmsdxfHisto, * cmsdyfHisto, * cmsdyfdc1Histo, * cmsdyfdc2Histo, * cmsdyfdc1qHisto, * cmsdyfdc2qHisto, * cmsdxfcHisto, * cmsdyfcHisto, * cmsdyfc1Histo, * cmsdyfc2Histo, * cmsdyfc3Histo, * cmsdyfc1cHisto, * cmsdyfc2cHisto, * cmsdyfc3cHisto, * cmsdyq0Histo, * cmsdyq1Histo, * cmsdyq2Histo, * cmsdyeta0Histo, * cmsdyeta1Histo, * cmsdxfctHisto, * cmsdyfctHisto, * cmsdyfcntHisto, * cmsdxfctqHisto, * cmsdyfctqHisto, * cmsdyfcntqHisto, * cmsdxfctq1Histo, * cmsdyfctq1Histo, * cmsdyfcntq1Histo, * cmsdyfctq1lHisto, * cmsdyfctq1rHisto, * cmsdxfctq2Histo, * cmsdyfctq2Histo, * cmsdxfctq3Histo, * cmsdyfctq3Histo, * cmsdy0fctq3Histo, * cmsdyfctqdotHisto, * cmsdyfctq3dHisto, * cmsdy0fctq4Histo, * cmsdyfctq4Histo, * cmsdyfctq4dc1Histo, * cmsdyfctq4dc2Histo, * cmsdyfctq4dc1qHisto, * cmsdyfctq4dc2qHisto, * cmsdyfctq4dHisto, * cmsdy0fctq4dHisto, * cmscolHisto, * cmsrowHisto, * cmsqHisto, * cmsq0Histo, * cmsqfctq4Histo, * cmslq0Histo, * trixlkHisto, * triylkHisto, * cmsqseedfHisto,* cmsdyfctOnePixelHisto, * cmsdyfctLowChargeHisto, * cmsdyfctHighChargeHisto, * cmsdyfctOnePixelLowChargeHisto,* cmsdyfctOnePixelHighChargeHisto,* cmsdxfctLowChargeHisto;

    AIDA::IHistogram2D * trixylkHisto;
    AIDA::IHistogram1D * cmstimingcut;
    AIDA::IHistogram1D * twoClusterDistanceHisto;
    AIDA::IHistogram1D * twoClusterXDistanceHisto;
    AIDA::IHistogram1D * twoClusterYDistanceHisto;
    AIDA::IHistogram1D * twoClusterDistanceLostSeedHisto;
    AIDA::IHistogram1D * twoClusterXDistanceLostSeedHisto;
    AIDA::IHistogram1D * twoClusterYDistanceLostSeedHisto;
    AIDA::IHistogram1D * twoClusterDistanceLinkedTrackHisto;
    AIDA::IHistogram1D * twoClusterXDistanceLinkedTrackHisto;
    AIDA::IHistogram1D * twoClusterYDistanceLinkedTrackHisto;
    AIDA::IHistogram1D * twoClusterDistanceLostSeedLinkedTrackHisto;
    AIDA::IHistogram1D * twoClusterXDistanceLostSeedLinkedTrackHisto;
    AIDA::IHistogram1D * twoClusterYDistanceLostSeedLinkedTrackHisto;
    AIDA::IProfile1D * cmsdxvsx;
    AIDA::IProfile1D * cmsdyvsx;
    AIDA::IProfile1D * cmsdxvsy;
    AIDA::IProfile1D * cmsdyvsy;
    AIDA::IProfile1D * cmsdxvstx;
    AIDA::IProfile1D * cmsdyvsty;
    AIDA::IHistogram2D * cmsdyvsxHisto;

    AIDA::IHistogram1D * cmsnpxHisto;
    AIDA::IHistogram1D * cmsnpxLowChargeHisto;
    AIDA::IHistogram1D * cmsnpx0Histo;
    AIDA::IHistogram1D * cmsnpx1Histo;
    AIDA::IHistogram1D * cmsnpx2Histo;
    AIDA::IHistogram1D * cmsnpx3Histo;
    AIDA::IHistogram1D * cmsncolHisto;
    AIDA::IHistogram1D * cmsncolLowChargeHisto;
    AIDA::IHistogram1D * cmsnrowHisto;
    AIDA::IHistogram1D * cmsnrowLowChargeHisto;
    AIDA::IHistogram1D * cmsnrowqHisto;
    AIDA::IProfile1D * cmsnrowvst1;
    AIDA::IHistogram1D * cmsetaHisto, *cmsetaevenHisto, * cmsetaoddHisto;
    AIDA::IHistogram1D * cmsqfHisto;
    AIDA::IHistogram1D * cmsqfOnePixeldyCutHisto;
    AIDA::IHistogram1D * cmsqfNotOnePixeldyCutHisto;
    AIDA::IHistogram1D * cmsqfcl1Histo;
    AIDA::IHistogram1D * cmsqfcl2Histo;
    AIDA::IHistogram1D * cmsqfrow1Histo;
    AIDA::IHistogram1D * cmsqfrow2Histo;
    AIDA::IHistogram1D * cmsq0fHisto;
    TH1D * cmsq0fHistoRoot;
    AIDA::IHistogram1D * cmsqf0Histo;
    AIDA::IHistogram1D * cmsqf1Histo;
    AIDA::IHistogram1D * cmsqf2Histo;
    AIDA::IHistogram1D * cmsqf3Histo;

    AIDA::IProfile1D * cmsdyvsxm, * cmsdy0vsxm, * cmsdxvsxm;
    AIDA::IProfile1D * cmsdyvsym, * cmsdxvsym;
    AIDA::IHistogram2D * cmspixvsxmym;
    AIDA::IHistogram2D * cmspix1vsxmym;
    AIDA::IHistogram2D * cmspixvsxmymLowCharge;

    // KIT: added for efficiency analysis - occupancy
    AIDA::IHistogram1D * cmspixvsxm50;
    AIDA::IHistogram1D * cmspixvsxm100; // Dot
    AIDA::IHistogram1D * cmspixvsxm150; // Pixel border
    AIDA::IHistogram1D * cmspixvsxm200; // Dot
    AIDA::IHistogram1D * cmspixvsxm250;
    AIDA::IHistogram1D * cmspixvsym25;
    AIDA::IHistogram1D * cmspixvsym50; // Dot
    AIDA::IHistogram1D * cmspixvsym75;
    AIDA::IHistogram1D * cmspixvsym100; // Pixel border
    AIDA::IHistogram1D * cmspixvsym125;
    AIDA::IHistogram1D * cmspixvsym150; // Dot
    AIDA::IHistogram1D * cmspixvsym175;

    AIDA::IHistogram2D * cmsxyHitMap;
    AIDA::IHistogram2D * cmsxyHitMapLowCharge;
    AIDA::IProfile1D * cmsqvsx;
    AIDA::IProfile1D * cmsqvsy;
    AIDA::IProfile1D * cmsqvsxm;
    AIDA::IProfile1D * cmsqvsym, * cmspxqsvsym;
    AIDA::IProfile2D * cmsqvsxmym, * cmspxqvsxmym;
    AIDA::IProfile2D * cmsqvsxmymdot;
    AIDA::IProfile1D * cmsskwvsym;
    AIDA::IProfile1D * cmsskwvsq, * cmsskwuvsq, * cmsskwvsqn;
    AIDA::IProfile1D * cmsqvsskw, * cmsqvsskwu;
    AIDA::IProfile1D * cmsskwvsxm;
    AIDA::IProfile1D * cmsdyvsskw, * cmsdyvsskwfctq4;
    AIDA::IProfile1D * cmsdy0vsskw, * cmsskwvsdy0, * cmsdy0vsskwfctq4, * cmsskwvsdy0fctq4;
    AIDA::IProfile1D * cmsqvsddt;
    AIDA::IProfile1D * cmsqvst1;
    AIDA::IProfile1D * cmsqvst2;
    AIDA::IProfile1D * cmsqvst3;
    AIDA::IProfile1D * cmsqvst4;
    // KIT: added for efficiency analysis - charge
    AIDA::IProfile1D * cmsqvsxm50;
    AIDA::IProfile1D * cmsqvsxm100; // Dot
    AIDA::IProfile1D * cmsqvsxm150; // Pixel border
    AIDA::IProfile1D * cmsqvsxm200; // Dot
    AIDA::IProfile1D * cmsqvsxm250;
    AIDA::IProfile1D * cmsqvsym25;
    AIDA::IProfile1D * cmsqvsym50; // Dot
    AIDA::IProfile1D * cmsqvsym75;
    AIDA::IProfile1D * cmsqvsym100; // Pixel border
    AIDA::IProfile1D * cmsqvsym125;
    AIDA::IProfile1D * cmsqvsym150; // Dot
    AIDA::IProfile1D * cmsqvsym175;

    AIDA::IProfile2D * cmsrmsxvsxmym, * cmsrmsyvsxmym, * cmsrmsxyvsxmym;
    AIDA::IProfile1D * cmsrmsxvsx;
    AIDA::IProfile1D * cmsrmsyvsx;
    AIDA::IProfile1D * cmsrmsxvsy;
    AIDA::IProfile1D * cmsrmsyvsy;
    AIDA::IProfile1D * cmsrmsxvsxm;
    AIDA::IProfile1D * cmsrmsyvsxm, * cmsrmsy0vsxm;
    AIDA::IProfile1D * cmsncolvsxm;
    AIDA::IProfile1D * cmsnrowvsxm;
    AIDA::IProfile1D * cmsrmsxvsym;
    AIDA::IProfile1D * cmsrmsyvsym;
    AIDA::IProfile1D * cmsrmsyvsym3;
    AIDA::IProfile1D * cmsrmsyvsym6;
    AIDA::IProfile1D * cmsrmsyvst;
    AIDA::IProfile1D * cmsrmsyvsddt;
    AIDA::IProfile1D * cmsrmsxvsq;
    AIDA::IProfile1D * cmsrmsyvsq;

    AIDA::IProfile1D * cmsdyvseta;
    AIDA::IProfile1D * cmsrmsyvseta;

    AIDA::IProfile1D * cmspMoyalvsq;
    AIDA::IHistogram1D * cmspMoyalHisto;
    AIDA::IProfile1D * cmsrmsyvsp;

    AIDA::IProfile2D * cmsnpxvsxmym;
    AIDA::IProfile2D * cmsnpx1vsxmym, * cmsnpx2vsxmym, * cmsnpx3vsxmym, * cmsnpx4vsxmym;
    AIDA::IProfile1D * cmsncolvsym;
    AIDA::IProfile1D * cmsnrowvsym;
    AIDA::IProfile1D * cmsetavsym;
    AIDA::IProfile1D * cmsetavsym3;
    AIDA::IProfile1D * cmsetavsym2;
    AIDA::IHistogram1D * cmsym1Histo;
    AIDA::IHistogram1D * cmsym2Histo;
    // KIT: added for efficiency analysis - cluster size
    AIDA::IProfile1D * cmsnpxvsxm50;
    AIDA::IProfile1D * cmsnpxvsxm100; // Dot
    AIDA::IProfile1D * cmsnpxvsxm150; // Pixel border
    AIDA::IProfile1D * cmsnpxvsxm200; // Dot
    AIDA::IProfile1D * cmsnpxvsxm250;
    AIDA::IProfile1D * cmsnpxvsym25;
    AIDA::IProfile1D * cmsnpxvsym50; // Dot
    AIDA::IProfile1D * cmsnpxvsym75;
    AIDA::IProfile1D * cmsnpxvsym100; // Pixel border
    AIDA::IProfile1D * cmsnpxvsym125;
    AIDA::IProfile1D * cmsnpxvsym150; // Dot
    AIDA::IProfile1D * cmsnpxvsym175;

    AIDA::IHistogram2D * effxyHisto;
    AIDA::IProfile2D * effvsxy;
    AIDA::IProfile1D * effvsx;
    AIDA::IProfile1D * effvsndri;
    AIDA::IProfile1D * effvsxg;
    AIDA::IProfile1D * effvsy;
    AIDA::IProfile1D * eff300;
    AIDA::IProfile1D * eff600;
    AIDA::IProfile1D * eff1200;
    AIDA::IProfile1D * eff1800;
    AIDA::IProfile1D * eff3600;
    AIDA::IProfile1D * effvsddt;
    AIDA::IProfile2D * effvsxmym;
    AIDA::IProfile1D * effd600;
    AIDA::IProfile1D * effn600;
    AIDA::IProfile1D * effm600;
    // KIT: added for efficiency analysis - efficiency
    AIDA::IProfile1D * effvsxm50;
    AIDA::IProfile1D * effvsxm100; // Dot
    AIDA::IProfile1D * effvsxm150; // Pixel border
    AIDA::IProfile1D * effvsxm200; // Dot
    AIDA::IProfile1D * effvsxm250;
    AIDA::IProfile1D * effvsym25;
    AIDA::IProfile1D * effvsym50; // Dot
    AIDA::IProfile1D * effvsym75;
    AIDA::IProfile1D * effvsym100; // Pixel border
    AIDA::IProfile1D * effvsym125;
    AIDA::IProfile1D * effvsym150; // Dot
    AIDA::IProfile1D * effvsym175;

    AIDA::IProfile2D * rffvsxy;
    AIDA::IProfile1D * rffvsx;

    AIDA::IHistogram1D * nTripClus;
    AIDA::IHistogram1D * nTripClusLostSeed;
    AIDA::IHistogram1D * nTripPixels;
    AIDA::IHistogram1D * nTripPixelsLostSeed;

    AIDA::IHistogram1D * nLinkedTripClus;

    AIDA::IHistogram1D * nTripClusLinkedTrack;
    AIDA::IHistogram1D * nTripClusLostSeedLinkedTrack;
    AIDA::IHistogram1D * nTripPixelsLinkedTrack;
    AIDA::IHistogram1D * nTripPixelsLostSeedLinkedTrack;

    AIDA::IHistogram1D * ntriHisto;
    AIDA::IProfile1D * lkAvst;

    // triplet eff w.r.t. CMS:

    AIDA::IHistogram1D * cmsxeHisto;
    AIDA::IHistogram1D * cmsyeHisto;

    AIDA::IHistogram1D * cmsdxeHisto;
    AIDA::IHistogram1D * cmsdyeHisto;

    AIDA::IHistogram1D * cmsnmHisto;
    AIDA::IProfile2D * trieffvsxy;

    // driplets 3-4-5:

    AIDA::IHistogram1D * dx35Histo;
    AIDA::IHistogram1D * dy35Histo;

    AIDA::IHistogram1D * dridxHisto;
    AIDA::IHistogram1D * dridyHisto;
    AIDA::IHistogram1D * drixHisto;
    AIDA::IHistogram1D * driyHisto;
    AIDA::IHistogram2D * drixyHisto;
    AIDA::IHistogram1D * dritxHisto;
    AIDA::IHistogram1D * drityHisto;

    AIDA::IProfile1D * dridxvsx;
    AIDA::IProfile1D * dridxvsy;
    AIDA::IProfile1D * dridxvstx;
    AIDA::IProfile1D * dridxvsty;
    AIDA::IProfile1D * dridyvsx;
    AIDA::IProfile1D * dridyvsy;
    AIDA::IProfile1D * dridyvstx;
    AIDA::IProfile1D * dridyvsty;

    AIDA::IHistogram1D * drixrefHisto; // at REF
    AIDA::IHistogram1D * driyrefHisto;
    AIDA::IHistogram2D * drixyrefHisto;

    AIDA::IHistogram1D * drixlkHisto;
    AIDA::IHistogram1D * driylkHisto;
    AIDA::IHistogram2D * drixylkHisto;
    AIDA::IHistogram2D * refpixvsxmym;

    AIDA::IHistogram1D * refqHisto;
    AIDA::IProfile2D * refqvsxmym;

    AIDA::IHistogram2D * refxxHisto; //REF vs driplet
    AIDA::IHistogram2D * refyyHisto;

    AIDA::IHistogram1D * refsxaHisto;
    AIDA::IHistogram1D * refdxaHisto;
    AIDA::IHistogram1D * refsyaHisto;
    AIDA::IHistogram1D * refdyaHisto;
    AIDA::IHistogram1D * refsxHisto;
    AIDA::IHistogram1D * refdyHisto;
    AIDA::IHistogram1D * refsxcHisto;
    AIDA::IHistogram1D * refdycHisto;

    AIDA::IProfile1D * refdyvsx;
    AIDA::IProfile1D * refdyvsy;
    AIDA::IProfile1D * refdyvsty;

    AIDA::IHistogram1D * reflkcolHisto;
    AIDA::IHistogram1D * reflkrowHisto;

    AIDA::IHistogram1D * bacsxaHisto;
    AIDA::IHistogram1D * bacdyaHisto;
    AIDA::IHistogram1D * bacsxcHisto;
    AIDA::IHistogram1D * bacdycHisto;
    AIDA::IHistogram1D * bacsxcqHisto;
    AIDA::IHistogram1D * bacdycqHisto;

    AIDA::IHistogram1D * ndriHisto;
    AIDA::IHistogram1D * ndrirefHisto;
    AIDA::IProfile1D * lkBvst;

    AIDA::IHistogram1D * nsixHisto;
    AIDA::IHistogram1D * sixkxHisto; //driplet-triplet
    AIDA::IHistogram1D * sixkyHisto;
    AIDA::IHistogram1D * sixdxHisto;
    AIDA::IHistogram1D * sixdyHisto;
    AIDA::IHistogram1D * sixdxcHisto;
    AIDA::IHistogram1D * sixdycHisto;

    AIDA::IHistogram1D * sixkxcHisto;
    AIDA::IHistogram1D * sixkycHisto;
    AIDA::IHistogram1D * sixxHisto;
    AIDA::IHistogram1D * sixyHisto;
    AIDA::IHistogram2D * sixxyHisto;
    AIDA::IHistogram2D * sixxycHisto;
    AIDA::IProfile2D * kinkvsxy;

    AIDA::IHistogram1D * sixx0Histo;
    AIDA::IHistogram1D * sixy0Histo;
    AIDA::IHistogram1D * sixx1Histo;
    AIDA::IHistogram1D * sixy1Histo;
    AIDA::IHistogram1D * sixx2Histo;
    AIDA::IHistogram1D * sixy2Histo;
    AIDA::IHistogram1D * sixx3Histo;
    AIDA::IHistogram1D * sixy3Histo;
    AIDA::IHistogram1D * sixx4Histo;
    AIDA::IHistogram1D * sixy4Histo;
    AIDA::IHistogram1D * sixx5Histo;
    AIDA::IHistogram1D * sixy5Histo;

    AIDA::IHistogram2D * sixxylkHisto;

    AIDA::IHistogram1D * derxtiltHisto;
    AIDA::IHistogram1D * derytiltHisto;
    AIDA::IHistogram1D * derxturnHisto;
    AIDA::IHistogram1D * deryturnHisto;

    AIDA::IHistogram1D * selxHisto;
    AIDA::IHistogram1D * selyHisto;
    AIDA::IHistogram1D * selaxHisto;
    AIDA::IHistogram1D * selayHisto;
    AIDA::IHistogram1D * seldxHisto;
    AIDA::IHistogram1D * seldyHisto;
    AIDA::IHistogram1D * selkxHisto;
    AIDA::IHistogram1D * selkyHisto;

    AIDA::IHistogram1D * seldx1Histo;
    AIDA::IHistogram1D * seldy1Histo;
    AIDA::IHistogram1D * seldx3Histo;
    AIDA::IHistogram1D * seldy3Histo;
    AIDA::IHistogram1D * seldx4Histo;
    AIDA::IHistogram1D * seldy4Histo;
    AIDA::IHistogram1D * seldx5Histo;
    AIDA::IHistogram1D * seldy5Histo;
    AIDA::IHistogram1D * seldx6Histo;
    AIDA::IHistogram1D * seldy6Histo;

    AIDA::IHistogram1D * gblndfHisto;
    AIDA::IHistogram1D * gblchi2aHisto;
    AIDA::IHistogram1D * gblchi2bHisto;
    AIDA::IHistogram1D * gblprbHisto;

    AIDA::IHistogram1D * badxHisto;
    AIDA::IHistogram1D * badyHisto;
    AIDA::IHistogram1D * badaxHisto;
    AIDA::IHistogram1D * badayHisto;
    AIDA::IHistogram1D * baddxHisto;
    AIDA::IHistogram1D * baddyHisto;
    AIDA::IHistogram1D * badkxHisto;
    AIDA::IHistogram1D * badkyHisto;

    AIDA::IHistogram1D * baddx1Histo;
    AIDA::IHistogram1D * baddy1Histo;
    AIDA::IHistogram1D * baddx3Histo;
    AIDA::IHistogram1D * baddy3Histo;
    AIDA::IHistogram1D * baddx4Histo;
    AIDA::IHistogram1D * baddy4Histo;
    AIDA::IHistogram1D * baddx5Histo;
    AIDA::IHistogram1D * baddy5Histo;
    AIDA::IHistogram1D * baddx6Histo;
    AIDA::IHistogram1D * baddy6Histo;

    AIDA::IHistogram1D * goodx1Histo;
    AIDA::IHistogram1D * goody1Histo;
    AIDA::IHistogram1D * goodx6Histo;
    AIDA::IHistogram1D * goody6Histo;

    AIDA::IHistogram1D * gblax0Histo;
    AIDA::IHistogram1D * gbldx0Histo;
    AIDA::IHistogram1D * gblrx0Histo;
    AIDA::IHistogram1D * gblpx0Histo;
    AIDA::IHistogram1D * gblqx0Histo;

    AIDA::IHistogram1D * gblax1Histo;
    AIDA::IHistogram1D * gbldx1Histo;
    AIDA::IHistogram1D * gblrx1Histo;
    AIDA::IHistogram1D * gblpx1Histo;
    AIDA::IHistogram1D * gblqx1Histo;
    AIDA::IHistogram1D * gblsx1Histo;
    AIDA::IHistogram1D * gbltx1Histo;

    AIDA::IHistogram1D * gblax2Histo;
    AIDA::IHistogram1D * gbldx2Histo;
    AIDA::IHistogram1D * gblrx2Histo;
    AIDA::IHistogram1D * gblpx2Histo;
    AIDA::IHistogram1D * gblqx2Histo;

    AIDA::IHistogram1D * gblax3Histo;
    AIDA::IHistogram1D * gbldx3Histo;
    AIDA::IHistogram1D * gblrx3Histo;
    AIDA::IHistogram1D * gblpx3Histo;
    AIDA::IHistogram1D * gblqx3Histo;

    AIDA::IHistogram1D * gblax4Histo;
    AIDA::IHistogram1D * gbldx4Histo;
    AIDA::IHistogram1D * gblrx4Histo;
    AIDA::IHistogram1D * gblpx4Histo;
    AIDA::IHistogram1D * gblqx4Histo;

    AIDA::IHistogram1D * gblax5Histo;
    AIDA::IHistogram1D * gbldx5Histo;
    AIDA::IHistogram1D * gblrx5Histo;
    AIDA::IHistogram1D * gblpx5Histo;
    AIDA::IHistogram1D * gblqx5Histo;

    AIDA::IHistogram1D * gblax6Histo;
    AIDA::IHistogram1D * gbldx6Histo;
    AIDA::IHistogram1D * gbldy6Histo;
    AIDA::IHistogram1D * gblrx6Histo;
    AIDA::IHistogram1D * gblry6Histo;
    AIDA::IHistogram1D * gblpx6Histo;
    AIDA::IHistogram1D * gblpy6Histo;
    AIDA::IHistogram1D * gblqx6Histo;
    AIDA::IHistogram1D * gblsx6Histo;
    AIDA::IHistogram1D * gbltx6Histo;

    AIDA::IHistogram1D * gblkx1Histo;
    AIDA::IHistogram1D * gblkx2Histo;
    AIDA::IHistogram1D * gblkx3Histo;
    AIDA::IHistogram1D * gblkx4Histo;
    AIDA::IHistogram1D * gblkx5Histo;
    AIDA::IHistogram1D * gblkx6Histo;

    AIDA::IHistogram1D * sixzx3Histo;
    AIDA::IHistogram1D * sixzy3Histo;
    AIDA::IHistogram1D * sixzx2Histo;
    AIDA::IHistogram1D * sixzy2Histo;
    AIDA::IHistogram1D * sixzx1Histo;
    AIDA::IHistogram1D * sixzy1Histo;

    AIDA::IHistogram1D * sixkxzyHisto;
    AIDA::IHistogram1D * sixkyzxHisto;
    AIDA::IHistogram1D * sixkxzxHisto;
    AIDA::IHistogram1D * sixkyzyHisto;

    AIDA::IHistogram1D * correvt100Histo;
    AIDA::IHistogram1D * correvt300Histo;
    AIDA::IHistogram1D * correvt1000Histo;
    AIDA::IHistogram1D * correvt4000Histo;
#endif


  };

  //! A global instance of the processor:
  //inline EUTelAnalysisCMSPixel aEUTelAnalysisCMSPixel;

}

#endif
