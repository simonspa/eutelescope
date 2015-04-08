#include "CMSPixelSingle.h"

namespace eutelescope {

  namespace geo {

    CMSPixelSingle::CMSPixelSingle(): EUTelGenericPixGeoDescr(8.1, 8.1, 0.285, //size X, Y, Z
							      0, 51, 0, 79, //min max X,Y
							      93.660734) //rad length					
    {
      //Create the material for the sensor
      matSi = new TGeoMaterial( "Si", 28.0855 , 14.0, 2.33, _radLength, 45.753206 );
      Si = new TGeoMedium("CMSPixelSilicon",1, matSi);

      /* Make a box for the sensitive area
	 Size is: x=2*300+50*150=8100 microns and y=1*200 + 79*100=8000 microns
	 MakeBox takes the half of those values in mm as arguments */
      plane = _tGeoManager->MakeBox("sensarea_cmspix", Si, 4.05, 4.05, 0.1425);

      //Create volumes for the centre and edge regions (again: half dimensions...)
      // The Fiducial volume of 50x79 pixels:
      TGeoVolume* centreregion = _tGeoManager->MakeBox("cmspixcentreregion" ,Si, 3.75, 3.95, 0.1425);
      // The left and right edge pixels (X = 0,51):
      TGeoVolume* xedgeregion   = _tGeoManager->MakeBox("cmspixedgeregion"   ,Si, 0.15 , 3.95, 0.1425);
      // The upper egde pixels (Y = 79):
      TGeoVolume* yedgeregion   = _tGeoManager->MakeBox("cmspixyedgeregion"  ,Si, 4.05 , 0.1, 0.1425);

      //Divide the regions to create pixels
      // The left and right pixels: divide along Y axis (2) in 79 blocks (the upper row
      // belongs to the yedge.
      xedgeregion->Divide("cmspixedgepixel",   2, 79, 0, 1, 0, "N");
      // The upper edge: divide along X (1) in 52 pixels.
      yedgeregion->Divide("cmspixyedgepixel",  1, 52, 0, 1, 0, "N");
      // This somehow neglects the size*4 pixels in the corners. FIXME

      // Divide fiducial into 79 rows of 50 pixels:
      TGeoVolume* centrerow = centreregion->Divide("cmspixcentrerow", 2, 79, 0, 1, 0, "N");
      centrerow ->Divide("cmspixcentrepixel", 1,  50, 0, 1, 0, "N"); 

      //And place them to make a singlechip
      // Move them by half the center width + half their own width:
      plane->AddNode(centreregion, 1, new TGeoTranslation( 0.0, 0, 0));
      plane->AddNode(xedgeregion,   2, new TGeoTranslation(-3.9, 0, 0));
      plane->AddNode(xedgeregion,   3, new TGeoTranslation( 3.9, 0, 0));
      // TGeo starts with (0,0) in upper left corner, so move to negative values:
      plane->AddNode(yedgeregion,   4, new TGeoTranslation(0, -4.05, 0));

    }

    CMSPixelSingle::~CMSPixelSingle()
    {
      //delete matSi;
      //delete Si;
    }

    void  CMSPixelSingle::createRootDescr(char const * planeVolume)
    {
      //Get the plane as provided by the EUTelGeometryTelescopeGeoDescription
      TGeoVolume* topplane =_tGeoManager->GetVolume(planeVolume);
      //Finaly add the sensitive area to the plane
      topplane->AddNode(plane, 1);
    }

    std::string CMSPixelSingle::getPixName(int x , int y)
    {
      char buffer [100];

      //since pixel 0|0 is located on the upper left corner we have to correct y by 335-y+1 
      //(one for the offset in TGeo which starts counting at 1)

      // Edge regions:
      if(y == 79) {
	snprintf(buffer, 100, "/sensarea_cmspix_1/cmspixyedgeregion_3/cmspixedgepixel_%d", x+1);
      }

      if(x == 0) {
	snprintf(buffer, 100, "/sensarea_cmspix_1/cmspixedgeregion_1/cmspixedgepixel_%d", 80-y);
      }
      else if(x == 51) {
	snprintf(buffer, 100, "/sensarea_cmspix_1/cmspixedgeregion_2/cmspixedgepixel_%d", 80-y);
      }

      // Fiducial volume (center region):
      if(x > 0 && x < 51) {
	snprintf( buffer, 100, "/sensarea_cmspix_1/cmspixcentreregion_1/cmspixcentrerow_%d/cmspixcentrepixel_%d", 80-y, x);
      }

      //Return the full path
      return std::string( buffer ); 
    }

    //TODO: parse the path to a pixel number!
    std::pair<int, int>  CMSPixelSingle::getPixIndex(char const*){return std::make_pair(0,0); }

    EUTelGenericPixGeoDescr* maker() {
      CMSPixelSingle* mPixGeoDescr = new CMSPixelSingle();
      return dynamic_cast<EUTelGenericPixGeoDescr*>(mPixGeoDescr);
    }

  } //namespace geo
} //namespace eutelescope

