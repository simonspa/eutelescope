#ifndef CMSPIXELSINGLE_H
#define	CMSPIXELSINGLE_H

  /** @class CMSPixelSingle
	* This class is the implementation of  @class EUTelGenericPixGeoDescr
	* for a CMSPixel layout with edge pixels which have double width, the
	* matrix has 52 x 80 pixels, 150 x 100 microns**2 size
	* with exception of the edge pixels (X=0,51 Y=79) which are longer.
    */

//EUTELESCOPE
#include "EUTelGenericPixGeoDescr.h"

//ROOT
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"

namespace eutelescope {
namespace geo {

class CMSPixelSingle : public EUTelGenericPixGeoDescr {
	
	public:
		CMSPixelSingle();
		~CMSPixelSingle();

		void createRootDescr(char const *);
		std::string getPixName(int, int);
		std::pair<int, int> getPixIndex(char const *);

	protected:
		TGeoMaterial* matSi;
		TGeoMedium* Si;
		TGeoVolume* plane;

};

} //namespace geo
} //namespace eutelescope

#endif	//CMSPIXELSINGLE_H
