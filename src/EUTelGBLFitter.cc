/* 
 * File:   EUTelGBLFitter.cc
 * Contact: alexander.morton975@gmail.com
 * 
 */

#ifdef USE_GBL

// its own header 
#include "EUTelGBLFitter.h"

// eutelescope includes ".h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTELESCOPE.h"
#include "EUTelNav.h"

// marlin util includes
#include "mille/Mille.h"

// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TVector3.h"
#else
#error *** You need ROOT to compile this code.  *** 
#endif

// GBL
#include "include/GblTrajectory.h"
#include "include/GblPoint.h"
#include "include/GblData.h"
#include "include/BorderedBandMatrix.h"
#include "include/MilleBinary.h"
#include "include/VMatrix.h"


// system includes <>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <iomanip>
#include <iterator>
#include <algorithm>

namespace eutelescope {

	EUTelGBLFitter::EUTelGBLFitter() :
	_alignmentMode(0),
	_beamQ(-1),
	_eBeam(4.),
	_mEstimatorType(),
	_mille(0),
	_parameterIdXShiftsMap(),
	_parameterIdYShiftsMap(),
	_parameterIdZShiftsMap(),
	_parameterIdXRotationsMap(),
	_parameterIdYRotationsMap(),
	_parameterIdZRotationsMap(),
	_counter_num_pointer(1),
	_kinkAngleEstimation(false) //This used to determine if the correction matrix from the GBL fit is 5 or 7 elements long. 
	{}

	EUTelGBLFitter::~EUTelGBLFitter() {
	}
	//THIS IS THE SETTERS. Not simple get function but the action of all these functions it in the end to set a member variable to something.	
	//TO DO: This is the iterative alignment part that varies the precision of the measurement to allow the GBL tracks to be fitted. The precision matrix itself is an input by us. In the long run can we not do proper error calculations? 
	//TO DO: This does not have to be done for each state since it only varies per plane. So could input this data another way. However in the long run this may not be favourable since we could have correct error analysis eventually. 
	void EUTelGBLFitter::setMeasurementCov(EUTelState& state){
		double hitcov[]= {0,0,0,0};
		int izPlane = state.getLocation();
		if( _parameterIdXResolutionVec.size() > 0 && _parameterIdYResolutionVec.size() > 0 ){
			hitcov[0] = _parameterIdXResolutionVec[izPlane];
			hitcov[3] = _parameterIdYResolutionVec[izPlane];

			hitcov[0] *= hitcov[0]; 
			hitcov[3] *= hitcov[3];
		}
		else{//Default in case you have not specified a variance  
			throw(lcio::Exception("There is no measurement variance specified.")); 	
		}
		state.setCombinedHitAndStateCovMatrixInLocalFrame(hitcov);
	}
	//TO DO: Kinks are set to zero as inital guess this is ol for low radiation length enviroments.
	//Note that we take the planes themselfs at scatters and also add scatterers to simulate the medium inbetween. 
	void EUTelGBLFitter::setScattererGBL(gbl::GblPoint& point, EUTelState & state ) {
		streamlog_out(DEBUG1) << " setScattererGBL ------------- BEGIN --------------  " << std::endl;
		TVectorD scat(2); 
		scat[0] = 0.0; scat[1]=0.0;//TO DO: This will depend on if the initial track guess has kinks. Only important if we reiterate GBL track in high radiation length enviroments  
		TMatrixDSym precisionMatrix =  state.getScatteringVarianceInLocalFrame();
		streamlog_out(MESSAGE1) << "The precision matrix being used for the sensor  "<<state.getLocation()<<":" << std::endl;
		streamlog_message( DEBUG0, precisionMatrix.Print();, std::endl; );
		point.addScatterer(scat, precisionMatrix);
		streamlog_out(DEBUG1) << "  setScattererGBL  ------------- END ----------------- " << std::endl;
	}
		//This is used when the we know the radiation length already
		void EUTelGBLFitter::setScattererGBL(gbl::GblPoint& point,EUTelState & state, float variance) {
		streamlog_out(MESSAGE1) << " setScattererGBL ------------- BEGIN --------------  " << std::endl;
		TVectorD scat(2); 
		scat[0] = 0.0; scat[1]=0.0;  
		TMatrixDSym precisionMatrix =  state.getScatteringVarianceInLocalFrame(variance);
		streamlog_out(MESSAGE1) << "The precision matrix being used for the scatter:  " << std::endl;
		streamlog_message( DEBUG0, precisionMatrix.Print();, std::endl; );
		point.addScatterer(scat, precisionMatrix);
		streamlog_out(MESSAGE1) << "  setScattererGBL  ------------- END ----------------- " << std::endl;
	}
	void EUTelGBLFitter::setLocalDerivativesToPoint(gbl::GblPoint& point, EUTelState & state, float distanceFromKinkTargetToNextPlane){
		TMatrixD derivatives(2,2);
		derivatives.Zero();
		//The derivative is the distance since the change in the measurements is  deltaX = distanceFromkink*(Angle of kink)
		derivatives[0][0] = distanceFromKinkTargetToNextPlane; 
		derivatives[1][1] = distanceFromKinkTargetToNextPlane; 
		point.addLocals(derivatives);
	}
	//This will add measurement information to the GBL point
	//Note that if we have a strip sensor then y will be ignored using projection matrix.
	void EUTelGBLFitter::setMeasurementGBL(gbl::GblPoint& point, const double *hitPos,  double statePos[3], double combinedCov[4], TMatrixD projection){
		streamlog_out(DEBUG1) << " setMeasurementGBL ------------- BEGIN --------------- " << std::endl;
		TVectorD meas(2);//Remember we need to pass the same 5 since gbl expects this due to jacobian
		meas.Zero();
		meas[0] = hitPos[0] - statePos[0];//This is how residuals are shown as output after fit (measurement-prediction)
		meas[1] = hitPos[1] - statePos[1];
		TVectorD measPrec(2); 
		//TO DO: Can make the variance vector a matrix to explore relationships between x/y and variance.  
		measPrec[0] = 1. / combinedCov[0];	// cov(x,x)
		measPrec[1] = 1. / combinedCov[3];	// cov(y,y)
		streamlog_out(DEBUG4) << "This is what we add to the measured point:" << std::endl;
		streamlog_out(DEBUG4) << "Residuals and covariant matrix for the hit:" << std::endl;
		streamlog_out(DEBUG4) << "X:" << std::setw(20) << meas[0] << std::setw(20) << measPrec[0] <<"," << std::endl;
		streamlog_out(DEBUG4) << "Y:" << std::setw(20) << meas[1] << std::setw(20)  <<"," << measPrec[1] << std::endl;
		streamlog_out(DEBUG4) << "This H matrix:" << std::endl;
		streamlog_message( DEBUG0, projection.Print();, std::endl; );
		//The gbl library creates 5 measurement vector and 5x5 propagation matrix automatically. If  
		point.addMeasurement(projection, meas, measPrec, 0);//The last zero is the minimum precision before this is set to 0. TO DO:Remove this magic number
		streamlog_out(DEBUG1) << " setMeasurementsGBL ------------- END ----------------- " << std::endl;
	}

	//This function calculates the alignment jacobian and labels and attaches this to the point
	void EUTelGBLFitter::setAlignmentToMeasurementJacobian(EUTelTrack& track, std::vector< gbl::GblPoint >& pointList ){
		for (size_t i = 0;i<_vectorOfPairsMeasurementStatesAndLabels.size() ;++i ){
			if(_vectorOfPairsMeasurementStatesAndLabels.at(i).first.getTrackerHits().empty()){
				throw(lcio::Exception("One of the points on the list of measurements states has no hit."));
			}
			for(size_t j = 0; j < pointList.size(); ++j){
				if( _vectorOfPairsMeasurementStatesAndLabels.at(i).second == pointList.at(j).getLabel()){
					EUTelState state = _vectorOfPairsMeasurementStatesAndLabels.at(i).first;
					streamlog_out(DEBUG0)<<"Point needs global parameters. It has  "<<pointList.at(j).hasMeasurement()<<" measurements."<<std::endl;
					if(getLabelToPoint(pointList,_vectorOfPairsMeasurementStatesAndLabels.at(i).second).hasMeasurement() == 0){
						throw(lcio::Exception("This point does not contain a measurements. Labeling of the state must be wrong "));
					} 
					_MilleInterface->computeAlignmentToMeasurementJacobian(state);//We calculate the jacobian. 
					_MilleInterface->setGlobalLabels(state); //Get the correct label for the sensors x,y,z shift and rotations. Depending on alignment mode and sensor the plane is on 
					TMatrixD&  alignmentJacobian = _MilleInterface->getAlignmentJacobian();//Return what was calculated by computeAlignmentToMeasurementJacobian
					std::vector<int> labels =  _MilleInterface->getGlobalParameters();//Return what was set by setGlobalLabels
					streamlog_out(DEBUG0)<<"The state associated with this alignment jacobian:  "<<std::endl;
					streamlog_message( DEBUG0, state.print() ;, std::endl; );
					if(!state.getIsThereAHit()){
						throw(lcio::Exception("You are just about to use a state with no hit to determine alignment jacobian."));
					}
					streamlog_out(DEBUG0)<<"This is the alignment jacobian we are about to add to the point at location "<<state.getLocation()<<std::endl;
					streamlog_message( DEBUG0, alignmentJacobian.Print();, std::endl; );
					streamlog_out(DEBUG0)<<"Here are the labels for these alignment parameters for the state at location: "<<state.getLocation()<<std::endl;
					for(size_t k=0; k<labels.size();++k){
						streamlog_out(DEBUG0)<<"Label just before adding: "<<labels.at(k) <<std::endl;
					}
					pointList.at(j).addGlobals(labels, alignmentJacobian);
					streamlog_out(DEBUG0)<<"The number of global parameters for this point are "<<pointList[i].getNumGlobals()<<std::endl;
					for(size_t k=0; k<pointList.at(j).getGlobalLabels().size();++k){
						streamlog_out(DEBUG0)<<"Label just after adding: "<<pointList.at(j).getGlobalLabels().at(k) <<std::endl;
					}
					streamlog_out(DEBUG0)<<"The alignment matrix after adding to point: "<<std::endl;
					streamlog_message( DEBUG0, pointList.at(j).getGlobalDerivatives().Print() ;, std::endl; );
					break;
				}
				//streamlog_out(DEBUG0)<<"This point has "<<pointList.at(j).hasMeasurement()<<" measurements. It has label "<<pointList.at(j).getLabel() <<std::endl;
				if(j == (pointList.size()-1)){
					throw(lcio::Exception("There is no point associated with this state."));
				}
			}
		}
	}
	//This creates the scatters point to simulate the passage through air. 
	void EUTelGBLFitter::setPointListWithNewScatterers(std::vector< gbl::GblPoint >& pointList,EUTelState & state, std::vector<float>  variance){
		if(_scattererJacobians.size() != _scattererPositions.size()){
			throw(lcio::Exception("The size of the scattering positions and jacobians is different.")); 	
		}
		for(size_t i = 0 ;i < _scattererJacobians.size()-1;++i){//The last jacobain is used to get to the plane! So only loop over to (_scatterJacobians-1)
			gbl::GblPoint point(_scattererJacobians[i]);
			point.setLabel(_counter_num_pointer);
			_counter_num_pointer++;
			if(variance.at(i) == 0){
				throw(lcio::Exception("Variance is 0 for the scattering plane."));
			}

			setScattererGBL(point,state, variance.at(i));//TO DO:This will only work for homogeneous distributions.
			pointList.push_back(point);
		}
	}
	std::vector<float> EUTelGBLFitter::computeVarianceForEachScatterer(EUTelState & state , float percentageRadiationLength){
		const double scatteringVariance  = pow(Utility::getThetaRMSHighland(state.getBeamEnergy(), percentageRadiationLength),2);
		std::vector<float> variance;
		float powMeanStart=pow((_normalMean - _start),2);
		float denominator=_normalVariance + powMeanStart;
		variance.push_back(scatteringVariance*(_normalVariance/denominator));
		variance.push_back(scatteringVariance*(powMeanStart/denominator));
		return variance;
	}
	//This set the estimate resolution for each plane in the X direction.
	void EUTelGBLFitter::setParamterIdXResolutionVec( const std::vector<float>& vector)
	{
		//We have a similar check after this to see that number of planes and elements in resolution vector are the same. We need this here since if they are different then it will just give an exception from the vector tryign to access a element that does not exist.
		if ( geo::gGeometry().sensorZOrdertoIDs().size() != vector.size() ){
			streamlog_out( ERROR5 ) << "The number of planes: "<< geo::gGeometry().sensorZOrdertoIDs().size()<<"  The size of input resolution vector: "<<vector.size()  << std::endl;
			throw(lcio::Exception("The size of the resolution vector and the total number of planes is different for x axis."));
		}

		for(size_t i=0; i < geo::gGeometry().sensorZOrdertoIDs().size(); ++i){
			_parameterIdXResolutionVec[ geo::gGeometry().sensorZOrderToID(i)] = vector.at(i);
		}
	}
	//This sets the estimated resolution for each plane in the Y direction.
	void EUTelGBLFitter::setParamterIdYResolutionVec( const std::vector<float>& vector)
	{
		if ( geo::gGeometry().sensorZOrdertoIDs().size() != vector.size() ){
			streamlog_out( ERROR5 ) << "The number of planes: "<< geo::gGeometry().sensorZOrdertoIDs().size()<<"  The size of input resolution vector: "<<vector.size()  << std::endl;
			throw(lcio::Exception("The size of the resolution vector and the total number of planes is different for y axis."));
		}
		for(size_t i=0; i < geo::gGeometry().sensorZOrdertoIDs().size(); ++i){
			_parameterIdYResolutionVec[ geo::gGeometry().sensorZOrderToID(i)] = vector.at(i);
		}
	}

	//This is used to deal with downweighting of ouliers. You must provide as input t,h or c. This specifies the function that will be used to do the downweigting.
	//TO DO: Check that this function works as expected since it has never been used. 
	void EUTelGBLFitter::setMEstimatorType( const std::string& mEstimatorType ) {
		std::string mEstimatorTypeLowerCase = mEstimatorType;
		std::transform( mEstimatorType.begin(), mEstimatorType.end(), mEstimatorTypeLowerCase.begin(), ::tolower);//Make the character lower case
		if ( mEstimatorType.size() != 1 ) {
			streamlog_out( WARNING1 ) << "More than one character supplied as M-estimator option" << std::endl;
			streamlog_out( WARNING1 ) << "No M-estimator downweighting will be used" << std::endl;
			return;
		}
		//Compare character to the ones that are accepted and if one is the same then set out member variable equal to it.
		if ( mEstimatorTypeLowerCase.compare("t") == 0 ||
			mEstimatorTypeLowerCase.compare("h") == 0 ||
			mEstimatorTypeLowerCase.compare("c") == 0   ) this->_mEstimatorType = _mEstimatorType;
		else {
			streamlog_out( WARNING1 ) << "M-estimator option " << mEstimatorType << " was not recognized" << std::endl;
			streamlog_out( WARNING1 ) << "No M-estimator downweighting will be used" << std::endl;
		}
	}
	//Not all points are states so we need a way to do:
	//state->label->point. This is what this function does.  
	void EUTelGBLFitter::setPointVec( std::vector< gbl::GblPoint >& pointList, gbl::GblPoint& point) {
	point.setLabel(_counter_num_pointer);
	_counter_num_pointer++;
	pointList.push_back(point);
	streamlog_out(DEBUG0) << std::endl << "pushBackPoint size: " << pointList.size() <<  std::endl;
	}
	//We use the labels for trajectory and equate them to the pointList index.
	//We create them when creating the points. However trajectory will overwrite these. However what we write should be the same as trajectory.
	//Note this is used by track fitting so we need to associate ANY state to the correct GBL point.
	//It is important to note the different between  setPairAnyStateAndPointLabelVec and setPairMeasurementStateAndPointLabelVec. One is for any state even if it has no hit and the other is used only for states that have a hit.
	void EUTelGBLFitter::setPairAnyStateAndPointLabelVec(std::vector< gbl::GblPoint >& pointList, gbl::GblTrajectory* traj){
		streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::setPairAnyStateAndPointLabelVec-- BEGIN " << std::endl;
		_vectorOfPairsStatesAndLabels.clear();
		streamlog_out(DEBUG5)<<"The number of states is: "<< _statesInOrder.size()<<std::endl;
		std::vector<unsigned int> labels;
		traj->getLabels(labels);
		for(size_t i=0; i<_statesInOrder.size();++i){
			//We don't add +1 since this is taken into account by the labels.
			int threei = 3*i; //This is done since we always have two scatters between states
			streamlog_out(DEBUG0)<<"Pair (state,label)  ("<<  &(_statesInOrder.at(i))<<","<<labels.at(threei)<<")"<<std::endl; 
			_vectorOfPairsStatesAndLabels.push_back(std::make_pair(_statesInOrder.at(i), labels.at(threei)));	
		}
		streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::setPairAnyStateAndPointLabelVec- END " << std::endl;
	}
	//This code must be repeated since we need to create the ling between states and labels before trajectory in alignment
	//Note here we create the link between MEASUREMENT states and label. 
	void EUTelGBLFitter::setPairMeasurementStateAndPointLabelVec(std::vector< gbl::GblPoint >& pointList){
		streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::setPairMeasurementStateAndPointLabelVec-- BEGIN " << std::endl;
		_vectorOfPairsMeasurementStatesAndLabels.clear();
		streamlog_out(DEBUG5)<<"The number of measurement states is: "<< _measurementStatesInOrder.size()<<std::endl;
		size_t counter = 0;
		for(size_t i=0; i<pointList.size();++i){
			if(pointList.at(i).hasMeasurement()>0){//Here we assume that traj has ordered the pointList the same.
				streamlog_out(DEBUG0)<<"Measurement found! Pair (state,label)  ("<<  &(_measurementStatesInOrder.at(counter))<<","<<pointList.at(i).getLabel()<<")"<<std::endl; 
				_vectorOfPairsMeasurementStatesAndLabels.push_back(std::make_pair(_measurementStatesInOrder.at(counter), pointList.at(i).getLabel()));	
				counter++;
			}
			if(counter == (_measurementStatesInOrder.size()+1)){//Add one since you will make an extra loop
				throw(lcio::Exception("The counter is larger than the number of states saved as measurements."));
			}
		}
		if(counter !=  (_measurementStatesInOrder.size())){//Since we make an extra loop and counter++ again
			throw(lcio::Exception("We did not add all the states."));
		}
		streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::setPairMeasurementStateAndPointLabelVec------------ END " << std::endl;
	}

	// Convert input TrackCandidates and TrackStates into a GBL Trajectory
	// This is done using the geometry setup, the scattering and the hits + predicted states.
	void EUTelGBLFitter::setInformationForGBLPointList(EUTelTrack& track, std::vector< gbl::GblPoint >& pointList){
		streamlog_out(DEBUG4)<<"EUTelGBLFitter::setInformationForGBLPointList-------------------------------------BEGIN"<<std::endl;
		TMatrixD jacPointToPoint(5, 5);
		jacPointToPoint.UnitMatrix();
		//We place this variable here since we want to set it every new track to false and then true again after we get to the scattering plane.
		bool kinkAnglePlaneEstimationAddedNow=false;
		float distanceFromKinkTargetToNextPlane=0;
		for(size_t i=0;i < track.getStates().size(); i++){		
			streamlog_out(DEBUG3) << "The jacobian to get to this state jacobian on state number: " << i<<" Out of a total of states "<<track.getStates().size() << std::endl;
			streamlog_message( DEBUG0, jacPointToPoint.Print();, std::endl; );
			gbl::GblPoint point(jacPointToPoint);
			EUTelState state = track.getStates().at(i);
			EUTelState nextState;
			if(i != (track.getStates().size()-1)){//Since we don't want to propagate from the last state.
				nextState = track.getStates().at(i+1);
			}
			if(state.getLocation() == 271){
				kinkAnglePlaneEstimationAddedNow=true;
				//We set the distance to the next plane to use to set the local derivatives for the next point.
				distanceFromKinkTargetToNextPlane=state.getArcLengthToNextState();
			}else{
				setScattererGBL(point,state);//Every sensor will have scattering due to itself. 
			}
			_statesInOrder.push_back(state);//This is list of measurements states in the correct order. This is used later to associate ANY states with point labels
			if(state.getTrackerHits().size() == 0 ){
				streamlog_out(DEBUG3)  << "This state does not have a hit."<<std::endl;
				setPointVec(pointList, point);//This creates the vector of points and keeps a link between states and the points they created
			}else{
				double localPositionForState[]	= {state.getPosition()[0], state.getPosition()[1],state.getPosition()[2]};//Need this since geometry works with const doubles not floats 
				setMeasurementCov(state);
				double cov[4] ;
				state.getCombinedHitAndStateCovMatrixInLocalFrame(cov);
				setMeasurementGBL(point, state.getTrackerHits()[0]->getPosition(),  localPositionForState,  cov, state.getProjectionMatrix());
				//Here we set if we want to add local derivatives to points after we have found the plane we want to determine it's kink angle.
				if(kinkAnglePlaneEstimationAddedNow){
					setLocalDerivativesToPoint(point,state,distanceFromKinkTargetToNextPlane); 
					//We add the distance to the next plane so we get the total distance to the target.
					distanceFromKinkTargetToNextPlane=distanceFromKinkTargetToNextPlane + state.getArcLengthToNextState();
				}
				_measurementStatesInOrder.push_back(state);//This is list of measurements states in the correct order. This is used later to associate MEASUREMENT states with point labels in alignment
				setPointVec(pointList, point);
			}//End of else statement if there is a hit.

			if(i != (track.getStates().size()-1)){//We do not produce scatterers after the last plane
				//Note here to determine the scattering we use a straight line approximation between the two points the particle will travel through. However to determine were to place the scatterer we use the exact arc length. We do this since to change the TGeo radiation length would be a lot of work for a very small change. 
				const double stateReferencePoint[] = {state.getPosition()[0], state.getPosition()[1],state.getPosition()[2]};
				double globalPosSensor1[3];
				geo::gGeometry().local2Master(state.getLocation(),stateReferencePoint , globalPosSensor1 );
				const double nextStateReferencePoint[] = {nextState.getPosition()[0], nextState.getPosition()[1],nextState.getPosition()[2]};
				double globalPosSensor2[3];
				geo::gGeometry().local2Master(nextState.getLocation(),nextStateReferencePoint , globalPosSensor2 );
				testDistanceBetweenPoints(globalPosSensor1,globalPosSensor2);
				float percentageRadiationLength  = geo::gGeometry().findRadLengthIntegral(globalPosSensor1,globalPosSensor2, false );//TO DO: This adds the radiation length of the plane again. If you chose true the some times it returns 0. The could be the reason for the slightly small residuals. 
				if(percentageRadiationLength == 0){
					streamlog_out(MESSAGE9)<<"The positions between the scatterers are: "<<std::endl;
					streamlog_out(MESSAGE9)<<"Start: "<<globalPosSensor1[0]<<" "<<globalPosSensor1[1]<<" "<<globalPosSensor1[2]<<std::endl;
					streamlog_out(MESSAGE9)<<"End: "<<globalPosSensor2[0]<<" "<<globalPosSensor2[1]<<" "<<globalPosSensor2[2]<<std::endl;
					//We will not use events that have a radiation length of zero. 
					throw std::string("The radiation length was precisely zero. This is bad!");
				} 
				setMomentsAndStartEndScattering(state);
				findScattersZPositionBetweenTwoStates();//We use the exact arc length between the two states to place the scatterers. 
				jacPointToPoint=findScattersJacobians(state,nextState);
				std::vector<float> variance =  computeVarianceForEachScatterer(state ,percentageRadiationLength);
				setPointListWithNewScatterers(pointList,state, variance);//We assume that on all scattering planes the incidence angle is the same as on the last measurement state. Not a terrible approximation and will be corrected by GBL anyway.
			}else{
				streamlog_out(DEBUG3)<<"We have reached the last plane"<<std::endl;
			}
		} 
		streamlog_out(DEBUG4)<<"EUTelGBLFitter::setInformationForGBLPointList-------------------------------------END"<<std::endl;

	}

	//THIS IS THE GETTERS
	gbl::GblPoint EUTelGBLFitter::getLabelToPoint(std::vector<gbl::GblPoint> & pointList, int label)
	{
		for(size_t i = 0; i< pointList.size();++i)
		{
			streamlog_out(DEBUG0)<<"This point has label : "<<pointList.at(i).getLabel()<<". The label number: "<<label<<std::endl;
			if(pointList.at(i).getLabel() == label){
				return pointList.at(i);
			}
		}

		//if no match after loop this label must belong to no point 
		throw(lcio::Exception("There is no point with this label"));
	}
	//This used after trackfit will fill a map between (sensor ID and residualx/y). 
	void EUTelGBLFitter::getResidualOfTrackandHits(gbl::GblTrajectory* traj, std::vector< gbl::GblPoint > pointList,EUTelTrack& track, std::map< int, std::map< float, float > > &  SensorResidual, std::map< int, std::map< float, float > >& sensorResidualError ){
		for(size_t j=0 ; j< _vectorOfPairsMeasurementStatesAndLabels.size();j++){
			EUTelState state = _vectorOfPairsMeasurementStatesAndLabels.at(j).first;
			if(getLabelToPoint(pointList,_vectorOfPairsMeasurementStatesAndLabels.at(j).second).hasMeasurement() == 0){
				throw(lcio::Exception("This point does not contain a measurements. Labeling of the state must be wrong "));
			} 
			streamlog_out(DEBUG0) << std::endl << "There is a hit on the state. Hit pointer: "<< state.getTrackerHits()[0]<<" Find updated Residuals!" << std::endl;
			unsigned int numData; //Not sure what this is used for??????
			TVectorD aResiduals(2);
			TVectorD aMeasErrors(2);
			TVectorD aResErrors(2);
			TVectorD aDownWeights(2); 
			streamlog_out(DEBUG0)<<"To get residual of states we use label: "<<_vectorOfPairsMeasurementStatesAndLabels.at(j).second<<std::endl; 
			traj->getMeasResults(_vectorOfPairsMeasurementStatesAndLabels.at(j).second, numData, aResiduals, aMeasErrors, aResErrors, aDownWeights);
			streamlog_out(DEBUG0) <<"State location: "<<state.getLocation()<<" The residual x " <<aResiduals[0]<<" The residual y " <<aResiduals[1]<<std::endl;
			std::map<float, float> res; //This is create on the stack but will pass thisa by value to the new map so it 
			res.insert(std::make_pair(aResiduals[0],aResiduals[1]));
			SensorResidual.insert(std::make_pair(state.getLocation(), res));		
			std::map<float, float> resError; //This is create on the stack but will pass thisa by value to the new map so it 
			resError.insert(std::make_pair(aResErrors[0],aResErrors[1]));
			sensorResidualError.insert(std::make_pair(state.getLocation(), resError));	
		}
	}
	std::string EUTelGBLFitter::getMEstimatorType( ) const {
			return _mEstimatorType;
	}
	//COMPUTE
	void EUTelGBLFitter::computeTrajectoryAndFit(std::vector< gbl::GblPoint >& pointList,  gbl::GblTrajectory* traj, double* chi2, int* ndf, int & ierr){
		streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::computeTrajectoryAndFit-- BEGIN " << std::endl;
		double loss = 0.;
		streamlog_out ( DEBUG0 ) << "This is the trajectory we are just about to fit: " << std::endl;
		streamlog_message( DEBUG0, traj->printTrajectory(10);, std::endl; );
		streamlog_out ( DEBUG0 ) << "This is the points in that trajectory " << std::endl;
		streamlog_message( DEBUG0, traj->printPoints(10);, std::endl; );


		if ( !_mEstimatorType.empty( ) ) ierr = traj->fit( *chi2, *ndf, loss, _mEstimatorType );
		else ierr = traj->fit( *chi2, *ndf, loss );

		if( ierr != 0 ){
			streamlog_out(MESSAGE0) << "Fit failed!" << " Track error: "<< ierr << " and chi2: " << *chi2 << std::endl;
		}
		else{
		streamlog_out(MESSAGE0) << "Fit Successful!" << " Track error; "<< ierr << " and chi2: " << *chi2 << std::endl;
		}
		streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::computeTrajectoryAndFit -- END " << std::endl;
	}
	//TEST
	void EUTelGBLFitter::testUserInput(){
		if(_parameterIdXResolutionVec.size() != _parameterIdYResolutionVec.size()){
				throw(lcio::Exception("The vector for resolutions for X and Y are different sizes."));
		}
		if(_parameterIdXResolutionVec.size() != geo::gGeometry().nPlanes() ){
				throw(lcio::Exception("The total number of planes and the resolution of the planes vector are different sizes."));
		}
	}
	void EUTelGBLFitter::testTrack(EUTelTrack& track){
		streamlog_out(DEBUG4)<<"EUTelGBLFitter::testTrack------------------------------------BEGIN"<<std::endl;
		if(track.getStates().size() == 0 ){
			throw(lcio::Exception("The number of states is zero."));
		}
		///Note we do not use excluded planes here. This should be dealt with in pattern recognition.
		if (track.getNumberOfHitsOnTrack() > geo::gGeometry().nPlanes() ){
			throw(lcio::Exception("The number of hits on the track is greater than the number of planes.")); 	
		}
		streamlog_out(DEBUG5)<< "Input track passed tests!" <<std::endl;
		streamlog_out(DEBUG4)<<"EUTelGBLFitter::testTrack------------------------------------END"<<std::endl;

	} 
	void EUTelGBLFitter::testDistanceBetweenPoints(double* position1,double* position2){
		TVectorD displacement(3);
		displacement(0) = position1[0] - position2[0];
		displacement(1) = position1[1] - position2[1];
		displacement(2) = position1[2] - position2[2];
		float distance= sqrt(displacement.Norm2Sqr());
		streamlog_out(DEBUG0)<<"The distance between the points to calculate radiation length is: "<<distance<<std::endl;
		if(distance<1){
			throw(lcio::Exception("The distance between the states is too small.")); 	
		}
		if(distance>1000){//The planes should not be over a 1m in width
			throw(lcio::Exception("The distance between the planes is over 1m.")); 	
		}
	}

	//OTHER FUNCTIONS
	//We want to create a jacobain from (Plane1 -> scatterer1) then (scatterer1->scatterer2) then (scatter2->plane2). We return the last jacobain
	TMatrixD EUTelGBLFitter::findScattersJacobians(EUTelState state, EUTelState nextState){
		_scattererJacobians.clear();
		TVector3 position = state.getPositionGlobal();
		TVector3 momentum = state.computeCartesianMomentum();
		TVector3 newMomentum;
		int location = state.getLocation();
		int locationEnd = 314; //This will create a scatterer with normal in beam direction
		const gear::BField&   Bfield = geo::gGeometry().getMagneticField();
		gear::Vector3D vectorGlobal(position[0],position[1],position[1]);//Since field is homogeneous this seems silly but we need to specify a position to geometry to get B-field.
		const double Bx = (Bfield.at( vectorGlobal ).x());
		const double By = (Bfield.at( vectorGlobal ).y());
		const double Bz = (Bfield.at( vectorGlobal ).z());
		TVector3 B;
		B[0]=Bx; B[1]=By; B[2]=Bz;
		for(size_t i=0;i<_scattererPositions.size();i++){
			newMomentum = EUTelNav::getXYZMomentumfromArcLength(momentum, position,state.getBeamCharge(), _scattererPositions[i] );
			TMatrixD curvilinearJacobian = EUTelNav::getPropagationJacobianCurvilinear(_scattererPositions[i], state.getOmega(), momentum.Unit(),newMomentum.Unit());
			streamlog_out(DEBUG0)<<"This is the curvilinear jacobian at sensor : " << location << " or scatter: "<< i << std::endl; 
			streamlog_message( DEBUG0, curvilinearJacobian.Print();, std::endl; );
			TMatrixD localToCurvilinearJacobianStart =  EUTelNav::getLocalToCurvilinearTransformMatrix(momentum, location ,state.getBeamCharge() );
			streamlog_out(DEBUG0)<<"This is the local to curvilinear jacobian at sensor : " << location << " or scatter: "<< i << std::endl; 
			streamlog_message( DEBUG0, localToCurvilinearJacobianStart.Print();, std::endl; );
			TMatrixD localToCurvilinearJacobianEnd =  EUTelNav::getLocalToCurvilinearTransformMatrix(newMomentum,locationEnd ,state.getBeamCharge() );
			streamlog_out(DEBUG0)<<"This is the local to curvilinear jacobian at sensor : " << locationEnd << " or scatter: "<< i << std::endl; 
			streamlog_message( DEBUG0, localToCurvilinearJacobianEnd.Print();, std::endl; );
			TMatrixD curvilinearToLocalJacobianEnd = localToCurvilinearJacobianEnd.Invert();
			streamlog_out(DEBUG0)<<"This is the curvilinear to local jacobian at sensor : " << locationEnd << " or scatter: "<< i << std::endl; 
			streamlog_message( DEBUG0, curvilinearToLocalJacobianEnd.Print();, std::endl; );
			TMatrixD localToNextLocalJacobian = curvilinearToLocalJacobianEnd*curvilinearJacobian*localToCurvilinearJacobianStart;
			streamlog_out(DEBUG0)<<"This is the full jacobian : " << locationEnd << " or scatter: "<< i << std::endl; 
			streamlog_message( DEBUG0, localToNextLocalJacobian.Print();, std::endl; );
			_scattererJacobians.push_back(localToNextLocalJacobian);//To DO if scatter then plane is always parallel to z axis
			momentum[0]=newMomentum[0]; momentum[1]=newMomentum[1];	momentum[2]=newMomentum[2];
			location = 314;//location will always be a scatter after first loop.  
			if(i == (_scattererPositions.size()-2)){//On the last loop we want to create the jacobain to the next plane
				locationEnd = nextState.getLocation();
			}
		}
		if(_scattererJacobians.size() != 3){
			throw(lcio::Exception("There are not 3 jacobians produced by scatterers!")); 	
		}
		return _scattererJacobians.back();//return the last jacobian so the next state can use this
	}
	//The distance from the first state to the next scatterer and then from that scatterer to the next all the way to the next state. 
	//TO DO: This uses the optimum positions as described by Claus.  However for non homogeneous material distribution this might not be the case.
	void EUTelGBLFitter::findScattersZPositionBetweenTwoStates(){
		streamlog_out(DEBUG1) << "  findScattersZPositionBetweenTwoStates------------- BEGIN --------------  " << std::endl;
		_scattererPositions.clear();	
		streamlog_out(DEBUG1) << "The arc length to the next state is: " << _end << std::endl;
		//We place the first scatter to model the air just after the plane
		_scattererPositions.push_back(_start);//Z position of 1st scatterer	
		float secondScatterPosition = _normalMean +_normalVariance/(_normalMean-_start);
//		cout<<"Scat position " << secondScatterPosition <<endl;
		if(secondScatterPosition < _start){
			throw(lcio::Exception("The distance of the second scatterer is smaller than the start. "));
		}
		_scattererPositions.push_back(secondScatterPosition);//Z position of 2nd scatterer
		if(secondScatterPosition > _end){
			streamlog_out(MESSAGE5) << "The second scatter distance: "<< secondScatterPosition <<". The distance of the arc length: " << _end  << std::endl;
			throw(lcio::Exception("The distance of the second scatterer is larger than the next plane. "));
		}
		_scattererPositions.push_back(_end-secondScatterPosition); 
		streamlog_out(DEBUG1) << "  findScattersZPositionBetweenTwoStates------------- END --------------  " << std::endl;
	}

	void EUTelGBLFitter::setMomentsAndStartEndScattering(EUTelState& state){
		_start = 0.6; //This is 600 micron from the centre of the sensor. 
		_end = state.getArcLengthToNextState();
		if(_end == 0){
			throw(lcio::Exception("The size of arc length to the next plane is 0"));
		}
		_normalMean = 0.5*(_end-_start);
//		cout<<"normal mean: " << _normalMean <<endl;
		if(_normalMean == 0){
			throw(lcio::Exception("The mean of the scattering integral is zero. "));
		}
		_normalVariance = ((1.0/3.0)*(pow(_end,3)-pow(_start,3))-_normalMean*(pow(_end,2)-pow(_start,2))+pow(_normalMean,2)*(_end-_start))/(_end-_start);
//		cout<<"normal variance: " << _normalVariance <<endl;
		if(_normalVariance == 0){
			throw(lcio::Exception("The variance of the scattering integral is zero. "));
		}
	}

	void EUTelGBLFitter::resetPerTrack() {
		_counter_num_pointer=1;
		_measurementStatesInOrder.clear();
		_statesInOrder.clear();

	}
	//This function will take the estimate track from pattern recognition and add a correction to it. This estimated track + correction is you final GBL track.
	void EUTelGBLFitter::updateTrackFromGBLTrajectory (gbl::GblTrajectory* traj, std::vector< gbl::GblPoint >& pointList,EUTelTrack &track, std::map<int, std::vector<double> > &  mapSensorIDToCorrectionVec){
		streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::UpdateTrackFromGBLTrajectory-- BEGIN " << std::endl;
		
		for(size_t i=0;i < track.getStatesPointers().size(); i++){//We get the pointers no since we want to change the track state contents		
			EUTelState* state = track.getStatesPointers().at(i);
			TVectorD corrections(5);
			TMatrixDSym correctionsCov(5);
			if(_kinkAngleEstimation){
				corrections.ResizeTo(7);
				correctionsCov.ResizeTo(7,7);
			}
			for(size_t j=0 ; j< _vectorOfPairsStatesAndLabels.size();++j){
				if(_vectorOfPairsStatesAndLabels.at(j).first == *state){
					streamlog_out(DEBUG0)<<"The loop number for states with measurements is: " << j << ". The label is: " << _vectorOfPairsStatesAndLabels.at(j).second <<std::endl; 
					streamlog_out(DEBUG0)<<"To update track we use label: "<<_vectorOfPairsStatesAndLabels.at(j).second<<std::endl; 
					//This part gets the track parameters.//////////////////////////////////////////////////////////////////////////
					traj->getResults(_vectorOfPairsStatesAndLabels.at(j).second, corrections, correctionsCov );
					streamlog_out(DEBUG3) << std::endl << "State before we have added corrections: " << std::endl;
					state->print();
					TVectorD newStateVec(5);
					streamlog_out(DEBUG0)<<"-----------------------Corrections----------------------"<<std::endl;
					streamlog_out(DEBUG0)<< std::setw(20)<<"State number "<< i <<" of "  << track.getStatesPointers().size() << std::endl;  
					streamlog_out(DEBUG0)<< std::setw(20)<<"Track Parameter"<<std::setw(20)<<"| Intial value |" <<std::setw(20)<<"| correction to add |"  << std::endl;  
					streamlog_out(DEBUG0)<<std::setw(20)<<"Omega(Charge/Energy) "<<std::setw(15)<<state->getOmega()<<std::setw(20)<< corrections[0]<<std::endl;
					streamlog_out(DEBUG0)<<std::setw(20)<<"Intersection XY "<<std::setw(15)<<state->getIntersectionLocalXZ() <<std::setw(20)<< corrections[1]<<std::endl;
					streamlog_out(DEBUG0)<<std::setw(20)<<"Intersection YZ "<<std::setw(15)<<state->getIntersectionLocalYZ() <<std::setw(20)<< corrections[2]<<std::endl;
					streamlog_out(DEBUG0)<<std::setw(20)<<"PositionX "<<std::setw(15)<<state->getPosition()[0] <<std::setw(20)<< corrections[3]<<std::endl;
					streamlog_out(DEBUG0)<<std::setw(20)<<"PositionY "<<std::setw(15)<<state->getPosition()[1] <<std::setw(20)<< corrections[4]<<std::endl;

					newStateVec[0] = state->getOmega() + corrections[0];
					newStateVec[1] = state->getIntersectionLocalXZ()+corrections[1];
					newStateVec[2] = state->getIntersectionLocalYZ()+corrections[2];
					newStateVec[3] = state->getPosition()[0]+corrections[3];
					newStateVec[4] = state->getPosition()[1]+corrections[4]; 
					//Here we collect the total correction for every track and state. This we use to work out an average to make sure that on each iteration the track corrections are reducing
					_omegaCorrections = _omegaCorrections + corrections[0];	
					_intersectionLocalXZCorrections= _intersectionLocalXZCorrections+ corrections[1];	
					_intersectionLocalYZCorrections = _intersectionLocalYZCorrections + corrections[2];	
					_localPosXCorrections = _localPosXCorrections + corrections[3];
					_localPosYCorrections = _localPosYCorrections + corrections[4];

					std::vector<double> correctionVec;
					correctionVec.push_back(corrections[0]);	
					correctionVec.push_back(corrections[1]);	
					correctionVec.push_back(corrections[2]);	
					correctionVec.push_back(corrections[3]);	
					correctionVec.push_back(corrections[4]);	

					mapSensorIDToCorrectionVec[state->getLocation()] = correctionVec;
					////////////////////////////////////////////////////////////////////////////////////////////////
					state->setStateVec(newStateVec);
					//If the state is just a normal scatterer then we get the scattering results from the getScatResults() function. Otherwise we use the local derivative to determine the kink angles.
					if(state->getLocation() != 271){
						//Note that this says meas but it simple means that the states you access must be a scatterer. Since all our planes are scatterers then we can access all of them.
						//However the only planes we are interest in are the ones that don't simply model the air.
						unsigned int numData;
						TVectorD aResidualsKink(2);//Measurement - Prediction
						TVectorD aMeasErrorsKink(2);
						TVectorD aResErrorsKink(2);
						TVectorD aDownWeightsKink(2); 
						traj->getScatResults( _vectorOfPairsStatesAndLabels.at(j).second, numData, aResidualsKink, aMeasErrorsKink, aResErrorsKink, aDownWeightsKink);
						state->setKinks(aResidualsKink);
						streamlog_out(DEBUG3) << std::endl << "State after we have added corrections: " << std::endl;
						state->print();
					}else{
						TVectorD correctionsKinks(5);
						TMatrixDSym correctionsCovKinks(5);
						if(_kinkAngleEstimation){
							corrections.ResizeTo(7);
							correctionsCov.ResizeTo(7,7);
						}
						//TO DO: This will only work for 3 planes in the forward region for scattering measurements.
						traj->getResults(_vectorOfPairsStatesAndLabels.at(j+1).second, correctionsKinks, correctionsCovKinks );
						TVectorD kinks(2);//Measurement - Prediction
						kinks[0] = correctionsKinks[5];
						kinks[1] = correctionsKinks[6];

						state->setKinks(kinks);
					}
					break;
				}
			}//END of loop of all states with hits	
		}//END of loop of all states
	}

} // namespace eutelescope


#endif
