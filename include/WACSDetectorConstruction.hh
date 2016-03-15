// ********************************************************************
// $Id: HallCBeamLineConstruction.hh,v 1.00, 2016/01/12 WACS Exp $
// --------------------------------------------------------------
//

#ifndef WACSDetectorConstruction_h
#define WACSDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "HRSVisAttribute.hh"
#include "HRSMaterial.hh"


class HRSEMFieldSetup;

class G4VPhysicalVolume;
class G4VSensitiveDetector;
class G4LogicalVolume;

class WACSDetectorConstruction : public G4VUserDetectorConstruction, public HRSVisAttribute
{
public:
	WACSDetectorConstruction(G4LogicalVolume *mother=0);
	virtual ~WACSDetectorConstruction();

public:
	virtual G4VPhysicalVolume* Construct();
	G4VPhysicalVolume* ConstructWACSScatChamber(G4LogicalVolume* motherLogical);
	G4VPhysicalVolume* ConstructWACSTarget(G4LogicalVolume* motherLogical);
	G4VPhysicalVolume* ConstructWACSChicane(G4LogicalVolume* motherLogical);
	G4VPhysicalVolume* ConstructWACSPlatform(G4LogicalVolume* motherLogical);
	

private:
	void ConstructMaterial();
	void GetConfig();

private:
	
	HRSMaterial* mMaterialManager;
	G4LogicalVolume *mMotherLogVol;

	G4Material *WACS_NH3He;


private:

	//////////////////////////
	//the following can be found in Detector_G2P.ini
	double mScatChamberRin,mScatChamberRout,mScatChamberL;
	double mScatChamberExitWindowThick;

	double mShieldLN2Rin,mShieldLN2Rout,mShieldLN2WindowThick,mShieldLN2L;
	double mShieldLHeRin,mShieldLHeRout,mShieldLHeL;

	int    mTargetType,mSetupG2PTarget;
	double mTargetL;

	int    mSetupG2PScatChamber,mSetupCoil;

	int    mSetupBeamDump;
	double mBeamDumpWidth,mBeamDumpHeight,mBeamDumpThick;
        double mPivot2BeamDumpX, mPivot2BeamDumpY, mPivot2BeamDumpZ;
        double mCollimatorDiameter;

	int    mSetupChicane,mSetupChicaneVD;

	int    mSetupPlatform;

	
	//////////////////////////
	//the following can be found in Detector_WACS.ini or Detecotor.ini
	double mPivotXOffset,mPivotYOffset,mPivotZOffset;
	double mScatChamberXOffset,mScatChamberYOffset,mScatChamberZOffset;
	double mTargetXOffset,mTargetYOffset,mTargetZOffset;

};

#endif

