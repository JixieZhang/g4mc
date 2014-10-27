// ********************************************************************
//
// $Id: LERDDetectorConstruction.hh,v 1.00, 2014/10/24 LERD Exp $
// --------------------------------------------------------------
//
//Construct LERD detector
//This is not an independent class, it has to be invoke by class HRSDetectorConstruction

#ifndef LERDDetectorConstruction_H_
#define LERDDetectorConstruction_H_ 1

#include "globals.hh"	//for units and g4types and g4io
#include "HRSVisAttribute.hh"
#include "G4VUserDetectorConstruction.hh"
#include "HRSMaterial.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

class LERDDetectorConstruction : public G4VUserDetectorConstruction, public HRSVisAttribute
{
public:
	LERDDetectorConstruction(G4LogicalVolume *mother=0);
	~LERDDetectorConstruction();

	// the function which builds everything
	G4VPhysicalVolume* Construct();
	G4VPhysicalVolume* ConstructSolenoid(G4LogicalVolume *pMotherLogVol);

private:
	void GetConfig();
	void ConstructMaterials();

private:
	HRSMaterial* mMaterialManager;
	G4LogicalVolume* mMotherLogVol;
	double mTargetXOffset,mTargetYOffset,mTargetZOffset;

private:

	G4Material* heliumGas;
	G4Material* WCGas;	//C4H10
	G4Material* deuteriumGas;

	G4Material* H2TgGas;
	G4Material* D2TgGas;
	G4Material* He3TgGas;
	G4Material* He4TgGas;
	

	G4Material* targetMaterial;
	G4Material* targetWallMaterial;

	G4double mD2GasD,mHeGasD,mWCGasD;

	G4double mD2GasL,mD2GasR,mD2GasT,mD2GasP;
	G4double mHeGasT,mHeGasP;
	G4double mWCGasT,mWCGasP;

	G4int mTargetType;  		//1=H2, 2=D2,3=He3,4=He4
	G4int mTgWallMaterialType;  	//1 is kapton, 2 is aluminum
	G4double mLERDLength,mTgWallThick;
	G4double m1stWCR,m1stWCWindowThick,m1stWCThick;
	G4double m2ndWCR,m2ndWCWindowThick,m2ndWCThick;
	
	G4double mSiXR,mSiYR,mSCR,mLERDChamberR;
	G4double mSiXThick,mSiYThick,mSCThick,mLERDChamberThick;

	G4double mBStepLimit,mDCStepLimit;


	G4int mSetupSolenoid;
	G4double mSolenoidPosX,mSolenoidPosY,mSolenoidPosZ;

	G4int mSetupEntranceNExitCap,mSetupEndPlateNCover;

};

#endif  //LERDDetectorConstruction_H_

