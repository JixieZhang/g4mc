// ********************************************************************
// $Id: HallCBeamLineConstruction.cc,v 1.00, 2016/01/12 WACS Exp $
//
// ********************************************************************
//
#include <stdio.h>
#include <math.h>
#include <iostream>
using namespace std;

#include "WACSDetectorConstruction.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trap.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "G4PVReplica.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Polycone.hh"
#include "G4AssemblyVolume.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4UserLimits.hh"

#include "HRSStdSD.hh"
#include "HRSDCSD.hh"
#include "UsageManager.hh"
#include "HRSMaterial.hh"
#include "HRSEMFieldSetup.hh"

//To verify some geometry, I use this flag to place some unit
//just for debugging
//#define G4DEBUG_GEOMETRY 0

//The G4RotationMatrix rotates the whole coordinate system, 
//Looking from top, it always rotate clockwise  

extern UsageManager* gConfig;	

/////////////////////////////////////////////////////////////////////
WACSDetectorConstruction::WACSDetectorConstruction(G4LogicalVolume *mother) : 
mMotherLogVol(mother) 
{
	GetConfig();
	mMaterialManager=HRSMaterial::GetHRSMaterialManager();
	
	ConstructMaterial();
	G4cout<<"Contrstruct WACS geometry ... done! "<<G4endl;
}

WACSDetectorConstruction::~WACSDetectorConstruction()
{
	if (WACS_NH3He)  delete WACS_NH3He;
	G4cout<<"Delete WACS geometry ... done! "<<G4endl;
}

/////////////////////////////////////////////////////////////////////
void WACSDetectorConstruction::GetConfig()
{
	gConfig->ReadFile("Detector_WACS.ini");
	//////////////////////////////////////////////////////////////////
	//the following is just for wacs

	gConfig->GetParameter("ScatChamberRin",mScatChamberRin);
	mScatChamberRin*=mm;
	gConfig->GetParameter("ScatChamberRout",mScatChamberRout);
	mScatChamberRout*=mm;
	gConfig->GetParameter("ScatChamberL",mScatChamberL);
	mScatChamberL*=mm;

	gConfig->GetParameter("ScatChamberExitWindowThick",mScatChamberExitWindowThick);
	mScatChamberExitWindowThick*=mm;
	gConfig->GetParameter("ShieldLN2WindowThick",mShieldLN2WindowThick);
	mShieldLN2WindowThick*=mm;

	gConfig->GetParameter("ShieldLN2Rin",mShieldLN2Rin);
	mShieldLN2Rin*=mm;
	gConfig->GetParameter("ShieldLN2Rout",mShieldLN2Rout);
	mShieldLN2Rout*=mm;
	gConfig->GetParameter("ShieldLN2L",mShieldLN2L);
	mShieldLN2L*=mm;
	gConfig->GetParameter("ShieldLHeRin",mShieldLHeRin);
	mShieldLHeRin*=mm;
	gConfig->GetParameter("ShieldLHeRout",mShieldLHeRout);
	mShieldLHeRout*=mm;
	gConfig->GetParameter("ShieldLHeL",mShieldLHeL);
	mShieldLHeL*=mm;

	gConfig->GetParameter("TargetType",mTargetType);
	mSetupG2PTarget=1;
	gConfig->GetParameter("SetupG2PTarget",mSetupG2PTarget);

	gConfig->GetParameter("TargetL",mTargetL);
	mTargetL*=mm;

	gConfig->GetParameter("SetupG2PScatChamber",mSetupG2PScatChamber);
	gConfig->GetParameter("SetupCoil",mSetupCoil);

	gConfig->GetParameter("SetupBeamDump",mSetupBeamDump);
	gConfig->GetParameter("BeamDumpWidth",mBeamDumpWidth);
	mBeamDumpWidth*=mm;
	gConfig->GetParameter("BeamDumpHeight",mBeamDumpHeight);
	mBeamDumpHeight*=mm;
	gConfig->GetParameter("BeamDumpThick",mBeamDumpThick);
	mBeamDumpThick*=mm;
	gConfig->GetParameter("Pivot2BeamDumpZ",mPivot2BeamDumpZ);
	mPivot2BeamDumpZ*=mm;
	gConfig->GetParameter("Pivot2BeamDumpY",mPivot2BeamDumpY);
	mPivot2BeamDumpY*=mm;


	gConfig->GetParameter("SetupChicane",mSetupChicane);
	gConfig->GetParameter("SetupChicaneVD",mSetupChicaneVD);

	mSetupPlatform=0;
	gConfig->GetParameter("SetupPlatform",mSetupPlatform);
	
	////////////////////////////////////////////////////////////////
	gConfig->GetParameter("PivotXOffset",mPivotXOffset);
	mPivotXOffset*=mm;
	gConfig->GetParameter("PivotYOffset",mPivotYOffset);
	mPivotYOffset*=mm;
	gConfig->GetParameter("PivotZOffset",mPivotZOffset);
	mPivotZOffset*=mm;

	gConfig->GetParameter("ScatChamberXOffset",mScatChamberXOffset);
	mScatChamberXOffset*=mm;
	gConfig->GetParameter("ScatChamberYOffset",mScatChamberYOffset);
	mScatChamberYOffset*=mm;
	gConfig->GetParameter("ScatChamberZOffset",mScatChamberZOffset);
	mScatChamberZOffset*=mm;
	gConfig->GetParameter("TargetXOffset",mTargetXOffset);
	mTargetXOffset*=mm;
	gConfig->GetParameter("TargetYOffset",mTargetYOffset);
	mTargetYOffset*=mm;
	gConfig->GetParameter("TargetZOffset",mTargetZOffset);
	mTargetZOffset*=mm;

	G4cout<<"\n****Load WACS detector config parameters done!***"<<G4endl;

	return ;
}

void WACSDetectorConstruction::ConstructMaterial()
{
	//Due to the fact that the HRSMaterial is built before reading Detector_WACS.ini
	//The mMaterialManager->NH3He was built with 55% NH3 volumn ratio. To allow user 
	//to change the packing ratio, I have to create the material here.

	double pSolidNH3D;	
	gConfig->GetParameter("SolidNH3D",pSolidNH3D);
	pSolidNH3D*=mg/cm3;

	double pLiquidHeD;	
	gConfig->GetParameter("LiquidHeD",pLiquidHeD);
	pLiquidHeD*=mg/cm3;

	double pNH3WeightRatio;	
	gConfig->GetParameter("NH3WeightRatio",pNH3WeightRatio);
	
	double density;
	density = 1.0/((1.0-pNH3WeightRatio)/pLiquidHeD+pNH3WeightRatio/pSolidNH3D);
	WACS_NH3He = new G4Material("WACS_NH3He", density, 2);
	WACS_NH3He->AddMaterial(mMaterialManager->solidNH3, pNH3WeightRatio);
	WACS_NH3He->AddMaterial(mMaterialManager->liquidHe, 1.0-pNH3WeightRatio);
}

/////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* WACSDetectorConstruction::Construct()
{	
	G4VPhysicalVolume* thePhysVol=0;

	/////////////////////////
	// WACS chicane magnets
	/////////////////////////
	if(mSetupChicane) this->ConstructWACSChicane(mMotherLogVol);

	/////////////////////////
	//g2p scatter chamber 
	/////////////////////////
	//since the coil will always rotate with the chamber, I use a container to hold them.
	//the ScatChamberContainer start from mShieldLHeRin to 
	//mScatChamberRout+pScatChamberEntranceWindow2UnionL+10*mm
	//The advantage is that all staff inside the containner have the same rotation 
	if(mSetupG2PScatChamber) this->ConstructWACSScatChamber(mMotherLogVol);


	/////////////////////////
	//g2p target container and target vessel
	/////////////////////////	
	if(mSetupG2PTarget) this->ConstructWACSTarget(mMotherLogVol);

	
	/////////////////////////
	// the platform
	/////////////////////////
	if(mSetupPlatform)  this->ConstructWACSPlatform(mMotherLogVol);

	return thePhysVol;
}

/////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* WACSDetectorConstruction::ConstructWACSScatChamber(G4LogicalVolume* motherLogical)
{
	const double inch=2.54*cm;
	double startphi,deltaphi;

	G4RotationMatrix* pRotX90deg = new G4RotationMatrix();
	pRotX90deg->rotateX(90.*deg);
	G4RotationMatrix* pRotX180deg = new G4RotationMatrix();
	pRotX180deg->rotateX(180.*deg);
	G4RotationMatrix* pRotX270deg = new G4RotationMatrix();
	pRotX270deg->rotateX(270.*deg);

	G4RotationMatrix* pRotY90deg = new G4RotationMatrix();
	pRotY90deg->rotateY(90.*deg);
	G4RotationMatrix* pRotY180deg = new G4RotationMatrix();
	pRotY180deg->rotateY(180.*deg);
	G4RotationMatrix* pRotY270deg = new G4RotationMatrix();
	pRotY270deg->rotateY(270.*deg);

	/////////////////////////////////////////////////////////
	//By default, the scatter chamber and field are in longitudinal status,
	//For transverse configuration, we must rotate the scatter chamber and coil
	//about Y axis by -90 degrees, minus here means anti-clockwise in the overlook view 

	//Note that the scatter chamber is put in a vertical orientation
	//By Jixie: I have verified that rotate by matrix A then by Matrix B should be written
	//as Rotate=A*B, always put the first rotation at the left

#ifdef G4DEBUG_GEOMETRY
	if(G4DEBUG_GEOMETRY>=1)
	{
		G4RotationMatrix *pRotBField=new G4RotationMatrix();
		BField_Helm *pHelmField=BField_Helm::GetInstance();
		pRotBField=(G4RotationMatrix*)(pHelmField->GetRotation_L2F());

		G4RotationMatrix* pRotX90deg = new G4RotationMatrix();
		pRotX90deg->rotateX(90.*deg);
		G4RotationMatrix *pRotScatChamber=new G4RotationMatrix();
		*pRotScatChamber = (*pRotX90deg) * (*pRotBField) ; //rotate with the coil	

		cout<<"\n pRotBField ==> EulerAngles: phi="<<pRotBField->getPhi()/deg
			<<"  theta="<<pRotBField->getTheta()/deg
			<<"  psi="<<pRotBField->getPsi()/deg<<endl;

		cout<<"\n pRotScatChamber ==> EulerAngles: phi="<<pRotScatChamber->getPhi()/deg
			<<"  theta="<<pRotScatChamber->getTheta()/deg
			<<"  psi="<<pRotScatChamber->getPsi()/deg<<endl;
	}
#endif

	// Magnetic field Rotation----------------------------------------------------------
	BField_Helm *pHelmField=BField_Helm::GetInstance();
	G4RotationMatrix *pRotBField=(G4RotationMatrix*)(pHelmField->GetRotation_L2F());
	G4RotationMatrix *pRotScatInHall=new G4RotationMatrix();
	*pRotScatInHall = (*pRotX90deg) * (*pRotBField) ; //rotate with the coil

	/////////////////////////////////////////////////////////
#ifdef G4DEBUG_GEOMETRY
	//this part is used to compare with pRotScatInHall, this is the way I figure out
	//how pRotScatInHall should be constructed

	G4RotationMatrix *pRotX90Z270deg=new G4RotationMatrix();
	pRotX90Z270deg->rotateX(90*deg);
	pRotX90Z270deg->rotateZ(270*deg); 

	G4RotationMatrix *pRotX90Z353deg=new G4RotationMatrix();
	pRotX90Z353deg->rotateX(90*deg);
	pRotX90Z353deg->rotateZ(353*deg); 

	if(G4DEBUG_GEOMETRY>=1)
	{
		cout<<"\n pRotBField ==> EulerAngles: phi="<<pRotBField->getPhi()/deg
			<<"  theta="<<pRotBField->getTheta()/deg
			<<"  psi="<<pRotBField->getPsi()/deg<<endl;

		cout<<"\n pRotX90deg ==> EulerAngles: phi="<<pRotX90deg->getPhi()/deg
			<<"  theta="<<pRotX90deg->getTheta()/deg
			<<"  psi="<<pRotX90deg->getPsi()/deg<<endl;

		cout<<"\n pRotX90Z270deg ==> EulerAngles: phi="<<pRotX90Z270deg->getPhi()/deg
			<<"  theta="<<pRotX90Z270deg->getTheta()/deg
			<<"  psi="<<pRotX90Z270deg->getPsi()/deg<<endl;

		cout<<"\n pRotX90Z353deg ==> EulerAngles: phi="<<pRotX90Z353deg->getPhi()/deg
			<<"  theta="<<pRotX90Z353deg->getTheta()/deg
			<<"  psi="<<pRotX90Z353deg->getPsi()/deg<<endl;

		//The following is used to verify that rotate by matrix A then rotate by Matrix B should 
		//be written as Rotate=A*B, always put the first rotation at the left
		cout<<"\n (*pRotX90deg) * (*pRotBField) ==> EulerAngles: phi="<<pRotScatInHall->getPhi()/deg
			<<"  theta="<<pRotScatInHall->getTheta()/deg
			<<"  psi="<<pRotScatInHall->getPsi()/deg<<endl;
		//to compare AB with BA	
		G4RotationMatrix pRot_BA=(*pRotBField) * (*pRotX90deg);
		cout<<"\n (*pRotBField) * (*pRotX90deg) ==> EulerAngles: phi="<<pRot_BA.getPhi()/deg
			<<"  theta="<<pRot_BA.getTheta()/deg
			<<"  psi="<<pRot_BA.getPsi()/deg<<endl;

		cout<<"\n pRotScatInHall ==> EulerAngles: phi="<<pRotScatInHall->getPhi()/deg
			<<"  theta="<<pRotScatInHall->getTheta()/deg
			<<"  psi="<<pRotScatInHall->getPsi()/deg<<endl;
	}
#endif
	/////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////
	//scattering chamber container
	/////////////////////////////////////////////////////////
	//This is just a container to enclose the taraget chamber, 10 mm larger than the 
	//scattering chamber itself.
	//The scattering chamber containner is made of helium gas, 
	//will put vacuum inside using innerScatChamber
	//double pScatChamberContainerR=(479.425+0.508+40+10.0)*mm, mScatChamberL=53.5*inch;
	double pSCEntranceWindowLongFlangeL=40.0*mm;		//the thickness of the flange
	double pScatChamberContainerRin=mShieldLHeRin;      //double mShieldLHeRin=38.1*mm
	double pScatChamberContainerRout=mScatChamberRout+mScatChamberExitWindowThick+
		pSCEntranceWindowLongFlangeL+10.0*mm;
	double pScatChamberContainerL=mScatChamberL+(3.50+17.0+1.25)*inch*2+10.0*mm;
	G4VSolid* scatChamberContainerExtendedSolid = new G4Tubs("scatChamberContainerExtendedTubs",
		pScatChamberContainerRin,pScatChamberContainerRout,
		pScatChamberContainerL/2.0,0.,360.*deg);
	G4VSolid* scatChamberContainerExtraSolid = new G4Tubs("scatChamberContainerExtraTubs",
		0,pScatChamberContainerRout+1*mm,17.25*inch/2.0,0.,360.*deg);
	G4SubtractionSolid* scatChamberContainerSolid=new G4SubtractionSolid("scatChamberContainerSolid",
		scatChamberContainerExtendedSolid,scatChamberContainerExtraSolid,
		0,G4ThreeVector(0,0,-mScatChamberL/2-17.25*inch/2.0-4.5*inch-10.1*mm));

	G4LogicalVolume* scatChamberContainerLogical = new G4LogicalVolume(scatChamberContainerSolid,
		mMaterialManager->heliumGas,"scatChamberContainerLogical",0,0,0);
	G4VPhysicalVolume* scatChamberContainerPhys=new G4PVPlacement(pRotScatInHall,
		G4ThreeVector(mScatChamberXOffset,mScatChamberYOffset,mScatChamberZOffset),
		scatChamberContainerLogical,"scatChamberContainerPhys",motherLogical,0,0);
	scatChamberContainerLogical->SetVisAttributes(HallVisAtt); 

#ifdef G4DEBUG_GEOMETRY
	if(G4DEBUG_GEOMETRY==1)
	{	
		//debug only
		new G4PVPlacement(pRotX90Z270deg,G4ThreeVector(),scatChamberContainerLogical,
			"scatChamberContainerPhys_Tran",motherLogical,0,0);
	}
	else if(G4DEBUG_GEOMETRY==2)
	{	
		new G4PVPlacement(pRotX90deg,G4ThreeVector(),scatChamberContainerLogical,
			"scatChamberContainerPhys_Long",motherLogical,0,0);
	}
	else if(G4DEBUG_GEOMETRY==3)
	{	
		new G4PVPlacement(pRotX90Z353deg,G4ThreeVector(),scatChamberContainerLogical,
			"scatChamberContainerPhys_Gep",motherLogical,0,0);
	}
#endif


	/////////////////////////
	// scatter chamber 
	/////////////////////////
	//Build the scatter chamber,it contains 5 windows,  4 rectangles and 1 circle
	//target chamber container || scattering chamber container
	//double mScatChamberRin=17.875*inch,mScatChamberRout=18.875*inch,mScatChamberL=27.25*inch;

	G4VSolid* scatChamberWholeSolid=0;
	//If mSetupG2PScatChamber==1, setup the body only, 
	//If mSetupG2PScatChamber==2, setup the body plus top flange and bottom flange, this
	//will make the program slower
	//The bottom flange is like the following
	//                       II                         II
	//                       II                         II
	//                       II                         II     __start here: -mScatChamberL/2
	//                      /II                         II\                                      0
	//                     /III                         III\                 -mScatChamberL/2-1.0"
	//                     IIII                         IIII                
	//                     IIII                         IIII               
	//                     IIII                         IIII                 -mScatChamberL/2-3.25"
	//                     IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
	//                     IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII     __end here: -mScatChamberL/2-4.5"
	//the top flange is similar in shape but upside down
	if(mSetupG2PScatChamber==1)
	{
		scatChamberWholeSolid = new G4Tubs("scatChamberWholeTubs",
			mScatChamberRin,mScatChamberRout,mScatChamberL/2.0,0.,360.*deg);
	}
	else if(mSetupG2PScatChamber>=2)
	{
		startphi=0.*deg; deltaphi=360.*deg;
		const int kNPlane_SC=11;
		double rInner_SC[] = {0,0,mScatChamberRin,
			mScatChamberRin,mScatChamberRin,mScatChamberRin,
			mScatChamberRin,mScatChamberRin,mScatChamberRin,
			0,0};
		double rOuter_SC[] = {
			mScatChamberRout+1.0*inch,mScatChamberRout+1.0*inch,mScatChamberRout+1.0*inch,
			mScatChamberRout,mScatChamberRout,mScatChamberRout+1.0*inch,
			mScatChamberRout+1.0*inch,mScatChamberRout,mScatChamberRout,
			mScatChamberRout+1*inch,mScatChamberRout+1*inch
		};
		double zPlane_SC[] = {
			-mScatChamberL/2-4.50*inch,-mScatChamberL/2-3.25*inch,-mScatChamberL/2-1.0*inch,
			-mScatChamberL/2,mScatChamberL/2+0.25*inch,mScatChamberL/2+1.25*inch,
			mScatChamberL/2+3.50*inch,mScatChamberL/2+3.50*inch,mScatChamberL/2+20.5*inch,
			mScatChamberL/2+20.5*inch,mScatChamberL/2+21.75*inch
		};

		G4Polycone* SCWholeSolid = new G4Polycone("SCPolycone",startphi,deltaphi,
			kNPlane_SC,zPlane_SC,rInner_SC,rOuter_SC);

		scatChamberWholeSolid = SCWholeSolid;
	}
	else if(mSetupG2PScatChamber>=30)
	{
		//build the body, then union with top and bottom flanges
		//this is the same as above
		startphi=0.*deg; deltaphi=360.*deg;
		G4VSolid* SCBodySolid = new G4Tubs("SCBodyTubs",
			mScatChamberRin,mScatChamberRout,mScatChamberL/2.0,startphi,deltaphi);

		const int kNPlane_SCBottom=4;
		double rInner_SCBottom[] = {0,0,mScatChamberRin,mScatChamberRin};
		double rOuter_SCBottom[] = {mScatChamberRout+1.0*inch,mScatChamberRout+1.0*inch,
			mScatChamberRout+1.0*inch,mScatChamberRout};
		double zPlane_SCBottom[] = {-mScatChamberL/2-4.50*inch,-mScatChamberL/2-3.25*inch,
			-mScatChamberL/2-1.0*inch,-mScatChamberL/2};
		G4Polycone* SCBottomFlangeSolid = new G4Polycone("SCBottomFlangePolycone",startphi,deltaphi,
			kNPlane_SCBottom,zPlane_SCBottom,rInner_SCBottom,rOuter_SCBottom);

		const int kNPlane_SCTop=8;
		double rInner_SCTop[] = {
			mScatChamberRin,mScatChamberRin,mScatChamberRin,
			mScatChamberRin,mScatChamberRin,mScatChamberRin,
			0,0
		};
		double rOuter_SCTop[] = {
			mScatChamberRout,mScatChamberRout,mScatChamberRout+1.0*inch,
			mScatChamberRout+1.0*inch,mScatChamberRout,mScatChamberRout,
			mScatChamberRout+1*inch,mScatChamberRout+1*inch 
		};
		double zPlane_SCTop[] = { 
			mScatChamberL/2,mScatChamberL/2+0.25*inch,mScatChamberL/2+1.25*inch,
			mScatChamberL/2+3.50*inch,mScatChamberL/2+3.50*inch,mScatChamberL/2+20.5*inch,
			mScatChamberL/2+20.5*inch,mScatChamberL/2+21.75*inch
		};

		G4Polycone* SCTopFlangeSolid = new G4Polycone("SCTopFlangePolycone",startphi,deltaphi,
			kNPlane_SCTop,zPlane_SCTop,rInner_SCTop,rOuter_SCTop);

		G4UnionSolid* SCBodyUnionBSolid = new G4UnionSolid("SCBodyUnionBSolid",
			SCBodySolid,SCBottomFlangeSolid);
		G4UnionSolid* SCBodyUnionBTSolid = new G4UnionSolid("SCBodyUnionBTSolid",
			SCBodyUnionBSolid,SCTopFlangeSolid);

		scatChamberWholeSolid = SCBodyUnionBTSolid;
	}


	//these are for the subtraction part, not the scatter chamber itself
	double pSCExitWindowH=15.0*inch;
	double pSCEntranceWindowH=6.44*inch;
	double pSCEntranceWindowH_GEP=6.44*inch;
	double pSCEntranceWindowVoffset=-1.766*inch;
	double pSCWindowRin=mShieldLN2Rin-1.0*cm;
	double pSCWindowRout=mScatChamberRout+1.0*cm;

	//In a bollean solid, if SolidA=SolidB-SolidC, it require solidC thicker than solidB 
	//the following angles are listed in the drawing
	//EntranceWindowTran (phi=72+36) for transverse, rec
	//EntranceWindowLong (d=5.75',phi=180) for logitidinal,  circle
	//EntranceWindowGEP (phi=156.45+7.1) for 20 degree GEP, rec
	//ExitWindowTran (phi=252+36), ExitWindowLong(phi=309+60), ExitWindowLeftSide (phi=19+32), 
	//Viewwindow1(d=3.75'), viewwindow2(d=3.75') are not included 

	//the following angles come from the angles in the drawing subtract 90 degrees, this is 
	//the G4 Lab convention
	//EntranceWindowTran (phi=342+36) for transverse, reduce to 7 degree for g2p 
	//EntranceWindowLong (d=5.75',phi=90) for logitidinal,  circle
	//EntranceWindowGEP (phi=66.5+7) for 20 degree GEP, rec
	//ExitWindowTran (phi=162+36), ExitWindowLong(phi=219+60), ExitWindowLeftSide (phi=289+32), 
	//Viewwindow1(d=3.75'), viewwindow2(d=3.75') are not included 

	//rectangle ExitWindowTran
	startphi=162.0*deg;deltaphi=36*deg;
	G4VSolid* SCExitWindowTranSolid = new G4Tubs("SCExitWindowTranTubs",
		pSCWindowRin,pSCWindowRout,pSCExitWindowH/2.0,startphi,deltaphi);

	//rectangle ExitWindowLong
	startphi=219.0*deg;deltaphi=60*deg;
	G4VSolid* SCExitWindowLongSolid = new G4Tubs("SCExitWindowLongTubs",
		pSCWindowRin,pSCWindowRout,pSCExitWindowH/2.0,startphi,deltaphi);

	//rectangle ExitWindowLeftSide
	startphi=289.0*deg;deltaphi=32*deg;
	G4VSolid* SCExitWindowLeftSideSolid = new G4Tubs("SCExitWindowLeftSideTubs",
		pSCWindowRin,pSCWindowRout,pSCExitWindowH/2.0,startphi,deltaphi);


	//////////////////////////////////////////////////////////////////////////
	//rectangle EntranceWindowTran
	//rectangle EntranceWindowTran, for SANE transverse, LN2 shielding will subtract this solid
	startphi=342.0*deg;deltaphi=36*deg;
	G4VSolid* SCEntranceWindowTranSANESolid = new G4Tubs("SCEntranceWindowTranSANETubs",
		pSCWindowRin,pSCWindowRout,pSCExitWindowH/2.0,startphi,deltaphi);

	//rectangle EntranceWindowTran, for G2P transverse
	//In G2P, Al patch this original window and dig a small window on the patch to use it as
	//the transverse entrance. The size of this entrance window is identical to the GEP entrace
	//in this program, I dig this small window from the scattering chamber
	//But I still dig the big window in LN2 shieding since Al did not patch it
	int pSetupHMS;
	gConfig->GetParameter("SetupHMS",pSetupHMS);    //in case this is for SANE or HALL-C WACS
	if(pSetupHMS)  {startphi=342.0*deg;deltaphi=36*deg;pSCEntranceWindowH=pSCExitWindowH;}
	else {startphi=356.25*deg;deltaphi=7.5*deg;}	
	G4VSolid* SCEntranceWindowTranSolid = new G4Tubs("SCEntranceWindowTranTubs",
		pSCWindowRin,pSCWindowRout,pSCEntranceWindowH/2.0,startphi,deltaphi);

	//////////////////////////////////////////////////////////////////////////

	//circle, located at phi=180 degree, EntranceWindowLong subtracted part
	double pSCEntranceWindowLongR=5.75*inch/2.0;
	double pSCEntranceWindowLongL_Half=pSCWindowRout-pSCWindowRin;
	double pSCEntranceWindowLongSubY=(pSCWindowRout+pSCWindowRin)/2.0;
	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* SCEntranceWindowLongSubSolid = new G4Tubs("SCEntranceWindowLongSubTubs",
		0.0,pSCEntranceWindowLongR,pSCEntranceWindowLongL_Half,startphi,deltaphi);

	//circle, located at phi=180 degree, EntranceWindowLong unioned part, 
	//the flange, outer diameter is 8 inch, opening is 5.6 inch,
	double pSCEntranceWindowLongFlangeRout=8.0*inch/2;
	double pSCEntranceWindowLongFlangeRin=5.6*inch/2;
	//double pSCEntranceWindowLongFlangeL=40.0*mm;  //the thickness of the flange, has been delared 
	//this flange attached to the openning, started position is
	//sqrt(SCRout^2 - Ropen^2) 
	double pSCEntranceWindowLongFlangeY = sqrt(pow(mScatChamberRout,2)-pow(pSCEntranceWindowLongR,2)) + 
		pSCEntranceWindowLongFlangeL/2.0;
	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* SCEntranceWindowLongFlangeSolid = new G4Tubs("SCEntranceWindowLongFlangeTubs",
		pSCEntranceWindowLongFlangeRin,pSCEntranceWindowLongFlangeRout,
		pSCEntranceWindowLongFlangeL/2.0,startphi,deltaphi);

	//rectangle EntranceWindowGEP 
	startphi=66.25*deg;deltaphi=7.5*deg;
	G4VSolid* SCEntranceWindowGEPSolid = new G4Tubs("SCEntranceWindowGEPTubs",
		pSCWindowRin,pSCWindowRout,pSCEntranceWindowH_GEP/2.0,startphi,deltaphi);

	//////////////////////////////////////////////////////////////////////
	//build the boolean geometry

	// subtract the transverse exit 
	G4SubtractionSolid* SCSubtractExitTSolid=new G4SubtractionSolid(
		"SCSubtractExitT",scatChamberWholeSolid,SCExitWindowTranSolid);

	// subtract the longitudinal exit 
	G4SubtractionSolid* SCSubtractExitTNLSolid=new G4SubtractionSolid(
		"SCSubtractExitTNL",SCSubtractExitTSolid,SCExitWindowLongSolid);

	// subtract the left side exit, the orignal size 
	G4SubtractionSolid* SCSubtractExitTNLNSSolid=new G4SubtractionSolid(
		"SCSubtractExitTNLNS",SCSubtractExitTNLSolid,SCExitWindowLeftSideSolid);

	// subtract the transverse entrance 
	// G2p entrance window has a vertical offset but SANE and WACS in Hall C will not  
	G4ThreeVector pV3EntranceWindowTran(0,0,pSCEntranceWindowVoffset);
	if(pSetupHMS) pV3EntranceWindowTran.set(0,0,0);

	G4SubtractionSolid* SCSubtractExitNEntranceTSolid=new G4SubtractionSolid(
		"SCSubtractExitNEntranceT",SCSubtractExitTNLNSSolid,SCEntranceWindowTranSolid,
		0,pV3EntranceWindowTran);

	// subtract the logitudinal entrance, which is a circle 
	//need to rotate about X (in container Coor, in Hall Coor it is Y) by 90 degrees then 
	//subtract it to make a circle hole
	G4SubtractionSolid* SCSubtractExitNEntranceTNLSolid=new G4SubtractionSolid(
		"SCSubtractExitNEntranceTNL",SCSubtractExitNEntranceTSolid,
		SCEntranceWindowLongSubSolid,pRotX90deg,G4ThreeVector(0,pSCEntranceWindowLongSubY,0));

	//subtract GEP entrance windows, 1.766 inch below
	G4SubtractionSolid* SCSubtractExitNEntranceTNLNGSolid=new G4SubtractionSolid(
		"SCSubtractExitNEntranceTNLNG",SCSubtractExitNEntranceTNLSolid,
		SCEntranceWindowGEPSolid,0,G4ThreeVector(0,0,pSCEntranceWindowVoffset));

	//union with the flange(the protrude part) of logitudinal entrance
	G4UnionSolid *SCSubtractExitNEntranceTNLNGUnionLSolid=new G4UnionSolid(
		"SCSubtractExitNEntranceTNLNGUnionL",SCSubtractExitNEntranceTNLNGSolid,
		SCEntranceWindowLongFlangeSolid,pRotX90deg,G4ThreeVector(0,pSCEntranceWindowLongFlangeY,0));

	//setup the scatter chamber
	G4LogicalVolume* scatChamberLogical = new G4LogicalVolume(
		SCSubtractExitNEntranceTNLNGUnionLSolid,mMaterialManager->aluminum,"scatChamberLogical",0,0,0);
	scatChamberLogical->SetVisAttributes(WhiteVisAtt); 

	new G4PVPlacement(0,G4ThreeVector(),scatChamberLogical,"scatChamberPhys",
		scatChamberContainerLogical,0,0);

	/////////////////////////
	// target chamber window covers 
	/////////////////////////
	//Covers for EntranceWindow Tran,Long,Side ,ExitWindowTran, ExitWindowLong,  
	//ExitWindowLeftSide, viewwindow1(d=3.75'), viewwindow2(d=3.75') are not included 
	/////////////////////////////////////////////////////////////////////
	//Al told that the exit window for 0, 20 and 90 degrees target field are all 20 mil of aluminum
	//double mScatChamberExitWindowThick=0.02*inch;
	double pSCExitWindowCoverH=pSCExitWindowH+0.8*inch;
	double pSCEntranceWindowCoverH=pSCEntranceWindowH+0.8*inch;
	double pSCExitWindowCoverRin=mScatChamberRout;
	double pSCExitWindowCoverRout=mScatChamberRout+mScatChamberExitWindowThick;

	//ExitWindowTranCover
	startphi=162.0*deg;deltaphi=36*deg;
	G4VSolid* SCExitWindowTranCoverSolid = new G4Tubs("SCExitWindowTranTubs",
		pSCExitWindowCoverRin,pSCExitWindowCoverRout,pSCExitWindowCoverH/2.0,startphi,deltaphi);
	G4LogicalVolume* SCExitWindowTranCoverLogical = new G4LogicalVolume(
		SCExitWindowTranCoverSolid,mMaterialManager->aluminum,"SCExitWindowTranCoverLogical",0,0,0);
	new G4PVPlacement(0,G4ThreeVector(),SCExitWindowTranCoverLogical,
		"SCExitWindowTranCoverPhys",scatChamberContainerLogical,false,0);
	SCExitWindowTranCoverLogical->SetVisAttributes(LightYellowVisAtt); 

	//ExitWindowLongCover
	startphi=219.0*deg;deltaphi=60*deg;
	G4VSolid* SCExitWindowLongCoverSolid = new G4Tubs("SCExitWindowLongTubs",
		pSCExitWindowCoverRin,pSCExitWindowCoverRout,pSCExitWindowCoverH/2.0,startphi,deltaphi);
	G4LogicalVolume* SCExitWindowLongCoverLogical = new G4LogicalVolume(
		SCExitWindowLongCoverSolid,mMaterialManager->aluminum,"SCExitWindowLongCoverLogical",0,0,0);
	SCExitWindowLongCoverLogical->SetVisAttributes(LightYellowVisAtt); 

	new G4PVPlacement(0,G4ThreeVector(),SCExitWindowLongCoverLogical,
		"SCExitWindowLongCoverPhys",scatChamberContainerLogical,false,0);

	//ExitWindowLeftSideCover
	startphi=289.0*deg;deltaphi=32*deg;
	G4VSolid* SCExitWindowLeftSideCoverSolid = new G4Tubs("SCExitWindowLeftSideTubs",
		pSCExitWindowCoverRin,pSCExitWindowCoverRout,pSCExitWindowCoverH/2.0,startphi,deltaphi);
	G4LogicalVolume* SCExitWindowLeftSideCoverLogical = new G4LogicalVolume(
		SCExitWindowLeftSideCoverSolid,mMaterialManager->aluminum,
		"SCExitWindowLeftSideCoverLogical",0,0,0);
	SCExitWindowLeftSideCoverLogical->SetVisAttributes(LightYellowVisAtt); 

	new G4PVPlacement(0,G4ThreeVector(),SCExitWindowLeftSideCoverLogical,
		"SCExitWindowLeftSideCoverPhys",scatChamberContainerLogical,false,0);

	/////////////////////////////////////////////////////////////////////
	//Al told that the entrance window for 20 and 90 degrees target field are 7 mil of aluminum
	double pSCEntranceWindowCoverRin=mScatChamberRout;
	double pSCEntranceWindowCoverRout=mScatChamberRout+0.007*inch;

	//EntranceWindowTranCover, G2P window
	if(pSetupHMS)  {startphi=342.0*deg;deltaphi=36*deg;}  //in case this is for SANE or HALL-C WACS
	else {startphi=356.25*deg;deltaphi=7.5*deg;}	//G2P
	G4VSolid* SCEntranceWindowTranCoverSolid = new G4Tubs("SCEntranceWindowTranCoverTubs",
		pSCEntranceWindowCoverRin,pSCEntranceWindowCoverRout,
		pSCEntranceWindowCoverH/2.0,startphi,deltaphi);

	G4LogicalVolume* SCEntranceWindowTranCoverLogical = new G4LogicalVolume(
		SCEntranceWindowTranCoverSolid,mMaterialManager->aluminum,
		"SCEntranceWindowTranCoverLogical",0,0,0);
	SCEntranceWindowTranCoverLogical->SetVisAttributes(LightYellowVisAtt); 

	G4ThreeVector pV3EntranceWindowTranCover(0,0,pSCEntranceWindowVoffset);
	if(pSetupHMS) pV3EntranceWindowTranCover.set(0,0,0);
	new G4PVPlacement(0,pV3EntranceWindowTranCover,
		SCEntranceWindowTranCoverLogical,"SCEntranceWindowTranCoverPhys",
		scatChamberContainerLogical,false,0);


	////EntranceWindowLongCover, this rotation is very strange, but it works!!!
	//!!!Note, only work if the BField rotation is about Y axis
	////circle, located at 180 phi=degree, this entrance window is for 0 degree target field
	//Al told me that the material is 20 mil Beryllium
	double pBeThick=0.02*inch;
	double pSCEntranceWindowLongCoverY=pSCEntranceWindowLongFlangeY
		+pSCEntranceWindowLongFlangeL/2+pBeThick/2.0;
	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* SCEntranceWindowLongCoverSolid = new G4Tubs("SCEntranceWindowLongCoverTubs",
		0.0,pSCEntranceWindowLongFlangeRout,pBeThick/2.0,startphi,deltaphi);
	G4LogicalVolume* SCEntranceWindowLongCoverLogical = new G4LogicalVolume(
		SCEntranceWindowLongCoverSolid,mMaterialManager->beryllium,
		"SCEntranceWindowLongCoverLogical",0,0,0);
	SCEntranceWindowLongCoverLogical->SetVisAttributes(LightYellowVisAtt);

	G4RotationMatrix* pRotEntranceWindowLongCover=new G4RotationMatrix();
	*pRotEntranceWindowLongCover=(*pRotBField)*(*pRotBField)*(*pRotX90deg);
	/////////////////////////////////////////////////////////
#ifdef G4DEBUG_GEOMETRY
	if(G4DEBUG_GEOMETRY>=2)
	{
		cout<<"\n (*pRotBField)*(*pRotBField)*(*pRotX90deg) ==> EulerAngles: phi="
			<<pRotEntranceWindowLongCover->getPhi()/deg
			<<"  theta="<<pRotEntranceWindowLongCover->getTheta()/deg
			<<"  psi="<<pRotEntranceWindowLongCover->getPsi()/deg<<endl;

		G4RotationMatrix* pRotEn2Cover=new G4RotationMatrix();
		*pRotEn2Cover=(*pRotX270deg)*(*pRotBField)*(*pRotBField);
		cout<<"\n (*pRotX270deg)*(*pRotBField)*(*pRotBField) ==> EulerAngles: phi="
			<<pRotEn2Cover->getPhi()/deg
			<<"  theta="<<pRotEn2Cover->getTheta()/deg
			<<"  psi="<<pRotEn2Cover->getPsi()/deg<<endl;

		*pRotEn2Cover=(*pRotY90deg)*(*pRotBField);
		cout<<"\n (*pRotY90deg)*(*pRotBField) ==> EulerAngles: phi="<<pRotEn2Cover->getPhi()/deg
			<<"  theta="<<pRotEn2Cover->getTheta()/deg
			<<"  psi="<<pRotEn2Cover->getPsi()/deg<<endl;
	}
#endif
	/////////////////////////////////////////////////////////
	new G4PVPlacement(pRotEntranceWindowLongCover,G4ThreeVector(0,pSCEntranceWindowLongCoverY,0),
		SCEntranceWindowLongCoverLogical,"SCEntranceWindowLongCoverPhys",
		scatChamberContainerLogical,false,0); 

	//EntranceWindowGEPCover, for GEP
	startphi=66.25*deg;deltaphi=7.5*deg;
	G4VSolid* SCEntranceWindowGEPCoverSolid = new G4Tubs("SCEntranceWindowGEPCoverTubs",
		pSCEntranceWindowCoverRin,pSCEntranceWindowCoverRout,
		pSCEntranceWindowCoverH/2.0,startphi,deltaphi);
	G4LogicalVolume* SCEntranceWindowGEPCoverLogical = new G4LogicalVolume(
		SCEntranceWindowGEPCoverSolid,mMaterialManager->aluminum,
		"SCEntranceWindowGEPCoverLogical",0,0,0);
	SCEntranceWindowGEPCoverLogical->SetVisAttributes(LightYellowVisAtt); 

	new G4PVPlacement(0,G4ThreeVector(0,0,pSCEntranceWindowVoffset),SCEntranceWindowGEPCoverLogical,
		"SCEntranceWindowGEPCoverPhys",scatChamberContainerLogical,false,0);
	/////////////////////////////////////////////////////////////////////

	/////////////////////////
	// innerOfChamber 
	/////////////////////////
	//inside the scatter chamber, it is vacuum, all stuff inside the chamber will use 
	//the local coordination of this mother volumn, 
	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* innerOfChamberSolid = new G4Tubs("innerOfChamberTubs",mShieldLHeRin,mScatChamberRin,
		mScatChamberL/2.0,startphi,deltaphi);
	G4LogicalVolume* innerOfChamberLogical = new G4LogicalVolume(innerOfChamberSolid,
		mMaterialManager->vacuum,"innerOfChamberLogical",0,0,0);
	//By Jixie: Add this step limit can help in calculating the integrated BdL
	double pScatChamberStepLimit=10;
	gConfig->GetArgument("ScatChamberStepLimit",pScatChamberStepLimit);
	pScatChamberStepLimit*=mm;
	G4UserLimits* uSCStepLimits = new G4UserLimits(pScatChamberStepLimit);
	innerOfChamberLogical->SetUserLimits(uSCStepLimits);
	innerOfChamberLogical->SetVisAttributes(HallVisAtt);  //white invisible

	new G4PVPlacement(0,G4ThreeVector(),innerOfChamberLogical,"innerOfChamberPhys",
		scatChamberContainerLogical,0,0);


	/////////////////////////
	// LN2 shielding 
	/////////////////////////

	//LN2 Shield cylinder itself is 1/4 inch, but have windows similar to the scattering chamber
	//there are also covers for each opening. In order to simplify the geometry, I union a
	//thin aluminum foil (1.5 mil) onto the outer surface of it as the cover
	//LN2 Shield cynlider with Din=32.5", Dout=33" and L=33.25", the cover is 0.0015'=0.0381mm,
	//Drawing# 67504-E-0015 sheet 2
	//double mShieldLN2Rin=16.15*inch,mShieldLN2Rout=16.50*inch, mShieldLN2L=33.25*inch;
	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* shieldLN2WholeSolid = new G4Tubs("shieldLN2WholeTubs",mShieldLN2Rin,mShieldLN2Rout,
		mShieldLN2L/2.0,startphi,deltaphi);

	if(mShieldLN2WindowThick/mm < 1.0E-05) mShieldLN2WindowThick=1.0E-06*mm;
	G4VSolid* shieldLN2CoverSolid = new G4Tubs("shieldLN2CoverTubs",mShieldLN2Rout,
		mShieldLN2Rout+mShieldLN2WindowThick,mShieldLN2L/2.0,startphi,deltaphi);

	//////////////////////////////////////////////////////////////////////
	//build the boolean geometry

	// subtract the transverse exit 
	G4SubtractionSolid* ShieldLN2SubtractExitTSolid=new G4SubtractionSolid(
		"ShieldLN2SubtractExitT",shieldLN2WholeSolid,SCExitWindowTranSolid);

	// subtract the longitudinal exit 
	G4SubtractionSolid* ShieldLN2SubtractExitTNLSolid=new G4SubtractionSolid(
		"ShieldLN2SubtractExitTNL",ShieldLN2SubtractExitTSolid,SCExitWindowLongSolid);

	// subtract the left side exit 
	G4SubtractionSolid* ShieldLN2SubtractExitTNLNSSolid=new G4SubtractionSolid(
		"ShieldLN2SubtractExitTNLNS",ShieldLN2SubtractExitTNLSolid,SCExitWindowLeftSideSolid);

	// subtract the SANE transverse entrance 
	G4SubtractionSolid* ShieldLN2SubtractExitNEntranceTSolid=new G4SubtractionSolid(
		"ShieldLN2SubtractExitNEntranceT",
		ShieldLN2SubtractExitTNLNSSolid,SCEntranceWindowTranSANESolid);

	// subtract the logitudinal entrance, which is a circle 
	//need to rotate about X (in container Coor, in Hall Coor it is Y) by 90 degrees then 
	//subtract it to make a circle hole
	G4SubtractionSolid* ShieldLN2SubtractExitNEntranceTNLSolid=new G4SubtractionSolid(
		"ShieldLN2SubtractExitNEntranceTNL",ShieldLN2SubtractExitNEntranceTSolid,
		SCEntranceWindowLongSubSolid,pRotX90deg,G4ThreeVector(0,pSCEntranceWindowLongSubY,0));

	//subtract GEP entrance windows, 1.766 inch below
	// Al make it a circle, since GEP angle is not finalize yet, I do not want to change it now
	// need to change it later
	G4SubtractionSolid* ShieldLN2SubtractExitNEntranceTNLNGSolid=new G4SubtractionSolid(
		"ShieldLN2SubtractExitNEntranceTNLNG",ShieldLN2SubtractExitNEntranceTNLSolid,
		SCEntranceWindowGEPSolid,0,G4ThreeVector(0,0,pSCEntranceWindowVoffset));

	//Union with the cover 
	G4VSolid* ShieldLN2Solid=(G4VSolid*)ShieldLN2SubtractExitNEntranceTNLNGSolid;
	if(mShieldLN2WindowThick/mm > 1.0E-05)
	{
		ShieldLN2Solid=new G4UnionSolid("shieldLN2Soild",
			ShieldLN2SubtractExitNEntranceTNLNGSolid,shieldLN2CoverSolid);
	}
	/////////////////////////////////////////////////////////////////////
	//now build and place the geometry

	G4LogicalVolume* shieldLN2Logical = new G4LogicalVolume(ShieldLN2Solid,
		mMaterialManager->aluminum,"shieldLN2Logical",0,0,0);
	shieldLN2Logical->SetVisAttributes(WhiteVisAtt);  //white

	new G4PVPlacement(0,G4ThreeVector(),shieldLN2Logical,"shieldLN2Phys",
		innerOfChamberLogical,0,0);


	/////////////////////////
	// target field coils 
	/////////////////////////
	//
	if (mSetupCoil==1)
	{
		// just 2 helm coils, downstream and upstream, connection part are not included
		//Downtream and upstream coil
		//this number here is for Bz along Z axis
		startphi=0.*deg; deltaphi=360.*deg;
		const int kNPlane_Coil=6;
		const double kInnerLineTheta=39.65*deg; //The inner surface is alone a line with this angle
		double rInner_Coil[] = {4.04*inch,4.04*inch,3.82*inch/tan(kInnerLineTheta),
			11.60*inch,12.06*inch,12.06*inch};
		double rOuter_Coil[] = {5.14*inch,11.60*inch,12.56*inch,12.56*inch,12.56*inch,12.56*inch};
		double zPlane_DownCoil[] = {1.57*inch,3.37*inch,3.82*inch,9.61*inch,9.611*inch,10.21*inch};
		double zPlane_UpCoil[kNPlane_Coil];
		for(int ii=0;ii<kNPlane_Coil;ii++) zPlane_UpCoil[ii]=-1*zPlane_DownCoil[ii];

		G4VSolid* downCoilSolid = new G4Polycone("DownCoilPcon",startphi,deltaphi,
			kNPlane_Coil,zPlane_DownCoil,rInner_Coil,rOuter_Coil);
		G4VSolid* upCoilSolid = new G4Polycone("UpCoilPcon",startphi,deltaphi,
			kNPlane_Coil,zPlane_UpCoil,rInner_Coil,rOuter_Coil);
		G4LogicalVolume* downCoilLogical = new G4LogicalVolume(downCoilSolid,
			mMaterialManager->stainlesssteel,"downCoilLogical",0,0,0);
		downCoilLogical->SetVisAttributes(DarkBlueVisAtt);
		G4LogicalVolume* upCoilLogical = new G4LogicalVolume(upCoilSolid,
			mMaterialManager->stainlesssteel,"upCoilLogical",0,0,0);
		upCoilLogical->SetVisAttributes(DarkBlueVisAtt);  //dark blue

		new G4PVPlacement(pRotX270deg,G4ThreeVector(),downCoilLogical,
			"DownCoilPhys",innerOfChamberLogical,false,0);
		new G4PVPlacement(pRotX270deg,G4ThreeVector(),upCoilLogical,
			"UpCoilPhys",innerOfChamberLogical,false,0);

	}
	else if (mSetupCoil>=2)
	{
		//Hall B Coil
		startphi=0.*deg;deltaphi=360.*deg;
		G4RotationMatrix* pRotDCoil = new G4RotationMatrix();
		pRotDCoil->rotateY(90.*deg);
		pRotDCoil->rotateX(-90.*deg);
		G4RotationMatrix* pRotUCoil = new G4RotationMatrix();
		pRotUCoil->rotateY(90.*deg);
		pRotUCoil->rotateX(90.*deg);

		//Coils
		const int kNTubs=2;
		double rInnerTubs[] = {116.5*mm,116.5*mm};
		double rOuterTubs[] = {153.0*mm,153.0*mm};
		double zPlaneTubs[] = {-16.0*mm,16.0*mm};
		G4Polycone* InnerCoilSolid = new G4Polycone("InnerCoilSolid",startphi,deltaphi,
			kNTubs,zPlaneTubs,rInnerTubs,rOuterTubs);
		//By Jixie: The tub is not shown correctly in the visulization, Chao use a polycone solid instead
		//This is a bug of HEPREP viewer itself
		//G4Tubs* InnerCoilSolid = new G4Tubs("InnerCoilSolid",116.5*mm,153.0*mm,32.0*mm/2.0,startphi,deltaphi);
		G4LogicalVolume* InnerCoilLogical = new G4LogicalVolume(InnerCoilSolid,
			mMaterialManager->copper,"InnerCoilLogical",0,0,0);
		InnerCoilLogical->SetVisAttributes(CuBrownVisAtt);
		new G4PVPlacement(pRotDCoil,G4ThreeVector(0,-68.0*mm,0),InnerCoilLogical,
			"InnerCoilPhys",innerOfChamberLogical,true,0,0);
		new G4PVPlacement(pRotUCoil,G4ThreeVector(0,68.0*mm,0),InnerCoilLogical,
			"InnerCoilPhys",innerOfChamberLogical,true,1,0);
		rInnerTubs[0]=rInnerTubs[1]=165.5*mm;
		rOuterTubs[0]=rOuterTubs[1]=222.5*mm;
		zPlaneTubs[0]=-24.0*mm;zPlaneTubs[1]=24.0*mm;
		G4Polycone* OuterCoilSolid = new G4Polycone("OuterCoilSolid",startphi,deltaphi,
			kNTubs,zPlaneTubs,rInnerTubs,rOuterTubs);
		//By Jixie: The tub is not shown correctly in the visulization, Chao use a polycone solid instead        
		//G4Tubs* OuterCoilSolid = new G4Tubs("OuterCoilSolid",165.5*mm,222.5*mm,48.0*mm/2.0,startphi,deltaphi);
		G4LogicalVolume* OuterCoilLogical = new G4LogicalVolume(OuterCoilSolid,
			mMaterialManager->copper,"OuterCoilLogical",0,0,0);
		OuterCoilLogical->SetVisAttributes(CuBrownVisAtt);
		new G4PVPlacement(pRotDCoil,G4ThreeVector(0,-99.0*mm,0),OuterCoilLogical,
			"OuterCoilPhys",innerOfChamberLogical,true,2,0);
		new G4PVPlacement(pRotUCoil,G4ThreeVector(0,99.0*mm,0),OuterCoilLogical,
			"OuterCoilPhys",innerOfChamberLogical,true,3,0);

		//Coil Former
		const int kNPlane_Coil_F=9;
		double rInner_Coil_F[] = {105.0*mm,105.0*mm,153.0*mm,
			153.0*mm,153.0*mm,153.0*mm,153.0*mm,222.5*mm,222.5*mm};
		double rOuter_Coil_F[] = {130.0*mm,130.0+34.5/tan(17.0*deg)*mm,130.0+34.5/tan(17.0*deg)*mm,
			250.0*mm,250.0*mm,244.5*mm,244.5*mm,244.5*mm,244.5*mm};
		double zPlane_Coil_F[] = {0.0*mm,34.5*mm,34.5*mm,
			120*tan(17.0*deg)*mm,62.0*mm,62.0*mm,74.0*mm,74.0*mm,76.0*mm};
		G4Polycone* CoilFPolycone = new G4Polycone("CoilFPcon",startphi,deltaphi,
			kNPlane_Coil_F,zPlane_Coil_F,rInner_Coil_F,rOuter_Coil_F);
		G4SubtractionSolid* CoilFSub1 = new G4SubtractionSolid("CoilFSub1",
			CoilFPolycone,InnerCoilSolid,0,G4ThreeVector(0,0,27.5*mm));
		G4SubtractionSolid* CoilFSolid = new G4SubtractionSolid("CoilFSolid",
			CoilFSub1,OuterCoilSolid,0,G4ThreeVector(0,0,58.5*mm));
		G4LogicalVolume* CoilFLogical = new G4LogicalVolume(CoilFSolid,
			mMaterialManager->stainlesssteel,"CoilFLogical",0,0,0);
		CoilFLogical->SetVisAttributes(DarkBlueVisAtt);
		new G4PVPlacement(pRotDCoil,G4ThreeVector(0,-40.5*mm,0),
			CoilFLogical,"CoilFormerPhys",innerOfChamberLogical,true,0,0);
		new G4PVPlacement(pRotUCoil,G4ThreeVector(0,40.5*mm,0),
			CoilFLogical,"CoilFormerPhys",innerOfChamberLogical,true,1,0);

		//Coil Retainer
		const int kNPlane_Coil_R=7;
		double rInner_Coil_R[] = {105.0*mm,105.0*mm,105.0*mm,105.0*mm,
			(150.0-9.0*tan(50.0*deg))*mm,(150.0-9.0*tan(50.0*deg))*mm,150.0*mm};
		double rOuter_Coil_R[] = {116.5*mm,116.5*mm,153.0*mm,153.0*mm,
			153.0*mm,165.5*mm,165.5*mm};
		double zPlane_Coil_R[] = {0.0*mm,9.0*mm,9.0*mm,(48.0-45.0/tan(50.0*deg))*mm,
			39.0*mm,39.0*mm,48.0*mm};
		G4Polycone* CoilRSolid = new G4Polycone("CoilRSolid",startphi,deltaphi,
			kNPlane_Coil_R,zPlane_Coil_R,rInner_Coil_R,rOuter_Coil_R);
		G4LogicalVolume* CoilRLogical = new G4LogicalVolume(CoilRSolid,
			mMaterialManager->stainlesssteel,"CoilRLogical",0,0,0);
		CoilRLogical->SetVisAttributes(DarkBlueVisAtt);
		new G4PVPlacement(pRotDCoil,G4ThreeVector(0,-75.0*mm,0),
			CoilRLogical,"CoilRetainerPhys",innerOfChamberLogical,true,0,0);
		new G4PVPlacement(pRotUCoil,G4ThreeVector(0,75.0*mm,0),
			CoilRLogical,"CoilRetainerPhys",innerOfChamberLogical,true,1,0);

		//Joint Flange
		const int kNPlane_Joint_F=5;
		double rInner_Joint_F[] = {170.0*mm,170.0*mm,170.0*mm,170.0*mm,177.0*mm};
		double rOuter_Joint_F[] = {221.0*mm,221.0*mm,240.0*mm,240.0*mm,240.0*mm};
		double zPlane_Joint_F[] = {0.0*mm,6.0*mm,6.0*mm,15.0*mm,22.0*mm};
		G4Polycone* JointFSolid = new G4Polycone("JointFSolid",startphi,deltaphi,
			kNPlane_Joint_F,zPlane_Joint_F,rInner_Joint_F,rOuter_Joint_F);
		G4LogicalVolume* JointFLogical = new G4LogicalVolume(JointFSolid,
			mMaterialManager->plastic,"ClampPLogical",0,0,0);
		JointFLogical->SetVisAttributes(YellowVisAtt);
		new G4PVPlacement(pRotDCoil,G4ThreeVector(0,-123.0*mm,0),
			JointFLogical,"JointFlangePhys",innerOfChamberLogical,true,0,0);
		new G4PVPlacement(pRotUCoil,G4ThreeVector(0,123.0*mm,0),
			JointFLogical,"JointFlangePhys",innerOfChamberLogical,true,1,0);

		//Left Side Surface
		const int kNPlane_Spinning=3;
		double rInner_Spinning[] = {102.0*mm,102.0*mm,241.5*mm};
		double rOuter_Spinning[] = {105.0*mm,105.0*mm,244.5*mm};
		double zPlane_Spinning[] = {40.5*mm,(123.0-45.0/tan(50.0*deg))*mm,205.7*mm};
		G4Polycone* SpinningSolid = new G4Polycone("SpinningSolid",startphi,deltaphi,
			kNPlane_Spinning,zPlane_Spinning,rInner_Spinning,rOuter_Spinning);
		rInnerTubs[0]=rInnerTubs[1]=244.5*mm;
		rOuterTubs[0]=rOuterTubs[1]=250.0*mm;
		zPlaneTubs[0]=-103.2/2.0*mm;zPlaneTubs[1]=103.2/2.0*mm;
		G4Polycone* LCTubSolid = new G4Polycone("LCTubSolid",startphi,deltaphi,
			kNTubs,zPlaneTubs,rInnerTubs,rOuterTubs);
		//G4Tubs* LCTubSolid = new G4Tubs("LCTubSolid",244.5*mm,250.0*mm,103.2/2.0*mm,startphi,deltaphi);
		G4LogicalVolume* SpinningLogical = new G4LogicalVolume(SpinningSolid,
			mMaterialManager->stainlesssteel,"SpinningLogical",0,0,0);
		G4LogicalVolume* LCTubLogical = new G4LogicalVolume(LCTubSolid,
			mMaterialManager->stainlesssteel,"LeftCoverTubLogical",0,0,0);
		SpinningLogical->SetVisAttributes(DarkBlueVisAtt);
		LCTubLogical->SetVisAttributes(DarkBlueVisAtt);
		new G4PVPlacement(pRotDCoil,G4ThreeVector(),
			SpinningLogical,"LeftSurfPhys",innerOfChamberLogical,true,0,0);
		new G4PVPlacement(pRotDCoil,G4ThreeVector(0,-154.1*mm,0),
			LCTubLogical,"LeftSurfPhys",innerOfChamberLogical,true,1,0);

		//Right Side Surface
		G4Tubs* FridgeTubsOuter = new G4Tubs("FridgeTubsOuter",0.0*mm,22.25*mm,500.0/2.0*mm,
			startphi,deltaphi);
		G4RotationMatrix* pRotTub = new G4RotationMatrix();
		pRotTub->rotateY(-25.*deg);
		G4Tubs* HeBaseFlgTubs1 = new G4Tubs("HeBaseFlgTubs1",38.1*mm,250.0*mm,19.5/2.0*mm,
			startphi,deltaphi);
		G4Tubs* HeBaseFlgTubs2 = new G4Tubs("HeBaseFlgTubs2",244.5*mm,250.0*mm,99.4/2.0*mm,
			startphi,deltaphi);
		G4UnionSolid* HeBaseFlgU1 = new G4UnionSolid("HeBaseFlgU1",
			HeBaseFlgTubs1,HeBaseFlgTubs2,0,G4ThreeVector(0,0,-39.95*mm));
		G4SubtractionSolid* HeBaseFlgSolid = new G4SubtractionSolid("HeBaseFlgSolid",
			HeBaseFlgU1,FridgeTubsOuter,
			pRotTub,G4ThreeVector((192.15+37.0/sin(25.0*deg))*tan(25.0*deg)*mm,0,0));
		G4LogicalVolume *HeBaseFlgLogical = new G4LogicalVolume(HeBaseFlgSolid,
			mMaterialManager->stainlesssteel,"HeBaseFlgLogical",0,0,0);
		HeBaseFlgLogical->SetVisAttributes(DarkBlueVisAtt);
		G4Tubs* HeCanBulkheadTubs1 = new G4Tubs("HeCanBulkheadTubs1",38.1*mm,105.0*mm,9.5/2.0*mm,
			startphi,deltaphi);
		G4Tubs* HeCanBulkheadTubs2 = new G4Tubs("HeCanBulkheadTubs2",102.0*mm,105.0*mm,24.4/2.0*mm,
			startphi,deltaphi);
		G4UnionSolid* HeCanBulkheadU1 = new G4UnionSolid("HeCBU1",
			HeCanBulkheadTubs1,HeCanBulkheadTubs2,0,G4ThreeVector(0,0,-7.45*mm));
		G4SubtractionSolid* HeCBSolid = new G4SubtractionSolid("HeCBSolid",
			HeCanBulkheadU1,FridgeTubsOuter,
			pRotTub,G4ThreeVector((60.15+37.0/sin(25.0*deg))*tan(25.0*deg)*mm,0,0));
		G4LogicalVolume* HeCBLogical = new G4LogicalVolume(HeCBSolid,
			mMaterialManager->stainlesssteel,"HeCBLogical",0,0,0);
		HeCBLogical->SetVisAttributes(DarkBlueVisAtt);
		//By Jixie: The tub is not shown correctly in the HEPREP visulization, 
		// Chao use an intersection solid to replace it, this is a bug of HEPREP viewer
		//G4Tubs* BeamTubSolid = new G4Tubs("BeamTubSolid",36.5*mm,36.6*mm,146.5/2.0*mm,startphi,deltaphi);
		G4Tubs* BeamTubs = new G4Tubs("BeamTubs",36.5*mm,36.6*mm,300.0/2.0*mm,startphi,deltaphi);
		G4Box* BeamBox = new G4Box("BeamBox",200.0/2.0*mm,200.0/2.0*mm,146.5/2.0*mm);
		G4IntersectionSolid* BeamTubSolid = new G4IntersectionSolid("BeamTubSolid",BeamBox,BeamTubs);
		G4LogicalVolume* BeamTubLogical = new G4LogicalVolume(BeamTubSolid,
			mMaterialManager->stainlesssteel,"BeamTubLogical",0,0,0);
		BeamTubLogical->SetVisAttributes(DarkBlueVisAtt);
		G4Tubs* FridgeTubs = new G4Tubs("FridgeTubs",21.85*mm,22.25*mm,300.0/2.0*mm,startphi,deltaphi);
		G4Box* FridgeBox = new G4Box("FridgeBox",200.0/2.0*mm,200.0/2.0*mm,146.5/2.0*mm);
		G4IntersectionSolid* FridgeTubSolid = new G4IntersectionSolid("FridgeTubSolid",
			FridgeBox,FridgeTubs,pRotTub,G4ThreeVector());
		G4LogicalVolume* FridgeTubLogical = new G4LogicalVolume(FridgeTubSolid,
			mMaterialManager->stainlesssteel,"FridgeTubLogical",0,0,0);
		FridgeTubLogical->SetVisAttributes(DarkBlueVisAtt);
		new G4PVPlacement(pRotUCoil,G4ThreeVector(0,192.15*mm,0),
			HeBaseFlgLogical,"RightSurfPhys",innerOfChamberLogical,true,0,0);
		new G4PVPlacement(pRotUCoil,G4ThreeVector(0,60.15*mm,0),
			HeCBLogical,"RightSurfPhys",innerOfChamberLogical,true,1,0);
		new G4PVPlacement(pRotUCoil,G4ThreeVector(0,128.65*mm,0),
			BeamTubLogical,"RightSurfPhys",innerOfChamberLogical,true,2,0);
		new G4PVPlacement(pRotUCoil,
			G4ThreeVector(0,128.65*mm,(128.65+37.0/sin(25.0*deg))*tan(25.0*deg)*mm),
			FridgeTubLogical,"RightSurfPhys",innerOfChamberLogical,true,3,0);

		//Left Side LHe Container
		const int kNPlane_LeftHe=4;
		double rInner_LeftHe[] = {222.5*mm,222.5*mm,150.0*mm,244.5*mm};
		double rOuter_LeftHe[] = {244.5*mm,244.5*mm,244.5*mm,244.5*mm};
		double zPlane_LeftHe[] = {116.5*mm,123.0*mm,123.0*mm,205.7*mm};
		double zPlane_LeftHe_S[] = {-1.0*mm,6.0*mm,6.0*mm,15.0*mm,22.0*mm};
		G4Polycone* LeftHeSPolycone = new G4Polycone("LeftHeSPcon",startphi,deltaphi,
			kNPlane_Joint_F,zPlane_LeftHe_S,rInner_Joint_F,rOuter_Joint_F);
		G4Polycone* LeftHePolycone = new G4Polycone("LeftHePcon",startphi,deltaphi,
			kNPlane_LeftHe,zPlane_LeftHe,rInner_LeftHe,rOuter_LeftHe);
		G4SubtractionSolid* LeftHeSolid = new G4SubtractionSolid("LeftHeSolid",
			LeftHePolycone,LeftHeSPolycone,0,G4ThreeVector(0,0,123.0*mm));
		G4LogicalVolume* LeftHeLogical = new G4LogicalVolume(LeftHeSolid,
			mMaterialManager->liquidHe,"LeftHeLogical",0,0,0);
		LeftHeLogical->SetVisAttributes(WhiteVisAtt);
		new G4PVPlacement(pRotDCoil,G4ThreeVector(),
			LeftHeLogical,"LeftHePhys",innerOfChamberLogical,false,0);

		//Right Side LHe Container
		const int kNPlane_RightHe=5;
		double rInner_RightHe[] = {38.1*mm,38.1*mm,38.1*mm,38.1*mm,38.1*mm};
		double rOuter_RightHe[] = {105.0*mm,105.0*mm,150.0*mm,244.5*mm,244.5*mm};
		double zPlane_RightHe[] = {64.9*mm,(123.0-45.0/tan(50.0*deg))*mm,123.0*mm,123.0*mm,182.4*mm};
		G4Polycone* RightHePolycone = new G4Polycone("RightHePcon",startphi,deltaphi,
			kNPlane_RightHe,zPlane_RightHe,rInner_RightHe,rOuter_RightHe);
		G4SubtractionSolid* RightHeS1 = new G4SubtractionSolid("RightHeS1",
			RightHePolycone,FridgeTubsOuter,pRotTub,G4ThreeVector(37.0/cos(25.0*deg)*mm,0,0));
		G4SubtractionSolid* RightHeS2 = new G4SubtractionSolid("RightHeS2",
			RightHeS1,LeftHeSPolycone,0,G4ThreeVector(0,0,123.0*mm));
		G4Tubs* RightHeTubs = new G4Tubs("RightHeTubs",222.5*mm,244.5*mm,(6.5+1.0)/2.0*mm,
			startphi,deltaphi);
		G4UnionSolid* RightHeSolid = new G4UnionSolid("RightHeSolid",
			RightHeS2,RightHeTubs,0,G4ThreeVector(0,0,120.25));
		G4LogicalVolume* RightHeLogical = new G4LogicalVolume(RightHeSolid,
			mMaterialManager->liquidHe,"RightHeLogical",0,0,0);
		RightHeLogical->SetVisAttributes(WhiteVisAtt);
		new G4PVPlacement(pRotUCoil,G4ThreeVector(),
			RightHeLogical,"RightHePhys",innerOfChamberLogical,false,0);

		//Large Wedge
		const int kNPlane_LWedge=10;
		double rInner_LWedge[] = {250.0*mm,(130.0+9.5/tan(17.0*deg))*mm,(130.0+9.5/tan(17.0*deg))*mm,
			130.0*mm,112.0*mm,112.0*mm,130.0*mm,(130.0+9.5/tan(17.0*deg))*mm,
			(130.0+9.5/tan(17.0*deg))*mm,250.0*mm};
		double rOuter_LWedge[] = {250.0*mm,250.0*mm,220.0*mm,
			220.0*mm,220.0*mm,220.0*mm,220.0*mm,220.0*mm,
			250.0*mm,250.0*mm};
		double zPlane_LWedge[] = {-(40.5+120*tan(17.0*deg))*mm,-50.0*mm,-50.0*mm,
			-40.5*mm,-40.5*mm,40.5*mm,40.5*mm,50.0*mm,
			50.0*mm,(40.5+120*tan(17.0*deg))*mm};
		G4Polycone* LWedgePolycone = new G4Polycone("LWedgePcon",-20.0*deg,40.0*deg,
			kNPlane_LWedge,zPlane_LWedge,rInner_LWedge,rOuter_LWedge);
		G4Tubs* LWedgeTubs = new G4Tubs("LWedgeTubs",0.0*mm,38.2*mm,1000/2.0*mm,
			startphi,deltaphi);//The orininal radius of this hole is 31.5 mm, change it to 38.2 mm
		G4Tubs* LWedgeTubs2 = new G4Tubs("LWedgeTubs2",0.0*mm,130.0*mm,40.0/2.0*mm,startphi,deltaphi);
		G4RotationMatrix* pRotWedgeTubs = new G4RotationMatrix();
		pRotWedgeTubs->rotateY(90.*deg);
		G4SubtractionSolid* LWedgeS1 = new G4SubtractionSolid("LWedgeS1",
			LWedgePolycone,LWedgeTubs,pRotWedgeTubs,G4ThreeVector());
		G4SubtractionSolid* LWedgeSolid = new G4SubtractionSolid("LWedgeSolid",LWedgeS1,LWedgeTubs2);
		G4LogicalVolume* LWedgeLogical = new G4LogicalVolume(LWedgeSolid,
			mMaterialManager->aluminum,"LWedgeLogical",0,0,0);
		LWedgeLogical->SetVisAttributes(SilverVisAtt);
		G4RotationMatrix* pRotWedge1 = new G4RotationMatrix();
		pRotWedge1->rotateY(90.*deg);
		pRotWedge1->rotateX(90.*deg);
		new G4PVPlacement(pRotWedge1,G4ThreeVector(),
			LWedgeLogical,"LargeWedgePhys",innerOfChamberLogical,true,0,0);
		G4RotationMatrix* pRotWedge2 = new G4RotationMatrix();
		pRotWedge2->rotateY(-90.*deg);
		pRotWedge2->rotateX(90.*deg);
		new G4PVPlacement(pRotWedge2,G4ThreeVector(),
			LWedgeLogical,"LargeWedgePhys",innerOfChamberLogical,true,1,0);

		//Small Wedge
		const int kNPlane_SWedge=10;
		double rInner_SWedge[] = {250.0*mm,130.0*mm,112.0*mm,112.0*mm,
			130.0*mm,130.0*mm,112.0*mm,112.0*mm,130.0*mm,250.0*mm};
		double rOuter_SWedge[] = {250.0*mm,250.0*mm,250.0*mm,250.0*mm,
			250.0*mm,250.0*mm,250.0*mm,250.0*mm,250.0*mm,250.0*mm};
		double zPlane_SWedge[] = {-(40.5+120*tan(17.0*deg))*mm,-40.5*mm,-40.5*mm,-20.0*mm,
			-20.0*mm,20.0*mm,20.0*mm,40.5*mm,40.5*mm,(40.5+120*tan(17.0*deg))*mm};
		G4Polycone* SWedgeSolid = new G4Polycone("SWedgeSolid",-7.5/2.0*deg,7.5*deg,
			kNPlane_SWedge,zPlane_SWedge,rInner_SWedge,rOuter_SWedge);
		G4LogicalVolume* SWedgeLogical= new G4LogicalVolume(SWedgeSolid,
			mMaterialManager->stainlesssteel,"SWedgeLogical",0,0,0);
		SWedgeLogical->SetVisAttributes(GrayVisAtt);
		G4RotationMatrix* pRotWedge3 = new G4RotationMatrix();
		pRotWedge3->rotateX(90.*deg);
		pRotWedge3->rotateZ(30.*deg);
		new G4PVPlacement(pRotWedge3,G4ThreeVector(),
			SWedgeLogical,"SmallWedgePhys",innerOfChamberLogical,true,0,0);
		G4RotationMatrix* pRotWedge4 = new G4RotationMatrix();
		pRotWedge4->rotateX(90.*deg);
		pRotWedge4->rotateZ(150.*deg);
		new G4PVPlacement(pRotWedge4,G4ThreeVector(),
			SWedgeLogical,"SmallWedgePhys",innerOfChamberLogical,true,1,0);
		G4RotationMatrix* pRotWedge5 = new G4RotationMatrix();
		pRotWedge5->rotateX(90.*deg);
		pRotWedge5->rotateZ(-30.*deg);
		new G4PVPlacement(pRotWedge5,G4ThreeVector(),
			SWedgeLogical,"SmallWedgePhys",innerOfChamberLogical,true,2,0);
		G4RotationMatrix* pRotWedge6 = new G4RotationMatrix();
		pRotWedge6->rotateX(90.*deg);
		pRotWedge6->rotateZ(-150.*deg);
		new G4PVPlacement(pRotWedge6,G4ThreeVector(),
			SWedgeLogical,"SmallWedgePhys",innerOfChamberLogical,true,3,0);

	}

	/////////////////////////
	// LHe Shielding (4k shielding)
	/////////////////////////
	//LHe Shield cylider with r=38.1mm and L=700mm, thickness=0.0015'=0.0381mm
	//double mShieldLHeRin=38.1*mm,mShieldLHeRout=38.1381*mm, mShieldLHeL=700.0*mm;
	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* shieldLHeSolid = new G4Tubs("shieldLHeTubs",mShieldLHeRin,mShieldLHeRout,
		mShieldLHeL/2.0,startphi,deltaphi);
	G4LogicalVolume* shieldLHeLogical = new G4LogicalVolume(shieldLHeSolid,
		mMaterialManager->aluminum,"shieldLHeLogical",0,0,0);
	shieldLHeLogical->SetVisAttributes(WhiteVisAtt);  //white

	new G4PVPlacement(0,G4ThreeVector(),shieldLHeLogical,"shieldLHePhys",innerOfChamberLogical,0,0);

	return scatChamberContainerPhys;
}


/////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* WACSDetectorConstruction::ConstructWACSTarget(G4LogicalVolume* motherLogical)
{
	G4String SDname;

	G4VSensitiveDetector* targetSD=new HRSDCSD(SDname="targetMaterial");
	G4VSensitiveDetector* targetNoseSD=new HRSDCSD(SDname="targetNose");
	G4VSensitiveDetector* targetEndCapSD=new HRSDCSD(SDname="targetEndCap");
	G4VSensitiveDetector* targetWallSD=new HRSDCSD(SDname="targetWall");

	G4SDManager* SDman = G4SDManager::GetSDMpointer();

	const double inch=2.54*cm;
	double startphi,deltaphi;

	G4RotationMatrix* pRotX90deg = new G4RotationMatrix();
	pRotX90deg->rotateX(90.*deg);
	G4RotationMatrix* pRotX270deg = new G4RotationMatrix();
	pRotX270deg->rotateX(270.*deg);

	//////////////////////////////////
	//target container
	//////////////////////////////////

	//The target containner should stay in the center of Scatter chamber, otherwise
	//there will be some overlapping
	//Note that the target offsets are in the hall coordinate system, I have to shift it 
	//to the target containner system(scat chamber system) 
	//after pRotX90deg, (x,y,z)_ScatChamber==>(x,-z,y)_Hall
	double pTargetOriginX=mTargetXOffset-mScatChamberXOffset;
	double pTargetOriginY=-(mTargetZOffset-mScatChamberZOffset);
	double pTargetOriginZ=mTargetYOffset-mScatChamberYOffset;

	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* targetContainerSolid = new G4Tubs("targetContainerTubs",0,mShieldLHeRin,
		mShieldLHeL/2.0,startphi,deltaphi);

	G4LogicalVolume* targetContainerLogical = new G4LogicalVolume(targetContainerSolid,
		mMaterialManager->vacuum,"targetContainerLogical",0,0,0);
	targetContainerLogical->SetVisAttributes(HallVisAtt);  //white invisible

	G4VPhysicalVolume* targetContainerPhys=new G4PVPlacement(pRotX90deg,
		G4ThreeVector(mScatChamberXOffset,mScatChamberYOffset,mScatChamberZOffset),
		targetContainerLogical,"targetContainerPhys",motherLogical,0,0);


	//////////////////////////////////
	//tail|nose shiled
	//////////////////////////////////
	//By Jixie: told by Josh that the thickness is 4 mil
	//Toby get the number from Josh: Inner diameter = 1.654", Outter diameter=1.662"

	double pShieldNoseRin=1.654/2*inch,pShieldNoseRout=pShieldNoseRin+0.004*inch; 
	double pShieldNoseL=mShieldLHeL-10.0*mm;
	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* shieldNoseSolid = new G4Tubs("shieldNoseTubs",pShieldNoseRin,pShieldNoseRout,
		pShieldNoseL/2.0,startphi,deltaphi);

	G4LogicalVolume* shieldNoseLogical = new G4LogicalVolume(shieldNoseSolid,
		mMaterialManager->aluminum,"shieldNoseLogical",0,0,0);
	shieldNoseLogical->SetVisAttributes(WhiteVisAtt);  //white

	//since I place the target inside the nose, I will put the center of the nose at the
	//target position, therefore only the nose should use the target offset
	new G4PVPlacement(0,G4ThreeVector(pTargetOriginX,pTargetOriginY,pTargetOriginZ),
		shieldNoseLogical,"shieldNosePhys",targetContainerLogical,0,0);


	//////////////////////////////////
	//inside the tail|nose 
	//////////////////////////////////

	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* innerOfNoseSolid = new G4Tubs("innerOfNoseTubs",0,pShieldNoseRin,
		pShieldNoseL/2.0,startphi,deltaphi);
	G4Material* theCoolant=mMaterialManager->liquidHe;
	//for optics, we do not put liquid helium inside the tail|nose as high as normal runs
	//I will use heliumGas4k instead
	if (mSetupG2PTarget>=10) theCoolant=mMaterialManager->heliumGas4k;

	G4LogicalVolume* innerOfNoseLogical = new G4LogicalVolume(innerOfNoseSolid,
		theCoolant,"innerOfNoseLogical",0,0,0);
	innerOfNoseLogical->SetVisAttributes(LightPurpleVisAtt);  //light purple
	
	SDman->AddNewDetector(targetNoseSD);
	innerOfNoseLogical->SetSensitiveDetector(targetNoseSD);

	//By Jixie: Add this step limit can help to calculate the integrated BdL
	double pTargetStepLimit=10;
	gConfig->GetArgument("TargetStepLimit",pTargetStepLimit);
	pTargetStepLimit*=mm;
	G4UserLimits* uTNStepLimits = new G4UserLimits(pTargetStepLimit);
	innerOfNoseLogical->SetUserLimits(uTNStepLimits);

	//since I place the target inside the nose, I will put the center of the nose at the
	//target position, therefore only the nose should use the target offset
	new G4PVPlacement(0,G4ThreeVector(pTargetOriginX,pTargetOriginY,pTargetOriginZ),
		innerOfNoseLogical,"innerOfNosePhys",targetContainerLogical,0,0);

	//////////////////////////////////
	//target insert, target wall and target cap
	//////////////////////////////////

	//#SetupG2PTarget is used to determine whether to setup target wall and also determine the
	//#the radius of the target cell and the target position w.r.t the target nose. 
	//#It can be 0, 1, 2, 3, 6, 10, 20, 30, 60 or 99
	//#0 means do not set up the target, 
	//#1 means it is the 1st target cell(R_out=1.142"/2 and R_in=TargetR=1.072"/2, centered at zero), 
	//#  for NH3He target. Can have target wall. The 4th and 5th target cell are identical to the 1st cell. 
	//#2 is the middle size hole(D=0.397", centered at -10.95mm) for C12. No target wall 
	//#3 is the smallest hole(D=0.312", centered at -10.95mm) for CH2. No target wall 
	//#6 is the 6th target cell(R_out=1.142"/2 and R_in=TargetR=0.960"/2, with the front face aligned to
	//#  the upstream end, z_face=-1.113"/2=-14.135mm) for C12. With target wall, No end caps. 
	//#10, 20, 30, and 60 are the the same as 1,2,3 and 6,respectively, but NO LHe inside the target nose.
	//#99 means use a pvc pipe and a 10-mil C12 target foil located at z=-896.0 mm, which is the situation 
	//#   of commissioning. In this setup, no target wall, no target cell, TargetType, TargetR, TargetL above 
	//#   will be ignored 

	//The target wall and end caps need to be rotated in its mother volumn 
	//in target nose coordinate system, 
	//x'=x_lab, y'=-z_lab, z'=y_lab

	//The 2nd and 3rd hole are smaller, their radii will be set later according to mSetupG2PTarget
	//R_2nd=0.397"/2   R_3rd=0.312"/2
	//Based on the drawing, the hole diameter of the insert is 1.098" 
	//But the outer diameter of the target cell is 1.142"
	//In order to simplify the geometry and avoid overlappings, I use 1.142" here
	double pTargetInsertHoleR=1.142*inch/2; 
	if(mSetupG2PTarget==2 || mSetupG2PTarget==20)  pTargetInsertHoleR=0.397*inch; 
	else if(mSetupG2PTarget==3 || mSetupG2PTarget==30) pTargetInsertHoleR=0.312*inch; 

	//////////////////////////////////////////////////////////////////////////
	//Note that if mSetupG2PTarget==2|3 then no need to put target wall and caps
	
	//#Setup the target wall, will setup target wall only if SetupG2PTarget==1,6 or 10,60
	int pSetupTargetWall=1; 
	if(mSetupG2PTarget!=1 && mSetupG2PTarget!=10 && mSetupG2PTarget!=6 && mSetupG2PTarget!=60) 
	{
		pSetupTargetWall=0;
	}
	int pSetupEndCap=1;
	if(mSetupG2PTarget!=1 && mSetupG2PTarget!=10) pSetupEndCap=0;

	
	//#The target wall is PCTFE, 1st, 3rd and 5th cells are 35 mil or 0.899 mm thick
	//#But 6th cell (for C12) is 90 mil
	//TargetWallThick=0.889;  #2.286 for cell 6 or 60 for C12
	double pTargetWallL=1.113*inch;
	double pTargetWallThick=0.035*inch;
	if(mSetupG2PTarget==6 || mSetupG2PTarget==60) pTargetWallThick=0.090*inch;

	double pTargetR=pTargetInsertHoleR;	
	if(pSetupTargetWall)  pTargetR=pTargetInsertHoleR-pTargetWallThick;
	//////////////////////////////////////////////////////////////////////////


	//////////////////////////////////
	//the insert, just build one block, not 6 blocks
	/////////////////////////////////
	G4VSolid* targetInsertHoleSolid = new G4Tubs("targetInsertHoleTubs",
		0,pTargetInsertHoleR, 0.13*inch/2,0.,360.*deg);

	double pTargetInsertX_Lab=1.156*inch; //in lab system
	double pTargetInsertY_Lab=1.230*inch; //in lab system
	double pTargetInsertZ_Lab=0.125*inch; //in lab system 
	G4VSolid* targetInsertWholeSolid = new G4Box("targetInsertWholeBox",
		pTargetInsertX_Lab/2,pTargetInsertY_Lab/2,pTargetInsertZ_Lab/2);

	//now dig the hole
	G4SubtractionSolid*	targetInsertSolid=new G4SubtractionSolid("targetInsertSolid",
		targetInsertWholeSolid,targetInsertHoleSolid);

	G4LogicalVolume* targetInsertLogical = new G4LogicalVolume(targetInsertSolid,
		mMaterialManager->aluminum, "targetInsertLogical",0,0,0);
	targetInsertLogical->SetVisAttributes(GrayVisAtt);

	//Insert position w.r.t nose, should be positive before 270 deg rotation
	//based on the drawing, it is 1.113/2-0.188+0.125/2 = 0.4310"=10.95mm
	double pTargetInsertInNoseCS=0.4310*inch;  //in nose system, 
											   
	new G4PVPlacement(pRotX270deg,
		G4ThreeVector(0.0,pTargetInsertInNoseCS,0.0),
		targetInsertLogical,"targetInsertPhys",innerOfNoseLogical,0,0);


	//////////////////////////////////////////////////////////////////////////

	if(pSetupTargetWall)
	{
		////////////////////////////////////////
		//target cell holder
		////////////////////////////////////////
		double pTargetCellHolderZ_Lab=0.10*inch; //in lab system 
		G4VSolid* targetCellHolderWholeSolid = new G4Box("targetCellHolderWholeBox",
			pTargetInsertX_Lab/2,pTargetInsertY_Lab/2,pTargetCellHolderZ_Lab/2);

		G4SubtractionSolid*	targetCellHolderSolid=new G4SubtractionSolid("targetCellHolderSolid",
			targetCellHolderWholeSolid,targetInsertHoleSolid);

		G4LogicalVolume* targetCellHolderLogical = new G4LogicalVolume(targetCellHolderSolid,
			mMaterialManager->PCTFE,"targetCellHolderLogical",0,0,0);
		targetCellHolderLogical->SetVisAttributes(LightYellowVisAtt);  

		double pTargetCellHolderInNoseCS=pTargetInsertInNoseCS-
			pTargetInsertZ_Lab/2-pTargetCellHolderZ_Lab/2;
		new G4PVPlacement(pRotX270deg,
			G4ThreeVector(0.0,pTargetCellHolderInNoseCS,0.0),
			targetCellHolderLogical,"targetCellHolderPhys",innerOfNoseLogical,0,0);

		////////////////////////////////////////
		//target wall
		////////////////////////////////////////
		G4VSolid* targetWallSolid = new G4Tubs("targetWallTubs",pTargetR,
			pTargetInsertHoleR,pTargetWallL/2.0,0.,360.*deg);

		G4LogicalVolume* targetWallLogical = new G4LogicalVolume(targetWallSolid,
			mMaterialManager->PCTFE,"targetWallLogical",0,0,0);
		targetWallLogical->SetVisAttributes(LightYellowVisAtt);
	
		SDman->AddNewDetector(targetWallSD);
		targetWallLogical->SetSensitiveDetector(targetWallSD);


		new G4PVPlacement(pRotX270deg,G4ThreeVector(),
			targetWallLogical,"targetWallPhys",innerOfNoseLogical,0,0);
	}

	////////////////////////////////////////
	//target entrance cap and exit cap
	////////////////////////////////////////
	//The end cap is 0.7 mil aluminum foil, told by Josh 
	double pTargetCapThick=0.0007*inch;
	if(pSetupEndCap) 
	{
		//double mTargetCapThick=0.0007*inch;
		double pTargetCapZ=pTargetWallL/2.0+pTargetCapThick/2.0;
		G4VSolid* targetCapSolid = new G4Tubs("targetCapTubs",0,pTargetR,
			pTargetCapThick/2.0,0.,360.*deg);

		G4LogicalVolume* targetCapLogical = new G4LogicalVolume(targetCapSolid,
			mMaterialManager->aluminum,"targetCapLogical",0,0,0);
		targetCapLogical->SetVisAttributes(GrayVisAtt);  //Grey invisable

		SDman->AddNewDetector(targetEndCapSD);
		targetCapLogical->SetSensitiveDetector(targetEndCapSD);

		//Entrance cap, 
		new G4PVPlacement(pRotX270deg,
			G4ThreeVector(0.0,pTargetCapZ,0.0),
			targetCapLogical,"targetEntranceCapPhys",innerOfNoseLogical,true,0);
		//Exit cap
		new G4PVPlacement(pRotX270deg,
			G4ThreeVector(0.0,-pTargetCapZ,0.0),
			targetCapLogical,"targetExitCapPhys",innerOfNoseLogical,true,1);
	}
	

	//////////////////////////////////
	//the target,
	//////////////////////////////////

	//major target is 55% solid NH3 mixed with LHe cylinder with r=13.6144mm and z=28.27mm 
	//double pTargetR=13.6144*mm, mTargetL=28.27*mm;
	G4VSolid* targetSolid = new G4Tubs("targetTubs",0.,pTargetR,mTargetL/2.0,0.,360.*deg);

	//for optics, we do not put liquid helium inside the tail|nose as high as normal run
	//I will use vacuum instead 
	G4Material* theTargetMaterial=WACS_NH3He;
	if (mTargetType==0) theTargetMaterial=mMaterialManager->vacuum;
	else if (mTargetType==1) theTargetMaterial=this->WACS_NH3He;
	else if (mTargetType==2) theTargetMaterial=mMaterialManager->CH2;
	else if (mTargetType==3) theTargetMaterial=mMaterialManager->carbon;
	else if (mTargetType==4) theTargetMaterial=mMaterialManager->tantalum;
	else if (mTargetType==5) theTargetMaterial=mMaterialManager->liquidH2;
	else if (mTargetType==6) theTargetMaterial=mMaterialManager->liquidD2;
	else if (mTargetType==7) theTargetMaterial=mMaterialManager->liquidHe3;
	else if (mTargetType==8) theTargetMaterial=mMaterialManager->liquidHe;
	else if (mTargetType==9) theTargetMaterial=mMaterialManager->aluminum;
	else if (mTargetType==10) theTargetMaterial=mMaterialManager->copper;
	else if (mTargetType==11) theTargetMaterial=mMaterialManager->lead;
	else if (mTargetType==12) theTargetMaterial=mMaterialManager->tungsten;
	else if (mTargetType==13) theTargetMaterial=mMaterialManager->stainlesssteel;
	else if (mTargetType==14) theTargetMaterial=mMaterialManager->kapton;
	else if (mTargetType==15) theTargetMaterial=mMaterialManager->air;       //room temperature
	else if (mTargetType==16) theTargetMaterial=mMaterialManager->heliumGas; //room temperature
	else if (mTargetType==17) theTargetMaterial=mMaterialManager->calcium;  
	else if (mTargetType==18) theTargetMaterial=mMaterialManager->NH3He;

	G4LogicalVolume* targetLogical = new G4LogicalVolume(targetSolid,theTargetMaterial,
		"targetLogical",0,0,0);
	targetLogical->SetVisAttributes(VioletVisAtt); 
	targetLogical->SetUserLimits(uTNStepLimits);

	//By Jixie: to study the radiator, I make the target sensitive
	SDman->AddNewDetector(targetSD);
	targetLogical->SetSensitiveDetector(targetSD);

	//since I place the target inside the nose, I will put the center of the nose at the
	//target position, therefore only the nose should use the target offset

	//this number is in target nose coordinate system, 
	//x'=x, y'=-z, z'=y
	double pTargetInNoseCS=0;
	if(mSetupG2PTarget==2 || mSetupG2PTarget==20 || mSetupG2PTarget==3 || mSetupG2PTarget==30) 
	{
		pTargetInNoseCS=pTargetInsertInNoseCS; //0.431 *inch;
	}
	else if(mSetupG2PTarget==6 || mSetupG2PTarget==60) 
	{
		pTargetInNoseCS=pTargetWallL/2-mTargetL/2; 
	}
	new G4PVPlacement(pRotX270deg,G4ThreeVector(0,pTargetInNoseCS,0),
		targetLogical,"targetPhys",innerOfNoseLogical,0,0);

	/////////////////////////////////////////////////////////////////////

	return targetContainerPhys;
}

/////////////////////////////////////////////////////////////////////
//By Jixie:
//There are 2 ways to place the FZ fields, A) using global field manager + field map
//B) placing the uniform field into a field logical volumn
//If there is overlaps, the field might now work.
//Note that one logical can take only onen field, therefore we have to build 2 identical logical
//if we use method B, while only one logical + 2 placements is needed in mentod A.
//One can update the field during run time in method A, while method B can only change it before run
//In this routine I choose method A
G4VPhysicalVolume* WACSDetectorConstruction::ConstructWACSChicane(G4LogicalVolume* motherLogical)
{
	const double inch=2.54*cm;
	// All G4VSensitiveDetector will be managed (deleted) by SDManager
	//therefore no need to delete it at the end of this subroutine
	G4String SDname;
	G4VSensitiveDetector* FZB2SD=new HRSStdSD(SDname="FZB2VD");
	G4VSensitiveDetector* FZDumpSD=new HRSStdSD(SDname="FZDump");
	G4SDManager* SDman = G4SDManager::GetSDMpointer();


	int pUseDefaultFZB2=1;
	gConfig->GetArgument("UseDefaultFZB2",pUseDefaultFZB2); 

	G4VPhysicalVolume* FZB2FieldContainerPhys=0;
	///////////////////////////////////////////////////////
	//Dipole FZ2 will be placed, and use B=1.1T all the time 

	double pFZB2TiltedAngle=0*deg;
	double pFZB2PosX=   -0.0*mm+mPivotXOffset;
	double pFZB2PosY= -100.0*mm+mPivotYOffset;
	double pFZB2PosZ=-2300.0*mm+mPivotZOffset;
	double pFZB2Bx=1.10,pFZB2By=0,pFZB2Bz=0;

	////////////////////////////////////////////////////////
	if(pUseDefaultFZB2)
	{
		//I have to updated these parameters in the map so they can be written into the config tree
		gConfig->SetArgument("FZB2TiltedAngle",pFZB2TiltedAngle/rad); 
		gConfig->SetArgument("FZB2PosX",pFZB2PosX/mm); 
		gConfig->SetArgument("FZB2PosY",pFZB2PosY/mm); 
		gConfig->SetArgument("FZB2PosZ",pFZB2PosZ/mm); 
		gConfig->SetArgument("FZB2Bx",pFZB2Bx/tesla); 
		gConfig->SetArgument("FZB2By",pFZB2By/tesla); 
		gConfig->SetArgument("FZB2Bz",pFZB2Bz/tesla); 
	}
	else
	{
		//reading the value from input arguments of option -FZB2 or -ChicaneMagnet2
		gConfig->GetArgument("FZB2TiltedAngle",pFZB2TiltedAngle); pFZB2TiltedAngle*=deg;
		gConfig->GetArgument("FZB2PosX",pFZB2PosX); pFZB2PosX*=mm;
		gConfig->GetArgument("FZB2PosY",pFZB2PosY); pFZB2PosY*=mm;
		gConfig->GetArgument("FZB2PosZ",pFZB2PosZ); pFZB2PosZ*=mm;
		gConfig->GetArgument("FZB2Bx",pFZB2Bx); pFZB2Bx*=tesla;
		gConfig->GetArgument("FZB2By",pFZB2By); pFZB2By*=tesla;
		gConfig->GetArgument("FZB2Bz",pFZB2Bz); pFZB2Bz*=tesla;
	}
	////////////////////////////////////////////////////////

	//now build the vectors
	G4ThreeVector pFZB2Pos3V(pFZB2PosX,pFZB2PosY,pFZB2PosZ);
	G4ThreeVector pFZB2Field3V(pFZB2Bx,pFZB2By,pFZB2Bz);


	G4RotationMatrix *pRotFZB2=new G4RotationMatrix();
	pRotFZB2->rotateX(pFZB2TiltedAngle); 

	//built the field for the chicane 
	HRSEMFieldSetup* mEMFieldSetup=HRSEMFieldSetup::GetHRSEMFieldSetup();

	mEMFieldSetup->SetBField3VFZB2(pFZB2Field3V); 
	G4FieldManager* FZB2FieldManager=mEMFieldSetup->GetFieldManagerFZB2();

	//set the step limit
	double pFZBStepLimit=10;
	gConfig->GetArgument("FZBStepLimit",pFZBStepLimit); 
	pFZBStepLimit*=mm;
	G4UserLimits* uFZBStepLimits = new G4UserLimits(pFZBStepLimit);

	///////////////////////////////////////////////////////////
	//#mSetupChicane
	//#setup the bender, 0 means none, 1 for FZB2, 2 for 1-m long Dipole, 3 for AV-magnet, 
	//4: C-type magnet, no coil, just have iron and field containner  

	/////////////////////////
	// FZB magnet Container
	/////////////////////////

	double pFZBX=17.66*inch, pFZBY=15.37*2*inch, pFZBZ=76.34*inch;
	//use half length dipole if (mSetupChicane==2
	if(mSetupChicane==2)  pFZBZ=39.37*inch;   //1.0 meter

	double pFZBContainerX=pFZBX+2*cm, pFZBContainerY=pFZBY+2*cm, pFZBContainerZ=pFZBZ+15*inch;
	G4VSolid* FZBContainerSolid = new G4Box("FZBContainerWholeBox",
		pFZBContainerX/2.0,pFZBContainerY/2.0,pFZBContainerZ/2.0);

	G4LogicalVolume* FZBContainerLogical = new G4LogicalVolume(FZBContainerSolid,
		mMaterialManager->vacuum,"FZBContainerLogical",0,0,0);
	FZBContainerLogical->SetVisAttributes(HallVisAtt);  

	if(mSetupChicane!=0)
	{
		new G4PVPlacement(pRotFZB2,pFZB2Pos3V,
			FZBContainerLogical,"FZB2ContainerPhys",motherLogical,0,0,0);
	}

	/////////////////////////////////////
	//FZB vacuum pipe and flanges
	/////////////////////////////////////

	//stainless steel rectangle pipe, thickness 0.188", there is one flange attached to each end 
	double pFZBVacuumXout=1.65*inch;
	double pFZBVacuumYout=10.0*inch;
	double pFZBVacuumZ=pFZBZ+13.17*inch;   //it is 89.51*inch; I want the variable pFZBZ can be configured
	double pFZBVacuumXin=pFZBVacuumXout-2*0.188*inch;
	double pFZBVacuumYin=pFZBVacuumYout-2*0.188*inch;

	//the outer box
	G4VSolid* FZBVacuumOutSolid = new G4Box("FZBVacuumOutBox",
		pFZBVacuumXout/2.0,pFZBVacuumYout/2.0,pFZBVacuumZ/2.0);

	//the inner box
	G4VSolid* FZBVacuumInSolid = new G4Box("FZBVacuumInBox",
		pFZBVacuumXin/2.0,pFZBVacuumYin/2.0,pFZBVacuumZ/2.0+1.0*mm);

	//the whole vacuum pipe = outer - inner
	G4SubtractionSolid *FZBVacuumPipeSolid = new G4SubtractionSolid("FZBVacuumPipeSolid",
		FZBVacuumOutSolid,FZBVacuumInSolid);

	//place the 2 ends part of the vacuum pipe into the FZ containner
	G4LogicalVolume* FZBVacuumPipeLogical = new G4LogicalVolume(FZBVacuumPipeSolid,
		mMaterialManager->stainlesssteel,"FZBVacuumPipeLogical",0,0,0);
	FZBVacuumPipeLogical->SetVisAttributes(SilverVisAtt);  

	//do not place the vacumn pipe for AV-magnet
	if(mSetupChicane==1 || mSetupChicane==2) {
	  new G4PVPlacement(0,G4ThreeVector(),FZBVacuumPipeLogical,
		"FZBVacuumPipePhys",FZBContainerLogical,false,0,0);
	}

	////////////////////////////////////
	//the end disk (flange)
	////////////////////////////////////
	double pFZBVacuumFlangeR=6.0*inch;
	double pFZBVacuumFlangeZ=0.6*inch;  //TODO: verify this thickness, not important ...
	G4VSolid* FZBVacuumFlangeSolid = new G4Tubs("FZBVacuumFlangeTubs",0,
		pFZBVacuumFlangeR,pFZBVacuumFlangeZ/2,0*deg,360*deg);

	//this flange will be subtracted by a rectange as large as the the vacuum pipe, 
	//then the vacuum pipe will be attached with 2 flanges to make the total length to 90 inch 
	G4SubtractionSolid *FZBVacuumFlangeSubRecSolid = new G4SubtractionSolid(
		"FZBVacuumFlangeSubRecSolid",FZBVacuumFlangeSolid,FZBVacuumOutSolid);



	//place 2 flanges, one at each end into the magnet container
	G4LogicalVolume* FZBVacuumFlangeLogical = new G4LogicalVolume(FZBVacuumFlangeSubRecSolid,
		mMaterialManager->stainlesssteel,"FZBVacuumFlangeLogical",0,0,0);
	FZBVacuumFlangeLogical->SetVisAttributes(GrayVisAtt);  

	//double pFZBVacuumFlangePosZ=90*inch/2-pFZBVacuumFlangeZ/2;
	double pFZBVacuumFlangePosZ=(pFZBZ+13.56*inch)/2-pFZBVacuumFlangeZ/2;

	//do not place the vacumn pipe for AV-magnet
	if(mSetupChicane==1 || mSetupChicane==2) {
	  new G4PVPlacement(0,G4ThreeVector(0,0,pFZBVacuumFlangePosZ),
		FZBVacuumFlangeLogical,"FZBVacuumFlangeDownPhys",FZBContainerLogical,true,0,0);
	  new G4PVPlacement(0,G4ThreeVector(0,0,-pFZBVacuumFlangePosZ),
		FZBVacuumFlangeLogical,"FZBVacuumFlangeUpPhys",FZBContainerLogical,true,1,0);
	}

	////////////////////////////
	//The field containner
	////////////////////////////

	////stainless steel rectangle pipe, thickness 0.188", there is one flange attached to each end	
	//has been declared above
	//double pFZBVacuumXout=1.65*inch;
	//double pFZBVacuumYout=10.0*inch;
	//double pFZBVacuumZ=pFZBZ+13.17*inch;   //it is 89.51*inch; I want the variable pFZBZ can be configured
	//double pFZBVacuumXin=pFZBVacuumXout-2*0.188*inch;
	//double pFZBVacuumYin=pFZBVacuumYout-2*0.188*inch;

	//in real life, the field container is the area in between the mid side plates, 1.66" x 11.5" x 76.34"
	//there will be a stainless steel rectangle pipe, thickness 0.188", with a flange disk attached to each end
	//In order to make life easy,  I use the inner volumn of this stainless steel rectangle pipe to be field containner
	//I will put this field containner inside FZ containner, so I do not have subtract this volumn from FZ containner

	double pFZBFieldContainerX=pFZBVacuumXin;
	double pFZBFieldContainerY=pFZBVacuumYin;
	double pFZBFieldContainerZ=pFZBZ;
	G4VSolid* FZBFieldContainerSolid = new G4Box("FZBFieldContainerBox",
		pFZBFieldContainerX/2.0,pFZBFieldContainerY/2.0,pFZBFieldContainerZ/2.0);

	G4LogicalVolume* FZB2FieldContainerLogical = new G4LogicalVolume(FZBFieldContainerSolid,
		mMaterialManager->vacuum,"FZB2FieldContainerLogical",FZB2FieldManager,0,uFZBStepLimits);
	FZB2FieldContainerLogical->SetVisAttributes(HallVisAtt);  

	if(mSetupChicane!=0)
	{
	  FZB2FieldContainerPhys = new G4PVPlacement(0,G4ThreeVector(0,0,0),
			FZB2FieldContainerLogical,"FZB2FieldContainerPhys",FZBContainerLogical,0,0,0);
	}


	/////////////////////////////////////
	//FZB iron yoke, the silicon steel part
	/////////////////////////////////////
	//The FZ magnet = whole - center plate + 2 middle side plates

	//The yoke of AV-magnet is cut to 1.3 m, but using original coils
	//so I temperarily keep the original length, will set it back when this block finish
	double oldFZBZ=pFZBZ;
	if(mSetupChicane==3)  pFZBZ=1.3*m;

	double pFZBSidePlateX=5.85*inch;
	double pFZBUpDownPlateY=5.98*inch;
	G4VSolid* FZBWholeSolid = new G4Box("FZBWholeBox",pFZBX/2.0,pFZBY/2.0,pFZBZ/2.0);

	//the top and the bottom box, make it longer for subtraction
	double pFZBSubPlateX=pFZBX-2*pFZBSidePlateX;
	double pFZBSubPlateY=pFZBY-2*pFZBUpDownPlateY;
	G4VSolid* FZBSubPlateSolid = new G4Box("FZBSubPlateBox",
		pFZBSubPlateX/2.0,pFZBSubPlateY/2.0,pFZBZ/2.0+1*mm);

	//the middle side plate
	double pFZBPlateMidSideX=pFZBX/2-pFZBVacuumXout/2-pFZBSidePlateX;
	double pFZBPlateMidSideY=11.5*inch;
	G4VSolid* FZBPlateMidSideSolid = new G4Box("FZBPlateMidSideBox",
		pFZBPlateMidSideX/2.0,pFZBPlateMidSideY/2.0,pFZBZ/2.0);

	//subtract the middle rectangle
	G4SubtractionSolid* FZBSubRecSolid=new G4SubtractionSolid("FZBSubRecSolid",
		FZBWholeSolid,FZBSubPlateSolid);
	
	//I was asked to built C-type magnet.....
	//so I subtract the bottom part
	if(mSetupChicane==4)
	{	  
	  G4SubtractionSolid* FZBCTypeSolid=new G4SubtractionSolid("FZBCTypeSolid",
	      FZBSubRecSolid,FZBSubPlateSolid,0,G4ThreeVector(0,-pFZBUpDownPlateY,0));
	  FZBSubRecSolid = FZBCTypeSolid;
	}

	//Union 2 middle side rectangles
	double pFZBPlateMidSidePosX=pFZBVacuumXout/2+pFZBPlateMidSideX/2;
	G4ThreeVector pFZBMidSideLPos3V=G4ThreeVector( pFZBPlateMidSidePosX,0,0);
	G4ThreeVector pFZBMidSideRPos3V=G4ThreeVector(-pFZBPlateMidSidePosX,0,0);
	G4UnionSolid* FZBSubRecUnionMidLSolid=new G4UnionSolid("FZBSubRecUnionMidLSolid",
		FZBSubRecSolid,FZBPlateMidSideSolid,0,pFZBMidSideLPos3V);
	G4UnionSolid* FZBSubRecUnionMidLRSolid=new G4UnionSolid("FZBSubRecUnionMidLRSolid",
		FZBSubRecUnionMidLSolid,FZBPlateMidSideSolid,0,pFZBMidSideRPos3V);

	G4LogicalVolume* FZBSteelLogical = new G4LogicalVolume(FZBSubRecUnionMidLRSolid,
		mMaterialManager->siliconsteel,"FZBSteelLogical",0,0,0);
	FZBSteelLogical->SetVisAttributes(DarkBlueVisAtt);  

	//place one copy of this silicon steel into each magnet
	new G4PVPlacement(0,G4ThreeVector(0,0,(oldFZBZ-pFZBZ)/2),
		FZBSteelLogical,"FZBSteelPhys",FZBContainerLogical,false,0,0);

	//now set the length back
	if(mSetupChicane==3)  pFZBZ=oldFZBZ;
	/////////////////////////////////////
	//FZB copper coils
	/////////////////////////////////////
	//this part is hard to build,  I simplified it as a big box subtract the siliconsteel block
	//made of copper

	//the whole box
	//each group of coil is 0.925" wide, place 2 of them with 1mm gap in between
	double pFZBCoilX=2*0.925*inch+1*mm; 
	double pFZBCoilY=pFZBY-2*5.98*inch; 
	double pFZBCoilZ=pFZBZ+2.62*inch+3.35*2*inch;  	
	G4VSolid* FZBCoilWholeSolid = new G4Box("FZBCoilWholeBox",
		pFZBCoilX/2.0,pFZBCoilY/2.0,pFZBCoilZ/2.0);

	//the small box need to be subtracted
	double pFZBCoilSmallRecX=pFZBCoilX+1*mm;
	double pFZBCoilSmallRecY=pFZBPlateMidSideY+1*mm;
	double pFZBCoilSmallRecZ=pFZBZ+2.62*inch;  
	G4VSolid* FZBCoilSmallRecSolid = new G4Box("FZBCoilSmallRecBox",
		pFZBCoilSmallRecX/2.0,pFZBCoilSmallRecY/2.0,pFZBCoilSmallRecZ/2.0);

	//coil = whole - small_rectangle
	G4SubtractionSolid* FZBCoilSolid=new G4SubtractionSolid("FZBCoilSolid",
		FZBCoilWholeSolid,FZBCoilSmallRecSolid);

	//TODO: cut off the corner ......

	//build the logical volumn
	G4LogicalVolume* FZBCoilLogical = new G4LogicalVolume(FZBCoilSolid,
		mMaterialManager->copper,"FZBCoilLogical",0,0,0);
	FZBCoilLogical->SetVisAttributes(CuBrownVisAtt);  

	//place 2 coils in each FZ magnet
	double pFZBCoilPosX=pFZBX/2-pFZBSidePlateX-pFZBCoilX/2;  //attached to the sife plate
	new G4PVPlacement(0,G4ThreeVector(pFZBCoilPosX,0,0),
		FZBCoilLogical,"FZBCoilLeftPhys",FZBContainerLogical,true,0,0);
	new G4PVPlacement(0,G4ThreeVector(-pFZBCoilPosX,0,0),
		FZBCoilLogical,"FZBCoilRightPhys",FZBContainerLogical,true,1,0);

	/////////////////////////////////////
	//FZB copper coil spacers
	/////////////////////////////////////
	//there are 4 spacers are used to support|fixed the vacuum
	//They thickness and density are of important since they are in the way that bremstrulung
	//photons come out
	//In the engineering drawing, the material is aluminum alloy bond with contact cemont on the surface
	// I am not going to build the cement here

	double pFZBCoilSpacerX=1.75*inch;
	//pFZBCoilSpacerY is 3.25" in the drawing, it should be about 3.64" in order to fill up the gap, 
	//I think the rest is cement, Here I fill it up
	//double pFZBCoilSpacerY=3.25*inch;   
	double pFZBCoilSpacerY=pFZBCoilY/2-pFZBPlateMidSideY/2;
	double pFZBCoilSpacerZ=8.00*inch;
	G4VSolid* FZBCoilSpacerSolid = new G4Box("FZBCoilSpacerBox",
		pFZBCoilSpacerX/2.0,pFZBCoilSpacerY/2.0,pFZBCoilSpacerZ/2.0);

	G4LogicalVolume* FZBCoilSpacerLogical = new G4LogicalVolume(FZBCoilSpacerSolid,
		mMaterialManager->aluminum,"FZBCoilSpacerLogical",0,0,0);
	FZBCoilSpacerLogical->SetVisAttributes(LightBlueVisAtt);  

	//place 4 copies of this spacer at the corner in each FZ magnet
	double pFZBCoilPosY=pFZBPlateMidSideY/2+pFZBCoilSpacerY/2;
	double pFZBCoilPosZ=pFZBZ/2-pFZBCoilSpacerZ/2;
	G4ThreeVector pSpacer1Pos3V=G4ThreeVector(0, pFZBCoilPosY, pFZBCoilPosZ);
	G4ThreeVector pSpacer2Pos3V=G4ThreeVector(0,-pFZBCoilPosY, pFZBCoilPosZ);
	G4ThreeVector pSpacer3Pos3V=G4ThreeVector(0,-pFZBCoilPosY,-pFZBCoilPosZ);
	G4ThreeVector pSpacer4Pos3V=G4ThreeVector(0, pFZBCoilPosY,-pFZBCoilPosZ);
	new G4PVPlacement(0,pSpacer1Pos3V,
		FZBCoilSpacerLogical,"FZBCoilSpacerPhys",FZBContainerLogical,true,0,0);
	new G4PVPlacement(0,pSpacer2Pos3V,
		FZBCoilSpacerLogical,"FZBCoilSpacerPhys",FZBContainerLogical,true,1,0);
	new G4PVPlacement(0,pSpacer3Pos3V,
		FZBCoilSpacerLogical,"FZBCoilSpacerPhys",FZBContainerLogical,true,2,0);
	new G4PVPlacement(0,pSpacer4Pos3V,
		FZBCoilSpacerLogical,"FZBCoilSpacerPhys",FZBContainerLogical,true,3,0);


	/////////////////////////////////////
	//setup a local dump
	if(mSetupBeamDump)
	  {
	    //double pFZDumpWidth=20.0*cm;
	    //double pFZDumpHeight=60.0*cm;
	    //double pFZDumpThick=15.0*cm;
	    //double pBDY_pos = -20.0*cm;
	    //double pBDZ_pos = mPivotZOffset - 80.0*cm;
	    
	    double pFZDumpWidth=mBeamDumpWidth;
	    double pFZDumpHeight=mBeamDumpHeight;
	    double pFZDumpThick=mBeamDumpThick;
	    double pBDY_pos = mPivot2BeamDumpY;
	    double pBDZ_pos = mPivotZOffset + mPivot2BeamDumpZ;
	    
	    
	    G4VSolid* FZDumpBox = new G4Box("FZDumpBox",
					    pFZDumpWidth/2.0, pFZDumpHeight/2.0, pFZDumpThick/2.0);
	    
	    //drill a 3mm diameter hole throught it
	    G4VSolid* FZDumpCollimatorTub = new G4Tubs("FZDumpCollimatorTub",
						       0,3*mm,pFZDumpThick/2+1*mm,0*deg,360.0*deg);
	    ///FZDumpBox  subtract the CollimatorTub
	    G4SubtractionSolid* FZDumpSolid = new G4SubtractionSolid("FZDumpSolid",
								     FZDumpBox, FZDumpCollimatorTub, 0, G4ThreeVector(0,-pBDY_pos,0));
	    
	    G4LogicalVolume* FZDumpLogical = new G4LogicalVolume(FZDumpSolid, 
								 mMaterialManager->tungsten,"FZDumpLogical",0,0,0);
	    SDman->AddNewDetector(FZDumpSD);
	    FZDumpLogical->SetSensitiveDetector(FZDumpSD);
	    FZDumpLogical->SetVisAttributes(PurpleVisAtt); 
	    
	    new G4PVPlacement(0,G4ThreeVector(0,pBDY_pos,pBDZ_pos),
			      FZDumpLogical,"FZDumpPhys",motherLogical,0,0,0);
	  }

	//virtual detector
	/////////////////////////////////////

	//mSetupChicaneVD = 0 menas none,
	if(mSetupChicaneVD)
	{
	       double pFZBVDWidth=60.0*cm;
	       double pFZBVDHeight=80.0*cm;
	       double pFZBVDThick=20.0*cm;
	       G4VSolid* FZB2VDSolid = new G4Box("FZB2VDBox",
						 pFZBVDWidth/2.0, pFZBVDHeight/2.0, pFZBVDThick/2.0);
	       
	       
	       G4LogicalVolume* FZB2VDLogical = new G4LogicalVolume(FZB2VDSolid, 
								    mMaterialManager->air,"FZB2VDLogical",0,0,0); 
	       SDman->AddNewDetector(FZB2SD);
	       FZB2VDLogical->SetSensitiveDetector(FZB2SD);
	       FZB2VDLogical->SetVisAttributes(LightYellowVisAtt);
	       
	       if(mSetupChicaneVD!=0)
		 {
		   new G4PVPlacement(0,G4ThreeVector(0,0*cm,120*cm+pFZB2PosZ),
				     FZB2VDLogical,"FZB2VDPhys",motherLogical,0,0,0);
		 }
	}

	return FZB2FieldContainerPhys;
}


G4VPhysicalVolume* WACSDetectorConstruction::ConstructWACSPlatform(G4LogicalVolume* motherLogical)
{
	const double inch=2.54*cm;
	G4RotationMatrix *pRotX90deg=new G4RotationMatrix();
	pRotX90deg->rotateX(-90*deg); 
	G4RotationMatrix *pRotY90deg=new G4RotationMatrix();
	pRotY90deg->rotateY(-90*deg); 

	/////////////////////////////////////
	//poles to support 3rd floor
	/////////////////////////////////////
	double PFPoleWidth=4.0*inch;	//x
	double PFPoleHeight=4.0*inch;   //y
	double PFPoleLength=48.0*inch;  //z
	G4VSolid* PFPoleSolid = new G4Box("PFPoleBox",PFPoleWidth/2.0,
		PFPoleHeight/2.0,PFPoleLength/2.0);

	G4LogicalVolume* PFPoleLogical = new G4LogicalVolume(PFPoleSolid,
		mMaterialManager->aluminum,"PFPoleLogical",0,0,0);
	PFPoleLogical->SetVisAttributes(LightGreenVisAtt); 

	double pPFPolePos_X=42.0*inch+PFPoleWidth/2.0;
	double pPFPolePos_Y=-1.0*inch;
	double pPFPoleUpPos_Z=-70.53*inch+PFPoleHeight/2.0;
	double pPFPoleDownPos_Z=32.0*inch+PFPoleHeight/2.0;
	G4VPhysicalVolume *PFPolePhys = new G4PVPlacement(pRotX90deg,
		G4ThreeVector(pPFPolePos_X,pPFPolePos_Y,pPFPoleUpPos_Z),
		PFPoleLogical,"PFPolePhys",motherLogical,true,0,0);	//up left
	new G4PVPlacement(pRotX90deg,
		G4ThreeVector(-pPFPolePos_X,pPFPolePos_Y,pPFPoleUpPos_Z),
		PFPoleLogical,"PFPolePhys",motherLogical,true,1,0);	//up right
	new G4PVPlacement(pRotX90deg,
		G4ThreeVector(-pPFPolePos_X,pPFPolePos_Y,pPFPoleDownPos_Z),
		PFPoleLogical,"PFPolePhys",motherLogical,true,0,0);	//down right
	new G4PVPlacement(pRotX90deg,
		G4ThreeVector(pPFPolePos_X,pPFPolePos_Y,pPFPoleDownPos_Z),
		PFPoleLogical,"PFPolePhys",motherLogical,true,0,0); //down left


	double pSideBoardWidth=0.75*inch;
	double pSideBoardHeight=4.0*inch;


	int pSetupPlatform=0;
	gConfig->GetParameter("SetupPlatform",pSetupPlatform);
	if(pSetupPlatform==1 || pSetupPlatform==3)
	{
		/////////////////////////////////////
		//bars to support 3rd floor
		/////////////////////////////////////
		double PF3FBarWidth=4.0*inch;	//x
		double PF3FBarHeight=4.0*inch;   //y
		double PF3FBarLRLength=98.53*inch;  //z
		double PF3FBarUDLength=170.0*inch;  //z
		G4VSolid* PF3FBarLRSolid = new G4Box("PF3FBarLRBox",PF3FBarWidth/2.0,
			PF3FBarHeight/2.0,PF3FBarLRLength/2.0);
		G4VSolid* PF3FBarUDSolid = new G4Box("PF3FBarUDBox",PF3FBarWidth/2.0,
			PF3FBarHeight/2.0,PF3FBarUDLength/2.0);

		G4LogicalVolume* PF3FBarLRLogical = new G4LogicalVolume(PF3FBarLRSolid,
			mMaterialManager->aluminum,"PF3FBarLRLogical",0,0,0);
		PF3FBarLRLogical->SetVisAttributes(LightGreenVisAtt); 
		G4LogicalVolume* PF3FBarUDLogical = new G4LogicalVolume(PF3FBarUDSolid,
			mMaterialManager->aluminum,"PF3FBarUDLogical",0,0,0);
		PF3FBarUDLogical->SetVisAttributes(LightGreenVisAtt); 

		double pPF3FBarUDPos_X=0.0;
		double pPF3FBarPos_Y=25.5*inch;
		double pPF3FBarUpPos_Z=-70.53*inch+PF3FBarHeight/2.0;
		double pPF3FBarDownPos_Z=32.0*inch+PF3FBarHeight/2.0;
		double pPF3FBarLRPos_X=42.0*inch+PF3FBarWidth/2.0;
		double pPF3FBarLRPos_Z=(pPF3FBarUpPos_Z+pPF3FBarDownPos_Z)/2;

		new G4PVPlacement(0,
			G4ThreeVector(pPF3FBarLRPos_X,pPF3FBarPos_Y,pPF3FBarLRPos_Z),
			PF3FBarLRLogical,"PF3FBarLPhys",motherLogical,true,0,0);	// left
		new G4PVPlacement(0,
			G4ThreeVector(-pPF3FBarLRPos_X,pPF3FBarPos_Y,pPF3FBarLRPos_Z),
			PF3FBarLRLogical,"PF3FBarRPhys",motherLogical,true,1,0);	// right
		new G4PVPlacement(pRotY90deg,
			G4ThreeVector(pPF3FBarUDPos_X,pPF3FBarPos_Y,pPF3FBarUpPos_Z),
			PF3FBarUDLogical,"PF3FBarUpPhys",motherLogical,true,0,0);	//up 
		new G4PVPlacement(pRotY90deg,
			G4ThreeVector(pPF3FBarUDPos_X,pPF3FBarPos_Y,pPF3FBarDownPos_Z),
			PF3FBarUDLogical,"PF3FBarDownPhys",motherLogical,true,0,0); //down 


		/////////////////////////////////////
		// 3rd floor
		/////////////////////////////////////

		double pPF3FX=170.0*inch;
		double pPF3FY=0.5*inch;
		double pPF3FZ=106.53*inch;
		//3rd floor, need to dig hole (R=31") for target chamber
		G4VSolid* PF3rdFloorWholeSolid = new G4Box("PF3rdFloorWholeBox",pPF3FX/2.0,
			pPF3FY/2.0,pPF3FZ/2.0);
		G4VSolid* PF3rdFloorSubCircleSolid = new G4Tubs("PF3rdFloorSubCircleSolid",0,
			39.0*inch/2,0.55*inch/2,0*deg,360*deg);
		G4SubtractionSolid *PF3rdFloorSolid = new G4SubtractionSolid(
			"PF3rdFloorSolid",PF3rdFloorWholeSolid,PF3rdFloorSubCircleSolid,pRotX90deg,
			G4ThreeVector(0,0,mScatChamberZOffset-pPF3FBarLRPos_Z));
		G4LogicalVolume* PF3rdFloorLogical = new G4LogicalVolume(PF3rdFloorSolid,
			mMaterialManager->aluminum,"PF3rdFloorLogical",0,0,0);
		PF3rdFloorLogical->SetVisAttributes(LightGreenVisAtt); 


		new G4PVPlacement(0,G4ThreeVector(0,27.5*inch,pPF3FBarLRPos_Z),
			PF3rdFloorLogical,"PF3rdFloorPhys",motherLogical,0,0,0);	


		/////////////////////////////////////
		// 3rd floor side boards
		/////////////////////////////////////

		G4VSolid* PF3FSideLRSolid = new G4Box("PF3FSideLRSolid",pSideBoardWidth/2.0,
			pSideBoardHeight/2.0,pPF3FZ/2.0);
		G4VSolid* PF3FSideUDSolid = new G4Box("PF3FBarUDBox",pSideBoardWidth/2.0,
			pSideBoardHeight/2.0,pPF3FX/2.0-pSideBoardWidth);

		G4LogicalVolume* PF3FSideLRLogical = new G4LogicalVolume(PF3FSideLRSolid,
			mMaterialManager->aluminum,"PF3FSideLRLogical",0,0,0);
		PF3FSideLRLogical->SetVisAttributes(LightGreenVisAtt); 
		G4LogicalVolume* PF3FSideUDLogical = new G4LogicalVolume(PF3FSideUDSolid,
			mMaterialManager->aluminum,"PF3FSideUDLogical",0,0,0);
		PF3FSideUDLogical->SetVisAttributes(LightGreenVisAtt); 


		double pPF3FSideUDPos_X=0.0;
		double pPF3FSidePos_Y=27.5*inch+pSideBoardHeight/2.0;
		double pPF3FSideUpPos_Z=-70.53*inch+pSideBoardWidth/2.0;
		double pPF3FSideDownPos_Z=36.0*inch-pSideBoardWidth/2.0;
		double pPF3FSideLRPos_X=pPF3FX/2-pSideBoardWidth/2.0;
		double pPF3FSideLRPos_Z=(pPF3FSideUpPos_Z+pPF3FSideDownPos_Z)/2;

		new G4PVPlacement(0,
			G4ThreeVector(pPF3FSideLRPos_X,pPF3FSidePos_Y,pPF3FSideLRPos_Z),
			PF3FSideLRLogical,"PF3FSideLPhys",motherLogical,true,0,0);	// left
		new G4PVPlacement(0,
			G4ThreeVector(-pPF3FSideLRPos_X,pPF3FSidePos_Y,pPF3FSideLRPos_Z),
			PF3FSideLRLogical,"PF3FSideRPhys",motherLogical,true,1,0);	// right
		new G4PVPlacement(pRotY90deg,
			G4ThreeVector(pPF3FSideUDPos_X,pPF3FSidePos_Y,pPF3FSideUpPos_Z),
			PF3FSideUDLogical,"PF3FSideUpPhys",motherLogical,true,0,0);	//up 
		new G4PVPlacement(pRotY90deg,
			G4ThreeVector(pPF3FSideUDPos_X,pPF3FSidePos_Y,pPF3FSideDownPos_Z),
			PF3FSideUDLogical,"PF3FSideDownPhys",motherLogical,true,0,0); //down 
	}

	/////////////////////////////////////
	//rails to support septum
	/////////////////////////////////////

	double PF2FRailWidth=4.0*inch;	//x
	double PF2FRailHeight=4.0*inch;   //y
	double PF2FRailLRLength=98.53*inch;  //z
	double PF2FRailUDLength=92.0*inch;  //z
	G4VSolid* PF2FRailLRSolid = new G4Box("PF2FRailLRBox",PF2FRailWidth/2.0,
		PF2FRailHeight/2.0,PF2FRailLRLength/2.0);
	G4VSolid* PF2FRailUDSolid = new G4Box("PF2FRailUDBox",PF2FRailWidth/2.0,
		PF2FRailHeight/2.0,PF2FRailUDLength/2.0);

	G4LogicalVolume* PF2FRailLRLogical = new G4LogicalVolume(PF2FRailLRSolid,
		mMaterialManager->aluminum,"PF2FRailLRLogical",0,0,0);
	PF2FRailLRLogical->SetVisAttributes(SteelVisAtt); 
	G4LogicalVolume* PF2FRailUDLogical = new G4LogicalVolume(PF2FRailUDSolid,
		mMaterialManager->stainlesssteel,"PF2FRailUDLogical",0,0,0);
	PF2FRailUDLogical->SetVisAttributes(SteelVisAtt); 

	double pPF2FRailUDPos_X=0.0;
	double pPF2FRailPos_Y=-27.5*inch;
	double pPF2FRailUpPos_Z=-70.53*inch+PF2FRailHeight/2.0;
	double pPF2FRailDownPos_Z=32.0*inch+PF2FRailHeight/2.0;
	double pPF2FRailLRPos_X=42.0*inch+PF2FRailWidth/2.0;
	double pPF2FRailLRPos_Z=(pPF2FRailUpPos_Z+pPF2FRailDownPos_Z)/2;

	new G4PVPlacement(0,
		G4ThreeVector(pPF2FRailLRPos_X,pPF2FRailPos_Y,pPF2FRailLRPos_Z),
		PF2FRailLRLogical,"PF2FRailLPhys",motherLogical,true,0,0);	// left
	new G4PVPlacement(0,
		G4ThreeVector(-pPF2FRailLRPos_X,pPF2FRailPos_Y,pPF2FRailLRPos_Z),
		PF2FRailLRLogical,"PF2FRailRPhys",motherLogical,true,1,0);	// right
	new G4PVPlacement(pRotY90deg,
		G4ThreeVector(pPF2FRailUDPos_X,pPF2FRailPos_Y,pPF2FRailUpPos_Z),
		PF2FRailUDLogical,"PF2FRailUpPhys",motherLogical,true,0,0);	//up 
	new G4PVPlacement(pRotY90deg,
		G4ThreeVector(pPF2FRailUDPos_X,pPF2FRailPos_Y,pPF2FRailDownPos_Z),
		PF2FRailUDLogical,"PF2FRailDownPhys",motherLogical,true,0,0); //down 


	if(pSetupPlatform==1 || pSetupPlatform==2)
	{
		/////////////////////////////////////
		// 2nd floor
		/////////////////////////////////////

		double pPF2FX=242.5*inch;
		double pPF2FY=0.5*inch;
		double pPF2FZ=106.53*inch;
		//3rd floor, need to dig hole (R=31") for target chamber
		G4VSolid* PF2ndFloorSolid = new G4Box("PF2ndFloorBox",pPF2FX/2.0,
			pPF2FY/2.0,pPF2FZ/2.0);

		G4LogicalVolume* PF2ndFloorLogical = new G4LogicalVolume(PF2ndFloorSolid,
			mMaterialManager->aluminum,"PF2ndFloorLogical",0,0,0);
		PF2ndFloorLogical->SetVisAttributes(LightGreenVisAtt); 


		new G4PVPlacement(0,G4ThreeVector(0,-38.38*inch,pPF2FRailLRPos_Z),
			PF2ndFloorLogical,"PF2ndFloorPhys",motherLogical,0,0,0);	


		/////////////////////////////////////
		// 2nd floor side boards
		/////////////////////////////////////

		G4VSolid* PF2FSideLRSolid = new G4Box("PF2FSideLRSolid",pSideBoardWidth/2.0,
			pSideBoardHeight/2.0,pPF2FZ/2.0);
		G4VSolid* PF2FSideUDSolid = new G4Box("PF2FBarUDBox",pSideBoardWidth/2.0,
			pSideBoardHeight/2.0,pPF2FX/2.0-pSideBoardWidth);

		G4LogicalVolume* PF2FSideLRLogical = new G4LogicalVolume(PF2FSideLRSolid,
			mMaterialManager->aluminum,"PF2FSideLRLogical",0,0,0);
		PF2FSideLRLogical->SetVisAttributes(LightGreenVisAtt); 
		G4LogicalVolume* PF2FSideUDLogical = new G4LogicalVolume(PF2FSideUDSolid,
			mMaterialManager->aluminum,"PF2FSideUDLogical",0,0,0);
		PF2FSideUDLogical->SetVisAttributes(LightGreenVisAtt); 


		double pPF2FSideUDPos_X=0.0;
		double pPF2FSidePos_Y=-38.13*inch+pSideBoardHeight/2.0;
		double pPF2FSideUpPos_Z=-70.53*inch+pSideBoardWidth/2.0;
		double pPF2FSideDownPos_Z=36.0*inch-pSideBoardWidth/2.0;
		double pPF2FSideLRPos_X=pPF2FX/2-pSideBoardWidth/2.0;
		double pPF2FSideLRPos_Z=(pPF2FSideUpPos_Z+pPF2FSideDownPos_Z)/2;

		new G4PVPlacement(0,
			G4ThreeVector(pPF2FSideLRPos_X,pPF2FSidePos_Y,pPF2FSideLRPos_Z),
			PF2FSideLRLogical,"PF2FSideLPhys",motherLogical,true,0,0);	// left
		new G4PVPlacement(0,
			G4ThreeVector(-pPF2FSideLRPos_X,pPF2FSidePos_Y,pPF2FSideLRPos_Z),
			PF2FSideLRLogical,"PF2FSideRPhys",motherLogical,true,1,0);	// right
		new G4PVPlacement(pRotY90deg,
			G4ThreeVector(pPF2FSideUDPos_X,pPF2FSidePos_Y,pPF2FSideUpPos_Z),
			PF2FSideUDLogical,"PF2FSideUpPhys",motherLogical,true,0,0);	//up 
		new G4PVPlacement(pRotY90deg,
			G4ThreeVector(pPF2FSideUDPos_X,pPF2FSidePos_Y,pPF2FSideDownPos_Z),
			PF2FSideUDLogical,"PF2FSideDownPhys",motherLogical,true,0,0); //down 
	}

	return PFPolePhys;
}

