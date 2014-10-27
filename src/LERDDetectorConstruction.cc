// ********************************************************************
//
// $Id: LERDDetectorConstruction.cc,v 1.00, 2014/10/24 LERD Exp $
// --------------------------------------------------------------
//
// ********************************************************************
//TODO: make a container for all layers, which subtract the target part
#include "LERDDetectorConstruction.hh"

#include "G4FieldManager.hh"
#include "G4MagneticField.hh"
#include "G4ChordFinder.hh"
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
#include "G4Polyhedra.hh"
#include "G4AssemblyVolume.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4UserLimits.hh"
#include "G4ios.hh"

#include "HRSStdSD.hh"
#include "HRSDCSD.hh"
#include "UsageManager.hh"
#include "HRSMaterial.hh"


extern UsageManager* gConfig;

////////////////////////////////////////////////////////////////////////////////////////////////
LERDDetectorConstruction::LERDDetectorConstruction(G4LogicalVolume *mother): 
mMotherLogVol(mother) 
{
	GetConfig();
	mMaterialManager=HRSMaterial::GetHRSMaterialManager();
	ConstructMaterials();

	G4cout<<"Contrstruct LERD geometry ... done! "<<G4endl;
}

LERDDetectorConstruction::~LERDDetectorConstruction()
{
	//I might need to delete the materials
	//But it does not hurt if I don't, since this class will have only one instance
	//in the program
	G4cout<<"Delete LERD geometry ... done! "<<G4endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////
void LERDDetectorConstruction::GetConfig()
{
	gConfig->ReadFile("Detector_LERD.ini");
	//////////////////////////////////////////////////////////////////

	const double mH2GasD_STP=0.08988*mg/cm3;
	const double mD2GasD_STP=0.180*mg/cm3;
	const double mHe3GasD_STP=0.1777*mg/cm3;
	const double mHeGasD_STP=0.1786*mg/cm3;
	const double mC4H10GasD_STP=2.52*mg/cm3*288.15/273.15;

	gConfig->GetParameter("TargetXOffset",mTargetXOffset);
	mTargetXOffset*=mm;
	gConfig->GetParameter("TargetYOffset",mTargetYOffset);
	mTargetYOffset*=mm;
	gConfig->GetParameter("TargetZOffset",mTargetZOffset);
	mTargetZOffset*=mm;

	gConfig->GetParameter("LERDLength",mLERDLength);mLERDLength*=mm;

	gConfig->GetParameter("D2GasL",mD2GasL);mD2GasL*=mm;
	gConfig->GetParameter("D2GasR",mD2GasR);mD2GasR*=mm;

	gConfig->GetParameter("TargetType",mTargetType);
	gConfig->GetParameter("D2GasT",mD2GasT);mD2GasT*=kelvin;
	gConfig->GetParameter("D2GasP",mD2GasP);mD2GasP*=atmosphere;
	double pTgGasD_STP = mD2GasD_STP;
	if(mTargetType==1) pTgGasD_STP = mH2GasD_STP;
	else if(mTargetType==3) pTgGasD_STP = mHe3GasD_STP;
	else if(mTargetType==4) pTgGasD_STP = mHeGasD_STP;
	mD2GasD=pTgGasD_STP*mD2GasP/atmosphere*273.15/mD2GasT;

	gConfig->GetParameter("HeGasT",mHeGasT);mHeGasT*=kelvin;
	gConfig->GetParameter("HeGasP",mHeGasP);mHeGasP*=atmosphere;
	mHeGasD=mHeGasD_STP*mHeGasP/atmosphere*273.15/mHeGasT;

	gConfig->GetParameter("WCGasT",mWCGasT);mWCGasT*=kelvin;
	gConfig->GetParameter("WCGasP",mWCGasP);mWCGasP*=atmosphere;
	mWCGasD=mC4H10GasD_STP*mWCGasP/atmosphere*273.15/mWCGasT;


	gConfig->GetParameter("TgWallMaterialType",mTgWallMaterialType);
	gConfig->GetParameter("TgWallThick",mTgWallThick);mTgWallThick*=mm;

	gConfig->GetParameter("1stWCR",m1stWCR);m1stWCR*=mm;
	gConfig->GetParameter("1stWCWindowThick",m1stWCWindowThick);m1stWCWindowThick*=mm;
	gConfig->GetParameter("1stWCThick",m1stWCThick);m1stWCThick*=mm;

	gConfig->GetParameter("2ndWCR",m2ndWCR);m2ndWCR*=mm;
	gConfig->GetParameter("2ndWCWindowThick",m2ndWCWindowThick);m2ndWCWindowThick*=mm;
	gConfig->GetParameter("2ndWCThick",m2ndWCThick);m2ndWCThick*=mm;


	gConfig->GetParameter("SiXR",mSiXR);mSiXR*=mm;
	gConfig->GetParameter("SiXThick",mSiXThick);mSiXThick*=mm;
	gConfig->GetParameter("SiYR",mSiYR);mSiYR*=mm;
	gConfig->GetParameter("SiYThick",mSiYThick);mSiYThick*=mm;
	gConfig->GetParameter("SCR",mSCR);mSCR*=mm;
	gConfig->GetParameter("SCThick",mSCThick);mSCThick*=mm;
	gConfig->GetParameter("LERDChamberR",mLERDChamberR);mLERDChamberR*=mm;
	gConfig->GetParameter("LERDChamberThick",mLERDChamberThick);mLERDChamberThick*=mm;

	gConfig->GetParameter("SetupSolenoid",mSetupSolenoid);
	gConfig->GetParameter("SolenoidPosX",mSolenoidPosX);mSolenoidPosX*=mm;
	gConfig->GetParameter("SolenoidPosY",mSolenoidPosY);mSolenoidPosY*=mm;
	gConfig->GetParameter("SolenoidPosZ",mSolenoidPosZ);mSolenoidPosZ*=mm;

	gConfig->GetParameter("BStepLimit",mBStepLimit);mBStepLimit*=mm;
	gConfig->GetParameter("DCStepLimit",mDCStepLimit);mDCStepLimit*=mm;


	gConfig->GetParameter("SetupEntranceNExitCap",mSetupEntranceNExitCap);
	gConfig->GetParameter("SetupEndPlateNCover",mSetupEndPlateNCover);

}


////////////////////////////////////////////////////////////////////////////////////////////////
void LERDDetectorConstruction::ConstructMaterials()
{
	//This is a stand alone version, most of this materials exist in HRSMaterial 
	G4double a;
	G4double z;
	G4double density;
	G4String name;
	G4String symbol;
	G4double pressure;
	G4double temperature;

	// elements for mixtures and compounds
	a = 1.01*g/mole;
	G4Element* elH = new G4Element(name="Hydrogen", symbol="H", z=1., a);
	a = 12.01*g/mole;
	G4Element* elC = new G4Element(name="Carbon", symbol="C", z=6., a);


	/////////////////////////////////////////////////////////
	//deuteriumGas 
	density = mD2GasD;
	a=2.014*g/mole;
	pressure=mD2GasP;
	temperature=mD2GasT;
	deuteriumGas = new G4Material(name="DeuteriumGas", z=1., a, density,
		kStateGas, temperature, pressure);


	/////////////////////////////////////////////////////////

	//Target Gas 1=H2, 2=D2,3=He3,4=He4

	a = 1.0081*g/mole;
	H2TgGas = new G4Material(name="H2TgGas", z=1., a, density=mD2GasD,
		kStateGas, temperature=mD2GasT, pressure=mD2GasP);
	a = 2.01410178*g/mole;
	D2TgGas = new G4Material(name="D2TgGas", z=1., a, density=mD2GasD,
		kStateGas, temperature=mD2GasT, pressure=mD2GasP);
	a = 3.0160293*g/mole;
	He3TgGas = new G4Material(name="He3TgGas", z=2., a, density=mD2GasD,
		kStateGas, temperature=mD2GasT, pressure=mD2GasP);
	a = 4.0026*g/mole;
	He4TgGas = new G4Material(name="He4TgGas", z=2., a, density=mD2GasD,
		kStateGas, temperature=mD2GasT, pressure=mD2GasP);

	//Helium Gas at 10 torr, or 1/76 ATM
	a = 4.0026*g/mole;
	density = mHeGasD ;
	pressure= mHeGasP ;
	temperature = mHeGasT ;
	heliumGas = new G4Material(name="HeliumGas", z=2., a, density,
		kStateGas, temperature, pressure);

	//////////////////////////////////////////////////////////////////////

	//WCGas isobutane at 10 torr, or 1/76 ATM
	a = 58.12*g/mole;
	pressure = mWCGasP;
	temperature = mWCGasT;
	density = mWCGasD;
	WCGas = new G4Material(name="WCGas", density, 2, kStateGas, temperature, pressure);
	WCGas -> AddElement(elH,10);
	WCGas -> AddElement(elC,4);
}


////////////////////////////////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* LERDDetectorConstruction::Construct()
{
	G4SDManager* SDMan = G4SDManager::GetSDMpointer();
	G4String SDname;

	G4VSensitiveDetector* LERDTargetSD = new HRSDCSD(SDname="LERDTargetGas");
	SDMan->AddNewDetector(LERDTargetSD);

	G4VSensitiveDetector* SDTargetWall = new HRSStdSD(SDname="LERDTargetWall");
	SDMan->AddNewDetector(SDTargetWall);

	G4VSensitiveDetector* SDWC1 = new HRSStdSD(SDname="LERDWC1");
	SDMan->AddNewDetector(SDWC1);

	G4VSensitiveDetector* SDWC2 = new HRSStdSD(SDname="LERDWC2");
	SDMan->AddNewDetector(SDWC2);

	G4VSensitiveDetector* SDSiX = new HRSStdSD(SDname="LERDSiX");
	SDMan->AddNewDetector(SDSiX);
	G4VSensitiveDetector* SDSiY = new HRSStdSD(SDname="LERDSiY");
	SDMan->AddNewDetector(SDSiY);
	G4VSensitiveDetector* SDSC = new HRSStdSD(SDname="LERDSC");
	SDMan->AddNewDetector(SDSC);

	G4VSensitiveDetector* SDLERDWall = new HRSStdSD(SDname="SDLERDWall");
	SDMan->AddNewDetector(SDLERDWall);

	G4VSensitiveDetector* virtualBoundarySD=new HRSStdSD("virtualBoundary");
	//############################################

	//////////////////////////
	//The mother box
	//////////////////////////

	//In this way we can place everything into the box as if placed them in the hall 
	// Size of this tub, large enough but not too big
	double LERDContainerL = mD2GasL+120*mm;       
	double LERDContainerRin = 0.0*m; 
	double LERDContainerRout = mLERDChamberR+50*mm; 

	if(mSetupSolenoid==1)  {LERDContainerRout=0.50*m; LERDContainerL=1.2*m;}
	else if(mSetupSolenoid>=2)  {LERDContainerRout=1.00*m; LERDContainerL=1.6*m;} 

	G4Tubs* LERDContainer = new G4Tubs("LERDContainer",LERDContainerRin,
		LERDContainerRout,LERDContainerL/2,0.0*deg,360.*deg);

	G4LogicalVolume* LERDContainerLogical = new G4LogicalVolume(LERDContainer, 
		mMaterialManager->air, "logLERDContainer", 0, 0, 0);
	LERDContainerLogical->SetVisAttributes(HallVisAtt);

	//the position at the hall
	double LERDContainerPosX=mTargetXOffset;
	double LERDContainerPosY=mTargetYOffset;
	double LERDContainerPosZ=mTargetZOffset;
	G4PVPlacement* phyLERDContainer= new G4PVPlacement(0,
		G4ThreeVector(LERDContainerPosX,LERDContainerPosY,LERDContainerPosZ),
		LERDContainerLogical,"LERDContainner",mMotherLogVol,0,0);


	//////////////////////////
	//The LERD ......//prepare some variables
	//////////////////////////
	//all radii in the ini file are inner radii, not center radii
	G4double Z_Half=mLERDLength/2;		
	G4double R40=m1stWCR;	
	G4double R80=m2ndWCR; 


	G4double phi_min,phi_max; //start and end phi angle
	G4double r_min,r_max;
	G4double z_down;

	G4RotationMatrix* RotZ90deg = new G4RotationMatrix();
	RotZ90deg->rotateZ(90.*deg);
	G4RotationMatrix* RotZ270deg = new G4RotationMatrix();
	RotZ270deg->rotateZ(270.*deg);	

	// set "user limits" for drawing smooth curve
	G4UserLimits* uDCStepLimits = new G4UserLimits(mDCStepLimit);
	G4UserLimits* uBStepLimits = new G4UserLimits(mBStepLimit);


	//////////////////////////
	//The solenoid
	//////////////////////////
	if(mSetupSolenoid)
	{
		ConstructSolenoid(LERDContainerLogical);
	}


	/////////////////////
	//Target 
	/////////////////////
	if(mTargetType==1) targetMaterial=H2TgGas;
	else if(mTargetType==3) targetMaterial=He3TgGas;
	else if(mTargetType==4) targetMaterial=He4TgGas;
	else targetMaterial=D2TgGas;
	G4VSolid* targetVesselSolid = new G4Tubs("targetVesselTubs",
		0.,mD2GasR,mD2GasL/2.0,0.,360.*deg);
	G4LogicalVolume* targetVesselLogical = new G4LogicalVolume(targetVesselSolid,
		targetMaterial,"targetVesselLogical",0,0,0);
	new G4PVPlacement(0,G4ThreeVector(),targetVesselLogical,
		"targetPhys",LERDContainerLogical,0,0);

	targetVesselLogical->SetSensitiveDetector(LERDTargetSD);

	/////////////////////
	//Target wall 
	/////////////////////
	//mTgWallMaterialType;  //1 is kapton, 2 is aluminum
	if(mTgWallMaterialType==2) targetWallMaterial=mMaterialManager->aluminum;
	else targetWallMaterial=mMaterialManager->kapton;
	G4VSolid* targetWallSolid = new G4Tubs("targetWallTubs",mD2GasR,
		mTgWallThick+mD2GasR,mD2GasL/2.0,0.,360.*deg);

	G4LogicalVolume* targetWallLogical = new G4LogicalVolume(targetWallSolid,
		targetWallMaterial,"targetWallLogical",0,0,0);

	new G4PVPlacement(0,G4ThreeVector(),targetWallLogical,
		"targetWallPhys",LERDContainerLogical,0,0);

	targetWallLogical->SetSensitiveDetector(SDTargetWall);
	if(mTgWallMaterialType==2)
		targetWallLogical->SetVisAttributes(GrayVisAtt);  //grey 
	else 
		targetWallLogical->SetVisAttributes(DarkRedVisAtt);  //dark red


	//////////////////////////
	//aluminum end caps 
	//////////////////////////
	//setup entrance cap and exit cap
	G4VSolid* entranceCoverSolid=0;
	if(mSetupEntranceNExitCap)
	{
		//////////////////////////
		//Downtream exit window aluminum cap 
		//////////////////////////
		//a cup of 3.1 mm height, including 10 mil thick bottom
		phi_min=0.*deg;
		phi_max=360.*deg;
		double pDownEndCapThick=10*0.0254*mm;
		const G4double  zPlane_cap[] ={-3*mm+mD2GasL/2,mD2GasL/2,mD2GasL/2,mD2GasL/2+pDownEndCapThick};
		const G4double  rInner_cap[] ={mD2GasR+mTgWallThick,mD2GasR+mTgWallThick,0.0*mm,0.0*mm};
		const G4double rOutner_cap[] ={mD2GasR+0.6*mm,mD2GasR+0.6*mm,mD2GasR+0.6*mm,mD2GasR+0.6*mm};
		G4VSolid* exitCapSolid = new G4Polycone("exitCapPcon",
			phi_min,phi_max,4,zPlane_cap,rInner_cap,rOutner_cap);
		G4LogicalVolume* exitCapLogical = new G4LogicalVolume(exitCapSolid,
			mMaterialManager->aluminum,"exitCapLogical",0,0,0);
		new G4PVPlacement(0,G4ThreeVector(),exitCapLogical,
			"exitCapPhys",LERDContainerLogical,0,0);
		exitCapLogical->SetVisAttributes(WhiteVisAtt);

		//////////////////////////
		//Uptream entrance aluminum cover  
		//////////////////////////
		phi_min=0.*deg;
		phi_max=360.*deg;
		double pUpEndCapThick=10*0.0254*mm;
		z_down=-mLERDLength/2+51.3*mm;  //Its down stream edge is 51.3*mm from LERD up edge
		r_min=mD2GasR+0.5*mm;  //3.5*mm
		const G4double  zPlane_u[] ={
			z_down,z_down-12*mm,z_down-12*mm,z_down-14*mm,z_down-14*mm,
			z_down-16*mm,z_down-21.75*mm,z_down-21.75*mm,z_down-55.1936*mm,z_down-55.1936*mm,
			z_down-65.1936*mm,z_down-65.1936*mm,z_down-86.3*mm,z_down-86.3*mm,
			z_down-86.3*mm-pUpEndCapThick
		};
		const G4double  rInner_u[] ={
			r_min-0.375*mm,r_min-0.375*mm,r_min-0.450*mm,r_min-0.450*mm,r_min,
			r_min,r_min,r_min,r_min,r_min,
			r_min,r_min,r_min,0,0
		};
		const G4double rOutner_u[] ={
			r_min-0.375*mm,r_min+3.226*mm,r_min+3.226*mm,r_min+4.993*mm,r_min+4.993*mm,
			r_min+5.760*mm,r_min+5.760*mm,r_min+8.50*mm,r_min+8.50*mm,r_min+11.50*mm,
			r_min+11.50*mm,r_min+8.50*mm,r_min+8.50*mm,r_min+8.50*mm,r_min+8.50*mm
		};
		entranceCoverSolid = new G4Polycone("entranceCoverPcon",
			phi_min,phi_max,15,zPlane_u,rInner_u,rOutner_u);
		G4LogicalVolume* entranceCoverLogical = new G4LogicalVolume(entranceCoverSolid,
			mMaterialManager->aluminum,"entranceCoverLogical",0,0,0);
		new G4PVPlacement(0,G4ThreeVector(),entranceCoverLogical,
			"entranceCoverPhys",LERDContainerLogical,0,0);

		entranceCoverLogical->SetVisAttributes(WhiteVisAtt); //white

		//place 10 mil Be as the beam exit window
		double pBeamExitWindowThick=10*0.0254*mm;
		G4VSolid* beamExitWindowSolid = new G4Tubs("beamExitWindowTubs",0,
			12.7*mm,pBeamExitWindowThick/2.0,0.,360.*deg);

		G4LogicalVolume* beamExitWindowLogical = new G4LogicalVolume(beamExitWindowSolid,
			mMaterialManager->beryllium,"beamExitWindowLogical",0,0,0);

		new G4PVPlacement(0,G4ThreeVector(0,0,z_down-86.3*mm-pUpEndCapThick-10*mm),
			beamExitWindowLogical, "beamExitWindowPgys", LERDContainerLogical,0,0);

		beamExitWindowLogical->SetVisAttributes(GrayVisAtt);  //grey 

	}

	
	//////////////////////////
	//LERD Low Pressure Chamber 
	//////////////////////////
	r_min=mLERDChamberR;
	r_max=mLERDChamberR+mLERDChamberThick;	
	phi_min=0*deg;
	phi_max=360.0*deg;

	G4VSolid* LERDChamberSolid
		= new G4Tubs("LERDChamberTubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* LERDChamberLogical
		= new G4LogicalVolume(LERDChamberSolid,mMaterialManager->aluminum,
		"LERDChamberLogical",0,0,0);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(),LERDChamberLogical,
		"LERDChamberPhys",LERDContainerLogical,0,0);

	LERDChamberLogical->SetSensitiveDetector(SDLERDWall);
	LERDChamberLogical->SetVisAttributes(PcbGreenVisAtt); //green, color for pcb


	//////////////////////////
	//Inner of the LERD Low Pressure Chamber 
	//////////////////////////
	//This is the containner that hold low presure helium gas
	//but not include the target straw
	r_min=mD2GasR+mTgWallThick;
	r_max=mLERDChamberR;	
	phi_min=0*deg;
	phi_max=360.0*deg;

	G4VSolid* innerLERDSolid = 0;
	G4VSolid* innerLERDSolid_whole = new G4Tubs("innerLERDSolid_whole",
		r_min,r_max,Z_Half,phi_min,phi_max);
	if(entranceCoverSolid)
	{
		innerLERDSolid = new G4SubtractionSolid("innerLERDSolid",
			innerLERDSolid_whole,entranceCoverSolid);
	}
	else
	{
		innerLERDSolid = innerLERDSolid_whole;
	}

	G4LogicalVolume* innerLERDLogical
		= new G4LogicalVolume(innerLERDSolid,heliumGas,"innerLERDLogical",0,0,0);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(),innerLERDLogical,
		"innerLERDPhys",LERDContainerLogical,0,0);

	innerLERDLogical->SetUserLimits(uBStepLimits);
	innerLERDLogical->SetVisAttributes(HallVisAtt); 


	/////////////////////
	//fisrt wire chamber
	/////////////////////
	//*****at r=40mm, 5mm thick drift gas (C4H10, 58.12 mg/cm^3 at 288k 1ATM, here it is  
	//10 torr or 0.0131 ATM at 300k. Window is 50 um Polypropylene([C3H6]n)

	phi_min=0*deg;
	phi_max=360*deg;

	r_min=R40;	
	r_max=R40+m1stWCWindowThick;	
	G4VSolid* firstWCWindow1Solid
		= new G4Tubs("firstWCWindow1Tubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* firstWCWindow1Logical
		= new G4LogicalVolume(firstWCWindow1Solid,mMaterialManager->polypropylene,
		"firstWCWindow1Logical",0,0,0);
	new G4PVPlacement(0,G4ThreeVector(),firstWCWindow1Logical,
		"firstWCWindow1Phys",innerLERDLogical,0,0);


	r_min=R40+m1stWCWindowThick;	
	r_max=R40+m1stWCThick+m1stWCWindowThick; 
	G4VSolid* firstWCSolid
		= new G4Tubs("firstWCTubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* firstWCLogical
		= new G4LogicalVolume(firstWCSolid,WCGas,"firstWCLogical",0,0,0);
	new G4PVPlacement(0,G4ThreeVector(),firstWCLogical,
		"firstWCPhys",innerLERDLogical,0,0);


	r_min=R40+m1stWCThick+m1stWCWindowThick;
	r_max=R40+m1stWCThick+2.0*m1stWCWindowThick;
	G4VSolid* firstWCWindow2Solid
		= new G4Tubs("firstWCWindow2Tubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* firstWCWindow2Logical
		= new G4LogicalVolume(firstWCWindow2Solid,mMaterialManager->polypropylene,
		"firstWCWindow2Logical",0,0,0);
	new G4PVPlacement(0,G4ThreeVector(),firstWCWindow2Logical,
		"firstWCWindow2Phys",innerLERDLogical,0,0);

	firstWCLogical->SetUserLimits(uDCStepLimits);
	firstWCLogical->SetSensitiveDetector(SDWC1);
	firstWCWindow1Logical->SetVisAttributes(WhiteVisAtt);   	//white
	firstWCLogical->SetVisAttributes(YellowVisAtt); 			//yellow
	firstWCWindow2Logical->SetVisAttributes(WhiteVisAtt);   	//white



	///////////////////// 
	//2nd wire chamber
	///////////////////// 
	//second wire chamber locate at r=80mm
	phi_min=0*deg;
	phi_max=360*deg;

	r_min=R80;	
	r_max=R80+m2ndWCWindowThick;	
	G4VSolid* secondWCWindow1Solid
		= new G4Tubs("secondWCWindow1Tubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* secondWCWindow1Logical
		= new G4LogicalVolume(secondWCWindow1Solid,
		mMaterialManager->polypropylene,"secondWCWindow1Logical",0,0,0);
	new G4PVPlacement(0,G4ThreeVector(),secondWCWindow1Logical,
		"secondWCWindow1Phys",innerLERDLogical,0,0);

	r_min=R80+m2ndWCWindowThick;	
	r_max=R80+m2ndWCWindowThick+m2ndWCThick; 
	G4VSolid* secondWCSolid
		= new G4Tubs("secondWCTubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* secondWCLogical
		= new G4LogicalVolume(secondWCSolid,WCGas,"secondWCLogical",0,0,0);
	new G4PVPlacement(0,G4ThreeVector(),secondWCLogical,
		"secondWCPhys",innerLERDLogical,0,0);

	r_min=R80+m2ndWCWindowThick+m2ndWCThick;
	r_max=R80+2.*m2ndWCWindowThick+m2ndWCThick;
	G4VSolid* secondWCWindow2Solid
		= new G4Tubs("secondWCWindow2Tubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* secondWCWindow2Logical
		= new G4LogicalVolume(secondWCWindow2Solid,
		mMaterialManager->polypropylene,"secondWCWindow2Logical",0,0,0);
	new G4PVPlacement(0,G4ThreeVector(),secondWCWindow2Logical,
		"secondWCWindow2Phys",innerLERDLogical,0,0);

	secondWCLogical->SetUserLimits(uDCStepLimits);
	secondWCLogical->SetSensitiveDetector(SDWC2);
	secondWCWindow1Logical->SetVisAttributes(WhiteVisAtt);  //white
	secondWCLogical->SetVisAttributes(YellowVisAtt);		//yellow
	secondWCWindow2Logical->SetVisAttributes(WhiteVisAtt);  //white


	//////////////////////////
	//Silicon layer 1
	//////////////////////////
	//0.1mm silicon located at r=100mm
	r_min=mSiXR;
	r_max=mSiXR+mSiXThick;
	phi_min=0*deg;
	phi_max=360.0*deg;
	//I want to use PGON for silicon other than tubs
	//G4VSolid* SiXLayerSolid
	//	= new G4Tubs("SiXLayerTubs",r_min,r_max,Z_Half,phi_min,phi_max);
	int nZplane_SiX=2;
	G4double zPlane_SiX[]={-Z_Half,Z_Half};
	G4double rInner_SiX[]={r_min,r_min};
	G4double rOuter_SiX[]={r_max,r_max};	
	G4VSolid* SiXLayerSolid = new G4Polyhedra("SiXLayer",
		phi_min,phi_max,6,nZplane_SiX,zPlane_SiX,rInner_SiX,rOuter_SiX);

	G4LogicalVolume* SiXLayerLogical = new G4LogicalVolume(SiXLayerSolid,
		mMaterialManager->silicon,"SiXLayerLogical",0,0,0);
	new G4PVPlacement(0,G4ThreeVector(),SiXLayerLogical,
		"SiXLayerPhys",innerLERDLogical,0,0);

	SiXLayerLogical->SetSensitiveDetector(SDSiX);
	SiXLayerLogical->SetVisAttributes(SkyBlueVisAtt);

	
	//////////////////////////
	//Silicon layer 2
	//////////////////////////
	//0.4mm silicon layer located at r=105mm
	r_min=mSiYR;
	r_max=mSiYR+mSiYThick;
	phi_min=0*deg;
	phi_max=360.0*deg;
	//I want to use PGON for silicon other than tubs
	//G4VSolid* SiYLayerSolid
	//	= new G4Tubs("SiYLayerTubs",r_min,r_max,Z_Half,phi_min,phi_max);
	int nZplane_SiY=2;
	G4double zPlane_SiY[]={-Z_Half,Z_Half};
	G4double rInner_SiY[]={r_min,r_min};
	G4double rOuter_SiY[]={r_max,r_max};	
	G4VSolid* SiYLayerSolid = new G4Polyhedra("SiYLayer",
		phi_min,phi_max,6,nZplane_SiY,zPlane_SiY,rInner_SiY,rOuter_SiY);

	G4LogicalVolume* SiYLayerLogical = new G4LogicalVolume(SiYLayerSolid,
		mMaterialManager->silicon,"SiYLayerLogical",0,0,0);
	new G4PVPlacement(0,G4ThreeVector(),SiYLayerLogical,
		"SiYLayerPhys",innerLERDLogical,0,0);

	SiYLayerLogical->SetSensitiveDetector(SDSiY);
	SiYLayerLogical->SetVisAttributes(SkyBlueVisAtt);


	//////////////////////////
	//Scintilator layer 
	//////////////////////////
	//SC located at r=100mm, 20 mm thick
	r_min=mSCR;
	r_max=mSCR+mSCThick;
	phi_min=0*deg;
	phi_max=360.0*deg;
	G4VSolid* SCLayerSolid
		= new G4Tubs("SCLayerTubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* SCLayerLogical
		= new G4LogicalVolume(SCLayerSolid,mMaterialManager->silicon,"SCLayerLogical",0,0,0);
	new G4PVPlacement(0,G4ThreeVector(),SCLayerLogical,
		"SCLayerPhys",innerLERDLogical,0,0);

	SCLayerLogical->SetSensitiveDetector(SDSC);
	SCLayerLogical->SetVisAttributes(YellowVisAtt);



	//############################################
	//////////////////////////
	// virtual boundary foil
	//////////////////////////
	//To speed up,
	//place virtual boundary foil here to kill all particles which penetrate the readout pcb
	//I put 7 z planes here such that I can choose the first 4 or the whole 7 planes
	phi_min=0.*deg;
	phi_max=360.*deg;
	//good for BoNuS
	double pDZ = 50.0*mm+mLERDLength/2; 
	double pRmax=mLERDChamberR+45.8*mm;  
	G4double zPlane_vb[]={-1.0*mm-pDZ,-pDZ,-pDZ,pDZ-32.0*mm,pDZ,pDZ,pDZ+1.0*mm};
	G4double rInner_vb[]={15.0*mm,15.0*mm,pRmax-0.1*mm,pRmax-0.1*mm,pRmax-0.1*mm,15.0*mm,15.0*mm};
	G4double rOuter_vb[]={pRmax,pRmax,pRmax,pRmax,pRmax,pRmax,pRmax};	//

	int nZplane=7;
	bool DoNotCoverDownEndPlate=false;
	int pSetupSuperBigBite=0;
	gConfig->GetParameter("SetupSuperBigBite",pSetupSuperBigBite);
	int pSetupBigBite=0;
	gConfig->GetParameter("SetupBigBite",pSetupBigBite); 
	if(pSetupBigBite || pSetupSuperBigBite)  DoNotCoverDownEndPlate=true;
	if(DoNotCoverDownEndPlate)  nZplane=4;

	G4VSolid* virtualBoundarySolid = new G4Polycone("virtualBoundaryPcon",
		phi_min,phi_max,nZplane,zPlane_vb,rInner_vb,rOuter_vb);
	G4LogicalVolume* virtualBoundaryLogical
		= new G4LogicalVolume(virtualBoundarySolid,mMaterialManager->air,"virtualBoundaryLogical",0,0,0);
	new G4PVPlacement(0,G4ThreeVector(),virtualBoundaryLogical,
		"virtualBoundaryPhys_LERD",LERDContainerLogical,0,0);

	virtualBoundaryLogical->SetSensitiveDetector(virtualBoundarySD);
	if(DoNotCoverDownEndPlate)
		virtualBoundaryLogical->SetVisAttributes(LightYellowVisAtt);
	else 
		virtualBoundaryLogical->SetVisAttributes(HallVisAtt);

	//############################################


	if (mSetupEndPlateNCover)
	{
		//now known yet
	}



	return phyLERDContainer;
}




//setup solenoid geometry
//if mSetupSolenoid==1, build DVCS solenoid, which is used in CLAS DVCS and BoNuS
//if mSetupSolenoid==2, build UVA solenoid, provided by Gorden Gates
//if mSetupSolenoid==3, CLAS12 solenoid, 
G4VPhysicalVolume* LERDDetectorConstruction::ConstructSolenoid(G4LogicalVolume *pMotherLogVol)
{
	double startphi=0.*deg, deltaphi=360.*deg;

	G4VPhysicalVolume* thePhysVol=0;

	//build the DVCS solenoid, Rin=230mm and Rout=912 mm
	//I build the whole body, not just the coil. Those numbers are measured from
	//pdf file: http://wwwold.jlab.org/Hall-B/secure/e1-dvcs/michel/sol/drawings/Ensemble.pdf
	//
	//|<--------- 510 mm -------->|<--228.8-->|	      
	//-----------------------------------------        --
	//|                                       |         132.7 mm in vertical
	//|                                       |        __
	//|                                      /                             p
	//|                                    /                               p
	//|                                  /                                 p 
	//-----------------------------------
	//|
	//|                           +(center)
	//|                        -->|     |<--106.2 mm 
	//-----------------------------------
	//|                                  \                                 p
	//|                                    \                               p 
	//|                                      \                             p
	//|                                       |
	//|                                       |  
	//-----------------------------------------

	G4VSolid* solenoidSolid = 0;
	G4LogicalVolume* solenoidLogical = 0;

	if(mSetupSolenoid==1)
	{
		//DVCS solenoid
		startphi=0.*deg; deltaphi=360.*deg;
		const int kNPlane_DVCSCoil=3;
		double rInner_DVCSCoil[] = {115*mm,115*mm,313.3*mm};
		double rOuter_DVCSCoil[] = {456*mm,456*mm,456*mm};
		double zPlane_DVCSCoil[] = {-510*mm,106.2*mm,228.8*mm};

		solenoidSolid = new G4Polycone("DVCSCoilPolycone",startphi,deltaphi,
			kNPlane_DVCSCoil,zPlane_DVCSCoil,rInner_DVCSCoil,rOuter_DVCSCoil);

	}
	else if(mSetupSolenoid==2)
	{
		//UVA Gorden Cates's solenoid
		startphi=0.*deg; deltaphi=360.*deg;
		//assuming it is a cylinder of Rin=200mm and Rout=851mm, L=1527mm
		solenoidSolid = new G4Tubs("UVASolenoidTubs",200.*mm,851.*mm,
			1527.*mm/2.0,startphi,deltaphi);
	}
	else if(mSetupSolenoid==3)
	{
		//CLAS12 solenoid, coil: 1.10m OD x 0.78m ID x 1.055m long
		//yoke: 1.96m OD x 1.10m ID x 1.18m long
		startphi=0.*deg; deltaphi=360.*deg;
		const int kNPlane_CLAS12Sol=3;
		double rInner_CLAS12Sol[] = {390*mm,390*mm,458.5*mm};
		double rOuter_CLAS12Sol[] = {1100*mm,1100*mm,1100*mm};
		double zPlane_CLAS12Sol[] = {-527.5*mm,448.7*mm,527.5*mm};

		solenoidSolid = new G4Polycone("CLAS12SolPolycone",startphi,deltaphi,
			kNPlane_CLAS12Sol,zPlane_CLAS12Sol,rInner_CLAS12Sol,rOuter_CLAS12Sol);
	}

	solenoidLogical = new G4LogicalVolume(solenoidSolid,
		mMaterialManager->stainlesssteel,"LERDCoilLogical",0,0,0);
	//solenoidLogical->SetVisAttributes(SteelVisAtt); 
	solenoidLogical->SetVisAttributes(DarkBlueVisAtt);   //dark blue looks better in 3D view

	//note that mSolenoidPosX,mSolenoidPosY,mSolenoidPosZ are in Hall coordinate,
	//if the pMotherLogVol is not the worldLogVol, one need to do the transfromation
	thePhysVol=new G4PVPlacement(0,
		G4ThreeVector(mSolenoidPosX-mTargetXOffset,mSolenoidPosY-mTargetYOffset,
		mSolenoidPosZ-mTargetZOffset),
		solenoidLogical,"SolenoidPhys",pMotherLogVol,false,0);

	return thePhysVol;
}

