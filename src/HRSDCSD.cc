// ********************************************************************
//
// $Id: HRSDCSD.cc,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//
//By Jixie: The DC SD is just for position detector, good for DC and SC. 
//It will do the following: 
//1) record the inpos and inmom from the prestep of current hit
//2) record the time, outpos,outmom from the poststep of current hit 
//3) sign deposited energy to non-ionized energy if it is a transportation or msc   
//In the event action, these hits will be stored into the root ntuple


#include "HRSDCSD.hh"
#include "HRSDCHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"

HRSDCSD::HRSDCSD(G4String name):G4VSensitiveDetector(name)
{
	G4String HCname;
	collectionName.insert(HCname="DCColl");
	HCID = -1;

	//The given 'name' will be stored into G4VSensitiveDetector::SensitiveDetectorName
	//you can give '/RTPCGEM' or 'RTPCGEM', the Geant4 system will both set 'RTPCGEM'
	//as the SensitiveDetectorName and '/RTPCGEM' as fullPathName
	//
	//For current hit colection, the collectionName is "SensitiveDetectorName/HCname"
	//To access this HC in event action, one need to get HCID through
	//G4HCtable::GetCollectionID(G4String collectionName) 

	hitsCollection = 0;
}

HRSDCSD::~HRSDCSD()
{
	;
}

void HRSDCSD::Initialize(G4HCofThisEvent* HCE)
{
	hitsCollection = new HRSDCHitsCollection(SensitiveDetectorName,collectionName[0]);
	if(HCID<0)
	{ 
		HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection); 
	}
	HCE->AddHitsCollection(HCID,hitsCollection);
}


G4bool HRSDCSD::ProcessHits(G4Step* aStep,G4TouchableHistory* /*aROHist*/)
{	
	if(!hitsCollection) return false;

	G4double edep = aStep->GetTotalEnergyDeposit();
	if(edep/keV<=0.)  return true;


	G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
	G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
	G4TouchableHistory* theTouchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
	G4VPhysicalVolume* thePhysical = theTouchable->GetVolume();
	
	G4int copyNo = thePhysical->GetCopyNo();
	G4double hitTime = preStepPoint->GetGlobalTime();

	HRSDCHit* aHit = new HRSDCHit(copyNo,hitTime);

	G4double edep_NonIon = 0;
	G4String process=postStepPoint->GetProcessDefinedStep()->GetProcessName();
	if(process=="Transportation" || process=="msc")
		edep_NonIon=edep;	
	//if(edep_NonIon) G4cout<<"edep_NonIon="<<edep_NonIon<<G4endl;

	G4int pdgid=aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
	if(pdgid>3000)
	{
		//G4cout<<"Wrong PDGCCoding pdgid="<<pdgid
		//	<<" name="<<aStep->GetTrack()->GetParticleDefinition()->GetParticleName()<<G4endl;
		pdgid=(pdgid%1000000)/10;
	}

	aHit->SetPhysV(theTouchable->GetVolume());
	aHit->SetPdgid(pdgid);
	aHit->SetEdep(edep);
	aHit->SetNonIonEdep(edep_NonIon);
	aHit->SetInPos(preStepPoint->GetPosition());
	aHit->SetInMom(preStepPoint->GetMomentum());
	aHit->SetOutPos(postStepPoint->GetPosition());
	aHit->SetOutMom(postStepPoint->GetMomentum());
	aHit->SetTrackId(aStep->GetTrack()->GetTrackID());
	aHit->SetParentTrackId(aStep->GetTrack()->GetParentID());

	hitsCollection->insert(aHit);

	return true;
}

void HRSDCSD::EndOfEvent(G4HCofThisEvent* HCE)
{
	if(!HCE) return;  //no hits collection found 
	//the above line just to avoid warning of not use HCE

	int nhitC = hitsCollection->GetSize();
	if(!nhitC) return;
	
	int HIT_VERBOSITY=0;
	if(HIT_VERBOSITY >= 2 && nhitC)
	{
		G4cout<<"<<<End of Hit Collections <" << SensitiveDetectorName << "/"<<collectionName[0]
		<<">: " << nhitC << " hits." << G4endl; 
		for (int i=0; i<nhitC; i++)
		{
			(*hitsCollection)[i]->Print();
		}
	}

	return;
}

