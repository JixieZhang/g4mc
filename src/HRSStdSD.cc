// ********************************************************************
//
// $Id: HRSStdSD.cc,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//
#include "HRSStdSD.hh"
#include "HRSStdHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"
//By Jixie: The std SD is just a flux detector or calorimeter detector. 
//It will not pay attention to the time window. It will do the following: 
//1) record the inpos and inmom from the first hit
//2) record the time, outpos,outmom from the last hit
//3) accumulate deposited energy (ionized or non-ionized) from all hits which 
//   are belonging to the same track.
//In the event action, these hits will be stored into the root ntuple


HRSStdSD::HRSStdSD(G4String name):G4VSensitiveDetector(name)
{
	G4String HCname;
	collectionName.insert(HCname="StdColl");
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

HRSStdSD::~HRSStdSD()
{
	;
}

void HRSStdSD::Initialize(G4HCofThisEvent* HCE)
{
	hitsCollection = new HRSStdHitsCollection(SensitiveDetectorName,collectionName[0]);
	if(HCID<0)
	{ 
		HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection); 
	}
	HCE->AddHitsCollection(HCID,hitsCollection);
}


G4bool HRSStdSD::ProcessHits(G4Step* aStep,G4TouchableHistory* /*aROHist*/)
{	
	if(!hitsCollection) return false;

	G4double edep = aStep->GetTotalEnergyDeposit();
	if(edep/keV<=0.)  return true;

	//By Jixie; it turns out that this is always zero due to the fact that G4 physics process never set 
	//a value to it
	G4double edep_NonIon = aStep->GetNonIonizingEnergyDeposit();

	G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
	G4TouchableHistory* theTouchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
	G4int copyNo = theTouchable->GetVolume()->GetCopyNo();
	G4int trackId = aStep->GetTrack()->GetTrackID();
	G4double hitTime = preStepPoint->GetGlobalTime();
	
	G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
	G4ThreeVector outpos = postStepPoint->GetPosition();       
	G4ThreeVector outmom = postStepPoint->GetMomentum();

	G4String process=postStepPoint->GetProcessDefinedStep()->GetProcessName();
	if(process=="Transportation" || process=="msc")
		edep_NonIon=edep;	
	//if(edep_NonIon) G4cout<<"edep_NonIon="<<edep_NonIon<<G4endl;

	// check if this finger already has a hit
	// generally a hit is the total energy deposit in a sensitive detector within a time window
	// it does not care which particle deposite these energy
	// In order to study the radiation, I need to isolate the hit source, therefore I require
	// comparing the trackid here

	HRSStdHit* aHit = 0;
	for(G4int i=0;i<hitsCollection->entries();i++)
	{
		aHit = (*hitsCollection)[i];
		if(aHit->GetID()==copyNo && aHit->GetTrackId()==trackId )
		{
			//found an exist hit 
			break;
		}
	}

	// if it has, then take the earlier time, accumulated the deposited energy
	// if not, create a new hit and set it to the collection
	if(aHit)
	{
		if(aHit->GetTime() > hitTime) aHit->SetTime(hitTime); 
		aHit->AddEdep(edep); 
		aHit->AddNonIonEdep(edep_NonIon);
		aHit->SetOutPos(outpos);
		aHit->SetOutMom(outmom);
	}
	else
	{
		//this is a new hit
		aHit = new HRSStdHit(copyNo,hitTime);

		G4int parentTrackId=aStep->GetTrack()->GetParentID();
		G4ThreeVector inpos = preStepPoint->GetPosition();       
		G4ThreeVector inmom = preStepPoint->GetMomentum();

		G4int pdgid=aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
		if(pdgid>3000)
		{
			//G4cout<<"Wrong PDGCCoding pdgid="<<pdgid
			//	<<" name="<<aStep->GetTrack()->GetParticleDefinition()->GetParticleName()<<G4endl;
			pdgid=(pdgid%1000000)/10;
		}

		aHit->SetPhysV(theTouchable->GetVolume());
		aHit->SetEdep(edep);
		aHit->SetNonIonEdep(edep_NonIon);
		aHit->SetTrackId(trackId);
		aHit->SetOutPos(outpos);
		aHit->SetOutMom(outmom);
		aHit->SetInPos(inpos);
		aHit->SetInMom(inmom);
		aHit->SetParentTrackId(parentTrackId);
		aHit->SetPdgid(pdgid);

		hitsCollection->insert(aHit);
	}

	return true;
}

void HRSStdSD::EndOfEvent(G4HCofThisEvent* HCE)
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

