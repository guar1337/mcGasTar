#include "PrimaryGeneratorAction.hh"
#include "TFile.h"


#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "ParticleInfo.hh"


#include <fstream>
#include <iostream>
#include <iomanip> 

	

	PrimaryGeneratorAction::PrimaryGeneratorAction()
	: G4VUserPrimaryGeneratorAction()
	{
		geometryID = 5;

		G4String beamFilePath;
		tar_thick = 5*mm;
		tar_angle = 33.0*TMath::DegToRad();
		beamFilePath = "/home/zalewski/dataTmp/beamSource/beamSource5.root";
		tar_pos_Z = 0.0;

		in_lvBeam = new TLorentzVector();
		TFile *inBeamF = new TFile{beamFilePath.data(),"READ"};
		if (!inBeamF->IsOpen()) G4Exception("PrimaryGeneratorAction::PrimaryGeneratorAction","beam sourceF not found", FatalException, ".");
		inBeamTree = (TTree*)inBeamF->Get("beamSource");
		inBeamTree->SetMakeClass(0);
		inBeamTree->SetBranchAddress("lvBeam", &in_lvBeam);

		inBeamTree->SetBranchAddress("MWPC_1_X", &in_MWPC_1_X);
		inBeamTree->SetBranchAddress("MWPC_2_X", &in_MWPC_2_X);
		inBeamTree->SetBranchAddress("MWPC_1_Y", &in_MWPC_1_Y);
		inBeamTree->SetBranchAddress("MWPC_2_Y", &in_MWPC_2_Y);

		
		// default particle kinematic
		particletable = G4ParticleTable::GetParticleTable();
		iontable = G4IonTable::GetIonTable();

		defNeut=particletable->FindParticle("neutron");

		ELC= new G4EmCalculator();
		G4NistManager* man = G4NistManager::Instance();
		//man->SetVerbose(0);
		silicon_material = man->FindOrBuildMaterial("G4_Si");

		excitedStateEnergy_6He=1797*keV;
		MWPC_equivalent_of_Si = 660*um;

		G4RotationMatrix *zeroRotMatrix = new G4RotationMatrix();
		deutDiscTube = new G4Tubs("deutDiscTube", 	0.0*mm, 12.5*mm,			//inner & outer radius
																			20.0*mm,						//length
																			0.0*rad, CLHEP::twopi);	//starting & ending angle

		deutSphere = new G4Sphere("deutSphere", 	0.0*mm, 78.63*mm,			//inner & outer radius
																			0.0*rad, 4*CLHEP::pi,	//starting & ending phi
																			0.0*rad, 12.95*2*deg);	//starting & ending theta

		sphereCutoff = new G4Box("sphereCutoff", 100.0*mm, 100.0*mm, 77.63*mm);
		deutCap = new G4SubtractionSolid("deutSphere-sphereCutoff", deutSphere, sphereCutoff);
		G4ThreeVector tempUnionDeutVector(0.0,0.0,(-76.63+1.0)*mm);
		tempDeuterUnion = new G4UnionSolid("deutDiscTube+deutCap", deutDiscTube, deutCap, zeroRotMatrix, tempUnionDeutVector);
		G4RotationMatrix *deutCapRot = new G4RotationMatrix();
		deutCapRot->rotateY(CLHEP::pi);
		G4ThreeVector secondCapShiftVect(0.0,0.0,(76.63-1.0)*mm);
		zAxis = G4ThreeVector(0.0,0.0,1.0);
		gasCellSolid = new G4UnionSolid("deutCapS+deutDiscTube", tempDeuterUnion, deutCap, deutCapRot, secondCapShiftVect);
		delete zeroRotMatrix;
		delete deutCapRot;

	}

		
	PrimaryGeneratorAction::~PrimaryGeneratorAction()
	{
		delete ELC;
	}

void
PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	inputEventNo = anEvent->GetEventID() - entryShift;	
	if (inputEventNo == 5942259)
	{
		inputTreeLoopCounter++;
		entryShift = inputTreeLoopCounter * 5942259;
		//printf("Loop number %d for 1st geo\nTree entry numer: %d\n",inputTreeLoopCounter, anEvent->GetEventID());
	}
	inBeamTree->GetEntry(inputEventNo);

	def6He = iontable->GetIon(2,6);
	def4He = iontable->GetIon(2,4);
	defTar = iontable->GetIon(1,cs::tarMass);

	mass6He=def6He->GetPDGMass();
	mass4He=def4He->GetPDGMass();
	massTar=defTar->GetPDGMass();
	massNeut=defNeut->GetPDGMass();
	
	//generate vertex position
	MWPC_1_X = in_MWPC_1_X + CLHEP::RandFlat::shoot(0.0,1.25)-0.625;
	MWPC_1_Y = in_MWPC_1_Y + CLHEP::RandFlat::shoot(0.0,1.25)-0.625;

	MWPC_2_X = in_MWPC_2_X + CLHEP::RandFlat::shoot(0.0,1.25)-0.625;
	MWPC_2_Y = in_MWPC_2_Y + CLHEP::RandFlat::shoot(0.0,1.25)-0.625;
	

	MWPC_1_Z = -816.0;
	MWPC_2_Z = -270.0;

	beam_T = (in_lvBeam->E()-in_lvBeam->M());
	dX=MWPC_2_X-MWPC_1_X;
	dY=MWPC_2_Y-MWPC_1_Y;
	dZ=MWPC_2_Z-MWPC_1_Z;

	double ene_beam = cs::mass6He + beam_T;
	double mom_beam = sqrt(ene_beam*ene_beam - cs::mass6He*cs::mass6He);
	G4ThreeVector *vBeam = new G4ThreeVector{dX, dY, dZ};
	vBeam->setMag(mom_beam);
	G4LorentzVector *lvBeam = new G4LorentzVector{*vBeam, in_lvBeam->E()};

	double tmpTarPos = tar_pos_Z + CLHEP::RandFlat::shoot(-tar_thick/2, tar_thick/2);
	Tcoef=(cos(tar_angle)*tmpTarPos-sin(tar_angle)*MWPC_1_X - cos(tar_angle)*MWPC_1_Z)/(sin(tar_angle)*dX+cos(tar_angle)*dZ);

	evX = MWPC_1_X + dX*Tcoef;
	evY = MWPC_1_Y + dY*Tcoef;
	evZ = MWPC_1_Z + dZ*Tcoef;

	//printf("evX: %f\tevY: %f\tevZ: %f\n", evX, evY, evZ);
	
	G4ThreeVector VertexPosition(evX,evY,evZ);
	G4ThreeVector minusvectBeam(-dX,-dY,-dZ);
	G4ThreeVector vectBeam(dX, dY, dZ);




	if (gasCellSolid->Inside(VertexPosition)==2)
	{
		maxDist = gasCellSolid->DistanceToOut(VertexPosition, vectBeam.unit());
		minDist = gasCellSolid->DistanceToOut(VertexPosition, minusvectBeam.unit());
		maxZ = maxDist*cos(zAxis.angle(vectBeam))*0.99;
		minZ = -minDist*cos(zAxis.angle(vectBeam))*0.99;
		G4double evZprym = CLHEP::RandFlat::shoot(minZ,maxZ);
		
		do
		{
			Tcoef = (evZprym-MWPC_1_Z)/dZ;

			evX = MWPC_1_X + dX*Tcoef;
			evY = MWPC_1_Y + dY*Tcoef;
			evZ = MWPC_1_Z + dZ*Tcoef;
		} while (gasCellSolid->Inside(VertexPosition)!=2);
		//printf("inside?: %d\n", gasCellSolid->Inside(VertexPosition));
	}

	
	//printf("evX: %f\tevY: %f\tevZ: %f\n", evX, evY, evZ);

	//Eloss estimation - CHECK
	//I can only use material which are in my setup, but how to get deut gas?
	//That's why I will use silicon equivalent of deuterium gas
	E_tar_loss = get_E(beam_T, MWPC_equivalent_of_Si, silicon_material);
	beam_T= E_tar_loss;

	E_tar_loss = get_E(beam_T, (tar_thick/2)*um*(1/cos(tar_angle)), silicon_material);
	beam_T= E_tar_loss;

	
	const G4LorentzVector lvTarget(0.0, 0.0, 0.0, massTar);


//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxx	START OF elastic PART	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	//ELASTIC SCATTERING VECTORS
	G4LorentzVector lv6He_EL(*lvBeam);
	G4LorentzVector lv2H_EL(lvTarget);
	G4LorentzVector lvCM_EL = lv6He_EL+lv2H_EL;
	G4ThreeVector boostVect_EL = lvCM_EL.boostVector();

	lv6He_EL.boost(-boostVect_EL);
	lv2H_EL.boost(-boostVect_EL);
	//CM PART
	double theta_CM = acos(CLHEP::RandFlat::shoot(-1.0, 1.0));
	double phi_CM = CLHEP::RandFlat::shoot(0.0,2*double(CLHEP::pi));

	lv6He_EL.setTheta(CLHEP::pi-theta_CM);
	lv2H_EL.setTheta(theta_CM);

	lv6He_EL.setPhi(phi_CM);
	lv2H_EL.setPhi(CLHEP::pi-phi_CM);
	
	G4LorentzVector lv2H_CM_EL(lv2H_EL);
	G4LorentzVector lv6He_CM_EL(lv6He_EL);

	lv6He_EL.boost(boostVect_EL);
	lv2H_EL.boost(boostVect_EL);


	//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	//xxxxxxxxxxxxxx	ELASTIC PRIMARIES PART	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	//deuterium
	G4PrimaryParticle *PrimaryParticle_2H_EL = new G4PrimaryParticle(defTar); //2 PP
	PrimaryParticle_2H_EL->SetKineticEnergy(lv2H_EL.e()-lv2H_EL.m());
	PrimaryParticle_2H_EL->SetMomentumDirection(lv2H_EL.vect().unit());
	//I put all particles, that are not in PrimaryVertex into ParticleInfo of 2H
	//That is: lv2H_CM, lv6He_CM and lvBeam
	ParticleInfo *ElasticINFO = new ParticleInfo;
	PrimaryParticle_2H_EL->SetUserInformation(ElasticINFO);
	ElasticINFO->Set_LV_Beam(*lvBeam);
	ElasticINFO->Set_LV_2H_CM(lv2H_CM_EL);
	ElasticINFO->Set_LV_6He_CM(lv6He_CM_EL);

	ElasticINFO->Set_MWPC_1_X(in_MWPC_1_X);
	ElasticINFO->Set_MWPC_2_X(in_MWPC_2_X);
	ElasticINFO->Set_MWPC_1_Y(in_MWPC_1_Y);
	ElasticINFO->Set_MWPC_2_Y(in_MWPC_2_Y);

	ElasticINFO->Set_nx1(in_nx1);
	ElasticINFO->Set_nx2(in_nx2);
	ElasticINFO->Set_ny1(in_ny1);
	ElasticINFO->Set_ny2(in_ny2);

	//Helium 6
	G4PrimaryParticle *PrimaryParticle_6He_EL = new G4PrimaryParticle(def6He);//3 PP
	PrimaryParticle_6He_EL->SetKineticEnergy(lv6He_EL.e()-mass6He);	
	PrimaryParticle_6He_EL->SetMomentumDirection(lv6He_EL.vect().unit() );

	elasticVertex = new G4PrimaryVertex(VertexPosition,0);
	elasticVertex->SetPrimary(PrimaryParticle_2H_EL);
	elasticVertex->SetPrimary(PrimaryParticle_6He_EL);



//
//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxx	START OF INELastic PART	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	G4LorentzVector lv6He_IN=*lvBeam;

	G4LorentzVector lv2H_IN=lvTarget;
	G4LorentzVector lvCM_IN = lv6He_IN+lv2H_IN;

	G4ThreeVector boostVect_IN = lvCM_IN.boostVector();
	lv6He_IN.boost(-boostVect_IN);
	lv2H_IN.boost(-boostVect_IN);
	//Let's go inelastic now!
	G4double newEcm	= lv2H_IN.e() + lv6He_IN.e() - excitedStateEnergy_6He;
	G4double newE2H	= (1/(2*newEcm))*(newEcm*newEcm-mass6He*mass6He+lv2H_IN.m()*lv2H_IN.m());
	G4double newE6He	= (1/(2*newEcm))*(newEcm*newEcm+mass6He*mass6He-lv2H_IN.m()*lv2H_IN.m())+excitedStateEnergy_6He;
	G4double newMom2H	= sqrt(newE2H*newE2H-lv2H_IN.m()*lv2H_IN.m());
	G4double newMom6He = sqrt(newE6He*newE6He-mass6He*mass6He);

	lv2H_IN.setE(newE2H);
	lv6He_IN.setE(newE6He+excitedStateEnergy_6He);
	
	lv2H_IN.setRho(newMom2H);
	lv6He_IN.setRho(newMom6He);

	//CM PART
	lv6He_IN.rotate(0,	theta_CM,	0);
	lv6He_IN.rotate(0,	0,				phi_CM	);
	lv2H_IN.rotate(0,		-theta_CM,	0	);
	lv2H_IN.rotate(0,		0,				-phi_CM	);
		
	G4LorentzVector lv2H_CM_IN(lv2H_IN);
	G4LorentzVector lv6He_CM_IN(lv6He_IN);

	lv6He_IN.boost(boostVect_IN);
	lv2H_IN.boost(boostVect_IN);


	//NEUTRON EMISSION - negligible lifetime of excited 6He allows for immidiate decay
	G4LorentzVector lv6He_exc_IN = lv6He_IN;
	G4ThreeVector boostVect_6He = lv6He_exc_IN.boostVector();
	//neutrons are sitting and waiting for some action
	G4double decayE = mass6He + excitedStateEnergy_6He -(mass4He + 2*massNeut);
	G4double momRatio = CLHEP::RandFlat::shoot(0.1,1);

	//Setting neutron vectors
	//4He + 2n decay of 6He is classical - differences are negligible
	G4double neut1Theta_CM = acos(CLHEP::RandFlat::shoot(-1.0,1.0));
	G4double neut1Phi_CM = CLHEP::RandFlat::shoot(0.0,2*G4double(CLHEP::pi));
	G4double neut2Theta_CM = acos(CLHEP::RandFlat::shoot(-1.0,1.0));
	G4double neut2Phi_CM = CLHEP::RandFlat::shoot(0.0,2*G4double(CLHEP::pi));

	G4ThreeVector vneu1 (0,0,1);
	G4ThreeVector vneu2 (0,0,1);
	vneu1.rotate(0,	neut1Theta_CM,	0	);
	vneu1.rotate(0,	0,	neut1Phi_CM	);
	vneu2.rotate(0,	neut2Theta_CM,	0	);
	vneu2.rotate(0,	0,	neut2Phi_CM	);

	G4ThreeVector v4He (0,0,1);
	G4double cos_angle_between_neutrons = cos(vneu1.angle(vneu2));
	G4double L_fac = momRatio * momRatio + 2*momRatio * cos_angle_between_neutrons + 1;
	G4double momNeu1=sqrt(2*decayE*mass4He*massNeut/(mass4He*(1+momRatio*momRatio)+massNeut*L_fac));
	G4double momNeu2 = momRatio*momNeu1;
	G4double mom4He	= momNeu1*sqrt(L_fac);
	G4double eneNeu1 = momNeu1*momNeu1/(2*massNeut);
	G4double eneNeu2 = momNeu2*momNeu2/(2*massNeut);
	G4double ene4He	= mom4He*mom4He/(2*mass4He);
	
	G4LorentzVector lvNeut1(vneu1,massNeut+eneNeu1);
	G4LorentzVector lvNeut2(vneu2,massNeut+eneNeu2);
	G4LorentzVector lv4He_IN=(-lvNeut1-lvNeut2);
	lv4He_IN.setE(ene4He+mass4He);

	lv4He_IN.boost(boostVect_6He);
	lvNeut1.boost(boostVect_6He);
	lvNeut2.boost(boostVect_6He);
//
//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxx	INELASTIC PRIMARIES PART	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/*
	G4PrimaryParticle *PrimaryParticle_2H_IN = new G4PrimaryParticle(defTar); //0 PP
	PrimaryParticle_2H_IN->SetKineticEnergy(lv2H_IN.e()-lv2H_IN.m());
	PrimaryParticle_2H_IN->SetMomentumDirection(lv2H_IN.vect().unit());
	
	ParticleInfo *INelasticINFO = new ParticleInfo;
	PrimaryParticle_2H_IN->SetUserInformation(INelasticINFO);
	INelasticINFO->Set_LV_Beam(*lvBeam);
	INelasticINFO->Set_LV_2H_CM(lv2H_CM_IN);
	INelasticINFO->Set_LV_6He_CM(lv6He_CM_IN);

	G4PrimaryParticle *PrimaryParticle_4He_IN = new G4PrimaryParticle(def4He);
	PrimaryParticle_4He_IN->SetKineticEnergy(lv4He_IN.e()-mass4He);
	PrimaryParticle_4He_IN->SetMomentumDirection(lv4He_IN.vect().unit());
	
	G4PrimaryParticle *PrimaryParticle_Neutron1_IN=new G4PrimaryParticle(defNeut);
	PrimaryParticle_Neutron1_IN->SetKineticEnergy(lvNeut1.e()-massNeut);
	PrimaryParticle_Neutron1_IN->SetMomentumDirection(lvNeut1.vect().unit());

	G4PrimaryParticle * PrimaryParticle_Neutron2_IN=new G4PrimaryParticle(defNeut);
	PrimaryParticle_Neutron2_IN->SetKineticEnergy(lvNeut2.e()-massNeut);
	PrimaryParticle_Neutron2_IN->SetMomentumDirection(lvNeut2.vect().unit());
	
	inelasticVertex = new G4PrimaryVertex(VertexPosition,0);
		
	inelasticVertex->SetPrimary(PrimaryParticle_2H_IN);				//0
	inelasticVertex->SetPrimary(PrimaryParticle_4He_IN);				//1
	//inelasticVertex->SetPrimary(PrimaryParticle_Neutron1_IN);		//2
	//inelasticVertex->SetPrimary(PrimaryParticle_Neutron2_IN);		//3

	inelasticVertex->SetWeight(0.1);*/
	elasticVertex->SetWeight(1.0);
	anEvent->AddPrimaryVertex(elasticVertex);
	//anEvent->AddPrimaryVertex(inelasticVertex);
//
//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxx	END OF INELASTIC PART	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
/*
delete PrimaryParticle_2H_EL;
delete ElasticINFO;
delete PrimaryParticle_6He_EL;
delete elasticVertex;
*/


}
	
G4double PrimaryGeneratorAction::get_E(G4double E, G4double r, G4Material* mat)
{
	G4double R0=ELC->GetRange(E,def6He,mat);
	G4double R;
	R = R0 - r;
	G4double E1=ELC->GetKinEnergy(R,def6He,mat);
	return E1;
}

G4double PrimaryGeneratorAction::get_R(G4double E, G4Material* mat)
{
	G4double R0=ELC->GetRange(E,def6He,mat);
 
	return R0;
}
