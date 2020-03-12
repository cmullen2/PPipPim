#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "TDatabasePDG.h"
#include "clas12reader.h"
#include <Math/Vector4D.h>
#include <Math/Point3D.h>
#include <Math/DisplacementVector3D.h>
#include <Math/VectorUtil.h> //for boosts etc.




using namespace clas12;
using ROOT::Math::VectorUtil::boost;



typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<Double32_t> > HSLorentzVector;
typedef ROOT::Math::PositionVector3D< ROOT::Math::Cartesian3D< Double32_t >, ROOT::Math::DefaultCoordinateSystemTag > HSPosition;

typedef ROOT::Math::DisplacementVector3D< ROOT::Math::Cartesian3D< Double_t >, ROOT::Math::DefaultCoordinateSystemTag > HSMomentum;

HSLorentzVector fElin;
HSLorentzVector fElsc;
HSLorentzVector fGamma;
HSLorentzVector fCM;
HSLorentzVector fTar;
HSLorentzVector fMes;
HSLorentzVector fBar;
HSMomentum fCMBoost;

Double_t fCosTh=0;
Double_t fPhi=0;
Double_t W=0;
Double_t Q2=0;
Double_t t=0;
Double_t Pol=0;





void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	     rp->par()->getPz(),p4.M());

}           


Double_t StartTime(Double_t ptime, Double_t rftime){
  Double_t fBeamBucket = 4.0080160;
  //Double_t fSTimePeak = 44.125;
  Double_t fSTimePeak = 124.25;
  Double_t fStartTime = fSTimePeak - fBeamBucket*( (Int_t)(std::round( ((fSTimePeak-(ptime-rftime))/fBeamBucket)  ) )) + rftime;
  return fStartTime;
}


Double_t StartTime(clas12::region_part_ptr rp3, Double_t EDeltaTime, Double_t RF){
  if (rp3->getRegion()<2000){
    Double_t FTShift = 0; //fstimepeak - ? //See Haspect function FindTimeOffSetFT to determine number
    return StartTime(EDeltaTime + FTShift,RF ); 
  }
  else{
    return StartTime(EDeltaTime,RF);
  }

}

Double_t DeltaTime( clas12::region_part_ptr rp2){
  Double_t fTime= rp2->getTime();
  Double_t fPath;
  if (rp2->getRegion()<2000 ){ 
    Double_t ftx = rp2->ft(FTCAL)->getX();
    Double_t fty = rp2->ft(FTCAL)->getY();
    Double_t ftz = rp2->ft(FTCAL)->getZ();
    fPath = (sqrt(ftx*ftx + fty*fty + ftz*ftz))/100;
  }
  else{
    fPath = (rp2->getPath())/100;  //Path is in cm change to m
  }
  Double_t code=rp2->getPid();
  TParticlePDG* part=TDatabasePDG::Instance()->GetParticle(code);
  Double_t fPDGMass =part->Mass() ;
  Double_t rp2P = rp2->getP();
  Double_t HypBeta= rp2P/sqrt(rp2P*rp2P + fPDGMass*fPDGMass);
  Double_t HypTime = fPath/HypBeta/2.99792e+08*1E9; //Path m, Time ns
  Double_t SomeTime = fTime - HypTime;
  return SomeTime;

}


Double_t DeltaTime( clas12::region_part_ptr rp2, Double_t fTime){
  Double_t fPath;
  if (rp2->getRegion()<2000 ){ 
    Double_t ftx = rp2->ft(FTCAL)->getX();
    Double_t fty = rp2->ft(FTCAL)->getY();
    Double_t ftz = rp2->ft(FTCAL)->getZ();
    fPath = (sqrt(ftx*ftx + fty*fty + ftz*ftz))/100;
  }
  else{
    fPath = (rp2->getPath())/100;  //Path is in cm change to m
  }
  Double_t code=rp2->getPid();
  TParticlePDG* part=TDatabasePDG::Instance()->GetParticle(code);
  Double_t fPDGMass =part->Mass() ;
  Double_t rp2P = rp2->getP();
  Double_t HypBeta= rp2P/sqrt(rp2P*rp2P + fPDGMass*fPDGMass);
  Double_t HypTime = fPath/HypBeta/2.99792e+08*1E9; //Path m, Time ns
  Double_t SomeTime = fTime - HypTime;
  return SomeTime;

}

void CLAS12PPipPimReader(){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();


  /////////////////////////////////////
  TFile *outputFile = new TFile("/project/Gruppo3/fiber7/cmullen/ServiceWork/output/FirstRun.root","recreate");
  /////////////////////////////////////

  /////////////////////////////////////////////
  TTree *treeVars = new TTree("ParticleVars","SomeTree");

  ////ELECTRON
  Double_t ElectronP;
  Double_t ElectronTheta;
  Double_t ElectronPhi;
  Double_t ElectronTime;
  Double_t ElectronEdep;
  Double_t ElectronDeltaE;
  Double_t ElectronRegion;
  Double_t ElectronSector;
  //Double_t ElectronStatus=0;
  Double_t ElectronDeltaTime;
  Double_t ElectronDeltaTime2;
  Double_t ElectronDeltaTimeFromRecStartTime;
  Double_t EStartTime;
  Double_t ElectronStartTimeFromREC;
  Double_t RFTime;
  Double_t ElectronShiftedTime;
  Double_t NElectrons;
  Double_t ElectronShiftedTimeByRecStart;

  treeVars->Branch("ElectronP",&ElectronP);
  treeVars->Branch("ElectronTheta",&ElectronTheta);
  treeVars->Branch("ElectronPhi",&ElectronPhi);
  treeVars->Branch("ElectronTime",&ElectronTime);
  treeVars->Branch("ElectronEdep",&ElectronEdep);
  treeVars->Branch("ElectronDeltaE",&ElectronDeltaE);
  treeVars->Branch("ElectronRegion",&ElectronRegion);
  treeVars->Branch("ElectronSector",&ElectronSector);
  //treeVars->Branch("ElectronStatus",&ElectronStatus);
  treeVars->Branch("ElectronDeltaTime",&ElectronDeltaTime); //Peak at 110
  treeVars->Branch("ElectronDeltaTime2",&ElectronDeltaTime2);
  treeVars->Branch("ElectronDeltaTimeFromRecStartTime",&ElectronDeltaTimeFromRecStartTime);
  treeVars->Branch("RFTime",&RFTime);
  treeVars->Branch("EStartTime",&EStartTime); //Peak at 110 and 83 (27diff)
  treeVars->Branch("ElectronStartTimeFromREC",&ElectronStartTimeFromREC); //Peak at 110 and 83 (27diff)
  treeVars->Branch("ElectronShiftedTime",&ElectronShiftedTime);
  treeVars->Branch("NElectrons",&NElectrons);

  ////PROTON
  Double_t ProtonP;
  Double_t ProtonTheta;
  Double_t ProtonPhi;
  Double_t ProtonTime;
  Double_t ProtonEdep;
  Double_t ProtonDeltaE;
  Double_t ProtonRegion;
  Double_t ProtonSector;
  //Double_t ProtonStatus=0;
  Double_t ProtonDeltaTime;
  Double_t ProtonDeltaTimeFromRecStartTime;
  Double_t ProtonShiftedTime;
  Double_t ProtonShiftedTimeByRecStart;
  Double_t NProtons;
  

  treeVars->Branch("ProtonP",&ProtonP);
  treeVars->Branch("ProtonTheta",&ProtonTheta);
  treeVars->Branch("ProtonPhi",&ProtonPhi);
  treeVars->Branch("ProtonTime",&ProtonTime);
  treeVars->Branch("ProtonEdep",&ProtonEdep);
  treeVars->Branch("ProtonDeltaE",&ProtonDeltaE);
  treeVars->Branch("ProtonRegion",&ProtonRegion);
  treeVars->Branch("ProtonSector",&ProtonSector);
  //treeVars->Branch("ProtonStatus",&ProtonStatus);
  treeVars->Branch("ProtonDeltaTime",&ProtonDeltaTime);  //Peak at -24 (7times greater than nearby peaks)
  treeVars->Branch("ProtonDeltaTimeFromRecStartTime",&ProtonDeltaTimeFromRecStartTime);  //Peak at -24 (7times greater than nearby peaks)
  treeVars->Branch("NProtons",&NProtons);



  Double_t MissMass2;
  Double_t MissMass;
  Double_t MissMom;
  Double_t MissEnergy;
  Double_t PionLabPhi;
  Double_t PionLabTheta;
  Double_t ReconProtonPhi;
  Double_t ReconProtonTheta;
  Double_t Coplanarity;
  Double_t ConeAngle;
  Double_t MissMassEP;
  Double_t MissMassEPi;
  Double_t MissMassEP2;
  Double_t MissMassEPi2;
  Double_t InvMass;
  //Double_t UID;
  Double_t ScElE;



  treeVars->Branch("MissMass2",&MissMass2);
  treeVars->Branch("MissMass",&MissMass);
  treeVars->Branch("MissMom",&MissMom);
  treeVars->Branch("MissEnergy",&MissEnergy);
  treeVars->Branch("PionLabPhi",&PionLabPhi);
  treeVars->Branch("PionLabTheta",&PionLabTheta);
  treeVars->Branch("ReconProtonPhi",&ReconProtonPhi);
  treeVars->Branch("ReconProtonTheta",&ReconProtonTheta);
  //  treeVars->Branch("Coplanarity",&Coplanarity);
  //  treeVars->Branch("ConeAngle",&ConeAngle);
  treeVars->Branch("MissMassEP",&MissMassEP);
  treeVars->Branch("MissMassEPi",&MissMassEPi);
  treeVars->Branch("MissMassEP2",&MissMassEP2);
  treeVars->Branch("MissMassEPi2",&MissMassEPi2);
  treeVars->Branch("InvMass",&InvMass);
  treeVars->Branch("Costh",&fCosTh);
  treeVars->Branch("CMPhi",&fPhi);
  treeVars->Branch("W",&W);
  treeVars->Branch("Q2",&Q2);
  treeVars->Branch("t",&t);
  treeVars->Branch("Pol",&Pol);
  //  treeVars->Branch("ScElE",&ScElE);//Scattered Electron Energy
  //  treeVars->Branch("UID",&UID);//Unique ID for (sub)event



  ////////////////////////////////////////////

  TChain chain;
  //  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_*");
//  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_5038*");
  //chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_50*");
  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_5046*");
  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_5051*");
  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_5117*");
  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_5124*");
  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_5424*");
  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_5425*");
  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_5428*");
  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_5429*");
  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_5430*");
//  chain.Add("/project/Gruppo3/fiber7/cmullen/pass1Initial/skim2_*");
  //get the hipo data
  //   reader.open(inputFile.Data());
  auto files=chain.GetListOfFiles();

  //some particles
  auto db=TDatabasePDG::Instance();
  TLorentzVector beam(0,0,10.6,10.6);
  TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector pr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector pip(0,0,0,db->GetParticle(211)->Mass());
  TLorentzVector pim(0,0,0,db->GetParticle(-211)->Mass());


  auto* hmiss=new TH1F("missM","missM",200,-2,3);
  auto* hm2g=new TH1F("m2g","m2g",200,0,1);
  auto* hm2gCut=new TH1F("m2gCut","m2g",200,0,1);
   
  gBenchmark->Start("timer");
  int counter=0;
 
   
  for(Int_t i=0;i<files->GetEntries();i++){
    //create the event reader
    clas12reader c12(files->At(i)->GetTitle());
      
    while(c12.next()==true){

      // get particles by type
      auto electrons=c12.getByID(11); //vector of region_part_ptr
      auto protons=c12.getByID(2212);
      auto pips=c12.getByID(211);
      auto pims=c12.getByID(-211);

      if(electrons.size()==1 && pips.size()==1 && pims.size()==1 && protons.size()>0){
      //if(electrons.size()>0 && gammas.size()>1 && protons.size()>0){
  	for(Int_t j=0; j<protons.size(), j++) {    

//	if(protons.region[i]->getRegion()) {}

	// set the particle momentum
	SetLorentzVector(el,electrons[0]);
	SetLorentzVector(pr,protons[0]);  //Need loop for these
	SetLorentzVector(pip,pips[0]);
	SetLorentzVector(pim,pims[0]);
	




	ElectronP = el.P();
	ElectronTheta = el.Theta();
	ElectronPhi = el.Phi();
	ElectronTime = electrons[0]->getTime();
	ElectronEdep = electrons[0]->getDetEnergy();
	ElectronDeltaE = electrons[0]->getDeltaEnergy();
	ElectronRegion = electrons[0]->getRegion();
	ElectronSector = electrons[0]->getSector();
	//	ElectronStatus = el.P();
	ElectronDeltaTime = DeltaTime(electrons[0]); 
        RFTime = c12.event()->getRFTime() ;
	EStartTime = StartTime(electrons[0],ElectronDeltaTime,RFTime);
        ElectronShiftedTime = ElectronTime - EStartTime;
	ElectronDeltaTime2 = DeltaTime(electrons[0],ElectronShiftedTime); //Recompute Delta Time with timing adjusted by the start time
	//	ElectronDeltaTime = DeltaTime(electrons[0],electrons[0]->getTime());
	ElectronStartTimeFromREC = c12.event()->getStartTime();
	ElectronShiftedTimeByRecStart = ElectronTime - ElectronStartTimeFromREC;
	ElectronDeltaTimeFromRecStartTime = DeltaTime(electrons[0],ElectronShiftedTimeByRecStart); 
	NElectrons = electrons.size();


	ProtonP = pr.P();
	ProtonTheta = pr.Theta();
	ProtonPhi = pr.Phi();
	ProtonTime = protons[0]->getTime();
	ProtonEdep = protons[0]->getDetEnergy();
	ProtonDeltaE = protons[0]->getDeltaEnergy();
	ProtonRegion = protons[0]->getRegion();
	ProtonSector = protons[0]->getSector();
	//	ProtonStatus = pr.P();
	//	ProtonDeltaTime = DeltaTime(protons[0]);
        ProtonShiftedTime = ProtonTime - EStartTime;
       	ProtonDeltaTime = DeltaTime(protons[0],ProtonShiftedTime);
        ProtonShiftedTimeByRecStart = ProtonTime - ElectronStartTimeFromREC;
    	ProtonDeltaTimeFromRecStartTime = DeltaTime(protons[0],ProtonShiftedTimeByRecStart);
        NProtons = protons.size();




	//Discrim Vars
	TLorentzVector miss=beam+target-el-pr-g1-g2;
	MissMass2 = miss.M2();
	MissMass = miss.M();
	MissMom = miss.P();
	MissEnergy = miss.E();

	TLorentzVector pi0=g1+g2;//Fix the mass to pi0 mass at some point
	TLorentzVector ReconProton = beam + target - el - pi0;
	PionLabPhi = pi0.Phi();
	PionLabTheta = pi0.Theta();
	ReconProtonPhi= ReconProton.Phi();
	ReconProtonTheta=ReconProton.Theta();

	Coplanarity = (ROOT::Math::VectorUtil::DeltaPhi(pi0, -(pr) ) )*TMath::RadToDeg();
	ConeAngle = ROOT::Math::VectorUtil::Angle(ReconProton,pr.Vect());  //Meaningless without rotation?

	TLorentzVector missEP = beam + target -el - pr;
	TLorentzVector missEPi = beam + target -el - pi0;
	MissMassEP = missEP.M();
	MissMassEPi = missEPi.M();
	MissMassEP2 = missEP.M2();
	MissMassEPi2 = missEPi.M2();
	TLorentzVector inv = g1 + g2;
	InvMass = inv.M();

	//HSLorentz Vector comes from https://root.cern.ch/root/html/MATH_GENVECTOR_Index.html with definition of
	//typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<Double32_t> > HSLorentzVector;
	//https://root.cern.ch/root/html/ROOT__Math__LorentzVector_-p1PxPyPzE4D_double___.html

	//See https://github.com/dglazier/HASPECT6/blob/hsfarm/HaSpect/THSKinematics.h for more details on below

	fElin=HSLorentzVector(beam.X(),beam.Y(),beam.Z(),beam.T()); //conversion to genvector
	fElsc=HSLorentzVector(el.X(),el.Y(),el.Z(),el.T());
	fGamma=fElin-fElsc;
	fTar=HSLorentzVector(target.X(),target.Y(),target.Z(),target.T());
	fCM=fGamma+fTar;
	fCMBoost=fCM.BoostToCM();
	fMes=HSLorentzVector(pi0.X(),pi0.Y(),pi0.Z(),pi0.T());
	fBar=HSLorentzVector(pr.X(),pr.Y(),pr.Z(),pr.T());

	HSLorentzVector CMBeam=boost(fElin,fCMBoost); 	//Some issue here with unresolved while linking
	HSLorentzVector CMScat=boost(fElsc,fCMBoost);
	HSLorentzVector CMMes=boost(fMes,fCMBoost);
	HSLorentzVector CMGamma=boost(fGamma,fCMBoost);

	HSMomentum zV=CMGamma.Vect().Unit();
	HSMomentum yV=CMGamma.Vect().Cross(CMBeam.Vect()).Unit();
	HSMomentum xV=yV.Cross(zV).Unit();

	HSMomentum angles(CMMes.Vect().Dot(xV),CMMes.Vect().Dot(yV),CMMes.Vect().Dot(zV));
	fCosTh=TMath::Cos(angles.Theta());  //CosThetaVector(&angles);
	fPhi=angles.Phi();
	W =  fCM.M();
	Q2 = -fGamma.M2();
	t = (fGamma-fMes).M2();
	Pol = 1./(1+2*(1+fGamma.E()*fGamma.E()/Q2)*TMath::Tan(fElsc.Theta()/2)*TMath::Tan(fElsc.Theta()/2));

                               
	///////////////////////////////////////////////////////////////////////////////////////////////

	treeVars->Fill();

	//	TLorentzVector miss=beam+target-el-pr-g1-g2;
	hmiss->Fill(miss.M2());
	//	TLorentzVector pi0 = g1+g2;
	hm2g->Fill(pi0.M());
	if(TMath::Abs(miss.M2())<0.5)hm2gCut->Fill(pi0.M());

	//could also get particle time etc. here too
	//Double_t eTime=electrons[0]->sci(FTOF1A)->getTime();
      }
    
       
      counter++;
    }
  }
  gBenchmark->Stop("timer");
  gBenchmark->Print("timer");
  TCanvas* can=new TCanvas();
  can->Divide(2,1);
  can->cd(1);
  hmiss->DrawCopy();
  can->cd(2);
  hm2g->DrawCopy();
  hm2gCut->SetLineColor(2);
  hm2gCut->DrawCopy("same");
  
  hmiss->Write();
  hm2g->Write();
  hm2gCut->Write();

  treeVars->Write();
  treeVars->Reset();
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";

  outputFile->Close();

}
