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


void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	     rp->par()->getPz(),p4.M());

}           



void CLAS12PPipPimReader(){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();


  /////////////////////////////////////
  TFile *outputFile = new TFile("/project/Gruppo3/fiber7/cmullen/ServiceWork/output/First7RecFilesWithPionsInFDorFT.root","recreate");
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
  Double_t ElectronStartTimeFromREC;
  Double_t RFTime;
  Double_t ElectronShiftedTime;
  Double_t NElectrons;

  treeVars->Branch("ElectronP",&ElectronP);
  treeVars->Branch("ElectronTheta",&ElectronTheta);
  treeVars->Branch("ElectronPhi",&ElectronPhi);
  treeVars->Branch("ElectronTime",&ElectronTime);
  treeVars->Branch("ElectronEdep",&ElectronEdep);
  treeVars->Branch("ElectronDeltaE",&ElectronDeltaE);
  treeVars->Branch("ElectronRegion",&ElectronRegion);
  treeVars->Branch("ElectronSector",&ElectronSector);
  //treeVars->Branch("ElectronStatus",&ElectronStatus);
  treeVars->Branch("RFTime",&RFTime);
  treeVars->Branch("ElectronStartTimeFromREC",&ElectronStartTimeFromREC); 
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
  treeVars->Branch("NProtons",&NProtons);




  ////Pip
  Double_t PipP;
  Double_t PipTheta;
  Double_t PipPhi;
  Double_t PipTime;
  Double_t PipEdep;
  Double_t PipDeltaE;
  Double_t PipRegion;
  Double_t PipSector;
  //Double_t PipStatus=0;
  Double_t NPips;
  

  treeVars->Branch("PipP",&PipP);
  treeVars->Branch("PipTheta",&PipTheta);
  treeVars->Branch("PipPhi",&PipPhi);
  treeVars->Branch("PipTime",&PipTime);
  treeVars->Branch("PipEdep",&PipEdep);
  treeVars->Branch("PipDeltaE",&PipDeltaE);
  treeVars->Branch("PipRegion",&PipRegion);
  treeVars->Branch("PipSector",&PipSector);
  //treeVars->Branch("PipStatus",&PipStatus);
  treeVars->Branch("NPips",&NPips);

  ////Pim
  Double_t PimP;
  Double_t PimTheta;
  Double_t PimPhi;
  Double_t PimTime;
  Double_t PimEdep;
  Double_t PimDeltaE;
  Double_t PimRegion;
  Double_t PimSector;
  //Double_t PimStatus=0;
  Double_t NPims;
  

  treeVars->Branch("PimP",&PimP);
  treeVars->Branch("PimTheta",&PimTheta);
  treeVars->Branch("PimPhi",&PimPhi);
  treeVars->Branch("PimTime",&PimTime);
  treeVars->Branch("PimEdep",&PimEdep);
  treeVars->Branch("PimDeltaE",&PimDeltaE);
  treeVars->Branch("PimRegion",&PimRegion);
  treeVars->Branch("PimSector",&PimSector);
  //treeVars->Branch("PimStatus",&PimStatus);
  treeVars->Branch("NPims",&NPims);


  Double_t MissMom;
  Double_t MissEnergy;


  Double_t ReconProtonPhi;
  Double_t ReconProtonMom;
  Double_t ReconProtonTheta;
  Double_t Coplanarity;
  Double_t ConeAngle;
  Double_t ScElE;

  Double_t  MissMass;
  Double_t  MissMass2;
  Double_t  MissMassEPipPimX;
  Double_t  MissMassEPipPimX2;
  Double_t  MissMassEPrX;
  Double_t  MissMassEPrX2;
  Double_t  MissMassEPimX ;
  Double_t  MissMassEPimX2;
  Double_t  MissMassEPipX ;
  Double_t  MissMassEPipX2;
  Double_t  MissMassEPrPipX ;
  Double_t  MissMassEPrPipX2;
  Double_t  MissMassEPrPimX;
  Double_t  MissMassEPrPimX2;
 
  Double_t  PhiDifferenceProton;
  Double_t  ThetaDifferenceProton;
  Double_t  MomentumDifferenceProton;


  treeVars->Branch("MissMass2",&MissMass2);
  treeVars->Branch("MissMass",&MissMass);
  treeVars->Branch("MissMassEPipPimX",&MissMassEPipPimX);
  treeVars->Branch("MissMassEPipPimX2",&MissMassEPipPimX2);
  treeVars->Branch("MissMassEPrX",&MissMassEPrX);
  treeVars->Branch("MissMassEPrX2",&MissMassEPrX2);
  treeVars->Branch("MissMassEPimX",&MissMassEPimX);
  treeVars->Branch("MissMassEPimX2",&MissMassEPimX2);
  treeVars->Branch("MissMassEPipX",&MissMassEPipX);
  treeVars->Branch("MissMassEPipX2",&MissMassEPipX2);
  treeVars->Branch("MissMassEPrPipX",&MissMassEPrPipX);
  treeVars->Branch("MissMassEPrPipX2",&MissMassEPrPipX2);
  treeVars->Branch("MissMassEPrPimX2",&MissMassEPrPimX);
  treeVars->Branch("MissMassEPrPimX2",&MissMassEPrPimX2);
  treeVars->Branch("MissMom",&MissMom);
  treeVars->Branch("MissEnergy",&MissEnergy);

  treeVars->Branch("ReconProtonPhi",&ReconProtonPhi);
  treeVars->Branch("ReconProtonTheta",&ReconProtonTheta);
  treeVars->Branch("ReconProtonMom",&ReconProtonMom);
  treeVars->Branch("PhiDifferenceProton",&PhiDifferenceProton);
  treeVars->Branch("ThetaDifferenceProton",&ThetaDifferenceProton);
  treeVars->Branch("MomentumDifferenceProton",&MomentumDifferenceProton);
  //  treeVars->Branch("ScElE",&ScElE);//Scattered Electron Energy
  //  treeVars->Branch("UID",&UID);//Unique ID for (sub)event

  ////////////////////////////////////////////

  TChain chain;
  chain.Add("/project/Gruppo3/fiber7/cmullen/ServiceWork/pass0v1_1.1.9/rec*");
  //  chain.Add("/project/Gruppo3/fiber7/cmullen/ServiceWork/pass0v1_1.1.9/skim*");
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

      // if(electrons.size()==1 && pips.size()==1 && pims.size()==1 && protons.size()>0){
      if(electrons.size()>0 && pips.size()>0 && pims.size()>0&& protons.size()>0){
	//if(electrons.size()>0 && gammas.size()>1 && protons.size()>0){
  	for(Int_t j=0; j<protons.size(); j++) {    
	  if(protons[j]->getRegion()==CD ) {  
	    if(electrons[0]->getRegion()==FT ) {  
	      if(pips[0]->getRegion()==FD ) {  
		if(pims[0]->getRegion()==FD ) {  

		  // set the particle momentum
		  SetLorentzVector(el,electrons[0]);
		  SetLorentzVector(pr,protons[j]);  //Need loop for these
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
		  RFTime = c12.event()->getRFTime() ;
		  ElectronStartTimeFromREC = c12.event()->getStartTime();
		  NElectrons = electrons.size();


		  ProtonP = pr.P();
		  ProtonTheta = pr.Theta();
		  ProtonPhi = pr.Phi();
		  ProtonTime = protons[j]->getTime();
		  ProtonEdep = protons[j]->getDetEnergy();
		  ProtonDeltaE = protons[j]->getDeltaEnergy();
		  ProtonRegion = protons[j]->getRegion();
		  ProtonSector = protons[j]->getSector();
		  //	ProtonStatus = pr.P();
		  NProtons = protons.size();

		  PipP = pip.P();
		  PipTheta = pip.Theta();
		  PipPhi = pip.Phi();
		  PipTime = pips[0]->getTime();
		  PipEdep = pips[0]->getDetEnergy();
		  PipDeltaE = pips[0]->getDeltaEnergy();
		  PipRegion = pips[0]->getRegion();
		  PipSector = pips[0]->getSector();
		  //	ProtonStatus = pr.P();
		  NPips = pips.size();


		  PimP = pim.P();
		  PimTheta = pim.Theta();
		  PimPhi = pim.Phi();
		  PimTime = pims[0]->getTime();
		  PimEdep = pims[0]->getDetEnergy();
		  PimDeltaE = pims[0]->getDeltaEnergy();
		  PimRegion = pims[0]->getRegion();
		  PimSector = pims[0]->getSector();
		  //	ProtonStatus = pr.P();
		  NPims = pims.size();


		  TLorentzVector Tmiss = beam + target - el - pr - pip -pim;
		  TLorentzVector TReconProton = beam + target -el -pip -pim;
		  TLorentzVector TmissEPrX = beam + target -el -pr;             
		  TLorentzVector TmissEPimX = beam + target - el -pim; 
		  TLorentzVector TmissEPipX = beam+target-el-pip;
		  TLorentzVector TmissEPrPipX = beam+ target -el -pr - pip;
		  TLorentzVector TmissEPrPimX = beam +target - el -pr - pim;


		  MissMom = Tmiss.P();
		  MissEnergy = Tmiss.E();
		  MissMass =  Tmiss.M();
		  MissMass2 = Tmiss.M2();
		  MissMassEPipPimX = TReconProton.M();
		  MissMassEPipPimX2 = TReconProton.M2();
		  MissMassEPrX = TmissEPrX.M();
		  MissMassEPrX2 = TmissEPrX.M2() ;
		  MissMassEPimX = TmissEPimX.M();
		  MissMassEPimX2 = TmissEPimX.M2();
		  MissMassEPipX = TmissEPipX.M();
		  MissMassEPipX2 = TmissEPipX.M2();
		  MissMassEPrPipX = TmissEPrPipX.M();
		  MissMassEPrPipX2 = TmissEPrPipX.M2();
		  MissMassEPrPimX = TmissEPrPimX.M();
		  MissMassEPrPimX2 = TmissEPrPimX.M2();

		  ReconProtonPhi = TReconProton.Phi();
		  ReconProtonTheta = TReconProton.Theta();
		  ReconProtonMom = TReconProton.P();

		  ThetaDifferenceProton = ReconProtonTheta - ProtonTheta;
		  PhiDifferenceProton = ReconProtonPhi - ProtonPhi;
		  MomentumDifferenceProton = ReconProtonMom - ProtonP;


		  ///////////////////////////////////////////////////////////////////////////////////////////////

		  treeVars->Fill();
		}
	      }
	    }

	  }

	}

      }
    
       
      counter++;
    }
  }
  gBenchmark->Stop("timer");
  gBenchmark->Print("timer");

  treeVars->Write();
  treeVars->Reset();
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";

  outputFile->Close();

}
