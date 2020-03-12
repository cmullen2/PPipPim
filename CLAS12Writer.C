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
#include "clas12reader.h"
#include "clas12writer.h"

using namespace clas12;


void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	      rp->par()->getPz(),p4.M());

}

void CLAS12Pi0Writer(){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();


  /////////////////////////////////////
  std::string outputFile;
 outputFile= "FirstPi0.root"; //Want this to be a root file so look at other examples,Probs need to do what haspect does or other example, Test first when chain is set
  /////////////////////////////////////



  HipoChain chain;
  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_*");

  //initialising clas12writer with path to output file
  clas12writer c12writer(outputFile.c_str());

  //can as writer not to write certain banks
  //c12writer.skipBank("REC::Cherenkov");
  //c12writer.skipBank("REC::Scintillator");

  cout<<"Analysing hipo file "<<inFile<<endl;

  //some particles
  auto db=TDatabasePDG::Instance();
  TLorentzVector beam(0,0,10.6,10.6);
  TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector pr(0,0,0,db->GetParticle(2212)->Mass());
   
  gBenchmark->Start("timer");
 
  int counter = 0;
  int writeCounter = 0;
   //Do I want the chain here or in a separate file like the HAspect code

  //create the event reader
  clas12reader c12(inFile.c_str());

  //assign a reader to the writer
  c12writer.assignReader(c12);
     
      
  while(c12.next()==true){

    // get particles by type
    auto electrons=c12.getByID(11);
    auto protons=c12.getByID(2212);
    auto pips=c12.getByID(211);
    auto pims=c12.getByID(-211);
       
//Need electron pi- pip+ in FD or FT and any number of protons in the central d.
//Need to access things like chi2/ndf and know the resolutions of cd. Also which data to use


    if(electrons.size()==1 && protons.size()>0 &&  pips.size()==1 &&pims.size() == 1){
//	for(const int& i : protons){
//		cout << i <<endl;
//}
	for (Int_t i=0;)

       //Loop for proton size and check each are in the central, calculate differences between each of resolutions
      // set the particle momentum
      SetLorentzVector(el,electrons[0]);
      SetLorentzVector(pr,protons[0]);
      SetLorentzVector(pip,pips[0]);
      SetLorentzVector(pim,pims[0]);
	
      TLorentzVector miss=beam+target-el-pr-pip-pim;
      if(TMath::Abs(miss.M2())<0.5){
	//write out an event
	c12writer.writeEvent(); 
	writeCounter++;
      }
      counter++;
    }
  }
  

  //close writer
  c12writer.closeWriter();

  gBenchmark->Stop("timer");
  gBenchmark->Print("timer");
  
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count()<< " read events = "<<counter<<" wrote events = "<<writeCounter<<" s\n";

}
