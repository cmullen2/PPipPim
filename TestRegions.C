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

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	     rp->par()->getPz(),p4.M());

}           



void TestRegions(){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();


  /////////////////////////////////////
  TFile *outputFile = new TFile("/project/Gruppo3/fiber7/cmullen/ServiceWork/output/FirstRun.root","recreate");
  /////////////////////////////////////

  /////////////////////////////////////////////
  TTree *treeVars = new TTree("ParticleVars","SomeTree");

  ////////////////////////////////////////////

  TChain chain;
  //  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_*");
//  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_5038*");
  //chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_50*");
//  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_5046*");
//  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_5051*");
  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_5117*");
//  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_5124*");
//  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_5424*");
//  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_5425*");
//  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_5428*");
//  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_5429*");
//  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_5430*");
//  chain.Add("/project/Gruppo3/fiber7/cmullen/pass1Initial/skim2_*");
//  chain.Add("/project/Gruppo3/fiber7/cmullen/ServiceWork/pass0v1_1.1.9/rec*");
  //get the hipo data
  //   reader.open(inputFile.Data());
  auto files=chain.GetListOfFiles();

   
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
      if(electrons.size()>0 && protons.size()>0 && pips.size()>0 &&pims.size()>0){
      //if(electrons.size()>0 && gammas.size()>1 && protons.size()>0){
  	for(Int_t j=0; j<electrons.size(); j++) {
	if(protons[0]->getRegion()==CD){    
	if(pips[0]->getRegion()<3000){    
	if(pims[0]->getRegion()<3000){    

//	if(protons[j]->getRegion()>3999 ) {
//	if(protons[j]->getRegion()==CD ) {  //gives region of 3000 which should be FD+FT
//	if(protons[j]->getRegion()==FD ) {  //gives region of 2000 which should be FD so fine (also gives majority of events)
	if(electrons[j]->getRegion()==FT ) {  
//	if(protons[j]->getStatus()>3999 ) {
//	if(electrons[0]->getRegion()<3000 && pips[0]->getRegion()<3000 && pims[0]->getRegion()<3000){

	
cout <<"      "<< electrons[j]->getRegion()<<endl;
cout <<"asdasdasdas" << endl;
	///////////////////////////////////////////////////////////////////////////////////////////////

	treeVars->Fill();
}
}
}
				}

			}

		}

  //    }
    
       
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
