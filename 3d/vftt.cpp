/* Kernel Independent Fast Multipole Method

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.  */

/*!
 This file provides an example test for the sequential fmm3d.
 To generate the executable, type "make tt" and then run tt options_file
*/

#include "Vfmm3d_FMM.hpp"
#include <omp.h>

using namespace std;

/*! parse the options file provided as input to tt
 * OptionsCreate takes a string and a character map and maps options file
 * names to their values to be parsed later
 */
int optionsCreate(const char* optionsFile, map<string,string>& options)
{
  options.clear();
  ifstream fileIn(optionsFile);
  if(fileIn.good()==false) {
	 cerr<<"wrong option file"<<endl;	 exit(1);
  }
  string name;  fileIn>>name;
  while(fileIn.good()) {
	 char content[100];	 fileIn.getline(content, 99);
	 options[name] = string(content);
	 fileIn>>name;
  }
  fileIn.close();
  return 0;
}

/*! Main function operates in several steps:
 * First, the options file as provided to optionsCreate is parsed.
 * From this, we get the number of sources and target points (numSrc, numTrg),
 * the kernel type (kt - see kernel3d.hpp and kernel3d.cpp), etc.  Random locations
 * are generated for source and target locations.  Random sources densities are allocated based
 * on the degrees of freedom of the specified kernel, and target value space is allocated as well.
 * Then, space for the fmm is allocated and variables are set.  Fmm is then this run and checked/clocked.
 */
int main(int argc, char** argv)
{
  srand48( (long)time(NULL) );  //srand48( 0 );
  iA(argc==2);
  map<string,string> optionsMap;
  optionsCreate(argv[1], optionsMap);
  
  ///1. allocate random data
  map<string,string>::iterator mapIdx;
  mapIdx = optionsMap.find("-numSrc"); assert(mapIdx!=optionsMap.end());
  /*! numSrc = Number of source points */
  int numSrc;  { istringstream ss((*mapIdx).second);  ss >> numSrc; }
  mapIdx = optionsMap.find("-numTrg"); assert(mapIdx!=optionsMap.end());
  /*! numTrg = Number of target points */
  int numTrg;  { istringstream ss((*mapIdx).second);  ss >> numTrg; }
  mapIdx = optionsMap.find("-kt"); assert(mapIdx!=optionsMap.end());
  /* ht = Kernel Type.  See kernel3d.hpp */
  int kernelType;  { istringstream ss((*mapIdx).second);  ss>>kernelType; }
    
  vector<double> tmp(2);	 tmp[0] = 1;	 tmp[1] = 0.25; //coefs in the kernel, work for all examples
  /*! Declare kernel as knl of type kt and coefficients temp */
  Kernel3d knl(kernelType, tmp);

  /*! Allocate random data for source positions (srcPos) and target positions (trgPos) */
  DblNumMat srcPos(3, numSrc);
  for(int i=0; i<numSrc; i++) {
    if (i < 0) { i++; }
    Point3 pnt(2.0*drand48()-1.0, 2.0*drand48()-1.0, 2.0*drand48()-1.0);
    //if (pnt.l2() >= 0.4 && pnt.l2() <= 0.9){
    //if (pnt.l2() <= 0.9){
	 if (1){
      srcPos(0,i) = pnt(0);
      srcPos(1,i) = pnt(1);
      srcPos(2,i) = pnt(2);
    }
    else {
      i--;
    }
  }
  DblNumMat trgPos(3, numTrg);
  for(int i=0; i<numTrg; i++) {
    if (i < 0) { i++; }
    Point3 pnt(2.0*drand48()-1.0, 2.0*drand48()-1.0, 2.0*drand48()-1.0);
    //if (pnt.l2() >= 0.4 && pnt.l2() <= 0.9){
    //if (pnt.l2() <= 0.9){
	 if (1){
      trgPos(0,i) = pnt(0);
      trgPos(1,i) = pnt(1);
      trgPos(2,i) = pnt(2);
    }
    else {
      i--;
    }
  }

  /*! srcDOF = source degree of freedom.  See kernel3d.hpp */
  int srcDOF = knl.srcDOF();
  /*! trgDOF = target degree of freedom.  See kernel3d.hpp */
  int trgDOF = knl.trgDOF();
  /*! srcDen = source density values of size sDOF*numSrc */
  DblNumVec srcDen(srcDOF * numSrc);
  for(int i=0; i<numSrc; i++) {
	 for(int d=0; d<srcDOF; d++)
	   srcDen(d + i*srcDOF) = drand48(); //(2.0*drand48()-1.0);
  }
  /* trgVal = target values of size trgDOF *numTrg */
  DblNumVec trgVal(trgDOF * numTrg);
  
  ///2. allocate fmm 
  clock_t clockZero, clockOne;
  
  FMM3d<Node>* _fmm = new FMM3d<Node>("fmm3d_");
  
  _fmm->srcPos()=&srcPos;  _fmm->srcNor()=&srcPos;
  _fmm->srcDen()=&srcDen;
  _fmm->trgPos()=&trgPos;
  _fmm->trgVal()=&trgVal;
  _fmm->center() = Point3(0,0,0); // CENTER OF THE TOPLEVEL BOX
  _fmm->rootLevel() = 0;         // 2^(-rootlvl) is the RADIUS OF THE TOPLEVEL BOX
  _fmm->knl() = knl;
  
  mapIdx = optionsMap.find("-vkt"); assert(mapIdx!=optionsMap.end());
  /* kt = Kernel Type.  See kernel3d.hpp */
  int vKernelType;  { istringstream ss((*mapIdx).second);  ss>>vKernelType; }
  mapIdx = optionsMap.find("-rootlevel"); assert(mapIdx!=optionsMap.end());
  int rootlevel; { istringstream ss((*mapIdx).second);  ss>>rootlevel; }
  cerr << "ROOTLEVEL = " << rootlevel << endl;
  mapIdx = optionsMap.find("-vfmm3d_lambda"); assert(mapIdx!=optionsMap.end());
  /* lambda.  ********/
  //#warning lambda needs to be changed to alpha since lambda = sqrt(alpha) or sqrt(alpha/mu)
  double lambda;  { istringstream ss((*mapIdx).second);  ss>>lambda; }
  
  mapIdx = optionsMap.find("-vfmm3d_rval"); assert(mapIdx!=optionsMap.end());
  double rval;  { istringstream ss((*mapIdx).second);  ss>>rval; }
  
  vector<double> vTmp(2);	 vTmp[0] = 1.0;	 vTmp[1] = 0.25; //coefs in the kernel, work for all examples
  vTmp[1] = lambda;
  vTmp[0] = rval;

  mapIdx = optionsMap.find("-vfmm3d_eqnType"); assert(mapIdx!=optionsMap.end());
  int eqntype;  { istringstream ss((*mapIdx).second);  ss>>eqntype; }

  /*! Declare kernel as knl of type kt and coefficients temp */
  Kernel3d vKnl(vKernelType, vTmp);
  Exsol3d* exsol3d = new Exsol3d(vKnl.kernelType(), vKnl.coefs(), eqntype);

  cout << "kt = " << vKernelType << endl;

  VFMM3d_FMM* vfmm = new VFMM3d_FMM("vfmm3d_");
  vfmm->fmm() = _fmm;
  vfmm->center() = Point3(0.0,0.0,0.0); // CENTER OF THE TOPLEVEL BOX
  vfmm->rootLevel() = rootlevel;         // 2^(-rootlvl) is the RADIUS OF THE TOPLEVEL BOX
  vfmm->knl() =vKnl;
  vfmm->srcPos() = vfmm->fmm()->srcPos();
  vfmm->trgPos() = vfmm->fmm()->srcPos(); //Both set to srcPos
  vfmm->exsol3d() = exsol3d;

  iC( vfmm->fmm()->setFromOptions(optionsMap) );
  iC( vfmm->setFromOptions(optionsMap) );
  iC( vfmm->setup() );
  iC( vfmm->evaluate(srcDen, trgVal) );
  double rerr;
  iC(_fmm->check(srcDen, trgVal, 20, rerr));
  
  delete vfmm;
  
  return 0;
}
