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

#include "Vfmm3d.hpp"
#include "Efmm3d.hpp"

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

int main(int argc, char** argv)
{
  srand48( (long)time(NULL) );  //srand48( 0 );
  iA(argc==2);
  map<string,string> optionsMap;
  optionsCreate(argv[1], optionsMap);
  
  ///1. allocate random data
  map<string,string>::iterator mapIdx;

  mapIdx = optionsMap.find("-vkt"); assert(mapIdx!=optionsMap.end());
  /* kt = Kernel Type.  See kernel3d.hpp */
  int kernelType;  { istringstream ss((*mapIdx).second);  ss>>kernelType; }

  mapIdx = optionsMap.find("-rootlevel"); assert(mapIdx!=optionsMap.end());
  int rootlevel; { istringstream ss((*mapIdx).second);  ss>>rootlevel; }
  cerr << "ROOTLEVEL = " << rootlevel << endl;
    
  mapIdx = optionsMap.find("-vfmm3d_lambda"); assert(mapIdx!=optionsMap.end());
  /* lambda.  ********/
  //#warning lambda needs to be changed to alpha since lambda = sqrt(alpha) or sqrt(alpha/mu)
  double lambda;  { istringstream ss((*mapIdx).second);  ss>>lambda; }
  
  mapIdx = optionsMap.find("-vfmm3d_rval"); assert(mapIdx!=optionsMap.end());
  double rval;  { istringstream ss((*mapIdx).second);  ss>>rval; }
  
  vector<double> tmp(2);	 tmp[0] = 1.0;	 tmp[1] = 0.25; //coefs in the kernel, work for all examples
  tmp[1] = lambda;
  tmp[0] = rval;

  mapIdx = optionsMap.find("-vfmm3d_eqnType"); assert(mapIdx!=optionsMap.end());
  int eqntype;  { istringstream ss((*mapIdx).second);  ss>>eqntype; }

  /*! Declare kernel as knl of type kt and coefficients temp */
  Kernel3d knl(kernelType, tmp);
  cout << kernelType << endl;
  
  VFMM3d<VolNode>* vfmm = new VFMM3d<VolNode>("vfmm3d_");
  Exsol3d* exsol3d = new Exsol3d(knl.kernelType(), knl.coefs(), eqntype);

  /*! srcDOF = source degree of freedom.  See kernel3d.hpp */
  int srcDOF = knl.srcDOF();
  /*! trgDOF = target degree of freedom.  See kernel3d.hpp */
  int trgDOF = knl.trgDOF();
  /*! srcDen = source density values of size sDOF*numSrc */
  /*! Allocate random data for source positions (srcPos) and target positions (trgPos) */

  DblNumMat srcPos(3, 1);
  DblNumMat trgPos(3, 1);
  DblNumVec srcDen(srcDOF);
  DblNumVec trgVal(trgDOF);

  /* trgVal = target values of size trgDOF *numTrg */
  cout << trgVal.m() << endl;
  
  ///2. allocate fmm 
  clock_t clockZero, clockOne;	
 
  vfmm->srcPos()=&srcPos;  vfmm->srcNor()=&srcPos;
  vfmm->trgPos()=&trgPos;
  vfmm->srcDen()=&srcDen;
  vfmm->trgVal()=&trgVal;
  
  vfmm->center() = Point3(0.0,0.0,0.0); // CENTER OF THE TOPLEVEL BOX
  vfmm->rootLevel() = rootlevel;         // 2^(-rootlvl) is the RADIUS OF THE TOPLEVEL BOX
  vfmm->knl() = knl;

  vfmm->exsol3d() = exsol3d;

  clockZero = clock();
  iC( vfmm->setFromOptions(optionsMap) );
  iC( vfmm->setup() );
  
  clockOne = clock();  cout<<"fmm setup used "<<double(clockOne-clockZero)/CLOCKS_PER_SEC<<"secs "<<endl;

  
  ///3. run fmm
  for(int i=0; i<1; i++) {
	 clockZero = clock();
	 iC( vfmm->evaluate(srcDen, trgVal) );
	 clockOne = clock();  cout<<"fmm eval used "<<double(clockOne-clockZero)/CLOCKS_PER_SEC<<"secs "<<endl;
  }
  
  ///4. check
  clockZero = clock();
  double relativeError, infinityError;
  iC( vfmm->vcheck(relativeError, infinityError) );
  cout << "relative error: " << relativeError << endl; 
  //iC( vfmm->grdCheck(20, relativeError) );
  //cout << "relative error: " << relativeError << endl << endl << endl;
  //cout << "relative error: " << relativeError << endl;
  clockOne = clock();  cout<<"fmm check used "<<double(clockOne-clockZero)/CLOCKS_PER_SEC<<"sec "<<endl;
  
  return 0;
}
