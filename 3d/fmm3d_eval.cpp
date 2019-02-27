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
#include "fmm3d.hpp"
#include "common/vecmatop.hpp"

#include <omp.h>

using std::cerr;
using std::endl;

// ----------------------------------------------------------------------
template <class N>
int FMM3d<N>::evaluate(const DblNumVec& srcDen, DblNumVec& trgVal)
{
  _matmgnt->report();
  //-----------------------------------
  cerr << srcDen.m() << " " << (*_srcPos).n() << endl;
  iA(srcDen.m()==srcDOF()*(*_srcPos).n());
  cerr << trgVal.m() << " " << (*_trgPos).n() << endl;
  iA(trgVal.m()==trgDOF()*(*_trgPos).n());
  
  //cerr<<"fmm src and trg numbers "<<pglbnum(_srcPos)<<" "<<pglbnum(_trgPos)<<endl;
  int srcDOF = this->srcDOF();
  int trgDOF = this->trgDOF();
  
  //1. zero out Vecs
  setvalue(trgVal, 0.0);
  
  setvalue(_srcExaDen, 0.0);
  setvalue(_srcUpwEquDen, 0.0);
  
  setvalue(_trgExaVal, 0.0);
  setvalue(_trgDwnChkVal, 0.0);
  //CLOCKING;
  double fmmevaltime = 0.0;

  vector<int> ordVec; iC( _let->upwOrderCollect(ordVec) ); //BOTTOM UP

  //2. for cbtr, load ExaDen

  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if(_let->tag(gNodeIdx) & LET_SRCNODE) {
		if(_let->terminal(gNodeIdx)==true) {
		  DblNumVec srcExaDen(this->srcExaDen(gNodeIdx));
		  vector<int>& curVecIdxs = _let->node(gNodeIdx).srcOwnVecIdxs();
		  for(int k=0; k<curVecIdxs.size(); k++) {
			 int poff = curVecIdxs[k];
			 for(int d=0; d<srcDOF; d++) {
				srcExaDen(k*srcDOF+d) = srcDen(poff*srcDOF+d);
			 }
		  }
		}
	 }
  }

  map<int, vector<int> > lvlOrdVec; iC( let()->revLvlOrderCollect(lvlOrdVec) ); //BOTTOM UP
  vector<bool> lvsAtLev;
  for (int j = 0; j < lvlOrdVec.size(); j++){
	 lvsAtLev.push_back(false);
	 vector<int>& thisLevBoxes = lvlOrdVec[j];
	 for (int i = 0; i < thisLevBoxes.size(); i++){
		int gNodeIdx = thisLevBoxes[i];
		if (let()->terminal(gNodeIdx) == true) { lvsAtLev[j] = true; }
	 }
  }
  
  //3. up computation
  
  //Pre-build matrices for openmp
  int levnum = (_matmgnt->hom() ? 1 : lvlOrdVec.size());
  for (int j = 0; j < levnum; j++){
	 vector<int>& thisLevBoxes = lvlOrdVec[j];
	 DblNumVec chkVal(_matmgnt->plnDatSze(UC));
	 DblNumVec tmpVec(datSze(UE));
	 for(int a=0; a<2; a++) { for(int b=0; b<2; b++) { for(int c=0; c<2; c++) {
			 Index3 idx(a,b,c);
			 iC( _matmgnt->UpwEqu2UpwChk_dgemv((j+1)+_rootLevel, idx, tmpVec, chkVal, 2.0) );
		  } } }
	 iC( _matmgnt->UpwChk2UpwEqu_dgemv(j+_rootLevel, chkVal, tmpVec) );
  }

  double startTime = omp_get_wtime();
  for (int j = lvlOrdVec.size(); j >= 0; j--){
    vector<int>& thisLevBoxes = lvlOrdVec[j];
#pragma omp parallel for
    for(int i=0; i<thisLevBoxes.size(); i++) {
      int gNodeIdx = thisLevBoxes[i];
      if(_let->tag(gNodeIdx) & LET_SRCNODE) {		//GNTra gnt = _let->gNodeIdx2gnt(gNodeIdx);
	//if(_let->depth(gNodeIdx)>=2) {
	if(_let->depth(gNodeIdx)>=0) {
	  DblNumVec chkVal(_matmgnt->plnDatSze(UC));
	  setvalue(chkVal, 0.0);
	  DblNumVec srcUpwEquDengNodeIdx(srcUpwEquDen(gNodeIdx));
	  if(_let->terminal(gNodeIdx)==true) {
	    //S2M - Source -> Multipole Exapnsion
	    iC( SrcEqu2UpwChk_dgemv(srcExaPos(gNodeIdx), srcExaNor(gNodeIdx), _let->center(gNodeIdx), _let->radius(gNodeIdx), srcExaDen(gNodeIdx), chkVal) );
	    cerr << chkVal << endl;
	  } else {
	    //M2M - Multipole -> Multipole
	    for(int a=0; a<2; a++) {
	      for(int b=0; b<2; b++) {
		for(int c=0; c<2; c++) {
		  Index3 idx(a,b,c);
		  int chi = _let->child(gNodeIdx, idx);
		  if(_let->tag(chi) & LET_SRCNODE) {
		    iC( _matmgnt->UpwEqu2UpwChk_dgemv(_let->depth(chi)+_rootLevel, idx, srcUpwEquDen(chi), chkVal) );
		  }
		}
	      }
	    }
	  }
	  //M2M - Multipole -> Multipole
	  iC( _matmgnt->UpwChk2UpwEqu_dgemv(_let->depth(gNodeIdx)+_rootLevel, chkVal, srcUpwEquDengNodeIdx) );
	  cerr << gNodeIdx << " " << _let->terminal(gNodeIdx) << " " << srcUpwEquDengNodeIdx << endl;
	}
      }
    }
  }
  double endTime = omp_get_wtime();
  fmmevaltime += (endTime - startTime);
  
  std::cout << "S2M used " << (endTime - startTime) << std::endl;

  vector<int> trgTrmVec;
  for (int i = 0; i < ordVec.size(); i++){
	 int gni = ordVec[i];
	 if (let()->terminal(gni) && let()->tag(gni) & LET_TRGNODE) { //evaluator
		trgTrmVec.push_back(gni);
	 }
  }

  ordVec.clear();  iC( _let->dwnOrderCollect(ordVec) );
  //U - list contribution calculation
  startTime = omp_get_wtime();
#pragma omp parallel for
  for(int i=0; i<trgTrmVec.size(); i++) {
	 int gNodeIdx = trgTrmVec[i];
	 if(_let->tag(gNodeIdx) & LET_TRGNODE) { //evaluator
		if( _let->terminal(gNodeIdx)==true ) { //terminal	
		  N& curNode = _let->node(gNodeIdx);
		  DblNumVec trgExaValgNodeIdx(trgExaVal(gNodeIdx));
		  DblNumMat trgExaPosgNodeIdx(trgExaPos(gNodeIdx));
		  for(vector<int>::iterator vi=curNode.Unodes().begin(); vi!=curNode.Unodes().end(); vi++) {
		    //S2T - source -> target
		    iC( SrcEqu2TrgChk_dgemv(srcExaPos(*vi), srcExaNor(*vi), trgExaPosgNodeIdx, srcExaDen(*vi), trgExaValgNodeIdx) );
		  }
		}
	 }
  }
  endTime = omp_get_wtime();
  fmmevaltime += (endTime - startTime);
  std::cout << "NEAR used " << (endTime - startTime) << endl;

  for (int j = 0; j < levnum; j++){
	 DblNumVec effVal(_matmgnt->effDatSze(DC));
	 DblNumVec effDen(_matmgnt->effDatSze(UE));
	 
	 //M2L - multipole -> local
	 for(int a=-3; a<=3; a++) { for(int b=-3; b<=3; b++) { for(int c=-3; c<=3; c++) {
			 Index3 idx(a,b,c);
			 if (idx.linfty() > 1) {
				iC( _matmgnt->UpwEqu2DwnChk_dgemv(j+_rootLevel, idx, effDen, effVal, 2.0) );
			 } } } }
  }
  
  //V - list contribution calculation
  startTime = omp_get_wtime();
  for (int j = lvlOrdVec.size(); j >= 1; j--){
	 vector<int>& thisLevBoxes = lvlOrdVec[j];
#pragma omp parallel for
	 for(int i=0; i<thisLevBoxes.size(); i++) {
		int gNodeIdx = thisLevBoxes[i];
		if(1 || (let()->tag(gNodeIdx) & LET_SRCNODE)){
		  N& srcPtr = node(gNodeIdx);
		  /* Need to find a good way to keep from recomputing this for dirichlet */
		  /* Seems to be taken care of */
		  srcPtr.effDen().resize( _matmgnt->effDatSze(UE) ); setvalue(srcPtr.effDen(), 0.0);//1. resize effDen
		  iC( _matmgnt->plnDen2EffDen(_let->depth(gNodeIdx)+_rootLevel, srcUpwEquDen(gNodeIdx),  srcPtr.effDen(), 2.0) );			 //2. transform from upeDen to effDen
		}
	 }
#pragma omp parallel for
	 for(int i=0; i<thisLevBoxes.size(); i++) {
		int gNodeIdx = thisLevBoxes[i];
		if( _let->tag(gNodeIdx) & LET_TRGNODE) { //eValuator		//GNTra gnt = _let->gNodeIdx2gnt(gNodeIdx);
		  Point3 gNodeIdxCtr(_let->center(gNodeIdx));
		  double D = 2.0 * _let->radius(gNodeIdx);
		  DblNumVec trgDwnChkVal_gNodeIdx(this->trgDwnChkVal(gNodeIdx));

		  N& trgPtr = node(gNodeIdx);
		  DblNumVec effVal(_matmgnt->effDatSze(DC));
		
		  for(vector<int>::iterator vi=_let->node(gNodeIdx).Vnodes().begin(); vi!=_let->node(gNodeIdx).Vnodes().end(); vi++) {
			 if((let()->tag(*vi) & LET_SRCNODE)){
				N& srcPtr = node(*vi);
				Point3 viCtr(_let->center(*vi));
				Index3 idx;
				for(int d=0; d<dim(); d++){
				  idx(d) = int(round( (viCtr[d]-gNodeIdxCtr[d])/D ));
				}
				//M2L - multipole -> local
				iC( _matmgnt->UpwEqu2DwnChk_dgemv(_let->depth(gNodeIdx)+_rootLevel, idx, srcPtr.effDen(), effVal) );
			 }
		  }
		  iC( _matmgnt->effVal2PlnVal(effVal, trgDwnChkVal_gNodeIdx) );			 //1. transform from effVal to dncVal
		}
	 }
#pragma omp parallel for
	 for(int i=0; i<thisLevBoxes.size(); i++) {
		int gNodeIdx = thisLevBoxes[i];
		N& srcPtr = node(gNodeIdx);
		srcPtr.effDen().resize(0);
	 }
  }
  endTime = omp_get_wtime();
  fmmevaltime += (endTime - startTime);
  std::cout << "V used " << (endTime - startTime) << endl;

	 
  //W - list contrubtion calculation

  trgTrmVec.clear();
  for (int i = 0; i < ordVec.size(); i++){
	 int gni = ordVec[i];
	 if (let()->node(gni).Wnodes().size() > 0 && let()->tag(gni) & LET_TRGNODE) { //evaluator
		trgTrmVec.push_back(gni);
	 }
  }

  if (trgTrmVec.size() > 0){
    startTime = omp_get_wtime();
#pragma omp parallel for
    for(int i=0; i<trgTrmVec.size(); i++) {
      int gNodeIdx = trgTrmVec[i];
      if( _let->tag(gNodeIdx) & LET_TRGNODE ) {
	if( _let->terminal(gNodeIdx)==true ) {
	  DblNumVec trgExaVal_gNodeIdx(this->trgExaVal(gNodeIdx));
	  N& curNode = _let->node(gNodeIdx);
	  for(vector<int>::iterator vi=curNode.Wnodes().begin(); vi!=curNode.Wnodes().end(); vi++) {
	    if(_let->terminal(*vi) && _let->node(*vi).srcExaNum()*srcDOF<_matmgnt->plnDatSze(UE)) { //use Exa instead
	      //S2T - source -> target
	      iC( SrcEqu2TrgChk_dgemv(srcExaPos(*vi), srcExaNor(*vi), trgExaPos(gNodeIdx), srcExaDen(*vi), trgExaVal_gNodeIdx) );
	    } else {
	      //M2T - multipole -> target
	      int vni = *vi;				
	      iC( UpwEqu2TrgChk_dgemv(_let->center(vni), _let->radius(vni), trgExaPos(gNodeIdx), srcUpwEquDen(*vi), trgExaVal_gNodeIdx) );
	    }
	  }
	}
      }
    }
	 endTime = omp_get_wtime();
	 fmmevaltime += (endTime - startTime);
	 std::cout << "W used " << (endTime - startTime) << endl;
  }
  
  trgTrmVec.clear();
  for (int i = 0; i < ordVec.size(); i++){
	 int gni = ordVec[i];
	 if (let()->node(gni).Xnodes().size() > 0 && let()->tag(gni) & LET_TRGNODE) { //evaluator
		trgTrmVec.push_back(gni);
	 }
  }

  //X - list contrubtion calculation
  if (trgTrmVec.size() > 0){
    startTime = omp_get_wtime();
#pragma omp parallel for
    for(int i=0; i<trgTrmVec.size(); i++) {
      int gNodeIdx = trgTrmVec[i];
      if( _let->tag(gNodeIdx) & LET_TRGNODE) {	
	N& curNode = _let->node(gNodeIdx);
	DblNumVec trgExaVal_gNodeIdx(trgExaVal(gNodeIdx));
	DblNumVec trgDwnChkVal_gNodeIdx(trgDwnChkVal(gNodeIdx));
	for(vector<int>::iterator vi=curNode.Xnodes().begin(); vi!=curNode.Xnodes().end(); vi++) {
	  if(_let->terminal(gNodeIdx) && _let->node(gNodeIdx).trgExaNum()*trgDOF<_matmgnt->plnDatSze(DC)) { //use Exa instead
	    iC( SrcEqu2TrgChk_dgemv(srcExaPos(*vi), srcExaNor(*vi), trgExaPos(gNodeIdx), srcExaDen(*vi), trgExaVal_gNodeIdx) );
	  } else {
	    //S2L - source -> local
	    iC( SrcEqu2DwnChk_dgemv(srcExaPos(*vi), srcExaNor(*vi), _let->center(gNodeIdx), _let->radius(gNodeIdx), srcExaDen(*vi), trgDwnChkVal_gNodeIdx) );
	  }
	}
      }
    }	
	 endTime = omp_get_wtime();
	 fmmevaltime += (endTime - startTime);
	 std::cout << "X used " << (endTime - startTime) << endl;
  }
  
  //7. combine
  //Pre-build matrices for openmp
  for (int j = 0; j < levnum; j++){
	 DblNumVec chkVal(_matmgnt->plnDatSze(DC));
	 DblNumVec denVal(datSze(DE));
	 iC( _matmgnt->DwnChk2DwnEqu_dgemv(j+_rootLevel, chkVal, denVal) );
	 for(int a=0; a<2; a++) { for(int b=0; b<2; b++) { for(int c=0; c<2; c++) {
			 Index3 idx(a,b,c);
			 iC( _matmgnt->DwnEqu2DwnChk_dgemv(j+_rootLevel, idx, denVal, chkVal) );
		  } } }
  }
  
  int rt = 2; //Will need to change for particle periodicity
  startTime = omp_get_wtime();
  for (int j = 0; j < lvlOrdVec.size(); j ++){
	 vector<int>& thisLevBoxes = lvlOrdVec[j];
#pragma omp parallel for
	 for(int i=0; i<thisLevBoxes.size(); i++) {
		int gNodeIdx = thisLevBoxes[i];
		DblNumVec trgDwnEquDen_gNodeIdx(datSze(DE));
		if(_let->depth(gNodeIdx)>=rt) {
		  //L2L - local -> local
		  iC( _matmgnt->DwnChk2DwnEqu_dgemv(_let->depth(gNodeIdx)+_rootLevel, trgDwnChkVal(gNodeIdx), trgDwnEquDen_gNodeIdx) );
		}
		if(_let->depth(gNodeIdx)>=rt && !(_let->terminal(gNodeIdx))) {
		  for (int a = 0; a < 2; a++){ for (int b = 0; b < 2; b++){ for (int c = 0; c < 2; c++){
				  Index3 chdIdx(a,b,c); int chi = _let->child(gNodeIdx, chdIdx); iA( chi != -1 );
				  DblNumVec trgDwnChkVal_chi(trgDwnChkVal(chi));
				  iC( _matmgnt->DwnEqu2DwnChk_dgemv(_let->depth(gNodeIdx)+_rootLevel, chdIdx, trgDwnEquDen_gNodeIdx, trgDwnChkVal_chi) );
				} } }
		}
		if(_let->terminal(gNodeIdx)) {
		  //L2T - local -> target
		  DblNumVec trgExaVal_gNodeIdx(trgExaVal(gNodeIdx));
		  iC( DwnEqu2TrgChk_dgemv(_let->center(gNodeIdx), _let->radius(gNodeIdx), trgExaPos(gNodeIdx), trgDwnEquDen_gNodeIdx, trgExaVal_gNodeIdx) );
		}
	 }
  }
  endTime = omp_get_wtime();
  fmmevaltime += (endTime - startTime);
  std::cout << "L2L/L2T used " << (endTime - startTime) << endl;

  cerr << "TOTAL TIME = " << fmmevaltime << endl;
  
  trgTrmVec.clear();
  for (int i = 0; i < ordVec.size(); i++){
	 int gni = ordVec[i];
	 if (let()->terminal(gni) && let()->tag(gni) & LET_TRGNODE) { //evaluator
		trgTrmVec.push_back(gni);
	 }
  }
  
  //8. save trgExaVal
#pragma omp parallel for
  for(int i=0; i<trgTrmVec.size(); i++) {
	 int gNodeIdx = trgTrmVec[i];
	 DblNumVec trgExaVal(this->trgExaVal(gNodeIdx));
	 vector<int>& curVecIdxs = _let->node(gNodeIdx).trgOwnVecIdxs();
	 for(int k=0; k<curVecIdxs.size(); k++) {
		int poff = curVecIdxs[k];
		for(int d=0; d<trgDOF; d++) {
		  trgVal(poff*trgDOF+d) = trgExaVal(k*trgDOF+d);
		}
	 }
  }
  return (0);
}



