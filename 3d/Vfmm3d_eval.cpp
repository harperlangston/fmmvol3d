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

#define PER_SCALE 3.0

#include "Vfmm3d.hpp"
#include "common/vecmatop.hpp"

#include <omp.h>
#include <time.h>

using namespace std;

// ----------------------------------------------------------------------
template <class VF>
int VFMM3d<VF>::evaluate(const DblNumVec& srcDen, DblNumVec& trgVal)
{  
  cout << trgVal.m() << " " << (*(this->trgPos())).n() << endl;
  //DblNumVec& tval = *(_trgVal);
  //-----------------------------------
  //iA(srcDen.m()==srcDOF()*(*_srcPos).n());
  int srcDOF = this->srcDOF();
  int trgDOF = this->trgDOF();

  int NK = vlet()->srcNk();

  double fmmevaltime = 0.0;
  bool trmsort = false;
  
  //1. zero out Vecs
  setvalue(trgVal, 0.0);
  setvalue(_grdExaVal, 0.0);
  
  int srcNodeCnt = vlet()->srcNodeCnt();
  (this->srcUpwEquDen()).resize(srcNodeCnt * (this->datSze(UE)));
  setvalue(this->srcUpwEquDen(), 0.0);
  
  setvalue(this->trgExaVal(), 0.0);

  int NMTHRDS = 1;

  vector<int> ordVec; iC( vlet()->upwOrderCollect(ordVec) ); //BOTTOM UP

  /* Clear out the vector Idxs for non-leaves - saves memory */
  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if (!vlet()->terminal(gNodeIdx)){
		Node& curNode =  this->node(gNodeIdx);
		curNode.srcOwnVecIdxs().clear();
	 }
  }
  
  vector<int> trgTrmVec;
  for (int i = 0; i < ordVec.size(); i++){
	 int gni = ordVec[i];
	 if (vlet()->terminal(gni) && vlet()->tag(gni) & LET_TRGNODE && vlet()->tag(gni) & LET_SRCNODE) { //evaluator
		trgTrmVec.push_back(gni);
	 }
  }

  cerr << "Number of source terminals " << trgTrmVec.size() << endl;

  map<int, vector<int> > lvlOrdVec; iC( vlet()->revLvlOrderCollect(lvlOrdVec) ); //BOTTOM UP
  vector<bool> lvsAtLev;
  vector<bool> ndsAtLev;
  for (int j = 0; j < lvlOrdVec.size(); j++){
	 lvsAtLev.push_back(false);
	 ndsAtLev.push_back(false);
	 vector<int>& thisLevBoxes = lvlOrdVec[j];
	 for (int i = 0; i < thisLevBoxes.size(); i++){
		int gNodeIdx = thisLevBoxes[i];
		if (vlet()->terminal(gNodeIdx) == true) { lvsAtLev[j] = true; }
		/* If in thi loop, there is a node at this level */
		ndsAtLev[j] = true;
	 }
  }
  
  cout << (this->knl()).kernelType() << endl;
  //3. up computation
  if (1){
	 cout << "S2M and M2M calculations" << endl;
	 //Pre-build matrices for openmp
	 int levnum = ((this->matmgnt())->hom() ? 1 : lvlOrdVec.size());
	 for (int j = 0; j < levnum; j++){
		vector<int>& thisLevBoxes = lvlOrdVec[j];
		DblNumVec chkVal((this->matmgnt())->plnDatSze(UC));
		DblNumVec tmpVec((this->datSze(UE)));
		if (lvsAtLev[j] == true || (this->matmgnt())->hom()){
		  DblNumVec coeffs(vlet()->srcNk()); setvalue(coeffs,1.0);
		  iC( GrdEqu2UpwEquTbls_dgemv(j+(this->rootLevel()), coeffs, chkVal, tmpVec) );
		}
		if (ndsAtLev[j] == true || (this->matmgnt())->hom()){
		  for(int a=0; a<2; a++) { for(int b=0; b<2; b++) { for(int c=0; c<2; c++) {
				  Index3 idx(a,b,c);
				  iC( (this->matmgnt())->UpwEqu2UpwChk_dgemv((j+1)+(this->rootLevel()), idx, tmpVec, chkVal, 2.0) );
				} } }
		  iC( (this->matmgnt())->UpwChk2UpwEqu_dgemv(j + (this->rootLevel()), chkVal, tmpVec, 2.0) );
		}
	 }

	 double startTime = omp_get_wtime();
	 for (int j = lvlOrdVec.size(); j >= 0; j--){
		vector<int>& thisLevBoxes = lvlOrdVec[j];
		if (ndsAtLev[j] == true){
#pragma omp parallel for
		  for(int i=0; i<thisLevBoxes.size(); i++) {
			 int gNodeIdx = thisLevBoxes[i];
			 if (omp_get_thread_num()+1 > NMTHRDS) { NMTHRDS = omp_get_thread_num()+1; }
			 if(vlet()->tag(gNodeIdx) & LET_SRCNODE) {
				DblNumVec srcUpwEquDengNodeIdx(this->srcUpwEquDen(gNodeIdx));
				if(vlet()->terminal(gNodeIdx)) {
				  DblNumVec chkVal((this->matmgnt())->plnDatSze(UC));
				  setvalue(chkVal, 0.0);

				  iC( GrdEqu2UpwEquTbls_dgemv(vlet()->depth(gNodeIdx)+(this->rootLevel()), srcCoeffs(gNodeIdx), chkVal, srcUpwEquDengNodeIdx) );
				} else {				
				  //M2M - Multipole -> Multipole
				  DblNumVec chkVal((this->matmgnt())->plnDatSze(UC));
				  for(int a=0; a<2; a++) {
					 for(int b=0; b<2; b++) {
						for(int c=0; c<2; c++) {
						  Index3 idx(a,b,c);
						  int chi = vlet()->child(gNodeIdx, idx);
						  if(vlet()->tag(chi) & LET_SRCNODE) {
							 DblNumVec srcUpwEquDenChi(this->srcUpwEquDen(chi));
							 iC( (this->matmgnt())->UpwEqu2UpwChk_dgemv(vlet()->depth(chi)+(this->rootLevel()), idx, srcUpwEquDenChi, chkVal, 2.0) );
							 //iC( (this->matmgnt())->UpwEqu2UpwEqu_dgemv(vlet()->depth(chi)+(this->rootLevel()), idx, srcUpwEquDenChi, srcUpwEquDengNodeIdx, 2.0) );
						  }
						}		
					 }
				  }
				  iC( (this->matmgnt())->UpwChk2UpwEqu_dgemv(vlet()->depth(gNodeIdx) + (this->rootLevel()), chkVal, srcUpwEquDengNodeIdx, 2.0) );
				}		
			 }
		  }
		}
	 }
	 double endTime = omp_get_wtime();
	 cout << "S2M used " << (endTime - startTime) << endl;
	 fmmevaltime += (endTime - startTime);
  
  }

  trgTrmVec.clear();
  for (int i = 0; i < ordVec.size(); i++){
	 int gni = ordVec[i];
	 bool inchk = false;
	 if ((vlet()->terminal(gni) && vlet()->node(gni).Unodes().size() > 0) && vlet()->tag(gni) & LET_TRGNODE) {
		trgTrmVec.push_back(gni);
		inchk = true;
	 }
	 if (vlet()->periodic() || vlet()->dirichlet()){
		if (inchk == false){
		  if (vlet()->pernode(gni).bdryUnodes().size() > 0){
			 trgTrmVec.push_back(gni);
		  }
		}
	 }
  }

  if (trmsort == true && trgTrmVec.size() > 0){
	 vector<int> trgTrmSortedVec; trgTrmSortedVec.resize(trgTrmVec.size());
	 trgTrmSortedVec[0] = trgTrmVec[0];
	 for (int i = 1; i < trgTrmVec.size(); i++){
		trgTrmSortedVec[i] = trgTrmVec[i];
		for (int j = i; j > 0; j--){
		  int gnj = trgTrmSortedVec[j];
		  int gnjm1 = trgTrmSortedVec[j-1];
		  if ((vlet()->node(gnj)).Unodes().size() > (vlet()->node(gnjm1)).Unodes().size()){
			 // Swap them
			 trgTrmSortedVec[j] = gnjm1;
			 trgTrmSortedVec[j-1] = gnj;
		  }
		}
	 }
	 
	 int NT = NMTHRDS;
	 int nperp = trgTrmVec.size()/NT;
	 int ndif = trgTrmVec.size() - ((nperp)*NT);
	 
	 int wcnt = 0;
	 for (int j = 0; j < nperp; j++){
		for (int i = 0; i < NT; i++){
		  if (wcnt < trgTrmVec.size()){
			 trgTrmVec[i*(nperp) + j] = trgTrmSortedVec[wcnt];
			 //sumU[i] = vlet()->node(trgTrmVec[i*(nperp+1)+j)
			 wcnt++;
		  }
		}
	 }
	 
	 cerr << NT << " " << (nperp) << " " << ndif << endl;
	 
	 if (ndif > 0 && wcnt != trgTrmVec.size()){
		iA(ndif != 0);
		for (int i = 0; i < ndif; i++){
		  trgTrmVec[wcnt] = trgTrmSortedVec[wcnt];
		  wcnt++;
		}
	 }
  }
  
  ordVec.clear();  iC( vlet()->dwnOrderCollect(ordVec) ); //TOP DOWN
  //U - list contribution calculation
  if (1){
	 int levnum = ((this->matmgnt())->hom() ? 1 : lvlOrdVec.size());
	 for (int lev = 0; lev < levnum; lev++){
		vector<int>& thisLevBoxes = lvlOrdVec[lev];
		if (lvsAtLev[lev] == true || (this->matmgnt())->hom()){
		  DblNumVec cfs(vlet()->srcNk());
		  DblNumVec vals(vlet()->trgGrdSze()*trgDOF);
		  for (int i = 0; i < 3; i++){
			 int nbrSze = (i == 0 ? 27 : 56);
			 for (int j = 0; j < nbrSze; j++){
				GrdCoeffs2TrgVal_dgemv(lev+(this->rootLevel()), i, j, cfs, vals);
			 }
		  }
		}
		if (!(vlet()->balance())){
		  for(int i=0; i<trgTrmVec.size(); i++) {
			 int gNodeIdx = trgTrmVec[i];
			 Node& curNode = vlet()->node(gNodeIdx);
			 for(vector<int>::iterator vi=curNode.Unodes().begin(); vi!=curNode.Unodes().end(); vi++) {
				if (vlet()->tag(*vi) & LET_SRCNODE && vlet()->tag(gNodeIdx) * LET_TRGNODE){
				  int depdiff = abs(vlet()->depth(gNodeIdx) - vlet()->depth(*vi));
				  if (depdiff > 1){
					 DblNumVec cfs(vlet()->srcNk());
					 DblNumVec vals(vlet()->trgGrdSze()*trgDOF);
					 int nbrTyp = nbrType(*vi, gNodeIdx);
					 iC( GrdCoeffs2UnbalTrgVal_dgemv(nbrTyp, gNodeIdx, *vi, cfs, vals));
				  }
				}
			 }
		  }
		}
	 }
  }
  if (1){
    cout << "U - list contribution calculation" << endl;
	 double startTime = omp_get_wtime();
	 for (int tt = 0; tt < 56; tt++){
#pragma omp parallel for
		for(int i=0; i<trgTrmVec.size(); i++) {
		  int gNodeIdx = trgTrmVec[i];
		  if( vlet()->tag(gNodeIdx) & LET_TRGNODE) { 
			 Node& curNode = vlet()->node(gNodeIdx);
			 DblNumVec grdExaValgNodeIdx(grdExaVal(gNodeIdx));
			 for(vector<int>::iterator vi=curNode.Unodes().begin(); vi!=curNode.Unodes().end(); vi++) {
				if((vlet()->tag(*vi) & LET_SRCNODE)){
				  DblNumVec coeffsvi(srcCoeffs(*vi));
				  int nbrTyp = nbrType(*vi, gNodeIdx);
				  /* Apply *vi's coefficients to gNodeIdx */
				  if (abs(vlet()->depth(gNodeIdx) - vlet()->depth(*vi)) <= 1){
					 int tabNbr = tabNbrType(nbrTyp, *vi, gNodeIdx);
					 if (tabNbr == tt) {
						iC( GrdCoeffs2TrgVal_dgemv(vlet()->depth(*vi)+(this->rootLevel()), nbrTyp, tabNbr, coeffsvi, grdExaValgNodeIdx));
					 }
				  }
				  else {
					 if (tt == 0){
						iC( GrdCoeffs2UnbalTrgVal_dgemv(nbrTyp, gNodeIdx, *vi, coeffsvi, grdExaValgNodeIdx));
					 }
				  }
				}
			 }
		  }
		}
	 }
#pragma omp parallel for
	 for(int i=0; i<trgTrmVec.size(); i++) {
		int gNodeIdx = trgTrmVec[i];
		if( vlet()->tag(gNodeIdx) & LET_TRGNODE) { 
		  if (vlet()->periodic()){ evaluateUnodes_BCs(gNodeIdx); }
		}
	 }
	 double endTime = omp_get_wtime();
	 fmmevaltime += (endTime - startTime);
	 cout << "NEAR used " << (endTime - startTime) << endl;
  }

  if ( (this->trgDwnChkVal()).m() == 0) (this->trgDwnChkVal()).resize(vlet()->trgNodeCnt() * (this->datSze(DC)));
  setvalue(this->trgDwnChkVal(), 0.0);

  trgTrmVec.clear();
  for (int i = 0; i < ordVec.size(); i++){
	 int gni = ordVec[i];
	 bool inchk = false;
	 if (vlet()->node(gni).Xnodes().size() > 0 || vlet()->node(gni).Wnodes().size() > 0) { //evaluator) { //evaluator
	   trgTrmVec.push_back(gni);
		inchk = true;
	 }
	 if (vlet()->periodic() || vlet()->dirichlet()){
		if (inchk == false){
		  if (vlet()->pernode(gni).bdryXnodes().size() > 0 || vlet()->pernode(gni).bdryWnodes().size() > 0) {
			 trgTrmVec.push_back(gni);
		  }
		}
	 }
  }

  /* Build tables in full  - speedup for testing */
  int levnum = ((this->matmgnt())->hom() ? 1 : lvlOrdVec.size());
  for (int lev = 0; lev < levnum; lev++){
	 vector<int>& thisLevBoxes = lvlOrdVec[lev];
	 if (lvsAtLev[lev] == true || (this->matmgnt())->hom()){
		DblNumVec cfs(vlet()->srcNk());
		DblNumVec vals(vlet()->trgGrdSze()*trgDOF);
		for (int i = WNODE; i <= XNODE; i++){
		  int nbrSze = 152;
		  for (int j = 0; j < nbrSze; j++){
			 WXcoeffs2TrgVal(lev+(this->rootLevel()), i, j, cfs, vals, false);
		  }
		}
	 }
  }
  if (!(vlet()->balance()) || 1){
	 for(int i=0; i<trgTrmVec.size(); i++) {
		int gNodeIdx = trgTrmVec[i];
		Node& curNode = vlet()->node(gNodeIdx);
		for(vector<int>::iterator vi=curNode.Wnodes().begin(); vi!=curNode.Wnodes().end(); vi++) {
		  bool dep = (abs(vlet()->depth(gNodeIdx) - vlet()->depth(*vi)) == 1);
		  if(!(vlet()->terminal(*vi) && dep)){
			 int tabNode = _tbls->lookup(XNODE,vlet()->depth(gNodeIdx)+(this->rootLevel()), vlet()->radius(gNodeIdx), vlet()->center(gNodeIdx), vlet()->depth(*vi)+(this->rootLevel()), vlet()->radius(*vi), vlet()->center(*vi));
			 /* Apply *vi's coefficients to gNodeIdx - so using XNODE, not WNODE */
			 iC( WXUpwEqu2DwnChk_dgemv(XNODE, *vi, gNodeIdx, (tabNode), false) );
		  }
		}
		for(vector<int>::iterator vi=curNode.Xnodes().begin(); vi!=curNode.Xnodes().end(); vi++) {
		  if((vlet()->tag(*vi) & LET_SRCNODE)){
			 bool dep = (abs(vlet()->depth(gNodeIdx) - vlet()->depth(*vi)) == 1);
			 if(!(vlet()->terminal(gNodeIdx) && dep)){	
				int tabNode = _tbls->lookup(WNODE,vlet()->depth(*vi)+(this->rootLevel()), vlet()->radius(*vi), vlet()->center(*vi), vlet()->depth(gNodeIdx)+(this->rootLevel()), vlet()->radius(gNodeIdx), vlet()->center(gNodeIdx));
				iC( WXUpwEqu2DwnChk_dgemv(WNODE, *vi, gNodeIdx, tabNode, false) );
			 }
		  }
		}
	 }
  }

  trgTrmVec.clear();
  for (int i = 0; i < ordVec.size(); i++){
	 int gni = ordVec[i];
	 bool inchk = false;
	 if (vlet()->node(gni).Wnodes().size() > 0) { //evaluator) { //evaluator
		trgTrmVec.push_back(gni);
		inchk = true;
	 }
  }
	  
  if (trmsort == true && trgTrmVec.size() > 0){
	 vector<int> trgTrmSortedVec; trgTrmSortedVec.resize(trgTrmVec.size());
	 trgTrmSortedVec[0] = trgTrmVec[0];
	 for (int i = 1; i < trgTrmVec.size(); i++){
		trgTrmSortedVec[i] = trgTrmVec[i];
		for (int j = i; j > 0; j--){
		  int gnj = trgTrmSortedVec[j];
		  int gnjm1 = trgTrmSortedVec[j-1];
		  if ((vlet()->node(gnj)).Wnodes().size() > (vlet()->node(gnjm1)).Wnodes().size()){
			 // Swap them
			 trgTrmSortedVec[j] = gnjm1;
			 trgTrmSortedVec[j-1] = gnj;
		  }
		}
	 }
	 
	 int NT = NMTHRDS;
	 int nperp = trgTrmVec.size()/NT;
	 int ndif = trgTrmVec.size() - (nperp*NT);
	 
	 int wcnt = 0;
	 for (int j = 0; j < nperp; j++){
		for (int i = 0; i < NT; i++){
		  trgTrmVec[i*nperp + j] = trgTrmSortedVec[wcnt];
		  wcnt++;
		}
	 }
	 
	 if (wcnt != trgTrmVec.size()){
		iA(ndif != 0);
		for (int i = 0; i < ndif; i++){
		  trgTrmVec[wcnt] = trgTrmSortedVec[wcnt];
		  wcnt++;
		}
	 }
  }
	 
  if (1){
    cout << "W/X - list contribution calculation" << endl;
    //W - list contribution calculation
	 double startTime = omp_get_wtime();
	 for (int tt = 0; tt < 152; tt++){
#pragma omp parallel for
		for(int i=0; i<trgTrmVec.size(); i++) {
		  int gNodeIdx = trgTrmVec[i];
		  Node& curNode =  vlet()->node(gNodeIdx);
		  for(vector<int>::iterator vi=curNode.Wnodes().begin(); vi!=curNode.Wnodes().end(); vi++) {
			 if((vlet()->tag(*vi) & LET_SRCNODE)){
				bool dep = (abs(vlet()->depth(gNodeIdx) - vlet()->depth(*vi)) == 1);
				DblNumVec grdExaValgNodeIdx(this->grdExaVal(gNodeIdx));
				if(vlet()->terminal(*vi) && dep){
				  int tabNode = XNodeType(*vi,gNodeIdx);
				  if (tabNode == tt){
					 iC( WXcoeffs2TrgVal(vlet()->depth(*vi)+(this->rootLevel()), XNODE, (tabNode), srcCoeffs(*vi), grdExaValgNodeIdx));
				  }
				}
				else {
				  if (tt == 0){
					 int tabNode = _tbls->lookup(XNODE,vlet()->depth(gNodeIdx)+(this->rootLevel()), vlet()->radius(gNodeIdx), vlet()->center(gNodeIdx), vlet()->depth(*vi)+(this->rootLevel()), vlet()->radius(*vi), vlet()->center(*vi));
					 /* Apply *vi's coefficients to gNodeIdx - so using XNODE, not WNODE */
					 iC( WXUpwEqu2DwnChk_dgemv(XNODE, *vi, gNodeIdx, (tabNode)));
				  }
				}
			 }
		  }
		}
	 }
  trgTrmVec.clear();
  for (int i = 0; i < ordVec.size(); i++){
	 int gni = ordVec[i];
	 if (vlet()->periodic() || vlet()->dirichlet()){
		if (vlet()->pernode(gni).bdryWnodes().size() > 0) {
		  trgTrmVec.push_back(gni);
		}
	 }
  }
#pragma omp parallel for
	 for(int i=0; i<trgTrmVec.size(); i++) {
		int gNodeIdx = trgTrmVec[i];
		if (vlet()->periodic()){ evaluateWnodes_BCs(gNodeIdx); }
	 }
	 double endTime = omp_get_wtime();
	 fmmevaltime += (endTime - startTime);
  }
  cout << "W - list contribution calculation done" << endl;
  trgTrmVec.clear();
  for (int i = 0; i < ordVec.size(); i++){
	 int gni = ordVec[i];
	 bool inchk = false;
	 if (vlet()->node(gni).Xnodes().size() > 0) { //evaluator) { //evaluator
		trgTrmVec.push_back(gni);
		inchk = true;
	 }
  }
  
  if (trmsort == true && trgTrmVec.size() > 0){
	 vector<int> trgTrmSortedVec; trgTrmSortedVec.resize(trgTrmVec.size());
	 trgTrmSortedVec[0] = trgTrmVec[0];
	 for (int i = 1; i < trgTrmVec.size(); i++){
		trgTrmSortedVec[i] = trgTrmVec[i];
		for (int j = i; j > 0; j--){
		  int gnj = trgTrmSortedVec[j];
		  int gnjm1 = trgTrmSortedVec[j-1];
		  if ((vlet()->node(gnj)).Xnodes().size() > (vlet()->node(gnjm1)).Xnodes().size()){
			 // Swap them
			 trgTrmSortedVec[j] = gnjm1;
			 trgTrmSortedVec[j-1] = gnj;
		  }
		}
	 }
	 
	 int NT = NMTHRDS;
	 int nperp = trgTrmVec.size()/NT;
	 int ndif = trgTrmVec.size() - (nperp*NT);
	 
	 int wcnt = 0;
	 for (int j = 0; j < nperp; j++){
		for (int i = 0; i < NT; i++){
		  trgTrmVec[i*nperp + j] = trgTrmSortedVec[wcnt];
		  wcnt++;
		}
	 }
	 
	 if (wcnt != trgTrmVec.size()){
		iA(ndif != 0);
		for (int i = 0; i < ndif; i++){
		  trgTrmVec[wcnt] = trgTrmSortedVec[wcnt];
		  wcnt++;
		}
	 }
  }


  if (1) {
	 cout << "W/X - list contribution calculation" << endl;
    //W - list contribution calculation
	 double startTime = omp_get_wtime();
	 for (int tt = 0; tt < 152; tt++){
#pragma omp parallel for
		for(int i=0; i< trgTrmVec.size(); i++) {
		  int gNodeIdx = trgTrmVec[i];
		  Node& curNode =  vlet()->node(gNodeIdx);
		  for(vector<int>::iterator vi=curNode.Xnodes().begin(); vi!=curNode.Xnodes().end(); vi++) {
			 if((vlet()->tag(*vi) & LET_SRCNODE)){
				bool dep = (abs(vlet()->depth(gNodeIdx) - vlet()->depth(*vi)) == 1);
				iA(vlet()->terminal(*vi));
				if(vlet()->terminal(gNodeIdx) && dep){
				  int tabNode = WNodeType(*vi,gNodeIdx);
				  if (tabNode == tt){
					 DblNumVec grdExaValVi(this->grdExaVal(*vi));
					 DblNumVec grdExaValgNodeIdx(this->grdExaVal(gNodeIdx));
					 iC( WXcoeffs2TrgVal(vlet()->depth(*vi)+(this->rootLevel()), WNODE, (tabNode), srcCoeffs(*vi), grdExaValgNodeIdx));
				  }
				}
				else {
				  if (tt == 0){
					 int tabNode = _tbls->lookup(WNODE,vlet()->depth(*vi)+(this->rootLevel()), vlet()->radius(*vi), vlet()->center(*vi), vlet()->depth(gNodeIdx)+(this->rootLevel()), vlet()->radius(gNodeIdx), vlet()->center(gNodeIdx));
					 //int tabNode = _tbls->lookup(WNODE,vlet()->depth(gNodeIdx)+(this->rootLevel()), vlet()->radius(gNodeIdx), vlet()->center(gNodeIdx), vlet()->depth(*vi)+(this->rootLevel()), vlet()->radius(*vi), vlet()->center(*vi));
					 iC( WXUpwEqu2DwnChk_dgemv(WNODE, *vi, gNodeIdx, tabNode));
				  }
				}
			 }
		  }
		}
	 }
	 trgTrmVec.clear();
	 for (int i = 0; i < ordVec.size(); i++){
		int gni = ordVec[i];
		if (vlet()->periodic() || vlet()->dirichlet()){
		  if (vlet()->pernode(gni).bdryXnodes().size() > 0) {
			 trgTrmVec.push_back(gni);
		  }
		}
	 }
#pragma omp parallel for
	 for(int i=0; i< trgTrmVec.size(); i++) {
		int gNodeIdx = trgTrmVec[i];
		if (vlet()->periodic()){ evaluateXnodes_BCs(gNodeIdx); }
	 }
	 double endTime = omp_get_wtime();
	 fmmevaltime += (endTime - startTime);
  }

  trgTrmVec.clear();
  
  (*_coeffs).resize(0);
  cleanSrcEqu2UpwChkTbls();
  
  /* Done with _dcuhre - free some space */
  delete _tbls; _tbls = NULL;

  if (1) {
	 //Pre-build matrices for openmp
	 if (1){
		int levnum = ((this->matmgnt())->hom() ? 1 : lvlOrdVec.size());
		for (int j = 0; j < levnum; j++){
		  DblNumVec effVal((this->matmgnt())->effDatSze(DC));
		  DblNumVec effDen((this->matmgnt())->effDatSze(UE));
		
		  //M2L - multipole -> local
		  for(int a=-3; a<=3; a++) { for(int b=-3; b<=3; b++) { for(int c=-3; c<=3; c++) {
				  Index3 idx(a,b,c);
				  if (idx.linfty() > 1) {
				    iC( (this->matmgnt())->UpwEqu2DwnChk_dgemv(j+(this->rootLevel()), idx, effDen, effVal, 2.0) );
				  } } } }
		}
	 }
    cout << "V - list contribution calculation" << endl;
	 //V - list contribution calculation
	 double startTime = omp_get_wtime();
	 for (int j = lvlOrdVec.size(); j >= 1; j--){
		vector<int>& thisLevBoxes = lvlOrdVec[j];
		if (ndsAtLev[j] == true){
#pragma omp parallel for
		  for(int i=0; i<thisLevBoxes.size(); i++) {
		    int gNodeIdx = thisLevBoxes[i];
			 if(1 || (vlet()->tag(gNodeIdx) & LET_SRCNODE)){
				Node& srcPtr = this->node(gNodeIdx);
				/* Need to find a good way to keep from recomputing this for dirichlet */
				/* Seems to be taken care of */
				srcPtr.effDen().resize( (this->matmgnt())->effDatSze(UE) ); setvalue(srcPtr.effDen(), 0.0);//1. resize effDen
				iC( (this->matmgnt())->plnDen2EffDen(vlet()->depth(gNodeIdx)+(this->rootLevel()), this->srcUpwEquDen(gNodeIdx),  srcPtr.effDen(), 2.0) );			 //2. transform from upeDen to effDen
			 }
		  }	
		}
		if (ndsAtLev[j] == true){
#pragma omp parallel for
		  for(int i=0; i<thisLevBoxes.size(); i++) {
			 int gNodeIdx = thisLevBoxes[i];
			 if( vlet()->tag(gNodeIdx) & LET_TRGNODE) { //eValuator		//GNTra gnt = vlet()->gNodeIdx2gnt(gNodeIdx);
				Point3 gNodeIdxCtr(vlet()->center(gNodeIdx));
				double D = 2.0 * vlet()->radius(gNodeIdx);
				DblNumVec trgDwnChkVal_gNodeIdx(this->trgDwnChkVal(gNodeIdx));
			 
				Node& trgPtr = this->node(gNodeIdx);
				DblNumVec effVal((this->matmgnt())->effDatSze(DC));
			 
				/* If periodic, evaluate these interactions */
				if (vlet()->periodic()){ evaluateVnodes_BCs(gNodeIdx, effVal); } 
				for(vector<int>::iterator vi=vlet()->node(gNodeIdx).Vnodes().begin(); vi!=vlet()->node(gNodeIdx).Vnodes().end(); vi++) {
				  if((vlet()->tag(*vi) & LET_SRCNODE)){
					 Node& srcPtr = this->node(*vi);
					 iA (srcPtr.effDen().m() != 0);
					 Point3 viCtr(vlet()->center(*vi));
					 Index3 idx;
					 for(int d=0; d<(this->dim()); d++){
						idx(d) = int(round( (viCtr[d]-gNodeIdxCtr[d])/D ));
					 }
					 //M2L - multipole -> local
					 iC( (this->matmgnt())->UpwEqu2DwnChk_dgemv(vlet()->depth(gNodeIdx) + (this->rootLevel()), idx, srcPtr.effDen(), effVal, 2.0) );
				  }
				}
				( (this->matmgnt())->effVal2PlnVal(effVal, trgDwnChkVal_gNodeIdx) );			 //1. transform from effVal to dncVal
			 }
		  }
#pragma omp parallel for
		  for(int i=0; i<thisLevBoxes.size(); i++) {
			 int gNodeIdx = thisLevBoxes[i];
			 Node& srcPtr = this->node(gNodeIdx);
			 srcPtr.effDen().resize(0);
		  }
		}
	 }
	 
	 double endTime = omp_get_wtime();
	 fmmevaltime += (endTime - startTime);
	 cout << "FFF used " << (endTime - startTime) << endl;

  }
    
  /* Make sure everything is cleared */
  for(int i=1; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 Node& Ptr = this->node(gNodeIdx);
	 DblNumVec trgDwnChkVal_gNodeIdx(this->trgDwnChkVal(gNodeIdx));
	 //cerr << gNodeIdx << " " <<  trgDwnChkVal_gNodeIdx << endl;
	 if (vlet()->terminal(gNodeIdx)){
		DblNumVec grdExaVal_gNodeIdx(grdExaVal(gNodeIdx));
		//cerr << gNodeIdx << " " <<  grdExaVal_gNodeIdx << endl;
	 }
	 iA (Ptr.effDen().m() == 0);
  }


  if ((vlet()->periodic() || (vlet()->dirichlet()))) {
	 iA(vlet()->root(0));
	 //double startTime = omp_get_wtime();
	 evaluate_far_BCs(0, (this->matmgnt())->perdirmaps()->rtTrgDwnEquDen());
	 //double endTime = omp_get_wtime();
	 //fmmevaltime += (endTime - startTime);
  }

  /*
  iC( (this->matmgnt())->cleanPlans());
  if ((this->matmgnt())->hom())
  	iC( (this->matmgnt())->UpwEqu2DwnChkCleanUp()); 
  _srcUpwEquDen.resize(0);
  */
  

  //7. combine
  cout << "L2L and L2T calculations" << endl;
  int rt = ((vlet()->periodic() || (vlet()->dirichlet())) ? 0 : 2);
  //Pre build matrices
  {
    for (int j = 0; j < lvlOrdVec.size(); j ++){
      vector<int>& thisLevBoxes = lvlOrdVec[j];
      if (ndsAtLev[j] == true){
	DblNumVec chkVal((this->matmgnt())->plnDatSze(DC));
	DblNumVec denVal((this->datSze(DE)));
	iC( (this->matmgnt())->DwnChk2DwnEqu_dgemv(j+(this->rootLevel()), chkVal, denVal, 2.0) );
	for(int a=0; a<2; a++) { for(int b=0; b<2; b++) { for(int c=0; c<2; c++) {
				     Index3 idx(a,b,c);
				     iC( (this->matmgnt())->DwnEqu2DwnChk_dgemv(j+(this->rootLevel()), idx, denVal, chkVal, 2.0) );
				   } } }
      }
      if (lvsAtLev[j] == true){
	DblNumVec denVal((this->datSze(DE))); setvalue(denVal,1.0); 
	DblNumVec grdVal(vlet()->trgGrdSze()); setvalue(grdVal,1.0); 
	iC( DwnEqu2GrdChk_dgemv(j+(this->rootLevel()), denVal, grdVal, false) );
      }
    }
  }
  double startTime = omp_get_wtime();
  for (int j = 0; j < lvlOrdVec.size(); j ++){
    vector<int>& thisLevBoxes = lvlOrdVec[j];
    if (ndsAtLev[j] == true){
#pragma omp parallel for
      for(int i=0; i<thisLevBoxes.size(); i++) {
	int gNodeIdx = thisLevBoxes[i];
	DblNumVec trgDwnEquDen_gNodeIdx((this->datSze(DE)));
	if (vlet()->periodic() || vlet()->dirichlet()){
	  if (vlet()->root(gNodeIdx)) { dcopy((this->matmgnt())->perdirmaps()->rtTrgDwnEquDen(), trgDwnEquDen_gNodeIdx); }
	}
	if (vlet()->tag(gNodeIdx) & LET_TRGNODE){
	  if(vlet()->depth(gNodeIdx) >= rt) {
	    //L2L - local -> local
	    iC( (this->matmgnt())->DwnChk2DwnEqu_dgemv(vlet()->depth(gNodeIdx)+(this->rootLevel()), this->trgDwnChkVal(gNodeIdx), trgDwnEquDen_gNodeIdx, 2.0) );
	  }
	  if (!vlet()->terminal(gNodeIdx)){
	    
	    for (int a = 0; a < 2; a++){ for (int b = 0; b < 2; b++){ for (int c = 0; c < 2; c++){
					     Index3 idx(a,b,c); int chi = vlet()->child(gNodeIdx, idx); iA( chi != -1 );
					     DblNumVec trgDwnChkVal_chi(this->trgDwnChkVal(chi));
					     iC( (this->matmgnt())->DwnEqu2DwnChk_dgemv(vlet()->depth(gNodeIdx)+(this->rootLevel()), idx, trgDwnEquDen_gNodeIdx, trgDwnChkVal_chi, 2.0) );
					   } } }
	    
	  }
	  else {
	    
	    DblNumVec grdExaVal_gNodeIdx(grdExaVal(gNodeIdx));
	    iC( DwnEqu2GrdChk_dgemv(vlet()->depth(gNodeIdx) + (this->rootLevel()), trgDwnEquDen_gNodeIdx, grdExaVal_gNodeIdx) );
	    
	    //iC( DwnEqu2GrdChk_dgemv(gNodeIdx, trgDwnEquDen_gNodeIdx, grdExaVal_gNodeIdx) );
	    //cerr<< gNodeIdx << " " << srcCoeffs(*vi) << " " << grdExaValgNodeIdx << endl;
		 //cerr << gNodeIdx << " " << grdExaVal_gNodeIdx << endl;


	  }
	}
      }
    }
  }
  
  double endTime = omp_get_wtime();
  fmmevaltime += (endTime - startTime);
  cout << "L2L + L2G used " << (endTime - startTime) << endl;
  (this->trgDwnChkVal()).resize(0);
  
  cerr << "TOTAL TIME = " << fmmevaltime << endl;
  
  return(0);
}


//PERIODIC CODE Evaluation Routines

/* for periodic or (zero) Dirichlet BCs, need to account for
 * flipped neighbors */
template <class VF>
int VFMM3d<VF>::evaluateUnodes_BCs(const int gNodeIdx){
  Point3 gNodeIdxCtr = vlet()->center(gNodeIdx);
  DblNumVec grdExaValgNodeIdx(grdExaVal(gNodeIdx));
  PerNode& curPerNode = vlet()->pernode(gNodeIdx);
  for(int a = 0; a < curPerNode.bdryUnodes().size(); a++) {
	 int U = curPerNode.bdryUnodes()[a].gni(); iA(vlet()->terminal(U));
	 Point3 offset = curPerNode.bdryUnodes()[a].bdryOffset();
	 /* Type of reflection - right now only dirichlet or periodic */
	 int type = (vlet()->dirichlet() ? curPerNode.bdryUnodes()[a].type() : FLN );
	 /* Apply *vi's coefficients to gNodeIdx */
	 DblNumVec coeffsvi;	coeffsvi.resize(vlet()->srcNk()*(this->srcDOF()));
	 /* Reflect coeffs as appropriate */
	 //cerr << vlet()->depth(gNodeIdx) << " " << vlet()->depth(U) << endl;
	 //cerr <<gNodeIdx << " " << U << " " << type << " " << offset << " " << srcCoeffs(U) << endl;
	 iC( reflectCoeffs(type, srcCoeffs(U), coeffsvi) );

	 /* evaluating U at gNodeIdx, so gni is the nbr */
	 int nbrTyp = nbrType(U, gNodeIdx);
	 if (abs(vlet()->depth(gNodeIdx) - vlet()->depth(U)) <= 1){
		int tabNbr;
		if (nbrTyp == NORM) { tabNbr = perNrmNbrType(U, gNodeIdx, offset, type); }
		else if (nbrTyp == FINE) { tabNbr = perFnNbrType(U, gNodeIdx, offset, type); }
		else if (nbrTyp == CRSE) { tabNbr = perCrsNbrType(U, gNodeIdx, offset, type); }
		else iA(0);
		GrdCoeffs2TrgVal_dgemv(vlet()->depth(gNodeIdx) + (this->rootLevel()), nbrTyp, tabNbr, coeffsvi, grdExaValgNodeIdx);
	 }	
	 else {
		iC( GrdCoeffs2UnbalTrgVal_dgemv(nbrTyp, gNodeIdx, U, coeffsvi, offset, grdExaValgNodeIdx));
		//iA(0); /* This needs to be fixed for cmp on fly */
		if (vlet()->dirichlet()) { iA(0); } // Impossible
		//Coeffs2TrgValCmptOnFly_dgemv(nbrTyp, U, gNodeIdx, coeffsvi, grdExaValgNodeIdx, offset);
	 }
  }
  return(0);
}

template <class VF>
int VFMM3d<VF>::evaluateVnodes_BCs(const int gNodeIdx, DblNumVec& effVal){
  PerNode& curPerNode = vlet()->pernode(gNodeIdx);
  Node& trgPtr = this->node(gNodeIdx);
  Point3 gNodeIdxCtr(vlet()->center(gNodeIdx));
  double D = 2.0 * vlet()->radius(gNodeIdx);
  DblNumVec trgDwnChkVal_gNodeIdx(this->trgDwnChkVal(gNodeIdx));
  int dep = vlet()->depth(gNodeIdx) + this->rootLevel();
  for(int a = 0; a < curPerNode.bdryVnodes().size(); a++) {
	 int v = curPerNode.bdryVnodes()[a].gni();
	 Point3 offset = curPerNode.bdryVnodes()[a].bdryOffset();
	 /* 0 == FLN */
	 int type = (vlet()->dirichlet() ? curPerNode.bdryVnodes()[a].type() : FLN );

	 Point3 vCtr = (vlet()->center(v));
	 if (vlet()->dirichlet()) vlet()->flpCtrDirNbr(vCtr, type);
	 for (int d = 0; d < (this->dim()); d++) { vCtr(d) += offset(d); }
	 Index3 idx;
	 for(int d=0; d<(this->dim()); d++){
		idx(d) = int(round( (vCtr[d]-gNodeIdxCtr[d])/D ));
	 }
	 /* The below section of code is repeated, so need to put it into a function */
	 Node& srcPtr = this->node(v);
	 //if (srcPtr.VotCnt()==0 && !vlet()->dirichlet() ) {
	 //srcPtr.effDen().resize( (this->matmgnt())->effDatSze(UE) ); setvalue(srcPtr.effDen(), 0.0);//1. resize effDen
	 //iC( (this->matmgnt())->plnDen2EffDen(dep, srcUpwEquDen(v),  srcPtr.effDen(), 2.0) );			 //2. transform from upeDen to effDen
	 //}
	 if (vlet()->dirichlet()){
		//NO#pragma omp critical
		DblNumVec srcEquDen(this->srcUpwEquDen(v));
		DblNumVec srcEquDenFlp; srcEquDenFlp.resize(srcEquDen.m());
		DblNumVec flpEffDen( (this->matmgnt())->effDatSze(UE) );
		DblNumMat& _DEM = ((this->matmgnt())->perdirmaps())->dwnEquGrdMap()[type];
		dgemv(1.0, _DEM, srcEquDen, 0.0, srcEquDenFlp);
		iC( (this->matmgnt())->plnDen2EffDen(dep, srcEquDenFlp,  flpEffDen, 2.0) );
		iC( (this->matmgnt())->UpwEqu2DwnChk_dgemv(dep, idx, flpEffDen, effVal, 2.0) );
		/*
		if (srcPtr.effDen(type).m() == 0){
		  DblNumVec srcEquDen(srcUpwEquDen(v));
		  DblNumVec srcEquDenFlp; srcEquDenFlp.resize(srcEquDen.m());
		  srcPtr.effDen(type).resize( (this->matmgnt())->effDatSze(UE) ); setvalue(srcPtr.effDen(J), 0.0);//1. resize effDen
		  DblNumMat& _DEM = ((this->matmgnt())->perdirmaps())->dwnEquGrdMap()[type];
		  dgemv(1.0, _DEM, srcEquDen, 0.0, srcEquDenFlp);
		  iC( (this->matmgnt())->plnDen2EffDen(dep, srcEquDenFlp,  srcPtr.effDen(type), 2.0) );
		  //2. transform from upeDen to effDen
		}
		*/
	 }
	 else {
		iC( (this->matmgnt())->UpwEqu2DwnChk_dgemv(dep, idx, srcPtr.effDen(), effVal, 2.0) );
	 }
	 //if(trgPtr.VinCnt()==0) {
	 //effVal.resize( (this->matmgnt())->effDatSze(DC) );	setvalue(effVal, 0.0);			 //1. resize effVal
	 //}
	 //M2L - multipole -> local
	 //iC( (this->matmgnt())->UpwEqu2DwnChk_dgemv(dep, idx, srcPtr.effDen(type), effVal, 2.0) );
	 //srcPtr.VotCnt()++;
	 //trgPtr.VinCnt()++;
	 //if (srcPtr.VotCnt()==srcPtr.VotNum()) {
	 //srcPtr.effDen().resize(0);			 //1. resize effDen to 0
	 //srcPtr.VotCnt()=0;
	 //}
	 //if (trgPtr.VinCnt()==trgPtr.VinNum()) {
	 //( (this->matmgnt())->effVal2PlnVal(effVal, trgDwnChkVal_gNodeIdx) );			 //1. transform from effVal to dncVal
	 //effVal.resize(0); //2. resize effVal to 0
	 //trgPtr.VinCnt()=0;
	 //}
  }
  return(0);
}

template <class VF>
int VFMM3d<VF>::evaluateWnodes_BCs(const int gNodeIdx){
  PerNode& curPerNode =  vlet()->pernode(gNodeIdx);
  DblNumVec grdExaVal_gNodeIdx(this->grdExaVal(gNodeIdx));
  Point3 gNodeIdxCtr = vlet()->center(gNodeIdx);
  for(int a = 0; a < curPerNode.bdryWnodes().size(); a++) {
	 int W = curPerNode.bdryWnodes()[a].gni();
	 Point3 offset = curPerNode.bdryWnodes()[a].bdryOffset();
	 int fliptype = (vlet()->dirichlet() ? curPerNode.bdryWnodes()[a].type() : FLN );
	 // Apply W's coefficients to gNodeIdx
	 if (vlet()->terminal(W) && abs(vlet()->depth(gNodeIdx) - vlet()->depth(W)) == 1) {
		DblNumVec coeffs; coeffs.resize(vlet()->srcNk());
		/* reflect coeffs as needed */
		iC( reflectCoeffs(fliptype, srcCoeffs(W), coeffs) );
		/* What type of X node is gNodeIdx relative to W */	
		int tabNode = perXnodeType(W, gNodeIdx, offset, fliptype);
		WXcoeffs2TrgVal(vlet()->depth(W) +(this->rootLevel()), XNODE, tabNode, coeffs, grdExaVal_gNodeIdx);
	 }
	 else {
		/* Apply *vi's coefficients to gNodeIdx - so using XNODE, not WNODE */
		int tabNode = _tbls->lookup(XNODE,vlet()->depth(gNodeIdx)+(this->rootLevel()), vlet()->radius(gNodeIdx), vlet()->center(gNodeIdx), vlet()->depth(W)+(this->rootLevel()), vlet()->radius(W), vlet()->center(W), offset);
		
		/* Apply *vi's coefficients to gNodeIdx - so using XNODE, not WNODE */
		iC( WXUpwEqu2DwnChk_dgemv(XNODE, W, gNodeIdx, (tabNode), offset, true, fliptype));
	 }
  }

  return(0);
}

template <class VF>
int VFMM3d<VF>::evaluateXnodes_BCs(const int gNodeIdx){

  PerNode& curPerNode = vlet()->pernode(gNodeIdx);
  DblNumVec grdExaVal_gNodeIdx(grdExaVal(gNodeIdx));
  Point3 gNodeIdxCtr = vlet()->center(gNodeIdx);			 
  for(int a = 0; a < curPerNode.bdryXnodes().size(); a++) {
	 int X = curPerNode.bdryXnodes()[a].gni(); iA(vlet()->terminal(X));
	 Point3 offset = curPerNode.bdryXnodes()[a].bdryOffset();
	 int fliptype = (vlet()->dirichlet() ? curPerNode.bdryXnodes()[a].type() : FLN );
	 if (vlet()->terminal(gNodeIdx) && abs(vlet()->depth(gNodeIdx) - vlet()->depth(X)) == 1) {
		// Apply X's coefficients to gNodeIdx
		DblNumVec coeffs; coeffs.resize(vlet()->srcNk());
		/* reflect coeffs as needed */
		iC( reflectCoeffs(fliptype, srcCoeffs(X), coeffs) );
		/* What type of W node is gNodeIdx relative to X */	
		int tabNode = perWnodeType(X, gNodeIdx, offset, fliptype);
		//int tabNode = 151 - perXnodeType(gNodeIdx, X, offset, fliptype);
		//cerr << tabNode << " " << srcCoeffs(X) << " " << coeffs << " " << grdExaVal_gNodeIdx << endl;
		WXcoeffs2TrgVal(vlet()->depth(X) + (this->rootLevel()), WNODE, tabNode, coeffs, grdExaVal_gNodeIdx);
		//cerr << grdExaVal_gNodeIdx << endl;
		//exit(0);
	 }
	 else { /* gni is not a leaf */
		int tabNode = _tbls->lookup(WNODE,vlet()->depth(X)+(this->rootLevel()), vlet()->radius(X), vlet()->center(X), vlet()->depth(gNodeIdx)+(this->rootLevel()), vlet()->radius(gNodeIdx), vlet()->center(gNodeIdx), offset);
		iC( WXUpwEqu2DwnChk_dgemv(WNODE, X, gNodeIdx, tabNode, offset, true, fliptype));
	 }
  }
  
  return(0);
}

template <class VF>
int VFMM3d<VF>::evaluate_far_BCs(const int gni, DblNumVec& FFtrgDwnEquDen){
  iA( gni == 0 );
  (this->matmgnt())->compPerDir() = true;
  //double D = 2.0 * vlet()->radius(); //2*radius of root
  DblNumVec srcUpwEquDen_gNodeIdx(this->srcUpwEquDen(gni));
  if (FFtrgDwnEquDen.m() == 0) { FFtrgDwnEquDen.resize((this->datSze(DE))); }
  
  vector<Node> _perLetVec;
  vector<Node> _perFMMsrcVec;
  vector<Node> _perFMMtrgVec;

  int rtlvl = this->rootLevel();
  /* For Dirichlet solver, correct for periodic part of the solver by moving root up one
	* Then, make the 8 necessary copies of upw ewn den, flipped as needed
	*/
  if (vlet()->dirichlet()) {
	 this->rootLevel() += (-1);
	 DblNumVec dirChkVal;
	 DblNumVec dirEquDen;
	 dirChkVal.resize((this->datSze(UC)));
	 dirEquDen.resize((this->datSze(UE)));
	 for(int a=0; a<2; a++) {
		for(int b=0; b<2; b++) {	
		  for(int c=0; c<2; c++) {
			 Index3 idx(a,b,c);
			 int n = idx(0)*4 + idx(1)*2 + idx(2);
			 iA( n >= 0 && n < 8);
			 DblNumVec srcUpwEquDenChi; srcUpwEquDenChi.resize(srcUpwEquDen_gNodeIdx.m());
			 DblNumMat& _DEM = ((this->matmgnt())->perdirmaps())->dwnEquGrdMap()[n];
			 dgemv(1.0, _DEM, srcUpwEquDen_gNodeIdx, 0.0, srcUpwEquDenChi);
			 /* Compute equ den on surface of super box */
			 iC( (this->matmgnt())->UpwEqu2UpwEqu_dgemv(vlet()->depth(gni)+(this->rootLevel()) + 1, idx, srcUpwEquDenChi, dirEquDen, 2.0) );
		  }
		}
	 }
	 /* Upw equ den to be used in periodic solver encompsses all flipped densities */
	 srcUpwEquDen_gNodeIdx = dirEquDen;
  }
  

  /* Sum of all densities should be zero on surface if total is zero in the domain as necessary */
  double sum = 0.0;
  for (int i = 0; i < srcUpwEquDen_gNodeIdx.m()/2; i++){
	 srcUpwEquDen_gNodeIdx(i) = -srcUpwEquDen_gNodeIdx(srcUpwEquDen_gNodeIdx.m()-1-i);
	 sum = sum + srcUpwEquDen_gNodeIdx(i) + srcUpwEquDen_gNodeIdx(srcUpwEquDen_gNodeIdx.m()-1-i);
  }
  
  /* This limit seems to work - need to do analytic testing */
  int LIM = (this->np()) + 1;

  /* Build increasing levels of upw equ den sities to tile the infinte domain with */
  DblNumVec _farSrcUpwChkVal; _farSrcUpwChkVal.resize((LIM) * (this->datSze(UC)));  setvalue(_farSrcUpwChkVal, 0.0);
  DblNumVec _farSrcUpwEquDen; _farSrcUpwEquDen.resize((LIM) * (this->datSze(UE)));  setvalue(_farSrcUpwEquDen, 0.0);
  DblNumVec _farSrcDwnChkVal; _farSrcDwnChkVal.resize((LIM) * (this->datSze(DC)));  setvalue(_farSrcDwnChkVal, 0.0);
  DblNumVec _farSrcDwnEquDen; _farSrcDwnEquDen.resize((LIM) * (this->datSze(DE)));  setvalue(_farSrcDwnEquDen, 0.0);
    
  _perLetVec.push_back( Node(1,-1,Index3(0,0,0), this->rootLevel()));

  /* Copy root level's data into our data structure */
  DblNumVec srcEquDen = DblNumVec((this->datSze(UE)), false, _farSrcUpwEquDen.data() + 0*(this->datSze(UE)));
  dcopy(srcUpwEquDen_gNodeIdx, srcEquDen);

  for (int d = 1; d < LIM; d++){
    _perLetVec.push_back( Node(-1,-1,Index3(0,0,0),this->rootLevel() - (d)));
    DblNumVec tmpChkVal = DblNumVec((this->datSze(UC)),	 false, _farSrcUpwChkVal.data() + d*(this->datSze(UC)));
    DblNumVec tmpEquDen = DblNumVec((this->datSze(UE)), false, _farSrcUpwEquDen.data() + d*(this->datSze(UE)));
    if (d > 0){ 
      DblNumVec tmpChiEquDen = DblNumVec((this->datSze(UE)), false, _farSrcUpwEquDen.data() + (d-1)*(this->datSze(UE)));
		/* Incorporate all 27 children as we are expanding in equal directions */
      for (int a = -1; a <= 1; a++){
		  for (int b = -1; b <= 1; b++){
			 for (int c = -1; c <= 1; c++){	
				Index3 idx(a,b,c);
				//iC( PerUpwEqu2UpwChk_dgemv(this->rootLevel() - (d) + 1, idx, tmpChiEquDen, tmpChkVal));
				iC( (this->matmgnt())->UpwEqu2UpwChk_dgemv(this->rootLevel() - (d) + 1, idx, tmpChiEquDen, tmpChkVal, PER_SCALE));
			 }
		  }
      }
    }
	 /* Compute equ density only at parent box at each level */
	 iC( (this->matmgnt())->UpwChk2UpwEqu_dgemv(this->rootLevel() - (d), tmpChkVal, tmpEquDen, PER_SCALE) );
	 //iC( PerUpwChk2UpwEqu_dgemv(this->rootLevel() - (d), tmpChkVal, tmpEquDen) );
	 /* Ensure zero sum to avoid instabilities */
    for (int i = 0; i < tmpEquDen.m()/2; i++){
		tmpEquDen(i) = -tmpEquDen(tmpEquDen.m()-1-i);	
	 }
  }

  _perFMMsrcVec.resize(_perLetVec.size());
  _perFMMtrgVec.resize(_perLetVec.size());

  /* Everything for the far fields and near fields are built to here thus far for all of the
   * "parents" */
  for (int d = LIM-1; d >= 0; d--){
	 /* Go out as far to left and rigt as necessary - incorporating all 316 possible V nodes, not just 189
	  * since we are expanding symmetrically */
    int begr = -4;
    int endr = 4;
    
    /* Equiv den from this level of node */
    srcEquDen = DblNumVec((this->datSze(UE)), false, _farSrcUpwEquDen.data() + d*(this->datSze(UE)));
    
    DblNumVec trgChkVal = DblNumVec((this->datSze(DC)), false, _farSrcDwnChkVal.data() + d*(this->datSze(DC)));
    Node& trgPtr = _perFMMtrgVec[d];
    Node& srcPtr = _perFMMsrcVec[d];
	 DblNumVec effVal( (this->matmgnt())->effDatSze(DC) );
    setvalue(effVal, 0.0);			 //1. resize effVal
    srcPtr.effDen().resize((this->matmgnt())->effDatSze(UE));
    setvalue(srcPtr.effDen(), 0.0);    
    
    /* 2. transform from upeDen to effDen
     * Have it stored at trgPtr since they are the same densities
     */
	 iC( (this->matmgnt())->plnDen2EffDen(this->rootLevel() - (d), srcEquDen,  srcPtr.effDen(), PER_SCALE) );
    for (int i = begr; i <= endr; i++){
      for (int j = begr; j <= endr; j++){
		  for (int k = begr; k <= endr; k++){
			 Index3 idx(i,j,k);
			 if (idx.linfty() > 1){ /* This is a Vnode */
				/* Map EffDen from Vnode at this level */
				/* Can use regular upwewu2DwnChk */
				iC( (this->matmgnt())->UpwEqu2DwnChk_dgemv(this->rootLevel() - (d), idx, srcPtr.effDen(), effVal, PER_SCALE));
			 }
		  }
		}
	 }

    iC( (this->matmgnt())->effVal2PlnVal(effVal, trgChkVal) );			 //1. transform from effVal to dncVal
    srcPtr.effDen().resize(0);
  }

  /* L2L from the infinite domain, passing parent's local to child.  Only needs to be
	* done at centered boxes */
  for (int d = LIM-1; d >=0; d--){
    if (d < LIM-1){
      Index3 chdIdx( _perLetVec[d+1].path2Node() - 2 * _perLetVec[d].path2Node() );

      DblNumVec trgDwnChkVal = DblNumVec((this->datSze(DC)), false, _farSrcDwnChkVal.data() + (d)*(this->datSze(DC)));
      DblNumVec srcDwnEquDen = DblNumVec((this->datSze(DE)), false, _farSrcDwnEquDen.data() + (d+1)*(this->datSze(DE)));
      //iC( PerDwnEqu2DwnChk_dgemv(this->rootLevel() - (d) - 1, Index3(0,0,0), srcDwnEquDen, trgDwnChkVal) );
		iC( (this->matmgnt())->DwnEqu2DwnChk_dgemv(this->rootLevel() - (d) - 1, Index3(0,0,0), srcDwnEquDen, trgDwnChkVal, PER_SCALE) );
	 }
    if (d < LIM){
      DblNumVec trgDwnEquDen = DblNumVec((this->datSze(DE)), false, _farSrcDwnEquDen.data() + d*(this->datSze(DE)));
      DblNumVec srcDwnChkVal = DblNumVec((this->datSze(DC)), false, _farSrcDwnChkVal.data() + d*(this->datSze(DC)));
      //iC( PerDwnChk2DwnEqu_dgemv(this->rootLevel() - (d), srcDwnChkVal, trgDwnEquDen) );
		iC( (this->matmgnt())->DwnChk2DwnEqu_dgemv(this->rootLevel() - (d), srcDwnChkVal, trgDwnEquDen, PER_SCALE) );
		/* Ensure stability */
      for (int i = 0; i < trgDwnEquDen.m()/2; i++){
		  trgDwnEquDen(i) = -trgDwnEquDen(trgDwnEquDen.m()-1-i);
		}
	 }
  }
  /* For main box */
  {
	 if ( (this->trgDwnChkVal()).m() == 0) (this->trgDwnChkVal()).resize(vlet()->trgNodeCnt() * (this->datSze(DC)));
	 DblNumVec trgDwnChkVal_gni = (this->trgDwnChkVal(gni));
	 DblNumVec trgDwnEquDen = DblNumVec((this->datSze(DE)), false, _farSrcDwnEquDen.data() + 0*(this->datSze(DE)));
	 DblNumVec trgDwnChkVal = DblNumVec((this->datSze(DC)), false, _farSrcDwnChkVal.data() + 0*(this->datSze(DC)));
	 /* For just periodic version, copy final dwn equ density to root of domain */
	 if (!vlet()->dirichlet()) {
		dcopy(trgDwnChkVal, trgDwnChkVal_gni);
		/* The next line is now taken care of in Vfmm3d_eval for rt == 0 */
		dcopy(trgDwnEquDen, FFtrgDwnEquDen);
	 }
	 else { /* Dirichlet or otherwise */
		/* gni is stored at (0,0,0) with respect to the super "root" */
		/* adjust should be negative 1 for super "root" */
		iC( (this->matmgnt())->DwnEqu2DwnChk_dgemv(vlet()->depth(gni) + (this->rootLevel()), Index3(0,0,0), trgDwnEquDen, trgDwnChkVal_gni, 2.0) );
		/* Change rootlevel back */
		this->rootLevel() = rtlvl;

		/* Evaluate root (not super-root)'s V-list */
		iC(evalPerDirVNodesRoot(gni));
		/* The next line is now taken care of in Vfmm3d_eval for rt == 0 */
		iC( (this->matmgnt())->DwnChk2DwnEqu_dgemv(vlet()->depth(gni) + (this->rootLevel()), trgDwnChkVal_gni, FFtrgDwnEquDen, 2.0) );
	 }
	 /* Make sure this is empty - not anymore */
	 setvalue(trgDwnChkVal_gni, 0.0);
	 /* Fix cleanup */
	 //iC( PerCleanup());
  }
  (this->matmgnt())->compPerDir() = false;

  return(0);
}

/* evaluate periodic dirichlet V Nodes
 * Have to take "reflections" into account
 * This is only done for the "Root" node
 */
template <class VF>
int VFMM3d<VF>::evalPerDirVNodesRoot(const int gni){
  iA(vlet()->root(gni));
  iA((this->matmgnt())->compPerDir() );
  cout << "evalPerDirVNode" << endl;
  DblNumVec trgDwnChkVal_gni = (this->trgDwnChkVal(gni));
  DblNumVec srcEquDen(this->srcUpwEquDen(gni));
    
  DblNumVec effVal((this->matmgnt())->effDatSze(DC));
  Node& srcPtr = this->node(gni);
    
  for (int n = 0; n < 8; n++){
	 DblNumVec srcEquDenFlp; srcEquDenFlp.resize(srcEquDen.m());
	 DblNumMat& _DEM = ((this->matmgnt())->perdirmaps())->dwnEquGrdMap()[n];
	 dgemv(1.0, _DEM, srcEquDen, 0.0, srcEquDenFlp);
	 DblNumVec effDen((this->matmgnt())->effDatSze(UE));
	 //setvalue(srcPtr.effDen(), 0.0);
	 iC( (this->matmgnt())->plnDen2EffDen(vlet()->depth(gni) + this->rootLevel(), srcEquDenFlp,  effDen, 3.0) );
	 //iC( PerPlnDen2EffDen(vlet()->depth(gni) + this->rootLevel(), srcEquDenFlp,  effDen) );
	 int I, J, K;
	 I = ((n == 0 || n == 1 || n == 2 || n == 3) ?  -2 : -1);
	 J = ((n == 0 || n == 1 || n == 4 || n == 5) ?  -2 : -1);
	 K = ((n == 0 || n == 2 || n == 4 || n == 6) ?  -2 : -1);
	 
	 for (int k = K; k <= K+4; k+=2){
		for (int i = I; i <= I+4; i+=2){
		  for (int j = J; j <= J+4; j+=2){
			 Index3 idx(i,j,k);
			 if (idx.linfty() > 1){
				/* Map Eff Den from Vnode at this level */
				/* Can use regular upwewu2DwnChk */
				iC( (this->matmgnt())->UpwEqu2DwnChk_dgemv(vlet()->depth(gni) + this->rootLevel(), idx, effDen, effVal, 2.0));
			 }
		  }
		}
	 }
  }

  iC( (this->matmgnt())->effVal2PlnVal(effVal, trgDwnChkVal_gni) );			 //1. transform from effVal to dncVal  
  
  return(0);

}


