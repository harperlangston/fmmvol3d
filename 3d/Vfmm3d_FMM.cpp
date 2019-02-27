#include "common/vecmatop.hpp"
#include "Vfmm3d_FMM.hpp"

#include <omp.h>
#include <time.h>

using namespace std;

VFMM3d_FMM::VFMM3d_FMM(const string& p):
  VFMM3d<VolNode>(p)
{
  
}

VFMM3d_FMM::~VFMM3d_FMM()
{
  
}

int VFMM3d_FMM::setup(){

  std::cerr << "m=  " << srcPos()->m() << endl;
  vlet()->ptsMax() = _fmm->let()->ptsMax();

  std::cerr << "pts max = " << vlet()->ptsMax() << " " << _fmm->let()->ptsMax() << std::endl;
  iC( vlet()->srcBalData() ); //iC( srcData() );
  iC( vlet()->trgData() ); //iC( trgData() );
  iC( vlet()->tblsSrcTrgData() );
  iC( vlet()->build() );
  iC( vlet()->trgData() ); //iC( trgData() );
  std::cerr << vlet()->trmNodeCnt() << std::endl;

  vector<int> ordVec; iC( let()->dwnOrderCollect(ordVec) ); //BOTTOM UP
  for (int i = 0; i < ordVec.size(); i++){
	 int gNodeIdx = ordVec[i];
	 //std::cerr << i << " " << ordVec[i] << " " << vlet()->terminal(ordVec[i]) << " " << vlet()->termidx(ordVec[i]) << " " << vlet()->termidx(ordVec[i-1]) << std::endl;
	 if (gNodeIdx >0 && _let->terminal(gNodeIdx) ) {
		if( ( vlet()->termidx(gNodeIdx) - 1 ) != vlet()->termidx(ordVec[i-1])) {
		  //std::cerr << i << " " << ordVec[i] << " " << vlet()->terminal(ordVec[i]) << " " << vlet()->termidx(ordVec[i]) << " " << vlet()->termidx(ordVec[i-1]) << std::endl;
		  //iA(0);
		}
	 }
  }
  
  /*
  DblNumMat grdSrcPos_FMM(3, vlet()->trmNodeCnt()*vlet()->srcGrdSze());
  //_grdSrcDen_FMM.resize(vlet()->trmNodeCnt()*vlet()->srcGrdSze()*srcDOF());

  for (int i = 0; i < ordVec.size(); i++){
	 int gNodeIdx = ordVec[i];
	 if (_let->terminal(gNodeIdx)){
		int tidx = vlet()->termidx(gNodeIdx);
		DblNumMat gpos(vlet()->grdSrcExaPos(gNodeIdx));
		
		for (int j = 0; j < vlet()->srcGrdSze(); j++){
		  for (int d = 0; d < dim(); d++){
			 grdSrcPos_FMM(d, tidx*vlet()->srcGrdSze() + j) = gpos(d, j);
		  }
		}
	 }
  }
  */

  //_fmm->trgPos() = &grdSrcPos_FMM; //Evaluating srcpos and srcden at vfmm grd pos
  //_fmm->trgVal() = &_grdSrcDen_FMM; //trg vals will be at grd pos


  iC( _fmm->setup() );

  iC( genNbrTypLsts());  
  if (vlet()->dirichlet() || vlet()->periodic()){
    iA (srcDOF() == 1 && trgDOF() == 1); //Only sdof=tdof=1 for now
    _matmgnt->pdMapsAlloc();
    iC( bldDirMaps() );
  }
  
  int grdExaCnt = vlet()->grdExaCnt();
  _grdExaVal.resize(grdExaCnt * trgDOF());
  
  _coeffs = new DblNumVec((vlet()->srcNk() * this->srcDOF()) * vlet()->trmNodeCnt());
  bldCoeffs(FRC);
  
  evaluate(*_fmm->srcDen(), *_fmm->trgVal());
  double relativeError;
  iC( fmm()->check(*_fmm->srcDen(), *_fmm->trgVal(), 20, relativeError) );
  std::cout << "relative error: " << relativeError << endl;

  double rerr, inferr;
  iC( vcheck(rerr, inferr) );
  cout << "relative error: " << rerr << endl; 

  /*
  iC( fmm()->evaluate(*_fmm->srcDen(), *_fmm->trgVal()) );
  relativeError = 0.0;
  iC( fmm()->check(*_fmm->srcDen(), *_fmm->trgVal(), 20, relativeError) );
  std::cout << "relative error: " << relativeError << endl;
  */

  exit(0);

  std::cerr << "let nodevec size " << _fmm->let()->nodeVec().size() << " " << vlet()->nodeVec().size() << std::endl;

  
  return(0);
}

// ---------------------------------------------------------------------- 
int VFMM3d_FMM::evaluate(const DblNumVec& srcDen, DblNumVec& trgVal)
{  
  
  bool _fmmcomp = false;
  bool _vfmmcomp = true;
  _matmgnt->report();
  //-----------------------------------
  cerr << srcDen.m() << " " << (*_srcPos).n() << endl;
  iA(srcDen.m()==srcDOF()*(*_srcPos).n());
  cerr << trgVal.m() << " " << (*_trgPos).n() << endl;
  iA(trgVal.m()==trgDOF()*(*_trgPos).n());

  cout << trgVal.m() << " " << (*_trgPos).n() << endl;
  //DblNumVec& tval = *(_trgVal);
  //-----------------------------------
  //iA(srcDen.m()==srcDOF()*(*_srcPos).n());
  int srcDOF = this->srcDOF();
  int trgDOF = this->trgDOF();

  int NK = vlet()->srcNk();

  double fmmevaltime = 0.0;
  
  //1. zero out Vecs
  setvalue(trgVal, 0.0);
  setvalue(_grdExaVal, 0.0);
  setvalue(_srcExaDen, 0.0);
  
  int srcNodeCnt = _let->srcNodeCnt();
  _srcUpwEquDen.resize(0);
  
  setvalue(_trgExaVal, 0.0);
  

  vector<int> ordVec; iC( let()->upwOrderCollect(ordVec) ); //BOTTOM UP

  // Make sure everything is marked - makes it simpler
  for(int i=0; i<ordVec.size(); i++) {
    int gNodeIdx = ordVec[i];
    //iA( _fmm->let()->tag(gNodeIdx) & LET_SRCNODE);
    //iA( _fmm->let()->tag(gNodeIdx) & LET_TRGNODE);
    //iA( this->let()->tag(gNodeIdx) & LET_SRCNODE);
    //iA( this->let()->tag(gNodeIdx) & LET_TRGNODE);
    _fmm->let()->tag(gNodeIdx) | LET_SRCNODE;
    _fmm->let()->tag(gNodeIdx) | LET_TRGNODE;
    this->let()->tag(gNodeIdx) | LET_SRCNODE;
    this->let()->tag(gNodeIdx) | LET_TRGNODE;
  }


  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if(_fmm->let()->tag(gNodeIdx) & LET_SRCNODE) {
	   iA(_let->tag(gNodeIdx) & LET_SRCNODE);
	   if(_fmm->let()->terminal(gNodeIdx)==true) {
		  DblNumVec srcExaDen(_fmm->srcExaDen(gNodeIdx));
		  vector<int>& curVecIdxs = _fmm->let()->node(gNodeIdx).srcOwnVecIdxs();
		  for(int k=0; k<curVecIdxs.size(); k++) {
			 int poff = curVecIdxs[k];
			 for(int d=0; d<srcDOF; d++) {
				srcExaDen(k*srcDOF+d) = srcDen(poff*srcDOF+d);
			 }
		  }
		}
	 }
  }

  /* Clear out the vector Idxs for non-leaves - saves memory */
  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if (!_let->terminal(gNodeIdx)){
		Node& curNode =  node(gNodeIdx);
		curNode.srcOwnVecIdxs().clear();
	 }
  }
  
  vector<int> trgTrmVec;
  for (int i = 0; i < ordVec.size(); i++){
	 int gni = ordVec[i];
	 
	 if (let()->terminal(gni) && let()->tag(gni) & LET_TRGNODE && let()->tag(gni) & LET_SRCNODE) { //evaluator
		trgTrmVec.push_back(gni);
	 }
  }

  cerr << "Number of source terminals " << trgTrmVec.size() << endl;

  map<int, vector<int> > lvlOrdVec; iC( let()->revLvlOrderCollect(lvlOrdVec) ); //BOTTOM UP
  vector<bool> lvsAtLev;
  vector<bool> ndsAtLev;
  for (int j = 0; j < lvlOrdVec.size(); j++){
	 lvsAtLev.push_back(false);
	 ndsAtLev.push_back(false);
	 vector<int>& thisLevBoxes = lvlOrdVec[j];
	 for (int i = 0; i < thisLevBoxes.size(); i++){
		int gNodeIdx = thisLevBoxes[i];
		if (let()->terminal(gNodeIdx) == true) { lvsAtLev[j] = true; }
		/* If in thi loop, there is a node at this level */
		ndsAtLev[j] = true;
	 }
  }
  
  cout << _knl.kernelType() << endl;
  //3. up computation
  if (1){
    cout << "S2M and M2M calculations" << endl;
    //Pre-build matrices for openmp
    int levnum = (_matmgnt->hom() ? 1 : lvlOrdVec.size());
    for (int j = 0; j < levnum; j++){
      vector<int>& thisLevBoxes = lvlOrdVec[j];
      DblNumVec chkVal(_matmgnt->plnDatSze(UC));
      DblNumVec tmpVec(datSze(UE));
      if (lvsAtLev[j] == true || _matmgnt->hom()){
	DblNumVec coeffs(vlet()->srcNk()); setvalue(coeffs,1.0);
	iC( GrdEqu2UpwEquTbls_dgemv(j+_rootLevel, coeffs, chkVal, tmpVec) );
      }
      if (ndsAtLev[j] == true || _matmgnt->hom()){
	for(int a=0; a<2; a++) { for(int b=0; b<2; b++) { for(int c=0; c<2; c++) {
				     Index3 idx(a,b,c);
				     iC( _matmgnt->UpwEqu2UpwChk_dgemv((j+1)+_rootLevel, idx, tmpVec, chkVal, 2.0) );
				   } } }
	iC( _matmgnt->UpwChk2UpwEqu_dgemv(j + _rootLevel, chkVal, tmpVec, 2.0) );
      }
    }
    
    double startTime = omp_get_wtime();
    for (int j = lvlOrdVec.size(); j >= 0; j--){
      vector<int>& thisLevBoxes = lvlOrdVec[j];
      if (ndsAtLev[j] == true){
#pragma omp parallel for
	for(int i=0; i<thisLevBoxes.size(); i++) {
	  int gNodeIdx = thisLevBoxes[i];
	  if(let()->tag(gNodeIdx) & LET_SRCNODE) {
	    DblNumVec srcUpwEquDengNodeIdx(_fmm->srcUpwEquDen(gNodeIdx));
	    DblNumVec chkVal(_matmgnt->plnDatSze(UC));
	    setvalue(chkVal, 0.0);
	    
	    if(let()->terminal(gNodeIdx)==true) {
	      //FMM S2M - Source -> Multipole Exapnsion
	      if (_fmmcomp){ iC( _fmm->SrcEqu2UpwChk_dgemv(_fmm->srcExaPos(gNodeIdx), _fmm->srcExaNor(gNodeIdx), _let->center(gNodeIdx), _let->radius(gNodeIdx), _fmm->srcExaDen(gNodeIdx), chkVal) ); }
	      //VFMM
	      if (_vfmmcomp) { iC( GrdEqu2UpwChkTbls_dgemv(_let->depth(gNodeIdx)+_rootLevel, srcCoeffs(gNodeIdx), chkVal) ); }
	    } else {				
	      //M2M - Multipole -> Multipole
	      //DblNumVec chkVal(_matmgnt->plnDatSze(UC));
	      for(int a=0; a<2; a++) {
		for(int b=0; b<2; b++) {
		  for(int c=0; c<2; c++) {
		    Index3 idx(a,b,c);
		    int chi = let()->child(gNodeIdx, idx);
		    if(let()->tag(chi) & LET_SRCNODE) {
		      DblNumVec srcUpwEquDenChi(_fmm->srcUpwEquDen(chi));
		      iC( _matmgnt->UpwEqu2UpwChk_dgemv(let()->depth(chi)+_rootLevel, idx, srcUpwEquDenChi, chkVal) );
		      //iC( _matmgnt->UpwEqu2UpwEqu_dgemv(let()->depth(chi)+_rootLevel, idx, srcUpwEquDenChi, srcUpwEquDengNodeIdx, 2.0) );
		    }
		  }		
		}
	      }
	    }
	    iC( _matmgnt->UpwChk2UpwEqu_dgemv(let()->depth(gNodeIdx) + _rootLevel, chkVal, srcUpwEquDengNodeIdx, 2.0) );
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
    if ((let()->terminal(gni) && _let->node(gni).Unodes().size() > 0) && let()->tag(gni) & LET_TRGNODE) { //evaluator) { //evaluator
      trgTrmVec.push_back(gni);
    }
  }
  
  ordVec.clear();  iC( let()->dwnOrderCollect(ordVec) ); //TOP DOWN
  //U - list contribution calculation
  if (_vfmmcomp){
	 int levnum = (_matmgnt->hom() ? 1 : lvlOrdVec.size());
	 for (int lev = 0; lev < levnum; lev++){
		vector<int>& thisLevBoxes = lvlOrdVec[lev];
		if (lvsAtLev[lev] == true || _matmgnt->hom()){
		  DblNumVec cfs(vlet()->srcNk());
		  DblNumVec vals(vlet()->trgGrdSze()*trgDOF);
		  for (int i = 0; i < 3; i++){
			 int nbrSze = (i == 0 ? 27 : 56);
			 for (int j = 0; j < nbrSze; j++){
				GrdCoeffs2TrgVal_dgemv(lev+_rootLevel, i, j, cfs, vals);
			 }
		  }
		}
		if (!(vlet()->balance())){
		  for(int i=0; i<trgTrmVec.size(); i++) {
			 int gNodeIdx = trgTrmVec[i];
			 Node& curNode = let()->node(gNodeIdx);
			 for(vector<int>::iterator vi=curNode.Unodes().begin(); vi!=curNode.Unodes().end(); vi++) {
				if (_let->tag(*vi) & LET_SRCNODE && _let->tag(gNodeIdx) * LET_TRGNODE){
				  int depdiff = abs(_let->depth(gNodeIdx) - _let->depth(*vi));
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
#pragma omp parallel for
    for(int i=0; i<trgTrmVec.size(); i++) {
      int gNodeIdx = trgTrmVec[i];
      if( let()->tag(gNodeIdx) & LET_TRGNODE) { 
	Node& curNode = let()->node(gNodeIdx);
	DblNumVec grdExaValgNodeIdx(grdExaVal(gNodeIdx));
	for(vector<int>::iterator vi=curNode.Unodes().begin(); vi!=curNode.Unodes().end(); vi++) {
	  if((let()->tag(*vi) & LET_SRCNODE)){
	    //S2T - source -> target
	    if (_fmmcomp)
	    {
	      DblNumVec trgExaValgNodeIdx(_fmm->trgExaVal(gNodeIdx));
	      DblNumMat trgExaPosgNodeIdx(_fmm->trgExaPos(gNodeIdx));
	      if (trgExaPosgNodeIdx.n() > 0){
		iC( _fmm->SrcEqu2TrgChk_dgemv(_fmm->srcExaPos(*vi), _fmm->srcExaNor(*vi), trgExaPosgNodeIdx, _fmm->srcExaDen(*vi), trgExaValgNodeIdx) );
	      }
	    }
	    //VFMM PART
	    if (_vfmmcomp)
	    {
	      DblNumVec coeffsvi(srcCoeffs(*vi));
	      int nbrTyp = nbrType(*vi, gNodeIdx);
	      /* Apply *vi's coefficients to gNodeIdx */
	      if (abs(_let->depth(gNodeIdx) - _let->depth(*vi)) <= 1){
		int tabNbr = tabNbrType(nbrTyp, *vi, gNodeIdx);
		iC( GrdCoeffs2TrgVal_dgemv(_let->depth(*vi)+_rootLevel, nbrTyp, tabNbr, coeffsvi, grdExaValgNodeIdx));
	      }
	      else {
		iC( GrdCoeffs2UnbalTrgVal_dgemv(nbrTyp, gNodeIdx, *vi, coeffsvi, grdExaValgNodeIdx));
	      }
	    }
	  }
	}
	if (_vfmmcomp) { if (vlet()->periodic()){ evaluateUnodes_BCs(gNodeIdx); } }
      }
    }
    double endTime = omp_get_wtime();
    fmmevaltime += (endTime - startTime);
    cout << "NEAR used " << (endTime - startTime) << endl;
  }
  
  if ( _trgDwnChkVal.m() == 0) _trgDwnChkVal.resize(0);
  setvalue(_trgDwnChkVal, 0.0);
  
  trgTrmVec.clear();
  for (int i = 0; i < ordVec.size(); i++){
    int gni = ordVec[i];
    if (_let->node(gni).Xnodes().size() > 0 || _let->node(gni).Wnodes().size() > 0) { //evaluator) { //evaluator
      trgTrmVec.push_back(gni);
    }
  }
  
  if (_vfmmcomp) {
    /* Build tables in full  - speedup for testing */
    int levnum = (_matmgnt->hom() ? 1 : lvlOrdVec.size());
    for (int lev = 0; lev < levnum; lev++){
      vector<int>& thisLevBoxes = lvlOrdVec[lev];
      if (lvsAtLev[lev] == true || _matmgnt->hom()){
	DblNumVec cfs(vlet()->srcNk());
	DblNumVec vals(vlet()->trgGrdSze()*trgDOF);
	for (int i = WNODE; i <= XNODE; i++){
	  int nbrSze = 152;
	  for (int j = 0; j < nbrSze; j++){
	    WXcoeffs2TrgVal(lev+_rootLevel, i, j, cfs, vals, false);
	  }
	}
      }
    }

    if (!(vlet()->balance()) || 1){
      for(int i=0; i<trgTrmVec.size(); i++) {
	int gNodeIdx = trgTrmVec[i];
	Node& curNode = let()->node(gNodeIdx);
	for(vector<int>::iterator vi=curNode.Wnodes().begin(); vi!=curNode.Wnodes().end(); vi++) {
	  bool dep = (abs(_let->depth(gNodeIdx) - _let->depth(*vi)) == 1);
	  if(!(let()->terminal(*vi) && dep)){
	    int tabNode = _tbls->lookup(XNODE,_let->depth(gNodeIdx)+_rootLevel, _let->radius(gNodeIdx), _let->center(gNodeIdx), _let->depth(*vi)+_rootLevel, _let->radius(*vi), _let->center(*vi));
	    /* Apply *vi's coefficients to gNodeIdx - so using XNODE, not WNODE */
	    iC( WXUpwEqu2DwnChk_dgemv(XNODE, *vi, gNodeIdx, (tabNode), false) );
	  }
	}
	for(vector<int>::iterator vi=curNode.Xnodes().begin(); vi!=curNode.Xnodes().end(); vi++) {
	  if((let()->tag(*vi) & LET_SRCNODE)){
	    bool dep = (abs(_let->depth(gNodeIdx) - _let->depth(*vi)) == 1);
	    if(!(let()->terminal(gNodeIdx) && dep)){	
	      int tabNode = _tbls->lookup(WNODE,_let->depth(*vi)+_rootLevel, _let->radius(*vi), _let->center(*vi), _let->depth(gNodeIdx)+_rootLevel, _let->radius(gNodeIdx), _let->center(gNodeIdx));
	      iC( WXUpwEqu2DwnChk_dgemv(WNODE, *vi, gNodeIdx, tabNode, false) );
	    }
	  }
	}
      }
    }
  }
  if (1){
    cout << "W/X - list contribution calculation" << endl;
    //W - list contribution calculation
    double startTime = omp_get_wtime();
    //#pragma omp parallel for
    for(int i=0; i<trgTrmVec.size(); i++) {
      int gNodeIdx = trgTrmVec[i];
      Node& curNode =  let()->node(gNodeIdx);
      DblNumVec trgExaVal_gNodeIdx(_fmm->trgExaVal(gNodeIdx));
      DblNumVec grdExaValgNodeIdx(this->grdExaVal(gNodeIdx));

      for(vector<int>::iterator vi=curNode.Wnodes().begin(); vi!=curNode.Wnodes().end(); vi++) {
	if((let()->tag(*vi) & LET_SRCNODE)){
	  //FMM 
	  if (_fmmcomp)
	  {
	    DblNumMat trgExaPosgNodeIdx(_fmm->trgExaPos(gNodeIdx));
	    if (trgExaPosgNodeIdx.n() > 0){

	      if(_fmm->let()->terminal(*vi) && _fmm->let()->node(*vi).srcExaNum()*srcDOF< _fmm->matmgnt()->plnDatSze(UE)) { //use Exa instead
		//S2T - source -> target
		iC( _fmm->SrcEqu2TrgChk_dgemv(_fmm->srcExaPos(*vi), _fmm->srcExaNor(*vi), trgExaPosgNodeIdx, _fmm->srcExaDen(*vi), trgExaVal_gNodeIdx) );
	      } else {
		//M2T - multipole -> target
		int vni = *vi;				
		iC( _fmm->UpwEqu2TrgChk_dgemv(_let->center(vni), _let->radius(vni), trgExaPosgNodeIdx, _fmm->srcUpwEquDen(*vi), trgExaVal_gNodeIdx) );
	      }
	    }
	  }
	  //VFMM
	  if (_vfmmcomp)
	  {
	    bool dep = (abs(_let->depth(gNodeIdx) - _let->depth(*vi)) == 1);
	    if(let()->terminal(*vi) && dep){
	      int tabNode = XNodeType(*vi,gNodeIdx);
	      iC( WXcoeffs2TrgVal(_let->depth(*vi)+_rootLevel, XNODE, (tabNode), srcCoeffs(*vi), grdExaValgNodeIdx));
	    }
	    else {
	      int tabNode = _tbls->lookup(XNODE,_let->depth(gNodeIdx)+_rootLevel, _let->radius(gNodeIdx), _let->center(gNodeIdx), _let->depth(*vi)+_rootLevel, _let->radius(*vi), _let->center(*vi));
	      /* Apply *vi's coefficients to gNodeIdx - so using XNODE, not WNODE */
	      if (_matmgnt->samPos(DC).n() < vlet()->grdTrgSamPos().n()){
		DblNumVec tdval(_fmm->trgDwnChkVal(gNodeIdx));
		iC( WXUpwEqu2DwnChk_dgemv(XNODE, *vi, gNodeIdx, (tabNode), _fmm->srcUpwEquDen(*vi), tdval) );
	      }
	      else {
		iC( WXUpwEqu2DwnChk_dgemv(XNODE, *vi, gNodeIdx, (tabNode), _fmm->srcUpwEquDen(*vi), grdExaValgNodeIdx) );
	      }
	    }
	  }
	}
      }
      if (_vfmmcomp) { if (vlet()->periodic()){ evaluateWnodes_BCs(gNodeIdx); }}
    }
    double endTime = omp_get_wtime();
    fmmevaltime += (endTime - startTime);
  }

  if (1) {
    cout << "W/X - list contribution calculation" << endl;
    //W - list contribution calculation
    double startTime = omp_get_wtime();
    //#pragma omp parallel for
    for(int i=0; i< trgTrmVec.size(); i++) {
      int gNodeIdx = trgTrmVec[i];
      Node& curNode =  let()->node(gNodeIdx);
      DblNumVec trgExaVal_gNodeIdx(_fmm->trgExaVal(gNodeIdx));
      DblNumVec trgDwnChkVal_gNodeIdx(_fmm->trgDwnChkVal(gNodeIdx));

      for(vector<int>::iterator vi=curNode.Xnodes().begin(); vi!=curNode.Xnodes().end(); vi++) {
	if((let()->tag(*vi) & LET_SRCNODE)){
	  //FMM
	  if (_fmmcomp)
	  {
	    if(_let->terminal(gNodeIdx) && _fmm->let()->node(gNodeIdx).trgExaNum()*trgDOF<_fmm->matmgnt()->plnDatSze(DC)) { //use Exa instead
	      DblNumMat trgExaPosgNodeIdx(_fmm->trgExaPos(gNodeIdx));
	      if (trgExaPosgNodeIdx.n() > 0){
		iC( _fmm->SrcEqu2TrgChk_dgemv(_fmm->srcExaPos(*vi), _fmm->srcExaNor(*vi), trgExaPosgNodeIdx, _fmm->srcExaDen(*vi), trgExaVal_gNodeIdx) );
	      }
	    } else {
	      //S2L - source -> local
	      iC( _fmm->SrcEqu2DwnChk_dgemv(_fmm->srcExaPos(*vi), _fmm->srcExaNor(*vi), _let->center(gNodeIdx), _let->radius(gNodeIdx), _fmm->srcExaDen(*vi), trgDwnChkVal_gNodeIdx) );
	    }
	  }
	  //VFMM
	  if (_vfmmcomp)
	  {
	    bool dep = (abs(_let->depth(gNodeIdx) - _let->depth(*vi)) == 1);
	    iA(_let->terminal(*vi));
	    if(let()->terminal(gNodeIdx) && dep){
	      int tabNode = WNodeType(*vi,gNodeIdx);
	      DblNumVec grdExaValgNodeIdx(this->grdExaVal(gNodeIdx));
	      DblNumVec grdExaValVi(this->grdExaVal(*vi));
	      iC( WXcoeffs2TrgVal(_let->depth(*vi)+_rootLevel, WNODE, (tabNode), srcCoeffs(*vi), grdExaValgNodeIdx));
	    }
	    else {
	      int tabNode = _tbls->lookup(WNODE,_let->depth(*vi)+_rootLevel, _let->radius(*vi), _let->center(*vi), _let->depth(gNodeIdx)+_rootLevel, _let->radius(gNodeIdx), _let->center(gNodeIdx));
	      //int tabNode = _tbls->lookup(WNODE,_let->depth(gNodeIdx)+_rootLevel, _let->radius(gNodeIdx), _let->center(gNodeIdx), _let->depth(*vi)+_rootLevel, _let->radius(*vi), _let->center(*vi));
	      DblNumVec tdval(_fmm->trgDwnChkVal(gNodeIdx));
	      iC( WXUpwEqu2DwnChk_dgemv(XNODE, *vi, gNodeIdx, (tabNode), _fmm->srcUpwEquDen(*vi), tdval) );
	    }
	  }
	}
      }
      if (_vfmmcomp){ if (vlet()->periodic()){ evaluateXnodes_BCs(gNodeIdx); } }
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
      int levnum = (_matmgnt->hom() ? 1 : lvlOrdVec.size());
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
	  if(1 || (let()->tag(gNodeIdx) & LET_SRCNODE)){
	    Node& srcPtr = node(gNodeIdx);
	    /* Need to find a good way to keep from recomputing this for dirichlet */
	    /* Seems to be taken care of */
	    srcPtr.effDen().resize( _matmgnt->effDatSze(UE) ); setvalue(srcPtr.effDen(), 0.0);//1. resize effDen
	    iC( _matmgnt->plnDen2EffDen(_let->depth(gNodeIdx)+_rootLevel, _fmm->srcUpwEquDen(gNodeIdx),  srcPtr.effDen(), 2.0) );			 //2. transform from upeDen to effDen
	  }
	}	
      }
      if (ndsAtLev[j] == true){
#pragma omp parallel for
	for(int i=0; i<thisLevBoxes.size(); i++) {
	  int gNodeIdx = thisLevBoxes[i];
	  if( let()->tag(gNodeIdx) & LET_TRGNODE) { //eValuator		//GNTra gnt = vlet()->gNodeIdx2gnt(gNodeIdx);
	    Point3 gNodeIdxCtr(let()->center(gNodeIdx));
	    double D = 2.0 * let()->radius(gNodeIdx);
	    DblNumVec trgDwnChkVal_gNodeIdx(_fmm->trgDwnChkVal(gNodeIdx));
	    
	    Node& trgPtr = node(gNodeIdx);
	    DblNumVec effVal(_matmgnt->effDatSze(DC));
	    
	    /* If periodic, evaluate these interactions */
	    if (vlet()->periodic()){ evaluateVnodes_BCs(gNodeIdx, effVal); } 
	    for(vector<int>::iterator vi=let()->node(gNodeIdx).Vnodes().begin(); vi!=let()->node(gNodeIdx).Vnodes().end(); vi++) {
	      if((let()->tag(*vi) & LET_SRCNODE)){
		Node& srcPtr = node(*vi);
		iA (srcPtr.effDen().m() != 0);
		Point3 viCtr(let()->center(*vi));
		Index3 idx;
		for(int d=0; d<dim(); d++){
		  idx(d) = int(round( (viCtr[d]-gNodeIdxCtr[d])/D ));
		}
		//M2L - multipole -> local
		iC( _matmgnt->UpwEqu2DwnChk_dgemv(_let->depth(gNodeIdx) + _rootLevel, idx, srcPtr.effDen(), effVal, 2.0) );
	      }
	    }
	    ( _matmgnt->effVal2PlnVal(effVal, trgDwnChkVal_gNodeIdx) );			 //1. transform from effVal to dncVal
	  }
	}
#pragma omp parallel for
	for(int i=0; i<thisLevBoxes.size(); i++) {
	  int gNodeIdx = thisLevBoxes[i];
	  Node& srcPtr = node(gNodeIdx);
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
    Node& Ptr = node(gNodeIdx);
    DblNumVec trgDwnChkVal_gNodeIdx(_fmm->trgDwnChkVal(gNodeIdx));
    //cerr << gNodeIdx << " " <<  trgDwnChkVal_gNodeIdx << endl;
    if (_let->terminal(gNodeIdx)){
      DblNumVec grdExaVal_gNodeIdx(grdExaVal(gNodeIdx));
      //cerr << gNodeIdx << " " <<  grdExaVal_gNodeIdx << endl;
    }
    iA (Ptr.effDen().m() == 0);
  }
  
  
  if ((vlet()->periodic() || (vlet()->dirichlet()))) {
    iA(_let->root(0));
    evaluate_far_BCs(0, _matmgnt->perdirmaps()->rtTrgDwnEquDen());
  }
  
  /*
    iC( _matmgnt->cleanPlans());
    if (_matmgnt->hom())
    iC( _matmgnt->UpwEqu2DwnChkCleanUp()); 
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
	DblNumVec chkVal(_matmgnt->plnDatSze(DC));
	DblNumVec denVal(datSze(DE));
	iC( _matmgnt->DwnChk2DwnEqu_dgemv(j+_rootLevel, chkVal, denVal, 2.0) );
	for(int a=0; a<2; a++) { for(int b=0; b<2; b++) { for(int c=0; c<2; c++) {
				     Index3 idx(a,b,c);
				     iC( _matmgnt->DwnEqu2DwnChk_dgemv(j+_rootLevel, idx, denVal, chkVal, 2.0) );
				   } } }
      }
      if (lvsAtLev[j] == true){
	DblNumVec denVal(datSze(DE)); setvalue(denVal,1.0); 
	DblNumVec grdVal(vlet()->trgGrdSze()); setvalue(grdVal,1.0); 
	iC( DwnEqu2GrdChk_dgemv(j+_rootLevel, denVal, grdVal, false) );
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
	DblNumVec trgDwnEquDen_gNodeIdx(datSze(DE));
	if (vlet()->periodic() || vlet()->dirichlet()){
	  if (_let->root(gNodeIdx)) { dcopy(_matmgnt->perdirmaps()->rtTrgDwnEquDen(), trgDwnEquDen_gNodeIdx); }
	}
	if (let()->tag(gNodeIdx) & LET_TRGNODE){
	  if(let()->depth(gNodeIdx) >= rt) {
	    //L2L - local -> local
	    iC( _matmgnt->DwnChk2DwnEqu_dgemv(let()->depth(gNodeIdx)+_rootLevel, _fmm->trgDwnChkVal(gNodeIdx), trgDwnEquDen_gNodeIdx, 2.0) );
	  }
	  if (!_let->terminal(gNodeIdx)){
	    for (int a = 0; a < 2; a++){ for (int b = 0; b < 2; b++){ for (int c = 0; c < 2; c++){
					     Index3 idx(a,b,c); int chi = _let->child(gNodeIdx, idx); iA( chi != -1 );
					     DblNumVec trgDwnChkVal_chi(_fmm->trgDwnChkVal(chi));
					     iC( _matmgnt->DwnEqu2DwnChk_dgemv(_let->depth(gNodeIdx)+_rootLevel, idx, trgDwnEquDen_gNodeIdx, trgDwnChkVal_chi, 2.0) );
					   } } }
	    
	  }
	  else {
	    //FMM
	    //L2T - local -> target
	    if (_fmmcomp){
	      DblNumVec trgExaVal_gNodeIdx(_fmm->trgExaVal(gNodeIdx));
	      DblNumMat trgExaPosgNodeIdx(_fmm->trgExaPos(gNodeIdx));
	      if (trgExaPosgNodeIdx.n() > 0){
		iC( _fmm->DwnEqu2TrgChk_dgemv(_let->center(gNodeIdx), _let->radius(gNodeIdx), trgExaPosgNodeIdx, trgDwnEquDen_gNodeIdx, trgExaVal_gNodeIdx) );
	      }
	    }
	    if (_vfmmcomp && !_fmmcomp) {
	      //VFMM - NOT DONE SINCE FAR FIELD COMPUTED AT TGTS
	      DblNumVec grdExaVal_gNodeIdx(grdExaVal(gNodeIdx));
	      iC( DwnEqu2GrdChk_dgemv(let()->depth(gNodeIdx) + _rootLevel, trgDwnEquDen_gNodeIdx, grdExaVal_gNodeIdx) );
	    }
	  }
	}
      }
    }
  }
  
  double endTime = omp_get_wtime();
  fmmevaltime += (endTime - startTime);
  cout << "L2L + L2G used " << (endTime - startTime) << endl;
  _trgDwnChkVal.resize(0);

  
  trgTrmVec.clear();
  for (int i = 0; i < ordVec.size(); i++){
    int gni = ordVec[i];
    if (let()->terminal(gni) && let()->tag(gni) & LET_TRGNODE) { //evaluator
      trgTrmVec.push_back(gni);
    }
  }

  //8. save trgExaVal
  //#pragma omp parallel for
  for(int i=0; i<trgTrmVec.size(); i++) {
    int gNodeIdx = trgTrmVec[i];
    DblNumVec trgExaVal(_fmm->trgExaVal(gNodeIdx));
    vector<int>& curVecIdxs = _fmm->let()->node(gNodeIdx).trgOwnVecIdxs();
    for(int k=0; k<curVecIdxs.size(); k++) {
      int poff = curVecIdxs[k];
      for(int d=0; d<trgDOF; d++) {
	trgVal(poff*trgDOF+d) = trgExaVal(k*trgDOF+d);
	//cerr << trgExaVal(k*trgDOF+d) << endl;
      }
    }
  }

  
  cerr << "TOTAL TIME = " << fmmevaltime << endl;
  
  return(0);
}

int VFMM3d_FMM::bldSrcCoeffsGrd_notused(const bool over, const int gNodeIdx){
  iA(vlet()->terminal(gNodeIdx));
  int NK = vlet()->srcNk();
  iA(!over);
 
  /* If the sources come from the forces on the GRID, then we are computing
	* a polynomial approximation to the force distribution.  Otherwise, we are computing
	* an approximation to the potential values on the GRID
	*/
  DblNumVec coeffVec(srcCoeffs(gNodeIdx));   setvalue(coeffVec,0.0);
  DblNumMat srcPos(vlet()->grdSrcExaPos(gNodeIdx)); /* Use oversampling to get good coeff approx */
  DblNumVec srcDen(grdExaDen(gNodeIdx));

  /* Need to get srcden over */
  
  if (srcDen.linfty() > 0){
	 double scale = (_matmgnt->hom()) ? pow(pow(0.5, (vlet()->depth(gNodeIdx)+_rootLevel)),2.0) : 1.0;
	 map<int, map<int, DblNumMat> >& pin_LevMap = (_matmgnt->hom()) ? _pinv[0] : _pinv[vlet()->depth(gNodeIdx)];
	 map<int, DblNumMat>& pin_map =  pin_LevMap[1];
	 DblNumMat& pinDOF = pin_map[pinType(NK)];
	 if (pinDOF.m() == 0){
		DblNumMat polys(srcPos.n(), NK);
		DblNumMat pin(NK, srcPos.n());
		/* Buil off of the node at this depth */
		iC( bldBasPolyMat(_matmgnt->hom(), vlet()->cheb(vlet()->kSrcVal()), FRC, gNodeIdx, NK, srcPos, polys));
		iC( pinv(polys, 1e-12, pin));
		pinDOF.resize(NK*srcDOF(), srcPos.n() * trgDOF());
		iC( bldDOFMat(pin, pinDOF, srcDOF(), trgDOF()));
	 }
	 DblNumVec resid(NK*srcDOF());
	 iC( dgemv(scale, pinDOF, srcDen, 1.0, resid));
	 for (int i = 0; i < NK*srcDOF(); i++) { 		coeffVec(i) = resid(i); }

	 bool all_zero = true;
	 for (int i = 0; i < NK*srcDOF(); i++){ if (abs(coeffVec(i)) != 0.0) all_zero = false; }
	 iA(!all_zero);
  }
  return(0);
}

DblNumVec VFMM3d_FMM::grdExaDen(int gNodeIdx)
{
  int tidx = vlet()->termidx(gNodeIdx);
  int num = vlet()->srcGrdSze()*srcDOF();
  return DblNumVec(num, false, _grdSrcDen_FMM.data()+tidx*num);
}
