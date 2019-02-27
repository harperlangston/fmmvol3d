#include "Vlet3d.hpp"

using std::istringstream;
using namespace std;

template <class V>
VLet3d<V>::VLet3d(const string& p)
  : Let3d<V>(p)
{
  
}

template <class V>
VLet3d<V>::~VLet3d()
{

}

template <class V>
bool VLet3d<V>::cheb(const int kval){
  if (kval <= 6) return false;
  else return true;
}

// ----------------------------------------------------------------------
template <class V>
int VLet3d<V>::setFromOptions(map<string,string>& optionsMap)
{
  //------------
  map<string,string>::iterator mapindex;
  mapindex = optionsMap.find("-" + this->prefix() + "maxLevel"); iA(mapindex!=optionsMap.end());
  { istringstream ss((*mapindex).second); ss>>this->_maxLevel; }

  mapindex = optionsMap.find("-" + this->prefix() + "periodic"); iA(mapindex!=optionsMap.end());
  { istringstream ss((*mapindex).second); ss>>_periodic; }
  mapindex = optionsMap.find("-" + this->prefix() + "dirichlet"); iA(mapindex!=optionsMap.end());
  { istringstream ss((*mapindex).second); ss>>_dirichlet; }

  mapindex = optionsMap.find("-" + this->prefix() + "ksrcval"); iA(mapindex!=optionsMap.end());
  { istringstream ss((*mapindex).second); ss>>_kSrcVal; }
  mapindex = optionsMap.find("-" + this->prefix() + "ktrgval"); iA(mapindex!=optionsMap.end());
  { istringstream ss((*mapindex).second); ss>>_kTrgVal; }

  mapindex = optionsMap.find("-" + this->prefix() + "rhs"); iA(mapindex!=optionsMap.end());
  { istringstream ss((*mapindex).second); ss>>_rhs; }
  
  mapindex = optionsMap.find("-" + this->prefix() + "balance"); iA(mapindex!=optionsMap.end());
  { istringstream ss((*mapindex).second); ss>>_balance; }

  mapindex = optionsMap.find("-" + this->prefix() + "adaptive"); iA(mapindex!=optionsMap.end());
  { istringstream ss((*mapindex).second); ss>>_adaptive; }
  return (0);
}

template <class V>
int VLet3d<V>::setup()
{
  srcNk() = kSrcVal()*(kSrcVal()+1)*(kSrcVal()+2)/6;
  trgNk() = kTrgVal()*(kTrgVal()+1)*(kTrgVal()+2)/6;
  cout << "PTS MAX = " << this->_ptsMax << endl;
  return (0);
}

template <class V>
int VLet3d<V>::build(){
  cout << "Balancing Tree" << endl;
  iC( srcBalData() );
  int i = 0;
  if (_balance){
	 _numViolators = -1;
	 while (_numViolators != 0){
		cout << "Pass " << i+1 << endl;	 
		iC( nbrsBalBld());
		iC( bal() );
		iC( subVltrs() );
		cout << "Number of remaining possible violators = " << _numViolators << endl;
		i++;
	 }
  }
  else {
	 iC( nbrsBalBld());
	 _numViolators = 0;
  }
  iC(srcTrgBalData());
  cout << "Tree Re-Balanced" << endl;
  

  return(0);
}


// ---------------------------------------------------------------------- 
/* For prebuilt tables, we have src and trg data as the same
 * Most of this is similar to the normal src and trg data setup.
 * Build the tree based on given information, then build the grids
 * and replace src and trg positions based on grid positions
 */
template <class V>
int VLet3d<V>::tblsSrcTrgData() {
  
  vector<V>& nodeVec = this->nodeVec();
  cout << "Number of Nodes = " << nodeVec.size() << endl;

  //ordering of the boxes, in top-down or bottom-up fashion
  vector<int> orderBoxesVec; iC( this->dwnOrderCollect(orderBoxesVec	) );

  /* Create Periodic or Dirichlet Nodes */
  if ((periodic() || dirichlet())){
	 cout << "PERIODIC/DIRICHLET CODE - VLET" << endl;
	 vector<PerNode>& perNodeVec = this->perNodeVec(); perNodeVec.clear();
	 for (int i = 0; i < orderBoxesVec.size(); i++){
		perNodeVec.push_back( PerNode() );		
	 }
	 for (int i = 0; i < orderBoxesVec.size(); i++){
		if (i == 0) perNodeVec[i].bdry() = true;
		if (!(this->terminal(i))){
		  for(int a=0; a<2; a++){
			 for(int b=0; b<2; b++){
				for(int c=0; c<2; c++){
				  setBdryLoc(i, ((a*2 + b)*2 + c));
				}
			 }
		  }
		}
		/* For debugging pruposes */
		iA( perNodeVec[i].bdryNbrs().size() == 0);
		iA( perNodeVec[i].bdryUnodes().size() == 0);
		iA( perNodeVec[i].bdryWnodes().size() == 0);
		iA( perNodeVec[i].bdryXnodes().size() == 0);
		iA( perNodeVec[i].bdryVnodes().size() == 0);
	 }

	 
	 
	 /* Build periodic or dirichlet nbr lists */
	 if ((periodic() && !dirichlet())){
		//Build periodic nbrs lists 
		bldPerNbrs();
		/* Build periodic V, U, W, and X lists */
		bldPerVnodes();
		bldPerUWXnodes();
	 }
	 else if (periodic() && dirichlet()){
		bldDirNbrs();
		//exit(0);
	 }
	 else {iA(0);}
  }


  return(0);
}

template <class V>
int VLet3d<V>::initTree(){
  bool regular = (!cheb(kSrcVal()));  
  bool edges = regular;

  /* recalc later */
  iC( grdSamPosCal(regular, kSrcVal(), 1.0, _grdSrcSamPos, false));
  iC( grdSamPosCal(regular, kSrcVal()*3, 1.0, _grdOverSamPos, false)); /* When doing grid tests, use a grid oversampling for coeff building */
  iC( grdSamPosCal(regular, kSrcVal()*2, 1.0, _grdDblSrcSamPos, false));
  

  regular = (!cheb(kTrgVal()));
  iC( grdSamPosCal(regular, kTrgVal(), 1.0, _grdTrgSamPos, false));
  
  vector<V>& nodeVec = this->nodeVec(); nodeVec.clear();
  nodeVec.push_back( V(-1,-1, Index3(0,0,0), 0) );
  iA( nodeVec.size() == 1);
  this->_level = 1;
  return(0);
}

/* This is very similar to srcData, but shorter for balancing
 * purposes since not all of the work in srcData() needs to
 * be done if balancing is taking place
 */
template <class V>
int VLet3d<V>::srcBalData(){

  DblNumMat& pos = *(this->_srcPos);
  cerr << "m = " << (this->_srcPos)->m() << endl;

  iA( pos.m()==this->dim() );

  vector<V>& nodeVec = this->nodeVec(); nodeVec.clear();
  vector< vector<int> > vecIdxs;  
  //local src number, the number of src point in each box
  vector<int> lclSrcNumVec; //glb 

  nodeVec.push_back( V(-1,-1, Index3(0,0,0), 0) );
  iA( nodeVec.size() == 1);
  vecIdxs.push_back( vector<int>() );
  vector<int>& curVecIdxs = vecIdxs[0];

  Point3 bbmin(this->center()-Point3(this->radius()));
  Point3 bbmax(this->center()+Point3(this->radius()));
  for(int k=0; k<pos.n(); k++) {
	 Point3 tmp(pos.clmdata(k));
	 iA(tmp>=bbmin && tmp<=bbmax);
	 curVecIdxs.push_back(k);
  }
  lclSrcNumVec.push_back( curVecIdxs.size() );
  
  int level = 0;
  int arrBeg = 0;
  int arrEnd = 1;
  int arrCnt = 0;
  while(arrBeg < arrEnd) {
	 //1.
	 arrCnt = arrEnd;
	 for(int k=arrBeg; k < arrEnd; k++) {
		//---
		if( lclSrcNumVec[k]> (this->ptsMax()) && level< (this->maxLevel())-1 ) {
		  nodeVec[k].child() = arrCnt;
		  arrCnt = arrCnt + pow2(this->dim());
		  //children's ess		  
		  for(int a=0; a<2; a++) {
			 for(int b=0; b<2; b++) {
				for(int c=0; c<2; c++) {
				  nodeVec.push_back( V(k,-1, 2*nodeVec[k].path2Node()+Index3(a,b,c), nodeVec[k].depth()+1) ); //par, chd
				  vecIdxs.push_back( vector<int>() );
				  lclSrcNumVec.push_back( 0 );
				}	
			 }
		  }
		  //children's vector of indices
		  Point3 centerCurNode( this->center(k) ); //get center of current node
		  for(vector<int>::iterator vecIdxsIt=vecIdxs[k].begin(); vecIdxsIt!=vecIdxs[k].end(); vecIdxsIt++) {
			 Point3 tmp(pos.clmdata(*vecIdxsIt));
			 Index3 idx;
			 for(int j=0; j<this->dim(); j++){
				idx(j) = (tmp(j) >= centerCurNode(j));
			 }
			 int chdGNodeIdx = this->child(k, idx);
			 vecIdxs[chdGNodeIdx].push_back(*vecIdxsIt);
		  }
		  vecIdxs[k].clear(); //VERY IMPORTANT
		  //children's lsm		  
		  for(int a=0; a<2; a++) {
			 for(int b=0; b<2; b++) {
				for(int c=0; c<2; c++) {
				  int chdGNodeIdx = this->child( k, Index3(a,b,c) );
				  lclSrcNumVec[chdGNodeIdx] = vecIdxs[chdGNodeIdx].size();
				}
			 }
		  }
		}
	 }
	 level++;
	 arrBeg = arrEnd;
	 arrEnd = arrCnt;
  }
  (this->_level) = level; //SET LEVEL

    /* Set Bdrys for periodic balancing */
  if ((periodic() || dirichlet())){
	 vector<int> ordVec; iC( this->dwnOrderCollect(ordVec	) );
	 vector<PerNode>& perNodeVec = this->perNodeVec(); perNodeVec.clear();
	 for (int i = 0; i < ordVec.size(); i++){
		perNodeVec.push_back( PerNode() );		
	 }
	 for (int i = 0; i < ordVec.size(); i++){
		if (i == 0) perNodeVec[i].bdry() = true;
		if (!(this->terminal(i))){
		  for(int a=0; a<2; a++){
			 for(int b=0; b<2; b++){
				for(int c=0; c<2; c++){
				  setBdryLoc(i, ((a*2 + b)*2 + c));
				}
			 }
		  }
		}
	 }
  }
  
  cout << "Max level of Tree = " << (this->_level) << " nodevec size = " << nodeVec.size() << endl;
  return(0);
}

template <class V>
int VLet3d<V>::nbrsBalBld(){
  
  vector<int> orderBoxesVec;
  iC( this->dwnOrderCollect(orderBoxesVec	) );

  vector<V>& nodeVec = this->nodeVec();
  for(int i=0; i < orderBoxesVec	.size(); i++) {
	 int gNodeIdx = orderBoxesVec[i];
	 nodeVec[gNodeIdx].Unodes().resize(0);
	 nodeVec[gNodeIdx].Wnodes().resize(0);
	 nodeVec[gNodeIdx].Xnodes().resize(0);
	 nodeVec[gNodeIdx].Vnodes().resize(0);
  }
  
  for(int g=0; g < orderBoxesVec.size(); g++) {
	 int gNodeIdx = orderBoxesVec[g];
	 nbrsBalBld(gNodeIdx);
  }

  /* Build periodic nbrs if needed - called for all bdry types except free-space */
  if (periodic() || dirichlet()){
	 bldPerNbrs(); /* Bld nbrs */
	 bldPerUnodes(); /* Bld just U Nodes */
  }
  
  return(0);
}

template <class V>
int VLet3d<V>::nbrsBalBld(int gNodeIdx){
  
  /* Only care about U nodes here, so only look at terminals */
  vector<V>& nodeVec = this->nodeVec();
  set<int> Uset;
  if(this->terminal(gNodeIdx) && gNodeIdx != 0) {
	 int parGNodeIdx = this->parent(gNodeIdx);
	 
	 Index3 minIdx(0);
	 Index3 maxIdx(pow2(this->depth(gNodeIdx)));
	
	 for(int i=-2; i<4; i++) {
		for(int j=-2; j<4; j++)	{
		  for(int k=-2; k<4; k++) {
			 Index3 tryPath( 2*this->path2Node(parGNodeIdx) + Index3(i,j,k) );
			 if(tryPath >= minIdx && tryPath <  maxIdx && tryPath != this->path2Node(gNodeIdx)) {	
				int resGNodeIdx = findgnt(this->depth(gNodeIdx), tryPath);
				if (this->adjacent(resGNodeIdx, gNodeIdx)){
				  if( this->depth(resGNodeIdx) < this->depth(gNodeIdx)) {
					 Uset.insert(resGNodeIdx);
				  }
				  if( this->depth(resGNodeIdx)==this->depth(gNodeIdx) ) {
					 queue<int> rest;
					 rest.push(resGNodeIdx);
					 while(rest.empty()==false) {
						int fntGNodeIdx = rest.front(); rest.pop();					 //int fntgNodeIdx = fntgnt.gNodeIdx();
						if(this->adjacent(fntGNodeIdx, gNodeIdx)) {
						  if(this->terminal(fntGNodeIdx)) {
							 Uset.insert(fntGNodeIdx);
						  }	
						  else { 
							 for(int a=0; a<2; a++) for(int b=0; b<2; b++) for(int c=0; c<2; c++) 
								rest.push( this->child(fntGNodeIdx, Index3(a,b,c)) );
						  }
						}
					 }
				  }
				}
			 }
		  }
		}
	 }
	 for(set<int>::iterator si=Uset.begin(); si!=Uset.end(); si++){
		nodeVec[gNodeIdx].Unodes().push_back(*si);
	 }
  }
  return(0);
}

template <class V>
int VLet3d<V>::bal(){
  
  vector<V>& nodeVec = this->nodeVec();
  vector<PerNode>& perNodeVec = this->perNodeVec();
  vector<int> orderBoxesVec	;  iC( this->dwnOrderCollect(orderBoxesVec	) );

  /* Go through whole tree only once and then just descendants of violators */
  bool descend = (_numViolators == -1 ? true : false);
  descend = true; 
  _numViolators = 0;
  {
	 /* Checking for primary violators */
	 for(int i=0; i < orderBoxesVec.size(); i++) {
		int gNodeIdx = orderBoxesVec[i];
		/* Set to false initially */
		nodeVec[gNodeIdx].prmVltr() = false;
		if ((descend || nodeVec[gNodeIdx].dscPrmVltr()) && this->terminal(gNodeIdx)){ /* Only care about leaves */
		  nodeVec[gNodeIdx].dscPrmVltr() = false; /* reset this */
		  V& curNode = this->node(gNodeIdx);
		  for(vector<int>::iterator vi=curNode.Unodes().begin(); vi!=curNode.Unodes().end(); vi++) {
			 //cout << *vi << endl;
		  
			 /* If *vi is more than one level deeper */
			 if ((this->depth(*vi) - this->depth(gNodeIdx)) > 1){
				nodeVec[gNodeIdx].prmVltr() = true;
				_numViolators++;
				//cout << "prmVltr = " << gNodeIdx << " " << this->depth(*vi) << " " <<this->depth(gNodeIdx) << endl;
				/* Once it is set, does not need to be set again */
				break;
			 }
		  }
		  if ((periodic() || dirichlet())){
			 /* If primary tag not yet set, check periodic nbrs */
			 if (!nodeVec[gNodeIdx].prmVltr()){
				PerNode curPerNode = pernode(gNodeIdx);
				for(int a = 0; a < curPerNode.bdryUnodes().size(); a++) {
				  int U = curPerNode.bdryUnodes()[a].gni();
				  if ((this->depth(U) - this->depth(gNodeIdx)) > 1){
					 nodeVec[gNodeIdx].prmVltr() = true;
					 _numViolators++;
					 /* Once it is set, does not need to be set again */
					 break;
				  }		
				}	
			 }
		  }
		}
	 }
  }

  /* secondary violators exist only once */
  if(descend){
	 /* Checking for secondary violators */
	 for(int i=0; i < orderBoxesVec.size(); i++) {
		int gNodeIdx = orderBoxesVec[i];
		/* Initially set to false */
		nodeVec[gNodeIdx].scdVltr() = false;
		/* Only care about leaves */
		if (this->terminal(gNodeIdx)){
		  /* Make sure it is not a prim violator */
		  if (!nodeVec[gNodeIdx].prmVltr()){
			 V& curNode = this->node(gNodeIdx);
			 for(vector<int>::iterator vi=curNode.Unodes().begin(); vi!=curNode.Unodes().end(); vi++) {
				/* If *vi is a primary violator AND finer */
				if (nodeVec[*vi].prmVltr() == true && (this->depth(*vi) > this->depth(gNodeIdx))){
				  nodeVec[gNodeIdx].scdVltr() = true;
				  _numViolators++;
				  /* Once it is set, does not need to be set again */
				  break;
				}
			 }
			 if ((periodic() || dirichlet())){
				/* If scdry tag not yet set, check  */
				if (!nodeVec[gNodeIdx].scdVltr()){
				  PerNode curPerNode = pernode(gNodeIdx);
				  for(int a = 0; a < curPerNode.bdryUnodes().size(); a++) {
					 int U = curPerNode.bdryUnodes()[a].gni();
					 if (nodeVec[U].prmVltr() == true && (this->depth(U) > this->depth(gNodeIdx))){
						nodeVec[gNodeIdx].scdVltr() = true;
						_numViolators++;
						/* Once it is set, does not need to be set again */
						break;
					 }		
				  }	
				}
			 }
		  }
		}
	 }
  }
  return(0);
}

template <class V>
int VLet3d<V>::subVltrs(){
  cout << "Number of Violators = " <<_numViolators << endl;
  vector<V>& nodeVec = this->nodeVec();
  vector<PerNode>& perNodeVec = this->perNodeVec();
  vector<int> orderBoxesVec	;  iC( this->dwnOrderCollect(orderBoxesVec	) );
  if (_numViolators != 0){
	 _numViolators = 0;
	 int arrCnt = nodeVec.size();
	 for (int i=0; i <  orderBoxesVec.size(); i++) {
		int gNodeIdx = orderBoxesVec[i];
		if (this->depth(gNodeIdx) > (this->_level)) { (this->_level) = this->depth(gNodeIdx); }
		if (nodeVec.size() < 50000000){ /* Seems to be a limit for seq ?*/
		  //cout << gNodeIdx << " " << this->terminal(gNodeIdx) << " " << (nodeVec[gNodeIdx].prmVltr()) << endl;
		  if ( (this->depth(gNodeIdx) + this->_rootLevel < this->maxLevel() - 1) && this->terminal(gNodeIdx) && (nodeVec[gNodeIdx].prmVltr() == true || nodeVec[gNodeIdx].scdVltr() == true)){
			 nodeVec[gNodeIdx].child() = arrCnt;
			 arrCnt = arrCnt + pow2(this->dim());
			 for(int a=0; a<2; a++) {
				for(int b=0; b<2; b++) {
				  for(int c=0; c<2; c++) {
					 nodeVec.push_back( V(gNodeIdx,-1, 2*nodeVec[gNodeIdx].path2Node()+Index3(a,b,c), nodeVec[gNodeIdx].depth()+1) ); //par, chd
					 int chdGNodeIdx = this->child( gNodeIdx, Index3(a,b,c));
					 if ((periodic() || dirichlet())){
						perNodeVec.push_back( PerNode() );
						setBdryLoc(gNodeIdx, ((a*2 + b)*2 + c));
					 }
					 if (nodeVec[gNodeIdx].prmVltr() == true){
						//cout << chdGNodeIdx << " " << this->terminal(chdGNodeIdx) << endl;
						/* Need to go deeper */
						nodeVec[chdGNodeIdx].dscPrmVltr() = true;
						_numViolators++;
					 }
				  }
				}
			 }
		  }
		  else { //Nothing
		  }
		}
	 }
  }
  return 0;
}
  
template <class V>
int VLet3d<V>::srcTrgBalData(){
  
  (this->_level) = 0; 
  DblNumMat& spos = *(this->_srcPos);  iA( spos.m()==this->dim() );

  vector<V>& nodeVec = this->nodeVec();
  vector< vector<int> > vecIdxs; vecIdxs.resize(nodeVec.size() );
  vector<int> lclSrcNumVec; lclSrcNumVec.resize(nodeVec.size(), 0);
  for (int i=0; i < nodeVec.size(); i++) {
	 vecIdxs.push_back( vector<int>() );
	 lclSrcNumVec.push_back( 0 );
  }
  /* We already know that the points or in th ebounding box */
  vector<int>& curVecIdxs = vecIdxs[0];
  for(int k=0; k < spos.n(); k++) {
	 curVecIdxs.push_back(k);
  }
  lclSrcNumVec[0] = curVecIdxs.size();
  vector<int> orderBoxesVec;
  iC( this->dwnOrderCollect(orderBoxesVec	) );
  for(int i=0; i < orderBoxesVec	.size(); i++) {
	 int gNodeIdx = orderBoxesVec[i];
	 V& curNode = nodeVec[gNodeIdx];
	 vector<int>& curVecIdxs = vecIdxs[gNodeIdx];	 
	 if(curNode.child()!=-1) { //not terminal
		//children's vecIdxs
		Point3 curCenter( this->center(gNodeIdx) );
		for(vector<int>::iterator curVecIdxsIt=curVecIdxs.begin();curVecIdxsIt !=curVecIdxs.end(); curVecIdxsIt++) {
		  Point3 tmp(spos.clmdata(*curVecIdxsIt));
		  Index3 idx;
		  for(int j=0; j<this->dim(); j++)
			 idx(j) = (tmp(j)>=curCenter(j));
		  int chdGNodeIdx = this->child(gNodeIdx, idx);
		  vector<int>& chdVecIdxs = vecIdxs[chdGNodeIdx];
		  chdVecIdxs.push_back(*curVecIdxsIt);
		}
		for(int a=0; a<2; a++) {
		  for(int b=0; b<2; b++) {
			 for(int c=0; c<2; c++) {
				int chdGNodeIdx = this->child(gNodeIdx, Index3(a,b,c));
				lclSrcNumVec[chdGNodeIdx] = vecIdxs[chdGNodeIdx].size();
			 }
		  }
		}
	 }
  }

  //#warning for now, only have targets at nonzero source locations
  int scnt = 0;
  int stcnt = 0;
  int tcnt = 0;
  map<int, vector<int> > lvlOrdVec; iC( this->revLvlOrderCollect(lvlOrdVec) ); //BOTTOM UP
  for (int j = lvlOrdVec.size(); j >= 0; j--){
	 vector<int>& thisLevBoxes = lvlOrdVec[j];
	 for (int i = 0; i < thisLevBoxes.size(); i++){
		int gNodeIdx = thisLevBoxes[i];
		nodeVec[gNodeIdx].tag() = 0;
		
		if (this->terminal(gNodeIdx) == true){
		  DblNumMat srcPos(grdSrcExaPos(gNodeIdx));
		  DblNumVec srcDen(grdSrcSamPos().n() * _knl.srcDOF());
		  exsol3d()->quantity(QNT_RHS, srcPos, srcDen);
		  if (srcDen.linfty() > 0.0 || exsol3d()->ct() == CHS_EMPTY) {
		    nodeVec[gNodeIdx].tag() = nodeVec[gNodeIdx].tag() | LET_SRCNODE;
		    nodeVec[gNodeIdx].tag() = nodeVec[gNodeIdx].tag() | LET_TRGNODE;
		    scnt++;
		    stcnt++;
		  }
		}
		else {
		  if (!(nodeVec[gNodeIdx].tag() & LET_SRCNODE)){
			 for(int a=0; a<2; a++) {
				for(int b=0; b<2; b++) {
				  for(int c=0; c<2; c++) {
					 Index3 idx(a,b,c);
					 int chi = this->child(gNodeIdx, idx);
					 if (nodeVec[chi].tag() & LET_SRCNODE){
						nodeVec[gNodeIdx].tag() = nodeVec[gNodeIdx].tag() | LET_SRCNODE;
						nodeVec[gNodeIdx].tag() = nodeVec[gNodeIdx].tag() | LET_TRGNODE;
						scnt++;
					 }
				  }
				}
			 }
		  }
		  else { cerr << "WRONG" << endl; exit(0); }
		}
	 }	
  }
  this->_srcNodeCnt = scnt;
  //_trgNodeCnt = scnt;

    /* For now all terminals have grids, and assign terminal index */
  int t = 0;
  _trmnodecnt = 0;
  for (int i = 0; i < orderBoxesVec.size(); i++){
	 //if (this->terminal(i) && nodeVec[i].tag() & LET_SRCNODE){
	 if (this->terminal(i)){
      /* add index to list of terminals if tables used */
		_trmnodecnt++;
		nodeVec[i].termidx() = t; t++;
	 }
	 else {
		//
	 }
  }

  cout << "Number of Leaves = " << _trmnodecnt << endl;
  cout << "Number of SRC NODES = " << scnt << " " << stcnt << endl;

  //No longer needed
  _grdDblSrcSamPos.resize(0,0);
  
  for(int i=0; i < orderBoxesVec	.size(); i++) {
	 int gNodeIdx = orderBoxesVec[i];
	 nodeVec[gNodeIdx].Unodes().resize(0);
  }

  /*
  for(int i=0; i < orderBoxesVec.size(); i++) {
	 int gNodeIdx = orderBoxesVec[i];
2A	 if(nodeVec[gNodeIdx].tag() & LET_TRGNODE) { //a evtr		
		iC( calgnext(gNodeIdx) );
	 }
  }
  */
  
  (this->_level)++; //level is max this->depth + 1
  cout << "New Max level of Tree = " << (this->_level) << endl;

  return(0);
}

/*******************************************
 * The following are all for handling periodic boundary conditions,
 * dirichlet bdry condistions, etc.
 */
template <class V>
int VLet3d<V>::setBdryLoc(int gni, int chiNum){

  iA(periodic());
  vector<V>& nodeVec = this->nodeVec();
  vector<PerNode>& perNodeVec = this->perNodeVec();
  Point3 ctrRt = this->center();
  /* For a node gni, check to see if it children are also at the boundary */
  if (perNodeVec[gni].bdry()){
	 int chi = nodeVec[gni].child() + chiNum;
	 Point3 ctrChi(this->center(chi));
	 Point3 ctr(this->center(chi));
	 /* For level root + 1, all 8 children are at the boundary,
	  * so we set these manually since root's boundary location
	  * is all boundaries, and that doesn't work well with this
	  * data structure
	  */
	 if (nodeVec[gni].depth()+ this->rootLevel() + 1 == this->rootLevel()+1) {
		perNodeVec[chi].bdry() = true;
		Index3 pDif = this->path2Node(chi) - 2*this->path2Node(gni);
		Point3 bloc;
		for (int i = 0; i < 3; i++){
		  if (pDif(i) == 0) bloc(i) = ctrRt(i) - this->radius();
		  else bloc(i) = ctrRt(i) + this->radius();
		}
		perNodeVec[chi].bdryLoc() = bloc;
	 }
	 else {
		Point3 pbLoc = perNodeVec[gni].bdryLoc(); /* Which boundary gni touches */
		Point3 chibLoc = Point3(-100000.0,-100000.0,-100000.0); /* Which boundary child touches as of yet */
		Index3 pDif = this->path2Node(chi) - 2*this->path2Node(gni);
		Point3 bloc;
		for (int i = 0; i < 3; i++){
		  if (pDif(i) == 0) bloc(i) = ctrChi(i) - this->radius(chi);
		  else bloc(i) = ctrChi(i) + this->radius(chi);
		}
		for (int m = 0; m < 3; m++){
		  if (bloc(m) == pbLoc(m)) chibLoc(m) = pbLoc(m);
		}
		perNodeVec[chi].bdryLoc() = chibLoc;
		/* If the resulting child boundary location is non-zero anywhere,
		 * it touches the boundary, so set it */
		for (int i = 0; i < 3; i++){
		  if (chibLoc(i) == (ctrRt(i) - this->radius()) || chibLoc(i) == (ctrRt(i) + this->radius())){
			 perNodeVec[chi].bdry() = true;
		  }
		}
	 }
  }
  return(0);
}

/* Checking for same level neighbors.  When we build the lists, we only look
 * at a node's a parent's neighbors and their children.  By using inverse
 * relationships, the X and Crse lists will be built at the same time as the W
 * and Fine lists
 */
template <class V>
int VLet3d<V>::bldPerNbrs(){
  
  vector<int> orderBoxesVec;
  vector<V>& nodeVec = this->nodeVec();
  vector<PerNode>& perNodeVec = this->perNodeVec();
  int gni;
  iC( this->dwnOrderCollect(orderBoxesVec	) );
  for (int i = 1; i < orderBoxesVec.size(); i++){
	 gni = orderBoxesVec[i];
	 perNodeVec[gni].bdryNbrs().clear();
  }
  /* Root is at zero */
  gni =  orderBoxesVec[0];
  Point3 ctr(this->center());
  /* For the root node, create the boundary neighbors since
	* we know where they are */
  for (int a = -1; a <= 1; a++){
	 for (int b = -1; b <= 1; b++){
		for (int c = -1; c <= 1; c++){
		  if (!(a == 0 && b == 0 && c == 0)){
			 double drad = 2.0*this->radius();
			 Point3 off(drad*a, drad*b, drad*c);
			 perNodeVec[gni].bdryNbrs().push_back(BdryNode(gni, off));
		  }
		}
	 }
  }

  /* For nodes starting at the root's children and down, do the following:
	* - Check to see if the node is on the boundary (already set)
	* - Look at a node's parent's periodic neighbors
	* - If these periodic parental neighbors have children,
	* see if these children are adjacent once they are offset the
	* proper amount.  If so, they are periodic neighbors of the
	* current box with offset the same as their parents
	*/
  for (int i = 1; i < orderBoxesVec.size(); i++){
	 gni = orderBoxesVec[i];
	 if ((nodeVec[gni].depth() + this->rootLevel() > this->rootLevel()) && perNodeVec[gni].bdry()){
		int par = nodeVec[gni].parent(); /* This node's parent */
		double parRad = this->radius(par); /* This node's parent's radius */
		double gniRad = this->radius(gni); /* This node's radius */

		for (int j = 0; j < perNodeVec[par].bdryNbrs().size(); j++){
		  /* parNbrOff tells us what we would have to add to the
			* center of this node's parent's j'th periodic nbr
			* to put it where it belongs periodically
			*/
		  Point3 parNbrOff = perNodeVec[par].bdryNbrs()[j].bdryOffset();
		  int parNbr = perNodeVec[par].bdryNbrs()[j].gni(); /* parent's j'th periodic nbr's global index */
		  if (!nodeVec[parNbr].terminal()){
			 /* Go through parent's periodic neighbor's children */
			 for (int a = 0; a < 2; a++){
				for (int b = 0; b < 2; b++){
				  for (int c = 0; c < 2; c++){
					 int parNbrChi = this->child(parNbr, Index3(a,b,c));
					 Point3 parNbrChiCtr = this->center(parNbrChi); /* par's jth per nbr's center */
					 for (int d = 0; d < 3; d++) { parNbrChiCtr(d) += (parNbrOff(d)); } /* Move by offset */
					 Point3 dif = parNbrChiCtr - this->center(gni);
					 if (dif.linfty() == parRad){ /* kth child is adjacent */
						perNodeVec[gni].bdryNbrs().push_back(BdryNode(parNbrChi, parNbrOff));
					 }
				  }
				}
			 }
		  }
		  else {
			 //Point3 parNbrCtr = this->center(parNbr) + parNbrOff;
		  }
		}
	 }
  }
  //exit(0);
  return(0);
}

template <class V>
int VLet3d<V>::bldPerVnodes(){
  iA( !dirichlet() );
  vector<int> orderBoxesVec;
  vector<V>& nodeVec = this->nodeVec();
  vector<PerNode>& perNodeVec = this->perNodeVec();
  iC( this->dwnOrderCollect(orderBoxesVec	) );

  for (int i = 1; i < orderBoxesVec.size(); i++){
	 int gni = orderBoxesVec[i];
	 int par = nodeVec[gni].parent(); /* This node's parent */
	 /* If this node's parent touches the boundary */
	 if ((nodeVec[gni].depth() + this->rootLevel() > this->rootLevel()) && perNodeVec[par].bdry()){
		double parRad = this->radius(par); /* This node's parent's radius */
		//double gniRad = this->radius(gni); /* This node's radius */
		for (int j = 0; j < perNodeVec[par].bdryNbrs().size(); j++){
		  Point3 parNbrOff = perNodeVec[par].bdryNbrs()[j].bdryOffset();
		  /* parent's j'th periodic nbr's global index */
		  int parNbr = perNodeVec[par].bdryNbrs()[j].gni(); 
		  if (!nodeVec[parNbr].terminal()){
			 /* Go through parent's periodic neighbor's children */
			 for (int a = 0; a < 2; a++){
				for (int b = 0; b < 2; b++){
				  for (int c = 0; c < 2; c++){
					 int parNbrChi = this->child(parNbr, Index3(a,b,c)); /* kth child of par's per' jth nbr */
					 Point3 parNbrChiCtr(this->center(parNbrChi) + parNbrOff); /* par's jth per nbr's center */
					 Point3 dif = parNbrChiCtr - this->center(gni);
					 if (dif.linfty() > parRad){ /* kth child is NOT adjacent */
						iA( dif.linfty() >= 2.0*parRad);
						if (nodeVec[parNbrChi].tag() & LET_SRCNODE){	
						  perNodeVec[gni].bdryVnodes().push_back(BdryNode(parNbrChi, parNbrOff));
						}
						//cout << "VNODE " << gni << " " << j << " " << parNbrChi << " " <<  parNbrOff << parNbrChiCtr << endl;
					 }
				  }
				}				 
			 }
		  }
		}
	 }
  }
  return(0);
}

/* This function is for dirichlet and periodic */
template <class V>
int VLet3d<V>::bldPerUWXnodes(){
  
  vector<int> orderBoxesVec;
  vector<V>& nodeVec = this->nodeVec();
  vector<PerNode>& perNodeVec = this->perNodeVec();
  
  iC( this->dwnOrderCollect(orderBoxesVec	) );

  for (int i = 1; i < orderBoxesVec.size(); i++){
	 int gni = orderBoxesVec[i];
	 double gniRad = this->radius(gni); /* This node's radius */
	 /* If this node's parent touches the boundary */
	 if (nodeVec[gni].terminal() && perNodeVec[gni].bdry()){
		//double parRad = this->radius(par); /* This node's parent's radius */
		for (int j = 0; j < perNodeVec[gni].bdryNbrs().size(); j++){
		  Point3 nbrOff = perNodeVec[gni].bdryNbrs()[j].bdryOffset();
		  /* parent's j'th periodic nbr's global index */
		  int nbr = perNodeVec[gni].bdryNbrs()[j].gni();
		  int type = ( dirichlet() ? perNodeVec[gni].bdryNbrs()[j].type() : FLN);
		  Point3 nbrCtr = this->center(nbr);
		  Point3 gniCtr = this->center(gni);
		  if (dirichlet() && type != FLN) flpCtrDirNbr(nbrCtr, type); /* "Flip" as necessary */
		  for (int a = 0; a < 3; a++) { nbrCtr(a) += (nbrOff(a)); }
		  Point3 dif = nbrCtr - gniCtr;
		  if (dif.linfty() == 2.0*gniRad){ /* nbr is adjacent */
			 if (nodeVec[nbr].terminal()){ /* Nbr is a leaf at same level */
				iA(this->depth(gni) == this->depth(nbr));
				if (nodeVec[nbr].tag() & LET_SRCNODE){
				  perNodeVec[gni].bdryUnodes().push_back(BdryNode(nbr, nbrOff, type));
				}
			 }
			 else { /* Descendants of nbr Possibly a fine neighbor or W node */
				/* Need to buidl descendants and then check them */
				if (perNodeVec[nbr].dscList().size() == 0){
				  bldDscList(nbr);
				}
				for (int k = 0; k < perNodeVec[nbr].dscList().size(); k++){
				  int nbrDsc = perNodeVec[nbr].dscList()[k];
				  Point3 nbrDscCtr = this->center(nbrDsc);
				  //cerr << nbrDscCtr << endl;
				  /* "Flip" as necessary */
				  if (type != FLN) flpCtrDirNbr(nbrDscCtr, type);

				  Point3 difDsc = (nbrDscCtr + nbrOff) - this->center(gni);
				  int levDif = this->depth(nbrDsc) - this->depth(gni);
				  /*
					 cout << i << " " << nbr <<  " " << nbrDsc << " " << parent(nbrDsc) << endl;
					 cerr << this->depth(gni) << " " << this->depth(nbr) << " " << this->depth(nbrDsc) << endl;
					 cerr << this->center(gni) << " " << this->center(nbr) << " " << this->center(nbrDsc) << endl;
					 cerr << nbrCtr << " " << nbrDscCtr << endl;
					 cout << levDif << " " << dif << " " << type << endl;
					 cout << nbrOff << endl;
				  */
				  //if (dirichlet()) {  iA(levDif == 1); }
				  if (this->terminal(nbrDsc)){
					 if (difDsc.linfty() == gniRad*(1.0 + pow(0.5, levDif))){
						if (nodeVec[nbrDsc].tag() & LET_SRCNODE){
						  /* Fine Nbr */
						  perNodeVec[gni].bdryUnodes().push_back(BdryNode(nbrDsc, nbrOff, type));
						}
						/* Coarse Nbr */
						if (dirichlet()) {
						  if (nodeVec[gni].tag() & LET_SRCNODE){
							 perNodeVec[nbrDsc].bdryUnodes().push_back(BdryNode(gni, nbrOff, type));
						  }
						}
						else {
						  if (nodeVec[gni].tag() & LET_SRCNODE){
							 perNodeVec[nbrDsc].bdryUnodes().push_back(BdryNode(gni, -nbrOff, type));
						  }
						}
					 }
					 else if (difDsc.linfty() == gniRad*(1.0 + 3.0*pow(0.5, levDif))){
						/* W node */
						//if (nodeVec[nbrDsc].tag() & LET_SRCNODE){
						perNodeVec[gni].bdryWnodes().push_back(BdryNode(nbrDsc, nbrOff, type));
						//}
						/* X node */
						if (dirichlet()) {
						  //if (nodeVec[gni].tag() & LET_SRCNODE){
						  perNodeVec[nbrDsc].bdryXnodes().push_back(BdryNode(gni, nbrOff, type));
						  //}
						}
						else {
						  if (nodeVec[gni].tag() & LET_SRCNODE){	
							 perNodeVec[nbrDsc].bdryXnodes().push_back(BdryNode(gni, -nbrOff, type));
						  }
						}
					 }
				  }
				  else {
					 if (difDsc.linfty() == gniRad*(1.0 + 3.0*pow(0.5, levDif))){
						/* W node */
						//if (nodeVec[nbrDsc].tag() & LET_SRCNODE){	
						perNodeVec[gni].bdryWnodes().push_back(BdryNode(nbrDsc, nbrOff, type));
						//}
						/* X node  - gni must be a terminal */
						if (dirichlet()) {
						  //if (nodeVec[gni].tag() & LET_SRCNODE){	
						  perNodeVec[nbrDsc].bdryXnodes().push_back(BdryNode(gni, nbrOff, type));
						  //}
						}
						else {
						  if (nodeVec[gni].tag() & LET_SRCNODE){	
							 perNodeVec[nbrDsc].bdryXnodes().push_back(BdryNode(gni, -nbrOff, type));
						  }
						}
					 }
				  }
				}
			 } /* END ELSE */
		  }
		}
	 }
  }
  return(0);
}

/* Just U Nodes - For balancing purposes across periodic bdrys
 * This is a lot like above, absent the W and X nodes and removing
 * anything dirichlet related
 */
template <class V>
int VLet3d<V>::bldPerUnodes(){

  vector<int> orderBoxesVec;
  vector<V>& nodeVec = this->nodeVec();
  vector<PerNode>& perNodeVec = this->perNodeVec();
  
  iC( this->dwnOrderCollect(orderBoxesVec	) );
  for (int i = 1; i < orderBoxesVec.size(); i++){
	 int gni = orderBoxesVec[i];
	 perNodeVec[gni].bdryUnodes().clear(); /* Make sure is clear */
  }
  for (int i = 1; i < orderBoxesVec.size(); i++){
	 int gni = orderBoxesVec[i];
	 double gniRad = this->radius(gni); /* This node's radius */
	 /* If this node's parent touches the boundary */
	 if (nodeVec[gni].terminal() && perNodeVec[gni].bdry()){
		//double parRad = this->radius(par); /* This node's parent's radius */
		for (int j = 0; j < perNodeVec[gni].bdryNbrs().size(); j++){
		  Point3 nbrOff = perNodeVec[gni].bdryNbrs()[j].bdryOffset();
		  /* parent's j'th periodic nbr's global index */
		  int nbr = perNodeVec[gni].bdryNbrs()[j].gni();
		  int type = FLN;
		  Point3 nbrCtr = this->center(nbr);
		  Point3 gniCtr = this->center(gni);
		  for (int a = 0; a < 3; a++) { nbrCtr(a) += (nbrOff(a)); }
		  Point3 dif = nbrCtr - gniCtr;
		  if (dif.linfty() == 2.0*gniRad){ /* nbr is adjacent */
			 if (nodeVec[nbr].terminal()){ /* Nbr is a leaf at same level */
				iA(this->depth(gni) == this->depth(nbr));
				perNodeVec[gni].bdryUnodes().push_back(BdryNode(nbr, nbrOff, type));
			 }
			 else {
				/* Descendants of nbr Possibly a fine neighbor or W node */
				/* Need to buidl descendants and then check them */
				if (perNodeVec[nbr].dscList().size() == 0){
				  bldDscList(nbr);
				}

				for (int i = 0; i < perNodeVec[nbr].dscList().size(); i++){
				  int nbrDsc = perNodeVec[nbr].dscList()[i];
				  Point3 nbrDscCtr = this->center(nbrDsc);
				  Point3 difDsc = (nbrDscCtr + nbrOff) - this->center(gni);
				  int levDif = this->depth(nbrDsc) - this->depth(gni);
				  //cout << levDif << " " << dif << endl;
				  iA(levDif >= 1);
				  if (this->terminal(nbrDsc)){
					 if (difDsc.linfty() == gniRad*(1.0 + pow(0.5, levDif))){
						/* Fine Nbr */
						perNodeVec[gni].bdryUnodes().push_back(BdryNode(nbrDsc, nbrOff, type));
						/* Coarse Nbr */
						perNodeVec[nbrDsc].bdryUnodes().push_back(BdryNode(gni, -nbrOff, type));
					 }
				  }
				}
			 }
		  }
		}
	 }
  }
  //exit(0);
  return(0);
}

/* This is all very similar to the periodic code, so I may combine it, but
 * for now this is for Dirichlet code
 */

/* Checking for same level neighbors.  When we build the lists, we only look
 * at a node's a parent's neighbors and their children.  By using inverse
 * relationships, the X and Crse lists will be built at the same time as the W
 * and Fine lists
 */
template <class V>
int VLet3d<V>::bldDirNbrs(){
  
  vector<int> orderBoxesVec;
  vector<V>& nodeVec = this->nodeVec();
  vector<PerNode>& perNodeVec = this->perNodeVec();
  
  iC( this->dwnOrderCollect(orderBoxesVec	) );
  /* Root is at zero */
  int gni =  orderBoxesVec[0];
  Point3 ctr(this->center());
  int type;
  /* For the root node, create the boundary neighbors since
	* we know where they are */
  for (int a = -1; a <= 1; a++){
	 for (int b = -1; b <= 1; b++){
		for (int c = -1; c <= 1; c++){
		  if (!(a == 0 && b == 0 && c == 0)){
			 double drad = 2.0*this->radius();
			 Point3 off(drad*a, drad*b, drad*c);
			 if (c == -1 || c == 1){
				if (b == -1 || b == 1){
				  if (a == -1 || a == 1) type = FLXYZ;
				  else type = FLYZ;
				}
				else { /* b == 0 */
				  iA (b == 0);
				  if (a == -1 || a == 1) type = FLXZ;
				  else type = FLZ;
				}
			 }
			 else {
				iA(c == 0);
				if (b == -1 || b == 1){
				  if (a == -1 || a == 1) type = FLXY;
				  else type = FLY;
				}
				else { /* b == 0 */
				  iA (b == 0);
				  if (a == -1 || a == 1) type = FLX;
				  else type = FLN; /* Should never happen */
				}
			 }
			 iA (type != FLN);
			 perNodeVec[gni].bdryNbrs().push_back(BdryNode(gni, off, type));
			 //cout << gni << " " << off << " " << type << endl;
		  }
		}
	 }
  }

  /* For nodes starting at the root's children and down, do the following:
	* - Check to see if the node is on the boundary (already set)
	* - Look at a node's parent's periodic dirichlet neighbors
	* - If these parental neighbors have children,
	* see if these children are adjacent once they are offset the
	* proper amount as well as "flipped" as necessary
	*.  If so, they are dirichlet neighbors of the
	* current box with offset and flip typethe same as their parents
	*/
  for (int i = 1; i < orderBoxesVec.size(); i++){
	 gni = orderBoxesVec[i];
	 if ((nodeVec[gni].depth() + this->rootLevel()> this->rootLevel()) && perNodeVec[gni].bdry()){
		int par = nodeVec[gni].parent(); /* This node's parent */
		double parRad = this->radius(par); /* This node's parent's radius */
		double gniRad = this->radius(gni); /* This node's radius */

		for (int j = 0; j < perNodeVec[par].bdryNbrs().size(); j++){
		  /* parNbrOff tells us what we would have to add to the
			* center of this node's parent's j'th periodic nbr
			* to put it where it belongs periodically
			*/
		  Point3 parNbrOff = perNodeVec[par].bdryNbrs()[j].bdryOffset();
		  int parNbr = perNodeVec[par].bdryNbrs()[j].gni(); /* parent's j'th periodic nbr's global index */
		  int type = perNodeVec[par].bdryNbrs()[j].type();
		  //cout << this->center(parNbr) << " " << parNbr << " " << parNbrOff << endl; //exit(0);
		  //cout << type << " " << this->center(gni) << endl;
		  //exit(0);
		  if (!nodeVec[parNbr].terminal()){
			 /* Go through parent's neighbor's children */
			 for (int a = 0; a < 2; a++){
				for (int b = 0; b < 2; b++){
				  for (int c = 0; c < 2; c++){
					 int parNbrChi = this->child(parNbr, Index3(a,b,c));
					 Point3 parNbrChiCtr = this->center(parNbrChi); /* par's jth per nbr's center */
					 //cout << parNbrChiCtr << " ";
					 flpCtrDirNbr(parNbrChiCtr, type);
					 //cout << parNbrChiCtr << " " << " " ;
					 for (int d = 0; d < 3; d++) { parNbrChiCtr(d) += (parNbrOff(d)); } /* Move by offset */
					 //cout << parNbrChiCtr << " " << endl;
					 Point3 dif = parNbrChiCtr - this->center(gni);
					 if (dif.linfty() == parRad){ /* kth child is adjacent */
						perNodeVec[gni].bdryNbrs().push_back(BdryNode(parNbrChi, parNbrOff, type));
						//cout << " DIF " << dif << endl;
						//if (this->terminal(parNbrChi)) cout << gni << " " << this->depth(parNbrChi) << endl; 
						//iA( gni != parNbrChi);
					 }
					 else {  }
				  }
				}
			 }
		  }
		  else {
			 //Point3 parNbrCtr = this->center(parNbr) + parNbrOff;
			 /* Why did I put this here? */
		  }
		}
	 }
	 //cout << perNodeVec[gni].bdryNbrs().size() << endl;
  }
  /* Build periodic V, U, W, and X lists */
  bldDirVnodes();
  bldPerUWXnodes();
  //exit(0);
  return(0);
 
}

template <class V>
int VLet3d<V>::bldDirVnodes(){
  
  vector<int> orderBoxesVec;
  vector<V>& nodeVec = this->nodeVec();
  vector<PerNode>& perNodeVec = this->perNodeVec();
  iC( this->dwnOrderCollect(orderBoxesVec	) );
  for (int i = 1; i < orderBoxesVec.size(); i++){
  //for (int i = 1; i < 0; i++){
	 int gni = orderBoxesVec[i];
	 int par = nodeVec[gni].parent(); /* This node's parent */
	 /* If this node's parent touches the boundary */
	 if ((nodeVec[gni].depth() + this->rootLevel() > this->rootLevel()) && perNodeVec[par].bdry()){
		double parRad = this->radius(par); /* This node's parent's radius */
		//double gniRad = this->radius(gni); /* This node's radius */
		for (int j = 0; j < perNodeVec[par].bdryNbrs().size(); j++){
		  Point3 parNbrOff = perNodeVec[par].bdryNbrs()[j].bdryOffset();
		  /* parent's j'th periodic nbr's global index */
		  int parNbr = perNodeVec[par].bdryNbrs()[j].gni();
		  int type = perNodeVec[par].bdryNbrs()[j].type();
		  if (!nodeVec[parNbr].terminal()){
			 /* Go through parent's periodic neighbor's children */
			 for (int a = 0; a < 2; a++){
				for (int b = 0; b < 2; b++){
				  for (int c = 0; c < 2; c++){
					 int parNbrChi = this->child(parNbr, Index3(a,b,c)); /* kth child of par's per' jth nbr */
					 Point3 parNbrChiCtr(this->center(parNbrChi)); /* par's jth per nbr's center */
					 flpCtrDirNbr(parNbrChiCtr, type); /* "Flip" as necessary */
					 for (int d = 0; d < 3; d++) { parNbrChiCtr(d) += (parNbrOff(d)); } /* Move by offset */
					 Point3 dif = parNbrChiCtr - this->center(gni);
					 if (dif.linfty() > parRad){ /* kth child is NOT adjacent */
						iA( dif.linfty() >= 2.0*parRad);
						if (nodeVec[parNbrChi].tag() & LET_SRCNODE){	
						  perNodeVec[gni].bdryVnodes().push_back(BdryNode(parNbrChi, parNbrOff, type));
						}
						//cout << "VNODE " << gni << " " << j << " " << type << " " << parNbrChi << " " <<  parNbrOff << parNbrChiCtr << endl;
					 }
				  }
				}				 
			 }
		  }
		}
	 }
  }
  //exit(0);
  return(0);
}

template <class V>
int VLet3d<V>::flpCtrDirNbr(Point3 &ctrDirNbr, const int flpTyp){
  
  Point3 ctr(this->center());
  /* Should do this with a series of "Reflection" matrices */
  iA(flpTyp != FLN);
  if (flpTyp == FLN) { /* Don't do anything! */
	 return(0);
  }

  if (flpTyp == FLX || flpTyp == FLXZ || flpTyp == FLXY || flpTyp == FLXYZ){
	 ctrDirNbr(0) = ctrDirNbr(0) + 2.0*(ctr(0) - ctrDirNbr(0));
  }
  if (flpTyp == FLY || flpTyp == FLYZ || flpTyp == FLXY || flpTyp == FLXYZ){
	 ctrDirNbr(1) = ctrDirNbr(1) + 2.0*(ctr(1) - ctrDirNbr(1));
  }
  if (flpTyp == FLZ || flpTyp == FLYZ || flpTyp == FLXZ || flpTyp == FLXYZ){
	 ctrDirNbr(2) = ctrDirNbr(2) + 2.0*(ctr(2) - ctrDirNbr(2));
  }
  //cerr << ctrDirNbr << " " << ctr << endl;
  return(0);
}

template <class V>
int VLet3d<V>::flpPDif(Index3 &pDif, const int flpTyp){
  
  Point3 ctr(this->center());
  /* Should do this with a series of "Reflection" matrices */
  iA(flpTyp != FLN);
  if (flpTyp == FLN) { /* Don't do anything! */
	 
  }
  if (flpTyp == FLX || flpTyp == FLXZ || flpTyp == FLXY || flpTyp == FLXYZ){
	 if (pDif(0) == 1) pDif(0) = 0;
	 else pDif(0) = 1;
  }
  if (flpTyp == FLY || flpTyp == FLYZ || flpTyp == FLXY || flpTyp == FLXYZ){
	 if (pDif(1) == 1) pDif(1) = 0;
	 else pDif(1) = 1;
  }
  if (flpTyp == FLZ || flpTyp == FLYZ || flpTyp == FLXZ || flpTyp == FLXYZ){
	 if (pDif(2) == 1) pDif(2) = 0;
	 else pDif(2) = 1;
  }
  return(0);
}

template <class V>
int VLet3d<V>::bldDscList(int gni){
  iA( !(this->terminal(gni)));
  vector<PerNode>& perNodeVec = this->perNodeVec();
  int beg = 0;
  int end = 8;

  for (int a = 0; a < 2; a++){
	 for (int b = 0; b < 2; b++){
		for (int c = 0; c < 2; c++){
		  perNodeVec[gni].dscList().push_back(this->child(gni,Index3(a,b,c)));
		}
	 }
  }

  int numnew = 0;
  
  while(beg < end){
	 numnew = 0;
	 for (int i = beg; i < end; i++){
		int dsc = perNodeVec[gni].dscList()[i];
		if (!(this->terminal(dsc))){
		  for (int a = 0; a < 2; a++){
			 for (int b = 0; b < 2; b++){
				for (int c = 0; c < 2; c++){
				  perNodeVec[gni].dscList().push_back(this->child(dsc,Index3(a,b,c)));
				  numnew++;
				}
			 }
		  }
		}
	 }
	 beg = end;
	 end = end + numnew;
  }
  return(0);
}

template <class V>
DblNumMat VLet3d<V>::grdOverExaPos(int gNodeIdx, bool dep) 
{	
  /* If dep is true, gNodeIdx represents a dep level, not node id */
  double R = (dep ? 1.0/pow(2.0,gNodeIdx + this->_rootLevel) : 1.0/pow(2.0,this->depth(gNodeIdx) + this->_rootLevel));
  DblNumMat grdPos(this->dim(),grdOverSamPos().n()); clear(grdPos);
  Point3 ctr(this->center(gNodeIdx));
  for (int l = 0; l < grdPos.n(); l++){
	 for (int d = 0; d < this->dim(); d++){
		grdPos(d,l) = ctr(d);
	 }
  }
  
  daxpy(R, grdOverSamPos(), grdPos);
  return grdPos;
}

template <class V>
DblNumMat VLet3d<V>::grdDblSrcExaPos(int gNodeIdx, bool dep)
{
  /* If dep is true, gNodeIdx represents a dep level, not node id */
  double R = (dep ? 1.0/pow(2.0,gNodeIdx + this->_rootLevel) : 1.0/pow(2.0,this->depth(gNodeIdx) + this->_rootLevel));
  DblNumMat grdPos(this->dim(),grdDblSrcSamPos().n()); clear(grdPos);
  Point3 ctr(this->center(gNodeIdx));
  for (int l = 0; l < grdPos.n(); l++){
	 for (int d = 0; d < this->dim(); d++){
		grdPos(d,l) = ctr(d);
	 }
  }
  
  daxpy(R, grdDblSrcSamPos(), grdPos);
  return grdPos;
}

// ---------------------------------------------------------------------
template <class V>
DblNumMat VLet3d<V>::grdSrcExaPos(int gNodeIdx, bool dep)
{
  /* If dep is true, gNodeIdx represents a dep level, not node id */
  double R = (dep ? 1.0/pow(2.0,gNodeIdx + this->_rootLevel) : 1.0/pow(2.0,this->depth(gNodeIdx) + this->_rootLevel));
  DblNumMat grdPos(this->dim(),grdSrcSamPos().n()); clear(grdPos);
  Point3 ctr(this->center(gNodeIdx));
  for (int l = 0; l < grdPos.n(); l++){
	 for (int d = 0; d < this->dim(); d++){
		grdPos(d,l) = ctr(d);
	 }
  }
  
  daxpy(R, grdSrcSamPos(), grdPos);
  return grdPos;
}

template <class V>
DblNumMat VLet3d<V>::grdTrgExaPos(int gNodeIdx, bool dep)
{
  /* If dep is true, gNodeIdx represents a dep level, not node id */
  double R = (dep ? 1.0/pow(2.0,gNodeIdx + this->_rootLevel) : 1.0/pow(2.0,this->depth(gNodeIdx) + this->_rootLevel));
  DblNumMat grdPos(this->dim(),trgGrdSze()); clear(grdPos);
  Point3 ctr(this->center(gNodeIdx));
  for (int l = 0; l < grdPos.n(); l++){
	 for (int d = 0; d < this->dim(); d++){
		grdPos(d,l) = ctr(d);
	 }
  }
  
  daxpy(R, grdTrgSamPos(), grdPos);
  return grdPos;
}
