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
#ifndef _LET3D_HPP_
#define _LET3D_HPP_

#include "common/vec3t.hpp"
#include "common/nummat.hpp"
#include "common/offtns.hpp"
#include "common/comobject.hpp"

using std::vector;

enum {
  LET_SRCNODE = 1,
  LET_TRGNODE = 2
};

class Node {
protected:
  int _parent, _child;
  Index3 _path2Node;
  int _depth;
  int _tag; //s empty, t empty ...

  map<int, DblNumVec> _dirEffDen;
  DblNumVec _effDen;
  
  /*! source node index */
  int _srcNodeIdx;
  /*! source exact beginning index */
  int _srcExaBeg;
  /*! source exact number */
  int _srcExaNum;
  /*! source own vector of indices */
  vector<int> _srcOwnVecIdxs;
  
  /*! target node index */
  int _trgNodeIdx;
  /*! target exact beginning index */
  int _trgExaBeg;
  /*! target exact number */
  int _trgExaNum;
  /*! target own vector of indices */
  vector<int> _trgOwnVecIdxs;
  
  /*! Nodes in the U-List */
  vector<int> _Unodes;
  /*! Nodes in the V-List */
  vector<int> _Vnodes;
  /*! Nodes in the W-List */
  vector<int> _Wnodes;
  /*! Nodes in the X-List */
  vector<int> _Xnodes;
  
public:
  Node(int p, int c, Index3 t, int d):
	 _parent(p), _child(c), _path2Node(t), _depth(d), _tag(false),
	 _srcNodeIdx(0), _srcExaBeg(0), _srcExaNum(0),
	 _trgNodeIdx(0), _trgExaBeg(0), _trgExaNum(0) {;}

  Node():
	 _parent(-1), _child(-1), _path2Node(Index3(-1,-1,-1)), _depth(-1), _tag(false),
	 _srcNodeIdx(0), _srcExaBeg(0), _srcExaNum(0),
	 _trgNodeIdx(0), _trgExaBeg(0), _trgExaNum(0) {;}

  /*! Return effective values */
  //DblNumVec& effVal() { return _effVal; }
  
  /*! Return effective densities */
  DblNumVec& effDen() { return _effDen; } /* Normally return the first one */
  DblNumVec& effDen(int i) { if (i == 0) return _effDen; else return _dirEffDen[i]; } /* For mixed bdry conditions */

  int& parent()        { return _parent; }
  int& child()         { return _child; }
  Index3& path2Node()  { return _path2Node; }
  int& depth()         { return _depth; }
  int& tag()           { return _tag; }
  bool terminal()      { return child()==-1; } //no child
  
  int& srcNodeIdx()          { return _srcNodeIdx; }
  int& srcExaBeg()           { return _srcExaBeg; }
  int& srcExaNum()           { return _srcExaNum; }
  vector<int>& srcOwnVecIdxs() { return _srcOwnVecIdxs; }
  
  int& trgNodeIdx()          { return _trgNodeIdx; }
  int& trgExaBeg()           { return _trgExaBeg; }
  int& trgExaNum()           { return _trgExaNum; }
  vector<int>& trgOwnVecIdxs() { return _trgOwnVecIdxs; }
  
  vector<int>& Unodes() { return _Unodes; }
  vector<int>& Vnodes() { return _Vnodes; }
  vector<int>& Wnodes() { return _Wnodes; }
  vector<int>& Xnodes() { return _Xnodes; }
};

//---------------------------------------------------------------------------
// unique identifier: a set of points and its bounding box
template <class N>
class Let3d: public ComObject
{
public:
  //---------------------------------------
  
  //----------------------------------------------
protected:
  //PARAMS(REQ)
  DblNumMat* _srcPos;
  DblNumMat* _trgPos;
  DblNumVec* _srcDen;
  DblNumVec* _trgVal;

  vector<int> _trgPosTrmGni; /* Keep track of nodes that target positions are in */
  
  Point3 _center;
  int _rootLevel;
  
  //PARAMS(OPT)
  int _ptsMax;
  int _maxLevel; //AT MOST 16

  //COMPONENTS
  vector<N> _nodeVec;
  
  /*! total number of levels */
  int _level; 
  int _srcNodeCnt, _srcExaCnt;
  int _trgNodeCnt, _trgExaCnt;
private:

public:
  Let3d(const string& p);
  ~Let3d();
  //MEMBER ACCESS
  DblNumMat*& srcPos() { return _srcPos; }
  DblNumMat*& trgPos() { return _trgPos; }
  DblNumVec*& srcDen() { return _srcDen; }
  DblNumVec*& trgVal() { return _trgVal; }

  vector<int>& trgPosTrmGni() { return _trgPosTrmGni; }
  
  Point3& center() { return _center; }
  int& rootLevel() { return _rootLevel; }
  double radius()  { return pow(2.0, -_rootLevel); }
  int& ptsMax()    { return _ptsMax; }
  int& maxLevel()  { return _maxLevel; }

  
  //SETUP AND USE
  virtual int setFromOptions(map<string,string>& opts);
  virtual int setup();
  
  int print();
  //access
  int level() { return _level; }
  vector<N>&  nodeVec() { return  _nodeVec; }
  int& srcNodeCnt()   { return _srcNodeCnt; }
  int srcExaCnt()  { return _srcExaCnt; }
  int trgNodeCnt()   { return _trgNodeCnt; }
  int trgExaCnt()  { return _trgExaCnt; }
  //construction
  int srcData();
  int trgData();
  
  //---------------------------------------------
  //LOCAL
  /*! build U,V,W,X lists */
  int calgnext(int gNodeIdx);
  
  /*! top down ordering of the nodes */
  int dwnOrderCollect(vector<int>&);
  /*! bottom up ordering of the nodes */
  int upwOrderCollect(vector<int>&);

  /*! in reverse order by levels */
  int revLvlOrderCollect(map<int, vector<int> > &);

  int postOrderCollect(vector<int>& orderBoxesVec	);
  void postOrderCollect(vector<int>& orderBoxesVec, int gni);
  /*! node access */
  N& node(int gNodeIdx) { return _nodeVec[gNodeIdx]; }

  //tree traversal and informational retrieval
  bool    root(int gNodeIdx)     { return node(gNodeIdx).parent()==-1; } //no parent
  bool    terminal(int gNodeIdx) { return node(gNodeIdx).child()==-1; } //no child
  int     parent(int gNodeIdx)   { assert(node(gNodeIdx).parent()!=-1); return node(gNodeIdx).parent(); }
  int     child( int gNodeIdx, const Index3& idx);
  Index3  path2Node(int gNodeIdx) { return node(gNodeIdx).path2Node(); }
  int     depth(int gNodeIdx)      { return node(gNodeIdx).depth(); }
  int&    tag(int gNodeIdx)        { return node(gNodeIdx).tag(); }
  
  Point3  center(int gNodeIdx);  
  double  radius(int gNodeIdx);  
  
  int findgnt(int depth, const Index3& path2Node);  
  bool adjacent(int, int);                
  
  int dim() const { return 3; }

};
  
#endif
