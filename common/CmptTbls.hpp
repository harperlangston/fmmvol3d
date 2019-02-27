#include "Dcuhre.hpp"
#include "common/vecmatop.hpp"
#include "kernel3d.hpp"
#include "syms.hpp"

#ifndef _VFMM_CMPTTBLS_HPP_
#define _VFMM_CMPTTBLS_HPP_

class CmptTbls: public ComObject
{
private:
  int _ksrcval;
  int _ktrgval;
  int _np;
  int _kt;
  double _lambda;

  Dcuhre* _dcuhre;

public:
  CmptTbls(const string& p): ComObject(p) {;}
  ~CmptTbls();

  int setup(map<string,string>& opts);

  int &ksrcval() { return _ksrcval; }
  int &ktrgval() { return _ktrgval; }
  int &np() { return _np; }
  int &kt() { return _kt; }
  double &lambda() { return _lambda; }

  Dcuhre* dcuhre() { return _dcuhre; }

  int eval(int s2m, int nbrs, int wx, int level, DblNumMat& _S2M, map<int, DblNumMat>& _SC2TV_map, map<int, DblNumMat>& _WXC2TV_map);
  //int eval(int s2m, int nbrs, int wx) { return eval(s2m,nbrs,wx,0); }
  int eval_s2m(int level, DblNumMat& _S2M) { map<int,DblNumMat> mt; return eval(1,0,0,level,_S2M, mt, mt); }
  int eval_nbrs(int level, map<int, DblNumMat>& _SC2TV_map) { DblNumMat mt; return eval(0,1,0,level,mt,_SC2TV_map, _SC2TV_map); }
  int eval_wx(int level, map<int,DblNumMat>& _WXC2TV_map) { DblNumMat mt; return eval(0,0,1,level, mt,_WXC2TV_map, _WXC2TV_map); }

  int eval_unbal(const int type, const int id, const int levgni, const double radgni, const Point3 ctrgni, const int levnbr, const double radnbr, const Point3 ctrnbr, DblNumMat& _evalMat) { return eval_unbal(type, id, levgni, radgni, ctrgni, levnbr, radnbr, ctrnbr, Point3(0,0,0), _evalMat); }
  int eval_unbal(const int type, const int id, const int levgni, const double radgni, const Point3 ctrgni, const int levnbr, const double radnbr, const Point3 ctrnbr, const Point3 offset, DblNumMat& _evalMat);
  int lookup(const int type, const int levgni, const double radgni, const Point3 ctrgni, const int levnbr, const double radnbr, const Point3 ctrnbr) { return lookup(type, levgni, radgni, ctrgni, levnbr, radnbr, ctrnbr, Point3(0,0,0)); }
  int lookup(const int type, const int levgni, const double radgni, const Point3 ctrgni, const int levnbr, const double radnbr, const Point3 ctrnbr, const Point3 offset);

};

#endif
