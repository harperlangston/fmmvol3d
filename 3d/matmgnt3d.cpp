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
	Software Foundation, Inc., 9 Temple Place - Suite 330, Boston, MA
	02111-1307, USA.  */
#include "matmgnt3d.hpp"
#include "common/vecmatop.hpp"
#include <iostream>

using std::cerr;
using std::endl;

// ---------------------------------------------------------------------- 
//double MatMgnt3d::_wsbuf[65536];

// ---------------------------------------------------------------------- 
MatMgnt3d::MatMgnt3d(): _np(6){
#ifdef FFTW3
  _forplan = NULL;
  _invplan = NULL;
  _ue2dcplan = NULL;
#else
  _forplan = NULL;
  _invplan = NULL;
#endif
}

MatMgnt3d::MatMgnt3d(Kernel3d knl, int np): _knl(knl), _np(np){
#ifdef FFTW3
  _forplan = NULL;
  _invplan = NULL;
  _ue2dcplan = NULL;
#else
  _forplan = NULL;
  _invplan = NULL;
#endif
  setup();
}

// ---------------------------------------------------------------------- 
MatMgnt3d::~MatMgnt3d()
{
#ifdef FFTW3
  if(_forplan!=NULL) {	 fftw_destroy_plan(_forplan); _forplan = NULL;    }
  if(_invplan!=NULL) { 	 fftw_destroy_plan(_invplan); _invplan=NULL;  }
  if(_ue2dcplan!=NULL) { fftw_destroy_plan(_ue2dcplan); _ue2dcplan=NULL;  }
#else
  if(_forplan!=NULL) {	 rfftwnd_destroy_plan(_forplan); _forplan=NULL;  }
  if(_invplan!=NULL) { 	 rfftwnd_destroy_plan(_invplan); _invplan=NULL;  }
#endif
}

double MatMgnt3d::alt()
{
  double p = (double)(_np+1);
  if (p > 9) p = 9.00;
  return pow(0.1, p);
}
// ---------------------------------------------------------------------- 
// ----------------------------------------------------------------------

int MatMgnt3d::cleanPlans(){
#ifdef FFTW3
  if(_forplan!=NULL) {	 fftw_destroy_plan(_forplan); _forplan = NULL;    }
  if(_invplan!=NULL) { 	 fftw_destroy_plan(_invplan); _invplan=NULL;  }
  if(_ue2dcplan!=NULL) { fftw_destroy_plan(_ue2dcplan); _ue2dcplan=NULL;  }
#else 
  if(_forplan!=NULL) {	 rfftwnd_destroy_plan(_forplan); _forplan=NULL;  }
  if(_invplan!=NULL) { 	 rfftwnd_destroy_plan(_invplan); _invplan=NULL;  }
#endif
  return (0);
}

int MatMgnt3d::setup()
{
  //--------------------------------------------------------
  _hom = _knl.homogeneous();
  _compPerDir = false;
  if(_hom==true) {
    _knl.homogeneousDeg(_degVec); iA(_degVec.size()==srcDOF());
  }
  iC( samPosCal(_np,   1.0, _samPos[UE], UE) );
  iC( samPosCal(_np+2,   3.0, _samPos[UC], UC) );
  iC( samPosCal(_np,   3.0, _samPos[DE], DE) );
  iC( samPosCal(_np,   1.0, _samPos[DC], DC) );
  iC( regPosCal(_np,   1.0, _regPos    ) ); //only one regPos

#ifdef FFTW3
  int effNum = (2*_np+2)*(2*_np)*(2*_np);
  DblNumVec effVal(effDatSze(DC)); fftw_complex* valPtr  = (fftw_complex*)(effVal.data());
  DblNumVec regVal(regPos().n()*trgDOF());
  DblNumVec regDen(regPos().n()*srcDOF()); clear(regDen);
  DblNumVec effDen(effDatSze(UE)); fftw_complex* denPtr  = (fftw_complex*)(effDen.data());
  DblNumMat tmp(trgDOF(),regPos().n()*srcDOF());
  DblNumMat UE2DC(trgDOF()*srcDOF(), effNum); //_memused[4] += dof*dof*effNum(UE)*sizeof(double);
  int np = _np;
  int nnn[3];  nnn[0] = 2*np;  nnn[1] = 2*np;  nnn[2] = 2*np;
  _invplan = fftw_plan_many_dft_c2r(3,nnn,trgDOF(), valPtr, NULL, trgDOF(),1, regVal.data(),NULL, trgDOF(),1, FFTW_ESTIMATE);
  _forplan = fftw_plan_many_dft_r2c(3,nnn,srcDOF(), regDen.data(),NULL, srcDOF(),1, denPtr, NULL, srcDOF(),1, FFTW_ESTIMATE);
  _ue2dcplan = fftw_plan_many_dft_r2c(3,nnn,srcDOF()*trgDOF(), tmp.data(),NULL, srcDOF()*trgDOF(),1, (fftw_complex*)(UE2DC.data()), NULL, srcDOF()*trgDOF(),1, FFTW_ESTIMATE);
#else 
  _forplan = rfftw3d_create_plan(2*_np,2*_np,2*_np, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_OUT_OF_PLACE);
  _invplan = rfftw3d_create_plan(2*_np,2*_np,2*_np, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_OUT_OF_PLACE);
#endif  
  //--------------------------------------------------------
  return (0);
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::report()
{
  cerr << "matrix map size"<<endl;
  cerr << _upwEqu2UpwChk.size() << " ";
  cerr << _upwChk2UpwEqu.size() << " ";
  cerr << _dwnChk2DwnEqu.size() << " ";
  cerr << _dwnEqu2DwnChk.size() << " ";
  cerr << _upwEqu2DwnChk.size() << endl;
  return (0);
}

// ---------------------------------------------------------------------- 
int MatMgnt3d::plnDatSze(int tp)
{
  return npPlnDatSze(tp,_np,srcDOF(),trgDOF());
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::effDatSze(int tp)
{
  return npEffDatSze(tp,_np,srcDOF(),trgDOF());  
}

// ---------------------------------------------------------------------- 
int MatMgnt3d::UpwChk2UpwEqu_dgemv(const int l, const DblNumVec& chk, DblNumVec& den, const double scale)
{
  iA(scale == 2.0 || scale == 3.0);
  DblNumMat& _UC2UE = ((scale == 2.0) ? ((_hom==true) ? _upwChk2UpwEqu[0] : _upwChk2UpwEqu[l]) : ((_hom==true) ? perdirmaps()->PerUpwChk2UpwEqu()[0] : perdirmaps()->PerUpwChk2UpwEqu()[l]));
  double R = (_hom==true) ? 1 : 1.0/pow(scale, l);
  //---------compute matrix
  if(_UC2UE.m()==0) {	 //cerr<<"UpwChk2UpwEqu compute"<<endl;
	 //set matrix
	 {
		DblNumMat ud2c(plnDatSze(UC), plnDatSze(UE));
		DblNumMat chkPos(dim(),samPos(UC).n());	 clear(chkPos);	 iC( daxpy(R, samPos(UC), chkPos) ); //scale
		DblNumMat denPos(dim(),samPos(UE).n());	 clear(denPos);	 iC( daxpy(R, samPos(UE), denPos) ); //scale
	 
		iC( _knl.kernel(denPos, denPos, chkPos, ud2c) );	
		_UC2UE.resize(plnDatSze(UE), plnDatSze(UC));
	 
		iC( pinv(ud2c, alt(), _UC2UE) );
	 }
  }
  //---------matvec
  if(_hom==true) {
    //matvec
    int srcDOF = this->srcDOF();
    DblNumVec tmpDen(srcDOF*samPos(UE).n());
	 //cerr << chk << endl;  
	 //cerr << tmpDen << endl;
	 //cerr << _UC2UE << endl;
	 iC( dgemv(_UC2UE, chk, tmpDen) );
	 //scale
    vector<double> sclvec(srcDOF);	 
    for(int s=0; s<srcDOF; s++)	 sclvec[s] = pow(scale, - l*_degVec[s]);
    int cnt = 0;
    for(int i=0; i < samPos(UE).n(); i++)
      for(int s=0; s < srcDOF; s++) {
		  den(cnt) = den(cnt) + tmpDen(cnt) * sclvec[s];
		  cnt++;
      }
  } else {
    //iC( dgemv(1.0, _UC2UE, chk, 1.0, den) );
	 iC( dgemv(_UC2UE, chk, den) );
  }
  return (0);
}
/* For periodic/dirichlet: idx is of form 0,-1,1 for offset from center (so 0,0,0 is center node)
 * For free-space: idx of form -1,1, offset b -0.5, +0.5 from center */
// ---------------------------------------------------------------------- 
int MatMgnt3d::UpwEqu2UpwChk_dgemv(const int lev, Index3 idx, const DblNumVec& den, DblNumVec& chk, const double scale)
{
  iA(scale == 2.0 || scale == 3.0);
  NumTns<DblNumMat>& _UE2UC = ((scale == 2.0) ? ((_hom==true) ? _upwEqu2UpwChk[0] : _upwEqu2UpwChk[lev]) : ((_hom==true) ? perdirmaps()->PerUpwEqu2UpwChk()[0] : perdirmaps()->PerUpwEqu2UpwChk()[lev]));
  double R = (_hom==true) ? 1 : 1.0/pow(scale, lev);

  if(_UE2UC.m()==0){
	 {
		if (scale == 2.0){ _UE2UC.resize(2,2,2); }
		if (scale == 3.0){ _UE2UC.resize(3,3,3); }
	 }
  }
  DblNumMat& _UE2UCii = ((scale == 2.0) ? _UE2UC(idx(0), idx(1), idx(2)) : _UE2UC(idx(0)+1, idx(1)+1, idx(2)+1));
  //---------compute matrix
  if(_UE2UCii.m()==0) {	 //cerr<<"UpwEqu2UpwChk compute"<<endl;
	 {
		double shft = (scale == 2.0) ? 1.0 : 0.0;
		_UE2UCii.resize(plnDatSze(UC), plnDatSze(UE)); //_memused[1] += plnnum(UC)*dof()*plnnum(UE)*dof()*sizeof(double);
		DblNumMat chkPos(dim(),samPos(UC).n());	 clear(chkPos);
		iC( daxpy(scale*R, samPos(UC), chkPos) ); //scale
		DblNumMat denPos(dim(),samPos(UE).n());	 clear(denPos);
		iC( daxpy(R, samPos(UE), denPos) ); //scale
		for(int i=0; i<dim(); i++){
		  for(int j=0; j<samPos(UE).n(); j++){
			 denPos(i,j) = denPos(i,j) + (2.0*idx(i)-shft)*R;//shift
		  }
		}
		iC( _knl.kernel(denPos, denPos, chkPos, _UE2UCii) );
	 }
  }
  //---------matvec
  if(_hom==true) {
    int srcDOF = this->srcDOF();
    DblNumVec tmpDen(srcDOF*samPos(UE).n());	 clear(tmpDen);
    vector<double> sclvec(srcDOF);
    for(int s=0; s<srcDOF; s++){
      sclvec[s] = pow(scale, (lev)*_degVec[s]);
    }
    int cnt = 0;
    for(int i=0; i<samPos(UE).n(); i++) {
      for(int s=0; s<srcDOF; s++) {
		  tmpDen(cnt) = den(cnt) * sclvec[s];
		  cnt++;
      }
	 }
	 iC( dgemv(_UE2UCii, tmpDen, chk) );
  }
  else {
	 iC( dgemv(_UE2UCii, den, chk) );
  }
  return (0);
}

/* If  need to store the check values, we do away with them */
int MatMgnt3d::UpwEqu2UpwEqu_dgemv(int l, Index3 idx, const DblNumVec& chdden, DblNumVec& parden, const double scale){
  DblNumVec chkval(plnDatSze(UC));
  iC( UpwEqu2UpwChk_dgemv(l, idx, chdden, chkval, scale));
  iC( UpwChk2UpwEqu_dgemv(l-1, chkval, parden, scale));
  return(0);
}

// ---------------------------------------------------------------------- 
int MatMgnt3d::DwnChk2DwnEqu_dgemv(const int level, const DblNumVec& chk, DblNumVec& den, const double scale)
{
  iA(scale == 2.0 || scale == 3.0);
  DblNumMat& _DC2DE = ((scale == 2.0) ? ((_hom==true) ? _dwnChk2DwnEqu[0]: _dwnChk2DwnEqu[level]) : ((_hom==true) ? perdirmaps()->PerDwnChk2DwnEqu()[0]: perdirmaps()->PerDwnChk2DwnEqu()[level]));
  double R = (_hom==true) ? 1 : 1.0/pow(scale,level);
  //---------compute matrix
  if(_DC2DE.m()==0) {	 //cerr<<"DwnChk2DwnEqu compute"<<endl;
    DblNumMat dd2c(plnDatSze(DC), plnDatSze(DE));
    DblNumMat chkPos(dim(),samPos(DC).n());		clear(chkPos);	 iC( daxpy(R, samPos(DC), chkPos) ); //scale
    DblNumMat denPos(dim(),samPos(DE).n());		clear(denPos);	 iC( daxpy(R, samPos(DE), denPos) ); //scale
    
    iC( _knl.kernel(denPos, denPos, chkPos, dd2c) );//matrix
    _DC2DE.resize(plnDatSze(DE), plnDatSze(DC)); //_memused[2] += plndnenum()*dof()*plndncnum()*dof()*sizeof(double);
    //if (_hom == true && _np < 10) { iC( pinv(dd2c, 1e-10, _DC2DE) ); }
    //else { iC( pinv(dd2c, alt(), _DC2DE) ); }
	 iC( pinv(dd2c, alt(), _DC2DE) );
	 //iC( tReg(dd2c, alt(), _DC2DE) );
  }
  //---------matvec
  if(_hom==true) {
    int srcDOF = this->srcDOF();
    DblNumVec tmpDen(srcDOF*samPos(DE).n());
    //iC( dgemv(1.0, _DC2DE, chk, 1.0, tmpDen) );
	 iC( dgemv(_DC2DE, chk, tmpDen) );
    //scale
    vector<double> sclvec(srcDOF);
    for(int s=0; s<srcDOF; s++)	{
      sclvec[s] = pow(scale, - level*_degVec[s]);
    }
    int cnt = 0;	 for(int i=0; i<samPos(DE).n(); i++)
		for(int s=0; s<srcDOF; s++) {
		  den(cnt) = den(cnt) + tmpDen(cnt) * sclvec[s];
		  cnt++;
		}
  } else {
    //iC( dgemv(1.0, _DC2DE, chk, 1.0, den) );
	 iC( dgemv(_DC2DE, chk, den) );
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::DwnEqu2DwnChk_dgemv(const int level, Index3 idx, const DblNumVec& den, DblNumVec& chk, const double scale)
{
  iA( scale == 2.0 || scale == 3.0);
  NumTns<DblNumMat>& _DE2DC = ((scale == 2.0) ? ((_hom==true) ? _dwnEqu2DwnChk[0] : _dwnEqu2DwnChk[level]) : ((_hom==true) ? perdirmaps()->PerDwnEqu2DwnChk()[0] : perdirmaps()->PerDwnEqu2DwnChk()[level]));
  double R = (_hom==true) ? 1.0 : 1.0/pow(scale, level);
  if(_DE2DC.m()==0){
	 if (scale == 2.0){ _DE2DC.resize(2,2,2); }
	 else if (scale == 3.0){ _DE2DC.resize(1,1,1); }
	 else { iA(0); }
  }

  DblNumMat& _DE2DCii = (scale == 2.0) ? _DE2DC(idx[0], idx[1], idx[2]) : _DE2DC(0, 0, 0);
  //---------compute matrix
  if(_DE2DCii.m()==0) {	 //cerr<<"DwnEqu2DwnChk compute"<<endl;
	 /* shift/scale factors for regular vs. periodic/dirichlet */
	 double shft = (scale == 2.0) ? 0.5 : 0.0;
	 _DE2DCii.resize(plnDatSze(DC), plnDatSze(DE)); //_memused[3] += plndncnum()*dof()*plndnenum()*dof()*sizeof(double);
	 DblNumMat denPos(dim(),samPos(DE).n());		  clear(denPos);	 iC( daxpy(R, samPos(DE), denPos) ); //scale
	 DblNumMat chkPos(dim(),samPos(DC).n());		  clear(chkPos);	 iC( daxpy(R/scale, samPos(DC), chkPos) ); //scale
	 for(int i=0; i<dim(); i++){
		for(int j=0; j<samPos(DC).n(); j++){
		  chkPos(i,j) = chkPos(i,j) + (double(idx(i))-shft)*R;
		}
	 }
	 
	 iC( _knl.kernel(denPos, denPos, chkPos, _DE2DCii) );
  }
  //---------matvec
  if(_hom==true) {
	 int srcDOF = this->srcDOF();
	 DblNumVec tmpDen(srcDOF*samPos(DE).n());
	 vector<double> sclvec(srcDOF);
	 for(int s=0; s<srcDOF; s++)		sclvec[s] = pow(scale, level*_degVec[s]);
	 int cnt = 0;
	 for(int i=0; i<samPos(DE).n(); i++) {
		for(int s=0; s<srcDOF; s++) {
		  tmpDen(cnt) = den(cnt) * sclvec[s];
		  cnt++;
		}
	 }
	 iC( dgemv(_DE2DCii, tmpDen, chk) );
  } else {
	 iC( dgemv(_DE2DCii, den, chk) );
  }	
  return (0);
}

int MatMgnt3d::DwnEqu2DwnEqu_dgemv(const int l, Index3 idx, const DblNumVec& parden, DblNumVec& chdden, const double scale){
  NumTns<DblNumMat>& _DE2DE = (_hom==true) ? _dwnEqu2DwnEqu[0] : _dwnEqu2DwnEqu[l];
  double R = (_hom==true) ? 1 : 1.0/pow(scale, l);  
  if(_DE2DE.m()==0)	 _DE2DE.resize(2,2,2);
  DblNumMat& _DE2DEii = _DE2DE(idx(0), idx(1), idx(2));
  //---------compute matrix
  if(_DE2DEii.m()==0) {	 //cerr<<"UpwEqu2UpwChk compute"<<endl;
    _DE2DEii.resize(plnDatSze(DE), plnDatSze(DE)); //_memused[1] += plnnum(UC)*dof()*plnnum(UE)*dof()*sizeof(double);
	 
	 DblNumMat de2dc(plnDatSze(DC), plnDatSze(DE));
    DblNumMat denPos(dim(),samPos(DE).n());	 clear(denPos);  iC( daxpy(R, samPos(DE), denPos) ); //scale
	 DblNumMat chkPos(dim(),samPos(DC).n());	 clear(chkPos);  iC( daxpy(0.5*R, samPos(DC), chkPos) ); //scale
	 for(int i=0; i<dim(); i++){
		for(int j=0; j<samPos(DC).n(); j++){
		  chkPos(i,j) = chkPos(i,j) + (double(idx(i))-0.5)*R;
		}
	 }
	 iC( _knl.kernel(denPos, denPos, chkPos, de2dc) );

	 DblNumMat dd2c(plnDatSze(DC), plnDatSze(DE));
	 clear(chkPos);	 iC( daxpy(0.5*R, samPos(DC), chkPos) ); //scale
    clear(denPos);	 iC( daxpy(0.5*R, samPos(DE), denPos) ); //scale
	 iC( _knl.kernel(denPos, denPos, chkPos, dd2c) );
	 
	 if (_hom == true){
		int srcDOF = this->srcDOF();
		vector<double> sclvec(srcDOF);
		for(int s=0; s<srcDOF; s++){
		  sclvec[s] = pow(scale, l*_degVec[s]);
		}
		iA(plnDatSze(DC) == de2dc.m() && samPos(DE).n() * srcDOF == de2dc.n());
		for (int i=0; i < plnDatSze(DC); i++){
		  for (int j=0; j < samPos(DE).n(); j++){
			 for (int s=0; s < srcDOF; s++){
				de2dc(i,j*srcDOF + s) = de2dc(i,j*srcDOF + s)*sclvec[s];
				dd2c(i,j*srcDOF + s) = dd2c(i,j*srcDOF + s)*sclvec[s];
			 }
		  }
		}
	 }
	 

	 DblNumMat dc2de(plnDatSze(DE), plnDatSze(DC));
	 iC( pinv(dd2c, alt(), dc2de) );
	 dgemm(1.0, dc2de, de2dc, 1.0, _DE2DEii);
	 //std::cout << _UE2UEii << endl;
  }
  iC( dgemv(1.0, _DE2DEii, parden, 1.0, chdden) );
  return(0);
}

// ---------------------------------------------------------------------- 
int MatMgnt3d::plnDen2EffDen(int level, const DblNumVec& plnDen, DblNumVec& effDen, const double scale)
{
  DblNumVec regDen(regPos().n()*srcDOF()); clear(regDen);
  //iC( samDen2RegularDen(plnDen, regDen) );
  //rfftwnd_real_to_complex(_forplan, srcDOF(), regDen.datmata(), srcDOF(), 1, (fftw_complex*)(effDen.data()), srcDOF(), 1);
  if(_hom==true) {
	 int srcDOF = this->srcDOF();
	 DblNumVec tmpDen(srcDOF*samPos(UE).n());
	 vector<double> sclvec(srcDOF);
	 for(int s=0; s<srcDOF; s++)		sclvec[s] = pow(scale, level*_degVec[s]);
	 int cnt = 0;
	 for(int i=0; i<samPos(UE).n(); i++)
		for(int s=0; s<srcDOF; s++) {
		  tmpDen(cnt) = plnDen(cnt) * sclvec[s];
		  cnt++;
		}
	 iC( samDen2RegDen(tmpDen, regDen) );
  } else {
	 iC( samDen2RegDen(plnDen, regDen) );
  }
#ifdef FFTW3
  //int np = _np;
  //int nnn[3];  nnn[0] = 2*np;  nnn[1] = 2*np;  nnn[2] = 2*np;
  fftw_complex* denPtr  = (fftw_complex*)(effDen.data());
  //fftw_plan forplan = fftw_plan_many_dft_r2c(3,nnn,srcDOF(), regDen.data(),NULL, srcDOF(),1, denPtr, NULL, srcDOF(),1, FFTW_ESTIMATE);
  fftw_execute_dft_r2c(_forplan, regDen.data(), denPtr);
  //fftw_destroy_plan(forplan);
#else
  rfftwnd_real_to_complex(_forplan, srcDOF(), regDen.data(), srcDOF(), 1, (fftw_complex*)(effDen.data()), srcDOF(), 1);
#endif  
  return (0);
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::samDen2RegDen(const DblNumVec& samDen, DblNumVec& regDen)
{
  int np = _np;
  int rgnum = 2*np;
  int srcDOF = this->srcDOF();
  int cnt=0;
  //the order of iterating is the same as SampleGrid
  for(int i=0; i<np; i++)
	 for(int j=0; j<np; j++)
		for(int k=0; k<np; k++) {
		  if(i==0 || i==np-1 || j==0 || j==np-1 || k==0 || k==np-1) {
			 //the position is fortran style
			 int rgoff = (k+np/2)*rgnum*rgnum + (j+np/2)*rgnum + (i+np/2);
			 for(int f=0; f<srcDOF; f++) {
				regDen(srcDOF*rgoff + f) += samDen(srcDOF*cnt + f);
			 }
			 cnt++;
		  }
		}  //iC( PetscLogFlops(np*np*np*dof) );
  return 0;
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::effVal2PlnVal(DblNumVec& effVal, DblNumVec& plnVal)
{
  DblNumVec regVal(regPos().n()*trgDOF());
  /* scale effective values */
  dscal(1.0/((double)(regPos().n())), effVal);
  int np = _np;
  int nnn[3];  nnn[0] = 2*np;  nnn[1] = 2*np;  nnn[2] = 2*np;
#ifdef FFTW3
  fftw_complex* valPtr  = (fftw_complex*)(effVal.data());
  //fftw_plan invplan = fftw_plan_many_dft_c2r(3,nnn,trgDOF(), valPtr, NULL, trgDOF(),1, regVal.data(),NULL, trgDOF(),1, FFTW_ESTIMATE);
  fftw_execute_dft_c2r(_invplan, valPtr, regVal.data());
  //fftw_destroy_plan(invplan);
  //invplan = NULL;
#else
  rfftwnd_complex_to_real(_invplan, trgDOF(), (fftw_complex*)(effVal.data()), trgDOF(), 1, regVal.data(), trgDOF(), 1);
#endif  
  iC( regVal2SamVal(regVal, plnVal) );
  return (0);
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::regVal2SamVal(const DblNumVec& regVal, DblNumVec& samVal)
{
  int np = _np;
  int rgnum = 2*np;
  int trgDOF = this->trgDOF();
  int cnt=0;
  //the order of iterating is the same as SampleGrid
  for(int i=0; i<np; i++)
	 for(int j=0; j<np; j++)
		for(int k=0; k<np; k++) {
		  if(i==0 || i==np-1 || j==0 || j==np-1 || k==0 || k==np-1) {
			 //the position is fortran style
			 int rgoff = (k+np/2)*rgnum*rgnum + (j+np/2)*rgnum + (i+np/2);
			 for(int f=0; f<trgDOF; f++) {
				samVal(trgDOF*cnt + f) += regVal(trgDOF*rgoff + f);
			 }
			 cnt++;
		  }
		}
  return 0;
}
// ---------------------------------------------------------------------- 
int MatMgnt3d::UpwEqu2DwnChk_dgemv(int lev, Index3 idx, const DblNumVec& effDen, DblNumVec& effVal, const double scale)
{
  OffTns<DblNumMat>& _UpwEqu2DwnChk = (_hom==true) ? _upwEqu2DwnChk[0] : _upwEqu2DwnChk[lev];
  double R = (_hom==true) ? 1.0 : 1.0/pow(scale, lev); //OffTns< DblNumMat >& _UpwEqu2DwnChk = _upwEqu2DwnChk[l];  double R       = 1.0/pow(2.0, l);
  /* Can make this 7 7 7 -3 -3 -3 if no periodicity - change later HARPER */
  if(_UpwEqu2DwnChk.m()==0){	 _UpwEqu2DwnChk.resize(9,9,9,-4,-4,-4); }
  //if(_UpwEqu2DwnChk.m()==0)	 _UpwEqu2DwnChk.resize(7,7,7,-3,-3,-3);
  DblNumMat& _UpwEqu2DwnChkii = _UpwEqu2DwnChk(idx[0], idx[1], idx[2]);
  //std::cout << idx << endl;	
  int np = _np;
  int effNum = (2*np+2)*(2*np)*(2*np);
  int srcDOF = this->srcDOF();
  int trgDOF = this->trgDOF();
  
  //---------compute matrix
  if(_UpwEqu2DwnChkii.m()==0) { //compute it if necessary
	 //-----------------------
	 //cerr<<"UpwEqu2DwnChk(FFT) compute"<<endl;	 //COMPUTE FFT	 //Index3 idx = ii.idx();	 
	 iA( idx.linfty()>1 );
	 DblNumMat denPos(dim(),1);	 for(int i=0; i<dim(); i++)		denPos(i,0) = double(idx(i))*2.0*R; //shift
	 DblNumMat chkPos(dim(),regPos().n());	 clear(chkPos);	 iC( daxpy(R, regPos(), chkPos) );
	 
	 DblNumMat tt(regPos().n()*trgDOF, srcDOF);
	 iC( _knl.kernel(denPos, denPos, chkPos, tt) );
	 // move data to tmp
	 
	 DblNumMat tmp(trgDOF,regPos().n()*srcDOF);
	 for(int k=0; k<regPos().n();k++) {
		for(int i=0; i<trgDOF; i++) {
		  for(int j=0; j<srcDOF; j++) {
			 tmp(i,j+k*srcDOF) = tt(i+k*trgDOF,j);
		  }
		}
	 }
	 _UpwEqu2DwnChkii.resize(trgDOF*srcDOF, effNum); //_memused[4] += dof*dof*effNum(UE)*sizeof(double);
	 //forward FFT from tmp to _UpwEqu2DwnChkii;
#ifdef FFTW3
	 int nnn[3];  nnn[0] = 2*np;  nnn[1] = 2*np;  nnn[2] = 2*np;
	 //fftw_execute_dft_r2c(_ue2dcplan, tmp.data(), (fftw_complex*)(_UpwEqu2DwnChkii.data()));
	 fftw_plan forplan = fftw_plan_many_dft_r2c(3,nnn,srcDOF*trgDOF, tmp.data(),NULL, srcDOF*trgDOF,1, (fftw_complex*)(_UpwEqu2DwnChkii.data()),NULL, srcDOF*trgDOF,1, FFTW_ESTIMATE);
	 fftw_execute(forplan);
	 fftw_destroy_plan(forplan);
#else
	 rfftwnd_real_to_complex(_forplan, srcDOF*trgDOF, tmp.data(), srcDOF*trgDOF, 1, (fftw_complex*)(_UpwEqu2DwnChkii.data()), srcDOF*trgDOF, 1);
#endif
  }
  //---------matvec
  //std::cout << effVal << endl;

  //double nrmfc = 1.0/double(regPos().n());
  fftw_complex* matPtr  = (fftw_complex*)(_UpwEqu2DwnChkii.data());
  fftw_complex* denPtr = (fftw_complex*)(effDen.data());
  fftw_complex* chkPtr   = (fftw_complex*)(effVal.data());
  int matStp  = srcDOF*trgDOF;
  int denStp = srcDOF;
  int chkStp   = trgDOF;

  for(int i=0; i<trgDOF; i++) {
	 for(int j=0; j<srcDOF; j++) {
		int matOff = j*trgDOF + i;
		int denOff = j;
		int chkOff = i;
		iC( cptwvv(effNum/2, 0.0, matPtr+matOff, matStp, denPtr+denOff, denStp, chkPtr+chkOff, chkStp) );
	 }
  }
  denPtr = NULL;
  chkPtr = NULL;
  matPtr = NULL;
  return (0);
}

// ---------------------------------------------------------------------- 
int MatMgnt3d::localPos(int tp, Point3 center, double radius, DblNumMat& positions)
{
  const DblNumMat& bas = samPos(tp);
  positions.resize(dim(), bas.n());
  for(int i=0; i<dim(); i++)
	 for(int j=0; j<positions.n(); j++)
		positions(i,j) = center(i) + radius * bas(i,j);
  return (0);
}

// ---------------------------------------------------------------------- 
int MatMgnt3d::cptwvv(int n, double alpha, fftw_complex* x, int incrementX, fftw_complex* y, int incrementY, fftw_complex* z, int incrementZ)
{
  for(int i=0; i<n; i++) {
#ifdef FFTW3
    	 int iz = i*incrementZ;
	 int iy = i*incrementY;
	 int ix = i*incrementX;
    //iA(0); // Needs to be fixed
	 (*(z+iz))[0] += ( (*(x+ix))[0] * (*(y+iy))[0] - (*(x+ix))[1] * (*(y+iy))[1]);
	 (*(z+iz))[1] += ( (*(x+ix))[0] * (*(y+iy))[1] + (*(x+ix))[1] * (*(y+iy))[0]);
	 //double k1 = ((*x)[0] != 0.0) ? (*x)[0] * ((*y)[0] + (*y)[1]) : 0.0;
	 //double k2 = ((*y)[1] != 0.0) ? (*y)[1] * ((*x)[0] + (*x)[1]) : 0.0;
	 //double k3 = ((*y)[0] != 0.0) ? (*y)[0] * ((*x)[1] - (*x)[0]) : 0.0;
	 //(*z)[0] += (k1 - k2);
	 //(*z)[1] += (k1 + k3);
#else
	 int iz = i*incrementZ;
	 int iy = i*incrementY;
	 int ix = i*incrementX;
	 (*(z+iz)).re += ( (*(x+ix)).re * (*(y+iy)).re - (*(x+ix)).im * (*(y+iy)).im);
	 (*(z+iz)).im += ( (*(x+ix)).re * (*(y+iy)).im + (*(x+ix)).im * (*(y+iy)).re);
#endif
	 
  }  //iC( PetscLogFlops( 10*n ) );
  
  return 0;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
vector<MatMgnt3d> MatMgnt3d::_mmvec;

MatMgnt3d* MatMgnt3d::getmmptr(Kernel3d knl, int np)
{
  for(int i=0; i<_mmvec.size(); i++)
	 if(_mmvec[i].knl()==knl && _mmvec[i].np()==np)
		return &(_mmvec[i]);
  
  _mmvec.push_back( MatMgnt3d() );
  int last = _mmvec.size()-1;
  MatMgnt3d* tmp = &(_mmvec[last]); //get the last one
  tmp->knl() = knl;  tmp->np() = np;
  tmp->setup();
  return tmp;
}

int MatMgnt3d::UpwEqu2DwnChkCleanUp(){
  int l = 0; /* For now - hom only */
  iA(_hom);

  int begr = -3;
  int endr = 3;

  OffTns<DblNumMat>& _UpwEqu2DwnChk = (hom()==true) ? _upwEqu2DwnChk[0] : _upwEqu2DwnChk[l];

  if(_UpwEqu2DwnChk.m()==0)	 _UpwEqu2DwnChk.resize(7,7,7,-3,-3,-3);

  for (int i = begr; i <= endr; i++){
	 for (int j = begr; j <= endr; j++){
		for (int k = begr; k <= endr; k++){
		  if ((i == begr || i == endr) || (j == begr || j == endr) || (k == begr || k == endr)){
			 DblNumMat& _UpwEqu2DwnChkii = _UpwEqu2DwnChk(i,j,k);
			 _UpwEqu2DwnChkii.resize(0,0);
		  }
		}
	 }
  }
  
  return(0);
}
	 
int MatMgnt3d::UpwChk2UpwEquCleanup(int maxlev){
  /* one level for now - fix later by iterating */
  for (int l = 0; l < maxlev; l++){
	 DblNumMat& _UC2UE = (_hom==true) ? _upwChk2UpwEqu[0] : _upwChk2UpwEqu[l];
	 _UC2UE.resize(0,0);

	 NumTns<DblNumMat>& _UE2UC = (_hom==true) ? _upwEqu2UpwChk[0] : _upwEqu2UpwChk[l];
	 if (_UE2UC.m() != 0){
		for (int i = 0; i < 2; i++){
		  for (int j = 0; j < 2; j++){
			 for (int k = 0; k < 2; k++){
				DblNumMat& _UE2UCii = _UE2UC(i,j,k);
				_UE2UCii.resize(0,0);
			 }
		  }
		}
		_UE2UC.resize(0,0,0);
	 }
  }
  return(0);
}
