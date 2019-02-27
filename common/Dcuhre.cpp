#include "Dcuhre.hpp"
#include "common/memAlloc.hpp"

extern "C" {
  void func(long *ndim, double *pnt, long *Xnumfun, double *funvls);
}

double glb_xtarg, glb_ytarg, glb_ztarg;
int glb_kt;
double glb_lambda;
int glb_nk;
int glb_k;
int glb_sdof;
int glb_tdof;

using namespace std;

int Dcuhre::setup(map<string,string>& opts){
  map<string,string>::iterator mapIdx;
  std::cout << "DCUHRE SETUP" << std::endl;
  long wksze, maxpts, minpts, key;
  /* For the dynamic computation using dcuhre */
  mapIdx = opts.find("-" + prefix() + "minpts"); assert(mapIdx!=opts.end());
  { istringstream ss((*mapIdx).second);  ss>> minpts; }
  std::cerr << "minpts = " << minpts << std::endl;
  mapIdx = opts.find("-" + prefix() + "maxpts"); assert(mapIdx!=opts.end());
  { istringstream ss((*mapIdx).second);  ss>> maxpts; }
  std::cerr << "maxpts = " << maxpts << std::endl;
  mapIdx = opts.find("-" + prefix() + "epsabs"); assert(mapIdx!=opts.end());
  double epsabs; { istringstream ss((*mapIdx).second);  ss >>epsabs; }
  mapIdx = opts.find("-" + prefix() + "epsrel"); assert(mapIdx!=opts.end());
  double epsrel; { istringstream ss((*mapIdx).second);  ss >>epsrel; }
  mapIdx = opts.find("-" + prefix() + "worksize"); assert(mapIdx!=opts.end());
  { istringstream ss((*mapIdx).second);  ss>> wksze; }
  mapIdx = opts.find("-" + prefix() + "key"); assert(mapIdx!=opts.end());
  { istringstream ss((*mapIdx).second);  ss>> key; }

  /* This will need to change if kval != 4 and Nk != 20 at some point */
  /* Should probably only allocate if needed */

  //double eps = 1e-12; /* Slows things down but seems necessary */
  _ndim = 3;
  _numfun = 0;
  _minpts = minpts;
  _maxpts = maxpts;
  _epsabs = epsabs;
  _epsrel = epsrel;
  _Alim=vecalloc((short int) _ndim); /*allocate memory for ndim double vector*/
  _Blim=vecalloc((short int) _ndim);
  _Alim[0] = -1.0; _Alim[1] = -1.0; _Alim[2] = -1.0;
  _Blim[0] = 1.0; _Blim[1] = 1.0; _Blim[2] = 1.0;
  _key = key;
  _restar = 0;
  _worksize = wksze;
  _isalloc = false;
  return(0);
}

Dcuhre::~Dcuhre(){
  if (_isalloc){
	 free(_Alim); free(_Blim);
	 free(_result);
	 free(_abserr);
	 free(_work);
  }
}

void Dcuhre::setLims(double a){
  setLims(a,-a);
}

void Dcuhre::setLims(double a, double b){
  setLims(a,a,a,b,b,b);
}

void Dcuhre::setLims(double ax, double ay, double az, double bx, double by, double bz){
  _Alim[0] = ax; _Alim[1] = ay; _Alim[2] = az;
  _Blim[0] = bx; _Blim[1] = by; _Blim[2] = bz;
}

void Dcuhre::printLims(){
  for (int i = 0; i < 3; i++){
	 std::cout << "A(" << i << ") = " << _Alim[i] << "  B(" << i << ") = " << _Blim[i] << std::endl;
  }
}

void Dcuhre::alloc(){
  iA(_numfun != (long)0);
  iA(_numfun == _nk*_sdof*_tdof);
  iA(_nk == _k*(_k+1)*(_k+2)/6);
  _result=vecalloc((int) _numfun);  /* resultant approximation to the integral(s) */
  _abserr=vecalloc((int) _numfun);  /* An array of estimates for absolute errors */
  _work=vecalloc((int) _worksize);
  _isalloc = true;
} 

int Dcuhre::eval(double x, double y, double z){
  if (!_isalloc) alloc();
  glb_xtarg = x; glb_ytarg = y; glb_ztarg = z;
  glb_kt = _kt;
  glb_lambda = _lambda;
  glb_nk = _nk;
  glb_k = _k;
  glb_sdof = _sdof;
  glb_tdof = _tdof;

  DDCUHRE(&_ndim,&_numfun,_Alim,_Blim,&_minpts,&_maxpts,(void *)func, &_epsabs, &_epsrel,&_key,&_worksize,&_restar,_result,_abserr,&_neval, &_ifail,_work);

  /*
  for (int i = 0; i < 3; i++){ 
	 std::cout << _Alim[i] << " " << _Blim[i] << " ";
  }
  std::cout << std::endl;
  std::cout << _ifail << std::endl;
  */
  return (int)_ifail;
}
