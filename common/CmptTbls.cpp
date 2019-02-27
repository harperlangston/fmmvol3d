#include "CmptTbls.hpp"

using namespace std;

CmptTbls::~CmptTbls(){
  delete _dcuhre;
}

int CmptTbls::setup(map<string,string>& optionsMap){
  map<string,string>::iterator mapIdx;
  std::cout << "CMPTTBLS SETUP" << std::endl;

  mapIdx = optionsMap.find("-" + prefix() + "ksrcval"); assert(mapIdx!=optionsMap.end());
  { istringstream ss((*mapIdx).second);  ss>> _ksrcval; }
  mapIdx = optionsMap.find("-" + prefix() + "ktrgval"); assert(mapIdx!=optionsMap.end());
  { istringstream ss((*mapIdx).second);  ss>> _ktrgval; }
  mapIdx = optionsMap.find("-" + prefix() + "np"); assert(mapIdx!=optionsMap.end());
  { istringstream ss((*mapIdx).second);  ss >>_np; }
  mapIdx = optionsMap.find("-" + prefix() + "kt"); assert(mapIdx!=optionsMap.end());
  { istringstream ss((*mapIdx).second);  ss >>_kt; }
  mapIdx = optionsMap.find("-" + prefix() + "lambda"); assert(mapIdx!=optionsMap.end());
  { istringstream ss((*mapIdx).second);  ss >>_lambda; }
  
  vector<double> coefs(2); coefs[0] = 1.0; coefs[1] = 1.0;
  Kernel3d knl(_kt, coefs);

  if (!(knl.homogeneous())){ knl.coefs()[1] = _lambda; }
	   
  _dcuhre = new Dcuhre(prefix()+"dcuhre_");
  iC( _dcuhre->setup(optionsMap) );	
  _dcuhre->kt() = knl.kernelType();

  int KSRC = _ksrcval; int srcNk = (KSRC+2)*(KSRC+1)*(KSRC)/6;
  _dcuhre->k() = KSRC; 
  _dcuhre->nk() = srcNk;
  _dcuhre->lambda() = _lambda;
  _dcuhre->numfun() = srcNk*knl.srcDOF()*knl.trgDOF();
  _dcuhre->sdof() = knl.srcDOF();
  _dcuhre->tdof() = knl.trgDOF();
 
  return(0);
}

int CmptTbls::eval(int s2m, int nbrs, int wx, int level, DblNumMat& _S2M, map<int, DblNumMat>& _SC2TV_map, map<int, DblNumMat>& _WXC2TV_map){
  int KSRC = _ksrcval; int srcNk = (KSRC+2)*(KSRC+1)*(KSRC)/6;
  int KTRG = _ktrgval; int trgNk = (KTRG+2)*(KTRG+1)*(KTRG)/6;

  //#warning FORCE KSRC=KTRG for now - alt. hardly used
  iA(KSRC == KTRG);
  
  DblNumMat _grdTrgSamPos;
  bool regular = (KTRG <= 6);
  iC( grdSamPosCal(regular, KTRG, 1.0, _grdTrgSamPos, false));

  int trgGrdSze = KTRG*KTRG*KTRG;
  int srcGrdSze = KSRC*KSRC*KSRC;

  map<int, map<int, DblNumMat> > _nrNbrPreCompF;
  map<int, map<int, DblNumMat> > _wxNbrPreCompF;

  vector<double> coefs(2); coefs[0] = 1.0; coefs[1] = 1.0;
  Kernel3d knl(_kt, coefs);

  /* Doesn't matter that done after the fact here */
  if (!(knl.homogeneous())){ knl.coefs()[1] = _lambda; }
  
  if (s2m){
	 double num;
	 string str;
	 stringstream converter;
	 ofstream myFileOut;

	 converter << "../include/" << knl.kernelType() << "/";
	 {
		ifstream myFileIn; myFileIn.open((converter.str()).c_str(), ios::in | ios::binary);
		if (!myFileIn) {
		  system((string("mkdir ").append(converter.str())).c_str());
		}
		myFileIn.close();
	 }
	 converter << KSRC << "/";
	 {
		ifstream myFileIn; myFileIn.open((converter.str()).c_str(), ios::in | ios::binary);
		if (!myFileIn) {
		  system((string("mkdir ").append(converter.str())).c_str());
		}
		myFileIn.close();
	 }
	 if (knl.homogeneous()) {
		converter <<"S2M/";
		{
		  ifstream myFileIn; myFileIn.open((converter.str()).c_str(), ios::in | ios::binary);
		  if (!myFileIn) {
			 system((string("mkdir ").append(converter.str())).c_str());
		  }
		  myFileIn.close();
		}
		converter <<"S2M_NP" << _np;
	 }
	 else {
		converter << _lambda << "/";
		{
		  ifstream myFileIn; myFileIn.open((converter.str()).c_str(), ios::in | ios::binary);
		  if (!myFileIn) {
			 system((string("mkdir ").append(converter.str())).c_str());
		  }
		  myFileIn.close();
		}
		converter << level << "/";
		{
		  ifstream myFileIn; myFileIn.open((converter.str()).c_str(), ios::in | ios::binary);
		  if (!myFileIn) {
			 system((string("mkdir ").append(converter.str())).c_str());
		  }
		  myFileIn.close();
		}
		converter << "S2M/";
		{
		  ifstream myFileIn; myFileIn.open((converter.str()).c_str(), ios::in | ios::binary);
		  if (!myFileIn) {
			 system((string("mkdir ").append(converter.str())).c_str());
		  }
		  myFileIn.close();
		}
		converter << "S2M_NP" << _np;
	 }
	 str.append(converter.str());
    ifstream myFileIn; myFileIn.open(str.c_str(), ios::in | ios::binary);
	 int NK = srcNk;
	 _S2M.resize(npPlnDatSze(UC,_np,knl.srcDOF(),knl.trgDOF()), NK*knl.trgDOF());
	 DblNumMat untPos(knl.dim(),0); clear(untPos);
	 double R = 1.0/pow(2.0,level);
	 samPosCal(_np+2,3.0,untPos,UC);
	 DblNumMat grdPos(untPos.m(),untPos.n()); clear(grdPos);
	 iC( daxpy(R, untPos, grdPos) ); //scale
	 iA(grdPos.n() == npPlnDatSze(UC,_np,knl.srcDOF(),knl.trgDOF())/knl.srcDOF());
	 cerr << "S2M LVL # " << level << " " << str.c_str() <<  endl;
	 if (!myFileIn) {
		//Cannot open file to read in, so create out
		myFileOut.open(str.c_str(), ios::out | ios::binary);	
		if (!myFileOut) { std::cerr << "S2M cannot open file " << str.c_str() << std::endl; exit(-1); }
		double scale = 1.0/pow(2.0,level);
		double init = -1.0*scale;
		_dcuhre->setLims(init); /* see Dcuhre.cpp */
		//_dcuhre->printLims();
		for (int j = 0; j < npPlnDatSze(UC,_np,knl.srcDOF(),knl.trgDOF())/knl.srcDOF(); j++){
		  double xt = grdPos(0,j); double yt = grdPos(1,j); double zt = grdPos(2,j);
		  _dcuhre->eval(xt,yt,zt);
		  int cnt = 0;
		  for (int k = 0; k < NK; k++){ for (int ds = 0; ds < knl.srcDOF(); ds++){ for (int dt = 0; dt < knl.trgDOF(); dt++){
				  /* Computed values go from src coeffs to chk vals */
				  _S2M(j*knl.srcDOF() + ds,k*knl.trgDOF() + dt) = _dcuhre->res(cnt);
				  //cerr << xt << " " << yt << " " << zt << " " << k << " " << _dcuhre->res(cnt) << endl;
				  cnt++;
				} } }
		}
	 
		myFileOut.write((char *) &(*(_S2M).data()), sizeof(*(_S2M).data())*grdPos.n()*NK*knl.srcDOF()*knl.trgDOF());
		//std::cerr << _S2M << endl;
		myFileOut.close();
	 }
	 else {
		myFileIn.read((char *)& *(_S2M).data(), sizeof(*(_S2M).data())*npPlnDatSze(UC,_np,knl.srcDOF(),knl.trgDOF())*NK*knl.trgDOF());
	 }
  }
  if (nbrs){ /* Compute Norm, Fine and Coarse Nbr Tables */
	 string str; stringstream converter;
	 converter << "../include/" << knl.kernelType() << "/";
	 {
		ifstream myFileIn; myFileIn.open((converter.str()).c_str(), ios::in | ios::binary);
		if (!myFileIn) {
		  system((string("mkdir ").append(converter.str())).c_str());
		}
		myFileIn.close();
	 }
	 converter << KSRC << "/";
	 {
		ifstream myFileIn; myFileIn.open((converter.str()).c_str(), ios::in | ios::binary);
		if (!myFileIn) {
		  system((string("mkdir ").append(converter.str())).c_str());
		}
		myFileIn.close();
	 }
	 if (knl.homogeneous()) {
		converter << "NBRFE16_" << KTRG;
	 }
	 else {
		converter << _lambda << "/";
		{
		  ifstream myFileIn; myFileIn.open((converter.str()).c_str(), ios::in | ios::binary);
		  if (!myFileIn) {
			 system((string("mkdir ").append(converter.str())).c_str());
		  }
		  myFileIn.close();
		}
		converter << level << "/";
		{
		  ifstream myFileIn; myFileIn.open((converter.str()).c_str(), ios::in | ios::binary);
		  if (!myFileIn) {
			 system((string("mkdir ").append(converter.str())).c_str());
		  }
		  myFileIn.close();
		}
		converter << "NBRFE16_" << KTRG;
	 }
	 str.append(converter.str());

	 ifstream myFileIn; myFileIn.open(str.c_str(), ios::in | ios::binary);
	 cerr << str.c_str() << endl;
	 
	 int NK = srcNk;
	 for (int type = NORM; type <= CRSE; type++){
		int star = ((type == NORM) ? N_SYM_CORNER : ((type == FINE) ? F_SYM_CORNER : C_SYM_CORNER));
		int stop = ((type == NORM) ? N_SYM_SELF : ((type == FINE) ? F_SYM_FACE : C_SYM_FACE));
		double scale = 1.0/pow(2.0,level);
		double init = -1.0*scale;
		_dcuhre->setLims(init); /* see Dcuhre.cpp */
		double rad = scale;
		/* Once init pts and radius set, change scales */
		if (type == FINE) scale *= 0.5;
		if (type == CRSE) scale *= 2.0;
		for (int nbr = star; nbr <= stop; nbr++){ /* See include/Syms.hpp*/
		  //cout << "NBR (lev,type,nb) " << level << " " << type << " " << nbr << endl;
		  DblNumMat& _SC2TV = _SC2TV_map[nbr]; iA (_SC2TV.m() == 0); /* If not, something wrong */
		  _SC2TV.resize(trgGrdSze * knl.srcDOF(),NK * knl.trgDOF());

		  Point3 ctr(0.0,0.0,0.0);
		  nbrSymMapSclVrt(type,nbr,rad,ctr);
		  DblNumMat grdPos(knl.dim(),trgGrdSze); clear(grdPos);

		  if (!myFileIn){
			 //Cannot open file, so create it
			 ofstream myFileOut; myFileOut.open(str.c_str(), ios::out | ios::app | ios::binary);	
			 if (!myFileOut) { std::cerr << "NBR cannot open file " << str.c_str() << std::endl; exit(-1); }
			 for (int l = 0; l < grdPos.n(); l++){ for (int d = 0; d < knl.dim(); d++){ grdPos(d,l) = ctr(d); }	 }
			 daxpy(scale, _grdTrgSamPos, grdPos);

			 for (int g = 0; g < grdPos.n(); g++){
				double xt = grdPos(0,g); double yt = grdPos(1,g); double zt = grdPos(2,g);
				_dcuhre->eval(xt,yt,zt);
				int cnt = 0;
				for (int c = 0; c < NK; c++){ for (int s = 0; s < knl.srcDOF(); s++){ for (int t = 0; t < knl.trgDOF(); t++){
						_SC2TV(g*knl.srcDOF() + s, c*knl.trgDOF() + t) = _dcuhre->res(cnt);
						cnt++;
					 } } }
			 }

			 myFileOut.write((char *) &(*(_SC2TV).data()), sizeof(*(_SC2TV).data())*grdPos.n()*NK*knl.srcDOF()*knl.trgDOF());
			 myFileOut.close();
		  }
		  else {
			 myFileIn.read((char *)& *(_SC2TV).data(), sizeof(*(_SC2TV).data())*grdPos.n()*NK*knl.srcDOF()*knl.trgDOF());
		  }
		}
	 }
	 myFileIn.close();
  }
  if (wx){
	 string str; stringstream converter;
	 converter << "../include/" << knl.kernelType() << "/";
	 {
		ifstream myFileIn; myFileIn.open((converter.str()).c_str(), ios::in | ios::binary);
		if (!myFileIn) {
		  system((string("mkdir ").append(converter.str())).c_str());
		}
		myFileIn.close();
	 }
	 converter << KSRC << "/";
	 {
		ifstream myFileIn; myFileIn.open((converter.str()).c_str(), ios::in | ios::binary);
		if (!myFileIn) {
		  system((string("mkdir ").append(converter.str())).c_str());
		}
		myFileIn.close();
	 }
	 if (knl.homogeneous()) {
		converter << "WXFE16_" << KTRG;
	 }
	 else {
		converter << _lambda << "/";
		{
		  ifstream myFileIn; myFileIn.open((converter.str()).c_str(), ios::in | ios::binary);
		  if (!myFileIn) {
			 system((string("mkdir ").append(converter.str())).c_str());
		  }
		  myFileIn.close();
		}
		converter << level << "/";
		{
		  ifstream myFileIn; myFileIn.open((converter.str()).c_str(), ios::in | ios::binary);
		  if (!myFileIn) {
			 system((string("mkdir ").append(converter.str())).c_str());
		  }
		  myFileIn.close();
		}
		converter << "WXFE16_" << KTRG;
	 }
	 str.append(converter.str());

	 ifstream myFileIn; myFileIn.open(str.c_str(), ios::in | ios::binary);
	 cerr << str.c_str() << endl;

	 //map<int, DblNumMat>& _WXC2TV_map  = _wxNbrPreCompF[level];
	 int NK = srcNk;
	 for (int type = WNODE; type <= XNODE; type++){
		int star = WX_SYM_C0;
		int stp = WX_SYM_F0;
		double scale = 1.0/pow(2.0,level);
		double init = -1.0*scale;
		_dcuhre->setLims(init); /* see Dcuhre.cpp */
		double rad = scale;
		/* Once init pts and radius set, change scales */
		if (type == WNODE) scale *= 0.5;
		if (type == XNODE) scale *= 2.0;
		for (int nbr = star; nbr <= stp; nbr++){ /* See include/Syms.hpp*/
		  int nbr_adj = (type == WNODE) ? nbr : nbr + 6;
		  //cout << "WX (lev,type,nb) " << level << " " << type << " " << nbr << endl;
		  DblNumMat& _WXC2TV = _WXC2TV_map[nbr_adj]; iA (_WXC2TV.m() == 0); /* If not, something wrong */
		  _WXC2TV.resize(trgGrdSze * knl.srcDOF(),NK * knl.trgDOF());
		  Point3 ctr(0.0,0.0,0.0);
		  nbrSymMapSclVrt(type,nbr,rad,ctr);
		  DblNumMat grdPos(knl.dim(),trgGrdSze); clear(grdPos);
		  
		  if (!myFileIn){
			 //Cannot open file, so create it
			 ofstream myFileOut; myFileOut.open(str.c_str(), ios::out | ios::app | ios::binary);	
			 if (!myFileOut) { std::cerr << "WX cannot open file " << str.c_str() << std::endl; exit(-1); }
			 for (int l = 0; l < grdPos.n(); l++){ for (int d = 0; d < knl.dim(); d++){ grdPos(d,l) = ctr(d); }	 }
			 daxpy(scale, _grdTrgSamPos, grdPos);
			 for (int g = 0; g < grdPos.n(); g++){
				double xt = grdPos(0,g); double yt = grdPos(1,g); double zt = grdPos(2,g);
				_dcuhre->eval(xt,yt,zt);
				int cnt = 0;
				for (int c = 0; c < NK; c++){ for (int s = 0; s < knl.srcDOF(); s++){ for (int t = 0; t < knl.trgDOF(); t++){
						_WXC2TV(g*knl.srcDOF() + s, c*knl.trgDOF() + t) = _dcuhre->res(cnt);
						cnt++;
					 } } }
			 }
			 myFileOut.write((char *) &(*(_WXC2TV).data()), sizeof(*(_WXC2TV).data())*grdPos.n()*NK*knl.srcDOF()*knl.trgDOF());
			 myFileOut.close();
		  }
		  else {
			 myFileIn.read((char *)& *(_WXC2TV).data(), sizeof(*(_WXC2TV).data())*grdPos.n()*NK*knl.srcDOF()*knl.trgDOF());
		  }
		}
	 }
	 myFileIn.close();
  }
  
  return(0);
}

int CmptTbls::eval_unbal(const int type, const int id, const int levgni, const double radgni, const Point3 ctrgni, const int levnbr, const double radnbr, const Point3 ctrnbr, const Point3 offset, DblNumMat& _evalMat){
  int KSRC = _ksrcval; int srcNk = (KSRC+2)*(KSRC+1)*(KSRC)/6;
  int KTRG = _ktrgval; int trgNk = (KTRG+2)*(KTRG+1)*(KTRG)/6;
  int NK = srcNk;
  
  //#warning FORCE KSRC=KTRG for now - alt. hardly used
  iA(KSRC == KTRG);
  
  DblNumMat _grdTrgSamPos;
  bool regular = (KTRG <= 6);
  iC( grdSamPosCal(regular, KTRG, 1.0, _grdTrgSamPos, false));

  int trgGrdSze = KTRG*KTRG*KTRG;
  int srcGrdSze = KSRC*KSRC*KSRC;
  vector<double> coefs(2); coefs[0] = 1.0; coefs[1] = 1.0;
  Kernel3d knl(_kt, coefs);
  /* Doesn't matter that done after the fact here */
  if (!(knl.homogeneous())){ knl.coefs()[1] = _lambda; }

  string str; stringstream converter;
  converter << "../include/" << knl.kernelType() << "/";
  {
	 ifstream myFileIn; myFileIn.open((converter.str()).c_str(), ios::in | ios::binary);
	 if (!myFileIn) {
		system((string("mkdir ").append(converter.str())).c_str());
	 }
	 myFileIn.close();
  }
  converter << KSRC << "/";
  if (!knl.homogeneous()) {
	 converter << _lambda << "/";
	 {
		ifstream myFileIn; myFileIn.open((converter.str()).c_str(), ios::in | ios::binary);
		if (!myFileIn) {
		  system((string("mkdir ").append(converter.str())).c_str());
		}
		myFileIn.close();
	 }
  }
  {
	 ifstream myFileIn; myFileIn.open((converter.str()).c_str(), ios::in | ios::binary);
	 if (!myFileIn) {
		system((string("mkdir ").append(converter.str())).c_str());
	 }
	 myFileIn.close();
  }
  converter << type << "/";
  {
	 ifstream myFileIn; myFileIn.open((converter.str()).c_str(), ios::in | ios::binary);
	 if (!myFileIn) {
		system((string("mkdir ").append(converter.str())).c_str());
	 }
	 myFileIn.close();
  }
  converter << abs(levgni-levnbr) << "/";
  {
	 ifstream myFileIn; myFileIn.open((converter.str()).c_str(), ios::in | ios::binary);
	 if (!myFileIn) {
		system((string("mkdir ").append(converter.str())).c_str());
	 }
	 myFileIn.close();
  }
  if (!knl.homogeneous()) {
	 converter << levgni << "/";
	 {
		ifstream myFileIn; myFileIn.open((converter.str()).c_str(), ios::in | ios::binary);
		if (!myFileIn) {
		  system((string("mkdir ").append(converter.str())).c_str());
		}
		myFileIn.close();
	 }
	 converter << levnbr << "/";
	 {
		ifstream myFileIn; myFileIn.open((converter.str()).c_str(), ios::in | ios::binary);
		if (!myFileIn) {
		  system((string("mkdir ").append(converter.str())).c_str());
		}
		myFileIn.close();
	 }
  }

  /* Reverse gni and nbr - way the code is written to generate unique number */
  //int look_up = lookup(type, levnbr, radnbr, ctrnbr, levgni, radgni, ctrgni, offset);
  int look_up = id;
  converter << "NBRFE16_" << KTRG << "_" << look_up;

  str.append(converter.str());

  ifstream myFileIn; myFileIn.open(str.c_str(), ios::in | ios::binary);
  cerr << str.c_str() << endl;
  int level = (knl.homogeneous() ? 0 : levnbr);
  double scale = 1.0/pow(2.0,level);
  double init = -1.0*scale;
  _dcuhre->setLims(init); /* see Dcuhre.cpp */
  double rad = scale;
  /* Once init pts and radius set, change scales */
  if (type == FINE || type == WNODE) scale *= pow(0.5,abs(levgni-levnbr));
  if (type == CRSE || type == XNODE) scale *= pow(2.0,abs(levgni-levnbr));
  
  _evalMat.resize(trgGrdSze * knl.srcDOF(),NK * knl.trgDOF());

  DblNumMat grdPos(knl.dim(),trgGrdSze); clear(grdPos);
  Point3 ctr((ctrgni-(ctrnbr+offset))/radnbr);
  
  if (!myFileIn){
	 //Cannot open file, so create it
	 ofstream myFileOut; myFileOut.open(str.c_str(), ios::out | ios::app | ios::binary);	
	 if (!myFileOut) { std::cerr << "NBR cannot open file " << str.c_str() << std::endl; exit(-1); }
	 for (int l = 0; l < grdPos.n(); l++){ for (int d = 0; d < knl.dim(); d++){ grdPos(d,l) = ctr(d); }	 }
	 daxpy(scale, _grdTrgSamPos, grdPos);

	 for (int g = 0; g < grdPos.n(); g++){
		double xt = grdPos(0,g); double yt = grdPos(1,g); double zt = grdPos(2,g);
		_dcuhre->eval(xt,yt,zt);
		int cnt = 0;
		for (int c = 0; c < NK; c++){ for (int s = 0; s < knl.srcDOF(); s++){ for (int t = 0; t < knl.trgDOF(); t++){
				_evalMat(g*knl.srcDOF() + s, c*knl.trgDOF() + t) = _dcuhre->res(cnt);
				cnt++;
			 } } }
	 }
	 myFileOut.write((char *) &(*(_evalMat).data()), sizeof(*(_evalMat).data())*grdPos.n()*NK*knl.srcDOF()*knl.trgDOF());
	 myFileOut.close();
  }
  else {
	 myFileIn.read((char *)& *(_evalMat).data(), sizeof(*(_evalMat).data())*grdPos.n()*NK*knl.srcDOF()*knl.trgDOF());
  }
  return(0);
}

int CmptTbls::lookup(const int type, const int levgni, const double radgni, const Point3 ctrgni, const int levnbr, const double radnbr, const Point3 ctrnbr, const Point3 offset){
  int look_up;
  
  if (type == FINE){
	 double radnb = radnbr; double radgn = radgni;
	 Point3 ctrnb(ctrnbr);  Point3 ctrgn(ctrgni + offset);
	 Point3 init(ctrgn);  init -= (radgn + radnb);
	 int n = (int)((radgn + radnb)/radnb) + 1;
	 iA( n == ((int)(pow(2.0, abs(levgni - levnbr))) + 2));
	 Index3 idx; for (int d = 0; d < 3; d++)  idx(d) = (int)((ctrnb(d) - init(d))/(2*radnb));
	 look_up = (idx(2)*n + idx(1))*n + idx(0);
  }
  else if (type == CRSE){
	 int n = (int)(pow(2.0, abs(levgni - levnbr))) + 2;
	 look_up  = ((n*n*n)-1) - lookup(FINE, levnbr, radnbr, ctrnbr, levgni, radgni, ctrgni, -offset);
  }
  else if (type == WNODE){
	 double radnb = radnbr; double radgn = radgni;
	 Point3 ctrnb(ctrnbr);  Point3 ctrgn(ctrgni + offset);
	 Point3 init(ctrgn);  init -= (radgn + 3.0*radnb);

	 int n = (int)((radgn + 3.0*radnb)/radnb) + 1;
	 iA( n == ((int)(pow(2.0, abs(levgni - levnbr))) + 4));
	 
	 Index3 idx; for (int d = 0; d < 3; d++)  idx(d) = (int)((ctrnb(d) - init(d))/(2*radnb));
	 look_up = (idx(2)*n + idx(1))*n + idx(0);
  }
  else if (type == XNODE){
	 int n = (int)(pow(2.0, abs(levgni - levnbr))) + 4;
	 look_up  = ((n*n*n)-1) - lookup(WNODE, levgni, radgni, ctrgni, levnbr, radnbr, ctrnbr, -offset);
  }
  else {
	 iA(0);
  }

  return look_up;
}
