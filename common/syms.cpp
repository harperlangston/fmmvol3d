#include "syms.hpp"

int npPlnDatSze(int tp, int np, int sdof, int tdof)
{
  int nn = (tp == UC ? np+2 : np);
  int tmpr = nn*nn*nn - (nn-2)*(nn-2)*(nn-2);
  if(tp==UE || tp ==DE)
	 return tmpr*sdof;
  else
	 return tmpr*tdof;
}

int npEffDatSze(int tp, int np, int sdof, int tdof)
{
  int effNum = (2*(np)+2)*(2*(np))*(2*(np));
  if(tp==UE || tp==DE)
	 return effNum*sdof;
  else
	 return effNum*tdof;
}

int samPosCal(int np, double R, DblNumMat& positions, int type)
{
  int n = np*np*np - (np-2)*(np-2)*(np-2);
  positions.resize(3,n);
  double step = 2.0/(np-1);
  double init = -1.0;
  int cnt = 0;
  for(int i=0; i<np; i++)
    for(int j=0; j<np; j++)
      for(int k=0; k<np; k++) {
		  if(i==0 || i==np-1 || j==0 || j==np-1 || k==0 || k==np-1) {
			 double x = init + i*step;
			 double y = init + j*step;
			 double z = init + k*step;

			 positions(0,cnt) = R*x;
			 positions(1,cnt) = R*y;
			 positions(2,cnt) = R*z;
			 cnt++;
		  }
      }

  iA(cnt==n);
  return 0;
}

// ---------------------------------------------------------------------- 
int regPosCal(int np, double R, DblNumMat& positions)
{
  int n = 2*np*2*np*2*np;
  positions.resize(3, n);
  double step = 2.0/(np-1);
  int cnt = 0;
  for(int k=0; k<2*np; k++)
	 for(int j=0; j<2*np; j++)
		for(int i=0; i<2*np; i++) {
		  int gi = (i<np) ? i : i-2*np;
		  int gj = (j<np) ? j : j-2*np;
		  int gk = (k<np) ? k : k-2*np;
		  positions(0, cnt) = R * gi*step;
		  positions(1, cnt) = R * gj*step;
		  positions(2, cnt) = R * gk*step;
		  cnt ++;
		}
  iA(cnt==n);
  return 0;
}

int grdSamPosCal(const bool regular, const int kval, const double R, DblNumMat& positions, const bool edges){
  int K = kval;
  double h = 2.0/((double)(K)); double hh = h/2.0; double init = -1.0+hh;
  if (edges) {
	 K = kval+2;
	 h = 2.0/((double)(K-1)); hh = 0.0; init = -1.0;
  }
  int n = K*K*K;
  positions.resize(3,n);


  int cnt = 0;
  if (regular){
	 for(int i = 1; i <= K; i++){ double z = init +(i-1)*h;
		for(int j = 1; j <= K; j++){ double y = init + (j-1)*h;
		  for(int k = 1;  k <= K; k++) { double x = init + (k-1)*h;
			 positions(0,cnt) = R*x; positions(1,cnt) = R*y; positions(2,cnt) = R*z; cnt++;
		  } } }
  }
  else {
	 double scale = 1.0; 
	 if (edges){
		scale = (0.0 + cos(((double)(2*(1) - 1))/((double)(2*K)) * M_PI));
	 }
	 for (int i = 1; i <= K; i++){ double z = -(0.0 + cos(((double)(2*(i) - 1))/((double)(2*K)) * M_PI))/scale;
		for (int j = 1; j <= K; j++){ double y = -(0.0 + cos(((double)(2*(j) - 1))/((double)(2*K)) * M_PI))/scale;
		  for (int k = 1; k <= K; k++){ double x = -(0.0 + cos(((double)(2*(k) - 1))/((double)(2*K)) * M_PI))/scale;
			 positions(0,cnt) = R*x; positions(1,cnt) = R*y; positions(2,cnt) = R*z; cnt++;
		  } } }
  }
  iA(cnt==n);
  return(0);
}

int symType(int nbrTyp, int tabNbr){
  int rtyp = -1;
  int nrmNbrSym[27] = {N_CORNER, N_EDGE, N_CORNER, N_EDGE, N_FACE, N_EDGE,
							  N_CORNER, N_EDGE, N_CORNER, N_EDGE, N_FACE, N_EDGE,
							  N_FACE, N_SELF, N_FACE, N_EDGE, N_FACE, N_EDGE,
							  N_CORNER, N_EDGE, N_CORNER, N_EDGE, N_FACE, N_EDGE,
							  N_CORNER, N_EDGE, N_CORNER};
  /* Coarse neighbors defined in the same fasion as Fine neighbors */
  int fnCrsNbrSym[56] = {FC_CORNER, FC_EDGE, FC_EDGE, FC_CORNER, FC_EDGE, FC_FACE,
								 FC_FACE, FC_EDGE, FC_EDGE, FC_FACE, FC_FACE, FC_EDGE,
								 FC_CORNER, FC_EDGE, FC_EDGE, FC_CORNER, FC_EDGE, FC_FACE,
								 FC_FACE, FC_EDGE, FC_FACE, FC_FACE, FC_FACE, FC_FACE,
								 FC_EDGE, FC_FACE, FC_FACE, FC_EDGE, FC_EDGE, FC_FACE,
								 FC_FACE, FC_EDGE, FC_FACE, FC_FACE, FC_FACE, FC_FACE,
								 FC_EDGE, FC_FACE, FC_FACE, FC_EDGE, FC_CORNER, FC_EDGE,
								 FC_EDGE, FC_CORNER, FC_EDGE, FC_FACE, FC_FACE, FC_EDGE,
								 FC_EDGE, FC_FACE, FC_FACE, FC_EDGE, FC_CORNER, FC_EDGE,
								 FC_EDGE, FC_CORNER};
 if (nbrTyp == NORM) rtyp =  nrmNbrSym[tabNbr];
 else if (nbrTyp == FINE || nbrTyp == CRSE) rtyp = fnCrsNbrSym[tabNbr];
 else { iA(0); }
 return rtyp;
}

int WXsymRef(int symNum){
  int symRef;
  if (symNum == WX_CORNER_0) symRef = WX_SYM_C0;
  else if (symNum == WX_CORNER_1) symRef = WX_SYM_C1;
  else if (symNum == WX_CORNER_2) symRef = WX_SYM_C2;
  else if (symNum == WX_EDGE_0) symRef = WX_SYM_E0;
  else if (symNum == WX_EDGE_1) symRef = WX_SYM_E1;
  else if (symNum == WX_FACE_0) symRef = WX_SYM_F0;
  else { iA(0); }
  
  return symRef;
}

int nbrSymRef(int type, int symNum){
  if (type == FINE){
	 if (symNum == FC_CORNER) return F_SYM_CORNER;
	 else if (symNum == FC_EDGE) return F_SYM_EDGE;
	 else return F_SYM_FACE;
  }
  else if (type == CRSE){
	 if (symNum == FC_CORNER) return C_SYM_CORNER;
	 else if (symNum == FC_EDGE) return C_SYM_EDGE;
	 else return C_SYM_FACE;
  }
  else {
	 iA( type == NORM );
	 if (symNum == N_CORNER) return N_SYM_CORNER;
	 else if (symNum == N_EDGE) return N_SYM_EDGE;
	 else if (symNum == N_FACE) return N_SYM_FACE;
	 else return N_SYM_SELF;
  }
}
																	
int symRef(int type, int symNum){
  if (type == FINE || type == CRSE || type == NORM){ return nbrSymRef(type,symNum); }
  else { return WXsymRef(symNum); }
}

int posNeg(Index3 &bdx){
  if (bdx(0) == 1){
	 if (bdx(1) == 1){ return ((bdx(2) == 1) ? 0 : 1); }
	 else { return((bdx(2) == 1) ? 2 : 3); }
  }
  else {
	 if (bdx(1) == 1){ return ((bdx(2) == 1) ? 4 : 5); }
	 else { return ((bdx(2) == 1) ? 6 : 7); }
  }
}

double posNeg(int posneg, int coeff, int knl, int ds, int dt){
  /* Different combinations for Positive and Negative for BasisSyms */
  /* (+, +, +)
	* (+, +, -)
	* (+, -, +)
	* (+, -, -)
	* (-, +, +)
	* (-, +, -)
	* (-, -, +)
	* (-, -, -)
	*/
  int symsPosNeg[8][120] = {
	 {POS, POS, POS, POS,   POS, POS, POS, POS, POS, POS,   POS, POS, POS, POS, POS, POS, POS, POS, POS, POS,   POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS,   POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS,   POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS,   POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS},
	 {POS, POS, POS, NEG,   POS, POS, NEG, POS, NEG, POS,   POS, POS, NEG, POS, NEG, POS, POS, NEG, POS, NEG,   POS, POS, NEG, POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS,   POS, POS, NEG, POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG,   POS, POS, NEG, POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS,   POS, POS, NEG, POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS, NEG},
	 {POS, POS, NEG, POS,   POS, NEG, POS, POS, NEG, POS,   POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS,   POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS, POS, NEG, POS, NEG, POS,   POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS,   POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS,   POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS},
	 {POS, POS, NEG, NEG,   POS, NEG, NEG, POS, POS, POS,   POS, NEG, NEG, POS, POS, POS, NEG, NEG, NEG, NEG,   POS, NEG, NEG, POS, POS, POS, NEG, NEG, NEG, NEG, POS, POS, POS, POS, POS,   POS, NEG, NEG, POS, POS, POS, NEG, NEG, NEG, NEG, POS, POS, POS, POS, POS, NEG, NEG, NEG, NEG, NEG, NEG,   POS, NEG, NEG, POS, POS, POS, NEG, NEG, NEG, NEG, POS, POS, POS, POS, POS, NEG, NEG, NEG, NEG, NEG, NEG, POS, POS, POS, POS, POS, POS, POS,   POS, NEG, NEG, POS, POS, POS, NEG, NEG, NEG, NEG, POS, POS, POS, POS, POS, NEG, NEG, NEG, NEG, NEG, NEG, POS, POS, POS, POS, POS, POS, POS, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG},
	 {POS, NEG, POS, POS,   POS, NEG, NEG, POS, POS, POS,   NEG, POS, POS, NEG, NEG, NEG, POS, POS, POS, POS,   POS, NEG, NEG, POS, POS, POS, NEG, NEG, NEG, NEG, POS, POS, POS, POS, POS,   NEG, POS, POS, NEG, NEG, NEG, POS, POS, POS, POS, NEG, NEG, NEG, NEG, NEG, POS, POS, POS, POS, POS, POS,   POS, NEG, NEG, POS, POS, POS, NEG, NEG, NEG, NEG, POS, POS, POS, POS, POS, NEG, NEG, NEG, NEG, NEG, NEG, POS, POS, POS, POS, POS, POS, POS,   NEG, POS, POS, NEG, NEG, NEG, POS, POS, POS, POS, NEG, NEG, NEG, NEG, NEG, POS, POS, POS, POS, POS, POS, NEG, NEG, NEG, NEG, NEG, NEG, NEG, POS, POS, POS, POS, POS, POS, POS, POS},
	 {POS, NEG, POS, NEG,   POS, NEG, POS, POS, NEG, POS,   NEG, POS, NEG, NEG, POS, NEG, POS, NEG, POS, NEG,   POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS, POS, NEG, POS, NEG, POS,   NEG, POS, NEG, NEG, POS, NEG, POS, NEG, POS, NEG, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG,   POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS,   NEG, POS, NEG, NEG, POS, NEG, POS, NEG, POS, NEG, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG},
	 {POS, NEG, NEG, POS,   POS, POS, NEG, POS, NEG, POS,   NEG, NEG, POS, NEG, POS, NEG, NEG, POS, NEG, POS,   POS, POS, NEG, POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS,   NEG, NEG, POS, NEG, POS, NEG, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, NEG, POS, NEG, POS, NEG, POS,   POS, POS, NEG, POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS,   NEG, NEG, POS, NEG, POS, NEG, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, POS, NEG, NEG, POS, NEG, POS, NEG, POS, NEG, POS},
	 {POS, NEG, NEG, NEG,   POS, POS, POS, POS, POS, POS,   NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG,   POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS,   NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG,   POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS, POS,   NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG}
  };
  iA( posneg < 8 && coeff < 120);
  double val = symsPosNeg[posneg][coeff];
  if (knl == 111 || knl == 211){
	 iA (ds == 0 && dt == 0);
  }
  else if (knl == 311 || knl == 411){
	 iA ((ds >= 0 && dt >= 0) && (ds <= 2 && dt <= 2));
	 int stkSrcTrgSymsPosNeg[8][3] = {
		{POS, POS, POS},
		{POS, POS, NEG},
		{POS, NEG, POS},
		{POS, NEG, NEG},
		{NEG, POS, POS},
		{NEG, POS, NEG},
		{NEG, NEG, POS},
		{NEG, NEG, NEG}	
	 };
	 val *=  stkSrcTrgSymsPosNeg[posneg][ds];
	 val *=  stkSrcTrgSymsPosNeg[posneg][dt];
  }
  return val;
}

int basRef(Index3 &cIdx){
  if (cIdx(0) == POSX || cIdx(0) == NEGX) {
	 if (cIdx(1) == POSY || cIdx(1) == NEGY){ return 0; }
	 else { return 1; }
  }
  else if (cIdx(0) == POSY || cIdx(0) == NEGY) {
	 if (cIdx(1) == POSX || cIdx(1) == NEGX){ return 2; }
	 else { return 3; }
  }
  else {
	 if (cIdx(1) == POSX || cIdx(1) == NEGX){ return 4; }
	 else { return 5; }
  }
}

int symNbrGetRefs(int nbr, int type, int &bref, int &symNum, int &posneg, Index3 &idx, Index3 &cIdx){

  if (type == NORM){
	 int nrSymVals[27][6] = {
		{ 1, 2, 3, 1, 2, 3}, { 1, 2, 3, 1, 2, 3}, { 1, 2, -3, 1, 2, -3},
		{ 1, 3, 2, 1, 3, 2}, { 1, 2, 3, 1, 2, 3}, { 1, -3, -2, 1, -3, -2},
		{ 1, -2, 3, 1, -2, 3}, { 1, -2, 3, 1, -2, 3}, { 1, -2, -3, 1, -2, -3},
		{ 3, 2, 1, 3, 2, 1}, { 2, 1, 3, 2, 1, 3}, { -3, 2, -1, -3, 2, -1},
		{ 3, 2, 1, 3, 2, 1}, { 1, 2, 3, 1, 2, 3}, { -3, 2, -1, -3, 2, -1},
		{ 3, -2, 1, 3, -2, 1}, { -2, -1, 3, -2, -1, 3}, { -3, -2, -1, -3, -2, -1},
		{ -1, 2, 3, -1, 2, 3}, { -1, 2, 3, -1, 2, 3}, { -1, 2, -3, -1, 2, -3},
		{ -1, 3, 2, -1, 3, 2}, { -1, 2, 3, -1, 2, 3}, { -1, -3, -2, -1, -3, -2},
		{ -1, -2, 3, -1, -2, 3}, { -1, -2, 3, -1, -2, 3}, { -1, -2, -3, -1, -2, -3}
	 };
	 idx(0) = nrSymVals[nbr][0]; idx(1) = nrSymVals[nbr][1]; idx(2) = nrSymVals[nbr][2];
	 cIdx(0) = nrSymVals[nbr][3]; cIdx(1) = nrSymVals[nbr][4]; cIdx(2) = nrSymVals[nbr][5];
	 
  }
  else {
	 int nrSymVals[56][6] = {
		{ 1, 2, 3, 1, 2, 3},       { 1, 2, 3, 1, 2, 3},      { 1, 2, -3, 1, 2, -3},
		{ 1, 2, -3, 1, 2, -3},     { 1, 3, 2, 1, 3, 2},      { 1, 2, 3, 1, 2, 3},
		{ 1, 2, -3, 1, 2, -3},     { 1, -3, 2, 1, 3, -2},    { 1, 3, -2, 1, -3, 2},
		{ 1, -2, 3, 1, -2, 3},     { 1, -2, -3, 1, -2, -3},  { 1, -3, -2, 1, -3, -2},
		{ 1, -2, 3, 1, -2, 3},     { 1, -2, 3, 1, -2, 3},    { 1, -2, -3, 1, -2, -3},
		{ 1, -2, -3, 1, -2, -3},   { 3, 2, 1, 3, 2, 1},      { 2, 1, 3, 2, 1, 3},
		{ 2, 1, -3, 2, 1, -3},     { -3, 2, 1, 3, 2, -1},    { 3, 2, 1, 3, 2, 1},
		{ -3, 2, 1, 3, 2, -1},     { 3, -2, 1, 3, -2, 1},    { -3, -2, 1, 3, -2, -1},
		{ 3, -2, 1, 3, -2, 1},     { -2, 1, 3, 2, -1, 3},    { -2, 1, -3, 2, -1, -3},
		{ -3, -2, 1, 3, -2, -1},   { 3, 2, -1, -3, 2, 1},    { 2, -1, 3, -2, 1, 3},
		{ 2, -1, -3, -2, 1, -3},   { -3, 2, -1, -3, 2, -1},  { 3, 2, -1, -3, 2, 1},
		{ -3, 2, -1, -3, 2, -1},   { 3, -2, -1, -3, -2, 1},  { -3, -2, -1, -3, -2, -1},
		{ 3, -2, -1, -3, -2, 1},   { -2, -1, 3, -2, -1, 3},  { -2, -1, -3, -2, -1, -3},
		{ -3, -2, -1, -3, -2, -1}, { -1, 2, 3, -1, 2, 3},    { -1, 2, 3, -1, 2, 3},
		{ -1, 2, -3, -1, 2, -3},   { -1, 2, -3, -1, 2, -3},  { -1, 3, 2, -1, 3, 2},
		{ -1, 2, 3, -1, 2, 3},     { -1, 2, -3, -1, 2, -3},  { -1, -3, 2, -1, 3, -2},
		{ -1, 3, -2, -1, -3, 2},   { -1, -2, 3, -1, -2, 3},  { -1, -2, -3, -1, -2, -3},
		{ -1, -3, -2, -1, -3, -2}, { -1, -2, 3, -1, -2, 3},  { -1, -2, 3, -1, -2, 3},
		{ -1, -2, -3, -1, -2, -3}, { -1, -2, -3, -1, -2, -3} };
	 idx(0) = nrSymVals[nbr][0]; idx(1) = nrSymVals[nbr][1]; idx(2) = nrSymVals[nbr][2];
	 cIdx(0) = nrSymVals[nbr][3]; cIdx(1) = nrSymVals[nbr][4]; cIdx(2) = nrSymVals[nbr][5];
  }

  symNum = symType(type,nbr);
  bref = basRef(cIdx); 
    
  Index3 bdx;
  bdx(0) = (cIdx(0) > 0 ? 1 : -1);
  bdx(1) = (cIdx(1) > 0 ? 1 : -1);
  bdx(2) = (cIdx(2) > 0 ? 1 : -1);

  posneg = posNeg(bdx);
  return(0);
}

int rvrsIdx(int i, int base){
  iA(base >= i);
  return (base-1) - i;
}

Index3 pntRef(Index3 &idx, int i, int j, int k, int base){
  int x, y, z;
  
  if (idx(0) == POSX) { x = i; }
  else if (idx(0) == NEGX) { x = rvrsIdx(i, base); }
  else if (idx(0) == POSY) { x = j; }
  else if (idx(0) == NEGY) { x = rvrsIdx(j, base); }
  else if (idx(0) == POSZ) { x = k; }
  else { x = rvrsIdx(k, base); }

  if (idx(1) == POSY) { y = j; }
  else if (idx(1) == NEGY) { y = rvrsIdx(j, base); }
  else if (idx(1) == POSX) { y = i; }
  else if (idx(1) == NEGX) { y = rvrsIdx(i, base); }
  else if (idx(1) == POSZ) { y = k; }
  else { y = rvrsIdx(k, base); }

  if (idx(2) == POSZ) { z = k; }
  else if (idx(2) == NEGZ) { z = rvrsIdx(k, base); }
  else if (idx(2) == POSX) { z = i; }
  else if (idx(2) == NEGX) { z = rvrsIdx(i, base); }
  else if (idx(2) == POSY) { z = j; }
  else { z = rvrsIdx(j, base); }

  return Index3(x,y,z);
}

int nbrSymMapSclVrt (int type, int nbr, double rad, Point3 &ctr){
  int nbrLoc;
  int sclVarNbrSymMap[10] = {N_CORNER, N_EDGE, N_FACE, N_SELF, FC_CORNER, FC_EDGE, FC_FACE, FC_CORNER, FC_EDGE, FC_FACE};
  int sclVarWXSymMap[12] = {WX_CORNER_0, WX_CORNER_1, WX_CORNER_2, WX_EDGE_0, WX_EDGE_1, WX_FACE_0, WX_CORNER_0, WX_CORNER_1, WX_CORNER_2, WX_EDGE_0, WX_EDGE_1, WX_FACE_0};
  
  if (type == NORM){
	 nbrLoc = sclVarNbrSymMap[nbr];
	 if      (nbrLoc == N_CORNER) { ctr(0) = 0.0 - 2.0*rad; ctr(1) = 0.0 - 2.0*rad; ctr(2) = 0.0 - 2.0*rad; }
	 else if (nbrLoc == N_EDGE)   { ctr(0) = 0.0 - 2.0*rad; ctr(1) = 0.0 - 2.0*rad; ctr(2) = 0.0; }
	 else if (nbrLoc == N_FACE)   { ctr(0) = 0.0 - 2.0*rad; ctr(1) = 0.0;           ctr(2) = 0.0; }
	 else if (nbrLoc == N_SELF)   { ctr(0) = 0.0;           ctr(1) = 0.0;           ctr(2) = 0.0; }
	 else { iA(0); }
  }
  else if (type == FINE){
	 nbrLoc = sclVarNbrSymMap[nbr];
	 if (nbrLoc == FC_CORNER)     { ctr(0) = 0.0 - 1.5*rad; ctr(1) = 0.0 - 1.5*rad; ctr(2) = 0.0 - 1.5*rad; }
	 else if (nbrLoc == FC_EDGE)  { ctr(0) = 0.0 - 1.5*rad; ctr(1) = 0.0 - 1.5*rad; ctr(2) = 0.0 - 0.5*rad; }
	 else if (nbrLoc == FC_FACE)  { ctr(0) = 0.0 - 1.5*rad; ctr(1) = 0.0 - 0.5*rad; ctr(2) = 0.0 - 0.5*rad; }
	 else { iA(0); }
  }
  else if (type == CRSE){
	 nbrLoc = sclVarNbrSymMap[nbr];
	 if (nbrLoc == FC_CORNER)     { ctr(0) = 0.0 - 3.0*rad; ctr(1) = 0.0 - 3.0*rad; ctr(2) = 0.0 - 3.0*rad; }
	 else if (nbrLoc == FC_EDGE)  { ctr(0) = 0.0 - 3.0*rad; ctr(1) = 0.0 - 3.0*rad; ctr(2) = 0.0 - 1.0*rad; }
	 else if (nbrLoc == FC_FACE)  { ctr(0) = 0.0 - 3.0*rad; ctr(1) = 0.0 - 1.0*rad; ctr(2) = 0.0 - 1.0*rad; }
	 else { iA(0); }
  }
  else if (type == WNODE){
	 nbrLoc = sclVarWXSymMap[nbr];
	 if (nbrLoc == WX_CORNER_0)     { ctr(0) = 0.0 - 2.5*rad; ctr(1) = 0.0 - 2.5*rad; ctr(2) = 0.0 - 2.5*rad; }
	 else if (nbrLoc == WX_CORNER_1){ ctr(0) = 0.0 - 2.5*rad; ctr(1) = 0.0 - 2.5*rad; ctr(2) = 0.0 - 1.5*rad; }
	 else if (nbrLoc == WX_CORNER_2){ ctr(0) = 0.0 - 2.5*rad; ctr(1) = 0.0 - 1.5*rad; ctr(2) = 0.0 - 1.5*rad; }
	 else if (nbrLoc == WX_EDGE_0)  { ctr(0) = 0.0 - 2.5*rad; ctr(1) = 0.0 - 2.5*rad; ctr(2) = 0.0 - 0.5*rad; }
	 else if (nbrLoc == WX_EDGE_1)  { ctr(0) = 0.0 - 2.5*rad; ctr(1) = 0.0 - 1.5*rad; ctr(2) = 0.0 - 0.5*rad; }
	 else if (nbrLoc == WX_FACE_0)  { ctr(0) = 0.0 - 2.5*rad; ctr(1) = 0.0 - 0.5*rad; ctr(2) = 0.0 - 0.5*rad; }
	 else { iA(0); }
  }
  else if (type == XNODE){
	 nbrLoc = sclVarWXSymMap[nbr];
	 if (nbrLoc == WX_CORNER_0)     { ctr(0) = 0.0 - 5.0*rad; ctr(1) = 0.0 - 5.0*rad; ctr(2) = 0.0 - 5.0*rad; }
	 else if (nbrLoc == WX_CORNER_1){ ctr(0) = 0.0 - 5.0*rad; ctr(1) = 0.0 - 5.0*rad; ctr(2) = 0.0 - 3.0*rad; }
	 else if (nbrLoc == WX_CORNER_2){ ctr(0) = 0.0 - 5.0*rad; ctr(1) = 0.0 - 3.0*rad; ctr(2) = 0.0 - 3.0*rad; }
	 else if (nbrLoc == WX_EDGE_0)  { ctr(0) = 0.0 - 5.0*rad; ctr(1) = 0.0 - 5.0*rad; ctr(2) = 0.0 - 1.0*rad; }
	 else if (nbrLoc == WX_EDGE_1)  { ctr(0) = 0.0 - 5.0*rad; ctr(1) = 0.0 - 3.0*rad; ctr(2) = 0.0 - 1.0*rad; }
	 else if (nbrLoc == WX_FACE_0)  { ctr(0) = 0.0 - 5.0*rad; ctr(1) = 0.0 - 1.0*rad; ctr(2) = 0.0 - 1.0*rad; }
	 else { iA(0); }
  }
  else { iA(0); }

  return(0);
}


int basSym(int bref, int c, int knl){
  /* Different Basis Combos */
  /* (X,Y,Z),
	* (X,Z,Y),
	* (Y,X,Z),
	* (Y,Z,X),
	* (Z,X,Y),
	* (Z,Y,X)
	*/
  int basisSyms[6][120] = {
	 {B0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15, B16, B17, B18, B19, B20, B21, B22, B23, B24, B25, B26, B27, B28, B29, B30, B31, B32, B33, B34, B35, B36, B37, B38, B39, B40, B41, B42, B43, B44, B45, B46, B47, B48, B49, B50, B51, B52, B53, B54, B55, B56, B57, B58, B59, B60, B61, B62, B63, B64, B65, B66, B67, B68, B69, B70, B71, B72, B73, B74, B75, B76, B77, B78, B79, B80, B81, B82, B83, B84, B85, B86, B87, B88, B89, B90, B91, B92, B93, B94, B95, B96, B97, B98, B99, B100, B101, B102, B103, B104, B105, B106, B107, B108, B109, B110, B111, B112, B113, B114, B115, B116, B117, B118, B119},
	 {B0, B1, B3, B2, B4, B6, B5, B9, B8, B7, B10, B12, B11, B15, B14, B13, B19, B18, B17, B16, B20, B22, B21, B25, B24, B23, B29, B28, B27, B26, B34, B33, B32, B31, B30, B35, B37, B36, B40, B39, B38, B44, B43, B42, B41, B49, B48, B47, B46, B45, B55, B54, B53, B52, B51, B50, B56, B58, B57, B61, B60, B59, B65, B64, B63, B62, B70, B69, B68, B67, B66, B76, B75, B74, B73, B72, B71, B83, B82, B81, B80, B79, B78, B77, B84, B86, B85, B89, B88, B87, B93, B92, B91, B90, B98, B97, B96, B95, B94, B104, B103, B102, B101, B100, B99, B111, B110, B109, B108, B107, B106, B105, B119, B118, B117, B116, B115, B114, B113, B112},
	 {B0, B2, B1, B3, B7, B5, B8, B4, B6, B9, B16, B13, B17, B11, B14, B18, B10, B12, B15, B19, B30, B26, B31, B23, B27, B32, B21, B24, B28, B33, B20, B22, B25, B29, B34, B50, B45, B51, B41, B46, B52, B38, B42, B47, B53, B36, B39, B43, B48, B54, B35, B37, B40, B44, B49, B55, B77, B71, B78, B66, B72, B79, B62, B67, B73, B80, B59, B63, B68, B74, B81, B57, B60, B64, B69, B75, B82, B56, B58, B61, B65, B70, B76, B83, B112, B105, B113, B99, B106, B114, B94, B100, B107, B115, B90, B95, B101, B108, B116, B87, B91,  B96,  B102, B109, B117, B85,  B88,  B92,  B97,  B103, B110, B118, B84,  B86,  B89,  B93,  B98,  B104, B111, B119},
	 {B0, B2, B3, B1, B7, B8, B5, B9, B6, B4, B16, B17, B13, B18, B14, B11, B19, B15, B12, B10, B30, B31, B26, B32, B27, B23, B33, B28, B24, B21, B34, B29, B25, B22, B20, B50, B51, B45, B52, B46, B41, B53, B47, B42, B38, B54, B48, B43, B39, B36, B55, B49, B44, B40, B37, B35, B77, B78, B71, B79, B72, B66, B80, B73, B67, B62, B81, B74, B68, B63, B59, B82, B75, B69, B64, B60, B57, B83, B76, B70, B65, B61, B58, B56, B112, B113, B105, B114, B106, B99, B115, B107, B100, B94, B116, B108, B101, B95, B90, B117, B109, B102, B96, B91,  B87,  B118, B110, B103, B97,  B92,  B88,  B85,  B119, B111, B104, B98,  B93,  B89,  B86,  B84},
	 {B0, B3, B1, B2, B9, B6, B8, B4, B5, B7, B19, B15, B18, B12, B14, B17, B10, B11, B13, B16, B34, B29, B33, B25, B28, B32, B22, B24, B27, B31, B20, B21, B23, B26, B30, B55, B49, B54, B44, B48, B53, B40, B43, B47, B52, B37, B39, B42, B46, B51, B35, B36, B38, B41, B45, B50, B83, B76, B82, B70, B75, B81, B65, B69, B74, B80, B61, B64, B68, B73, B79, B58, B60, B63, B67, B72, B78, B56, B57, B59, B62, B66, B71, B77, B119, B111, B118, B104, B110, B117, B98, B103, B109, B116, B93, B97, B102, B108, B115, B89, B92,  B96,  B101, B107, B114, B86,  B88,  B91,  B95,  B100, B106, B113, B84,  B85,  B87,  B90,  B94,  B99,  B105, B112},
	 {B0, B3, B2, B1, B9, B8, B6, B7, B5, B4, B19, B18, B15, B17, B14, B12, B16, B13, B11, B10, B34, B33, B29, B32, B28, B25, B31, B27, B24, B22, B30, B26, B23, B21, B20, B55, B54, B49, B53, B48, B44, B52, B47, B43, B40, B51, B46, B42, B39, B37, B50, B45, B41, B38, B36, B35, B83, B82, B76, B81, B75, B70, B80, B74, B69, B65, B79, B73, B68, B64, B61, B78, B72, B67, B63, B60, B58, B77, B71, B66, B62, B59, B57, B56, B119, B118, B111, B117, B110, B104, B116, B109, B103, B98, B115, B108, B102, B97, B93, B114, B107, B101, B96,  B92,  B89,  B113, B106, B100, B95,  B91,  B88,  B86,  B112, B105, B99,  B94,  B90,  B87,  B85,  B84}
  };
  int val = basisSyms[bref][c];
  if (knl == 111 || knl == 211){
	 //val = basisSyms[bref][c];
  }
  else if (knl == 311 || knl == 411){
	 //val = basisSyms[bref][c];
  }
  return val;
}

int srcTrgSym(int bref, int dts, int knl){
  int stkSrcTrgSyms[6][3] = {
	 {0, 1, 2},
	 {0, 2, 1},
	 {1, 0, 2},
	 {1, 2, 0},
	 {2, 0, 1},
	 {2, 1, 0}
  };
  int val;
  if (knl == 311 || knl == 411){
	 iA( dts >= 0 && dts <= 2 );
	 val = stkSrcTrgSyms[bref][dts];
  }
  else { iA(0); }
  return val;
}

/* For W and X list symmetries */

int WXsymType(int tabNde){
  int WXNodeSym[152] = {WX_CORNER_0, WX_CORNER_1, WX_EDGE_0, WX_EDGE_0, WX_CORNER_1, WX_CORNER_0,
								WX_CORNER_1, WX_CORNER_2, WX_EDGE_1, WX_EDGE_1, WX_CORNER_2, WX_CORNER_1,
								WX_EDGE_0,   WX_EDGE_1,   WX_FACE_0, WX_FACE_0, WX_EDGE_1,   WX_EDGE_0,
								WX_EDGE_0,   WX_EDGE_1,   WX_FACE_0, WX_FACE_0, WX_EDGE_1,   WX_EDGE_0,
								WX_CORNER_1, WX_CORNER_2, WX_EDGE_1, WX_EDGE_1, WX_CORNER_2, WX_CORNER_1,
								WX_CORNER_0, WX_CORNER_1, WX_EDGE_0, WX_EDGE_0, WX_CORNER_1, WX_CORNER_0,

								WX_CORNER_1, WX_CORNER_2, WX_EDGE_1, WX_EDGE_1, WX_CORNER_2, WX_CORNER_1,
								WX_CORNER_2,      						                          WX_CORNER_2,
								WX_EDGE_1,                             					        WX_EDGE_1,
								WX_EDGE_1,                            					        WX_EDGE_1,
								WX_CORNER_2,                             					     WX_CORNER_2,
								WX_CORNER_1, WX_CORNER_2, WX_EDGE_1, WX_EDGE_1, WX_CORNER_2, WX_CORNER_1,

								WX_EDGE_0, WX_EDGE_1, WX_FACE_0, WX_FACE_0, WX_EDGE_1, WX_EDGE_0,
								WX_EDGE_1,                             					  WX_EDGE_1,
								WX_FACE_0,                             					  WX_FACE_0,
								WX_FACE_0,															  WX_FACE_0,
								WX_EDGE_1,                             					  WX_EDGE_1,
								WX_EDGE_0, WX_EDGE_1, WX_FACE_0, WX_FACE_0, WX_EDGE_1, WX_EDGE_0,

								WX_EDGE_0, WX_EDGE_1, WX_FACE_0, WX_FACE_0, WX_EDGE_1, WX_EDGE_0,
								WX_EDGE_1,                             					  WX_EDGE_1,
								WX_FACE_0,                             					  WX_FACE_0,
								WX_FACE_0,															  WX_FACE_0,
								WX_EDGE_1,                             					  WX_EDGE_1,
								WX_EDGE_0, WX_EDGE_1, WX_FACE_0, WX_FACE_0, WX_EDGE_1, WX_EDGE_0,

								WX_CORNER_1, WX_CORNER_2, WX_EDGE_1, WX_EDGE_1, WX_CORNER_2, WX_CORNER_1,
								WX_CORNER_2,      						                          WX_CORNER_2,
								WX_EDGE_1,                             					        WX_EDGE_1,
								WX_EDGE_1,                            					        WX_EDGE_1,
								WX_CORNER_2,                             					     WX_CORNER_2,
								WX_CORNER_1, WX_CORNER_2, WX_EDGE_1, WX_EDGE_1, WX_CORNER_2, WX_CORNER_1,

								WX_CORNER_0, WX_CORNER_1, WX_EDGE_0, WX_EDGE_0, WX_CORNER_1, WX_CORNER_0,
								WX_CORNER_1, WX_CORNER_2, WX_EDGE_1, WX_EDGE_1, WX_CORNER_2, WX_CORNER_1,
								WX_EDGE_0,   WX_EDGE_1,   WX_FACE_0, WX_FACE_0, WX_EDGE_1,   WX_EDGE_0,
								WX_EDGE_0,   WX_EDGE_1,   WX_FACE_0, WX_FACE_0, WX_EDGE_1,   WX_EDGE_0,
								WX_CORNER_1, WX_CORNER_2, WX_EDGE_1, WX_EDGE_1, WX_CORNER_2, WX_CORNER_1,
								WX_CORNER_0, WX_CORNER_1, WX_EDGE_0, WX_EDGE_0, WX_CORNER_1, WX_CORNER_0,
  };
  return WXNodeSym[tabNde];
}

int WXsymNodeGetRefs(int WXnode, int &bref, int &symNum, int &posneg, Index3 &idx, Index3& cIdx){
  symNum = WXsymType(WXnode);
  
  int wxSymVals[152][6] = {
	 { 1, 2, 3, 1, 2, 3}, { 1, 2, 3, 1, 2, 3}, { 1, 2, 3, 1, 2, 3},
	 { 1, 2, -3, 1, 2, -3}, { 1, 2, -3, 1, 2, -3}, { 1, 2, -3, 1, 2, -3},
	 { 1, 3, 2, 1, 3, 2}, { 1, 2, 3, 1, 2, 3}, { 1, 2, 3, 1, 2, 3},
	 { 1, 2, -3, 1, 2, -3}, { 1, 2, -3, 1, 2, -3}, { 1, -3, 2, 1, 3, -2},
	 { 1, 3, 2, 1, 3, 2}, { 1, 3, 2, 1, 3, 2}, { 1, 2, 3, 1, 2, 2},
	 { 1, 2, -3, 1, 2, -3},	{ 1, -3, 2, 1, 3, -2}, { 1, -3, 2, 1, 3, -2},
	 { 1, 3, -2, 1, -3, 2}, { 1, 3, -2, 1, -3, 2}, { 1, -2, 3, 1, -2, 2},
	 { 1, -2, -3, 1, -2, -3}, { 1, -3, -2, 1, -3, -2}, { 1, -3, -2, 1, -3, -2},
	 { 1, 3, -2, 1, -3, 2}, { 1, -2, 3, 1, -2, 3}, { 1, -2, 3, 1, -2, 3},
	 { 1, -2, -3, 1, -2, -3}, { 1, -2, -3, 1, -2, -3}, { 1, -3, -2, 1, -3, -2},
	 { 1, -2, 3, 1, -2, 3}, { 1, -2, 3, 1, -2, 3}, { 1, -2, 3, 1, -2, 3},
	 { 1, -2, -3, 1, -2, -3}, { 1, -2, -3, 1, -2, -3}, { 1, -2, -3, 1, -2, -3},
	 { 3, 2, 1, 3, 2, 1}, { 2, 1, 3, 2, 1, 3}, { 2, 1, 3, 2, 1, 3},
	 { 2, 1, -3, 2, 1, -3}, { 2, 1, -3, 2, 1, -3}, { -3, 2, 1, 3, 2, -1},
	 { 3, 2, 1, 3, 2, 1}, { -3, 2, 1, 3, 2, -1}, { 3, 1, 2, 2, 3, 1},
	 { -3, 1, 2, 2, 3, -1}, { 3, 1, -2, 2, -3, 1}, { -3, 1, -2, 2, -3, -1},
	 { 3, -2, 1, 3, -2, 1}, { -3, -2, 1, 3, -2, -1}, { 3, -2, 1, 3, -2, 1},
	 { -2, 1, 3, 2, -1, 3}, { -2, 1, 3, 2, -1, 3}, { -2, 1, -3, 2, -1, -3},
	 { -2, 1, -3, 2, -1, -3}, { -3, -2, 1, 3, -2, -1}, { 3, 2, 1, 3, 2, 1},
	 { 2, 3, 1, 3, 1, 2}, { 2, 1, 3, 2, 1, 3}, { 2, 1, -3, 2, 1, -3},
	 { 2, -3, 1, 3, 1, -2}, { -3, 2, 1, 3, 2, -1}, { 3, 2, 1, 3, 2, 1},
	 { -3, 2, 1, 3, 2, -1}, { 3, 2, 1, 3, 2, 1}, { -3, 2, 1, 3, 2, -1},
	 { 3, -2, 1, 3, -2, 1}, { -3, -2, 1, 3, -2, -1}, { 3, -2, 1, 3, -2, 1},
	 { -3, -2, 1, 3, -2, -1}, { 3, -2, 1, 3, -2, 1}, { -2, 3, 1, 3, -1, 2},
	 { -2, 1, 3, 2, -1, 3}, { -2, 1, -3, 2, -1, -3}, { -2, -3, 1, 3, -1, -2},
	 { -3, -2, 1, 3, -2, -1}, { 3, 2, -1, -3, 2, 1}, { 2, 3, -1, -3, 1, 2},
	 { 2, -1, 3, -2, 1, 3}, { 2, -1, -3, -2, 1, -3}, { 2, -3, -1, -3, 1, -2},
	 { -3, 2, -1, -3, 2, -1}, { 3, 2, -1, -3, 2, 1}, { -3, 2, -1, -3, 2, -1},
	 { 3, 2, -1, -3, 2, 1}, { -3, 2, -1, -3, 2, -1}, { 3, -2, -1, -3, -2, 1},
	 { -3, -2, -1, -3, -2, -1}, { 3, -2, -1, -3, -2, 1}, { -3, -2, -1, -3, -2, -1},
	 { 3, -2, -1, -3, -2, 1}, { -2, 3, -1, -3, -1, 2}, { -2, -1, 3, -2, -1, 3},
	 { -2, -1, -3, -2, -1, -3}, { -2, -3, -1, -3, -1, -2}, { -3, -2, -1, -3, -2, -1},
	 { 3, 2, -1, -3, 2, 1}, { 2, -1, 3, -2, 1, 3}, { 2, -1, 3, -2, 1, 3},
	 { 2, -1, -3, -2, 1, -3}, { 2, -1, -3, -2, 1, -3}, { -3, 2, -1, -3, 2, -1},
	 { 3, 2, -1, -3, 2, 1}, { -3, 2, -1, -3, 2, -1}, { 3, -1, 2, -2, 3, 1},
	 { -3, -1, 2, -2, 3, -1}, { 3, -1, -2, -2, -3, 1},	{ -3, -1, -2, -2, -3, -1},
	 { 3, -2, -1, -3, -2, 1}, { -3, -2, -1, -3, -2, -1},	{ 3, -2, -1, -3, -2, 1},
	 { -2, -1, 3, -2, -1, 3}, { -2, -1, 3, -2, -1, 3},	{ -2, -1, -3, -2, -1, -3},
	 { -2, -1, -3, -2, -1, -3},	{ -3, -2, -1, -3, -2, -1},	{ -1, 2, 3, -1, 2, 3},
	 { -1, 2, 3, -1, 2, 3},	{ -1, 2, 3, -1, 2, 3},	{ -1, 2, -3, -1, 2, -3},
	 { -1, 2, -3, -1, 2, -3},	{ -1, 2, -3, -1, 2, -3},	{ -1, 3, 2, -1, 3, 2},
	 { -1, 2, 3, -1, 2, 3},	{ -1, 2, 3, -1, 2, 3},	{ -1, 2, -3, -1, 2, -3},
	 { -1, 2, -3, -1, 2, -3},	{ -1, -3, 2, -1, 3, -2},	{ -1, 3, 2, -1, 3, 2},
	 { -1, 3, 2, -1, 3, 2},	{ -1, 2, 3, -1, 2, 2},	{ -1, 2, -3, -1, 2, -3},
	 { -1, -3, 2, -1, 3, -2},	{ -1, -3, 2, -1, 3, -2},	{ -1, 3, -2, -1, -3, 2},
	 { -1, 3, -2, -1, -3, 2},	{ -1, -2, 3, -1, -2, 2},	{ -1, -2, -3, -1, -2, -3},
	 { -1, -3, -2, -1, -3, -2}, { -1, -3, -2, -1, -3, -2}, { -1, 3, -2, -1, -3, 2},
	 { -1, -2, 3, -1, -2, 3}, { -1, -2, 3, -1, -2, 3}, { -1, -2, -3, -1, -2, -3},
	 { -1, -2, -3, -1, -2, -3}, { -1, -3, -2, -1, -3, -2}, { -1, -2, 3, -1, -2, 3},
	 { -1, -2, 3, -1, -2, 3}, { -1, -2, 3, -1, -2, 3}, { -1, -2, -3, -1, -2, -3},
	 { -1, -2, -3, -1, -2, -3}, { -1, -2, -3, -1, -2, -3}
  };

  idx(0) = wxSymVals[WXnode][0]; idx(1) = wxSymVals[WXnode][1]; idx(2) = wxSymVals[WXnode][2];
  cIdx(0) = wxSymVals[WXnode][3]; cIdx(1) = wxSymVals[WXnode][4]; cIdx(2) = wxSymVals[WXnode][5];

  bref = basRef(cIdx); 
    
  Index3 bdx;
  bdx(0) = (cIdx(0) > 0 ? 1 : -1);
  bdx(1) = (cIdx(1) > 0 ? 1 : -1);
  bdx(2) = (cIdx(2) > 0 ? 1 : -1);

  posneg = posNeg(bdx);
  return(0);
}

int WXsymNodeGetRefsDwnChk(int WXnode, int &bref, int &symNum, int &posneg, Index3 &idx, Index3 &cIdx){
  symNum = WXsymType(WXnode);
  iA(0);

  int wxSymVals[152][6] = {
	 { 1, 2, 3, 1, 2, 3}, { 1, 2, 3, 1, 2, 3}, { 1, 2, 3, 1, 2, 3},
	 { 1, 2, -3, 1, 2, -3}, { 1, 2, -3, 1, 2, -3}, { 1, 2, -3, 1, 2, -3},
	 { 1, 3, 2, 1, 3, 2}, { 1, 2, 3, 1, 2, 3}, { 1, 2, 3, 1, 2, 3},
	 { 1, 2, -3, 1, 2, -3}, { 1, 2, -3, 1, 2, -3}, { 1, -3, 2, 1, 3, -2},
	 { 1, 3, 2, 1, 3, 2}, { 1, 3, 2, 1, 3, 2}, { 1, 2, 3, 1, 2, 2},
	 { 1, 2, -3, 1, 2, -3},	{ 1, -3, 2, 1, 3, -2}, { 1, -3, 2, 1, 3, -2},
	 { 1, 3, -2, 1, -3, 2}, { 1, 3, -2, 1, -3, 2}, { 1, -2, 3, 1, -2, 2},
	 { 1, -2, -3, 1, -2, -3}, { 1, -3, -2, 1, -3, -2}, { 1, -3, -2, 1, -3, -2},
	 { 1, 3, -2, 1, -3, 2}, { 1, -2, 3, 1, -2, 3}, { 1, -2, 3, 1, -2, 3},
	 { 1, -2, -3, 1, -2, -3}, { 1, -2, -3, 1, -2, -3}, { 1, -3, -2, 1, -3, -2},
	 { 1, -2, 3, 1, -2, 3}, { 1, -2, 3, 1, -2, 3}, { 1, -2, 3, 1, -2, 3},
	 { 1, -2, -3, 1, -2, -3}, { 1, -2, -3, 1, -2, -3}, { 1, -2, -3, 1, -2, -3},
	 { 3, 2, 1, 3, 2, 1}, { 2, 1, 3, 2, 1, 3}, { 2, 1, 3, 2, 1, 3},
	 { 2, 1, -3, 2, 1, -3}, { 2, 1, -3, 2, 1, -3}, { -3, 2, 1, 3, 2, -1},
	 { 3, 2, 1, 3, 2, 1}, { -3, 2, 1, 3, 2, -1}, { 3, 1, 2, 2, 3, 1},
	 { -3, 1, 2, 2, 3, -1}, { 3, 1, -2, 2, -3, 1}, { -3, 1, -2, 2, -3, -1},
	 { 3, -2, 1, 3, -2, 1}, { -3, -2, 1, 3, -2, -1}, { 3, -2, 1, 3, -2, 1},
	 { -2, 1, 3, 2, -1, 3}, { -2, 1, 3, 2, -1, 3}, { -2, 1, -3, 2, -1, -3},
	 { -2, 1, -3, 2, -1, -3}, { -3, -2, 1, 3, -2, -1}, { 3, 2, 1, 3, 2, 1},
	 { 2, 3, 1, 3, 1, 2}, { 2, 1, 3, 2, 1, 3}, { 2, 1, -3, 2, 1, -3},
	 { 2, -3, 1, 3, 1, -2}, { -3, 2, 1, 3, 2, -1}, { 3, 2, 1, 3, 2, 1},
	 { -3, 2, 1, 3, 2, -1}, { 3, 2, 1, 3, 2, 1}, { -3, 2, 1, 3, 2, -1},
	 { 3, -2, 1, 3, -2, 1}, { -3, -2, 1, 3, -2, -1}, { 3, -2, 1, 3, -2, 1},
	 { -3, -2, 1, 3, -2, -1}, { 3, -2, 1, 3, -2, 1}, { -2, 3, 1, 3, -1, 2},
	 { -2, 1, 3, 2, -1, 3}, { -2, 1, -3, 2, -1, -3}, { -2, -3, 1, 3, -1, -2},
	 { -3, -2, 1, 3, -2, -1}, { 3, 2, -1, -3, 2, 1}, { 2, 3, -1, -3, 1, 2},
	 { 2, -1, 3, -2, 1, 3}, { 2, -1, -3, -2, 1, -3}, { 2, -3, -1, -3, 1, -2},
	 { -3, 2, -1, -3, 2, -1}, { 3, 2, -1, -3, 2, 1}, { -3, 2, -1, -3, 2, -1},
	 { 3, 2, -1, -3, 2, 1}, { -3, 2, -1, -3, 2, -1}, { 3, -2, -1, -3, -2, 1},
	 { -3, -2, -1, -3, -2, -1}, { 3, -2, -1, -3, -2, 1}, { -3, -2, -1, -3, -2, -1},
	 { 3, -2, -1, -3, -2, 1}, { -2, 3, -1, -3, -1, 2}, { -2, -1, 3, -2, -1, 3},
	 { -2, -1, -3, -2, -1, -3}, { -2, -3, -1, -3, -1, -2}, { -3, -2, -1, -3, -2, -1},
	 { 3, 2, -1, -3, 2, 1}, { 2, -1, 3, -2, 1, 3}, { 2, -1, 3, -2, 1, 3},
	 { 2, -1, -3, -2, 1, -3}, { 2, -1, -3, -2, 1, -3}, { -3, 2, -1, -3, 2, -1},
	 { 3, 2, -1, -3, 2, 1}, { -3, 2, -1, -3, 2, -1}, { 3, -1, 2, -2, 3, 1},
	 { -3, -1, 2, -2, 3, -1}, { 3, -1, -2, -2, -3, 1},	{ -3, -1, -2, -2, -3, -1},
	 { 3, -2, -1, -3, -2, 1}, { -3, -2, -1, -3, -2, -1},	{ 3, -2, -1, -3, -2, 1},
	 { -2, -1, 3, -2, -1, 3}, { -2, -1, 3, -2, -1, 3},	{ -2, -1, -3, -2, -1, -3},
	 { -2, -1, -3, -2, -1, -3},	{ -3, -2, -1, -3, -2, -1},	{ -1, 2, 3, -1, 2, 3},
	 { -1, 2, 3, -1, 2, 3},	{ -1, 2, 3, -1, 2, 3},	{ -1, 2, -3, -1, 2, -3},
	 { -1, 2, -3, -1, 2, -3},	{ -1, 2, -3, -1, 2, -3},	{ -1, 3, 2, -1, 3, 2},
	 { -1, 2, 3, -1, 2, 3},	{ -1, 2, 3, -1, 2, 3},	{ -1, 2, -3, -1, 2, -3},
	 { -1, 2, -3, -1, 2, -3},	{ -1, -3, 2, -1, 3, -2},	{ -1, 3, 2, -1, 3, 2},
	 { -1, 3, 2, -1, 3, 2},	{ -1, 2, 3, -1, 2, 2},	{ -1, 2, -3, -1, 2, -3},
	 { -1, -3, 2, -1, 3, -2},	{ -1, -3, 2, -1, 3, -2},	{ -1, 3, -2, -1, -3, 2},
	 { -1, 3, -2, -1, -3, 2},	{ -1, -2, 3, -1, -2, 2},	{ -1, -2, -3, -1, -2, -3},
	 { -1, -3, -2, -1, -3, -2}, { -1, -3, -2, -1, -3, -2}, { -1, 3, -2, -1, -3, 2},
	 { -1, -2, 3, -1, -2, 3}, { -1, -2, 3, -1, -2, 3}, { -1, -2, -3, -1, -2, -3},
	 { -1, -2, -3, -1, -2, -3}, { -1, -3, -2, -1, -3, -2}, { -1, -2, 3, -1, -2, 3},
	 { -1, -2, 3, -1, -2, 3}, { -1, -2, 3, -1, -2, 3}, { -1, -2, -3, -1, -2, -3},
	 { -1, -2, -3, -1, -2, -3}, { -1, -2, -3, -1, -2, -3}
  };
  
  idx(0) = wxSymVals[WXnode][2]; idx(1) = wxSymVals[WXnode][1]; idx(2) = wxSymVals[WXnode][0];
  cIdx(0) = wxSymVals[WXnode][3]; cIdx(1) = wxSymVals[WXnode][4]; cIdx(2) = wxSymVals[WXnode][5];
  
  bref = basRef(cIdx); 
    
  Index3 bdx;
  bdx(0) = (cIdx(0) > 0 ? 1 : -1);
  bdx(1) = (cIdx(1) > 0 ? 1 : -1);
  bdx(2) = (cIdx(2) > 0 ? 1 : -1);

  if (cIdx(0) == POSX || cIdx(0) == NEGX){
	 if (cIdx(0) == POSX){
		if (cIdx(1) == POSY || cIdx(1) == NEGY){
		  if (cIdx(1) == POSY){
			 if (cIdx(2) == NEGX){
				bdx(0) *= -1;
				bdx(2) *= -1;
			 }
		  }
		  else { }
		}
		else {
		  if (cIdx(1) == POSZ){
			 if (cIdx(2) == POSY){
			 }
			 else {
				bdx(1) *= -1;
				bdx(2) *= -1;
			 }
		  }
		  else {
			 if (cIdx(2) == POSY){
				bdx(1) *= -1;
				bdx(2) *= -1;
			 }
		  }
		}
	 }
	 else {
		if (cIdx(1) == POSY || cIdx(1) == NEGY){
		  if (cIdx(1) == POSY){
			 if (cIdx(2) == NEGX){
				bdx(0) *= -1;
				bdx(2) *= -1;
			 }
		  }
		  else { }
		}
		else {
		  if (cIdx(1) == POSZ){
			 if (cIdx(2) == NEGY	){
				bdx(1) *= -1;
				bdx(2) *= -1;
			 }
		  }
		  else {
			 if (cIdx(2) == POSY){
				bdx(1) *= -1;
				bdx(2) *= -1;
			 }
			 else {
			
			 }
		  }
		}
	 }
  }
  else if (cIdx(0) == POSY || cIdx(0) == NEGY){
	 if (cIdx(0) == POSY){
		if (cIdx(1) == POSX || cIdx(1) == NEGX) {
		  if (cIdx(1) == POSX){

		  }
		  else {
			 if (cIdx(2) == POSZ){
				bdx(0) *= -1;
				bdx(1) *= -1;
			 }
			 else {
				bdx(0) *= -1;
				bdx(1) *= -1;
			 }
		  }
		}
		else {
		  if (cIdx(1) == POSZ){
			 if (cIdx(2) == POSX){
			 }
			 else {
				bdx(0) *= -1;
				bdx(2) *= -1;
			 }
		  }
		  else {
			 if (cIdx(2) == POSX){
				bdx(1) *= -1;
				bdx(2) *= -1;
			 }
			 else {
				bdx(0) *= -1;
				bdx(1) *= -1;
			 }
		  }
		}
	 }
	 else {
		if (cIdx(1) == POSX || cIdx(1) == NEGX){
		  if (cIdx(1) == POSX){
			 if (cIdx(2) == POSY){
			 }
			 else {
				bdx(0) *= -1;
				bdx(1) *= -1;
			 }
		  }
		}
		else {
		  if (cIdx(1) == POSZ){
			 if (cIdx(2) == POSX){
				bdx(0) *= -1;
				bdx(1) *= -1;
			 }
			 else {
				bdx(1) *= -1;
				bdx(2) *= -1;
			 }
		  }
		  else {
			 if (cIdx(2) == POSX){
				bdx(0) *= -1;
				bdx(2) *= -1;
			 }
		  }
		}
	 }
  }
  else {
	 if (cIdx(0) == POSZ){
		if (cIdx(1) == POSY || cIdx(1) == NEGY){
		  if (cIdx(1) == POSY){
			 if (cIdx(2) == POSX){
			 }
			 else{
				bdx(0) *= -1;
				bdx(2) *= -1;
			 }
		  }
		  else {
			 if (cIdx(2) == POSX){
			 }
			 else {
				bdx(0) *= -1;
				bdx(2) *= -1;
			 }
		  }
		}
		else {
		  if (cIdx(1) == POSX){
			 if (cIdx(2) == POSY){
			 }
			 else {
				bdx(1) *= -1;
				bdx(2) *= -1;
			 }
		  }
		  else {
			 if (cIdx(2) == POSY){
				bdx(0) *= -1;
				bdx(1) *= -1;
			 }
			 else {
				bdx(0) *= -1;
				bdx(2) *= -1;
			 }
		  }
		}
	 }
	 else {
		if (cIdx(1) == POSY || cIdx(1) == NEGY){
		  if (cIdx(1) == POSY){
			 if (cIdx(2) == POSX){
				bdx(0) *= -1;
				bdx(2) *= -1;
			 }
		  }
		  else {
			 if (cIdx(2) == POSX){
				bdx(0) *= -1;
				bdx(2) *= -1;
			 }
			 else {
				//bdx(1) *= -1;
				//bdx(2) *= -1;
			 }
		  }
		}
		else {
		  if (cIdx(1) == POSX){
			 if (cIdx(2) == POSY){
				bdx(0) *= -1;
				bdx(2) *= -1;
			 }
			 else {
				bdx(0) *= -1;
				bdx(1) *= -1;
			 }
		  }
		  else {
			 if (cIdx(2) == POSY){
				bdx(1) *= -1;
				bdx(2) *= -1;
			 }
			 else {
			 }
		  }
		}
	 }
  }

		
  
  int pn = posNeg(bdx);
  posneg = pn;
  //if (pn == 1) posneg = 2;
}

int nrmNbr(int WXnode){
  int WXNodeNormNbrLoc[152] = {NRM0, NRM0, NRM1, NRM1, NRM2, NRM2,
										 NRM0, NRM0, NRM1, NRM1, NRM2, NRM2,
										 NRM3, NRM3, NRM4, NRM4, NRM5, NRM5,
										 NRM3, NRM3, NRM4, NRM4, NRM5, NRM5,
										 NRM6, NRM6, NRM7, NRM7, NRM8, NRM8,
										 NRM6, NRM6, NRM7, NRM7, NRM8, NRM8,
										 NRM0, NRM0, NRM1, NRM1, NRM2, NRM2,
										 NRM0,                         NRM2,
										 NRM3,                         NRM5,
										 NRM3,                         NRM5,
										 NRM6,                         NRM8,
										 NRM6, NRM6, NRM7, NRM7, NRM8, NRM8,

										 NRM9, NRM9, NRM10, NRM10, NRM11,   NRM11,
										 NRM9,                              NRM11,
										 NRM12,                             NRM14,
										 NRM12,                             NRM14,
										 NRM15,                             NRM17,
										 NRM15, NRM15, NRM16, NRM16, NRM17, NRM17,
										 NRM9,  NRM9, NRM10, NRM10, NRM11,  NRM11,
										 NRM9,                              NRM11,
										 NRM12,                             NRM14,
										 NRM12,                             NRM14,
										 NRM15,                             NRM17,
										 NRM15, NRM15, NRM16, NRM16, NRM17, NRM17,

										 NRM18, NRM18, NRM19, NRM19, NRM20, NRM20,
										 NRM18,                             NRM20,
										 NRM21,                             NRM23,
										 NRM21,                             NRM23,
										 NRM24,                             NRM26,
										 NRM24, NRM24, NRM25, NRM25, NRM26, NRM26,

										 NRM18, NRM18, NRM19, NRM19, NRM20, NRM20,
										 NRM18, NRM18, NRM19, NRM19, NRM20, NRM20,
										 NRM21, NRM21, NRM22, NRM22, NRM23, NRM23,
										 NRM21, NRM21, NRM22, NRM22, NRM23, NRM23,
										 NRM24, NRM24, NRM25, NRM25, NRM26, NRM26,
										 NRM24, NRM24, NRM25, NRM25, NRM26, NRM26
  };
	
  return WXNodeNormNbrLoc[WXnode];

}

int nrmNbrByDif(Index3 dif){
  int val = -1;
  if (dif(0) == 1) { /* Normal Nbr 0-8 */
	 if (dif(1) == 1) { /* Normal Nbr 0-2 */
		if (dif(2) == 1) val = NRM0;
		else if (dif(2) == 0) val = NRM1;
		else if (dif(2) == -1) val = NRM2;
		else iA( 0 );
	 }
	 else if (dif(1) == 0) { /* Normal Nbr 3-5 */
		if (dif(2) == 1) val = NRM3;
		else if (dif(2) == 0) val = NRM4;
		else if (dif(2) == -1) val = NRM5;
		else iA( 0 );
	 }
	 else if (dif(1) == -1) { /* Normal Nbr 3-5 */
		if (dif(2) == 1) val = NRM6;
		else if (dif(2) == 0) val = NRM7;
		else if (dif(2) == -1) val = NRM8;
		else iA( 0 );
	 }
  }
  else if (dif(0) == 0) { /* Normal Nbr 9-17 */
	 if (dif(1) == 1) { /* Normal Nbr 9-11 */
		if (dif(2) == 1) val = NRM9;
		else if (dif(2) == 0) val = NRM10;
		else if (dif(2) == -1) val = NRM11;
		else iA( 0 );
	 }
	 else if (dif(1) == 0) { /* Normal Nbr 12-14 */
		if (dif(2) == 1) val = NRM12;
		else if (dif(2) == 0) val = NRM13;
		else if (dif(2) == -1) val = NRM14;
		else iA( 0 );
	 }
	 else if (dif(1) == -1) { /* Normal Nbr 15-17 */
		if (dif(2) == 1) val = NRM15;
		else if (dif(2) == 0) val = NRM16;
		else if (dif(2) == -1) val = NRM17;
		else iA( 0 );
	 }
  }
  else if (dif(0) == -1) { /* Normal Nbr 18-26 */
	 if (dif(1) == 1) { /* Normal Nbr 18-20 */
		if (dif(2) == 1) val = NRM18;
		else if (dif(2) == 0) val = NRM19;
		else if (dif(2) == -1) val = NRM20;
		else iA( 0 );
	 }
	 else if (dif(1) == 0) { /* Normal Nbr 21-23 */
		if (dif(2) == 1) val = NRM21;
		else if (dif(2) == 0) val = NRM22;
		else if (dif(2) == -1) val = NRM23;
		else iA( 0 );
	 }
	 else if (dif(1) == -1) { /* Normal Nbr 24-26 */
		if (dif(2) == 1) val = NRM24;
		else if (dif(2) == 0) val = NRM25;
		else if (dif(2) == -1) val = NRM26;
		else iA( 0 );
	 }
  }	
  else iA ( 0 );
  return val;
}

int fnNbrByDif(Index3 dif, Index3 pDif){
  int val = -1;
  
  if (dif(0) == 1) { /* Fine Nbr 0-15 */
	 if (dif(1) == 1) { /* Fine Nbr 0-3 */
		if (dif(2) == 1){
		  if ( !(pDif(0) == 1 && pDif(1) == 1 && pDif(2) == 1 )) val = -1;
		  else val =  FN0;
		}
		else if (dif(2) == 0){
		  if( !(pDif(0) == 1 && pDif(1) == 1 )) val = -1;
		  else {
			 if (pDif(2) == 0) val =  FN1;
			 else if (pDif(2) == 1) val =  FN2;
			 else iA(0);
		  }
		}
		else if (dif(2) == -1){
		  if ( !(pDif(0) == 1 && pDif(1) == 1 && pDif(2) == 0 )) val = -1;
		  else val =  FN3;
		}
		else iA(0);
	 }
	 else if (dif(1) == 0) { /* Fine Nbr 4-11 */
		if (dif(2) == 1) {
		  if( !(pDif(0) == 1 && pDif(2) == 1 )) val = -1;
		  else {
			 if (pDif(1) == 0) val =  FN4;
			 else if (pDif(1) == 1) val =  FN8;
			 else iA(0);
		  }
		}
		else if (dif(2) == 0){
		  if( !(pDif(0) == 1 )) val = -1;
		  else {
			 if (pDif(1) == 0){
				if (pDif(2) == 0) val =  FN5;
				else if (pDif(2) == 1) val =  FN6;
				else iA(0);
			 }
			 else if (pDif(1) == 1){
				if (pDif(2) == 0) val =  FN9;
				else if (pDif(2) == 1) val =  FN10;
				else iA(0);
			 }
		  }
		}
		else if (dif(2) == -1){
		  if( !(pDif(0) == 1 && pDif(2) == 0 )) val = -1;
		  else {
			 if (pDif(1) == 0) val =  FN7;
			 else if (pDif(1) == 1) val =  FN11;
			 else iA(0);
		  }
		}
	 }
	 else if (dif(1) == -1) { /* Fine Nbr 12-15 */
		if (dif(2) == 1) {
		  if ( !(pDif(0) == 1 && pDif(1) == 0 && pDif(2) == 1 )) val = -1;
		  else val =  FN12;
		}
		else if (dif(2) == 0){
		  if( !(pDif(0) == 1 && pDif(1) == 0 )) val = -1;
		  else {
			 if (pDif(2) == 0) val =  FN13;
			 else if (pDif(2) == 1) val =  FN14;
			 else iA(0);
		  }
		}
		else if (dif(2) == -1){
		  if ( !(pDif(0) == 1 && pDif(1) == 0 && pDif(2) == 0 )) val = -1;
		  else val =  FN15;
		}
		else iA(0);
	 }
	 else iA(0);
  }
  else if (dif(0) == 0){ /* Fine Nbr 16-39 */
	 if (dif(1) == 1) { /* Fine Nbr 16-19 and 28-31 */
		if (dif(2) == 1){
		  if( !(pDif(1) == 1 && pDif(2) == 1 )) val = -1;
		  else {
			 if (pDif(0) == 0) val =  FN16;
			 else if (pDif(0) == 1) val =  FN28;
			 else iA(0);
		  }
		}
		else if (dif(2) == 0) {
		  if( !(pDif(1) == 1 )) val = -1;
		  else {
			 if (pDif(2) == 0) {
				if (pDif(0) == 0) val =  FN17;
				else if (pDif(0) == 1) val =  FN29;
				else iA(0);
			 }
			 else if (pDif(2) == 1) {
				if (pDif(0) == 0) val =  FN18;
				else if (pDif(0) == 1) val =  FN30;
				else iA(0);
			 }
			 else iA(0);
		  }
		}
		else if (dif(2) == -1){
		  if( !(pDif(1) == 1 && pDif(2) == 0 )) val = -1;
		  else {	
			 if (pDif(0) == 0) val =  FN19;
			 else if (pDif(0) == 1) val =  FN31;
			 else iA(0);
		  }
		}
	 }
	 else if (dif(1) == 0) { /* Fine Nbr 20-23 and 32-35 */
		if( !(dif(2) != 0)) val = -1;
		else {
		  if (dif(2) == 1){ /* Fine Nbr 20,22,32,34 */
			 if( !(pDif(2) == 1)) val = -1;
			 else {
				if (pDif(1) == 0){
				  if (pDif(0) == 0) val =  FN20;
				  else if (pDif(0) == 1) val =  FN32;
				  else iA(0);
				}
				else if (pDif(1) == 1){
				  if (pDif(0) == 0) val =  FN22;
				  else if (pDif(0) == 1) val =  FN34;
				  else iA(0);
				}
				else iA(0);
			 }
		  }
		  else if (dif(2) == -1){ /* Fine Nbr 21,23,33,35 */
			 if( !(pDif(2) == 0)) val = -1;
			 else {
				if (pDif(1) == 0){
				  if (pDif(0) == 0) val =  FN21;
				  else if (pDif(0) == 1) val =  FN33;
				  else iA(0);
				}
				else if (pDif(1) == 1){
				  if (pDif(0) == 0) val =  FN23;
				  else if (pDif(0) == 1) val =  FN35;
				  else iA(0);
				}
				else iA(0);
			 }
		  }
		  else iA(0);
		}
	 }
	 else if (dif(1) == -1) { /* Fine Nbr 24-27 and 36-39 */
		  if (dif(2) == 1){
			 if( !(pDif(1) == 0 && pDif(2) == 1 )) val = -1;
			 else {
				if (pDif(0) == 0) val =  FN24;
				else if (pDif(0) == 1) val =  FN36;
				else iA(0);
			 }
		  }
		  else if (dif(2) == 0) {
			 if( !(pDif(1) == 0 )) val = -1;
			 else {
				if (pDif(2) == 0) {
				  if (pDif(0) == 0) val =  FN25;
				  else if (pDif(0) == 1) val =  FN37;
				  else iA(0);
				}
				else if (pDif(2) == 1) {
				  if (pDif(0) == 0) val =  FN26;
				  else if (pDif(0) == 1) val =  FN38;
				  else iA(0);
				}
				else iA(0);
			 }
		  }
		  else if (dif(2) == -1){
			 if( !(pDif(1) == 0 && pDif(2) == 0 )) val = -1;
			 else {
				if (pDif(0) == 0) val =  FN27;
				else if (pDif(0) == 1) val =  FN39;
				else iA(0);
			 }
		  }
		  else iA(0);
	 }
	 else iA(0);
  }
  else if (dif(0) == -1) { /* Fine Nbr 40-55 */
	 if (dif(1) == 1) { /* Fine Nbr 40-43 */
		if (dif(2)  == 1) {
		  if ( !(pDif(0) == 0 && pDif(1) == 1 && pDif(2) == 1 )) val = -1;
		  else val =  FN40;
		}
		else if (dif(2) == 0){
		  if( !(pDif(0) == 0 && pDif(1) == 1 )) val = -1;
		  else {
			 if (pDif(2) == 0) val =  FN41;
			 else if (pDif(2) == 1) val =  FN42;
			 else iA(0);
		  }
		}
		else if (dif(2) == -1){
		  if ( !(pDif(0) == 0 && pDif(1) == 1 && pDif(2) == 0 )) val = -1;
		  else val =  FN43;
		}
		else iA(0);
	 }
	 else if (dif(1) == 0) { /* Fine Nbr 44-51 */
		if (dif(2) == 1) {
		  if( !(pDif(0) == 0 && pDif(2) == 1 )) val = -1;
		  else {
			 if (pDif(1) == 0) val =  FN44;
			 else if (pDif(1) == 1) val =  FN48;
			 else iA(0);
		  }
		}
		else if (dif(2) == 0){
		  if( !(pDif(0) == 0 )) { val = -1; }
		  else {
			 if (pDif(1) == 0){
				if (pDif(2) == 0) { val =  FN45; }
				else if (pDif(2) == 1) { val =  FN46; }
				else {iA(0); }
			 }
			 else if (pDif(1) == 1){
				if (pDif(2) == 0) val =  FN49;
				else if (pDif(2) == 1) val =  FN50;
				else iA(0);
			 }
		  }
		}
		else if (dif(2) == -1){
		  if( !(pDif(0) == 0 && pDif(2) == 0 )) val = -1;
		  else {
			 if (pDif(1) == 0) val =  FN47;
			 else if (pDif(1) == 1) val =  FN51;
			 else iA(0);
		  }
		}
	 }
	 else if (dif(1) == -1) { /* Fine Nbr 52-55 */
		if (dif(2) == 1) {
		  if ( !(pDif(0) == 0 && pDif(1) == 0 && pDif(2) == 1 )) val = -1;
		  else val =  FN52;
		}
		else if (dif(2) == 0){
		  if( !(pDif(0) == 0 && pDif(1) == 0 )) val = -1;
		  else {
			 if (pDif(2) == 0) val =  FN53;
			 else if (pDif(2) == 1) val =  FN54;
			 else iA(0);
		  }
		}
		else if (dif(2) == -1){
		  if ( !(pDif(0) == 0 && pDif(1) == 0 && pDif(2) == 0)) val = -1;
		  else val =  FN55;
		}
		else iA(0);
	 }
	 else iA(0);
  }
  else iA(0);

  return val;
}

int WNodeByDif(int pNrmNum, Index3 pDif){
  //cout << dif << std::endl;;
  int val = -1;
  for (int i = 0; i < 3; i++) iA ( ( pDif(i) == 0 || pDif(i) == 1) );
  if (pNrmNum == NRM0) {/* Wnode 0, 1, 6, 7, 36, 37 or 42 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN0; }
		  else { val = WN1; }
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN6; }
		  else { val = WN7; }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN36; }
		  else { val = WN37; }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN42; }
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN0 as Wnode" << std::endl;; exit(-1);
		  }
		}
	 }
  }
  else if (pNrmNum == NRM1) {/* Wnode 2, 3, 8, 9, 38 or 39 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN2; }
		  else { val = WN3; }
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN8; }
		  else { val = WN9; }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN38; }
		  else { val = WN39; }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN1 as Wnode" << std::endl;; exit(-1);
		  }
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN2 as Wnode" << std::endl;; exit(-1);
		  }
		}
	 }
  }
  else if (pNrmNum == NRM2) {/* Wnode 4, 5, 10, 11, 40, 41 or 43 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN4; }
		  else { val = WN5; }
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN10; }
		  else { val = WN11; }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN40; }
		  else { val = WN41; }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN3 as Wnode" << std::endl;; exit(-1);
		  }
		  else { val = WN43; }
		}
	 }
  }
  else if (pNrmNum == NRM3) {/* Wnode 12, 13, 18, 19, 44, 46 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN12; }
		  else { val = WN13; }
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN18; }
		  else { val = WN19; }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN44; }
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN4 as Wnode" << std::endl;; exit(-1);
		  }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN46; }
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN8 as Wnode" << std::endl;; exit(-1);
		  }
		}
	 }
  }
  else if (pNrmNum == NRM4) {/* Wnode 14, 15, 18, 20, 21 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN14; }
		  else { val = WN15; }
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN20; }
		  else { val = WN21; }
		}
	 } 
	 else { /* pDif(0) == 1 */
		val = -1;
		//std::cerr << "Wnodes Error:  Trying to val = FN5, FN6, FN9 or FN10 as Wnode" << std::endl;; exit(-1);
	 }
  }
  else if (pNrmNum == NRM5) {/* Wnode 16, 17, 22, 23, 45, 47 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN16; }
		  else { val = WN17; }
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN22; }
		  else { val = WN23; }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN7 as Wnode" << std::endl;; exit(-1);
		  }
		  else { val = WN45; }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN11 as Wnode" << std::endl;; exit(-1);
		  } 
		  else { val = WN47; }
		}
	 }
  }
  else if (pNrmNum == NRM6) {/* Wnode 24, 25, 30, 31, 48, 50, 51 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN24; }
		  else { val = WN25; }
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN30; }
		  else { val = WN31; }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN48; }
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN12 as Wnode" << std::endl;; exit(-1);
		  }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN50; }
		  else { val = WN51; }
		}
	 }
  }
  else if (pNrmNum == NRM7) {/* Wnode 26, 27, 32, 33, 52, 53 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN26; }
		  else { val = WN27; }
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN32; }
		  else { val = WN33; }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN13 as Wnode" << std::endl;; exit(-1);
		  }
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN14 as Wnode" << std::endl;; exit(-1);
		  }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN52; }
		  else { val = WN53; }
		}
	 }
  }
  else if (pNrmNum == NRM7) {/* Wnode 26, 27, 32, 33, 52, 53 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN26; }
		  else { val = WN27; }
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN32; }
		  else { val = WN33; }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){
			 val = -1;
			 std::cerr << "Wnodes Error:  Trying to val = FN13 as Wnode" << std::endl;; exit(-1);
		  }
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN14 as Wnode" << std::endl;; exit(-1);
		  }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN52; }
		  else { val = WN53; }
		}
	 }
  }
  else if (pNrmNum == NRM8) {/* Wnode 28, 29, 34, 35, 49, 54, 55 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN28; }
		  else { val = WN29; }
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN34; }
		  else { val = WN35; }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN15 as Wnode" << std::endl;; exit(-1);
		  }
		  else { val = WN49; }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN54; }
		  else { val = WN55; }
		}
	 }
  }
  else if (pNrmNum == NRM9) {/* Wnode 56, 57, 62, 76, 77, 82 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN56; }
		  else { val = WN57; }
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN62; }
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN16 as Wnode" << std::endl;; exit(-1);
		  }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN76; }
		  else { val = WN77; }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN82; }
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN28 as Wnode" << std::endl;; exit(-1);
		  }
		}
	 }
  }
  else if (pNrmNum == NRM10) {/* Wnode 58, 59, 78, 79 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN58; }
		  else { val = WN59; }
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN17 as Wnode" << std::endl;; exit(-1);
		  }
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN18 as Wnode" << std::endl;; exit(-1);
		  }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN78; }
		  else { val = WN79; }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN29 as Wnode" << std::endl;; exit(-1);
		  }
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN30 as Wnode" << std::endl;; exit(-1);
		  }
		}
	 }
  }
  else if (pNrmNum == NRM11) {/* Wnode 60, 61, 63, 80, 81, 83 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN60; }
		  else { val = WN61; }
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN19 as Wnode" << std::endl;; exit(-1);
		  }
		  else { val = WN63; }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN80; }
		  else { val = WN81; }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN31 as Wnode" << std::endl;; exit(-1);
		  } 
		  else { val = WN83; }
		}
	 }
  }
  else if (pNrmNum == NRM12) {/* Wnode 64, 66, 84, 86 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN64; }
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN20 as Wnode" << std::endl;; exit(-1);
		  }
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN66; }
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN22 as Wnode" << std::endl;; exit(-1);
		  }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN84; } 
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN32 as Wnode" << std::endl;; exit(-1);
		  }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN86; }
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN34 as Wnode" << std::endl;; exit(-1);
		  }
		}
	 }
  }
  else if (pNrmNum == NRM13) {
	 val = -1;
	 //std::cerr << "Wnodes Error:  Trying to val = something inside of B as a Wnode" << std::endl;; exit(-1);  
  }
  else if (pNrmNum == NRM14) {/* Wnode 65, 67, 85, 87 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN21 as Wnode" << std::endl;; exit(-1);  }
		  }
		  else { val = WN65; }
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN23 as Wnode" << std::endl;; exit(-1);
		  } 
		  else { val = WN67; }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN33 as Wnode" << std::endl;; exit(-1);
		  } 
		  else { val = WN85; }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN35 as Wnode" << std::endl;; exit(-1);
		  }
		  else { val = WN87; }
		}
	 }
  }
  else if (pNrmNum == NRM15) {/* Wnode 68, 70, 71, 88, 90, 91 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN68; }
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN24 as Wnode" << std::endl;; exit(-1);
		  } 
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN70; }
		  else { val = WN71; }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN88; }
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN36 as Wnode" << std::endl;; exit(-1);
		  } 
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN90; }
		  else { val = WN91; }
		}
	 }
  }
  else if (pNrmNum == NRM15) {/* Wnode 68, 70, 71, 88, 90, 91 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN68; }
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN24 as Wnode" << std::endl;; exit(-1);
		  } 
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN70; }
		  else { val = WN71; }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN88; }
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN36 as Wnode" << std::endl;; exit(-1);
		  } 
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN90; }
		  else { val = WN91; }
		}
	 }
  }
  else if (pNrmNum == NRM16) {/* Wnode 72, 73, 92, 93 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN25 as Wnode" << std::endl;; exit(-1);
		  } 
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN26 as Wnode" << std::endl;; exit(-1);
		  } 
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN72; } 
		  else { val = WN73; }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN37 as Wnode" << std::endl;; exit(-1);
		  } 
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN38 as Wnode" << std::endl;; exit(-1);
		  } 
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN92; }
		  else { val = WN93; }
		}
	 }
  }
  else if (pNrmNum == NRM17) {/* Wnode 69, 74, 75, 89, 94, 95 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN27 as Wnode" << std::endl;; exit(-1);
		  } 
		  else { val = WN69; }
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN74; } 
		  else { val = WN75; }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN39 as Wnode" << std::endl;; exit(-1);
		  } 
		  else { val = WN89; }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN94; }
		  else { val = WN95; }
		}
	 }
  }
  else if (pNrmNum == NRM18) {/* Wnode 116, 117, 122, 123, 96, 97, 102 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN96; }
		  else { val = WN97; }
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN102; }
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN40 as Wnode" << std::endl;; exit(-1);
		  } 
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN116; }
		  else { val = WN117; }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN122; }
		  else { val = WN123; }
		}
	 }
  }
  else if (pNrmNum == NRM19) {/* Wnode 87, 98, 118, 119, 124, 125 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN98; }
		  else { val = WN99; }
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN41 as Wnode" << std::endl;; exit(-1);
		  } 
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN42 as Wnode" << std::endl;; exit(-1);
		  } 
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN118; }
		  else { val = WN119; }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN124; }
		  else { val = WN125; }
		}
	 }
  }
  else if (pNrmNum == NRM20) {/* Wnode 100, 101, 103, 120, 121, 126, 127 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN100; }
		  else { val = WN101; }
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN43 as Wnode" << std::endl;; exit(-1);
		  } 
		  else { val = WN103; }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN120; }
		  else { val = WN121; }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN126; }
		  else { val = WN127; }
		}
	 }
  }
  else if (pNrmNum == NRM21) {/* Wnode 104, 106, 128, 129, 134, 135 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN104; }
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN44 as Wnode" << std::endl;; exit(-1);
		  } 
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN106; }
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN48 as Wnode" << std::endl;; exit(-1);
		  } 
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN128; }
		  else { val = WN129; }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN134; }
		  else { val = WN135; }
		}
	 }
  }
  else if (pNrmNum == NRM22) {/* Wnode 130, 131, 136, 137 */
	 if (pDif(0) == 0){
		val = -1;
		//std::cerr << "Wnodes Error:  Trying to val = FN45, FN46, FN49 or FN50 as Wnode" << std::endl;; exit(-1);  
	 } 	
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN130; }
		  else { val = WN131; }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN136; }
		  else { val = WN137; }
		}
	 }
  }
  else if (pNrmNum == NRM23) {/* Wnode 105, 107, 132, 133, 138, 139 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN47 as Wnode" << std::endl;; exit(-1);
		  }  
		  else { val = WN105; }
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN51 as Wnode" << std::endl;; exit(-1);
		  }  
		  else { val = WN107; }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN132; }
		  else { val = WN133; }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN138; }
		  else { val = WN139; }
		}
	 }
  }
  else if (pNrmNum == NRM24) {/* Wnode 108, 110, 111, 140, 141, 146, 147 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN108; }
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN52 as Wnode" << std::endl;; exit(-1);
		  }  
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN110; }
		  else { val = WN111; }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN140; }
		  else { val = WN141; }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN146; }
		  else { val = WN147; }
		}
	 }
  }
  else if (pNrmNum == NRM25) {/* Wnode 112, 113, 142, 143, 148, 149 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN53 as Wnode" << std::endl;; exit(-1);
		  }  
		  else {
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN54 as Wnode" << std::endl;; exit(-1);
		  }  
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN112; }
		  else { val = WN113; }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN142; }
		  else { val = WN143; }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN148; }
		  else { val = WN149; }
		}
	 }
  }
  else if (pNrmNum == NRM26) {/* Wnode 109, 114, 115, 144, 145, 150, 151 */
	 if (pDif(0) == 0){
		if (pDif(1) == 0){
		  if (pDif(2) == 0){
			 val = -1;
			 //std::cerr << "Wnodes Error:  Trying to val = FN55 as Wnode" << std::endl;; exit(-1);
		  }  
		  else { val = WN109; }
		} 
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN114; }
		  else { val = WN115; }
		}
	 } 
	 else { /* pDif(0) == 1 */
		if (pDif(1) == 0){
		  if (pDif(2) == 0){ val = WN144; }
		  else { val = WN145; }
		}
		else { /* pDif(1) == 1 */
		  if (pDif(2) == 0){ val = WN150; }
		  else { val = WN151; }
		}
	 }
  }
  else { std::cerr << "Error with Wnodes" << std::endl;; exit(-1); }

  return val;
}


//  LocalWords:  stkSrcTrgSymsPosNeg wxSymVals
