#ifndef __SYMS_HPP__
#define __SYMS_HPP__

#include "common/nummat.hpp"
#include "common/vec3t.hpp"
#include "NbrDefs.hpp"

/*! UE = Upper Equivalent
	*  UC = Upper Check
	*  DE = Downward Equivalent
	*  DC = Downward Check
	*/
enum {	 UE=0,	 UC=1,	 DE=2,	 DC=3  };

int symType(int nbrTyp, int tabNbr);
int WXsymRef(int symNum);
int nbrSymRef(int type, int symNum);
int symRef(int type, int symNum);
int posNeg(Index3 &bdx);
double posNeg(int posneg, int coeff, int knl, int ds, int dt);
int basRef(Index3 &cIdx);
int symNbrGetRefs(int nbr, int type, int &bref, int &symNum, int &posneg, Index3 &idx, Index3 &cIdx);
int rvrsIdx(int i, int base);
Index3 pntRef(Index3 &idx, int i, int j, int k, int base);
int nbrSymMapSclVrt (int type, int nbr, double rad, Point3 &ctr);
int basSym(int bref, int c, int knl);
int srcTrgSym(int bref, int dts, int knl);
int WXsymType(int tabNde);
int WXsymNodeGetRefs(int WXnode, int &bref, int &symNum, int &posneg, Index3 &idx, Index3& cIdx);
int WXsymNodeGetRefsDwnChk(int WXnode, int &bref, int &symNum, int &posneg, Index3 &idx, Index3 &cIdx);
int nrmNbr(int WXnode);

int nrmNbrByDif(Index3 dif);
int fnNbrByDif(Index3 dif, Index3 pDif);
int WNodeByDif(int pNrmNum, Index3 pDif);

int npEffDatSze(int tp, int np, int sdof, int tdof);
int npPlnDatSze(int tp, int np, int sdof, int tdof);
int samPosCal(int np, double R, DblNumMat& positions, int type);
int regPosCal(int np, double R, DblNumMat& positions);
int grdSamPosCal(const bool regular, const int kval, const double R, DblNumMat& positions, const bool edges);

#define NEG -1
#define POS +1


/* For Stokes solver */
#define RXRX 0
#define RXRY 1
#define RXRZ 2
#define RYRY 3
#define RYRZ 4
#define RZRZ 5

/* Information for Basis Functions for Laplacian and Modified Laplacian */

#define B0 0 /* 1 */
#define B1 1 /* X */
#define B2 2 /* Y */
#define B3 3 /* Z */
#define B4 4 /* X^2 */
#define B5 5 /* X Y */
#define B6 6 /* X Z */
#define B7 7 /* Y^2 */
#define B8 8 /* Y Z */
#define B9 9 /* Z^2 */
#define B10 10 /* X^3 */
#define B11 11 /* X^2 Y */
#define B12 12 /* X^2 Z */
#define B13 13 /* X Y^2 */ 
#define B14 14 /* X Y Z */
#define B15 15 /* X Z^2 */
#define B16 16 /* Y^3 */
#define B17 17 /* Y^2 Z */
#define B18 18 /* Y Z^2 */
#define B19 19 /* Z^3 */

#define B20 20 /* X^4 */
#define B21 21 /* X^3 Y */
#define B22 22 /* X^3 Z	*/
#define B23 23 /* X^2 Y^2 */
#define B24 24 /* X^2 Y Z */
#define B25 25 /* X^2 Z^2 */
#define B26 26 /* X Y^3 */
#define B27 27 /* X Y^2 Z */
#define B28 28 /* X Y Z^2 */
#define B29 29 /* X Z^3 */ 
#define B30 30 /* Y^4 */
#define B31 31 /* Y^3 Z */
#define B32 32 /* Y^2 Z^2 */
#define B33 33 /* Y Z^3 */
#define B34 34 /* Z^4 */

#define B35 35 /* X^5 */
#define B36 36 /* X^4 Y */
#define B37 37 /* X^4 Z */
#define B38 38 /* X^3 Y^2 */
#define B39 39 /* X^3 Y Z */
#define B40 40 /* X^3 Z^2 */
#define B41 41 /* X^2 Y^3 */
#define B42 42 /* X^2 Y^2 Z */
#define B43 43 /* X^2 Y Z^2 */
#define B44 44 /* X^2 Z^3 */
#define B45 45 /* X Y^4 */
#define B46 46 /* X Y^3 Z */
#define B47 47 /* X Y^2 Z^2 */
#define B48 48 /* X Y Z^3 */
#define B49 49 /* X Z^4 */
#define B50 50 /* Y^5 */
#define B51 51 /* Y^4 Z */
#define B52 52 /* Y^3 Z^2 */
#define B53 53 /* Y^2 Z^3 */
#define B54 54 /* Y Z^4 */
#define B55 55 /* Z^5 */

/* Actually uses Chebyshev polys, so comments not correct but give sense of order */
#define B56 56 /* X^6 */
#define B57 57 /* X^5 Y */
#define B58 58 /* X^5 Z */
#define B59 59 /* X^4 Y^2 */
#define B60 60 /* X^4 Y Z */
#define B61 61 /* X^4 Z^2 */
#define B62 62 /* X^3 Y^3 */
#define B63 63 /* X^3 Y^2 Z */
#define B64 64 /* X^3 Y Z^2 */
#define B65 65 /* X^3 Z^3 */
#define B66 66 /* X^2 Y^4 */
#define B67 67 /* X^2 Y^3 Z */
#define B68 68 /* X^2 Y^2 Z^2 */
#define B69 69 /* X^2 Y Z^3 */
#define B70 70 /* X^2 Z^4 */
#define B71 71 /* X Y^5 */
#define B72 72 /* X Y^4 Z */
#define B73 73 /* X Y^3 Z^2 */
#define B74 74 /* X Y^2 Z^3 */
#define B75 75 /* X Y Z^4 */
#define B76 76 /* X Z^5 */
#define B77 77 /* Y^6 */
#define B78 78 /* Y^5 Z */
#define B79 79 /* Y^4 Z^2 */
#define B80 80 /* Y^3 Z^3 */
#define B81 81 /* Y^2  Z^4 */
#define B82 82 /* Y Z^5 */
#define B83 83 /* Z^6 */
/* Actually uses Chebyshev polys, so comments not correct but give sense of order */
#define B84 84 /* X^7 */
#define B85 85 /* X^6 Y */
#define B86 86 /* X^6 Z */
#define B87 87 /* X^5 Y^2 */
#define B88 88 /* X^5 Y Z */
#define B89 89 /* X^5 Z^2 */
#define B90 90 /* X^4 Y^3 */
#define B91 91 /* X^4 Y^2 Z */
#define B92 92 /* X^4 Y Z^2 */
#define B93 93 /* X^4 Z^3 */
#define B94 94 /* X^3 Y^4 */
#define B95 95 /* X^3 Y^3 Z */
#define B96 96 /* X^3 Y^2 Z^2 */
#define B97 97 /* X^3 Y Z^3 */
#define B98 98 /* X^3 Z^4 */
#define B99 99 /* X^2 Y^5 */
#define B100 100 /* X^2 Y^4 Z */
#define B101 101 /* X^2 Y^3 Z^2 */
#define B102 102 /* X^2 Y^2 Z^3 */
#define B103 103 /* X^2 Y Z^4 */
#define B104 104 /* X^2 Z^5 */
#define B105 105 /* X Y^6 */
#define B106 106 /* X Y^5 Z */
#define B107 107 /* X Y^4 Z^2 */
#define B108 108 /* X Y^3 Z^3 */
#define B109 109 /* X Y^2 Z^4 */
#define B110 110 /* X Y Z^5 */
#define B111 111 /* X Z^6 */
#define B112 112 /* Y^7 */
#define B113 113 /* Y^6 Z */
#define B114 114 /* Y^5 Z^2 */
#define B115 115 /* Y^4 Z^3 */
#define B116 116 /* Y^3  Z^4 */
#define B117 117 /* Y^2 Z^5 */
#define B118 118 /* Y Z^6 */
#define B119 119 /* Z^7 */

#endif
