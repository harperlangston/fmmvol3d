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

#ifndef _EXSOL3D_HPP_
#define _EXSOL3D_HPP_

#include "common/kernel3d.hpp"

using std::vector;

enum {
  QNT_U = 0,
  QNT_P = 1,
  QNT_RHS = 2,
  QNT_MAX_U = 3,
  QNT_MAX_RHS = 4,
  QNT_ERR = -1
};

enum {
  CHS_EMPTY = 0,
  CHS_LAP_POLY_TEST = 1,
  CHS_LAP_FREE_SPACE = 10,
  CHS_LAP_SPH = 11,
  CHS_LAP_COL = 12,
  CHS_LAP_PER = 13,
  CHS_LAP_DIR = 14,
  CHS_LAP_DIR_UNB = 15,
  CHS_MODHEL_FREE_SPACE = 20,
  CHS_STK_FREE_SPACE = 30,
  
  CHS_NAV_CST = 51,
  CHS_NAV_NAV = 52,
  CHS_NAV_ROT = 53,
  CHS_NAV_PB3 = 54,
  //error
  CHS_ERR = -1
};

class Exsol3d
{
protected:
  int _et; //equation type
  vector<double> _coefs; //coefs
  int _ct; //choice type  //int _qt; //quantity type
  
public:
  Exsol3d() {;}
  Exsol3d(int et, const vector<double>& coefs, int ct): _et(et), _coefs(coefs), _ct(ct) {;}
  int& et() { return _et; }
  vector<double>& coefs() { return _coefs; }
  int& ct() { return _ct; }
  int tdof(int qt);
  int quantity(int qt, const DblNumMat& trgpos, DblNumVec& trgqnt);
};

#endif
