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
#include "Vfmm3d.hpp"
#include "common/vecmatop.hpp"

using std::cerr;
using std::endl;

template <class VF>
int VFMM3d<VF>::vcheck(double& rerr, double& inferr) {

  int trgDOF = this->trgDOF();
  vector<int> ordVec; iC( vlet()->dwnOrderCollect(ordVec) ); //BOTTOM UP
  DblNumVec diffs(trgDOF); DblNumVec sums(trgDOF); DblNumVec infs(trgDOF); DblNumVec exainf(trgDOF); DblNumVec L2diff(trgDOF); DblNumVec L2ediff(trgDOF);
  for (int d = 0; d < trgDOF; d++) {
    diffs(d) = 0.0; sums(d) = 0.0; infs(d) = 0.0; exainf(d) = 0.0; L2diff(d) = 0.0; L2ediff(d) = 0.0;
  }
  DblNumMat ptmp(3,1);
  if ((this->exsol3d())->ct() == CHS_LAP_COL){ (this->exsol3d())->quantity(QNT_MAX_U, ptmp, exainf);
    cerr << "COL = " << exainf(0) << endl;
  }

  int nonznodes = 0;

  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];

	 VF& curNode = this->node(gNodeIdx);
	 curNode.errcheck() = false;
	 bool chkhere = false;
	 //No real reason to separate out periodic/dirichlet vs. other tests except for debugging purposes if needed
	 if (vlet()->periodic()){
		PerNode& curPerNode = vlet()->pernode(gNodeIdx);
		//cerr << gNodeIdx << " " << curPerNode.bdryUnodes().size() << " " << curPerNode.bdryVnodes().size() << " " << curPerNode.bdryWnodes().size() << " " << curPerNode.bdryXnodes().size() << endl;
		  if (vlet()->terminal(gNodeIdx) == true){
		  chkhere = true;
		}
	 }
	 else {
		if ( vlet()->terminal(gNodeIdx)==true ) {
		  chkhere = true;
		}
	 }
	 if (chkhere){
	   if(vlet()->tag(gNodeIdx) & LET_TRGNODE && vlet()->tag(gNodeIdx) && LET_SRCNODE) { //evaluator
	     DblNumMat srcPos(vlet()->grdSrcExaPos(gNodeIdx));
	     DblNumVec srcDen(vlet()->grdSrcSamPos().n() * (this->knl()).srcDOF());
	     (this->exsol3d())->quantity(QNT_RHS, srcPos, srcDen);
	     if (srcDen.linfty() > 0.0) {
	       nonznodes++;
		  DblNumVec grdExaValgNodeIdx(grdExaVal(gNodeIdx));
		  DblNumMat grdExaPosgNodeIdx(vlet()->grdTrgExaPos(gNodeIdx));
		  DblNumVec grdRealValgNodeIdx(grdExaValgNodeIdx.m());
		  DblNumVec grdSrcDengNodeIdx(grdExaValgNodeIdx.m());

		  (this->exsol3d())->quantity(QNT_U, grdExaPosgNodeIdx, grdRealValgNodeIdx);
		  (this->exsol3d())->quantity(QNT_RHS, grdExaPosgNodeIdx, grdSrcDengNodeIdx);
		  double h = vlet()->radius(gNodeIdx)*2.0/vlet()->kTrgVal();
		  
		  curNode.errcheck() = true;
		  for (int j = 0; j < grdExaPosgNodeIdx.n(); j++){
			 for (int d = 0; d < trgDOF; d++){
				double rval = grdRealValgNodeIdx(j*trgDOF + d);
				double sval = grdSrcDengNodeIdx(j*trgDOF + d);
				if (fabs(rval) > 10e-12){
				  double rval = grdRealValgNodeIdx(j*trgDOF + d);
				  double eval = grdExaValgNodeIdx(j*trgDOF + d);
				  double tmp = abs(rval - eval);
				  //cerr << Point3(grdExaPosgNodeIdx(0,j),grdExaPosgNodeIdx(1,j),grdExaPosgNodeIdx(2,j)) << " " << grdRealValgNodeIdx(j*trgDOF + d) << " " <<  grdExaValgNodeIdx(j*trgDOF + d) << " " << tmp << " " << grdRealValgNodeIdx(j*trgDOF + d)/grdExaValgNodeIdx(j*trgDOF + d) << endl;
				  infs(d) = max(infs(d), tmp);
				  
				  exainf(d) = max(exainf(d), rval);

				  diffs(d) += tmp*tmp;
				  L2diff(d) += tmp*tmp*(h*h*h);
				  L2ediff(d) += (tmp/exainf(d))*(tmp/exainf(d))*(pow(h,3.0));
				  sums(d) += rval*rval;
				}	
			 }
		  }
	     }
		}	
	 }
  }
  std::cout << "Num points = " << vlet()->trmNodeCnt()*vlet()->trgGrdSze() << " " << nonznodes*vlet()->trgGrdSze() << endl;
  for (int d = 0; d < trgDOF; d++){
	 
	 diffs(d) = sqrt(diffs(d));
	 sums(d) = sqrt(sums(d));
	 rerr = diffs(d)/sums(d);
	 double L2err = sqrt(L2diff(d));
	 double L2Cerr = sqrt(L2ediff(d));
	 cerr << exainf(d) << " " <<diffs(d) << " " << sqrt(diffs(d)) << " " << infs(d) << " " << rerr << " " << sqrt(L2diff(d)) << " " << sqrt(L2diff(d))/exainf(d) << " " << L2err << " " << L2Cerr << " " << infs(d) << " " << infs(d)/exainf(d) << endl;
	 cerr << "EPS = " << vlet()->eps_rhs() << endl;;
  }
  return(0);
  }
