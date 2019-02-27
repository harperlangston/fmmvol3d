#include "common/vecmatop.hpp"
#include "common/vec3t.hpp"
#include "exsol3d.hpp"

using std::cerr;

int Exsol3d::tdof(int qt)
{
  int dof = 0;
  if(       _et==KNL_LAP_S_U) {
	 switch(qt) {
	 case QNT_U: dof = 1; break;
	 case QNT_RHS: dof = 1; break;
	 case QNT_MAX_U: dof = 1; break;
	 case QNT_MAX_RHS: dof = 1; break;
	 }
  } else if(       _et==KNL_MODHEL_S_U) {
	 switch(qt) {
	 case QNT_U: dof = 1; break;
	 case QNT_RHS: dof = 1; break;
	 case QNT_MAX_U: dof = 1; break;
	 case QNT_MAX_RHS: dof = 1; break;
	 }
  } else if(_et==KNL_STK_S_U) {
	 switch(qt) {
	 case QNT_U: dof = 3; break;
	 case QNT_P: dof = 1; break;
	 case QNT_RHS: dof = 3; break;
	 case QNT_MAX_U: dof = 3; break;
	 case QNT_MAX_RHS: dof = 3; break;
	 }
  } else if(_et==KNL_NAV_S_U) {
	 switch(qt) {
	 case QNT_U: dof = 3; break;
	 case QNT_RHS: dof = 3; break;
	 }
  } else {
	 cerr<<"error"<<endl;
  }
  return dof;
}

// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "Exsol3d::quantity"
int Exsol3d::quantity(int qt, const DblNumMat& trgpos, DblNumVec& trgval)
{

  if (_ct != CHS_EMPTY) { iA(trgpos.n()*tdof(qt)==trgval.m()); }
  
  double L = 250.0;
  if(       _et==KNL_LAP_S_U) {
	 //----------------------------------------------
	 if(       _ct==CHS_LAP_FREE_SPACE) {
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double r2 = x*x + y*y + z*z;
			 trgval(i) = exp(-L*r2);
			 for (int a = -1; a < 2; a++) {
				for (int b = -1; b < 2; b++) {
				  for (int c = -1; c < 2; c++) {
					 double x1 = x +  a*(0.075);	
					 double y1 = y +  b*(0.075);
					 double z1 = z +  c*(0.075);
					 r2 = x1*x1 + y1*y1 + z1*z1;
					 trgval(i) += exp(-L*r2);
				  }
				}
			 }
			 //cerr << i << " " << trgval(i) << endl;
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double r2 = x*x + y*y + z*z;
			 trgval(i) = -exp(-r2*L)*((4.0*L*L*r2) - 6.0*L);
			 for (int a = -1; a < 2; a++) {
				for (int b = -1; b < 2; b++) {
				  for (int c = -1; c < 2; c++) {
					 double x1 = x +  a*(0.075);	
					 double y1 = y +  b*(0.075);
					 double z1 = z +  c*(0.075);
					 r2 = x1*x1 + y1*y1 + z1*z1;
					 trgval(i) += -exp(-r2*L)*((4.0*L*L*r2) - 6.0*L);
				  }
				}
			 }	
			 //trgval(i) = 1 + x + y + z + x*x + x*y + x*z + y*y + y*z + z*z + x*x*x + x*x*y + x*x*z + x*y*y + x*y*z + x*z*z + y*y*y + y*y*z + y*z*z + z*z*z;
			 //trgval(i) = 1.0 + x + y + z;
			 //if (abs(trgval(i)) < 10e-64) trgval(i) = 0.0;
		  }
		} else if(qt== QNT_MAX_U){
		  trgval(0) = exp(0.0);
		  for (int a = -1; a < 2; a++) {
			 for (int b = -1; b < 2; b++) {
				for (int c = -1; c < 2; c++) {
				  double x1 = a*(0.075);	
				  double y1 = b*(0.075);
				  double z1 = c*(0.075);
				  double r2 = x1*x1 + y1*y1 + z1*z1;
				  trgval(0) += exp(-L*r2);
				}
			 }
		  }
		  trgval(0) = abs(trgval(0));
		} else if(qt== QNT_MAX_RHS){
		  trgval(0) = -exp(0)*((0.0) - 6.0*L);
		  for (int a = -1; a < 2; a++) {
			 for (int b = -1; b < 2; b++) {
				for (int c = -1; c < 2; c++) {
				  double x1 = a*(0.075);	
				  double y1 = b*(0.075);
				  double z1 = c*(0.075);
				  double r2 = x1*x1 + y1*y1 + z1*z1;
				  trgval(0) += -exp(-r2*L)*((4.0*L*L*r2) - 6.0*L);
				}
			 }
		  }
		  trgval(0) = abs(trgval(0));
		}
	 } else if(       _ct==CHS_LAP_SPH) {
		double RSPH = coefs()[0];
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double r2 = x*x + y*y + z*z;
			 if (sqrt(r2) < RSPH) trgval(i) = 1.0/(6.0) * (RSPH*RSPH - r2) + RSPH*RSPH/(3.0);
			 else trgval(i) = RSPH*RSPH*RSPH/(3.0*sqrt(r2));
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double r2 = x*x + y*y + z*z;
			 if (sqrt(r2) < RSPH) trgval(i) = 1.0;
			 else trgval(i) = 0.0;
		  }
		} else if(qt== QNT_MAX_U){
		  trgval(0) = abs(1.0/(6.0) * (RSPH*RSPH) + RSPH*RSPH/(3.0));
		} else if(qt== QNT_MAX_RHS){
		  trgval(0) = 1.0;
		}
	 }
	 else if(       _ct==CHS_LAP_COL) {
		Point3 C1(-0.3125, -0.0625, 0.3125);
		Point3 C2(-0.0625, 0.3125, -0.3125);
		Point3 C3(0.3125, -0.3125, -0.0625);
		double RVAL = coefs()[0];
		double MVAL = coefs()[1];
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double al = 4.0*M_PI*MVAL;

			 double x1 = x-C1(0); double y1 = y-C1(1); double z1 = z-C1(2);
			 double r_1= sqrt(x1*x1 + y1*y1 + z1*z1);	 
			 x1 = x-C2(0); y1 = y-C2(1); z1 = z-C2(2);
			 double r_2= sqrt(x1*x1 + y1*y1 + z1*z1);
			 x1 = x-C3(0); y1 = y-C3(1); z1 = z-C3(2);
			 double r_3= sqrt(x1*x1 + y1*y1 + z1*z1);

			 double al2 = al*al;
			 double al3 = al*al2;
			 double al4 = al*al3;
			 double al5 = al*al4;
			 double al6 = al*al5;
			 double al7 = al*al6;
			 double tmp = 0.0;
			 double RS1 = r_1/RVAL;
			 double RS2 = r_2/RVAL;
			 double RS3 = r_3/RVAL;
			 double val = 0.0;
			 
			 if (RS1 < 1.0) {
				double RS = RS1;
				tmp = pow(RS,6.0)/84.0 - pow(RS,5.0)/30.0 + pow(RS,4.0)/40.0;
				tmp = tmp + 60.0/(al6) - 9.0/(al4) - 1.0/120.0 + 120.0/(al6*RS);
				
				tmp = tmp + (-120.0/(al6*RS) - 9.0/(al4) + 300.0/(al6))*cos(al*RS);
				tmp = tmp + (36.0*RS/(al4) + pow(RS,2.0)/(2.0*al2))*cos(al*RS);
				tmp = tmp + (-30.0*pow(RS,2.0)/(al4) - pow(RS,3.0)/(al2))*cos(al*RS);
				tmp = tmp + (pow(RS,4.0)/(2.0*al2))*(cos(al*RS));

				tmp = tmp + (12.0/(al5*RS) - 360.0/(al7*RS) - 96.0/(al5) + 120.0*RS/(al5))*(sin(al*RS));
				tmp = tmp + (-3.0*RS/(al3) + 8.0*pow(RS,2.0)/al3 - 5.0*pow(RS,3.0)/al3)*(sin(al*RS));
		
				//Add in RS2
				tmp = tmp + (-1.0/210.0 - 12.0/(al4) + 360.0/(al6))/(RS2);
				tmp = tmp + (-1.0/210.0 - 12.0/(al4) + 360.0/(al6))/(RS3);
				val = -tmp/RVAL;
			 }
			 else if (RS2 < 1.0) {
				double RS = RS2;
				tmp = pow(RS,6.0)/84.0 - pow(RS,5.0)/30.0 + pow(RS,4.0)/40.0;
				tmp = tmp + 60.0/(al6) - 9.0/(al4) - 1.0/120.0 + 120.0/(al6*RS);

				tmp = tmp + (-120.0/(al6*RS) - 9.0/(al4) + 300.00/(al6))*cos(al*RS);
				tmp = tmp + (36.0*RS/(al4) + pow(RS,2.0)/(2.0*al2))*cos(al*RS);
				tmp = tmp + (-30.0*pow(RS,2.0)/(al4) - pow(RS,3.0)/(al2))*cos(al*RS);
				tmp = tmp + (pow(RS,4.0)/(2.0*al2))*(cos(al*RS));

				tmp = tmp + (12.0/(al5*RS) - 360.0/(al7*RS) - 96.0/(al5) + 120.0*RS/(al5))*(sin(al*RS));
				tmp = tmp + (-3.0*RS/(al3) + 8.0*pow(RS,2.0)/al3 - 5.0*pow(RS,3.0)/al3)*(sin(al*RS));

				//Add in RS2
				tmp = tmp + (-1.0/210.0 - 12.0/(al4) + 360.0/(al6))/(RS1);
				tmp = tmp + (-1.0/210.0 - 12.0/(al4) + 360.0/(al6))/(RS3);
				val = -tmp/RVAL;
			 }	 
			 else if (RS3 < 1.0) {
				double RS = RS3;
				tmp = pow(RS,6.0)/84.0 - pow(RS,5.0)/30.0 + pow(RS,4.0)/40.0;
				tmp = tmp + 60.0/(al6) - 9.0/(al4) - 1.0/120.0 + 120.0/(al6*RS);
				
				tmp = tmp + (-120.0/(al6*RS) - 9.0/(al4) + 300.00/(al6))*cos(al*RS);
				tmp = tmp + (36.0*RS/(al4) + pow(RS,2.0)/(2.0*al2))*cos(al*RS);
				tmp = tmp + (-30.0*pow(RS,2.0)/(al4) - pow(RS,3.0)/(al2))*cos(al*RS);
				tmp = tmp + (pow(RS,4.0)/(2.0*al2))*(cos(al*RS));

				tmp = tmp + (12.0/(al5*RS) - 360.0/(al7*RS) - 96.0/(al5) + 120.0*RS/(al5))*(sin(al*RS));
				tmp = tmp + (-3.0*RS/(al3) + 8.0*pow(RS,2.0)/al3 - 5.0*pow(RS,3.0)/al3)*(sin(al*RS));

				//Add in RS2
				tmp = tmp + (-1.0/210.0 - 12.0/(al4) + 360.0/(al6))/(RS1);
				tmp = tmp + (-1.0/210.0 - 12.0/(al4) + 360.0/(al6))/(RS2);
		
				val = -tmp/RVAL;
			 }
			 else{
				tmp = (-1.0/210.0 - 12.0/(al4) + 360.0/(al6))/(RS1);
				tmp = tmp + (-1.0/210.0 - 12.0/(al4) + 360.0/(al6))/(RS2);
				tmp = tmp + (-1.0/210.0 - 12.0/(al4) + 360.0/(al6))/(RS3);
				val = -tmp/RVAL;
			 }
			 trgval(i) = val;
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double al = 4.0*M_PI*MVAL;

			 double x1 = x-C1(0); double y1 = y-C1(1); double z1 = z-C1(2);
			 double r2 = (x1*x1 + y1*y1 + z1*z1);
			 double r_1= sqrt(r2)/RVAL;
			 	 
			 x1 = x-C2(0); y1 = y-C2(1); z1 = z-C2(2);
			 r2 = (x1*x1 + y1*y1 + z1*z1);
			 double r_2= sqrt(r2)/RVAL;
			 
			 x1 = x-C3(0); y1 = y-C3(1); z1 = z-C3(2);
			 r2 = (x1*x1 + y1*y1 + z1*z1);
			 double r_3= sqrt(r2)/RVAL;

			 double val = 0.0;
			 if (r_1 < 1.0) {
				double r_12 = r_1*r_1;
				val = ((r_1 - r_12)*sin(al*r_1/2.0));
				val = val*val;
			 }
			 else if (r_2 < 1.0) {
				double r_22 = r_2*r_2;
				val = ((r_2 - r_22)*sin(al*r_2/2.0));
				val = val*val;
			 }
			 else if (r_3 < 1.0) {
				double r_32 = r_3*r_3;
				val = ((r_3 - r_32)*sin(al*r_3/2.0));
				val = val*val;
			 }
			 else val = 0.0;
			 val = val/(RVAL*RVAL*RVAL);
			 if (val != 0.0){
				//cerr << val << endl; 
			 }
			 trgval(i) = val;
		  }
		} else if(qt== QNT_MAX_U){
		  trgval(0) = abs((-1.0/120.0 - (6.0/(pow(4.0*M_PI*MVAL,4.0))))/RVAL + (-1.0/105.0 - 24.0/(pow(4.0*M_PI*MVAL, 4.0)) + 720.0/(pow(4.0*M_PI*MVAL, 6.0)))/(sqrt(0.59375)));
		} else if(qt== QNT_MAX_RHS){
		  trgval(0) = 1.0;
		}
	 } else if(       _ct==CHS_LAP_PER) {
		//------------------
		double CON = coefs()[0];
		double PK = coefs()[1];
		double adjval  = 0.0;
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 trgval(i) = -PK*sin(M_PI*x*CON)*sin(M_PI*y*CON)*sin(M_PI*z*CON);
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 trgval(i) = -PK*(3.0*(CON*CON)*M_PI*M_PI + adjval)*sin(M_PI*x*CON)*sin(M_PI*y*CON)*sin(M_PI*z*CON);
		  }
		  //#warning need to change max vals
		}  else if(qt== QNT_MAX_U){
		  trgval(0) = 1.0;
		} else if(qt== QNT_MAX_RHS){
		  trgval(0) = 1.0;
		}
	 } else if(       _ct==CHS_LAP_DIR) {
		//------------------
		double CON = coefs()[0];
		double PK = coefs()[1];
		double adjval  = 0.0;
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 trgval(i) = -PK*sin(CON*M_PI*(1.0+x))*sin(CON*M_PI*(1.0+y))*sin(CON*M_PI*(1.0+z));
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 trgval(i) = -PK*(3.0*(CON*CON)*M_PI*M_PI)*sin(CON*M_PI*(1.0+x))*sin(CON*M_PI*(1.0+y))*sin(CON*M_PI*(1.0+z));
		  }
		  //#warning need to change max vals
		}  else if(qt== QNT_MAX_U){
		  trgval(0) = 1.0;
		} else if(qt== QNT_MAX_RHS){
		  trgval(0) = 1.0;
		}
	 } else if(       _ct==CHS_LAP_DIR_UNB) {
		//------------------
		double CON = coefs()[0];
		double PK = coefs()[1];
		double L = sqrt(2.0); double L2 = L*L;
		double adjval  = 0.0;
		double cx = 0.75; double cy = 0.75; double cz = 0.75;
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double r2 = (x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz);
			 trgval(i) = PK*sin(M_PI*CON*(1.0+x))*sin(M_PI*CON*(1.0+y))*sin(M_PI*CON*(1.0+z))*exp(r2*L);
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double r2 = (x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz);
			 trgval(i) = 3*CON*CON*M_PI*M_PI*PK*exp(L*r2)*sin(CON*M_PI*(x + 1))*sin(CON*M_PI*(y + 1))*sin(CON*M_PI*(z + 1)) - 6*L*PK*exp(L*r2)*sin(CON*M_PI*(x + 1))*sin(CON*M_PI*(y + 1))*sin(CON*M_PI*(z + 1)) - L2*PK*exp(L*r2)*sin(CON*M_PI*(x + 1))*sin(CON*M_PI*(y + 1))*sin(CON*M_PI*(z + 1))*(2*cx - 2*x)*(2*cx - 2*x) - L2*PK*exp(L*r2)*sin(CON*M_PI*(x + 1))*sin(CON*M_PI*(y + 1))*sin(CON*M_PI*(z + 1))*(2*cy - 2*y)*(2*cy - 2*y) - L2*PK*exp(L*r2)*sin(CON*M_PI*(x + 1))*sin(CON*M_PI*(y + 1))*sin(CON*M_PI*(z + 1))*(2*cz - 2*z)*(2*cz - 2*z) + 2*L*CON*M_PI*PK*exp(L*r2)*cos(CON*M_PI*(x + 1))*sin(CON*M_PI*(y + 1))*sin(CON*M_PI*(z + 1))*(2*cx - 2*x) + 2*L*CON*M_PI*PK*exp(L*r2)*cos(CON*M_PI*(y + 1))*sin(CON*M_PI*(x + 1))*sin(CON*M_PI*(z + 1))*(2*cy - 2*y) + 2*L*CON*M_PI*PK*exp(L*r2)*cos(CON*M_PI*(z + 1))*sin(CON*M_PI*(x + 1))*sin(CON*M_PI*(y + 1))*(2*cz - 2*z);
		  }
		  //#warning need to change max vals
		}  else if(qt== QNT_MAX_U){
		  trgval(0) = 1.0;
		} else if(qt== QNT_MAX_RHS){
		  trgval(0) = 1.0;
		}
   } else if(       _ct==CHS_LAP_POLY_TEST) {
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 trgval(i) = 0.0;
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 trgval(i) = 1 + x + y + z + x*x + x*y + x*z + y*y + y*z + z*z + x*x*x + x*x*y + x*x*z + x*y*y + x*y*z + x*z*z + y*y*y + y*y*z + y*z*z + z*z*z + x*x*x*x;
		  }
		}
	 } else if ( _ct == CHS_EMPTY){
		//Do nothing
	 }
	 else {
		//------------------
		iA(0);
	 }
  }else if(       _et==KNL_MODHEL_S_U) {
	 double lambda = coefs()[1];
	 double alpha = lambda*lambda;
	 //----------------------------------------------
	 if(       _ct==CHS_MODHEL_FREE_SPACE) {
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double r2 = x*x + y*y + z*z;
			 trgval(i) = exp(-L*r2);
			 for (int a = -1; a < 2; a++) {
				for (int b = -1; b < 2; b++) {
				  for (int c = -1; c < 2; c++) {
					 double x1 = x +  a*(0.075);	
					 double y1 = y +  b*(0.075);
					 double z1 = z +  c*(0.075);
					 r2 = x1*x1 + y1*y1 + z1*z1;
					 trgval(i) += exp(-L*r2);
				  }
				}
			 }
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double r2 = x*x + y*y + z*z;
			 trgval(i) = -exp(-r2*L)*((4.0*L*L*r2) - 6.0*L - alpha);
			 for (int a = -1; a < 2; a++) {
				for (int b = -1; b < 2; b++) {
				  for (int c = -1; c < 2; c++) {
					 double x1 = x +  a*(0.075);	
					 double y1 = y +  b*(0.075);
					 double z1 = z +  c*(0.075);
					 r2 = x1*x1 + y1*y1 + z1*z1;
					 trgval(i) += -exp(-r2*L)*((4.0*L*L*r2) - 6.0*L - alpha);
				  }
				}
			 }
			 //trgval(i) = 1 + x + y + z + x*x + x*y + x*z + y*y + y*z + z*z + x*x*x + x*x*y + x*x*z + x*y*y + x*y*z + x*z*z + y*y*y + y*y*z + y*z*z + z*z*z + x*x*x*x + x*x*x*x*x + x*x*x*x*x*x + x*x*x*x*x*x*x + x*y*z*x*y*z*x*y*z + sin(7.0*x);
			 //trgval(i) = 1.0;
			 //if (abs(trgval(i)) < 10e-64) trgval(i) = 0.0;
		  }
		} else if(qt== QNT_MAX_U){
		  trgval(0) = exp(0.0);
		  for (int a = -1; a < 2; a++) {
			 for (int b = -1; b < 2; b++) {
				for (int c = -1; c < 2; c++) {
				  double x1 = a*(0.075);	
				  double y1 = b*(0.075);
				  double z1 = c*(0.075);
				  double r2 = x1*x1 + y1*y1 + z1*z1;
				  trgval(0) += exp(-L*r2);
				}
			 }
		  }
		  trgval(0) = abs(trgval(0));
		} else if(qt== QNT_MAX_RHS){
		  trgval(0) = -exp(0)*((0.0) - 6.0*L - alpha);
		  for (int a = -1; a < 2; a++) {
			 for (int b = -1; b < 2; b++) {
				for (int c = -1; c < 2; c++) {
				  double x1 = a*(0.075);	
				  double y1 = b*(0.075);
				  double z1 = c*(0.075);
				  double r2 = x1*x1 + y1*y1 + z1*z1;
				  trgval(0) += -exp(-r2*L)*((4.0*L*L*r2) - 6.0*L - alpha);
				}
			 }
		  }
		  trgval(0) = abs(trgval(0));
		}
	 }
	 else {
		iA(0);
	 }
  } else if(_et==KNL_STK_S_U) {
	 //----------------------------------------------
	 if(       _ct==CHS_STK_FREE_SPACE) {
		double l = L/2.0;
		double l2 = l*l;
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double r2 = x*x + y*y + z*z;
			 int trgdof = tdof(qt);
			 for (int d = 0; d < trgdof; d++){
				if (d == 0){      trgval(i*trgdof + d) = 2.0*l*exp(-l*r2)*(z - y); }
				else if (d == 1){ trgval(i*trgdof + d) = 2.0*l*exp(-l*r2)*(x - z); }
				else if (d == 2){ trgval(i*trgdof + d) = 2.0*l*exp(-l*r2)*(y - x); }
			 }
			 
			 for (int a = -1; a < 2; a++) {
				for (int b = -1; b < 2; b++) {
				  for (int c = -1; c < 2; c++) {
					 double x1 = x +  a*(0.075);	
					 double y1 = y +  b*(0.075);
					 double z1 = z +  c*(0.075);
					 r2 = x1*x1 + y1*y1 + z1*z1;
					 for (int d = 0; d < trgdof; d++){
						if (d == 0){      trgval(i*trgdof + d) += 2.0*l*exp(-l*r2)*(z1 - y1); }
						else if (d == 1){ trgval(i*trgdof + d) += 2.0*l*exp(-l*r2)*(x1 - z1); }
						else if (d == 2){ trgval(i*trgdof + d) += 2.0*l*exp(-l*r2)*(y1 - x1); }
					 }
				  }
				}
			 }
			 
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double r2 = x*x + y*y + z*z;
			 int trgdof = tdof(qt);
			 for (int d = 0; d < trgdof; d++){
				if (d == 0){
				  trgval(i*trgdof + d) = -4.0*l2*exp(-l*r2)*(z - y)*(2.0*l*r2 - 5.0);
				  //trgval(i*trgdof + d) = -(z - y);
				}
				else if (d == 1){
				  trgval(i*trgdof + d) = -4.0*l2*exp(-l*r2)*(x - z)*(2.0*l*r2 - 5.0);
				  //trgval(i*trgdof + d) = -(x - z);
				}
				else if (d == 2){
				  trgval(i*trgdof + d) = -4.0*l2*exp(-l*r2)*(y - x)*(2.0*l*r2 - 5.0);
				  //trgval(i*trgdof + d) = -(y - x);
				}
			 }
			 for (int a = -1; a < 2; a++) {
				for (int b = -1; b < 2; b++) {
				  for (int c = -1; c < 2; c++) {
					 double x1 = x +  a*(0.075);	
					 double y1 = y +  b*(0.075);
					 double z1 = z +  c*(0.075);
					 r2 = x1*x1 + y1*y1 + z1*z1;
					 for (int d = 0; d < trgdof; d++){
						if (d == 0){		
						  trgval(i*trgdof + d) += -4.0*l2*exp(-l*r2)*(z1 - y1)*(2.0*l*r2 - 5.0);
						}
						else if (d == 1){
						  trgval(i*trgdof + d) += -4.0*l2*exp(-l*r2)*(x1 - z1)*(2.0*l*r2 - 5.0);
						}
						else if (d == 2){
						  trgval(i*trgdof + d) += -4.0*l2*exp(-l*r2)*(y1 - x1)*(2.0*l*r2 - 5.0);
						}
					 }
				  }
				}
			 }
		  }
		} else if(qt== QNT_MAX_U){
		  //#warning need to fix
		  trgval(0) = 1.0;
		} else if(qt== QNT_MAX_RHS){
		  //#warning need to fix
		  trgval(0) = 1.0;
		} else {
		  //------------------
		  iA(0);
		}
	 }
  } else if(_et==KNL_NAV_S_U) {
	 //----------------------------------------------
	 //------------------
		iA(0);
  }

  return (0);
}
