all: libfmm3d

include ../makefile.opt

LIB_SRC =   fmm3d_impl.cpp \
				matmgnt3d.cpp \
				fmm3d.cpp \
				fmm3d_setup.cpp \
				fmm3d_eval.cpp \
				fmm3d_check.cpp \
				dense3d.cpp \
			   exsol3d.cpp \
				let3d.cpp \
				Vlet3d.cpp \
				Vfmm3d_setup.cpp \
				Vfmm3d.cpp \
				Vfmm3d_eval.cpp \
				Vfmm3d_check.cpp \
				Vfmm3d_Coeffs.cpp \
				Vfmm3d_FMM.cpp \
				Efmm3d.cpp



LIB_CSRC = 

NBRS_SRC = 

LIB_OBJ = 	$(LIB_SRC:.cpp=.o)
NBRS_OBJ =  $(NBRS_SRC:.cpp=.o)
C_OBJ = $(LIB_CSRC:.c=.o)

DEP     = 	$(LIB_SRC:.cpp=.d) $(LIB_CSRC:.c=.d) $(TST_SRC:.cpp=.d)

libfadgs3d: 	$(LIB_OBJ)
	$(AR) $(ARFLAGS) libfadgs3d.a $(LIB_OBJ)
	$(RANLIB) libfadgs3d.a

vtt: libfadgs3d $(C_OBJ) vtt.o
	${CXX} -o vtt vtt.o libfadgs3d.a ../common/libcommon.a  ${LDFLAGS}

vftt: libfadgs3d $(C_OBJ) vftt.o
	${CXX} -o vftt vftt.o libfadgs3d.a ../common/libcommon.a  ${LDFLAGS}


-include $(DEP)

#	${MEX} fmm3d_mex.cpp ${INCLUDES} libfmm3d.a ../common/libcommon.a librfftw.a libfftw.a liblapack.a libblas.a



# DO NOT DELETE
