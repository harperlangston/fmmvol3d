all: libcommon

include ../makefile.opt

LIB_SRC = 	vecmatop.cpp svdrep.cpp memAlloc.cpp Dcuhre.cpp CmptTbls.cpp kernel3d.cpp syms.cpp
LIB_FSRC = f_dcuhre.f f_qshep3d.f
LIB_CSRC = ffunc.c

LIB_OBJ = 	$(LIB_FSRC:.f=.o) $(LIB_SRC:.cpp=.o) $(LIB_CSRC:.c=.o)

DEP     = 	$(LIB_SRC:.cpp=.d) $(LIB_FSRC:.f=.d) $(LIB_CSRC:.c=.d)

libcommon: 	$(LIB_OBJ)
	$(AR) $(ARFLAGS) libcommon.a $(LIB_OBJ) 
	$(RANLIB) libcommon.a

-include $(DEP)

