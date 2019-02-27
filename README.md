# fmmvol3d
Sequential Fast Multipole Method (FMM) Volume Solver Code in Free-Space or with Periodic / Dirichlet Boundary Conditions in a Box

//-----------------------------------------------------

//-----------------------------------------------------

SEQUENTIAL CODE

//---------------------------------------
To compile the code

1. Declare FMM3D_DIR as the directory containing the code (e.g., fmmvol)

2. Make sure fftw3 (or fftw2.1.5), blas and lapack libraries are installed, and
change the library path variable LIBS in makefile.opt accordingly. If
not installed, they can be easily downloaded and installed.

3. From ${FMM3D_DIR} type:
        cd common
        make 

3. From ${FMM3D_DIR} type:
        cd 3d
        make vtt

//---------------------------------------
To run the test examples:

go to the 3d directory, and try the following:
(1) single layer laplace kernel, mid accuracy
        vtt 111.opt
(2) double layer stokes kernel, mid accuracy
        vtt 311.opt

        If the precomputed tables for a test are unavailable,
        the code will take a while to run the Adaptive Gaussian Quadrature
        for the first time, and results will be stored in
        ${FMM3D_DIR}/include for future use. To alter the precision of the
        tables (speeds up precomputation time), edit the dcuhre variables/

        Additional tests for various kernels as well as periodic and
        dirichlet tests are available.  Test cases can be viewed or
        augmented in 3d/exsol3d.hpp/cpp

        Currently, only precomputed tables are in the tarfile for the
        Laplacian.  For other tables, run the code as needed and wait for
        the tables to precompute.  If the tables are needed sooner, contact
        harperlangston@gmail.com

