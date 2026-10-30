
################################################################################
## Configuring the pkg-config .pc file.                                       ##
################################################################################

# Get include directories for dependencies.
get_filename_component(TBB_INCLUDE_DIR "-L${TBB}" DIRECTORY)

# Check whether we're using OpenMP.
if(OpenMP)
	set(OPENMP_FLAG "-fopenmp")
else()
	set(OPENMP_FLAG "")
endif()

set(PKG_CONFIG_LIBS "-lflint ${TBB_INCLUDE_DIR} -ltbb ${PKG_CONFIG_MEMALLOC}")
set(PKG_CONFIG_CFLAGS "-std=c++23 ${OPENMP_FLAG} -Ofast -march=native -flto -funroll-loops -fprefetch-loop-arrays -fomit-frame-pointer -fno-rtti")

configure_file(ATEAMS++.pc.in ATEAMS++.pc @ONLY)
