
################################################################################
## Configuring the pkg-config .pc file.                                       ##
################################################################################

# Get include directories for dependencies.
get_filename_component(TBB_INCLUDE_DIR "-L${TBB}" DIRECTORY)

set(PKG_CONFIG_LIBS "-lflint ${TBB_INCLUDE_DIR} -ltbb ${PKG_CONFIG_MEMALLOC}")
set(PKG_CONFIG_CFLAGS "-std=c++23 -fopenmp -Ofast -march=native -flto -funroll-loops -fprefetch-loop-arrays -fomit-frame-pointer -fno-rtti")

configure_file(ATEAMS++.pc.in ATEAMS++.pc @ONLY)
