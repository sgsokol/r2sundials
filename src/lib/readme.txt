# 2023-03-17 sokol@insa-toulouse.fr

# build internal lib for sundials
SUNTOP=/usr/local/src/cvodes-6.5.0
MYTOP=$HOME/dev/r/rcpp-pkgs/r2sundials
mkdir -p $MYTOP/src/lib/

# cmake build original CVODES to get sundials_config.h
cd $SUNTOP
mkdir instdir
mkdir build
cd build
ccmake ..
# c - as configure
# set:
# CMAKE_INSTALL_PREFIX = /usr/local/src/cvodes-6.5.0/instdir
# SUNDIALS_INDEX_SIZE = 32
# c - as configure
# g - as generate

make -j4

# copy sources
cp -a $SUNTOP/src/{cvodes,nvector,sundials,sunlinsol,sunmatrix,sunnonlinsol} $MYTOP/src/lib/
# remove fmod dirs
find $MYTOP/src/lib/ -type d -name fmod -exec rm -fr {} \;
find $MYTOP/src/lib/ -name CMakeLists.txt -exec rm -fr {} \;
# remove pthreads etc
rm -rf $MYTOP/src/lib/nvector/{cuda,manyvector,mpiplusx,openmp,openmpdev,parallel,parhyp,petsc,raja,trilinos,pthreads,sycl}
rm -rf $MYTOP/src/lib/sunnonlinsol/petscsnes
rm -rf $MYTOP/src/lib/sunlinsol/{klu,superludist,superlumt}
rm -rf $MYTOP/src/lib/sunmatrix/slunrloc
rm -rf $MYTOP/src/lib/sundials/sundials_xbraid.c

( cd $MYTOP/inst/include/nvector/ && rm -rf $(ls -1 | grep -v serial) )

# copy includes
cp -a $SUNTOP/src/sundials/*.h $MYTOP/inst/include/sundials/
cp -a $SUNTOP/include/{cvodes,nvector,sundials,sunlinsol,sunmatrix,sunnonlinsol} $MYTOP/inst/include
cp -a $SUNTOP/build/include/sundials/{sundials_config,sundials_export}.h $MYTOP/inst/include/sundials/
