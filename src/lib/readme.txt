# 2023-03-17 sokol@insa-toulouse.fr

# build internal lib for sundials
SUNTOP=/usr/local/src/cvodes-7.2.1
MYTOP=$HOME/dev/r/rcpp-pkgs/r2sundials
mkdir -p $MYTOP/src/lib/
vold=$(fgrep Version: "$MYTOP"/DESCRIPTION | cut -f2 -d " ")

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

# rename old version to .old
( cd $MYTOP/src/lib/ && rename 's/(.+)$/$1-'$vold'.old/'  $MYTOP/src/lib/{cvodes,nvector,sundials,sunlinsol,sunmatrix,sunnonlinsol} )

# copy sources
cp -a $SUNTOP/src/{cvodes,nvector,sundials,sunlinsol,sunmatrix,sunnonlinsol} $MYTOP/src/lib/
# remove fmod dirs
find $MYTOP/src/lib/ -type d -name fmod_'*' -exec rm -fr {} \;
find $MYTOP/src/lib/ -name CMakeLists.txt -exec rm -fr {} \;
# remove pthreads etc
rm -rf $MYTOP/src/lib/nvector/{cuda,manyvector,mpiplusx,openmp,openmpdev,parallel,parhyp,petsc,raja,trilinos,pthreads,sycl}
rm -rf $MYTOP/src/lib/sunnonlinsol/petscsnes
rm -rf $MYTOP/src/lib/sunlinsol/{klu,superludist,superlumt}
rm -rf $MYTOP/src/lib/sunmatrix/slunrloc
rm -rf $MYTOP/src/lib/sundials/sundials_xbraid.c

( cd $MYTOP/inst/include/nvector/ && rm -rf $(ls -1 | grep -v serial) )

# copy includes
( cd $MYTOP/inst/include/ && rename 's/(.+)$/$1-'$vold'.old/' {cvodes,nvector,sundials,sunlinsol,sunmatrix,sunnonlinsol} )
#cp -a $SUNTOP/src/sundials/*.h $MYTOP/inst/include/sundials/
cp -a $SUNTOP/include/{cvodes,nvector,sundials,sunlinsol,sunmatrix,sunnonlinsol} $MYTOP/inst/include
cp -a $SUNTOP/build/include/sundials/{sundials_config,sundials_export}.h $MYTOP/inst/include/sundials/

# replace sunsrc var in src/Makevars with the output of
find lib -type f -name '*'.c ! -path '*.old' ! -name '*_mpi_*' | tr $'\n' ' '

# patch nvector_serial.c to exclude stdout
