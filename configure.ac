AC_INIT([r2sundials], 4.1.0-1, [sokol@insa-toulouse.fr])

if test -z "${R_HOME}"; then
   : ${R_HOME=`R RHOME`}
fi
if test -z "${R_HOME}"; then
   echo "could not determine R_HOME"
   exit 1
fi
CXX11=`"${R_HOME}/bin/R" CMD config CXX11`
CXX11STD=`"${R_HOME}/bin/R" CMD config CXX11STD`
CXX="${CXX11} ${CXX11STD}"
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXX11FLAGS`
AC_LANG(C++)

# config cvodes header path
AC_ARG_WITH([cvodes-include],
            AC_HELP_STRING([--with-cvodes-include=CVODES_INCLUDE_DIR],
                           [the location of CVODES header files, where cvodes/cvodes.h resides (default /usr/local/include)]),
            [cvodes_include=$withval])
CVODES_INCLUDE=""
if test [ -n "$cvodes_include" ] ; then
   CVODES_INCLUDE="${cvodes_include}"
else
   if test [ -n "${CVODES_INCLUDE}" ] ; then
      CVODES_INCLUDE="${CVODES_INCLUDE}"
   else
      CVODES_INCLUDE="/usr/local/include"
   fi
fi

dnl Check the headers can be found
SAVE_CPPFLAGS=$CPPFLAGS
CPPFLAGS="$CPPFLAGS -I$CVODES_INCLUDE/cvodes"
AC_CHECK_HEADERS([cvodes.h])
CPPFLAGS=$SAVE_CPPFLAGS
if test "${ac_cv_header_cvodes_h}" = no; then
   AC_MSG_ERROR("CVODES header cvodes.h not found")
fi

# config cvodes lib path
AC_ARG_WITH([cvodes-libs],
            AC_HELP_STRING([--with-cvodes-libs=CVODES_LIBS_DIR],
                           [the location of CVODES libraries, where libsundials_cvodes.so resides (default /usr/local/lib64)]),
            [cvodes_libs=$withval])
CVODES_LIBS=""
if test [ -n "$cvodes_libs" ] ; then
   CVODES_LIBS="${cvodes_libs}"
else
   if test [ -n "${CVODES_LIBS}" ] ; then
      CVODES_LIBS="${CVODES_LIBS}"
   else
      CVODES_LIBS="/usr/local/lib64"
   fi
fi

AC_SUBST(CVODES_INCLUDE)
AC_SUBST(CVODES_LIBS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
