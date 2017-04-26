/* myfloat.h
 *
 * Copyright (C) 2012 Ido Tal <idotal@ieee.org> and Alexander Vardy <avardy@ucsd.edu>
 *
 * This file is part of polarList.
 *
 */

#ifndef MYFLOAT
#define MYFLOAT

//void myfloat_atStart();
//void myfloat_atEnd();

//#define MYFLOAT_FLOAT
#define MYFLOAT_DOUBLE
//#define MYFLOAT_LONG_DOUBLE
//#define MYFLOAT_MPREAL

#ifdef MYFLOAT_FLOAT
#include <limits>
#include <math.h> //for log
typedef float myfloat;

//const myfloat plusInfinity = std::numeric_limits<myfloat>::infinity();
//const myfloat normalizedMin = std::numeric_limits<myfloat>::min();
#endif // MYFLOAT_FLOAT

#ifdef MYFLOAT_DOUBLE
#include <limits>
#include <math.h> //for log
typedef double myfloat;

//const myfloat plusInfinity = std::numeric_limits<myfloat>::infinity();
//const myfloat normalizedMin = std::numeric_limits<myfloat>::min();
#endif // MYFLOAT_DOUBLE

#ifdef MYFLOAT_LONG_DOUBLE
#include <limits>
#include <math.h>
#define log logl // have log take long double as an argument
#define sqrt sqrtl // have sqrt take long double as an argument
typedef long double myfloat;

//const myfloat plusInfinity = std::numeric_limits<myfloat>::infinity();
//const myfloat normalizedMin = std::numeric_limits<myfloat>::min();
#endif // MYFLOAT_LONG_DOUBLE

#ifdef MYFLOAT_MPREAL
#include "mpreal.h"
using namespace mpfr;
typedef mpreal myfloat;

//const myfloat plusInfinity = const_infinity();
//const myfloat normalizedMin = 0.0; // no denormalized numbers in MPFR

const int myfloat_default_prec = 256;
//const int myfloat_default_prec = 53;
#endif // MYFLOAT_MPREAL

#endif // MYFLOAT
