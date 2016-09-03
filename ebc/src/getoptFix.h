/*! \file getoptFix.h
	\brief to have appropriate include`s for getopt in OS X, Sun and other (Linux) systems  
*/
#ifndef __GETOPTFIX_H__
#define __GETOPTFIX_H__

/*
  Apple (OS X) and Sun systems declare getopt in unistd.h,
  other systems (Linux) use getopt.h
*/
#if defined ( __APPLE__ ) || ( defined (__SVR4) && defined (__sun) )
#include <unistd.h>
#else
#include "getopt.h"
#endif

#endif
