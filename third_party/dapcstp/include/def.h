/**
 * \file   def.h
 * \brief  definitions and macros
 *
 * \author Martin Luipersbeck 
 * \date   2015-05-03
 */

#ifndef DEF_H_
#define DEF_H_

#include <limits>
#include <string>
#include <sstream>
#include <iostream>

#define PROGRAM_NAME "da"
#define PROGRAM_VERSION "1.0"

#define RED      "\x1b[31m"
#define GREEN    "\x1b[32m"
#define NORMAL   "\x1b[0m"
#define CYAN     "\x1b[36m"
#define YELLOW   "\x1b[33m"
#define BLUE     "\x1b[34m"
#define GRAY     "\033[1m\033[30m"
#define YELLOWBI "\033[1m\033[93m"
#define GREENBI  "\033[1m\033[92m"
#define CYANBI   "\033[1m\033[96m"
#define PURPLEBI "\033[1m\033[95m"
#define BLUEBI   "\033[1m\033[94m"

#define WMAX std::numeric_limits<weight_t>::max()
#define DMAX std::numeric_limits<double>::max()
#define IMAX std::numeric_limits<int>::max()
#define gap(lb,ub) ((double)abs(ub-lb)/ub)
#define gapP(lb,ub) ((double)abs(ub-lb)*100.0/ub)
#define EXIT(...) {fprintf(stderr, __VA_ARGS__);exit(1);}

typedef int64_t weight_t;
typedef char flag_t;

#endif // DEF_H_
