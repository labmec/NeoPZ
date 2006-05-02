//
// C++ Interface: pzlog
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PZLOGH
#define PZLOGH

#include <string>
#ifdef LOG4CXX

#ifdef DEBUG
  #define DEBUG2
#endif

#include <sstream>
#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"
#include <log4cxx/propertyconfigurator.h>
using namespace log4cxx;
using namespace log4cxx::helpers;

#define LOGPZ_DEBUG(A,B) LOG4CXX_DEBUG(A,B)
#define LOGPZ_INFO(A,B) LOG4CXX_INFO(A,B)
#define LOGPZ_WARN(A,B) LOG4CXX_WARN(A,B)
#define LOGPZ_ERROR(A,B) LOG4CXX_ERROR(A,B)
#define LOGPZ_FATAL(A,B) LOG4CXX_FATAL(A,B)



#else

#define LOGPZ_DEBUG(A,B) {}
#define LOGPZ_INFO(A,B) {}
#define LOGPZ_WARN(A,B) {}
#define LOGPZ_ERROR(A,B) {}
#define LOGPZ_FATAL(A,B) {}


#endif

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


inline void InitializePZLOG(std::string &configfile)
{
#ifdef LOG4CXX
  log4cxx::PropertyConfigurator::configure(configfile);
  {
    log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("pz.mesh.tpzgeoelrefpattern"));
    logger->setAdditivity(false);
    logger->setLevel(log4cxx::Level::DEBUG);
  }
 {
    log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("pz.mesh.refpattern"));
    logger->setAdditivity(false);
    logger->setLevel(log4cxx::Level::DEBUG);
  }
#endif
}

void InitializePZLOG();

#endif
