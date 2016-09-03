
/*
 * declarations for making populations for runs
*/

#ifndef __MCT_MAKE_RUN_POPS__
#define __MCT_MAKE_RUN_POPS__

#include <mscplusplus/population_structure.hpp>



#include <boost/smart_ptr.hpp>



boost::shared_ptr< hudson_ms::PopulationStructure > makePopSingle();

boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandSmall(size_t Nzero);
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandBasic(size_t Nzero);
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandComplex(size_t Nzero);

boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandBalanced(size_t Nzero);
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandBalLessComm(size_t Nzero);
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandMSstats(size_t Nzero);

boost::shared_ptr< hudson_ms::PopulationStructure > makePopRaazRuns(
				size_t subNsamSizeEach,
				double alpha, double beta);

#endif
