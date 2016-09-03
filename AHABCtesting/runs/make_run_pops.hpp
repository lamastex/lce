
/*
 * declarations for making populations for runs for Sharyn
*/

#ifndef __MAKE_S_RUN_POPS__
#define __MAKE_S_RUN_POPS__

#include <mscplusplus/population_structure.hpp>



#include <boost/smart_ptr.hpp>



boost::shared_ptr< hudson_ms::PopulationStructure > makePopScenario1(size_t nZero);

boost::shared_ptr< hudson_ms::PopulationStructure > makePopScenario2(size_t nZero);

boost::shared_ptr< hudson_ms::PopulationStructure > makePopScenario3(size_t nZero);

boost::shared_ptr< hudson_ms::PopulationStructure > makePopScenario4(size_t nZero);

boost::shared_ptr< hudson_ms::PopulationStructure > makePopScenario4alt(size_t nZero);

boost::shared_ptr< hudson_ms::PopulationStructure > makePopScenario4contrast(size_t nZero);

boost::shared_ptr< hudson_ms::PopulationStructure > makePopScenario5(size_t nZero);

#endif


