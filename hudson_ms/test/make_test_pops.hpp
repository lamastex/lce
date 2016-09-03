
/*
 * declarations for making test populations for testing
*/

#ifndef __MAKE_TEST_POPS__
#define __MAKE_TEST_POPS__

#include <mscplusplus/population_structure.hpp>



#include <boost/smart_ptr.hpp>




boost::shared_ptr< hudson_ms::PopulationStructure > makeTestPop1();
boost::shared_ptr< hudson_ms::PopulationStructure > makeTestPop2();
boost::shared_ptr< hudson_ms::PopulationStructure > makePopSingle();
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandSmall(size_t Nzero);
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandBasic(size_t Nzero);
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandComplex(size_t Nzero);
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandBad1(size_t Nzero);
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandBad2(size_t Nzero);
boost::shared_ptr< hudson_ms::PopulationStructure > makePopInstantSizeChange();
boost::shared_ptr< hudson_ms::PopulationStructure > makePopOutgroupSequence();
boost::shared_ptr< hudson_ms::PopulationStructure > makePopInstantSizeChangeAndExpGrowth();
boost::shared_ptr< hudson_ms::PopulationStructure > makePopTwoSpeciesWithPopSizeDiff();
boost::shared_ptr< hudson_ms::PopulationStructure > makePopBackwardsSplit();
boost::shared_ptr< hudson_ms::PopulationStructure > makePopSteppingStoneRecentBarrier1();
boost::shared_ptr< hudson_ms::PopulationStructure > makePopSteppingStoneRecentBarrier2();

boost::shared_ptr< hudson_ms::PopulationStructure > makeTestPopForSS1();




#endif
