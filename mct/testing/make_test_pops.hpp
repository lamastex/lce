
/*
 * declarations for making test populations for testing
*/

#ifndef __MCT_MAKE_TEST_POPS__
#define __MCT_MAKE_TEST_POPS__

#include <mscplusplus/population_structure.hpp>



#include <boost/smart_ptr.hpp>



boost::shared_ptr< hudson_ms::PopulationStructure > makePopSingle();
boost::shared_ptr< hudson_ms::PopulationStructure > makePopSingleZero();
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandOneZero(size_t Nzero);
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandThreeOneZero(size_t Nzero);
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandFourTwoZero(size_t Nzero);
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandSmall(size_t Nzero);
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandThreeSmall(size_t Nzero);
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandThreeSmallest(size_t Nzero);
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandBasic(size_t Nzero);
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandComplex(size_t Nzero);
boost::shared_ptr< hudson_ms::PopulationStructure > makePopInstantSizeChange();
boost::shared_ptr< hudson_ms::PopulationStructure > makePopOutgroupSequence();
boost::shared_ptr< hudson_ms::PopulationStructure > makePopInstantSizeChangeAndExpGrowth();
boost::shared_ptr< hudson_ms::PopulationStructure > makePopTwoSpeciesWithPopSizeDiff();
boost::shared_ptr< hudson_ms::PopulationStructure > makePopBackwardsSplit();
boost::shared_ptr< hudson_ms::PopulationStructure > makePopSteppingStoneRecentBarrier1();
boost::shared_ptr< hudson_ms::PopulationStructure > makePopSteppingStoneRecentBarrier2();

boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandBalanced(size_t Nzero);
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandBalLessComm(size_t Nzero);
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandMSstats(size_t Nzero);
boost::shared_ptr< hudson_ms::PopulationStructure > makePopIslandPTstats(size_t Nzero);


boost::shared_ptr< hudson_ms::PopulationStructure > makePopRaazRuns(
				size_t subNsamSizeEach,
				double alpha, double beta);

#endif
