
// declarations for testing emulation of ms args code


#ifndef __TEST_MS_ARGS__
#define __TEST_MS_ARGS__

#include <mscplusplus/population_structure.hpp>
#include <mscplusplus/recomb_params.hpp>


#include <boost/smart_ptr.hpp>


void testArgsSub(boost::shared_ptr< hudson_ms::PopulationStructure > pop,
		const hudson_ms::HudsonMSRecombParams& recombParams);


void testNewArgs1();
void testNewArgs2();
void testMSArgPrinting();


#endif
