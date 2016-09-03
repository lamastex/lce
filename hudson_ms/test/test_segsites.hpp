
// declarations for testing emulation of ms segsites code


#ifndef __TEST_SEGSITES__
#define __TEST_SEGSITES__

#include <mscplusplus/population_structure.hpp>
#include <mscplusplus/recomb_params.hpp>
#include <mscplusplus/segsite_params.hpp>


#include <boost/smart_ptr.hpp>


void testSegSitesSub(boost::shared_ptr< hudson_ms::PopulationStructure > pop,
	const hudson_ms::HudsonMSRecombParams& recombParams,
	const hudson_ms::HudsonMSSegSiteParams& ssParams);

void testMSSegsitesFramework(boost::shared_ptr< hudson_ms::PopulationStructure > pop);

void testMSSegsites1();
void testMSSegsites2();



#endif
