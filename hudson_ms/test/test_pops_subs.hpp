
/*
 * declarations for populations testing subroutines
*/

#ifndef __TEST_POPS_SUBS__
#define __TEST_POPS_SUBS__

#include <mscplusplus/population_structure.hpp>
#include <mscplusplus/recomb_params.hpp>


#include <boost/smart_ptr.hpp>

#include <vector>
#include <string>



std::string string_from_file(const std::string name);

std::string getDefaultOutputDir();

void makeOutputDir(const std::string& outputDir);

void testPopStructSub(const hudson_ms::HudsonMSRecombParams& recombParams,
		const boost::shared_ptr< hudson_ms::PopulationStructure >& pop,
		const std::string& filename,
		int seed);


void testPopConstructorSub(const hudson_ms::PopulationStructure& pop);

void testGetAndSetSimple(hudson_ms::PopulationStructure& pop, int order);
void testGetAndSetMigrationSimple(hudson_ms::PopulationStructure& pop);

void testPopOrderingSub(hudson_ms::PopulationStructure& pop,
						const std::vector < std::string >& newOrder);

void testPopOutputWithMigCheck(const hudson_ms::PopulationStructure& pop);

// expecting a failure in migCheck
void testPopOutputWithMigCheckFailure(const hudson_ms::PopulationStructure& pop);

void testPopMigSetSub(hudson_ms::PopulationStructure& pop,
			const std::vector < std::string >& labels,
			int order1, int order21, int order22,
			double m1, double m2);

void testPopMigSetSubByLabel(hudson_ms::PopulationStructure& pop,
			const std::vector < std::string >& labels,
			int order1, int order21, double m1);
			
void testPopMigSetSubByOrder(hudson_ms::PopulationStructure& pop,
			const std::vector < std::string >& labels,
			int order1, int order22, double m2);

void testPopEraseByOrder(hudson_ms::PopulationStructure& pop,
			const std::vector < std::string >& labels,
			int orderErase);

void testPopEraseByLabel(hudson_ms::PopulationStructure& pop,
			const std::vector < std::string >& labels,
			int orderErase);

void testPopRelabellingSub(hudson_ms::PopulationStructure& pop,
						const std::vector < std::string >& newlabels);

void testPopRelabellingSubByOrder(hudson_ms::PopulationStructure& pop,
						std::vector < std::string >& labels,
						int order, const std::string& newlabel);

void testPopRelabellingSubByLabel(hudson_ms::PopulationStructure& pop,
						std::vector < std::string >& labels,
						int order, const std::string& newlabel);

						
void testPopMigCheckSub(boost::shared_ptr< hudson_ms::PopulationStructure > pop);

void popOutOfRangeTestSub(const hudson_ms::PopulationStructure& pop, int index);

void popNotRecognisedTestSub(const hudson_ms::PopulationStructure& pop, 
							const std::string& label);


void testSubPopsInvalidArgSub(int nsam, double sr, double gr);

void testSubPopsSetInvalidArgSub(const hudson_ms::PopulationStructure& pop,
									int nsam, double sr, double gr);

void testPopRelabellingIndividualInvalidArgSub(const hudson_ms::PopulationStructure& pop,
					int index1, std::string oldLabel, std::string newLabel, bool duplicate);
						
void testPopRelabellingIndividualInvalidArg(const hudson_ms::PopulationStructure& pop,
										int index1, int index2);

void testPopRelabellingInvalidArg(const hudson_ms::PopulationStructure& pop);

void testPopBeforeAndAfterEventsSub(const hudson_ms::PopulationStructure& pop);

void testPopBeforeAndAfterEvents(const hudson_ms::PopulationStructure& pop);
										
void throwLogicErrorCheck(bool outOfRange);


#endif
