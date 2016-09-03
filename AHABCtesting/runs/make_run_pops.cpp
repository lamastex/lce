
/*
 * definitions for making populations for runs for Sharyn
*/


#include "config.h"

#include "make_run_pops.hpp"

#include <cmath>
#include <stdexcept>
#include <iostream>

using namespace std;



/* Two subpopulations, 30 individuals in each sub-pop sample;
 * Each has haploid population size 1000 = diploid N0 500;
 * No migration (island model);
 * Subpopulations separated 20 generations in the past;
 * No growth. 
 * */
boost::shared_ptr< hudson_ms::PopulationStructure > makePopScenario1(size_t nZero)
{
	try {
	
		boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
		
		{
			// add sub-pop 1
			size_t subNsamSize = 30;
			pop->add(subNsamSize);
		}
		{
			// add sub-pop 2
			size_t subNsamSize = 30;
			pop->add(subNsamSize);
		}
		
		// indices to make it easy to refer to subpops
		int pop1 = 0;
		int pop2 = 1;
		
		double T1gen = 20.0; // time back in generations
		double T1 = T1gen/(4*nZero); // scaled to units of 4N0 generations
		
		// past demographic events
		std::vector< hudson_ms::HudsonMSDemogEventPtr > _eList;
	
		{
			// lineages in second pop moved to first pop at time T1
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_j(T1, pop2, pop1)); 
			_eList.push_back(ev);
		}
		{
			// first pop instantaneous change in size at T1
			// this effectively adds total size of second pop to first at T1
			double size_ratio = 2.0; // ie N0 + N0
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_n(T1, pop1, size_ratio));
			_eList.push_back(ev);
		}
		
		pop->setEventList(_eList);
			
		return pop;	
	}
	catch (std::exception& e) {
	
		std::cout << "\nError in makePopScenario1:\n" + std::string(e.what()) << std::endl;
		throw;
	}
}

/* Two subpopulations, 30 individuals in each sub-pop sample;
 * Pop1 has haploid population size 1000 = diploid N0 500;
 * Pop2 has haploid population size 100 (1/10 of N0 for Pop1);
 * No migration (island model);
 * Subpopulations separated 20 generations in the past;
 * No growth. 
 * */
boost::shared_ptr< hudson_ms::PopulationStructure > makePopScenario2(size_t nZero)
{
	try {
	
		boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
		
		{
			// add sub-pop 1
			size_t subNsamSize = 30;
			pop->add(subNsamSize);
		}
		{
			// add sub-pop 2
			size_t subNsamSize = 30;
			double growth = 0.0;
			double size_ratio = 0.1;
			pop->add(subNsamSize, size_ratio, growth);
		}
		
		// indices to make it easy to refer to subpops
		int pop1 = 0;
		int pop2 = 1;
		
		double T1gen = 20.0; // time back in generations
		double T1 = T1gen/(4*nZero); // scaled to units of 4N0 generations
		
		// past demographic events
		std::vector< hudson_ms::HudsonMSDemogEventPtr > _eList;
	
		{
			// lineages in second pop moved to first pop at time T1
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_j(T1, pop2, pop1)); 
			_eList.push_back(ev);
		}
		{
			// first pop instantaneous change in size at T1
			// this effectively add total size of second pop to first at T1
			double size_ratio = 1.1; // ie N0 + 0.1N0
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_n(T1, pop1, size_ratio));
			_eList.push_back(ev);
		}
		
		pop->setEventList(_eList);
			
		return pop;	
	}
	catch (std::exception& e) {
	
		std::cout << "\nError in makePopScenario2:\n" + std::string(e.what()) << std::endl;
		throw;
	}
}


/* Two subpopulations, 30 individuals in each sub-pop sample;
 * Pop1 has haploid population size 1000 = diploid N0 500;
 * Pop2 has haploid population size 100 (1/10 of N0 for Pop1);
 * Migration into Pop2 from Pop1: fraction of Pop2 made of migrants from Pop1 = 0.01 ;
 * Subpopulations separated 20 generations in the past;
 * No growth. 
 * */
boost::shared_ptr< hudson_ms::PopulationStructure > makePopScenario3(size_t nZero)
{
	try {
	
		boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
		
		{
			// add sub-pop 1
			size_t subNsamSize = 30;
			pop->add(subNsamSize);
		}
		{
			// add sub-pop 2
			size_t subNsamSize = 30;
			double growth = 0.0;
			double size_ratio = 0.1;
			pop->add(subNsamSize, size_ratio, growth);
		}
		
		// indices to make it easy to refer to subpops
		int pop1 = 0;
		int pop2 = 1;
		
		// set migration matrix elements
		{
			double _m = 0.01; 
			pop->setMigration(pop2, pop1, _m);
		}
		
		double T1gen = 20.0; // time back in generations
		double T1 = T1gen/(4*nZero); // scaled to units of 4N0 generations
		// past demographic events
		std::vector< hudson_ms::HudsonMSDemogEventPtr > _eList;
	
		{
			// lineages in second pop moved to first pop at time T1
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_j(T1, pop2, pop1)); 
			_eList.push_back(ev);
		}
		{
			// first pop instantaneous change in size at T1
			// this effectively add total size of second pop to first at T1
			double size_ratio = 1.1; // ie N0 + 0.1N0
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_n(T1, pop1, size_ratio));
			_eList.push_back(ev);
		}
		
		pop->setEventList(_eList);
			
		return pop;	
	}
	catch (std::exception& e) {
	
		std::cout << "\nError in makePopScenario3:\n" + std::string(e.what()) << std::endl;
		throw;
	}
}

/* Three subpopulations (first two our pops 1 and 2, third 
 * is the big population they came from)
 * Pops 1 and 2 have 30 individuals in each sub-pop sample;
 * Each of pops 1 and 2 has current haploid population size 1000 = diploid N0 500;
 * They were each 0.01 of current size (ie 0.01N0, 1/100 of current size) 18 generations ago
 * and grew exponentially from then on.
 * Before 18 generations ago there was no growth in pops 1 and 2.
 * At 20 generations ago pops 1 and 2 join - consider as 2 joining 1
 * with combined size 0.02N0.  
 * At 21 generations ago pop 1 (which has subsumed pop2) joins pop 3 - the main pop
 * Total population size (for pop 3) 10N0
 * No migration (island model);
 * I think we can just let pop3 continue on its merry way from 21 generations ago
 * - we have no samples from it - ie treat as having no popsize change etc
 * */
boost::shared_ptr< hudson_ms::PopulationStructure > makePopScenario4(size_t nZero)
{
	try {
	
		boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
		
		double T1gen = 18.0; // time back in generations
		double T1 = T1gen/(4*nZero); // scaled to units of 4N0 generations
		double T2gen = 20.0; // time back in generations
		double T2 = T2gen/(4*nZero); // scaled to units of 4N0 generations
		double T3gen = 21.0; // time back in generations
		double T3 = T3gen/(4*nZero); // scaled to units of 4N0 generations
		
		/* work out a growth rate for pops 1 and 2 so that at T1 both were 
		 * 0.01 of current size 
		 * ie solve exp(-g * T1) = 0.1 or exp (g * T1) = 100
		 * g = 1/T1 * ln(100) */
		double g = log(100)/T1;
		
		{
			// add sub-pop 1
			size_t subNsamSize = 30;
			double growth = g;
			double size_ratio = 1.0;
			pop->add(subNsamSize, size_ratio, growth);
		}
		{
			// add sub-pop 2
			size_t subNsamSize = 30;
			double growth = g;
			double size_ratio = 1.0;
			pop->add(subNsamSize, size_ratio, growth);
		}
		
		{
			// add sub-pop 3 --- no samples from this one
			size_t subNsamSize = 0;
			double growth = 0.0;
			double size_ratio = 10.0 - 0.02;
			pop->add(subNsamSize, size_ratio, growth);
		}
		
		// indices to make it easy to refer to subpops
		int pop1 = 0;
		int pop2 = 1;
		int pop3 = 2;
		
		
		// past demographic events
		std::vector< hudson_ms::HudsonMSDemogEventPtr > _eList;
	
		{
			// zero growth before time T1 for all subpopulations
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_G(T1, 0.0)); 
			_eList.push_back(ev);
		}
		{
			// lineages in second pop moved to first pop at time T2
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_j(T2, pop2, pop1)); 
			_eList.push_back(ev);
		}
		{
			// first pop instantaneous change in size at T2
			// this effectively adds total size of second pop to first at T2
			double size_ratio = 0.02; // ie N0/100 + N0/100
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_n(T2, pop1, size_ratio));
			_eList.push_back(ev);
		}
		{
			// lineages in first pop moved to third pop at time T3
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_j(T3, pop1, pop3)); 
			_eList.push_back(ev);
		}
		{
			// third pop instantaneous change in size at T3
			// this effectively adds total size of first pop to third at T3
			double size_ratio = 10; 
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_n(T3, pop3, size_ratio));
			_eList.push_back(ev);
		}
		
		pop->setEventList(_eList);
			
		return pop;	
	}
	catch (std::exception& e) {
	
		std::cout << "\nError in makePopScenario4:\n" + std::string(e.what()) << std::endl;
		throw;
	}
}


/* Three subpopulations (first two our pops 1 and 2, third 
 * is the big population they came from)
 * Pops 1 and 2 have 30 individuals in each sub-pop sample;
 * Each of pops 1 and 2 has current haploid population size 1000 = diploid N0 500;
 * They were each 0.01 of current size (ie 0.01N0, 1/100 of current size) 18 generations ago
 * and grew exponentially from then on.
 * Before 18 generations ago there was no growth in pops 1 and 2.
 * At 20 generations ago pops 1 and 2 join - consider as 2 joining 1
 * with combined size 0.02N0.  
 * At 21 generations ago pop 1 (which has subsumed pop2) joins pop 3 - the main pop
 * Total population size (for pop 3) 10N0
 * No migration (island model);
 * 
 * This scenario just demonstrates that what happens to pop3 after the others
 * split off is irrelevant because we have no samples from it. 
 * */
boost::shared_ptr< hudson_ms::PopulationStructure > makePopScenario4alt(size_t nZero)
{
	try {
	
		boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
		
		double T1gen = 18.0; // time back in generations
		double T1 = T1gen/(4*nZero); // scaled to units of 4N0 generations
		double T2gen = 20.0; // time back in generations
		double T2 = T2gen/(4*nZero); // scaled to units of 4N0 generations
		double T3gen = 21.0; // time back in generations
		double T3 = T3gen/(4*nZero); // scaled to units of 4N0 generations
		
		/* work out a growth rate for pops 1 and 2 so that at T1 both were 
		 * 0.01 of current size 
		 * ie solve exp(-g * T1) = 0.1 or exp (g * T1) = 100
		 * g = 1/T1 * ln(100) */
		double g = log(100)/T1;
		
		/* work out a growth rate for pop3 so that at T3 it had size 
		 * 10 - 0.01 - 0.01  
		 * ie solve exp(-g * T3) = 10 -0.01 -0.01 or exp (g * T3) = 1/(10 -0.01 -0.01)
		 * g = 1/T1 * ln(1/(10 -0.01 -0.01)) */
		double g3 = log(1.0/(10-0.01-0.01))/T3;
		
		{
			// add sub-pop 1
			size_t subNsamSize = 30;
			double growth = g;
			double size_ratio = 1.0;
			pop->add(subNsamSize, size_ratio, growth);
		}
		{
			// add sub-pop 2
			size_t subNsamSize = 30;
			double growth = g;
			double size_ratio = 1.0;
			pop->add(subNsamSize, size_ratio, growth);
		}
		
		{
			// add sub-pop 3 --- no samples from this one
			size_t subNsamSize = 0;
			double growth = g3;
			double size_ratio = 1.0;
			pop->add(subNsamSize, size_ratio, growth);
		}
		
		// indices to make it easy to refer to subpops
		int pop1 = 0;
		int pop2 = 1;
		int pop3 = 2;
		
		
		// past demographic events
		std::vector< hudson_ms::HudsonMSDemogEventPtr > _eList;
	
		{
			// zero growth before time T1 for all subpopulations
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_G(T1, 0.0)); 
			_eList.push_back(ev);
		}
		{
			// lineages in second pop moved to first pop at time T2
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_j(T2, pop2, pop1)); 
			_eList.push_back(ev);
		}
		{
			// first pop instantaneous change in size at T2
			// this effectively adds total size of second pop to first at T2
			double size_ratio = 0.02; // ie N0/100 + N0/100
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_n(T2, pop1, size_ratio));
			_eList.push_back(ev);
		}
		{
			// lineages in first pop moved to third pop at time T3
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_j(T3, pop1, pop3)); 
			_eList.push_back(ev);
		}
		{
			// third pop instantaneous change in size at T3
			// this effectively adds total size of first pop to third at T3
			double size_ratio = 10; 
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_n(T3, pop3, size_ratio));
			_eList.push_back(ev);
		}
		
		pop->setEventList(_eList);
			
		return pop;	
	}
	catch (std::exception& e) {
	
		std::cout << "\nError in makePopScenario4alt:\n" + std::string(e.what()) << std::endl;
		throw;
	}
}


/* Two subpopulations (our pops 1 and 2)
 * Pops 1 and 2 have 30 individuals in each sub-pop sample;
 * Each of pops 1 and 2 has current haploid population size 1000 = diploid N0 500;
 * They were each 0.01 of current size (ie 0.01N0, 1/100 of current size) 18 generations ago
 * and grew exponentially from then on.
 * Before 18 generations ago there was no growth in pops 1 and 2.
 * At 20 generations ago pops 1 and 2 join - consider as 2 joining 1
 * with combined size 0.02N0.  
 * No migration (island model);
 * So we just ignore pop1+pop2 coming from some bigger pop
 * */
boost::shared_ptr< hudson_ms::PopulationStructure > makePopScenario4contrast(size_t nZero)
{
	try {
	
		boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
		
		double T1gen = 18.0; // time back in generations
		double T1 = T1gen/(4*nZero); // scaled to units of 4N0 generations
		double T2gen = 20.0; // time back in generations
		double T2 = T2gen/(4*nZero); // scaled to units of 4N0 generations
		
		/* work out a growth rate for pops 1 and 2 so that at T1 both were 
		 * 0.01 of current size 
		 * ie solve exp(-g * T1) = 0.1 or exp (g * T1) = 100
		 * g = 1/T1 * ln(100) */
		double g = log(100)/T1;
		
		{
			// add sub-pop 1
			size_t subNsamSize = 30;
			double growth = g;
			double size_ratio = 1.0;
			pop->add(subNsamSize, size_ratio, growth);
		}
		{
			// add sub-pop 2
			size_t subNsamSize = 30;
			double growth = g;
			double size_ratio = 1.0;
			pop->add(subNsamSize, size_ratio, growth);
		}
		
		
		// indices to make it easy to refer to subpops
		int pop1 = 0;
		int pop2 = 1;
		
		// past demographic events
		std::vector< hudson_ms::HudsonMSDemogEventPtr > _eList;
	
		{
			// zero growth before time T1 for all subpopulations
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_G(T1, 0.0)); 
			_eList.push_back(ev);
		}
		{
			// lineages in second pop moved to first pop at time T2
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_j(T2, pop2, pop1)); 
			_eList.push_back(ev);
		}
		{
			// first pop instantaneous change in size at T2
			// this effectively adds total size of second pop to first at T2
			double size_ratio = 0.02; // ie N0/100 + N0/100
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_n(T2, pop1, size_ratio));
			_eList.push_back(ev);
		}
		
		
		pop->setEventList(_eList);
			
		return pop;	
	}
	catch (std::exception& e) {
	
		std::cout << "\nError in makePopScenario4contrast:\n" + std::string(e.what()) << std::endl;
		throw;
	}
}

/* Three subpopulations (first two our pops 1 and 2, third 
 * is the big population they came from)
 * Pops 1 and 2 have 30 individuals in each sub-pop sample;
 * Each of pop 1 has current haploid population size 1000 = diploid N0 500;
 * Each of pop 2 has current haploid population size 0.1N0
 * They were each 0.01 of current size (ie 0.01N0, 1/100 of current size) 18 generations ago
 * and grew exponentially from then on.
 * Before 18 generations ago there was no growth in pops 1 and 2.
 * At 20 generations ago pops 1 and 2 join - consider as 2 joining 1
 * with combined size 0.011N0.  
 * At 21 generations ago pop 1 (which has subsumed pop2) joins pop 3 - the main pop
 * Total population size (for pop 3) 10N0
 * No migration (island model);
 * I think we can just let pop3 continue on its merry way from 21 generations ago
 * - we have no samples from it - ie treat as having no popsize change etc
 * */
boost::shared_ptr< hudson_ms::PopulationStructure > makePopScenario5(size_t nZero)
{
	try {
	
		boost::shared_ptr< hudson_ms::PopulationStructure > pop ( new hudson_ms::PopulationStructure );
		
		double T1gen = 18.0; // time back in generations
		double T1 = T1gen/(4*nZero); // scaled to units of 4N0 generations
		double T2gen = 20.0; // time back in generations
		double T2 = T2gen/(4*nZero); // scaled to units of 4N0 generations
		double T3gen = 21.0; // time back in generations
		double T3 = T3gen/(4*nZero); // scaled to units of 4N0 generations
		
		/* work out a growth rate for pops 1 and 2 so that at T1 both were 
		 * 0.01 of current size 
		 * ie solve exp(-g * T1) = 0.1 or exp (g * T1) = 100
		 * g = 1/T1 * ln(100) */
		double g = log(100)/T1;
		
		{
			// add sub-pop 1
			size_t subNsamSize = 30;
			double growth = g;
			double size_ratio = 1.0;
			pop->add(subNsamSize, size_ratio, growth);
		}
		{
			// add sub-pop 2
			size_t subNsamSize = 30;
			double growth = g;
			double size_ratio = 0.1;
			pop->add(subNsamSize, size_ratio, growth);
		}
		
		{
			// add sub-pop 3 --- no samples from this one
			size_t subNsamSize = 0;
			double growth = 0.0;
			double size_ratio = 10.0 - 0.01 - 0.001;
			pop->add(subNsamSize, size_ratio, growth);
		}
		
		// indices to make it easy to refer to subpops
		int pop1 = 0;
		int pop2 = 1;
		int pop3 = 2;
		
		
		// past demographic events
		std::vector< hudson_ms::HudsonMSDemogEventPtr > _eList;
	
		{
			// zero growth before time T1 for all subpopulations
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_G(T1, 0.0)); 
			_eList.push_back(ev);
		}
		{
			// lineages in second pop moved to first pop at time T2
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_j(T2, pop2, pop1)); 
			_eList.push_back(ev);
		}
		{
			// first pop instantaneous change in size at T2
			// this effectively adds total size of second pop to first at T2
			double size_ratio = 0.01 + 0.001; // ie N0/100 + 0.1N0/100
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_n(T2, pop1, size_ratio));
			_eList.push_back(ev);
		}
		{
			// lineages in first pop moved to third pop at time T3
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_j(T3, pop1, pop3)); 
			_eList.push_back(ev);
		}
		{
			// third pop instantaneous change in size at T3
			// this effectively adds total size of first pop to third at T3
			double size_ratio = 10; 
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > 
				ev(new hudson_ms::HudsonMSDemogEvent_n(T3, pop3, size_ratio));
			_eList.push_back(ev);
		}
		
		pop->setEventList(_eList);
			
		return pop;	
	}
	catch (std::exception& e) {
	
		std::cout << "\nError in makePopScenario5:\n" + std::string(e.what()) << std::endl;
		throw;
	}
}
