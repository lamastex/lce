
/*! \file     
\brief AutoMCMCFactory definitions.
*/

#include "automcmc_factory.hpp"

#include "MCMCGRDiagnosticPSRFLeaves.hpp"
#include "MCMCGRDiagnosticPSRFLogpost.hpp"
#include "MCMCGRDiagnosticPSRFCherries.hpp"
#include "MCMCGRDiagnosticPSRFAvDepth.hpp"
#include "MCMCGRDiagnosticIntervalLeaves.hpp"
#include "MCMCGRDiagnosticIntervalCherries.hpp"
#include "MCMCGRDiagnosticIntervalAvDepth.hpp"


#include <stdexcept>
#include <vector>
//#include <boost/lambda/lambda.hpp>
//#include <boost/lambda/construct.hpp>
//#include <boost/lambda/bind.hpp>


using namespace std;

namespace ahabc {

	// static members
	AutoMCMCFactory::AutoTypeMap AutoMCMCFactory::typemap;

	bool AutoMCMCFactory::_dummy = loadTypeMap(typemap);


	std::vector< subpavings::MCMCGRAuto::Diagnostic * >& 
		AutoMCMCFactory::makeAutoMCMCDiags(
			AUTOTYPE autotype,
			std::vector< subpavings::MCMCGRAuto::Diagnostic * >& diagPtrs,
			real tolerance,
			size_t samplingInterval,
			size_t maxStatesForCalcs )
	{
		/* parameters controlling the diagnostics */
		double percent = 0.8;
	
		return makeAutoMCMCDiags(
			autotype,
			diagPtrs,
			tolerance, samplingInterval, 
			maxStatesForCalcs, percent);
	}

	
	std::vector< subpavings::MCMCGRAuto::Diagnostic * >& 
		AutoMCMCFactory::makeAutoMCMCDiags(
			AUTOTYPE autotype,
			std::vector< subpavings::MCMCGRAuto::Diagnostic * >& diagPtrs,
			real tolerance)
	{
		/* parameters controlling the diagnostics */
		size_t samplingInterval = 100;
		size_t maxStatesForCalcs = 1000000;
		double percent = 0.8;
	
		return makeAutoMCMCDiags(
			autotype,
			diagPtrs,
			tolerance, samplingInterval, 
			maxStatesForCalcs, percent);
	}	

	std::vector< subpavings::MCMCGRAuto::Diagnostic * >& 
		AutoMCMCFactory::makeAutoMCMCDiags(
			AUTOTYPE autotype,
			std::vector< subpavings::MCMCGRAuto::Diagnostic * >& diagPtrs,
			real tolerance, size_t samplingInterval, 
			size_t maxStatesForCalcs, double percent)
	{
		
		/* load_model_info should have been called. */
		if (typemap.empty()) throw std::runtime_error("No typemap");
			//loadTypeMap(model_info);

		AutoTypeMap::iterator it = typemap.find(autotype);

		if (it == typemap.end()) {
			std::ostringstream stm;
				stm << ":\nUnknown type: " << autotype << std::endl;
				throw std::invalid_argument(
				std::string("AutoMCMCFactory::createAutoMCMC(AUTOTYPE, ")
				+ std::string("unsigned long int, cxsc::real)"));
		}
			
		

		return  it->second.createFunc(diagPtrs, tolerance, 
					samplingInterval, maxStatesForCalcs, percent);
		
	}
	
	std::string AutoMCMCFactory::getDiagnosticName(AUTOTYPE autotype)
	{
		AutoTypeMap::iterator it = typemap.find(autotype);

		if (it == typemap.end()) {
			std::ostringstream stm;
				stm << ":\nUnknown type: " << autotype << std::endl;
				throw std::invalid_argument(
				std::string("AutoMCMCFactory::getDiagnosticName(AUTOTYPE)"));
		}
			
		return  it->second.name;
	}
	
				


	// private

	bool AutoMCMCFactory::loadTypeMap(AutoTypeMap& m) {
		m.clear();
		
		{
			AutoInfo a;
			a.name = string("Leaves");
			a.createFunc = &getAutoMCMCDiagsLeaves;
			m[LEAVES] = a;
		}
		{
			AutoInfo a;
			a.name = string("Cherries");
			a.createFunc = &getAutoMCMCDiagsCherries;
			m[CHERRIES] = a;
		}
		{
			AutoInfo a;
			a.name = string("AvDepth");
			a.createFunc = &getAutoMCMCDiagsAvDepth;
			m[AVDEPTH] = a;
		}
		{
			AutoInfo a;
			a.name = string("LeavesCherriesAvDepth");
			a.createFunc = &getAutoMCMCDiagsAllThree;
			m[ALLTHREE] = a;
		}
		{
			AutoInfo a;
			a.name = string("LeavesCherries");
			a.createFunc = &getAutoMCMCDiagsLeavesAndCherries;
			m[LEAVESANDCHERRIES] = a;
		}
		return true;
		
	}

	std::vector< subpavings::MCMCGRAuto::Diagnostic * >&
		AutoMCMCFactory::getAutoMCMCDiagsLeaves(
			std::vector< subpavings::MCMCGRAuto::Diagnostic * >& diagPtrs,
			real tolerance, size_t samplingInterval, 
			size_t maxStatesForCalcs, double percent)
	{
		
		cout << " \nTesting convergence using leaves (interval) \n" << endl;
		
		subpavings::MCMCGRDiagnosticIntervalLeaves* diagPtr = 
			new subpavings::MCMCGRDiagnosticIntervalLeaves(tolerance, 
							samplingInterval, maxStatesForCalcs,
							percent);
			
		diagPtrs.push_back(diagPtr);
		
		return diagPtrs;
			
	}
	
	std::vector< subpavings::MCMCGRAuto::Diagnostic * >&
		AutoMCMCFactory::getAutoMCMCDiagsCherries(
			std::vector< subpavings::MCMCGRAuto::Diagnostic * >& diagPtrs,
			real tolerance, size_t samplingInterval, 
			size_t maxStatesForCalcs, double percent)
	{
		
		cout << " \nTesting convergence using cherries (interval) \n" << endl;
		
		subpavings::MCMCGRDiagnosticIntervalCherries* diagPtr = 
			new subpavings::MCMCGRDiagnosticIntervalCherries(tolerance, 
							samplingInterval, maxStatesForCalcs,
							percent);
							
		diagPtrs.push_back(diagPtr);
		
		return diagPtrs;
			
	}
	
	std::vector< subpavings::MCMCGRAuto::Diagnostic * >&
		AutoMCMCFactory::getAutoMCMCDiagsAvDepth(
			std::vector< subpavings::MCMCGRAuto::Diagnostic * >& diagPtrs,
			real tolerance, size_t samplingInterval, 
			size_t maxStatesForCalcs, double percent)
	{
		
		cout << " \nTesting convergence using av depth (interval) \n" << endl;
		
		subpavings::MCMCGRDiagnosticIntervalAvDepth* diagPtr = 
			new subpavings::MCMCGRDiagnosticIntervalAvDepth(tolerance, 
							samplingInterval, maxStatesForCalcs,
							percent);
			
		diagPtrs.push_back(diagPtr);
		
		return diagPtrs;
			
	}
	
	std::vector< subpavings::MCMCGRAuto::Diagnostic * >&
		AutoMCMCFactory::getAutoMCMCDiagsAllThree(
			std::vector< subpavings::MCMCGRAuto::Diagnostic * >& diagPtrs,
			real tolerance, size_t samplingInterval, 
			size_t maxStatesForCalcs, double percent)
	{
		getAutoMCMCDiagsLeaves( diagPtrs, tolerance, 
							samplingInterval, maxStatesForCalcs,
							percent);
		getAutoMCMCDiagsCherries( diagPtrs, tolerance, 
							samplingInterval, maxStatesForCalcs,
							percent);
		getAutoMCMCDiagsAvDepth( diagPtrs, tolerance, 
							samplingInterval, maxStatesForCalcs,
							percent);
		return diagPtrs;
			
	}
	
	std::vector< subpavings::MCMCGRAuto::Diagnostic * >&
		AutoMCMCFactory::getAutoMCMCDiagsLeavesAndCherries(
			std::vector< subpavings::MCMCGRAuto::Diagnostic * >& diagPtrs,
			real tolerance, size_t samplingInterval, 
			size_t maxStatesForCalcs, double percent)
	{
		getAutoMCMCDiagsLeaves( diagPtrs, tolerance, 
							samplingInterval, maxStatesForCalcs,
							percent);
		getAutoMCMCDiagsCherries( diagPtrs, tolerance, 
							samplingInterval, maxStatesForCalcs,
							percent);
		return diagPtrs;
			
	}


}
