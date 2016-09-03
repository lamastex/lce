
/*! \file     
\brief AutoMCMCFactory declarations.
*/

#ifndef __INC_AUTOMCMC_FACTORY__H
#define __INC_AUTOMCMC_FACTORY__H

#include "MCMCGRAuto.hpp"

#include "real.hpp"

//#include <vector>
#include <string>
#include <map>

#include <boost/function.hpp>
#include <boost/smart_ptr.hpp>


namespace ahabc {
	
	/* A class to make automcmc objects. */
	class AutoMCMCFactory {
		public:
		
			typedef enum{LEAVES = 1, CHERRIES = 10, AVDEPTH = 100,
						ALLTHREE = 111, LEAVESANDCHERRIES = 11}
							AUTOTYPE;
		
			//explicit AutoMCMCFactory(AUTO_TYPE a);

			static std::vector< subpavings::MCMCGRAuto::Diagnostic * >& 
			makeAutoMCMCDiags(
					AUTOTYPE autotype,
					std::vector< subpavings::MCMCGRAuto::Diagnostic * >& diagPtrs,
					real tolerance,
					size_t samplingInterval,
					size_t maxStatesForCalcs );

			static std::vector< subpavings::MCMCGRAuto::Diagnostic * >&
				makeAutoMCMCDiags(
					AUTOTYPE autotype,
					std::vector< subpavings::MCMCGRAuto::Diagnostic * >& diagPtrs,
					real tolerance);
			
			static std::vector< subpavings::MCMCGRAuto::Diagnostic * >& 
			makeAutoMCMCDiags(
				AUTOTYPE autotype,
				std::vector< subpavings::MCMCGRAuto::Diagnostic * >& diagPtrs,
				real tolerance, size_t samplingInterval, 
				size_t maxStatesForCalcs, double percent);
			
			static std::string getDiagnosticName(AUTOTYPE autotype);
			
		private:
		
			AutoMCMCFactory();
		
			struct AutoInfo {
				std::string name;
				boost::function< std::vector< subpavings::MCMCGRAuto::Diagnostic * >&
					(std::vector< subpavings::MCMCGRAuto::Diagnostic * >& diagPtrs,
					cxsc::real tolerance, size_t samplingInterval, 
					size_t maxStatesForCalcs, double percent)> createFunc;
				
			};
			
			typedef std::map < AUTOTYPE, AutoInfo > AutoTypeMap;
			
		
			static AutoTypeMap typemap;
			
			static bool _dummy;
					
			static bool loadTypeMap(AutoTypeMap& m);
			
			static std::vector< subpavings::MCMCGRAuto::Diagnostic * >&
				getAutoMCMCDiagsLeaves(
					std::vector< subpavings::MCMCGRAuto::Diagnostic * >& diagPtrs,
					real tolerance, size_t samplingInterval, 
					size_t maxStatesForCalcs, double percent);
			
			static std::vector< subpavings::MCMCGRAuto::Diagnostic * >&
				getAutoMCMCDiagsCherries(
					std::vector< subpavings::MCMCGRAuto::Diagnostic * >& diagPtrs,
					real tolerance, size_t samplingInterval, 
					size_t maxStatesForCalcs, double percent);
			
			static std::vector< subpavings::MCMCGRAuto::Diagnostic * >&
				getAutoMCMCDiagsAvDepth(
					std::vector< subpavings::MCMCGRAuto::Diagnostic * >& diagPtrs,
					real tolerance, size_t samplingInterval, 
					size_t maxStatesForCalcs, double percent);
					
			static std::vector< subpavings::MCMCGRAuto::Diagnostic * >&
				getAutoMCMCDiagsAllThree(
					std::vector< subpavings::MCMCGRAuto::Diagnostic * >& diagPtrs,
					real tolerance, size_t samplingInterval, 
					size_t maxStatesForCalcs, double percent);
			
			static std::vector< subpavings::MCMCGRAuto::Diagnostic * >&
				getAutoMCMCDiagsLeavesAndCherries(
					std::vector< subpavings::MCMCGRAuto::Diagnostic * >& diagPtrs,
					real tolerance, size_t samplingInterval, 
					size_t maxStatesForCalcs, double percent);		
			
			
	};

} // end namespace ahabc

#endif 

