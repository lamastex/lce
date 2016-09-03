/* mscplusplus - ms for c++.

   Copyright (C) 2011, 2012 Jenny Harlow
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */


/*! \file     
\brief HudsonMSDemogEvent declarations.
*/

#ifndef __HUDSON_MS_DEMOG_EVENT_HPP__
#define __HUDSON_MS_DEMOG_EVENT_HPP__

#include <boost/smart_ptr.hpp>

//#include <exception>
#include <vector>
//#include <string>
#include <sstream>

namespace hudson_ms {
	

	//forward declarations
	class HudsonMSArgFactory;

	/*! \brief An abstract class for events in the Hudson MS world.
	
	All HudsonMSEvents have an event time, accessed using getTime().
	Event time is counted backwards from the present, ie if an event 
	\f$ a \f$ has a shorter time than an event \f$ b \f$ then
	\f$ a \f$ happens closer (in the past) to the present than \f$ b \f$.
	
	\internal
	Events and the HudsonMSArgFactory together use the visitor pattern:
	the arg factory accepts the event (presumably at the event time). 
	Accepting the event will result in the event's callEvent method being
	called, passing in the instance of the arg factory.  The event will
	then call the factory's own event process, passing in the appropriate
	information (eg time, populations involved, parameters).  This is 
	all a bit overkill here but it potentially gives extra flexibility
	by some decoupling of events and arg factory.  This means that
	events can be set up outside of the arg factory and passed in, and
	(if we want to extend) the same basic
	factory machinery can be used to process
	different kinds of events (only by adding different accept methods)
	and similarly the same events could be used for different types of
	factories.   
	*/
    class HudsonMSEvent
    {
		public:

			HudsonMSEvent (double _time);

			/*! \brief Destructor */
			virtual ~HudsonMSEvent () {};

			/*! \brief Call the event processing machinery. */
			virtual void callEvent(const hudson_ms::HudsonMSArgFactory& builder) const = 0;
		
			/*! \brief Get the time at which the event takes place. */
			virtual double getTime() const;
			
		protected:
			const double time;
    };
	
	/*! \brief An abstract class for demographic events in the 
	Hudson MS world.
	
	Demographic events can be ordered.  If two events have different
	times associated with them, the event with the smaller time is 
	considered to be the lesser of the two, ie closer (in the past)
	to the present.  If two events have exactly
	the same times, the events must provide other criteria for determining
	their ordering.
	*/
	class HudsonMSDemogEvent : public HudsonMSEvent
    {
		public:
		
			HudsonMSDemogEvent (double _time) : HudsonMSEvent(_time) {};

			/*! \brief Destructor */
			virtual ~HudsonMSDemogEvent () {};
			
			virtual std::string msString() const;
			
			bool operator < (const hudson_ms::HudsonMSDemogEvent& rhs) const;
			
		protected:

			static int precEst(double param);
			
				
		private:
		
		static const int default_ordering = 0;
		
		virtual int getOrdering() const;
		
		virtual std::ostringstream& outputDescription(
									std::ostringstream& stm) const = 0;


    };

	/*! \brief An event to change (rescale by some multiplier)
	all subpopulation sizes (\f$ N0 \f$ for each subpopulation).
	
	Event ordering should ensure that this event is ordered less than
	(happens before, in backwards time) any event changing growth rates
	that have the same event time. 	*/
	class HudsonMSDemogEvent_N : public HudsonMSDemogEvent
	{
		public :
			/*! \brief Constructor,
			
			\param _time The event time,
			\param _multNzero The multiple to be applied to 
			\f$ N_0 \f$ for all subpopulations.*/
			HudsonMSDemogEvent_N(double _time, double _multNzero);
			
			void callEvent(const hudson_ms::HudsonMSArgFactory& builder) const;
			
					
		private :
			HudsonMSDemogEvent_N();

			std::ostringstream& outputDescription(
									std::ostringstream& stm ) const;
									
			static const int ordering = -1;

			int getOrdering() const;
			
			const double multNzero;
			
	};
	
	/*! \brief An event to change (rescale by some multiplier)
	the size of a specified subpopulation.
	
	Event ordering should ensure that this event is ordered less than
	(happens before, in backwards time) any event changing growth rates
	that have the same event time. 	*/
	class HudsonMSDemogEvent_n : public HudsonMSDemogEvent
	{
		public :
			// _pop is index of pop to change, ie 0, 1, ... npop-1
			HudsonMSDemogEvent_n(double _time, int _pop, double _multNzero);
			
			void callEvent(const hudson_ms::HudsonMSArgFactory& builder) const;
			
		private :
			HudsonMSDemogEvent_n();

			std::ostringstream& outputDescription(
									std::ostringstream& stm ) const;
			
			static const int ordering = -1;

			int getOrdering() const;

			size_t pop;
			const double multNzero;
	};
	
	/*! \brief An event to change (reset to a specified rate)
	the growth rates of all subpopulations. */
	class HudsonMSDemogEvent_G : public HudsonMSDemogEvent
	{
		public :
			HudsonMSDemogEvent_G(double _time, double _newGrowthRate);
			
			void callEvent(const hudson_ms::HudsonMSArgFactory& builder) const;
			
		private :
			HudsonMSDemogEvent_G();

			std::ostringstream& outputDescription(
									std::ostringstream& stm ) const;
			const double newGrowthRate;
	};

	/*! \brief An event to change (reset to a specified rate)
	the growth rate for a specified subpopulation. */
	class HudsonMSDemogEvent_g : public HudsonMSDemogEvent
	{
		public :
			// _pop is index of pop to change, ie 0, 1, ... npop-1
			HudsonMSDemogEvent_g(double _time, int _pop, double _newGrowthRate);
			
			void callEvent(const hudson_ms::HudsonMSArgFactory& builder) const;
			
		private :
			HudsonMSDemogEvent_g();

			std::ostringstream& outputDescription(
									std::ostringstream& stm ) const;
			size_t pop;
			const double newGrowthRate;
	};


	/*! \brief An event to join (in backwards time) two
	subpopulations (equivalent to a split in forwards time).
	
	Event ordering should ensure that this event is ordered less than
	(happens before, in backwards time) any event changing subpopulation
	sizes that have the same event time. 	*/
	class HudsonMSDemogEvent_j : public HudsonMSDemogEvent
	{
		public :
			// _pop_i, _popj are indices of pops to join, ie 0, 1, ... npop-1
			HudsonMSDemogEvent_j(double _time, int _pop_i, int _pop_j);
		
			void callEvent(const hudson_ms::HudsonMSArgFactory& builder) const;
			
		private :
			HudsonMSDemogEvent_j();

			std::ostringstream& outputDescription(
									std::ostringstream& stm ) const;
			
			static const int ordering = -2;
			int getOrdering() const;
			
			size_t pop_i;
			size_t pop_j;
	};

	/*! \brief An event to split (in backwards time) two
	subpopulations 
	(equivalent to a subpopulation admixture in forwards time).
	
	Event ordering should ensure that this event is ordered less than
	(happens before, in backwards time) any event changing subpopulation
	sizes, growth rates, or migration rates, 
	that has the same event time. 	*/
	class HudsonMSDemogEvent_s : public HudsonMSDemogEvent
	{
		public :
			// _pop is indices of pop to split, ie 0, 1, ... npop-1
			HudsonMSDemogEvent_s(double _time, int _pop, double _probStaying);
		
			void callEvent(const hudson_ms::HudsonMSArgFactory& builder) const;
			
		private :
			HudsonMSDemogEvent_s();

			std::ostringstream& outputDescription(
									std::ostringstream& stm ) const;
			
			static const int ordering = -2;
			int getOrdering() const;

			size_t pop;
			const double probStaying;
	};

	
	/*! \brief An event to set all the elements of the migration 
	matrix for all subpopulations to the same value.*/
	class HudsonMSDemogEvent_M : public HudsonMSDemogEvent
	{
		public :
		
			HudsonMSDemogEvent_M(double _time, double _newMigRate);

			void callEvent(const hudson_ms::HudsonMSArgFactory& builder) const;
			
		private :
			HudsonMSDemogEvent_M();

			std::ostringstream& outputDescription(
									std::ostringstream& stm ) const;
			const double newMigRate;
	};

	/*! \brief An event to set an element of the migration 
	matrix, ie set the migration parameter for migration
	into a specified subpopulation from a specified subpopulation.*/
	class HudsonMSDemogEvent_m : public HudsonMSDemogEvent
	{
		public :
			// _pop_i, _popj are indices of pops for mig, ie 0, 1, ... npop-1
			HudsonMSDemogEvent_m(double _time, int _pop_i, int _pop_j, double _newMigRate);

			void callEvent(const hudson_ms::HudsonMSArgFactory& builder) const;
			
		private :
			HudsonMSDemogEvent_m();

			std::ostringstream& outputDescription(
									std::ostringstream& stm ) const;
			size_t pop_i;
			size_t pop_j;
			const double newMigRate;
	};
	
	/*! \brief An event to set each element of the migration 
	matrix individually, ie set the migration parameter for migration
	into each subpopulation from each other subpopulation.*/
	class HudsonMSDemogEvent_ma : public HudsonMSDemogEvent
	{
		public :
		
			HudsonMSDemogEvent_ma(double _time, const std::vector< std::vector < double > >& _newMigMatrix);
			
			void callEvent(const hudson_ms::HudsonMSArgFactory& builder) const;
			
		private :
			HudsonMSDemogEvent_ma();

			std::ostringstream& outputDescription(
									std::ostringstream& stm ) const;
			
			std::string migrationMatrixString() const;

			const std::vector< std::vector < double > > newMigMatrix;
	};
	
	/*! \brief A class providing the same interface as 
	 HudsonMSDemogEvent implemented through a pointer.*/
	class HudsonMSDemogEventPtr
	{
		public:
		
		HudsonMSDemogEventPtr(
			boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > e);
		
		void callEvent(const hudson_ms::HudsonMSArgFactory& builder) const;
		
		double getTime() const;
		
		std::string msString() const;
		
		bool operator < (const hudson_ms::HudsonMSDemogEventPtr& rhs) const;
		
		private:
		HudsonMSDemogEventPtr();
		boost::shared_ptr < const hudson_ms::HudsonMSDemogEvent > eventPtr;
		
	};


} // end namespace hudson_ms

#endif


