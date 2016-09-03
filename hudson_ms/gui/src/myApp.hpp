/*
 *      GUI application for mscplusplus
 *      
 *      Copyright 2012 Jennifer A Harlow <jah217@mathxeon3>
 *      
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; either version 2 of the License, or
 *      (at your option) any later version.
 *      
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *      
 *      You should have received a copy of the GNU General Public License
 *      along with this program; if not, write to the Free Software
 *      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *      MA 02110-1301, USA.
 */

#ifndef INC_MYAPP_H
#define INC_MYAPP_H

#include <wx/wx.h>

#include <mscplusplus/population_structure.hpp>

#include <boost/smart_ptr.hpp>


namespace jah {

	class MyApp : public wxApp
	{
		public:
			virtual bool OnInit();
		
		private:
		
			void setupPop();
			
			boost::shared_ptr< hudson_ms::PopulationStructure > popPtr;
	};

}

#endif
