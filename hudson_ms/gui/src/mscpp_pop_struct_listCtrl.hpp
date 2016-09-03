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

#ifndef __INC_MSCPP_POPSTRUCT_LISTCTRL_H__
#define __INC_MSCPP_POPSTRUCT_LISTCTRL_H__

#include <mscplusplus/population_structure.hpp>

#include <boost/smart_ptr.hpp>

#include <wx/wx.h>
#include <wx/listctrl.h>


namespace jah {

	
	class MScppPopStructListCtrl : public wxListCtrl
	{
	public:
		MScppPopStructListCtrl(wxWindow *parent,
						   const wxWindowID id = wxID_ANY,
						   const wxPoint& pos = wxDefaultPosition,
						   const wxSize& size = wxDefaultSize,
						   long style = wxLC_REPORT | wxLC_SINGLE_SEL);
		
		void populateItems(
					const boost::shared_ptr< const hudson_ms::PopulationStructure > popPtr);
		
		void resetSize();	
		
		void deleteSelectedItem();
		
		bool haveItemSelected();
		
		long getSelectedItem();			   
		
	private:
	
		void OnSelected(wxListEvent& event);
		void OnFocused(wxListEvent& event);
	
		long currentItem;
	
		//This class handles events
		DECLARE_EVENT_TABLE()
	
		
	};
}
#endif
