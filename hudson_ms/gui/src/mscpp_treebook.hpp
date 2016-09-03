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

#ifndef __INC_MSCPP_TREEBOOK_H__
#define __INC_MSCPP_TREEBOOK_H__

#include "mscpp_pop_struct_listCtrl.hpp"
#include "mscpp_migmat_grid.hpp"

#include <wx/wx.h>
#include <wx/treebook.h>


namespace jah {

	
	class MScppTreeBook : public wxTreebook
	{
	public:
		MScppTreeBook(const boost::shared_ptr< hudson_ms::PopulationStructure > _popPtr,
							wxWindow *parent,
							wxWindowID id = wxID_ANY,
							const wxPoint& pos = wxDefaultPosition,
							const wxSize& size = wxDefaultSize,
							long style = wxBK_DEFAULT);
		
		int getPopStructIndex();
		int getMigmatIndex();
		
		void deleteSelectedPop();
		
		void editSelectedPop();
		
		void addPop();
		
		
	private:
	
		enum
		{
			POP_LIST = wxID_HIGHEST,
			MIG_GRID
			
		};
	
		void recreatePagePopStruct();
		
		void recreatePageMigmat();
		
		void createPagePopStruct();
		
		void createPageMigmat();
		
		void recreatePopStruct(wxWindow *window);
		
		void recreateMigmat(wxWindow *window);
		
		void OnPageChanging(wxTreebookEvent& event);
		
		void OnPageChangingFromMigMat(wxTreebookEvent& event);
		
		void OnPageChangingFromPopStruct(wxTreebookEvent& event);
		
		void OnDeletePop(wxListEvent& event);
		
		jah::MScppPopStructListCtrl* m_listCtrl;
		jah::MScppMigmatGrid* m_grid;
		
		static const int popStructIndex = 0;
		static const int migmatIndex = 1;
		
		const boost::shared_ptr< hudson_ms::PopulationStructure > popPtr;
		
		bool autoPageChangeEvent; // indicator for whether we are in automatic page changes
		
		//This class handles events
		DECLARE_EVENT_TABLE()
		
	};
}
#endif
