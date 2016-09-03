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

#ifndef INC_MYFRAME_H
#define INC_MYFRAME_H

#include "mscpp_treebook.hpp"

#include <mscplusplus/population_structure.hpp>

#include <boost/smart_ptr.hpp>

#include <wx/wx.h>
#include <wx/panel.h>



namespace jah {

	enum
	{
		LIST_ADD = wxID_HIGHEST,
		LIST_EDIT,
		LIST_DELETE,
		
		MIGMAT_RESTORE,
		MIGMAT_ALL_ZERO
		
	};

	class MyFrame : public wxFrame
	{
	public:
		MyFrame(
		const boost::shared_ptr< hudson_ms::PopulationStructure > _popPtr,
		const wxString& title, 
		const wxPoint& pos = wxDefaultPosition, 
		const wxSize& size = wxDefaultSize
		);
		
		~MyFrame();
		
	private:
	
		enum
		{
			POP_MENU_ADD = wxID_HIGHEST,
			POP_MENU_EDIT,
			POP_MENU_DELETE,
			
			MIGMAT_MENU_RESTORE,
			MIGMAT_MENU_ALL_ZERO
			
		};
	
		void createMenuBar();
		
		wxMenu* createMainMenu();
		
		void createPopListMenu();
		
		void createMigmatMenu();
		
		void createBook();
		
		void createStatusBar();
		
		void setMenu();
		
		// event handlers
		void OnQuit(wxCommandEvent& event);
		
		void OnAbout(wxCommandEvent& event);
		
		void OnPopMenuDelete(wxCommandEvent& event);
		void OnPopMenuAdd(wxCommandEvent& event);
		void OnPopMenuEdit(wxCommandEvent& event);
		
		void OnPageChanged(wxTreebookEvent& event);
		
		// data members
		jah::MScppTreeBook* m_treebook;
		
		int popStructIndex;
		int migmatIndex;
		
		wxString popMenuTitle;
		wxString migMenuTitle;
		
		wxBoxSizer* m_sizerFrame;
		
		wxPanel* m_panel;
		
		wxMenuBar* m_menuBar;
		
		wxMenu* m_menuPopStruct;
		
		wxMenu* m_menuMigmat;
		
		
		bool bookSetupDone;
				
		const boost::shared_ptr< hudson_ms::PopulationStructure > popPtr;
		
		//This class handles events
		DECLARE_EVENT_TABLE()
		
	};
}
#endif
