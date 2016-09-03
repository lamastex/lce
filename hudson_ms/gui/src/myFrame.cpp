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

#include "myFrame.hpp"

#define MYDEBUG

namespace jah {

	/*
	 * In the constructor of the main window (or later on)
	 * we create a menu with two menu items as well as a 
	 * status bar to be shown at the bottom of the main 
	 * window. Both have to be "announced" to 
	 * the frame with respective calls.
	 * */
	MyFrame::MyFrame(
		const boost::shared_ptr< hudson_ms::PopulationStructure > _popPtr,
		const wxString& title, 
		const wxPoint& pos, 
		const wxSize& size
		)
		: wxFrame( NULL, -1, title, pos, size ), m_treebook(NULL), 
			popStructIndex(0), migmatIndex(0),
			popMenuTitle(_T("&Populations")), migMenuTitle(_T("&Migration")),
			m_panel(NULL),
			m_menuBar(NULL), m_menuPopStruct(NULL), m_menuMigmat(NULL),
			bookSetupDone (false), popPtr(_popPtr)
		
	{
		m_panel = new wxPanel(this, wxID_ANY);
		
		createBook();
		
		// do the menu after the book
		createMenuBar();
		
		// set menu to match whatever book is showing
		setMenu();
		
		createStatusBar();
		
		//sizer
		m_sizerFrame = new wxBoxSizer(wxVERTICAL);
		m_sizerFrame->Add(m_treebook,
					wxSizerFlags(1).Expand().Border(wxALL, 10)
					);
							// vertically stretchable
							// horizontally stretchable
							// border all around
							// border width 10);
		
		//m_sizerFrame->Layout();
		
		m_panel->SetSizer(m_sizerFrame);
		
		m_sizerFrame->Fit(this);
		m_sizerFrame->SetSizeHints(this);
	}
	
	MyFrame::~MyFrame()
	{
		// only the currently attached menu will be destroyed
		
		// if each menu is not currently on the menubar, delete it
		if (m_menuBar->FindMenu(popMenuTitle) == wxNOT_FOUND) {
			std::cout << "deleting popmenu" << std::endl;
			delete m_menuPopStruct;
		}
		// if each menu is not currently on the menubar, delete it
		if (m_menuBar->FindMenu(migMenuTitle) == wxNOT_FOUND) {
			std::cout << "deleting migmenu" << std::endl;
			delete m_menuMigmat;
		}
		
	}
	
	void MyFrame::createMenuBar()
	{
		m_menuBar = new wxMenuBar;

		wxMenu* menuFile = createMainMenu();
		m_menuBar->Append(menuFile, _("&File") );
		
		//create both menus but only attach the 
		createPopListMenu();
		//menuBar->Append(m_menuPopStruct, _("&Sub-populations") );

		createMigmatMenu();
		//menuBar->Append(m_menuMigmat, _("&Sub-populations") );
	
		SetMenuBar(m_menuBar);
		
		
	}
	
	wxMenu * MyFrame::createMainMenu()
	{
		wxMenu *menuFile = new wxMenu;

		menuFile->Append( wxID_ABOUT, wxT("&About...") );
		menuFile->AppendSeparator();
		menuFile->Append( wxID_EXIT, _("E&xit\tCtrl-x") );
	
		return menuFile;
	}
	
	void MyFrame::createPopListMenu()
	{
		m_menuPopStruct = new wxMenu();
		m_menuPopStruct->Append(LIST_ADD, _T("&Add a sub-pop\tCtrl-a"));
		m_menuPopStruct->Append(LIST_EDIT, _T("&Edit the sub-pop\tCtrl-e"));
		m_menuPopStruct->Append(LIST_DELETE, _T("&Delete the sub-pop\tCtrl-d"));
		
	}	
	
	void MyFrame::createMigmatMenu()
	{
		m_menuMigmat = new wxMenu();
		m_menuMigmat->Append(MIGMAT_RESTORE, _T("&Restore matrix\tCtrl-r"));
		m_menuMigmat->Append(MIGMAT_ALL_ZERO, _T("Set entries to &Zero\tCtrl-z"));
		
	}	
	
	
		
	
	void MyFrame::createBook()
	{
		// book goes in panel
		m_treebook = new jah::MScppTreeBook(popPtr, m_panel);
		
		// no image list at the moment
		
		popStructIndex = m_treebook->getPopStructIndex();
		migmatIndex = m_treebook->getMigmatIndex();
		
		bookSetupDone = true;
	
	}
	
	

	void MyFrame::createStatusBar()
	{
	
		CreateStatusBar();
			
		SetStatusText( _("Welcome to MSCplusplus!") );
	}
	
	// set the right menu for the tree book selection
	void MyFrame::setMenu()
	{
		if (NULL != m_treebook && NULL != m_menuBar) { 
			
			if (m_menuBar->GetMenuCount() > 1) {
				m_menuBar->Remove(1);
				
				assert(m_menuBar->GetMenuCount() == 1);
			}
			int page =  m_treebook->GetSelection();
			
			if (page != wxNOT_FOUND) {
				
				if (page == popStructIndex) {
						m_menuBar->Append(m_menuPopStruct, popMenuTitle);
				}
				if (page == migmatIndex) {
						m_menuBar->Append(m_menuMigmat, migMenuTitle);
				}
			}
		}
		
		m_menuBar->Refresh();
	}
	
	/*
	 * Here are the actual event handlers.
	 */

	/*
	 * MyFrame::OnQuit() closes the main window by 
	 * calling Close(). The paramter TRUE indicates 
	 * that other windows have no veto power such 
	 * as after asking "Do you really want to close?".
	 * If there is no other main window left, 
	 * the applicatin will quit. 
	 */
	void MyFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
	{
		Close(true);
	}

	/* 
	 * MyFrame::OnAbout() will display a 
	 * small window with some text in it. 
	 * In this case a typical "About" 
	 * window with information about the program. 
	 
	 */
	void MyFrame::OnAbout(wxCommandEvent& WXUNUSED(event))
	{
		wxMessageBox( _("wxWidgets Hello World example."), 
					  _("About Hello World"),
					  wxOK|wxICON_INFORMATION, this );
	}
	
	
	void MyFrame::OnPopMenuDelete(wxCommandEvent& event)
	{
		// just tell the treebook to process
		#ifdef MYDEBUG
			std::cout << "Processing pop menu delete item event" << std::endl;
		#endif
		
		m_treebook->deleteSelectedPop();
			
		Layout();
	}
	
	void MyFrame::OnPopMenuAdd(wxCommandEvent& event)
	{
		// just tell the treebook to process
		#ifdef MYDEBUG
			std::cout << "Processing pop menu edit item event" << std::endl;
		#endif
		
		m_treebook->addPop();
		
		Layout();
	}
	
	void MyFrame::OnPopMenuEdit(wxCommandEvent& event)
	{
		// just tell the treebook to process
		#ifdef MYDEBUG
			std::cout << "Processing pop menu edit item event" << std::endl;
		#endif
		
		m_treebook->editSelectedPop();
	}
	
	
	void MyFrame::OnPageChanged(wxTreebookEvent& event)
	{
		const int id = event.GetSelection();
		
		// only process events on user-initiated changes, ie after set up
		if (bookSetupDone && (id != wxNOT_FOUND)) {
			
			#ifdef MYDEBUG
				std::cout << "Page change event in MyFrame: change to " << id << std::endl;
			#endif
					
			setMenu();
		}

	}
	
	//Event table for the frame
	BEGIN_EVENT_TABLE(MyFrame, wxFrame)
		
		// file menu events
		EVT_MENU(wxID_ABOUT, MyFrame::OnAbout)
		EVT_MENU(wxID_EXIT, MyFrame::OnQuit)
		
		EVT_MENU(MyFrame::POP_MENU_ADD, MyFrame::OnPopMenuAdd)
		EVT_MENU(MyFrame::POP_MENU_EDIT, MyFrame::OnPopMenuEdit)
		EVT_MENU(MyFrame::POP_MENU_DELETE, MyFrame::OnPopMenuDelete)
		
		EVT_TREEBOOK_PAGE_CHANGED(wxID_ANY, MyFrame::OnPageChanged)
		
	END_EVENT_TABLE()
	
	

}
