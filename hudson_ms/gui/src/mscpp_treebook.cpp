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

#include "mscpp_treebook.hpp"
#include "mscpp_pop_dialogue.hpp"

#include <stdexcept>

#define MYDEBUG

namespace jah {

	const int MScppTreeBook::popStructIndex;
	const int MScppTreeBook::migmatIndex;
	
	/*
	 * In the constructor of the main window (or later on)
	 * we create a menu with two menu items as well as a 
	 * status bar to be shown at the bottom of the main 
	 * window. Both have to be "announced" to 
	 * the frame with respective calls.
	 * */
	MScppTreeBook::MScppTreeBook(
							const boost::shared_ptr< hudson_ms::PopulationStructure > _popPtr,
							wxWindow *parent,
							wxWindowID id,
							const wxPoint& pos,
							const wxSize& size,
							long style)
        : wxTreebook(parent, id, pos, size, style), 
			m_listCtrl(NULL), m_grid(NULL),
			popPtr(_popPtr),
			autoPageChangeEvent(true)
	{
		
		// create windows to be the pages
		createPagePopStruct();
			
		createPageMigmat();
		
		int width, height;
		GetSize(& width, & height);
		std::cout << "treebook width is " << width << " << and height is " << height << std::endl;
		
		// set the opening page to the population structure one
		// use set because we do want a page change event
		SetSelection(popStructIndex);
		
		autoPageChangeEvent = false;
			
	}
	
	int MScppTreeBook::getPopStructIndex()
	{
		return popStructIndex;
	}
	
	int MScppTreeBook::getMigmatIndex()
	{
		return migmatIndex;
	}
	
	void MScppTreeBook::deleteSelectedPop()
	{
		// check that we are on pop page and there is a pop selected
		if (GetSelection() == popStructIndex) {
			
			if ( m_listCtrl->haveItemSelected()) {
				
				#ifdef MYDEBUG
					int item =  static_cast<int>(m_listCtrl->getSelectedItem());
				
					std::cout << "start of deleting pop, selected item is " << item << std::endl;
				#endif	
				
				if ( wxMessageBox
                     (
                      wxT("Are you sure you want to\n")
                      wxT("delete this sub-population?"),
                      wxT("Confirm deletion"),
                      wxICON_QUESTION | wxYES_NO |wxNO_DEFAULT,
                      this
                     ) == wxYES ) {
                   
					m_listCtrl->deleteSelectedItem();
                }
			}
		
			else {
				wxMessageBox
					 (
					  wxT("You need to select a sub-population to delete\n")
					  wxT(""),
					  wxT("No sub-population selected"),
					  wxICON_ERROR | wxOK,
					  this);
			}
			
		}
		// just ignore if not on popPage
	}
		
	void MScppTreeBook::editSelectedPop()
	{
		#ifdef MYDEBUG
			std::cout << "In editSelectedPop" << std::endl;
		#endif
		
		// check that we are on pop page and there is a pop selected
		if (GetSelection() == popStructIndex) {
			
			if ( m_listCtrl->haveItemSelected()) {
				
				//get the item
				int item =  static_cast<int>(m_listCtrl->getSelectedItem());
				#ifdef MYDEBUG
					std::cout << "selected pop index is " << item << std::endl;
				#endif
				
				bool changed = false;
				
				std::string label = popPtr->getLabelAt(item);
				size_t nsam = popPtr->getSubNsamAt(item);
				double sizeRatio = popPtr->getSizeRatioAt(item);
				double growthRate = popPtr->getGrowthRateAt(item);
				// collect the changes
				
				#ifdef MYDEBUG
					std::cout << "current label is " << label << ", nsam = " << nsam << ", sizeRatio = " << sizeRatio << ", growthRate = " << growthRate << std::endl;
				#endif
				
				MScppPopDialogue editDlg(wxString(label), static_cast<int>(nsam),
							sizeRatio, growthRate, this);
							
				if (editDlg.ShowModal() == wxID_OK) {
					
					std::string _label = (editDlg.getLabel()).ToStdString();
					size_t _nsam = static_cast<size_t>(editDlg.getNsam());
					double _sizeRatio = editDlg.getSizeRatio();
					double _growthRate = editDlg.getGrowthRate();
				
					#ifdef MYDEBUG
						std::cout << "new label is " << _label << ", _nsam = " << _nsam << ", _sizeRatio = " << _sizeRatio << ", _growthRate = " << _growthRate << std::endl;
					#endif
				
					
					try {
					
						//change the values first
						if (_nsam != nsam) {
							#ifdef MYDEBUG
								std::cout << "Changing nsam for label = \"" << label << "\" to " << _nsam << std::endl;
							#endif
							popPtr->setSubNsamAt(item, _nsam);
							changed = true;
						}
						if (_sizeRatio != sizeRatio) {
							#ifdef MYDEBUG
								std::cout << "Changing sizeRatio for label = \"" << label << "\" to " << _sizeRatio << std::endl;
							#endif
							popPtr->setSizeRatioAt(item, _sizeRatio);
							changed = true;
						}
						if (_growthRate != growthRate) {
							#ifdef MYDEBUG
								std::cout << "Changing growthRate for label = \"" << label << "\" to " << _growthRate << std::endl;
							#endif
							popPtr->setGrowthRateAt(item, _growthRate);
							changed = true;
						}
						
						// and then change the labels if necessary
						if (_label != label) {
							// changing label could throw an exception
							try {
								#ifdef MYDEBUG
									std::cout << "Changing label for index = " << item << " to " << _label << std::endl;
								#endif
							
								popPtr->relabel(item, _label);
								changed = true;
							}
							catch (std::invalid_argument& ia) {
								
								wxMessageBox
								 (
								  wxT("Sorry, could not change\n")
								  wxT("label to ") + wxString(_label),
								  wxT("Sub-population relabel failed"),
								  wxICON_ERROR | wxOK,
								  this);
								
							}

						}	
					
						if (changed) {
							#ifdef MYDEBUG
								std::cout << "recreating the page pop" << std::endl;
							#endif
							m_listCtrl->populateItems(popPtr);
						}
					}
					catch (std::invalid_argument& ia) {
						
						std::cout << "Invalid argument:\n" << ia.what() << std::endl;
								
								
								
						wxMessageBox
						 (
						  wxT("Sorry, there was an error processing your edits:\n")
						  wxT("please check your values and try again."),
						  wxT("Sub-population edit failed"),
						  wxICON_ERROR | wxOK,
						  this);
								
					}
					
				} // end ShowModal
				
			}
		
			else {
				wxMessageBox
					 (
					  wxT("You need to select a sub-population to edit\n")
					  wxT(""),
					  wxT("No sub-population selected"),
					  wxICON_ERROR | wxOK,
					  this);
			}
			
			
			// just ignore if not on popPage
			#ifdef MYDEBUG
				int item =  static_cast<int>(m_listCtrl->getSelectedItem());
			
				std::cout << "end of adding pop, selected item is " << item << std::endl;
			#endif	
			
			ChangeSelection(popStructIndex);
		}
	}
	
	void MScppTreeBook::addPop()
	{
		#ifdef MYDEBUG
			std::cout << "In addPop" << std::endl;
		#endif
		
		// check that we are on pop page
		if (GetSelection() == popStructIndex) {
			
			MScppPopDialogue editDlg(this);
						
			if (editDlg.ShowModal() == wxID_OK) {
				
				std::string _label = (editDlg.getLabel()).ToStdString();
				int _nsam = editDlg.getNsam();
				double _sizeRatio = editDlg.getSizeRatio();
				double _growthRate = editDlg.getGrowthRate();
			
				#ifdef MYDEBUG
					std::cout << "new label is " << _label << ", _nsam = " << _nsam << ", _sizeRatio = " << _sizeRatio << ", _growthRate = " << _growthRate << std::endl;
				#endif
			
				
				try {
				
					#ifdef MYDEBUG
							std::cout << "Adding new pop" << std::endl;
					#endif
					int old_nsps = static_cast<int> (popPtr->nSubPops());
					popPtr->add(_nsam, _sizeRatio, _growthRate);
					int nsps = static_cast<int> (popPtr->nSubPops());
					assert(nsps == old_nsps + 1);
					
					//change the default label if the user has supplied one
					if (!_label.empty()) {
						int nsps = static_cast<int> (popPtr->nSubPops());
						popPtr->relabel(nsps-1, _label);
					}
					
					recreatePagePopStruct();
				}
				catch (std::invalid_argument& ia) {
					
					std::cout << "Invalid argument:\n" << ia.what() << std::endl;
							
					wxMessageBox
					 (
					  wxT("Sorry, there was an error adding your new sub-population:\n")
					  wxT("please check your values and try again."),
					  wxT("Add sub-population failed"),
					  wxICON_ERROR | wxOK,
					  this);
							
				}
				
			} // end ShowModal
			
			
			#ifdef MYDEBUG
			int item =  static_cast<int>(m_listCtrl->getSelectedItem());
		
			std::cout << "end of editing pop, selected item is " << item << std::endl;
			#endif	
			
			ChangeSelection(popStructIndex);
		}		
		// just ignore if not on popPage
			
		
	}
	
	
	void MScppTreeBook::recreatePagePopStruct()
	{
		#ifdef MYDEBUG
			std::cout << "In recreatePagePopStruct, about to DeletePage(popStructIndex)" << std::endl;
		#endif
		
		autoPageChangeEvent = true;
		
		DeletePage(popStructIndex);
		
		m_listCtrl = NULL;
		
		#ifdef MYDEBUG
			std::cout << "In recreatePagePopStruct, about to createPagePopStruct()" << std::endl;
		#endif
		
		createPagePopStruct();
		
		autoPageChangeEvent = false;
		
	}
	
	void MScppTreeBook::recreatePageMigmat()
	{
		autoPageChangeEvent = true;
		
		#ifdef MYDEBUG
			std::cout << "In recreatePagePopStruct, about to DeletePage(migmatIndex)" << std::endl;
		#endif
		
		DeletePage(migmatIndex);
		
		m_grid = NULL;
		
		#ifdef MYDEBUG
			std::cout << "In recreatePagePopStruct, about to createPageMigmat()" << std::endl;
		#endif
		
		createPageMigmat();
		
		autoPageChangeEvent = false;
		
	}
	
	void MScppTreeBook::createPagePopStruct()
	{
		wxPanel* window = new wxPanel(this);
		
		/* should have our own list cntrl type 
		 * constructor would be passed the current population */
		/* if we add a new population, the frame will get the dialogue results
		 * and do things to the population, and call recreate list control
		 * and similarly if we use a dialogue owned by the frame to edit an item
		 * but what if we change the order?  then finish of the ordering must trigger
		 * an event for the frame which will give the pop a new order */
		
		recreatePopStruct(window);
		
		{
			int width, height;
			m_listCtrl->GetSize(& width, & height);
			std::cout << "listctrl width is " << width << " << and height is " << height << std::endl;
		}
		
		wxBoxSizer *sizerPanel = new wxBoxSizer(wxVERTICAL);
		
		sizerPanel->Add(m_listCtrl,
					wxSizerFlags(1).Expand().FixedMinSize()
					);
						// vertically stretchable	
						// horizontally stretchable
						// initial size of window is its size
		
		window->SetSizer(sizerPanel);
		sizerPanel->Fit(window);
	
		//sizerPanel->Layout();
		
		InsertPage(popStructIndex, window, wxT("Population Structure"), false);
	}
	
	void MScppTreeBook::createPageMigmat()
	{
		wxPanel* window = new wxPanel(this);
		
		recreateMigmat(window);
		
		{
			int width, height;
			m_grid->GetSize(& width, & height);
			std::cout << "grid width is " << width << " << and height is " << height << std::endl;
		}
		
		wxBoxSizer *sizerPanel = new wxBoxSizer(wxVERTICAL);
		
		sizerPanel->Add(m_grid,
					wxSizerFlags(1).Expand().FixedMinSize()
					);
						// vertically stretchable	
						// horizontally stretchable
						// initial size of window is its size
		
		window->SetSizer(sizerPanel);
		sizerPanel->Fit(window);
	
		//sizerPanel->Layout();
		
		InsertPage(migmatIndex, window, wxT("Migration matrix"), false);
		
		
	}
	
	
		/* should have our own list cntrl type 
		 * constructor would be passed the current population */
		/* if we add a new population, the frame will get the dialogue results
		 * and do things to the population, and call recreate list control
		 * and similarly if we use a dialogue owned by the frame to edit an item
		 * but what if we change the order?  then finish of the ordering must trigger
		 * an event for the frame which will give the pop a new order */
		


	void MScppTreeBook::recreatePopStruct(wxWindow *window)
	{
		#ifdef MYDEBUG
			std::cout << "In recreatePopStruct, about to delete m_listCtrl" << std::endl;
		#endif
		
		delete m_listCtrl;
		
		#ifdef MYDEBUG
			std::cout << " and set new m_listCtrl" << std::endl;
		#endif
		m_listCtrl = new MScppPopStructListCtrl(window, POP_LIST);
		
		m_listCtrl->populateItems(popPtr);
		
		m_listCtrl->resetSize();

	}
	
	void MScppTreeBook::recreateMigmat(wxWindow *window)
	{
		delete m_grid;
		
		m_grid = new MScppMigmatGrid(window, MIG_GRID);
		
		m_grid->populateMatrix(popPtr);
		
		//m_listCtrl->resetSize();

	}
	
	void MScppTreeBook::OnPageChanging(wxTreebookEvent& event)
	{
		const int idx = event.GetOldSelection();
		
		// only process events on user-initiated changes
		if (!autoPageChangeEvent && (idx != wxNOT_FOUND)) {
			
			#ifdef MYDEBUG
				std::cout << "Page change event: change from " << idx << std::endl;
			#endif
					
			if(idx == popStructIndex) {
				OnPageChangingFromPopStruct(event);
				
			}
			if (idx == migmatIndex) {
				
				OnPageChangingFromMigMat(event);
			}
		}

	}
	
	void MScppTreeBook::OnPageChangingFromMigMat(wxTreebookEvent& event)
	{
		#ifdef MYDEBUG
			std::cout << "Page change event: change from MigMat" << std::endl;
		#endif
		
		// we want to get the migration matrix back out
		std::vector < std::vector < double > > migmat = m_grid->getMatrix();
		
		if (!m_grid->getMatrixValidation()) {
			
			event.Veto();
			
			wxMessageBox
					 (
					  wxT("Error in your matrix\n")
					  wxT("Please check and reenter"),
					  wxT("Migration matrix warning"),
					  wxICON_ERROR | wxOK,
					  this);
		}
		
		else // have a valid matrix
		{
			#ifdef MYDEBUG
				std::cout << "setting new migmat" << std::endl;
			#endif
			popPtr->setMigrationMatrix(migmat);
			
			#ifdef MYDEBUG
				std::cout << "migration matrix is now" << std::endl;
				std::cout << popPtr->migrationMatrixToString() << std::endl;
			#endif
			
			//don't update the population structure page
			
		}
	}
	
	void MScppTreeBook::OnPageChangingFromPopStruct(wxTreebookEvent& event)
	{
		#ifdef MYDEBUG
			std::cout << "Page change event: change from PopStruct" << std::endl;
		#endif
		
		#ifdef MYDEBUG
				std::cout << "about to recreatePageMigMat()" << std::endl;
		#endif
		recreatePageMigmat();
			
	}
	
	void MScppTreeBook::OnDeletePop(wxListEvent& event)
	{
		long deleted = event.GetIndex();
		
		#ifdef MYDEBUG
			std::cout << "OnDelete event from the pop list ctrl" << std::endl;
			std::cout << "Index of item to be deleted is " << deleted << std::endl;
		#endif
		
		// delete from the population
		popPtr->erase(deleted);
			
	}
	
	
	//Event table for the frame
	BEGIN_EVENT_TABLE(MScppTreeBook, wxTreebook)
		
		// treebook events
		EVT_TREEBOOK_PAGE_CHANGING(wxID_ANY, MScppTreeBook::OnPageChanging)
		
		// events in the treebook controls
		EVT_LIST_DELETE_ITEM(POP_LIST, MScppTreeBook::OnDeletePop)
		
	END_EVENT_TABLE()
	

}
