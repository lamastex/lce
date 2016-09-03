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

#include "mscpp_pop_struct_listCtrl.hpp"

#include <iostream>
#include <sstream>

#define MYDEBUG

namespace jah {

	/*
	 * In the constructor of the main window (or later on)
	 * we create a menu with two menu items as well as a 
	 * status bar to be shown at the bottom of the main 
	 * window. Both have to be "announced" to 
	 * the frame with respective calls.
	 * */
	MScppPopStructListCtrl::MScppPopStructListCtrl(wxWindow *parent,
						   const wxWindowID id,
						   const wxPoint& pos,
						   const wxSize& size,
						   long style)
        : wxListCtrl(parent, id, pos, size, style),
		  currentItem(-1)
	{
		wxListItem itemCol;
		itemCol.SetText(_T("Label"));
		itemCol.SetImage(-1);
		InsertColumn(0, itemCol);

		itemCol.SetText(_T("Nsam"));
		//itemCol.SetAlign(wxLIST_FORMAT_CENTRE);
		InsertColumn(1, itemCol);

		itemCol.SetText(_T("Size ratio"));
		//itemCol.SetAlign(wxLIST_FORMAT_RIGHT);
		InsertColumn(2, itemCol);
		
		itemCol.SetText(_T("Growth rate"));
		//itemCol.SetAlign(wxLIST_FORMAT_RIGHT);
		InsertColumn(3, itemCol);
		
		
	}
	
	//this is where we would pass in the population
	void MScppPopStructListCtrl::populateItems(
			const boost::shared_ptr< const hudson_ms::PopulationStructure > popPtr)
	{

		// to speed up inserting, hide temporarily
		Hide();

		DeleteAllItems();

		size_t nsps =popPtr->nSubPops();
		
		for ( size_t i = 0; i < nsps; ++i )
		{
			wxString buf;
			long tmp = InsertItem(i, wxString(popPtr->getLabelAt(i)), 0);
			SetItemData(tmp, i);

			
			//buf.Printf(_T("%u"), popPtr->getSubNsamAt(i));
			{
				std::ostringstream stm;
				stm << popPtr->getSubNsamAt(i);
				
			//buf.Printf(_T("%u"), i);
				SetItem(tmp, 1, wxString(stm.str()));
			}
			{
				std::ostringstream stm;
				stm << popPtr->getSizeRatioAt(i);

			//buf.Printf(_T("%f"), popPtr->getSizeRatioAt(i));
				SetItem(tmp, 2, wxString(stm.str()));
			}
			{
				std::ostringstream stm;
				stm << popPtr->getGrowthRateAt(i);
				//buf.Printf(_T("%f"), popPtr->getGrowthRateAt(i));
				SetItem(tmp, 3, wxString(stm.str()));
			}
		}
		
		SetColumnWidth( 0, wxLIST_AUTOSIZE_USEHEADER );
		SetColumnWidth( 1, wxLIST_AUTOSIZE_USEHEADER );
		SetColumnWidth( 2, wxLIST_AUTOSIZE_USEHEADER );
		SetColumnWidth( 3, wxLIST_AUTOSIZE_USEHEADER );
		
		currentItem = -1;
		
		Show();
		
	}
	
	void MScppPopStructListCtrl::deleteSelectedItem() 
	{
		#ifdef MYDEBUG
			std::cout << "deleteSelectedItem: before deletion current item is " << currentItem << std::endl;	
		#endif
		if (haveItemSelected()) 
		{
			DeleteItem(currentItem);
		}
		#ifdef MYDEBUG
			std::cout << "deleteSelectedItem: after deletion, setting currentItem to -1" << std::endl;	
		#endif
		
		currentItem = -1;
	}
	
	bool MScppPopStructListCtrl::haveItemSelected()
	{
		return (currentItem >= 0 && currentItem < GetItemCount());
	}
	
	long MScppPopStructListCtrl::getSelectedItem()
	{
		return currentItem;
	}
	
	void MScppPopStructListCtrl::resetSize()
	{
		int width, height;
		GetSize(& width, & height);
		std::cout << "before resize, listctrl width is " << width << " << and height is " << height << std::endl;
		
		int reqWidth = 0;
		for (int i = 0; i < 4; i++) reqWidth += (GetColumnWidth(i)+2);
		
		int reqHeight = 0.0;
		int guessHeaderHeight = 15;
		SetSize(wxDefaultCoord, wxDefaultCoord, reqWidth, reqHeight+guessHeaderHeight);
		
		
		
		#ifdef MYDEBUG
			std::cout << "before getItemCount " << GetItemCount() << " << and GetCountPerPage() " << GetCountPerPage() << std::endl;
		#endif
		
		while (GetItemCount() > GetCountPerPage()) {
			
			reqHeight +=5;
		
			SetSize(wxDefaultCoord, wxDefaultCoord, reqWidth, reqHeight+guessHeaderHeight);
		}
		
	}
	
	void MScppPopStructListCtrl::OnSelected(wxListEvent& event)
	{
		currentItem = event.GetIndex();
		
		#ifdef MYDEBUG
			std::cout << "onSelected event: current item is " << currentItem << std::endl;	
		#endif
	}
	
	void MScppPopStructListCtrl::OnFocused(wxListEvent& event)
	{
		#ifdef MYDEBUG
			std::cout << "onFocused event: focused item is " << event.GetIndex() << std::endl;	
		#endif
	}	
	
	//Event table for the frame
	BEGIN_EVENT_TABLE(MScppPopStructListCtrl, wxListCtrl)	
		EVT_LIST_ITEM_SELECTED(wxID_ANY, MScppPopStructListCtrl::OnSelected)
		EVT_LIST_ITEM_FOCUSED(wxID_ANY, MScppPopStructListCtrl::OnFocused)
		
	END_EVENT_TABLE()
}
