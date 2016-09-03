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

#include "mscpp_pop_dialogue.hpp"

#include <wx/valnum.h>
#include <wx/statline.h>


#define MYDEBUG

namespace jah {
	
	MScppPopDialogue::MScppPopDialogue(wxWindow *parent,
							const wxWindowID id,
							const wxString& caption,
							const wxPoint& pos,
							const wxSize& size,
							long style )
		:	wxDialog(parent, id, caption, pos, size, style),
			edit(false), m_labelCtrl(NULL),
			label(wxEmptyString), nsam(0), sizeRatio(1.0), growthRate(0.0)
	{
		createControls();
	}
						   
	MScppPopDialogue::MScppPopDialogue(wxString _label,
							int _nsam,
							double _sizeRatio,
							double _growthRate,
							wxWindow *parent,
							const wxWindowID id,
							const wxString& caption,
							const wxPoint& pos,
							const wxSize& size,
							long style )
		:	wxDialog(parent, id, caption, pos, size, style),
			edit(true), m_labelCtrl(NULL),
			label(_label), nsam(_nsam), sizeRatio(_sizeRatio), growthRate(_growthRate)
	{
		#ifdef MYDEBUG
			std::cout << "edit dialogue constructor:" << std::endl;
			std::cout << "label is " << _label << ", nsam = " << _nsam << ", sizeRatio = " << _sizeRatio << ", growthRate = " << _growthRate << std::endl;
		#endif
				
		
		createControls();
	}
	
	wxString MScppPopDialogue::getLabel()
	{
		return label;
	}
	
	int MScppPopDialogue::getNsam()
	{
		return nsam;
	}
	
	double MScppPopDialogue::getSizeRatio()
	{
		return sizeRatio;
	}
	
	double MScppPopDialogue::getGrowthRate()
	{
		return growthRate;
	}
	
	void MScppPopDialogue::createControls()
	{
		#ifdef MYDEBUG
			std::cout << "creating controls:" << std::endl;
		#endif
		int border = 5;
		
		// top levels sizer and inner sizer
		wxBoxSizer* topSizer = new wxBoxSizer(wxVERTICAL);
		
		
		wxBoxSizer* boxSizer = new wxBoxSizer(wxVERTICAL);
		topSizer->Add(boxSizer, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, border);
		
		SetSizer(topSizer);
		
		//instruction
		wxString instructTxt = (edit ? wxT("Edit ") : wxT("Please enter "));
						
		instructTxt += 
			wxT("the label, the number of samples, size ratio and growth rate for this sub-population");
						
		
		wxStaticText* instruction = new wxStaticText( this,
						wxID_STATIC,
						instructTxt,
						wxDefaultPosition,
						wxDefaultSize,
						0);
		
		boxSizer->Add(instruction, 0, wxALIGN_LEFT|wxALL, border);

		#ifdef MYDEBUG
			std::cout << "added instruction:" << std::endl;
		#endif
		

		//spacer
		boxSizer->Add(border, border, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, border);
		
		//sizer for the controls
		int rows = 4;
		int cols = 3;
		wxFlexGridSizer* gridSizer = new wxFlexGridSizer(rows, cols, 0, 0);
		boxSizer->Add(gridSizer, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, border);
		
		#ifdef MYDEBUG
			std::cout << "added gridsizer" << std::endl;
		#endif
		
		//label
		{
			wxStaticText* labelLabel = new wxStaticText( this,
						wxID_STATIC,
						wxT("Label"),
						wxDefaultPosition,
						wxDefaultSize,
						0);
		
			gridSizer->Add(labelLabel, 0, wxALIGN_LEFT|wxALL, border);
			
			m_labelCtrl = new wxTextCtrl( this,
						ID_LABEL,
						wxEmptyString,
						wxDefaultPosition,
						wxDefaultSize,
						0,
						wxTextValidator(wxFILTER_ALPHANUMERIC, &label)
						);
		
			gridSizer->Add(m_labelCtrl, 0, wxALIGN_LEFT|wxALL, border);
			
			wxStaticText* labelHelp = new wxStaticText( this,
						wxID_STATIC,
						wxT("Labels can contain any alpha-numeric characters"),
						wxDefaultPosition,
						wxDefaultSize,
						0);
		
			
			gridSizer->Add(labelHelp, 0, wxALIGN_LEFT|wxALL, border);
			
			#ifdef MYDEBUG
				std::cout << "label controls" << std::endl;
			#endif
		}
		
		//nsam
		{
			wxStaticText* nsamLabel = new wxStaticText( this,
							wxID_STATIC,
							wxT("Nsam"),
							wxDefaultPosition,
							wxDefaultSize,
							0);
			
			gridSizer->Add(nsamLabel, 0, wxALIGN_LEFT|wxALL, border);
			
			//validator that will only take +ve ints
			wxIntegerValidator<int>
						val(&nsam);
			//validator that will only take +ve ints
			int min = 0; 
			val.SetMin(min);
			
			
			wxTextCtrl* nsamCtrl = new wxTextCtrl( this,
						ID_NSAM,
						wxEmptyString,
						wxDefaultPosition,
						wxDefaultSize,
						0,
						val);
		
			
		
			gridSizer->Add(nsamCtrl, 0, wxALIGN_LEFT|wxALL, border);
			
			wxStaticText* nsamHelp = new wxStaticText( this,
						wxID_STATIC,
						wxT("Nsam must be >= 0"),
						wxDefaultPosition,
						wxDefaultSize,
						0);
		
			gridSizer->Add(nsamHelp, 0, wxALIGN_LEFT|wxALL, border);
			
			#ifdef MYDEBUG
				std::cout << "added nsam controls" << std::endl;
			#endif
			
		}
			
		// precision for floating point numbers
		unsigned prec = 5;
		
		//size ratio
		{
			wxStaticText* srLabel = new wxStaticText( this,
							wxID_STATIC,
							wxT("Size ratio"),
							wxDefaultPosition,
							wxDefaultSize,
							0);
			
			gridSizer->Add(srLabel, 0, wxALIGN_LEFT|wxALL, border);
			
			wxFloatingPointValidator<double>
							val(prec, &sizeRatio, wxNUM_VAL_NO_TRAILING_ZEROES);
			// make sure the sr validator can only accept values >0.0
			double min = 0.00001;  // min given the precision we are using
			val.SetMin(min);
						
			wxTextCtrl* srCtrl = new wxTextCtrl( this,
						ID_SR,
						wxEmptyString,
						wxDefaultPosition,
						wxDefaultSize,
						0,
						val);
		
			gridSizer->Add(srCtrl, 0, wxALIGN_LEFT|wxALL, border);
			
			wxStaticText* srHelp = new wxStaticText( this,
						wxID_STATIC,
						wxT("Size ratio must be > 0.0;\nLeading 0's entered will not be shown"),
						wxDefaultPosition,
						wxDefaultSize,
						0);
			
			gridSizer->Add(srHelp, 0, wxALIGN_LEFT|wxALL, border);
		
			#ifdef MYDEBUG
				std::cout << "added sizeRatio controls" << std::endl;
			#endif
			
		}
			
		//growth rate
		{
			wxStaticText* grLabel = new wxStaticText( this,
							wxID_STATIC,
							wxT("Growth rate"),
							wxDefaultPosition,
							wxDefaultSize,
							0);
			
			gridSizer->Add(grLabel, 0, wxALIGN_LEFT|wxALL, border);
			
			wxFloatingPointValidator<double>
				val(prec, &growthRate, wxNUM_VAL_NO_TRAILING_ZEROES);
			// no restrictions
		
			wxTextCtrl* grCtrl = new wxTextCtrl( this,
						ID_GR,
						wxEmptyString,
						wxDefaultPosition,
						wxDefaultSize,
						0,
						val);
		
			gridSizer->Add(grCtrl, 0, wxALIGN_LEFT|wxALL, border);
			
			wxStaticText* grHelp = new wxStaticText( this,
						wxID_STATIC,
						wxT("Growth rate can be positive or negative"),
						wxDefaultPosition,
						wxDefaultSize,
						0);
			
			gridSizer->Add(grHelp, 0, wxALIGN_LEFT|wxALL, border);
			#ifdef MYDEBUG
				std::cout << "added growthRate controls" << std::endl;
			#endif
		}
	
		//dividing line
		
		wxStaticLine* line = new wxStaticLine(this, wxID_STATIC,
								wxDefaultPosition, wxDefaultSize,
								wxLI_HORIZONTAL);
		
		boxSizer->Add(line, 0, wxGROW|wxALL, border);
		
		
		//sizer for buttons
		wxBoxSizer* buttonSizer = new wxBoxSizer(wxHORIZONTAL);
		
		boxSizer->Add(buttonSizer, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, border);
		
		//ok button
		wxButton* okButton = new wxButton(this, wxID_OK, wxT("&OK"),
								wxDefaultPosition, wxDefaultSize, 0);
		buttonSizer->Add(okButton, 0, wxALIGN_CENTER_VERTICAL|wxALL, border);
		
		
		//cancel button
		
		//ok button
		wxButton* cancelButton = new wxButton(this, wxID_CANCEL, wxT("&Cancel"),
								wxDefaultPosition, wxDefaultSize, 0);
		buttonSizer->Add(cancelButton, 0, wxALIGN_CENTER_VERTICAL|wxALL, border);
		
		#ifdef MYDEBUG
			std::cout << "added buttons" << std::endl;
		#endif
		
		topSizer->Fit(this);
	}	


	bool MScppPopDialogue::TransferDataToWindow()
	{
		bool r = wxDialog::TransferDataToWindow();
		// These function calls have to be made here, after the
		// dialog has been created.
		if (m_labelCtrl != NULL) m_labelCtrl->SetFocus();
		return r;
	}

}
