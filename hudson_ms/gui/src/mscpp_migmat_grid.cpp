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

#include "mscpp_migmat_grid.hpp"

#include <vector>
#include <cassert>

#define MYDEBUG

namespace jah {
	
	MScppMigmatGrid::MScppMigmatGrid(wxWindow *parent,
						   const wxWindowID id,
						   const wxPoint& pos,
						   const wxSize& size,
						   long style,
						   const wxString &name)
        : wxGrid(parent, id, pos, size, style, name)
	{
		nonEditColour = wxTheColourDatabase->Find(wxT("LIGHT GREY"));
		labelColour = wxTheColourDatabase->Find(wxT("LIGHT BLUE"));
		problemColour = wxTheColourDatabase->Find(wxT("RED"));
		cellTextColour = GetDefaultCellTextColour();
		
		CreateGrid(0,0); // no rows, no columns
		SetDefaultRenderer(new wxGridCellFloatRenderer);
		SetDefaultEditor(new wxGridCellFloatEditor);
		
				
	}
	
	//this is where we would pass in the population
	void MScppMigmatGrid::populateMatrix(
			const boost::shared_ptr< const hudson_ms::PopulationStructure > popPtr)
	{
		
		#ifdef MYDEBUG
				std::cout << "in populateMatrix, pop migration matrix is" << std::endl;
				std::cout << popPtr->migrationMatrixToString() << std::endl;
			#endif

		// to speed up inserting, hide temporarily
		Hide();
		
		size_t popSize = popPtr->nSubPops();
		
		std::vector < std::vector < double > > migmat = popPtr->makeMigration_matrix();
		
		AppendRows(popSize);
		AppendCols(popSize);
		SetLabelBackgroundColour(labelColour);
		
		#ifdef MYDEBUG
			std::cout << "CanDragCell() = " << CanDragCell() << std::endl;
			std::cout << "CanDragColMove() = " << CanDragColMove() << std::endl;
			std::cout << "CanDragColSize(0) = " << CanDragColSize(0) << std::endl;
			std::cout << "CanDragGridSize() = " << CanDragGridSize() << std::endl;
			std::cout << "CanDragRowSize(0) = " << CanDragRowSize(0) << std::endl;
		#endif
		

		for ( size_t i = 0; i < popSize; ++i ) {
			
			{
				wxString lab = wxString(popPtr->getLabelAt(i));
				SetColLabelValue(i, lab);
				SetRowLabelValue(i, lab);
				SetColFormatFloat(i, 6, 1); // width 6, precision 1
			}
			
			for ( size_t j = 0; j < popSize; ++j ) {
				
				if (i != j) {
					wxString buf;
					buf.Printf(_T("%f"), migmat[i][j]);
					SetCellValue(i, j, buf);
				}
				else {
					SetReadOnly(i, j, true);
					SetCellBackgroundColour(i, j, nonEditColour);
				}
			}
				

		}
		
		AutoSize(); //autosize rows and columns
		Fit(); // set grid size to min to show content
		
		int width, height;
		GetSize(& width, & height);
		std::cout << "after autosize, grid width is " << width << " << and height is " << height << std::endl;
		
		setMatrixValidation(true);
		
		Show();
		
	}

	void MScppMigmatGrid::setMatrixValidation(bool setTo)
	{
		isValidMatrix = setTo;
	}
	
	bool MScppMigmatGrid::getMatrixValidation()
	{
		return isValidMatrix;
	}

	std::vector < std::vector < double > > MScppMigmatGrid::getMatrix()
	{
		setMatrixValidation(true);

		size_t rows = static_cast<size_t>(GetNumberRows());
		size_t cols = static_cast<size_t>(GetNumberCols());
		
		assert(rows == cols);
		
		std::vector < std::vector < double > > migmat;
		
		if (rows > 0) {
			
			//make sure any current change is 'written in'
			GoToCell(0,0);
		
			std::vector < double > tmp(cols, 0.0);
			
			while (migmat.size() < rows) {
				migmat.push_back(tmp);
			}
				
			for ( size_t i = 0; i < rows; ++i ) {
				
				for ( size_t j = 0; j < cols; ++j ) {
					
					if (i != j) {
						wxString val = GetCellValue(i, j);
						double m;
						if (val.ToDouble(&m) && m >= 0.0) {
							SetCellTextColour(cellTextColour, i, j);
							migmat[i][j] = m;
						}
						else {
							setMatrixValidation(false);
							SetCellTextColour(problemColour, i, j);
						}
					}
				}
			}	
			
			// repaint grid
			ForceRefresh();
			
			#ifdef MYDEBUG
				for ( size_t i = 0; i < rows; ++i ) {
				
					for ( size_t j = 0; j < cols; ++j ) {
						
						if (GetCellTextColour(i, j) == problemColour) std::cout << "Problem colour at (" << i << ", " << j << ")" << std::endl;
						
					}
					
				}
			#endif	
			
			
		}
		return migmat;
	}
}
