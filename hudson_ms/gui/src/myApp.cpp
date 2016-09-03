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

#include "myApp.hpp"

#include "myFrame.hpp"

namespace jah {
	
	/*
	 * wxApp::OnInit() is called upon startup and should be used 
	 * to initialize the program, maybe showing a "splash screen" 
	 * and creating the main window (or several). The 
	 * frame should get a title bar text ("Hello World") and a 
	 * position and start-up size. One frame can also be declared 
	 * to be the top window. Returning TRUE indicates a 
	 * successful intialization.
	 */
	bool MyApp::OnInit()
	{
		
		setupPop();
				
		/*
		 * When deleting a frame or dialog, use Destroy 
		 * rather than delete so that the wxWidgets delayed 
		 * deletion can take effect. This waits until 
		 * idle time (when all messages have been processed) 
		 * to actually delete the window, to avoid 
		 * problems associated with the GUI sending 
		 * events to deleted windows.
		 * 
		 * Don't create a window on the stack, 
		 * because this will interfere with 
		 * delayed deletion.
		 */
		MyFrame *frame = new MyFrame(popPtr, wxT("MSCplusplus"));
		
		int width, height;
		frame->GetSize(& width, & height);
		std::cout << "frame width is " << width << " << and height is " << height << std::endl;
		
		//frame->SetSize(wxDefaultCoord, wxDefaultCoord, width, height);

		frame->Show(true);
		SetTopWindow(frame);
		return true;
	}
	
	void MyApp::setupPop()
	{
		
		//set up the population structure
		popPtr = boost::shared_ptr< hudson_ms::PopulationStructure >
					( new hudson_ms::PopulationStructure() );
		
		int nsps = 4;
			
		for (int nsam = 1; nsam <= nsps; ++nsam) {
			
			if (nsam < 3) {
				popPtr->add(nsam);
			}	
			else if (nsam == 3) {
				double sizeRatio = 0.5;
				double growthRate = 10.0;
				popPtr->add(nsam, sizeRatio, growthRate);
			}
			else {
				double growthRate = 3.5;
				
				popPtr->add(nsam, growthRate);
			}
		}
		
		double allMigs = 5.1111;
		popPtr->setMigrationMatrix(allMigs);
		
			
	}	
}
