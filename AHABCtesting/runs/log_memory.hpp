/* See http://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-run-time-in-c
 * Answer from Don Wakefield - thanks Don!
 * 
 */

#ifndef _INC_LOG_MEMORY_H
#define _INC_LOG_MEMORY_H


#include <iostream>
#include <fstream>
#include <string>



void logUsage(std::ostream& os);

void logUsage(std::ostream& os, const std::string& logline);

void logUsage(const std::string& filename, const std::string& logline,
							bool append = true);

#endif
