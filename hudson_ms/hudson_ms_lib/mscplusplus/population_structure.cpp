/* mscplusplus - ms for c++.

   Copyright (C) 2011, 2012 Jenny Harlow
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/*! \file     
\brief SubPopulationStructure and PopulationStructure definitions.
*/

#include "population_structure.hpp"


#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>

#include <iterator> // for ostream_iterator

#include <algorithm>
#include <cmath>
#include <utility>
#include <stdexcept>

#include <cassert>

//#define MYDEBUG

using namespace hudson_ms;

PopulationStructure::SubPopulationStructure::SubPopulationStructure(
												size_t _sn,
												double _sr,
												double _gr,
												std::string _label,
												bool _defLab)
	: subNsam(_sn), sizeRatio(_sr), growthRate(_gr),
		label(_label), isDefaultLabel(_defLab)
{}

PopulationStructure::SubPopulationStructure::SubPopulationStructure(
												size_t _sn,
												double _gr,
												std::string _label,
												bool _defLab)
	: subNsam(_sn), sizeRatio(1.0), growthRate(_gr),
		label(_label), isDefaultLabel(_defLab)
{}

PopulationStructure::SubPopulationStructure::SubPopulationStructure(
												size_t _sn,
												std::string _label,
												bool _defLab)
: subNsam(_sn), sizeRatio(1.0), growthRate(0.0),
	label(_label), isDefaultLabel(_defLab)
{}
		
PopulationStructure::SubPopulationStructure::SubPopulationStructure(
												std::string _label,
												bool _defLab)
: subNsam(0), sizeRatio(1.0), growthRate(0.0),
	label(_label), isDefaultLabel(_defLab)
	{}

PopulationStructure::SubPopulationStructure& 
		PopulationStructure::SubPopulationStructure::operator=(
				PopulationStructure::SubPopulationStructure tmp)
{
	tmp.swap(*this);
	return *this;
}
		
void PopulationStructure::SubPopulationStructure::swap(
		PopulationStructure::SubPopulationStructure& other)
{
	std::swap(subNsam, other.subNsam);
	
	std::swap(sizeRatio, other.sizeRatio);
			
	std::swap(growthRate, other.growthRate);
	
	std::swap(label, other.label);
	
	std::swap(isDefaultLabel, other.isDefaultLabel);
	
	std::swap(inwardMigration, other.inwardMigration);
}

size_t PopulationStructure::SubPopulationStructure::getSubNsam() const
{
	return subNsam;
}

void PopulationStructure::SubPopulationStructure::setSubNsam(size_t _sn)
{
	subNsam = _sn;
}	

double PopulationStructure::SubPopulationStructure::getSizeRatio() const
{
	return sizeRatio;
}

void PopulationStructure::SubPopulationStructure::setSizeRatio(double _sr)
{
	sizeRatio = _sr;
}	

double PopulationStructure::SubPopulationStructure::getGrowthRate() const
{
	return growthRate;
}

void PopulationStructure::SubPopulationStructure::setGrowthRate(double _gr)
{
	growthRate = _gr;
}		

std::string PopulationStructure::SubPopulationStructure::getLabel() const
{
	return label;
}

// no change in isDefaultLabel
void PopulationStructure::SubPopulationStructure::setLabel(const std::string& _label)
{
	label = _label;
}


bool PopulationStructure::SubPopulationStructure::getIsDefaultLabel() const
{
	return isDefaultLabel;
}

void PopulationStructure::SubPopulationStructure::setIsDefaultLabel(bool newSetting)
{
	isDefaultLabel = newSetting;
}

// return 0.0 if there is no matching fromLabel
double PopulationStructure::SubPopulationStructure::getInwardMig(
									const std::string& fromLabel) const
{
	double ret = 0.0;
	// use find
	migMap::const_iterator it = inwardMigration.find(fromLabel);
	
	if (it != inwardMigration.end() ) { // fromLabel was not there
		ret = it->second; 
	}
	
	return ret;
	
}

// throw exception if new migration rate < 0.0				
// change if fromLabel is there already, otherwise add
void PopulationStructure::SubPopulationStructure::setInwardMig(
									const std::string& fromLabel, double mig)
{
	if (mig < 0.0 ) { 
		throw std::invalid_argument("mig < 0.0");
	}
	
	//try inserting as a new pair in the map if mig > 0.0
	if (mig > 0.0) {
		std:: pair< migMap::iterator, bool> ret
			= inwardMigration.insert (
					std::pair< std::string, double >(fromLabel, mig) ); 
	
		if (!(ret.second)) { // fromLabel was already a key in the map
			
			ret.first->second = mig;
		}
	}
	else { // mig == 0.0 so take it out if fromLabel is there
		eraseInwardMig(fromLabel);
	}
}

// do nothing if there is no matching fromLabel
void PopulationStructure::SubPopulationStructure::eraseInwardMig(
									const std::string& fromLabel)
{
	// use find
	migMap::iterator it = inwardMigration.find(fromLabel);
	
	if (it != inwardMigration.end() ) { // fromLabel was there
		
		inwardMigration.erase(it);
	} 
}


/* do nothing if oldFromLabel not there
 * else throw an exception if there is already a newFromLabel
 * and replace oldFromLabel with newFromLabel if there is not already a newFromLabel
 */
void PopulationStructure::SubPopulationStructure::relabelInwardMig(
						const std::string& oldFromLabel,
						const std::string& newFromLabel)
{
	// use find
	migMap::iterator it = inwardMigration.find(oldFromLabel);
	
	if (it != inwardMigration.end() ) { // oldFromLabel was there
	
		//check the new label
		if ( inwardMigration.count(newFromLabel) ) {
			throw std::invalid_argument(newFromLabel);
		}
	
		double storeMig = it->second;
		
		inwardMigration.erase(it);
		
		setInwardMig(newFromLabel, storeMig);
		
	}
}


std::string PopulationStructure::SubPopulationStructure::toString() const
{
	// label, size, size ratio, growth rate
	// \t migStuff
	
	std::ostringstream stm;
	
	stm << "\"" << label << "\" : ";
	stm << "\tnsam = " << subNsam;
	stm << "\tsize ratio = " << sizeRatio;
	stm << "\tgrowth rate = " << growthRate;
	
	for (migMap::const_iterator it = inwardMigration.begin();
			it != inwardMigration.end();
			++it) {
		
		stm << "\n\t<- from " << it->first << " = " << it->second;
				
	}
	return stm.str();
	
}

// ------------- PopulationStructure definitions
	

PopulationStructure::PopulationStructure() {}

PopulationStructure::PopulationStructure(int _sn)
{
	if (_sn < 0) {
		throw std::invalid_argument(
		"PopulationStructure::PopulationStructure(int) : _sn < 0)");
	}
	else add(_sn);
	
}
			
PopulationStructure::PopulationStructure(const PopulationStructure& other)
	: container(other.container), eList(other.eList) {}

PopulationStructure::~PopulationStructure() {}


PopulationStructure& PopulationStructure::operator=(PopulationStructure tmp)
{
	tmp.swap(*this);
	return *this;
	
}

void PopulationStructure::swap(PopulationStructure& other)
{
	std::swap(container, other.container);
	
	std::swap(eList, other.eList);
}
			
size_t PopulationStructure::getSubNsamAt(int index) const
{
	try {
		if ( index < 0 || static_cast<size_t>(index) >= nSubPops() ) {
			throw std::invalid_argument("index");
		}
		
		return (container.at(index)).getSubNsam();
		
	}
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::getSubNsamAt(int)");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
	
}

size_t PopulationStructure::getSubNsamAt(const std::string& label) const
{
	try {
		size_t order = getOrder(label);
		
		return (container.at(order)).getSubNsam();
	}
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::getSubNsamAt(const std::string&)");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
}



			
double PopulationStructure::getSizeRatioAt(int index) const
{
	try {
		if ( index < 0 || static_cast<size_t>(index) >= nSubPops() ) {
			throw std::invalid_argument("index");
		}
		
		return (container.at(index)).getSizeRatio();
		
	}
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::getSizeRatioAt(int)");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
	
}

double PopulationStructure::getSizeRatioAt(const std::string& label) const
{
	try {
		size_t order = getOrder(label);
		
		return (container.at(order)).getSizeRatio();
		
	}
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::getSizeRatioAt(int)");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
	
}


			
double PopulationStructure::getGrowthRateAt(int index) const
{
	try {
		if ( index < 0 || static_cast<size_t>(index) >= nSubPops() ) {
			throw std::invalid_argument("index");
		}
		
		return (container.at(index)).getGrowthRate();
		
	}
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::getGrowthRateAt(int)");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
	
}

double PopulationStructure::getGrowthRateAt(const std::string& label) const
{
	try {
		size_t order = getOrder(label);
		
		return (container.at(order)).getGrowthRate();
		
	}
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::getGrowthRateAt(int)");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
	
}


std::vector < std::string > PopulationStructure::getLabels() const
{
	std::vector < std::string > labels;
	
	for (const_iterator it = begin(); it < end(); ++it) {
		labels.push_back(it->getLabel());
	}
	
	return labels;
}

std::string PopulationStructure::getLabelAt(int index) const
{
	try {
		if ( index < 0 || static_cast<size_t>(index) >= nSubPops() ) {
			throw std::invalid_argument("index");
		}
		
		return (container.at(index)).getLabel();
		
	}
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::getLabelAt(int)");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
	
}

bool PopulationStructure::hasLabel(const std::string& label) const
{
	try {
		getOrder(label);
		
		return true;
		
	}
	catch (std::invalid_argument& ia) {
		return false;
	}
	
}

// adding to the population
			
void PopulationStructure::add(int _sn)
{
	if (_sn < 0) {
		throw std::invalid_argument(
		"PopulationStructure::add(...) : _sn < 0)");
	}
	
	std::vector < PopulationStructure::SubPopulationStructure > 
								tmp = container;
	
	std::string newLabel = getNewDefaultLabel(nSubPops());
	bool isDefaultLabel = true;

	size_t sn = static_cast<size_t>(_sn);
	tmp.push_back( PopulationStructure::SubPopulationStructure(sn,
												newLabel, isDefaultLabel) );
	
	container.swap(tmp);
}

void PopulationStructure::add(	int _sn, double _sr, double _gr)
{
	if (_sn < 0) {
		throw std::invalid_argument(
		"PopulationStructure::add(...) : _sn < 0)");
	}
	
	if (_sr <= 0.0) {
		throw std::invalid_argument(
			"PopulationStructure::add(...) : _sr <= 0.0)" );
	}
	
	std::vector < PopulationStructure::SubPopulationStructure > 
								tmp = container;
	
	std::string newLabel = getNewDefaultLabel(nSubPops());
	bool isDefaultLabel = true;
	
	size_t sn = static_cast<size_t>(_sn);
	tmp.push_back( PopulationStructure::SubPopulationStructure(sn, _sr, _gr,
												newLabel, isDefaultLabel) );
	container.swap(tmp);
}

void PopulationStructure::add(int _sn, double _gr)
{
	if (_sn < 0) {
		throw std::invalid_argument(
		"PopulationStructure::add(...) : _sn < 0)");
	}
	
	std::vector < PopulationStructure::SubPopulationStructure > 
								tmp = container;
	
	std::string newLabel = getNewDefaultLabel(nSubPops());
	bool isDefaultLabel = true;
	
	size_t sn = static_cast<size_t>(_sn);
	tmp.push_back( PopulationStructure::SubPopulationStructure(sn, _gr,
												newLabel, isDefaultLabel) );
	
	container.swap(tmp);
}

// erase a subpopulation
// must have eList empty
// vectors are not ideal for erasing but lets assume we don't do this often
void PopulationStructure::erase(const std::string& eraseLabel)
{
	try {
		
		if (!eList.empty()) {
			throw std::runtime_error(
			std::string(
			"PopulationStructure::erase(const std::string&)")
			+ std::string(" : Event list is not empty") );
		}
		
		size_t orderErase = getOrder(eraseLabel);
		
		completeErase(orderErase);
		
	}
	
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::erase(const std::string&)");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
}

void PopulationStructure::erase(int orderErase)
{
	try {
		if (!eList.empty()) {
			throw std::runtime_error(
			std::string(
			"PopulationStructure::erase(int)")
			+ std::string(" : Event list is not empty") );
		}
		
		if ( (orderErase < 0) || (static_cast<size_t>(orderErase) >= nSubPops()) ) {
			throw std::invalid_argument("orderErase");
		}
		
		completeErase(orderErase);
	}
	
	catch (std::invalid_argument& ia) {
		std::string msg(
				"PopulationStructure::erase(int)");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
}


// can only reorder if the event list is empty
// new ordering must be a permutation of the current set of labels
// put the stuff in the container in the order indicated by the the vector of labels
void PopulationStructure::resetOrdering(
		const std::vector < std::string >& newLabelOrder)
{
	try {
		
		if (!eList.empty()) {
			throw std::runtime_error(
			std::string(
			"PopulationStructure::resetOrdering(const std::vector < std::string >&)")
			+ std::string(" : Event list is not empty") );
		}
		
		std::set < std::string > newLabelSet(newLabelOrder.begin(),
											newLabelOrder.end());
		
		size_t nsps = nSubPops();
		
		if ( newLabelSet.size() != nsps ) {
			throw std::invalid_argument("");
		}
		
		// every current label should be in there somewhere
		for (const_iterator cit = begin(); cit < end(); ++cit) {
			if (!(newLabelSet.count(cit->getLabel())) ) {
				throw std::invalid_argument("");
			}
		} 
		
		std::vector < PopulationStructure::SubPopulationStructure > tmp;
		
		for (size_t newOrder = 0; newOrder < nsps; ++newOrder) {
			
			// current order of the one at this pos in the new label order
			size_t oldOrder = 	getOrder( newLabelOrder[newOrder] );
				
			tmp.push_back( container[oldOrder] );
			
		}
	
		assert(tmp.size() == nsps);
		
		// no change in the actual names so no change to subpops migration maps
		
		/* but also reset any default labels 
		 * do this once everything is in tmp so that we can change
		 * labels in inWardsMigration maps as well */
		for (size_t i = 0; i < nsps; i++) {
			
			if (tmp[i].isDefaultLabel) {
				
				std::string oldDefaultLabel = tmp[i].getLabel();
				
				std::string newDefaultLabel = getNewDefaultLabel ( i );
				
				if (newDefaultLabel != oldDefaultLabel) {
					tmp[i].setLabel(newDefaultLabel); // still a default label
					
					for (size_t j = 0; j < nsps; j++) {
						if (i != j) {
							tmp[j].relabelInwardMig(oldDefaultLabel, newDefaultLabel);
						}
					}
				}
			}
		}
				
		//finally swap tmp into container		
		tmp.swap(container);

	}
	catch (std::invalid_argument& ia) {
		
		std::string msg(
		"PopulationStructure::resetOrdering(const std::vector < std::string >&)");
		std::string addMsg = std::string(ia.what());
		
		if (!addMsg.empty()) {
			msg += std::string(" : ") + addMsg;
		}
		throw std::invalid_argument(msg);
	}
	

}



void PopulationStructure::relabel(const std::string& oldLabel,
									const std::string& newLabel)
{
	try {
				
		size_t order = getOrder(oldLabel);
		
		if (!checkUniqueLabel(order, newLabel)) {
			throw std::invalid_argument(newLabel);
		}
		
		completeRelabel(order, newLabel);
		
	}
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::relabel(const std::string&, const std::string&)");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);

	}
}


void PopulationStructure::relabel(int order,
									const std::string& newLabel)
{
	try {
		
		if ((order < 0) || (static_cast<size_t>(order) >= nSubPops())) {
			throw std::invalid_argument("order");
		}
		
		if (!checkUniqueLabel(order, newLabel)) {
			throw std::invalid_argument(newLabel);
		}
		
		completeRelabel(order, newLabel);
		
	}
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::relabel(int, const std::string&)");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
}


void PopulationStructure::relabel(
						const std::vector < std::string >& newLabels)
{
	try {
		
		size_t nsps = nSubPops();
		if (newLabels.size() != nsps) {
			throw std::invalid_argument("newLabels.size()");
		}
		
		std::set < std::string > labelSet(newLabels.begin(), newLabels.end());
		
		if (labelSet.size() != nsps ) {
			throw std::invalid_argument("newLabels");
		}
		
		for (size_t i = 0; i < nsps; ++i) {
			
			completeRelabel(i, newLabels[i]);
		}
	}
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::relabel(const std::vector < std::string >&)");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
}

void PopulationStructure::setSubNsamAt(int index, 
											int _sn)
{
	try {
		if ( index < 0 || static_cast<size_t>(index) >= nSubPops() ) {
			throw std::invalid_argument("into");
		}
		
		if ( _sn < 0 ) {
			throw std::invalid_argument("_sn < 0");
		}
		
		container.at(index).setSubNsam(_sn);
		
	}
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::setSubNsamAt(int, double)");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
	
}

void PopulationStructure::setSubNsamAt(const std::string& label, 
											int _sn)
{
	try {
		size_t order = getOrder(label);
		
		if ( _sn < 0 ) {
			throw std::invalid_argument("_sn < 0");
		}
		
		container.at(order).setSubNsam(_sn);
	}
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::setSubNsamAt(const std::string&, double)");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
}

void PopulationStructure::setSizeRatioAt(int index, 
											double _sr)
{
	try {
		if ( index < 0 || static_cast<size_t>(index) >= nSubPops() ) {
			throw std::invalid_argument("into");
		}
		
		if ( _sr <= 0.0 ) {
			throw std::invalid_argument("_sr <= 0.0");
		}
		
		container.at(index).setSizeRatio(_sr);
		
	}
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::setSizeRatioAt(int, double)");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
	
}

void PopulationStructure::setSizeRatioAt(const std::string& label, 
											double _sr)
{
	try {
		size_t order = getOrder(label);
		
		if ( _sr <= 0.0 ) {
			throw std::invalid_argument("_sr <= 0.0");
		}
		
		container.at(order).setSizeRatio(_sr);
	}
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::setSizeRatioAt(const std::string&, double)");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
}

			
void PopulationStructure::setGrowthRateAt(int index, 
											double _gr)
{
	try {
		if ( index < 0 || static_cast<size_t>(index) >= nSubPops() ) {
			throw std::invalid_argument("into");
		}
		
		container.at(index).setGrowthRate(_gr);
		
	}
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::setGrowthRateAt(int, double)");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
	
}

void PopulationStructure::setGrowthRateAt(const std::string& label, 
											double _gr)
{
	try {
		size_t order = getOrder(label);
		container.at(order).setGrowthRate(_gr);
	}
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::setGrowthRateAt(const std::string&, double)");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
}

void PopulationStructure::setAllGrowthRates(double _gr)
{
	for (std::vector< PopulationStructure::SubPopulationStructure >::iterator it
			= container.begin();
			it < container.end(); ++it ) {
		it->setGrowthRate(_gr);
	}
}

void PopulationStructure::setMigrationMatrix(
					const std::vector < std::vector <double> >& _mm)
{
	try {
	
		size_t nsps = nSubPops();
		
		if (_mm.size() != nsps) {
			throw std::invalid_argument(
				"_mm.size()");
		}
		
		for (std::vector < std::vector <double> >::const_iterator it = _mm.begin();
			it < _mm.end();
			++it) {
				
			if ( it->size() != nsps ) {
				throw std::invalid_argument(
					"_mm[row].size()");
			}
			if ( ( *std::min_element(it->begin(), it->end()) ) < 0.0 ) {
				throw std::invalid_argument(
					"some element < 0.0");
			}
		}
		
		/* Each inner vector 'row' i gives the inwards migration for subpop at index i
		 * For each of these we have to find the labels for the given elements
		 * and then tell the subpop at i to set these inward migration parameters*/
		
		// save time by storing names in a vector rather than repeating lookups
		std::vector < std::string > subPopLabels = getLabels();
				
		for (size_t i = 0; i < nsps; ++i) {
			
			for (size_t j = 0; j < nsps; ++j) {
				// only set for nondiagonals and if mig param is > 0.0
				if ( (i != j) && _mm[i][j] > 0.0) {
					container[i].setInwardMig( subPopLabels[j], _mm[i][j] );
				}
			}
		}
		// subpops are storing the migration parameters
	}
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::setMigrationMatrix(const std::vector < std::vector <double> >&)");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
}


void PopulationStructure::setMigrationMatrix(double _m)
{
	try {
	
		
		if (_m < 0.0) {
			throw std::invalid_argument(
				"_m < 0.0");
		}
		
		if (_m > 0.0) {  // only need to set if mig param is > 0.0 
			
			size_t nsps = nSubPops();
			
			// save time by storing names in a vector rather than repeating lookups
			std::vector < std::string > subPopLabels = getLabels();
						
			for (size_t i = 0; i < nsps; ++i) {
				
				for (size_t j = 0; j < nsps; ++j) {
					// only set for nondiagonals
					if ( (i != j) ) {
						container[i].setInwardMig( subPopLabels[j], _m );
					}
				}
			}
			// subpops are storing the migration parameters
		}
	}
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::setMigrationMatrix(double)");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
}

// throws exception if into or from outof range, or if _m < 0.0
void PopulationStructure::setMigration(int into, 
											int from,
											double _m)
{
	try {
		if (_m < 0.0) {
			throw std::invalid_argument("_m < 0.0");
		}
		
		if ( into < 0 || static_cast<size_t>(into) >= nSubPops() ) {
			throw std::invalid_argument("into");
		}
		
		if ( from < 0 || static_cast<size_t>(from) >= nSubPops() ) {
			throw std::invalid_argument("from");
		}
		
		container[into].setInwardMig( container[from].getLabel(), _m );
	}
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::setMigration(int, int, double");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
}

// throws exception if toLabel or fromLabel not valid, or if _m < 0.0
void PopulationStructure::setMigration(const std::string& toLabel, 
										const std::string& fromLabel,
										double _m)
{
	try {
		if (_m < 0.0) {
			throw std::invalid_argument("_m < 0.0");
		}
		
		size_t into = getOrder(toLabel);
		getOrder(fromLabel); // just to check label exists
		
		container[into].setInwardMig( fromLabel, _m );
		
	}
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::setMigration(const std::string&, const std::string&, double");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
}


double PopulationStructure::getMigration(int into, int from) const
{
	try {
		if ( into < 0 || static_cast<size_t>(into) >= nSubPops() ) {
			throw std::invalid_argument("into");
		}
		
		if ( from < 0 || static_cast<size_t>(from) >= nSubPops() ) {
			
			throw std::invalid_argument("from");
		}
		
		return ( container[into].getInwardMig( container[from].getLabel() ) );
	}
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::getMigration(int, int)");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
	
}


double PopulationStructure::getMigration(const std::string& toLabel, 
									const std::string& fromLabel) const
{
	try {
		
		size_t into = getOrder(toLabel);
		getOrder(fromLabel); // just to check label exists
		
		return ( container[into].getInwardMig( fromLabel ) );
	}
	catch (std::invalid_argument& ia) {
		std::string msg(
			"PopulationStructure::getMigration(const std::string&, const std::string&");
		msg = addToErrorMsg(msg, ia);
		throw std::invalid_argument(msg);
	}
}
			
size_t PopulationStructure::nSubPops() const
{
	return container.size();
		
}
			
bool PopulationStructure::empty() const
{
	return container.empty();
}


bool PopulationStructure::eListEmpty() const
{
	return eList.empty();
}

			
size_t PopulationStructure::totalNsam() const
{
	size_t result = 0;
	
	boost::function< size_t  
				(SubPopulationStructure) >
					f1 = boost::bind(&SubPopulationStructure::getSubNsam, _1);
		
	boost::lambda::placeholder1_type _p;
	
	std::for_each(container.begin(), container.end(), 
		result += boost::lambda::bind(f1, _p) );
		
	return result;
}


void PopulationStructure::setEventList(
	const std::vector< hudson_ms::HudsonMSDemogEventPtr >& _eList) const
{
	eList = _eList;
	// sort but keep relative ordering for equal events
	std::stable_sort( eList.begin(), eList.end() );
}
			
void PopulationStructure::addToEventList(const hudson_ms::HudsonMSDemogEventPtr& e) const
{
	eList.push_back(e);
	std::stable_sort( eList.begin(), eList.end() );
}			


std::vector < int > PopulationStructure::makePopnsams() const 
{
	std::vector < int > popnsams;
	for ( const_iterator it = begin(); it < end(); ++it ) {
		popnsams.push_back( it->getSubNsam() );
	}
	return popnsams;
}
			
std::vector < double > PopulationStructure::makePopsizeratios() const 
{
	std::vector < double > sizeratios;
	for ( const_iterator it = begin(); it < end(); ++it ) {
		sizeratios.push_back( it->getSizeRatio() );
	}
	return sizeratios;
}
			
std::vector < double > 
	PopulationStructure::makePopgrowthrates() const 
{
	std::vector < double > growthrates;
	for ( const_iterator it = begin(); it < end(); ++it ) {
		growthrates.push_back( it->getGrowthRate() );
	}
	return growthrates;
}

// actual pop growthrates are ignored and result just uses \a growthRate
std::vector < double > 
	PopulationStructure::makePopgrowthrates(double growthRate) const 
{
	return	std::vector < double > (nSubPops(), growthRate);
}
			
std::vector< std::vector < double > > 
						PopulationStructure::makeMigration_matrix() const 
{
	return makeMigMatrixFromSubPops();
}

// return the sorted eList;
std::vector< hudson_ms::HudsonMSDemogEventPtr > PopulationStructure::makeEventList() const
{
	return eList;
}			


std::string PopulationStructure::migrationMatrixToString() const
{
	std::vector < std::vector < double > > migration_matrix 
										= makeMigMatrixFromSubPops();
	
	int prec= 1; // use 1 dp precision
	std::ostringstream stm;
	stm.precision(prec);
	stm <<   std::fixed;
	
	if (!migration_matrix.empty()) {
		stm << "SubPop";
	
	
		for (size_t i = 0; i < migration_matrix.size(); ++i) {
			stm << "\t" << i;
		}
		
		for (size_t i = 0; i < migration_matrix.size(); ++i) {
		
			stm << "\n" << i ;
			
			for (size_t j = 0; j < (migration_matrix[i]).size(); ++j) {
			
				if (i != j) {
					stm << "\t" << (migration_matrix[i])[j];
				}
				else {
					stm << "\tn/a"; // don't print diagonals
				}
			}
		}
		stm << "\n";
	
	}
	
	return stm.str();
}

std::string PopulationStructure::makeMScommand(
								int seed, bool time) const 
{
	hudson_ms::HudsonMSRecombParams recombParams;
	hudson_ms::HudsonMSSegSiteParams ssParams;
	
	return makeMScommand(recombParams, ssParams, seed, time);
}



std::string PopulationStructure::makeMScommand(
								const hudson_ms::HudsonMSRecombParams& recombParams,
								int seed, bool time) const 
{
	hudson_ms::HudsonMSSegSiteParams ssParams;
	
	return makeMScommand(recombParams, ssParams, seed, time);
}

std::string PopulationStructure::makeMScommand(
					const hudson_ms::HudsonMSRecombParams& recombParams,
					const hudson_ms::HudsonMSSegSiteParams& ssParams,
					int seed, bool time) const 
{
	std::ostringstream stm;
	stm << std::fixed;
	size_t reps = 1;

	outputStart(stm, reps);
	
	stm << ssParams.msString();
	
	stm << recombParams.msString();
	
	outputSubpops(stm);
	
	stm << eventStringForMS(); 
	if (time) stm << " -L";
	stm << " -T -seeds " << seed;
	
	return stm.str();
}


std::string PopulationStructure::toString() const 
{
	std::ostringstream stm;
	stm << "nsam " << totalNsam();
	
		stm << " npops " << nSubPops(); // 
	
	if (nSubPops() > 0) {	
		for (const_iterator it = begin(); it < end(); ++it) {
			stm << "\n" << it->toString(); 
		
	}
		
	// add on the migration matrix
	stm << "\nM " << migrationMatrixStringForMS();
	stm << "\nEvents" << eventStringForMS(); 
	}
			
	return stm.str();
}



std::vector < std::string > PopulationStructure::makePopulationLabels() const
{
	std::vector < std::string > result;
	
	size_t p = 0; 
	for (PopulationStructure::const_iterator it = begin();
		it < end();
		++it) {
			std::string label = it->getLabel();
			for (size_t i = 0; i < getSubNsamAt(p); ++i) {
				result.push_back(label);
			}
			++p;
	}
	return result;
}


PopulationStructure::const_iterator PopulationStructure::begin() const
{
	return container.begin();
}

PopulationStructure::const_iterator PopulationStructure::end() const
{
	return container.end();
}



std::string PopulationStructure::migrationMatrixStringForMS() const
{
	std::vector < std::vector < double > > migration_matrix 
										= makeMigMatrixFromSubPops();
	
	std::ostringstream stm;
		
	if (!migration_matrix.empty()) {
		
		int prec = 1;
		
		// rough and ready way to find an appropriate precision for printing mig value for ms
		for (size_t i = 0; i < migration_matrix.size(); ++i) {
			for (size_t j = 0; j < (migration_matrix[i]).size(); ++j) {
			
				int altPrec = precEst( migration_matrix[i][j] );
				if (altPrec > prec) prec = altPrec;
			}
		}

		stm.precision(prec);
		stm << std::fixed;
		
		for (size_t i = 0; i < migration_matrix.size(); ++i) {
		
			for (size_t j = 0; j < (migration_matrix[i]).size(); ++j) {
			
				if (i != j) {
					stm << " " << (migration_matrix[i])[j];
				}
				else {
					stm << " x"; // x for diagonals
				}
			}
		}
	}
	
	return stm.str();
}

std::string PopulationStructure::eventStringForMS() const
{
	
	std::string result("");

	for (std::vector< hudson_ms::HudsonMSDemogEventPtr >::const_iterator it = eList.begin();
			it < eList.end(); ++it) {
		result += " " + (it->msString());
	}
	
	return result;
}



std::string PopulationStructure::doMigCheck() const
{
	typedef std::set < size_t > places;
	typedef places::iterator places_it;
	typedef std::map < size_t, places > routemap;
	typedef routemap::iterator routemap_it;
	typedef std::map < places, places > miscommap;
	typedef miscommap::iterator miscommap_it;
	 
	routemap routes; 
	places starts; 
	
	routemap miscommunications;
	
	std::vector < std::vector < double > > migration_matrix 
										= makeMigMatrixFromSubPops();
	
	
	for (size_t i =0; i < nSubPops(); ++i) {
		//send each i off with it's own copy of been
		
		#ifdef MYDEBUG
			std::cout << "\nPop " << i << std::endl;
		#endif
		
		starts.insert(i);
		
		places been;
		been = traceMig(been, migration_matrix, i); // all the places we can get to from i
		
		#ifdef MYDEBUG
		{
			std::cout << "Pop " << i << " can get to: ";
			for (places_it pit = been.begin(); pit != been.end(); ++pit) {
				std::cout << (*pit) << '\t';
			}
			std::cout << std::endl;
		}
		#endif
		
		routes.insert( std::pair < size_t, std::set < size_t > > 
				( i, been ) );
	}
	
	// go through the map
	for (size_t i =0; i < nSubPops(); ++i) {
		routemap_it mit = routes.find(i);
		
		bool okay = true;
		places cant_get_to;
			
		if (mit != routes.end() ) {
			//get any places you can't get to from i
			places can_get_to = mit->second;
			for (places_it sit = starts.begin(); sit != starts.end(); ++sit ) {
				places_it found_it = can_get_to.find(*sit);
				if (found_it == can_get_to.end() ) { // not there
					cant_get_to.insert(*sit);
				}
			} 
		}
		
		#ifdef MYDEBUG
		{
			std::cout << "Pop " << i << " cannot get to: ";
			for (places_it pit = cant_get_to.begin(); pit != cant_get_to.end(); ++pit) {
				std::cout << (*pit) << '\t';
			}
			std::cout << std::endl;
		}	
		#endif

		// for each place we can't get to from i, can that place get to i?
		for (places_it sit = cant_get_to.begin(); sit != cant_get_to.end(); ++sit ) {
			#ifdef MYDEBUG
			{
			std::cout << "\tCan " << (*sit) << " get to " << i << std::endl;
		}	
			#endif
			routemap_it cgtit = routes.find(*sit);
			if ( cgtit != routes.end() ) {
				places can_get_to = cgtit->second;
				places_it found_it = can_get_to.find(i);
				okay = ( found_it != can_get_to.end() );
			}
			if ( !okay ) {
				#ifdef MYDEBUG
				{
				std::cout << "\tNo" << std::endl;
				}	
				#endif
				//break; // break out of loop if can't get to i
				
				routemap_it misit_i = miscommunications.find(i);
				if (misit_i != miscommunications.end() ) {
					misit_i->second.insert(*sit);
				}
				else {
					routemap_it misit_j = miscommunications.find(*sit);
					if (misit_j != miscommunications.end() ) {
						misit_j->second.insert(i);
					}
					else {
						places p;
						p.insert(i);
						miscommunications.insert( 
							std:: pair < size_t, places >(*sit, p));
							
					}
				}
				
			}
		} // end loop through places we can't get to from i
		
	}// end loop through all places
	
	miscommap summary;
	
	for (routemap_it rit = miscommunications.begin();
			rit != miscommunications.end();
			++rit ) {
		miscommap_it it = summary.find( (rit->second) );
		if ( it != summary.end() ) {
			(it->second).insert(rit->first);
		}
		else {
			places p;
			p.insert(rit->first);
			summary.insert( std::pair < places, places > 
								( rit->second, p) );
		}
	}
	
	std::string result("");
	
	if ( !(summary.empty()) ) {
		std::ostringstream stm;
		
		for (miscommap_it mit = summary.begin();
				mit != summary.end();
				++mit ) {
			places from = mit->first;	
			places to = mit->second;		
			stm << "Subpops ("; 
				for (places_it pit = from.begin();
				pit != from.end();
				++pit ) {
					stm << "\t" << (*pit);
				}
			stm << "\t) do not communicate with subpops (";
				for (places_it pit = to.begin();
				pit != to.end();
				++pit ) {
					stm << "\t" << (*pit);
				}
				stm << "\t)" << std::endl;
		}
		
		result += (std::string("Failed PopulationStructure::doMigCheck() :") + stm.str());
	}
	
	return result;
}


// ------------------ non-public methods --------------------


std::vector < std::vector < double > >
				PopulationStructure::makeMigMatrixFromSubPops() const
{
	/* For each subpopulation (in the current order)*/
	
	int nsps = nSubPops();
	
	// set up something of the right size
	std::vector < double > inner(nsps, 0.0);
	std::vector < std::vector < double > > mm(nsps, inner);
	
	/* Each inner vector 'row' i should give the inwards migration into subpop at index i
	 * For each of these we have to find the labels for the given elements
	 * and then get from subpop at i these inward migration parameters*/
	
	// save time by storing names in a vector rather than repeating lookups
	std::vector < std::string > subPopLabels;
	for (const_iterator it = begin(); it < end(); ++it) {
		subPopLabels.push_back(it->getLabel());
	}
	
	for (int i = 0; i < nsps; ++i) {
		
		for (int j = 0; j < nsps; ++j) {
			// only set for nondiagonals and if mig param is > 0.0
			if ( (i != j) ) {
				mm[i][j] = container[i].getInwardMig( subPopLabels[j] );
			}
		}
	}
	
	return mm;
}

std::set < size_t > & PopulationStructure::traceMig(
			std::set < size_t > & been,
			const std::vector < std::vector < double > >& migration_matrix,
			size_t at) const
{
	// add at to been
	std::pair < std::set < size_t >:: iterator, bool >
		ins = been.insert(at);
	
	//if (ins.second) { // if it's a new place
		
	// collect where we can go to from at
	std::set < size_t > canGetTo;
	for (size_t j = 0; j < nSubPops(); ++j) {
		if (migration_matrix[at][j] > 0.0) {
			canGetTo.insert(j);
		}
	}
	//}
	
	// find any new places we have not been
	for (std::set < size_t >:: iterator it = canGetTo.begin();
				it != canGetTo.end();
				++it) {
		std::pair < std::set < size_t >:: iterator, bool >
		ins = been.insert(*it);
	
		if (ins.second) { // if we have not been there before
			been = traceMig( been, migration_matrix, *(ins.first) );
			
		}
	}
	
	return been;
}

int PopulationStructure::precEst() const
{
	int prec = 1;
	int maxprec = 10;
		
	// rough and ready way to find an appropriate precision for printing values for ms
	for (size_t i = 0; i < nSubPops(); ++i) {
		{
			double adjval = getGrowthRateAt(i)*pow(10.0,prec);
			while ((adjval - floor(adjval) > 0.0) && (prec < maxprec)) {
				prec++;
				adjval *= 10;
			}
		}
		{
			double adjval = getSizeRatioAt(i)*pow(10.0,prec);
			while ((adjval - floor(adjval) > 0.0) && (prec < maxprec)) {
				prec++;
				adjval *= 10;
			}
		}
	}
	return prec;
}

int PopulationStructure::precEst(double val) const
{
	int prec = 1;
	int maxprec = 10;
		
	// rough and ready way to find an appropriate precision for printing values for ms
	double adjval = val*pow(10.0,prec);
	while ( (std::abs(adjval - floor(adjval)) > 0.0) && (prec < maxprec) ) {
		prec++;
		adjval *= 10;
	}
	return prec;
}


//get a default label for a new subpop
// string(order + 1);
std::string PopulationStructure::getNewDefaultLabel(int order)
{
	std::ostringstream stm;
	stm << (order + 1);
	return stm.str();
}


// erase from the actual container and from everyone else's migration rates
// and also relabel any pops with default (order-related) labels
void PopulationStructure::completeErase(size_t orderErase)
{	
	std::vector< SubPopulationStructure > tmp = container;
	
	std::string eraseLabel = tmp.at(orderErase).getLabel();
		
	std::vector< SubPopulationStructure >::iterator eraseIt = tmp.begin();
	advance(eraseIt, orderErase);
	
	tmp.erase(eraseIt);
	
	// and tell everyone else to erase it from their mig matrices
	for (std::vector< SubPopulationStructure >::iterator it = tmp.begin();
			it < tmp.end();
			++it) {
		it->eraseInwardMig(eraseLabel);
		
	}
	
	// and reset any default labels 
	for (size_t i = orderErase; i < tmp.size(); i++) {
		
		if (tmp[i].isDefaultLabel) {
			
			std::string oldDefaultLabel = tmp[i].getLabel();
			
			std::string newDefaultLabel = getNewDefaultLabel(i);
			tmp[i].setLabel(newDefaultLabel); 
			
			for (size_t j = 0; j < tmp.size(); j++) {
				if (i != j) {
					tmp[j].relabelInwardMig(oldDefaultLabel, newDefaultLabel);
				}
			}
		}
	}
	
	tmp.swap(container);
}


void PopulationStructure::completeRelabel(size_t order,
									const std::string& newLabel)
{
	
	std::string oldLabel = container.at(order).getLabel();
	
	if (oldLabel != newLabel) {
		
		container.at(order).setLabel(newLabel);
		
		// if this is not what the default would be, reset isDefaultLabel
		if (newLabel != getNewDefaultLabel(order)) {
			container.at(order).setIsDefaultLabel(false);
		}
	
		for (size_t i = 0; i < nSubPops(); ++i) {
			
			if (i != order) {
				container[i].relabelInwardMig(oldLabel, newLabel);
				
			}
		}
	}
}

// throw "label" if label not found
size_t PopulationStructure::getOrder(const std::string& label) const
{
	size_t nsps = nSubPops();
	size_t order = nsps;
	
	for (size_t i = 0; i < container.size(); ++i) {
		if (container[i].getLabel() == label) {
			order = i;
			break;
		}
	}
	if (order>= nsps) {
		throw std::invalid_argument(label);
	}

	return order;
}

bool PopulationStructure::checkUniqueLabel(size_t order,
										const std::string& label) const
{
	bool result = true;
	for (size_t i = 0; i < nSubPops() && result; ++i) {
		if ( (i != order) 
			&& container[i].getLabel() == label ) {
				result = false;
		}
	}
	return result;
}

std::string& PopulationStructure::addToErrorMsg(std::string& msg, const std::exception& e)
{
	std::string addMsg(e.what());
	if ( !addMsg.empty() ) {
		msg += ( std::string(" : ") + addMsg );
	}
	return msg;
}


std::ostringstream& PopulationStructure::outputStart(
								std::ostringstream& stm,
								size_t reps)	const
{
	stm << "ms " << totalNsam() << " " << reps;
	
	return stm;
	
} 



std::ostringstream& PopulationStructure::outputSubpops(
								std::ostringstream& stm)	const
{
	if (nSubPops() == 1) {
		double growth = getGrowthRateAt(0);
		stm.precision(precEst(growth));
		stm << " -G " << growth;
	}
	
	else { // nSubPops() > 1
		stm.precision(precEst());
	
		stm << " -I " << nSubPops(); // -I and number of subpops
		
		for (size_t i = 0; i < nSubPops(); ++i) // sub population sizes
		{
			stm << " " << getSubNsamAt(i); 
		}
		
		for (size_t i = 0; i < nSubPops(); ++i) // sub population sizes
		{
			stm << " -n " << i + 1 << " " << getSizeRatioAt(i); 
		}
		
		for (size_t i = 0; i < nSubPops(); ++i) // sub population growth rates
		{
			stm << " -g " << i + 1 << " " << getGrowthRateAt(i); 
		}
		
		// add on the migration matrix
		stm << " -ma " << migrationMatrixStringForMS();
		
	}
	return stm;
}

//  non-member functions
// Full specializations of the templates in std namespace can be added in std namespace.
template <>
void std::swap(hudson_ms::PopulationStructure & p1, 
				hudson_ms::PopulationStructure & p2) // throw ()
{
	p1.swap(p2);
}
