// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Thomas Tschager $
// $Authors: Simon Rösch and Thomas Tschager $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_DENOVO_SYMDIFF_AMINOACIDS_H
#define OPENMS_ANALYSIS_DENOVO_SYMDIFF_AMINOACIDS_H

#include <string>
#include <map>
#include <cmath>
#include <set>
#include "data_types.h"
#include "config.h"

set<pair<double, string> > amino_masses;
set<peak,comp_peak> spectrum;

pair<double,string> aminoacids[19] = 
{
	mp(57.02147,"G"),
	mp(71.03712,"A"),
	mp(87.03203,"S"),
	mp(97.05277,"P"),
	mp(99.06842,"V"),
	mp(101.04768,"T"),
	mp(103.00919,"C"),
	mp(113.08407,"L"),
	mp(114.04293,"N"),
	mp(115.02695,"D"),
	mp(128.05858,"Q"),
	mp(128.09497,"K"),
	mp(129.04260,"E"),
	mp(131.04049,"M"),
	mp(137.05891,"H"),
	mp(147.06842,"F"),
	mp(156.10112,"R"),
	mp(163.06333,"Y"),
	mp(186.07932,"W")
};

/**
  @brief Initializes a table of all possible amino acid combinations up to length two. This table can then be used to find the edge labels.
*/
void init_mass_table() 
{
	if(amino_masses.size()==0)
	{
		// Insert amino acids and pairs of aminoacids in order to permit gaps in the list of prefixes.
	  	for (int i = 0; i < 19;++i)
		{
			amino_masses.insert(mp(aminoacids[i].first, aminoacids[i].second));
	  	}
	  	for (int i = 0; i < 19;++i) 
	  	{
		  	for (int j = 0; j < 19;++j) 
			{
				pair<double,string> p = mp(aminoacids[i].first + aminoacids[j].first,  "("+aminoacids[i].second + aminoacids[j].second+")");
				amino_masses.insert(p);
			}
	  	}
	}
}

/**
  @brief Returns a list of possible labels for a mass difference
*/
vector<string> get_label(double mass, config conf) 
{
	vector<string> rez;
	set<pair<double,string> >::iterator it = amino_masses.lower_bound(mp(mass - conf.EPS, ""));
	for (; it!= amino_masses.end() && (*it).first - mass <= conf.EPS ; ++it) 
	{
		if (fabs((*it).first - mass) <= conf.EPS) 
		{
			rez.push_back((*it).second);
		}
	}
	return rez;
}

/**
  @brief Returns the mass of the amino acid @em c.
*/
double get_aminoacid_mass(char c)
{
	for(int i=0; i<sizeof(aminoacids)/sizeof(aminoacids[0]); i++)
	{
		pair<double,string> p = aminoacids[i];
		if(p.second[0]==c)
		{
			return p.first;
		}
	}
	return 0;
}

/**
  @brief Returns the set of prefix masses of an amino acid sequence seq
*/
vector<double> get_prefix_masses( string seq ) 
{
	vector<double> prefix_masses;
	for ( int i = 0; i < seq.size(); i++ )
	{
		string aa = seq.substr(i,1);
		// substitute all L's by I's
		if ( aa == "L" )
		{
			aa == "I";
		}
		if ( prefix_masses.size() == 0 )
		{
			prefix_masses.push_back( get_aminoacid_mass(aa.at(0)) );
		}
		else
		{
			prefix_masses.push_back( prefix_masses.back() + get_aminoacid_mass(aa.at(0)) );
		}
	}
	return prefix_masses;
}


/**
  @brief Returns the recall of an amino acid sequence seq2 wrt. the true amino acid sequence seq1 and a given accuracy epsilon>0
  recall = number of common prefix masses of seq1 and seq2 divided by number of prefix masses of seq1
*/
double get_recall( string seq1, string seq2, const double eps)
{
	vector<double> prefix_masses_seq1 = get_prefix_masses(seq1);
	vector<double> prefix_masses_seq2 = get_prefix_masses(seq2);
	set<double> intersection( prefix_masses_seq1.begin(), prefix_masses_seq1.end() );

	bool found = true;
	double low,high;
	while (found)
	{
		found = false;
		set<double>::iterator iter;
		for (iter = intersection.begin(); iter != intersection.end(); ++iter) 
		{
			double m = *iter;
			vector<double>::iterator up = upper_bound(prefix_masses_seq2.begin(), prefix_masses_seq2.end(), m);
			if ( up != prefix_masses_seq2.begin() ) 
			{
				low = *(up - 1);
			}
			else
			{
				low = prefix_masses_seq2.front();
			}
			high = *up;
			if ( fabs(high - m) > eps && fabs(low-m) > eps )
			{
				intersection.erase(m);
				found = true;
			}
		}
	}
	
	double recall = (double)intersection.size() / prefix_masses_seq1.size();
	return recall;
}

#endif
