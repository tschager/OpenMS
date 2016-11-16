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

#ifndef OPENMS_ANALYSIS_DENOVO_SYMDIFF_SYMDIFFSCOREMULTYIONS_H
#define OPENMS_ANALYSIS_DENOVO_SYMDIFF_SYMDIFFSCOREMULTYIONS_H

#include <iostream>
#include <cmath>
#include "config.h"
#include "data_types.h"
#include "amino_acids.h"
#include <set>

///simple_edge. A holder for an edge an its label, as well as an id.
typedef struct{
	int from;
	int to;
	string label;
	int id;
} simple_edge;

const double NEG_INF = -std::numeric_limits<double>::infinity();

///mass_index. A pair of a double mass and its corresponding peak index. index -1 is used for no peak.
typedef struct{
  double mass;
  int index;
  bool punishType;
} mass_index;

struct comp_mass_index {
  bool operator() (const mass_index& lhs, const mass_index& rhs) const 
  {
    return lhs.mass < rhs.mass;
  }
};

/**
  @brief This method searches a peak representing the @em mass.
*/
mass_index getIndexForMass(vector<peak>& peaks, double mass, config conf, bool punishType)
{
	for(int i=0; i<peaks.size(); i++)
	{
		if(fabs(peaks[i].mz-mass)<conf.EPS)
		{
			return {peaks[i].mz, i};
		}
	}
	
	return  {mass, -1, punishType};
}

///Caching for the computed peaks for every vertex
vector<vector<mass_index> > explainedPeakCacheVector;
vector<set<mass_index, comp_mass_index> > explainedPeakCache;

/**
  @brief This method computes all ion peaks for a given b-ion peak.
*/
vector<mass_index> getExplainedPeaksForIndexVector(vector<peak>& peaks, int index, config conf)
{
	vector<mass_index> result;
	if(index>0 && (index>=explainedPeakCacheVector.size() || explainedPeakCacheVector[index].size()==0))
	{
		
		mass_index b;
		b.punishType = false;
		mass_index y;
		y.punishType = false;

		if(index%2==0)
		{
			b.mass = peaks[index/2].mz;
			b.index = index/2;
			y.mass = peaks[peaks.size()-1-index/2].mz;
			y.index = peaks.size()-1-index/2;
		}
		else
		{
			b.mass = peaks[peaks.size()-1-index/2].mz;
			b.index = peaks.size()-1-index/2;
			y.mass = peaks[index/2].mz;
			y.index = index/2;
		}

		result.push_back(b);
		result.push_back(y);
		
		for(int i=0; i<conf.b_ion_offsets.size(); i++)
		{
			result.push_back(getIndexForMass(peaks, b.mass+conf.b_ion_offsets[i], conf, true));
		}
		
		for(int i=0; i<conf.y_ion_offsets.size(); i++)
		{
			result.push_back(getIndexForMass(peaks, y.mass+conf.y_ion_offsets[i], conf, true));
		}

		explainedPeakCacheVector.resize(max((int)explainedPeakCache.size(),index+1));
		explainedPeakCacheVector[index] = result;

	}
	else if(index>0)
	{
		result = explainedPeakCacheVector[index];	
	}

	return result;
}

/**
  @brief This method does the same as getExplainedPeaksForIndexVector but returns the result in a set.
*/
set<mass_index, comp_mass_index> getExplainedPeaksForIndex(vector<peak>& peaks, int index, config conf)
{
	set<mass_index, comp_mass_index> result;
	if(index>0 && (index>=explainedPeakCache.size() || explainedPeakCache[index].size()==0))
	{
		
		vector<mass_index> v = getExplainedPeaksForIndexVector(peaks, index, conf);
		for(int i=0; i<v.size(); i++ )
		{
			result.insert(v[i]);
		}

		explainedPeakCache.resize(max((int)explainedPeakCache.size(),index+1));
		explainedPeakCache[index] = result;

	}
	else if(index>0)
	{
		result = explainedPeakCache[index];	
	}

	return result;
}

/**
  @brief This method computes the score for an edge, under the condition that @em e_given was already scored.

*/
double score(vector<peak>& peaks, simple_edge e, simple_edge e_given, config conf)
{
	double score = 0;

	//vertex score
	set<mass_index, comp_mass_index> given = getExplainedPeaksForIndex(peaks, e_given.to, conf);
	vector<mass_index> new_peaks = getExplainedPeaksForIndexVector(peaks, e.to, conf);

	for(int i=0; i<new_peaks.size(); i++)
	{
		mass_index p = new_peaks[i];
		if(p.index>=0)
		{
			if(given.find(p)==given.end())
			{
				if(peaks[p.index].h < conf.MISSING_PEAK_THRESHOLD_SCORE)
				{
					score += conf.MISSING_ION_PEAK_PUNISHMENT;
				}
				else
				{
					score += peaks[p.index].h;
				}
			}
		}
		else
		{
			bool add=true;
			mass_index low;
			low.mass = p.mass-conf.EPS;
			mass_index upp;
			upp.mass = p.mass+conf.EPS;
			for(set<mass_index,comp_mass_index>::iterator it=given.lower_bound(low); it!=given.upper_bound(upp); ++it)
			{
				add=false;
			}
			if(add && p.punishType)
			{
				score += conf.MISSING_ION_PEAK_PUNISHMENT;
			}
		}
	}

	//Punish half missing b/y peak
	if(e.to/2 != e_given.to/2 && e.to/2+e_given.to/2!=peaks.size()-1)
	{
		if(peaks[e.to/2].mz != peaks[peaks.size()-1-e.to/2].mz)
		{
			if(peaks[e.to/2].h < conf.MISSING_PEAK_THRESHOLD_SCORE)
			{
				score += conf.HALF_MISSING_PEAK_PUNISHMENT;
			}
			else
			{
				//score += peaks[e.to/2].h;
			}
			if(peaks[peaks.size()-1-e.to/2].h < conf.MISSING_PEAK_THRESHOLD_SCORE)
			{
				score += conf.HALF_MISSING_PEAK_PUNISHMENT;
			}
			else
			{
				//score += peaks[peaks.size()-1-e.to/2].h;
			}
		}
		else
		{
			//score += peaks[e.to/2].h;
		}
	}

	//missing b-/y-peaks score for multi amino acid label edges
	if(e.label.size()>1)
	{
		if(e_given.label.size()>1)
		{
			if(fabs( (peaks[e.from/2].mz+get_aminoacid_mass(e.label[1])) - (peaks[e_given.from/2].mz + get_aminoacid_mass(e_given.label[1])) )<conf.EPS)
			{
				//The missing peak was already accounted for
				score += 0;
			}
			else
			{
				score += conf.MISSING_PEAK_PUNISHMENT;
			}
		}
		else
		{
			score += conf.MISSING_PEAK_PUNISHMENT;
		}
	}

	return score;
}

/**
  @brief This method clears the explained peak cache.

  @Note This method has to be called before every new run of the algorithm.
*/
void scoreCacheClear()
{
	explainedPeakCache.clear();
	explainedPeakCacheVector.clear();
}

#endif
