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

#ifndef OPENMS_ANALYSIS_DENOVO_SYMDIFF_SOLVEWRAPPER_H
#define OPENMS_ANALYSIS_DENOVO_SYMDIFF_SOLVEWRAPPER_H

#include <string>
#include <iostream>
#include <map>
#include <algorithm>
#include "config.h"

#include <OpenMS/METADATA/PeptideIdentification.h>

using namespace OpenMS;

#include "symdiffscore_multyions.h"
#include "amino_acids.h"
#include "peptide_solver_symmetric_difference_general.h"

/// H_MASS. The mass of H molecule in Daltons.
const double H_MASS = 1.00794;
/// H20_MASS. The mass of H20 molecule in Daltons.
const double H20_MASS = 18.01048;

/// PeptideOrder. Is used to sort the peptides by their score.
bool peptideOrder(peptide a, peptide b) 
{ 
  return a.first > b.first || (a.first==b.first && a.second > b.second); 
}

/// Mergepeak. This represents a peak pair and is needed in the merge process.
typedef struct{
  double time;
  double mz;
  double h;
  double low_mass_h;
  double high_mass_h;
  bool deleted;
} mergepeak;


/**
  @brief This method merges similar peaks and then calls the search routine

  The merge process simultaniously merges also the complementary peaks. For this for every measured peak the complementary peak is generated. To do this we have to know the offset at the last character @em offset. Afterwards peaks from 0 to W/2 are merged with a greedy approach. Because every merged peak represents two complementary peak, they are restored after the merge process. At the end the search routine is called.

*/
void solve_prec_mass_with_offset(std::vector<peak>& spectrum, double mass, double offset, std::vector<peptide>& result, char last_letter, double rho, config conf)
{

  if(spectrum.size() == 0) 
  {
    return;
  }

  double W = mass - offset - H20_MASS - H_MASS; //Compute the mass of the amino acid sequence without the C-terminus and without the tag of the last amino acid

  std::vector<peak> v_spectrum;
  for (int i=0; i<spectrum.size(); i++)
  {
    if (spectrum[i].mz < mass)
    {
      v_spectrum.push_back(spectrum[i]);
    }
  }

  std::vector<mergepeak> potential_b_peaks;
  std::set<std::pair<double,int> > peakIntensityOrdering;
  int front=0;
  int back=v_spectrum.size()-1;
  //Generate all potential b-ion-masses (uncharged) in range 0 - (mass-H_MASS)/2
  while (front<=back) 
  {
    if(mass-v_spectrum[back].mz>v_spectrum[front].mz-H_MASS)
    {
      peak p = v_spectrum[front];
      mergepeak mpeak;
      mpeak.time = p.time;
      mpeak.mz = p.mz-H_MASS;
      mpeak.h = p.h;
      mpeak.low_mass_h = p.h;
      mpeak.high_mass_h = 0;
      mpeak.deleted = false;
      peakIntensityOrdering.insert(std::make_pair(-p.h,potential_b_peaks.size()));
      potential_b_peaks.push_back(mpeak);
      front++;
    }
    else
    {
      peak p = v_spectrum[back];
      mergepeak mpeak;
      mpeak.time = p.time;
      mpeak.mz = mass-p.mz;
      mpeak.h = p.h;
      mpeak.low_mass_h = 0;
      mpeak.high_mass_h = p.h;
      mpeak.deleted = false;
      peakIntensityOrdering.insert(std::make_pair(-p.h,potential_b_peaks.size()));
      potential_b_peaks.push_back(mpeak);
      back--;
    }
  }

  //Greedy algorithm for peak merging: go through peaks in decreasing intensity order. For every peak merge with closest neighbor and update position as weighted average, repeat this until no peak is closer than eps.
  for (std::set<std::pair<double,int> >::iterator iter=peakIntensityOrdering.begin(); iter!=peakIntensityOrdering.end(); iter++)
  {
    std::pair<double,int> peakPos = *iter;
    int pos = peakPos.second;
    if(potential_b_peaks[pos].deleted)
    {
      continue;
    }

    int left = pos-1;
    int right = pos+1;
    //merge with undeleted neighbors
    while ((left>=0 && fabs(potential_b_peaks[left].mz-potential_b_peaks[pos].mz)<conf.MERGING_EPS) || (right<potential_b_peaks.size() && fabs(potential_b_peaks[right].mz-potential_b_peaks[pos].mz)<conf.MERGING_EPS) )
    {

      double leftDiff=conf.MERGING_EPS;
      if(left>=0)
      {
        leftDiff = fabs(potential_b_peaks[left].mz-potential_b_peaks[pos].mz);
      }
      double rightDiff=conf.MERGING_EPS;
      if(right<potential_b_peaks.size())
      {
        rightDiff = fabs(potential_b_peaks[right].mz-potential_b_peaks[pos].mz);
      }

      if(leftDiff<rightDiff)
      {
        if(!potential_b_peaks[left].deleted)
        {
          potential_b_peaks[pos].mz = (potential_b_peaks[pos].mz*potential_b_peaks[pos].h + potential_b_peaks[left].mz*potential_b_peaks[left].h)/(potential_b_peaks[pos].h+potential_b_peaks[left].h);
          potential_b_peaks[pos].h += potential_b_peaks[left].h;
          potential_b_peaks[pos].low_mass_h += potential_b_peaks[left].low_mass_h;
          potential_b_peaks[pos].high_mass_h += potential_b_peaks[left].high_mass_h;
          potential_b_peaks[left].deleted = true;
        }
        left--;
      }
      else
      {
        if(!potential_b_peaks[right].deleted)
        {
          potential_b_peaks[pos].mz = (potential_b_peaks[pos].mz*potential_b_peaks[pos].h + potential_b_peaks[right].mz*potential_b_peaks[right].h)/(potential_b_peaks[pos].h+potential_b_peaks[right].h);
          potential_b_peaks[pos].h += potential_b_peaks[right].h;
          potential_b_peaks[pos].low_mass_h += potential_b_peaks[right].low_mass_h;
          potential_b_peaks[pos].high_mass_h += potential_b_peaks[right].high_mass_h;
          potential_b_peaks[right].deleted = true;
        }
        right++;
      }
    }
  }

  std::vector<peak> peaks;
  double cutoff = 0;

  do 
  {
    peaks.clear();

    peak startpeak = {0, 0, 0};
    peaks.push_back(startpeak);
    //Create list of peaks in whole range with h > NOISE_H_CUTOFF
    for(int i=0; i<potential_b_peaks.size(); i++)
    {
      mergepeak p = potential_b_peaks[i];
      if(!p.deleted && p.h > cutoff)
      {
        peak pp;
        pp.time = p.time;
        pp.mz = p.mz;
        pp.h = p.low_mass_h;
        peaks.push_back(pp);
      }
    }
    for(int i=potential_b_peaks.size()-1; i>=0; i--)
    {
      mergepeak p = potential_b_peaks[i];
      if(!p.deleted && p.h > cutoff)
      {
        peak pp;
        pp.time = p.time;
        pp.mz = (mass-p.mz)-1*H_MASS;
        pp.h = p.high_mass_h;
        peaks.push_back(pp);
      }
    }
    peak endpeak = {0, W, 0};
    peaks.push_back(endpeak);

    cutoff += conf.NOISE_H_CUTOFF;

  } while(peaks.size()>500);

  //Call the search routine that finds amino acid sequences with the peaks we generated now
  find_all_peptides(peaks, last_letter, rho, result, conf);

}

/**
  @brief This method invokes the search procedure for every last character

  This method invokes the search procedure for every @em lastChars. Afterwards all solutions are added to a PeptideIdentification which is then returned. Because for every last character an individual search is executed, we have to check afterwards that we only return the solutions with a score>=rho*max_score, where the max_score is the max of the maximum scores of the individual searches.
*/
PeptideIdentification solve_one_prec_mass(std::vector<peak>& spectrum, double mass, double rho, double time, config conf, StringList lastChars, DoubleList lastCharTags) 
{

  init_mass_table();
  std::vector<peptide> peptides;
  for(int i=0; i<lastChars.size(); i++)
  {
    solve_prec_mass_with_offset(spectrum, mass, lastCharTags[i], peptides, lastChars[i][0], rho, conf);
  }

  PeptideIdentification identification;
  identification.setRT(time);
  identification.setHigherScoreBetter(true);
    identification.setSignificanceThreshold(rho);

  if(peptides.size()>0)
  {

    double max_score=0;
    std::sort(peptides.begin(), peptides.end(), peptideOrder);
    max_score = std::max(max_score, peptides[0].first);

    for(int j=0; j<peptides.size(); j++)
    {
      peptide p = peptides[j];
      if(p.first >= rho*max_score)
      {
        //Generate peptide hit (remove brackets of multi amino acid edges)
        std::stringstream ss;
        for(int si=0; si<p.second.size(); si++)
        {
          for(int ci=0; ci<p.second[si].size(); ci++)
          {
            char c = p.second[si][ci];
            if(c!='(' && c!=')')
            {
              ss << c;
            }
          }
        }
        PeptideHit hit = PeptideHit(p.first, 0, 1, AASequence::fromString(ss.str()));
        identification.insertHit(hit);
      }
    }

    identification.assignRanks();
  }

  return identification;
}

#endif
