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

std::set<std::pair<double, std::string> > amino_masses;
std::set<peak, comp_peak> spectrum;

std::pair<double, std::string> aminoacids[19] = 
{
  std::mp(57.02147,"G"), 
  std::mp(71.03712,"A"), 
  std::mp(87.03203,"S"), 
  std::mp(97.05277,"P"), 
  std::mp(99.06842,"V"), 
  std::mp(101.04768,"T"), 
  std::mp(103.00919,"C"), 
  std::mp(113.08407,"L"), 
  std::mp(114.04293,"N"), 
  std::mp(115.02695,"D"), 
  std::mp(128.05858,"Q"), 
  std::mp(128.09497,"K"), 
  std::mp(129.04260,"E"), 
  std::mp(131.04049,"M"), 
  std::mp(137.05891,"H"), 
  std::mp(147.06842,"F"), 
  std::mp(156.10112,"R"), 
  std::mp(163.06333,"Y"), 
  std::mp(186.07932,"W")
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
        std::pair<double,std::string> p = std::mp(aminoacids[i].first + aminoacids[j].first,  "("+aminoacids[i].second + aminoacids[j].second+")");
        amino_masses.insert(p);
      }
    }
  }
}

/**
  @brief Returns a list of possible labels for a mass difference
*/
std::vector<std::string> get_label(double mass, config conf) 
{
  std::vector<std::string> rez;
  std::set<std::pair<double,std::string> >::iterator it = amino_masses.lower_bound(std::mp(mass - conf.EPS, ""));
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
    std::pair<double,std::string> p = aminoacids[i];
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
std::vector<double> get_prefix_masses( std::string seq ) 
{
  std::vector<double> prefix_masses;
  for ( int i = 0; i < seq.size(); i++ )
  {
    std::string aa = seq.substr(i, 1);
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
double get_recall( std::string seq1, std::string seq2, const double eps)
{
  std::vector<double> prefix_masses_seq1 = get_prefix_masses(seq1);
  std::vector<double> prefix_masses_seq2 = get_prefix_masses(seq2);
  std::set<double> intersection( prefix_masses_seq1.begin(), prefix_masses_seq1.end() );

  bool found = true;
  double low,high;
  while (found)
  {
    found = false;
    std::set<double>::iterator iter;
    for (iter = intersection.begin(); iter != intersection.end(); ++iter) 
    {
      double m = *iter;
      std::vector<double>::iterator up = upper_bound(prefix_masses_seq2.begin(), prefix_masses_seq2.end(), m);
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
