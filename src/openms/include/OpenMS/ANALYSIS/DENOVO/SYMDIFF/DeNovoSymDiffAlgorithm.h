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

#ifndef OPENMS_ANALYSIS_DENOVO_SYMDIFF_DENOVOSYMDIFFALGORITHM_H
#define OPENMS_ANALYSIS_DENOVO_SYMDIFF_DENOVOSYMDIFFALGORITHM_H

// OpenMS includes
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include "solve_wrapper.h"
#include "config.h"
#include "peptide_solver_symmetric_difference_general.h"
#include "symdiffscore_multyions.h"


namespace OpenMS
{
  /**
  @brief A class handling the parameters to run the symmetric difference de novo sequencing algorithm

  This class defines all parameters for the algorithm and it's default values.

  For exection the parameters are filled in a config struct and passed to the solver.

  The only parameter check that is conducted, is that the number of last_characters matches the number of last_character_tags. A basic parameter check is already done by the framework.
  */
  class OPENMS_DLLAPI DeNovoSymDiffAlgorithm :
    public DefaultParamHandler
  {


  public:

    /// Default constructor that fills the default params
    DeNovoSymDiffAlgorithm() :
      DefaultParamHandler("DeNovoSymDiffAlgorithm")
    {
      defaults_.setValue("rho", 0.9, "Relative score a solution must have to be listed in the results");

      StringList lastChars;
      lastChars.push_back("K");
      lastChars.push_back("R");
      defaults_.setValue("last_characters", lastChars, "The last character of a peptide can only be one of these characters");
      DoubleList lastCharTags;
      lastCharTags.push_back(8.0142);
      lastCharTags.push_back(10.008);
      defaults_.setValue("last_character_tag", lastCharTags, "For every last character a mass offset must be defined, use 0 if no mass offset is necessary");

      defaults_.setValue("accuracy", 0.02, "Accuracy used for equality of peaks (in Daltons)");
      defaults_.setValue("merge_accuracy", 0.04, "Accuracy used for merging peaks (in Daltons)");
      defaults_.setValue("noise_cutoff", 100.0, "Peaks with lower intensity are ignored");

      //SymDiff Parameters
      defaults_.setValue("count_peaks", "false", "Count number of peaks instead of summing up their intensities (if true, missing_peak_punishment should be -2 and half_missing_peak_punishment should be -1 )");
      defaults_.setValidStrings("count_peaks", ListUtils::create<String>("true,false"));
      defaults_.setValue("half_missing_peak_punishment", -2500.0, "Punishment if only one peak of a b-/y-ion peak pair is missing");
      defaults_.setValue("missing_peak_punishment", -5000.0, "Negative score for a missing b-/y-ion peak pair");
      defaults_.setValue("missing_peak_threshold", 1.0, "Maximum intensity a peak can have to be counted as missing");
      defaults_.setValue("missing_ion_peak_punishment", -1000.0, "Punishment if a ion other than b-/y-ion is missing");

      //Multiple IonTypes Parameters
      DoubleList b_ion_offsets;
      b_ion_offsets.push_back(-27.9949);//a-ion
      b_ion_offsets.push_back(-18.0106);//b-ion - H20
      defaults_.setValue("b_ion_offsets", b_ion_offsets, "Offsets for additional ion types to the b-ion");
      DoubleList y_ion_offsets;
      y_ion_offsets.push_back(-18.0106);//y-ion - H20
      y_ion_offsets.push_back(1);//y-ion with 13C isotope
      defaults_.setValue("y_ion_offsets", y_ion_offsets, "Offsets for additional ion types to the y-ion");

      defaults_.setValue("annotation_sequence", "", "A amino acid sequence from annotation data. This can be used for testing the algorithm with a test set where the correct sequence is already known. The tool then outputs performance values as well");

      // write defaults into Param object param_
      defaultsToParam_();
    }

    /**
  @brief This method is used to execute a search for amino acid sequences.

  The @em input is the spectrum on which the search as to be performed. All other parameters are passed by the framework over the @em param_ field.

  @exception Exception::InvalidParameter is thrown if the number of elements in the parameter last_characters does not match the number of elements in last_character_tag.

    */
    PeptideIdentification searchSequences(const MSSpectrum<Peak1D>& input, double mass)
    {

    std::vector<peak> spectrum;

      for(Size i=0; i<input.size(); i++)
      {
        peak p = {input.getRT(), input[i].getMZ(), input[i].getIntensity()};
        spectrum.push_back(p);
      }

      config conf;
      conf.EPS = param_.getValue("accuracy");
      conf.MERGING_EPS = param_.getValue("merge_accuracy");
      conf.NOISE_H_CUTOFF = param_.getValue("noise_cutoff");
      conf.COUNT_PEAKS = param_.getValue("count_peaks").toBool();
      conf.MISSING_PEAK_PUNISHMENT = param_.getValue("missing_peak_punishment");
      conf.HALF_MISSING_PEAK_PUNISHMENT = param_.getValue("half_missing_peak_punishment");
      conf.MISSING_PEAK_THRESHOLD_SCORE = param_.getValue("missing_peak_threshold");
      conf.MISSING_ION_PEAK_PUNISHMENT = param_.getValue("missing_ion_peak_punishment");
      conf.START_SCORE = -10*conf.MISSING_PEAK_PUNISHMENT;
      conf.b_ion_offsets = param_.getValue("b_ion_offsets");
      conf.y_ion_offsets = param_.getValue("y_ion_offsets");

      StringList lastChars = param_.getValue("last_characters");
      DoubleList lastCharTags = param_.getValue("last_character_tag");
      if(lastChars.size()!=lastCharTags.size())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "For every last character a corresponding tag is required! Use 0 if you don't need an offset.");
      }

      PeptideIdentification result = solve_one_prec_mass(spectrum, mass, param_.getValue("rho"), input.getRT(), conf, lastChars, lastCharTags);

      //If an annotation sequence is given we search this sequence and output maxscore and score and position of the sequence
      //Format: annotated sequence, retention time, position of annotated seq, score of annotated seq, best score,  best scoring sequence, recall of best scoring sequence
      std::string annotationSequence = param_.getValue("annotation_sequence");
      if(annotationSequence.size()>0)
      {
        replace(annotationSequence.begin(),annotationSequence.end(), 'I', 'L');//Replace all I with L because they have the same mass and our algorithm only generates L's
        std::cout << "RESULT," << annotationSequence << "," << input.getRT() << ",";
        AASequence annotation = AASequence::fromString(annotationSequence);
        if(result.getHits().size()==0)
        {
          std::cout << "Inf,Inf,Inf,Inf,0,***";
        }
        else
        {
          const double bestScore = result.getHits()[0].getScore();

          // find the position of the annotated sequence
          int solutionPos=-1;
          for(int i=0; i<result.getHits().size(); i++)
          {
            if(result.getHits()[i].getSequence()==annotation)
            {
              solutionPos = i;
              break;
            }
          }

          // find the position of the first sequence that has the same score as the annotated sequence
          while( solutionPos>0 && result.getHits()[solutionPos-1].getScore()==result.getHits()[solutionPos].getScore() )
          {
            solutionPos--;
          }

          // find the best recall of a best-scoring sequence
          int i = 0;
          double r;
          double bestRecall = 0;
          std::string bestSequence = result.getHits()[0].getSequence().toString();
          while ( result.getHits()[i].getScore() == bestScore )
          {
            r = get_recall( annotationSequence, result.getHits()[i].getSequence().toString(), conf.EPS );
            if (r > bestRecall) 
            {
              bestRecall = r;
              bestSequence = result.getHits()[i].getSequence().toString();
            }
            i++;
          }

          if (solutionPos >= 0)
          {
            std::cout << solutionPos << ",";
            std::cout << result.getHits()[solutionPos].getScore() << ",";
          }
          else
          {
            std::cout << "Inf,Inf,";
          }
          std::cout << bestScore << ",";
          std::cout << bestSequence << ",";
          std::cout << bestRecall << ",***" ;
        }
        std::cout << std::endl;
      }

      return result;
    }

  };

}

#endif