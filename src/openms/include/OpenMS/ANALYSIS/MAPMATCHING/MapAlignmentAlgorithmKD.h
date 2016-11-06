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
// $Maintainer: Johannes Veit $
// $Authors: Johannes Veit $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMKD_H
#define OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMKD_H

#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeData.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>

namespace OpenMS
{

/**
    @brief An efficient reference-free feature map alignment algorithm for unlabeled data

    This algorithm uses a kd-tree to efficiently compute conflict-free connected components (CCC)
    in a compatibility graph on feature data. This graph is comprised of nodes corresponding
    to features and edges connecting features f and f' iff both are within each other's tolerance
    windows (wrt. RT and m/z difference). CCCs are those CCs (i) that do not contain multiple features
    from the same input map, (ii) whose total diameter does not violate the user-specified m/z and RT
    tolerance (since CCs are contiguity clusters, they can suffer from chaining and easily exceed
    these thresholds, as these are defined only on pairs of features and not on clusters) and (iii)
    whose features all have the same charge state.

    All CCCs above a user-specified minimum size are considered true sets of corresponding features
    and based on these, LOWESS transformations are computed for each input map such that the average
    deviation from the mean retention time within each CCC is minimized.

    @htmlinclude OpenMS_MapAignmentAlgorithmKD.parameters

    @ingroup MapAlignment
*/

class OPENMS_DLLAPI MapAlignmentAlgorithmKD
{
public:

  /// Constructor
  MapAlignmentAlgorithmKD(Size num_maps, Size min_cc_size);

  /// Default destructor
  virtual ~MapAlignmentAlgorithmKD();

  /// Compute data points needed for RT transformation in the current @p kd_data, add to fit_data_
  void addRTFitData(const KDTreeData& kd_data);

  /// Fit LOWESS to fit_data_, store final models in transformations_
  void fitLOWESS();

  /// Transform RTs for @p kd_data
  void transform(KDTreeData& kd_data) const;

protected:

  /// Compute connected components, store CC indices in member cc_index. Return number of CCs.
  Size computeCCs_(const KDTreeData& kd_data, std::vector<Size>& cc_index) const;

  /// Return connected components
  void getCCs_(const KDTreeData& kd_data, std::map<Size, std::vector<Size> >& result) const;

  /// Filter connected components (return conflict-free CCs of sufficiently large size and small diameter)
  void filterCCs_(const KDTreeData& kd_data, std::map<Size, std::vector<Size> >& filtered_ccs, const std::map<Size, std::vector<Size> >& ccs) const;

private:

  /// Default constructor is not supposed to be used.
  MapAlignmentAlgorithmKD();

  /// RT data for fitting the LOWESS
  std::vector<TransformationModel::DataPoints> fit_data_;

  /// LOWESS transformations
  std::vector<TransformationModelLowess*> transformations_;

  /// Minimum size for a connected component to be used for alignment
  double min_cc_size_;

};

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMKD_H
