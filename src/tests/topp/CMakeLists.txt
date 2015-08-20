# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2015.
#
# This software is released under a three-clause BSD license:
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of any author or any participating institution
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
# For a full list of authors, refer to the file AUTHORS.
# --------------------------------------------------------------------------
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
# INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche, Chris Bielow $
# $Authors: Andreas Bertsch, Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------

project("OpenMS_TOPP_tests")
cmake_minimum_required(VERSION 2.8.3 FATAL_ERROR)

#
# README:
#
# This file contains the declarations of all TOPP/UTILS tests that are executed by CTest
#
# In order to add a test do the following:
# - add specific tests for your tool (this usually involves FuzzyDiff) - see below for examples
# - for UTILS we do not require any tests, but you are encouraged to provide them!
#
# Test names (the first argument of 'add_test') must be distinct.
# If not, the first definition of the test is overwritten!
#
# Therefore we use the following naming scheme for tests:
# TOPP_<tool>_<num>[_out<num2>]
#
# 'TOPP_' This prefix is used to distinguish TOPP and OpenMS tests in the nightly builds.
# <tool>  is the TOPP tool name used to identify failed TOPP tool tests during the nightly builds.
# <num>   Index that is increased for each call of a TOPP tool during a test.
#         For each new TOPP tool test this is reset to 1.
# '_out'  When the result file of a test is checked with FuzzyDiff, this suffix is appended.
# <num2>  If several files need to be checked in the same test, the checks are numbered from 1.
#
# When the result file of a test is checked with FuzzyDiff, we need to set the dependency to the test
# that generates the output file by
#
# set_tests_properties("TOPP_<tool>_<num>" PROPERTIES DEPENDS "TOPP_<tool>_<num>_out<num2>")
#
# Warning: When passing parameters avoid usage of "<" and ">" signs, as this will confuse CDASH when
#          interpreting the test results file (which happens to be XML)
#
# Hint: Most TOPP tools add DataProcessing information to their output.
#       The -test flag ensures that this output is independent of local path settings etc., so use it when *testing* TOPP tools!
# Hint: Provide individual names for each test output file e.g., the first test for the FileConverter should have FileConverter_1.tmp
#       as output file.

#------------------------------------------------------------------------------
# General settings
#------------------------------------------------------------------------------
set(TOPP_BIN_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
set(DATA_DIR_TOPP ${PROJECT_SOURCE_DIR})
set(DATA_DIR_TOPP_BIN ${PROJECT_BINARY_DIR})
set(DIFF ${TOPP_BIN_PATH}/FuzzyDiff -ini ${DATA_DIR_TOPP}/FuzzyDiff.ini)

## configure resource files (for TOPPAS/ExecutePipeline) - there is no other way to specify the input
configure_file(${DATA_DIR_TOPP}/ExecutePipeline_1.trf.in ${DATA_DIR_TOPP_BIN}/ExecutePipeline_1.trf)

#------------------------------------------------------------------------------
# test if -write_wdsl works
#------------------------------------------------------------------------------
add_test("TOPP_WSDL_1" ${TOPP_BIN_PATH}/AdditiveSeries -write_wsdl WSDL_out.tmp)
add_test("TOPP_WSDL_1_out1" ${DIFF} -in1 WSDL_out.tmp -in2 ${DATA_DIR_TOPP}/WSDL_out.wsdl)
set_tests_properties("TOPP_WSDL_1_out1" PROPERTIES DEPENDS "TOPP_WSDL_1")

#------------------------------------------------------------------------------
# WRITE_INI test for each TOPP/UTILS tool, checking -write_ini functionality
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# PeakPicker tests
add_test("TOPP_PeakPickerWavelet_1" ${TOPP_BIN_PATH}/PeakPickerWavelet  -test -ini ${DATA_DIR_TOPP}/PeakPickerWavelet_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerWavelet_input.mzML -out PeakPickerWavelet_1.tmp)
add_test("TOPP_PeakPickerWavelet_1_out1" ${DIFF} -in1 PeakPickerWavelet_1.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerWavelet_output.mzML )
set_tests_properties("TOPP_PeakPickerWavelet_1_out1" PROPERTIES DEPENDS "TOPP_PeakPickerWavelet_1")
add_test("TOPP_PeakPickerWavelet_2" ${TOPP_BIN_PATH}/PeakPickerWavelet  -test -ini ${DATA_DIR_TOPP}/PeakPickerWavelet_deconv_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerWavelet_deconv_input.mzML -out PeakPickerWavelet_2.tmp)
add_test("TOPP_PeakPickerWavelet_2_out1" ${DIFF} -in1 PeakPickerWavelet_2.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerWavelet_deconv_output.mzML)
set_tests_properties("TOPP_PeakPickerWavelet_2_out1" PROPERTIES DEPENDS "TOPP_PeakPickerWavelet_2")
add_test("TOPP_PeakPickerWavelet_3" ${TOPP_BIN_PATH}/PeakPickerWavelet  -test -ini ${DATA_DIR_TOPP}/PeakPickerWavelet_2Dopt_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerWavelet_2Dopt_input.mzML -out PeakPickerWavelet_3.tmp)
#  
# See https://github.com/OpenMS/OpenMS/issues/1202
# On Ubuntu 14.04 in Release mode only (not Debug mode), the output of this
# file is slightly different due to a small numerical error. Since FuzzyDiff is
# not able to compare the numbers in base64, we have to "whitelist" the
# corresponding base64 string.
#
#  Validation: 
#   mynr = base64.standard_b64decode("3bAxk01KjUAAAACgmVGNQIvf6JpEWo1AhaXzLTpijUA=")
#   [struct.unpack('<d', mynr[i:i+8]) for i in range(0, len(mynr), 8) ]
#   [(937.2878784067249,), (938.2000122070312,), (939.2834985917301,), (940.2784079584968,)]
#
#   mynr = base64.standard_b64decode("3bAxk01KjUAAAACgmVGNQIvf6JpEWo1AhqXzLTpijUA=") 
#   [struct.unpack('<d', mynr[i:i+8]) for i in range(0, len(mynr), 8) ]
#   [(937.2878784067249,), (938.2000122070312,), (939.2834985917301,), (940.2784079584969,)]
#
#   mynr = base64.standard_b64decode("bSj7DzRKjUDnZAujlFGNQINk5GQ2Wo1AeTS9qztijUA=") 
#   [struct.unpack('<d', mynr[i:i+8]) for i in range(0, len(mynr), 8) ]
#   [(937.2754211065018,), (938.1975766077202,), (939.2765596240209,), (940.2791361600531,)]
#
#   -- case 2
# 
#   mynr = base64.standard_b64decode("3bAxk01KjUAAAACgmVGNQIvf6JpEWo1AhqXzLTpijUAAAACgmWmNQA==") 
#   [struct.unpack('<d', mynr[i:i+8]) for i in range(0, len(mynr), 8) ]
#   [(937.2878784067249,), (938.2000122070312,), (939.2834985917301,), (940.2784079584969,), (941.2000122070312,)]
#
#   mynr = base64.standard_b64decode("3bAxk01KjUAAAACgmVGNQIvf6JpEWo1AhaXzLTpijUAAAACgmWmNQA==") 
#   [struct.unpack('<d', mynr[i:i+8]) for i in range(0, len(mynr), 8) ]
#   [(937.2878784067249,), (938.2000122070312,), (939.2834985917301,), (940.2784079584968,), (941.2000122070312,)]
#
add_test("TOPP_PeakPickerWavelet_3_out1" ${DIFF} -in1 PeakPickerWavelet_3.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerWavelet_2Dopt_output.mzML -matched_whitelist 3bAxk01KjUAAAACgmVGNQIvf6JpEWo1AhaXzLTpijUA:3bAxk01KjUAAAACgmVGNQIvf6JpEWo1AhqXzLTpijUAAAACgmWmNQA 3bAxk01KjUAAAACgmVGNQIvf6JpEWo1AhaXzLTpijUA:bSj7DzRKjUDnZAujlFGNQINk5GQ2Wo1AeTS9qztijUA 3bAxk01KjUAAAACgmVGNQIvf6JpEWo1AhqXzLTpijUAAAACgmWmNQA:3bAxk01KjUAAAACgmVGNQIvf6JpEWo1AhaXzLTpijUAAAACgmWmNQA)
set_tests_properties("TOPP_PeakPickerWavelet_3_out1" PROPERTIES DEPENDS "TOPP_PeakPickerWavelet_3")
add_test("TOPP_PeakPickerWavelet_4" ${TOPP_BIN_PATH}/PeakPickerWavelet  -test -ini ${DATA_DIR_TOPP}/PeakPickerWavelet_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerWavelet_input.mzML -out PeakPickerWavelet_4.tmp -threads 2)
add_test("TOPP_PeakPickerWavelet_4_out1" ${DIFF} -in1 PeakPickerWavelet_4.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerWavelet_output.mzML)
set_tests_properties("TOPP_PeakPickerWavelet_4_out1" PROPERTIES DEPENDS "TOPP_PeakPickerWavelet_4")
add_test("TOPP_PeakPickerWavelet_5" ${TOPP_BIN_PATH}/PeakPickerWavelet  -test -ini ${DATA_DIR_TOPP}/PeakPickerWavelet_parameters_noMetaData.ini -in ${DATA_DIR_TOPP}/PeakPickerWavelet_input.mzML -out PeakPickerWavelet_5.tmp -threads 2)
add_test("TOPP_PeakPickerWavelet_5_out1" ${DIFF} -in1 PeakPickerWavelet_5.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerWavelet_output_noMetaData.mzML)
set_tests_properties("TOPP_PeakPickerWavelet_5_out1" PROPERTIES DEPENDS "TOPP_PeakPickerWavelet_5")

# "high_res" algorithm with "ms_levels" option set to only 1:
add_test("TOPP_PeakPickerHiRes_1" ${TOPP_BIN_PATH}/PeakPickerHiRes -test -ini ${DATA_DIR_TOPP}/PeakPickerHiRes_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerHiRes_input.mzML -out PeakPickerHiRes_1.tmp)
add_test("TOPP_PeakPickerHiRes_1_out1" ${DIFF} -in1 PeakPickerHiRes_1.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerHiRes_output.mzML)
set_tests_properties("TOPP_PeakPickerHiRes_1_out1" PROPERTIES DEPENDS "TOPP_PeakPickerHiRes_1")

add_test("TOPP_PeakPickerHiRes_2" ${TOPP_BIN_PATH}/PeakPickerHiRes -test -ini ${DATA_DIR_TOPP}/PeakPickerHiRes_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerHiRes_2_input.mzML -out PeakPickerHiRes_2.tmp)
add_test("TOPP_PeakPickerHiRes_2_out1" ${DIFF} -in1 PeakPickerHiRes_2.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerHiRes_2_output.mzML)
set_tests_properties("TOPP_PeakPickerHiRes_2_out1" PROPERTIES DEPENDS "TOPP_PeakPickerHiRes_2")

# "high_res" peakpicker with the lowmemory option
add_test("TOPP_PeakPickerHiRes_3" ${TOPP_BIN_PATH}/PeakPickerHiRes -test -ini ${DATA_DIR_TOPP}/PeakPickerHiRes_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerHiRes_input.mzML -out PeakPickerHiRes_3.tmp -processOption lowmemory)
add_test("TOPP_PeakPickerHiRes_3_out1" ${DIFF} -in1 PeakPickerHiRes_3.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerHiRes_output.mzML)
set_tests_properties("TOPP_PeakPickerHiRes_3_out1" PROPERTIES DEPENDS "TOPP_PeakPickerHiRes_3")

add_test("TOPP_PeakPickerHiRes_4" ${TOPP_BIN_PATH}/PeakPickerHiRes -test -ini ${DATA_DIR_TOPP}/PeakPickerHiRes_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerHiRes_2_input.mzML -out PeakPickerHiRes_4.tmp -processOption lowmemory)
add_test("TOPP_PeakPickerHiRes_4_out1" ${DIFF} -in1 PeakPickerHiRes_4.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerHiRes_2_output.mzML)
set_tests_properties("TOPP_PeakPickerHiRes_4_out1" PROPERTIES DEPENDS "TOPP_PeakPickerHiRes_4")

# peak picking of a smoothed chromatogram (like in PeakPickerMRM; requires disabling of spacing constraints):
add_test("TOPP_PeakPickerHiRes_5" ${TOPP_BIN_PATH}/PeakPickerHiRes -test -in ${DATA_DIR_TOPP}/PeakPickerHiRes_5_input.mzML -out PeakPickerHiRes_5.tmp)
add_test("TOPP_PeakPickerHiRes_5_out1" ${DIFF} -in1 PeakPickerHiRes_5.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerHiRes_5_output.mzML)
set_tests_properties("TOPP_PeakPickerHiRes_5_out1" PROPERTIES DEPENDS "TOPP_PeakPickerHiRes_5")

