## the directory name
set(directory source/TEST)

set(concept_executables_list
	ClassTest_test
	Exception_Base_test
	FactoryBase_test
	Factory_test
	SingletonRegistry_test
	VersionInfo_test
	FuzzyStringComparator_test
)

set(datastructures_executables_list
	String_test
	StringList_test
	IntList_test
	DoubleList_test
	Date_test
	Map_test
	DPosition_test
	DIntervalBase_test
	DRange_test
	DBoundingBox_test
	DataValue_test
	DateTime_test
	RangeManager_test
	Matrix_test
	ConvexHull2D_test
	DefaultParamHandler_test
	SuffixArray_test
	SuffixArraySeqan_test
	SuffixArrayTrypticSeqan_test
	SuffixArrayTrypticCompressed_test
	SuffixArrayPeptideFinder_test
	SparseVector_test
	DistanceMatrix_test
	Adduct_test
	ChargePair_test
	Compomer_test
	MassExplainer_test
	CVMappingTerm_test
	CVMappingRule_test
	CVReference_test
	CVMappings_test
)

set(metadata_executables_list
  MetaInfoRegistry_test
  MetaInfo_test
  MetaInfoInterface_test
  ScanWindow_test
  IonSource_test
  IonDetector_test
  IDTagger_test
  MassAnalyzer_test
  Instrument_test
  ProteinHit_test
  PeptideHit_test
  SampleTreatment_test
  Digestion_test
  Modification_test
  Tagging_test
  Sample_test
  Acquisition_test
  ProteinIdentification_test
  PeptideIdentification_test
  AcquisitionInfo_test
  Precursor_test
  Product_test
  ContactPerson_test
  SourceFile_test
  Software_test
  MetaInfoDescription_test
  DataProcessing_test
  InstrumentSettings_test
  SpectrumSettings_test
  Gradient_test
  HPLC_test
  ExperimentalSettings_test
  DocumentIdentifier_test
)

set(system_executables_list
  StopWatch_test
  File_test
  FileWatcher_test
)

set(kernel_executables_list
  PeakIndex_test
  Peak1D_test
  Peak2D_test
  DPeak_test
  RichPeak1D_test
  RichPeak2D_test
  DRichPeak_test
  ConstRefVector_test
  ComparatorUtils_test
  StandardTypes_test
  Feature_test
  FeatureMap_test
  RangeUtils_test
  MSSpectrum_test
  MSExperiment_test
  ConsensusMap_test
  ConsensusFeature_test
  AreaIterator_test
)

set(visual_executables_list
  MultiGradient_test
  AxisTickCalculator_test
)

set(format_executables_list
  CsvFile_test
  XMLValidator_test
  TextFile_test
  XMLFile_test
  LibSVMEncoder_test
  DTAFile_test
  DTA2DFile_test
  MascotInfile_test
  MascotOutfile_test
  Base64_test
  PersistentObject_test
  FASTAFile_test
  FileHandler_test
  MzXMLFile_test
  MzMLFile_test
  MzDataFile_test
  Param_test
  FeatureXMLFile_test
  MascotXMLFile_test
  PeakFileOptions_test
  PeakTypeEstimator_test
  InspectInfile_test
  InspectOutfile_test
  SequestInfile_test
  SequestOutfile_test
  PepNovoInfile_test
  PepNovoOutfile_test
  PTMXMLFile_test
  ConsensusXMLFile_test
  OMSSAXMLFile_test
  IdXMLFile_test
  BigString_test
  XTandemXMLFile_test
  MSPFile_test
  OMSSACSVFile_test
  XTandemInfile_test
  UnimodXMLFile_test
  PepXMLFile_test
  PepXMLFileMascot_test
  ControlledVocabulary_test
  TransformationXMLFile_test
  SemanticValidator_test
  MzMLValidator_test
  MS2File_test
	MascotGenericFile_test
	MascotRemoteQuery_test
	CVMappingFile_test
	#MzIdentMLFile_test
	#MzIdentMLValidator_test
	#TraMLFile_test
)

if (DB_TEST)
  list(APPEND format_executables_list DBConnection_test)
  list(APPEND format_executables_list DBAdapter_test)
endif()

if (USE_ANDIMS)
  list(APPEND format_executables_list ANDIFile_test)
endif()

set(math_executables_list
  LinearInterpolation_test
  AveragePosition_test
  BasicStatistics_test
  AsymmetricStatistics_test
  StatisticFunctions_test
  ROCCurve_test
  Histogram_test
  MathFunctions_test
  LinearRegression_test
  BilinearInterpolation_test
  GammaDistributionFitter_test
  GaussFitter_test
  NonNegativeLeastSquaresSolver_test
)

set(filtering_executables_list
  LinearResampler_test
  MorphologicalFilter_test
  SignalToNoiseEstimator_test
  SignalToNoiseEstimatorMedian_test
  SignalToNoiseEstimatorMeanIterative_test
  GaussFilter_test
  PreprocessingFunctor_test
  NLargest_test
  Scaler_test
  Normalizer_test
  ParentPeakMower_test
  SqrtMower_test
  ThresholdMower_test
  WindowMower_test
  FilterFunctor_test
  ComplementFilter_test
  GoodDiffFilter_test
  IntensityBalanceFilter_test
  IsotopeDiffFilter_test
  NeutralLossDiffFilter_test
  TICFilter_test
  PeakMarker_test
  ComplementMarker_test
  IsotopeMarker_test
  NeutralLossMarker_test
  MarkerMower_test
  SavitzkyGolayFilter_test
  InternalCalibration_test
  IDFilter_test
  TOFCalibration_test
  BernNorm_test
  DataFilters_test
)

set(comparison_executables_list
  SpectrumCheapDPCorr_test
  PeakSpectrumCompareFunctor_test
  SpectrumPrecursorComparator_test
  ZhangSimilarityScore_test
  SpectrumAlignment_test
  SpectrumAlignmentScore_test
  SpectraSTSimilarityScore_test
  SteinScottImproveScore_test
  CompareFouriertransform_test
  BinnedSpectrum_test
  BinnedSpectrumCompareFunctor_test
  BinnedSumAgreeingIntensities_test
  BinnedSpectralContrastAngle_test
  BinnedSharedPeakCount_test
  ClusterFunctor_test
  SingleLinkage_test
  CompleteLinkage_test
  AverageLinkage_test
  PeakAlignment_test
  ClusterHierarchical_test
  ClusterAnalyzer_test
  EuclideanSimilarity_test
)

set(chemistry_executables_list
  IsotopeDistribution_test
  Element_test
  ElementDB_test
  EmpiricalFormula_test
  ResidueModification_test
  Residue_test
  ResidueDB_test
  AASequence_test
  AAIndex_test
  EnzymaticDigestion_test
  TheoreticalSpectrumGenerator_test
  AdvancedTheoreticalSpectrumGenerator_test
  ModifierRep_test
  FastaIterator_test
  FastaIteratorIntern_test
  PepIterator_test
  EdwardsLippertIterator_test
  EdwardsLippertIteratorTryptic_test
  TrypticIterator_test
  ModificationsDB_test
  ModificationDefinition_test
  ModificationDefinitionsSet_test
)


set(analysis_executables_list
  ItraqChannelExtractor_test
  ItraqQuantifier_test
	ItraqConstants_test
  SVMWrapper_test
  TransformationDescription_test
  FeatureHandle_test
  MapAlignmentAlgorithm_test
  MapAlignmentAlgorithmSpectrumAlignment_test
  MapAlignmentAlgorithmPoseClustering_test
  LabeledPairFinder_test
  BaseGroupFinder_test
  DelaunayPairFinder_test
  SimplePairFinder_test
  StablePairFinder_test
  BaseSuperimposer_test
  PoseClusteringAffineSuperimposer_test
  PoseClusteringShiftSuperimposer_test
  FeatureGroupingAlgorithm_test
  FeatureGroupingAlgorithmLabeled_test
  FeatureGroupingAlgorithmUnlabeled_test
	FeatureDeconvolution_test
  ConsensusID_test
  ProteinInference_test
  PILISScoring_test
  IDDecoyProbability_test
  FalseDiscoveryRate_test
  LocalLinearMap_test
  PeakIntensityPredictor_test
  IDMapper_test
  PrecursorIonSelectionPreprocessing_test
  PrecursorIonSelection_test
  MapAlignmentEvaluationAlgorithm_test
  MapAlignmentEvaluationAlgorithmPrecision_test
  MapAlignmentEvaluationAlgorithmRecall_test
  OfflinePrecursorIonSelection_test
	DeNovoAlgorithm_test
	DeNovoIdentification_test
	DeNovoIonScoring_test
	DeNovoPostScoring_test
	MassDecomposition_test
	MassDecompositionAlgorithm_test
	CompNovoIdentificationBase_test 
	CompNovoIdentification_test
	CompNovoIonScoringCID_test
	CompNovoIdentificationCID_test
	CompNovoIonScoringBase_test
	CompNovoIonScoring_test
  ILPWrapper_test
	MRMFragmentSelection_test
	ProtonDistributionModel_test
	HiddenMarkovModel_test
	PILISModel_test
	PILISModelGenerator_test
)

set(applications_executables_list
  TOPPBase_test
)

set(transformations_executables_list
  ModelDescription_test
  BaseModel_test
  EmgModel_test
  LmaGaussModel_test
  ExtendedIsotopeModel_test
  InterpolationModel_test
  GaussModel_test
  BiGaussModel_test
  IsotopeModel_test
  LmaIsotopeModel_test
  ProductModel_test
  FeatureFinder_test
  FeaFiModule_test
  FeatureFinderAlgorithm_test
  FeatureFinderAlgorithmPicked_test
  FeatureFinderAlgorithmSimplest_test
  FeatureFinderAlgorithmSimple_test
  FeatureFinderAlgorithmWavelet_test
  FeatureFinderAlgorithmIsotopeWavelet_test
  FeatureFinderAlgorithmMRM_test
  SimpleSeeder_test
  SimpleExtender_test
  GaussFitter1D_test
  BiGaussFitter1D_test
  IsotopeFitter1D_test
  EmgFitter1D_test
  LmaGaussFitter1D_test
  LmaIsotopeFitter1D_test
  ExtendedIsotopeFitter1D_test
  Fitter1D_test
  LevMarqFitter1D_test
  MaxLikeliFitter1D_test
  ModelFitter_test
  IsotopeWavelet_test
  IsotopeWaveletTransform_test
  ContinuousWaveletTransform_test
  PeakShape_test
  ContinuousWaveletTransformNumIntegration_test
  OptimizePick_test
  PeakPickerCWT_test
  PeakPickerHiRes_test
  OptimizePeakDeconvolution_test
  TwoDOptimization_test
)

set(simulation_executables_list
  MSSim_test
  DigestSimulation_test
  PTMSimulation_test
  DetectabilitySimulation_test
  RTSimulation_test
  IonizationSimulation_test
  RawMSSignalSimulation_test
  ElutionModel_test
  IsotopeModelGeneral_test
  MixtureModel_test
)

### collect test executables
set(TEST_executables
		${concept_executables_list}
		${system_executables_list}
		${datastructures_executables_list}
		${kernel_executables_list}
		${metadata_executables_list}
		${visual_executables_list}
		${format_executables_list}
		${math_executables_list}
		${filtering_executables_list}
		${comparison_executables_list}
		${chemistry_executables_list}
		${analysis_executables_list}
		${applications_executables_list}
		${transformations_executables_list}
		${simulation_executables_list})
		

### add filenames to Visual Studio solution tree
set(sources_VS)
foreach(i ${TEST_executables})
	list(APPEND sources_VS "${i}.C")
endforeach(i)
source_group("" FILES ${sources_VS})
