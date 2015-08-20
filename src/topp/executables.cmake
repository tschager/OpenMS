### the directory name
set(directory source/APPLICATIONS/TOPP)

### list all filenames of the directory here
set(TOPP_executables
NoiseFilterGaussian
NoiseFilterSGolay
OMSSAAdapter
PeakPickerHiRes
PeakPickerWavelet
)

## all targets requiring OpenMS_GUI
set(TOPP_executables_with_GUIlib
ExecutePipeline
Resampler
)

### add filenames to Visual Studio solution tree
set(sources_VS)
foreach(i ${TOPP_executables} ${TOPP_executables_with_GUIlib})
	list(APPEND sources_VS "${i}.cpp")
endforeach(i)

source_group("" FILES ${sources_VS})
