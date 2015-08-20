### the directory name
set(directory source/APPLICATIONS/UTILS)

### list all filenames of the directory here
set(UTILS_executables
DecoyDatabase
DeMeanderize
FuzzyDiff
)

## all targets requiring OpenMS_GUI
set(UTILS_executables_with_GUIlib
IDEvaluator
ImageCreator
INIUpdater
)

### add filenames to Visual Studio solution tree
set(sources_VS)
foreach(i ${UTILS_executables} ${UTILS_executables_with_GUIlib})
	list(APPEND sources_VS "${i}.cpp")
endforeach(i)

source_group("" FILES ${sources_VS})
