# the top-level README is used for describing this module, just
# re-used it for documentation here
get_filename_component(MY_CURRENT_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
file(READ "${MY_CURRENT_DIR}/README.md" DOCUMENTATION)

# itk_module() defines the module dependencies in StructurePreservingColorNormalization
# StructurePreservingColorNormalization depends on ITKCommon
# The testing module in StructurePreservingColorNormalization depends on ITKTestKernel
# and ITKMetaIO(besides StructurePreservingColorNormalization and ITKCore)
# By convention those modules outside of ITK are not prefixed with
# ITK.

# define the dependencies of the include module and the tests
itk_module(StructurePreservingColorNormalization
  DEPENDS
    ITKCommon
    ITKStatistics
  COMPILE_DEPENDS
    ITKImageSources
  TEST_DEPENDS
    ITKTestKernel
    ITKMetaIO
  DESCRIPTION
    "${DOCUMENTATION}"
  EXCLUDE_FROM_DEFAULT
  ENABLE_SHARED
)
