itk_module_test()

set(StructurePreservingColorNormalizationTests
  itkStructurePreservingColorNormalizationFilterTest.cxx
  )

CreateTestDriver(StructurePreservingColorNormalization "${StructurePreservingColorNormalization-Test_LIBRARIES}" "${StructurePreservingColorNormalizationTests}")

itk_add_test(NAME itkStructurePreservingColorNormalizationFilterTest
  COMMAND StructurePreservingColorNormalizationTestDriver
    --compare
    ${ITK_TEST_OUTPUT_DIR}/itkStructurePreservingColorNormalizationFilterTestOutput.png
    DATA{Baseline/itkStructurePreservingColorNormalizationFilterTestOutput.png}
    --compareNumberOfPixelsTolerance 2
  itkStructurePreservingColorNormalizationFilterTest
    DATA{Baseline/itkStructurePreservingColorNormalizationFilterTestInput0.png}
    DATA{Baseline/itkStructurePreservingColorNormalizationFilterTestInput1.png}
    ${ITK_TEST_OUTPUT_DIR}/itkStructurePreservingColorNormalizationFilterTestOutput.png
  )
