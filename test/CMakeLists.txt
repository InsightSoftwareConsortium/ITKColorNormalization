itk_module_test()

set(StructurePreservingColorNormalizationTests
  itkMinimalStandardRandomVariateGeneratorTest.cxx
  itkStructurePreservingColorNormalizationFilterTest.cxx
  )

CreateTestDriver(StructurePreservingColorNormalization "${StructurePreservingColorNormalization-Test_LIBRARIES}" "${StructurePreservingColorNormalizationTests}")

itk_add_test(NAME itkMinimalStandardRandomVariateGeneratorTest
  COMMAND StructurePreservingColorNormalizationTestDriver itkMinimalStandardRandomVariateGeneratorTest
  )

itk_add_test(NAME itkStructurePreservingColorNormalizationFilterTest
  COMMAND StructurePreservingColorNormalizationTestDriver
    --compare
    ${ITK_TEST_OUTPUT_DIR}/itkStructurePreservingColorNormalizationFilterTestOutput.mha
    DATA{Baseline/itkStructurePreservingColorNormalizationFilterTestOutput.mha}
  itkStructurePreservingColorNormalizationFilterTest
    DATA{Baseline/itkStructurePreservingColorNormalizationFilterTestInput0.png}
    DATA{Baseline/itkStructurePreservingColorNormalizationFilterTestInput1.png}
    ${ITK_TEST_OUTPUT_DIR}/itkStructurePreservingColorNormalizationFilterTestOutput.mha
  )
