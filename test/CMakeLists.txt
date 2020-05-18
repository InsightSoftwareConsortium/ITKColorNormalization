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
    ${ITK_TEST_OUTPUT_DIR}/itkStructurePreservingColorNormalizationFilterTestOutput.mha
  )
