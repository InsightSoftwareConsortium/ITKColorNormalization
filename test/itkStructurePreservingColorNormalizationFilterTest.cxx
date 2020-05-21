/*=========================================================================
 *
 *  Copyright NumFOCUS!!!
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "itkStructurePreservingColorNormalizationFilter.h"

#include "itkCommand.h"
#include "itkImageFileWriter.h"
#include "itkTestingMacros.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkNormalVariateGenerator.h"

namespace
{

class ShowProgress : public itk::Command
{
public:
  itkNewMacro( ShowProgress );

  void
  Execute( itk::Object * caller, const itk::EventObject & event ) override
  {
    Execute( ( const itk::Object * )caller, event );
  }

  void
  Execute( const itk::Object * caller, const itk::EventObject & event ) override
  {
    if( !itk::ProgressEvent().CheckEvent( &event ) )
    {
      return;
    }
    const auto * processObject = dynamic_cast< const itk::ProcessObject *>( caller );
    if( !processObject )
    {
      return;
    }
    std::cout << " " << processObject->GetProgress();
  }
};

} // namespace

int itkStructurePreservingColorNormalizationFilterTest( int argc, char * argv[] )
{
// Install correct output (to be compared to during test) to
// ~/git/ITKStructurePreservingColorNormalization-build/ExternalData/test/Baseline/itkStructurePreservingColorNormalizationFilterTestOutput.mha
// or similar location!!!

  if( argc < 2 )
  {
    std::cerr << "Missing parameters." << std::endl;
    std::cerr << "Usage: " << itkNameOfTestExecutableMacro( argv );
    std::cerr << " outputImage";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }
  const char * const outputImageFileName = argv[1];

  constexpr unsigned int Dimension = 2;
  using PixelType = typename itk::RGBPixel< unsigned char >;
  static constexpr unsigned int InputImageLength = PixelType::Length;
  using ImageType = typename itk::Image< PixelType, Dimension >;

  using FilterType = itk::StructurePreservingColorNormalizationFilter< ImageType, ImageType >;
  FilterType::Pointer filter = FilterType::New();

  EXERCISE_BASIC_OBJECT_METHODS( filter, StructurePreservingColorNormalizationFilter, ImageToImageFilter );

  // Create input images to avoid test dependencies.
  const ImageType::SizeValueType testSize = 128;
  ImageType::SizeType size;
  size.Fill( testSize );
  ImageType::Pointer input = ImageType::New();
  ImageType::Pointer refer = ImageType::New();

  // We will need some random number generators.
  using UniformGeneratorType = itk::Statistics::MersenneTwisterRandomVariateGenerator;
  UniformGeneratorType::Pointer uniformGenerator = UniformGeneratorType::New();
  uniformGenerator->Initialize( 20200519 );

  using NormalGeneratorType = itk::Statistics::NormalVariateGenerator;
  NormalGeneratorType::Pointer normalGenerator = NormalGeneratorType::New();
  normalGenerator->Initialize( 20200520 );

  // Define some useful colors
  PixelType white;              // For unstained / background pixels
  white.SetRed( 240 );
  white.SetGreen( 240 );
  white.SetBlue( 240 );
  PixelType hematoxylin;        // dominant effect is to suppress red
  hematoxylin.SetRed( 16 );
  hematoxylin.SetGreen( 67 );
  hematoxylin.SetBlue( 118 );
  PixelType eosin;              // dominant effect is to suppress green
  eosin.SetRed( 199 );
  eosin.SetGreen( 21 );
  eosin.SetBlue( 133 );

  using CalcElementType = double;
  using CalcVectorType = vnl_vector< CalcElementType >;
  CalcVectorType logWhite {InputImageLength};
  CalcVectorType logHematoxylin {InputImageLength};
  CalcVectorType logEosin {InputImageLength};
  for( int color {0}; color < InputImageLength; ++color )
    {
    logWhite.put( color, std::log( static_cast< CalcElementType >( white[color] ) ) );
    logHematoxylin.put( color, logWhite.get( color ) - std::log( static_cast< CalcElementType >( hematoxylin[color] ) ) );
    logEosin.put( color, logWhite.get( color ) - std::log( static_cast< CalcElementType >( eosin[color] ) ) );
    }

  // { std::ostringstream mesg; mesg << "logHematoxylin = " << logHematoxylin << ", logEosin = " << logEosin << std::endl; std::cout << mesg.str(); }

  // Randomly generate both input images
  for( ImageType::Pointer image : {input, refer} )
    {
    image->SetRegions( size );
    image->Allocate();
    image->FillBuffer( white );

    using InputRegionIterator = typename itk::ImageRegionIterator< ImageType >;
    InputRegionIterator iter {image, size};
    PixelType tmp;
    for( iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
      {
      const double hematoxylinContribution {0.1 * ( 1.0 / uniformGenerator->GetVariate() - 1.0 )};
      const double eosinContribution {0.1 * ( 1.0 / uniformGenerator->GetVariate() - 1.0 )};
      const double noise {5.0 * normalGenerator->GetVariate()};
      const CalcVectorType randomPixelValue {( logWhite - ( hematoxylinContribution * logHematoxylin ) - ( eosinContribution * logEosin ) ).apply( std::exp ) + noise};
      // std::cout << "hematoxylinContribution = " << hematoxylinContribution << ", eosinContribution = " << eosinContribution << ", randomPixelValue = " << randomPixelValue << std::endl;
      for( int color {0}; color < InputImageLength; ++color )
        {
        tmp[color] = std::max( 0.0, std::min( 255.0, randomPixelValue.get( color ) ) );
        }
      iter.Set( tmp );
      }
    }

  ShowProgress::Pointer showProgress = ShowProgress::New();
  filter->AddObserver( itk::ProgressEvent(), showProgress );
  filter->SetInput( 0, input );   // image to be normalized using ...
  filter->SetInput( 1, refer );   // reference image

  using WriterType = itk::ImageFileWriter< ImageType >;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputImageFileName );
  writer->SetInput( filter->GetOutput() );
  writer->SetUseCompression( true );

  TRY_EXPECT_NO_EXCEPTION( writer->Update() );


  std::cout << "Test finished." << std::endl;
  return EXIT_SUCCESS;
}
