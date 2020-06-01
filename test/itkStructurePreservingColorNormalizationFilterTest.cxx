/*=========================================================================
 *
 *  Copyright NumFOCUS
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
#include "itkImageFileReader.h"
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
    const auto * processObject = dynamic_cast< const itk::ProcessObject * >( caller );
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
  // At compile time, these examples should fail a static_assert and
  // report "Images need at least 3 colors".
  #if 0
  {
  using tmpImageType = typename itk::Image< float >;
  auto tmp = itk::StructurePreservingColorNormalizationFilter< tmpImageType >::New();
  }
  {
  using tmpImageType = typename itk::Image< typename itk::FixedArray< unsigned char, 2 > >;
  auto tmp = itk::StructurePreservingColorNormalizationFilter< tmpImageType >::New();
  }
  #endif
  // These examples should compile.
  {
  using tmpImageType = typename itk::Image< typename itk::RGBPixel< unsigned char > >;
  auto tmp = itk::StructurePreservingColorNormalizationFilter< tmpImageType >::New();
  }
  {
  using tmpImageType = typename itk::Image< typename itk::RGBPixel< double > >;
  auto tmp = itk::StructurePreservingColorNormalizationFilter< tmpImageType >::New();
  }
  {
  using tmpImageType = typename itk::Image< typename itk::RGBAPixel< unsigned char > >;
  auto tmp = itk::StructurePreservingColorNormalizationFilter< tmpImageType >::New();
  }
  {
  using tmpImageType = typename itk::Image< typename itk::RGBAPixel< double > >;
  auto tmp = itk::StructurePreservingColorNormalizationFilter< tmpImageType >::New();
  }
  {
  using tmpImageType = typename itk::Image< typename itk::Vector< unsigned char, 4 > >;
  auto tmp = itk::StructurePreservingColorNormalizationFilter< tmpImageType >::New();
  }
  {
  using tmpImageType = typename itk::Image< typename itk::Vector< double, 5 > >;
  auto tmp = itk::StructurePreservingColorNormalizationFilter< tmpImageType >::New();
  }

  // Run-time test
  if( argc < 4 )
  {
    std::cerr << "Missing parameters." << std::endl;
    std::cerr << "Usage: " << itkNameOfTestExecutableMacro( argv );
    std::cerr << " input0Image";
    std::cerr << " input1Image";
    std::cerr << " outputImage";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  const char * const input0ImageFileName = argv[1];
  const char * const input1ImageFileName = argv[2];
  const char * const outputImageFileName = argv[3];

  constexpr unsigned int Dimension = 2;
  using PixelType = typename itk::RGBPixel< unsigned char >;
  static constexpr unsigned int NumberOfColors = PixelType::Length;
  using ImageType = typename itk::Image< PixelType, Dimension >;

  using FilterType = itk::StructurePreservingColorNormalizationFilter< ImageType >;
  FilterType::Pointer filter = FilterType::New();

  EXERCISE_BASIC_OBJECT_METHODS( filter, StructurePreservingColorNormalizationFilter, ImageToImageFilter );

  ShowProgress::Pointer showProgress = ShowProgress::New();
  filter->AddObserver( itk::ProgressEvent(), showProgress );

  using ReaderType = itk::ImageFileReader< ImageType >;
  ReaderType::Pointer reader0 = ReaderType::New();
  reader0->SetFileName( input0ImageFileName );
  filter->SetInput( 0, reader0->GetOutput() );   // image to be normalized using ...

  ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( input1ImageFileName );
  filter->SetInput( 1, reader1->GetOutput() );   // reference image for normalization

#if 0
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

  using CalcElementType = typename FilterType::CalcElementType;
  using CalcRowVectorType = typename FilterType::CalcRowVectorType;
  using CalcUnaryFunctionPointer = typename FilterType::CalcUnaryFunctionPointer;
  CalcRowVectorType logWhite {CalcRowVectorType::Constant( 1, NumberOfColors, 1.0 )};
  CalcRowVectorType logHematoxylin {CalcRowVectorType::Constant( 1, NumberOfColors, 1.0 )};
  CalcRowVectorType logEosin {CalcRowVectorType::Constant( 1, NumberOfColors, 1.0 )};
  for( int color {0}; color < NumberOfColors; ++color )
    {
    logWhite( color ) = std::log( static_cast< CalcElementType >( white[color] ) );
    logHematoxylin( color ) = logWhite( color ) - std::log( static_cast< CalcElementType >( hematoxylin[color] ) );
    logEosin( color ) = logWhite( color ) - std::log( static_cast< CalcElementType >( eosin[color] ) );
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
      const CalcElementType hematoxylinContribution( 0.1 * ( 1.0 / uniformGenerator->GetVariate() - 1.0 ) );
      const CalcElementType eosinContribution( 0.1 * ( 1.0 / uniformGenerator->GetVariate() - 1.0 ) );
      const CalcElementType noise( 5.0 * normalGenerator->GetVariate() );
      const CalcRowVectorType randomPixelValue
        {( logWhite - ( hematoxylinContribution * logHematoxylin ) - ( eosinContribution * logEosin ) ).unaryExpr( CalcUnaryFunctionPointer( std::exp ) ).array() + noise};
      // std::cout << "hematoxylinContribution = " << hematoxylinContribution << ", eosinContribution = " << eosinContribution << ", randomPixelValue = " << randomPixelValue << std::endl;
      for( int color {0}; color < NumberOfColors; ++color )
        {
        tmp[color] = std::max( CalcElementType( 0.0 ), std::min( CalcElementType( 255.0 ), randomPixelValue( color ) ) );
        }
      iter.Set( tmp );
      }
    }

  filter->SetInput( 0, input );   // image to be normalized using ...
  filter->SetInput( 1, refer );   // reference image
#endif

  using WriterType = itk::ImageFileWriter< ImageType >;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputImageFileName );
  writer->SetInput( filter->GetOutput() );
  writer->SetUseCompression( true );

  TRY_EXPECT_NO_EXCEPTION( writer->Update() );


  std::cout << "Test finished." << std::endl;
  return EXIT_SUCCESS;
}
