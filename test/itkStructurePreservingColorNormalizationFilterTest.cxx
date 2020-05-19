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
#include "itkImageFileWriter.h"
#include "itkTestingMacros.h"

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
  using ImageType = typename itk::Image< PixelType, Dimension >;

  using FilterType = itk::StructurePreservingColorNormalizationFilter< ImageType, ImageType >;
  FilterType::Pointer filter = FilterType::New();

  EXERCISE_BASIC_OBJECT_METHODS( filter, StructurePreservingColorNormalizationFilter, ImageToImageFilter );

  // Create input image to avoid test dependencies.
  ImageType::SizeType size;
  size.Fill( 128 );
  ImageType::Pointer image = ImageType::New();
  image->SetRegions( size );
  image->Allocate();
  PixelType white;
  white.SetRed( 255 );
  white.SetGreen( 255 );
  white.SetBlue( 255 );
  PixelType hematoxylin;        // dominant effect is to suppress red
  hematoxylin.SetRed( 16 );
  hematoxylin.SetGreen( 67 );
  hematoxylin.SetBlue( 118 );
  PixelType eosin;              // dominant effect is to suppress green
  eosin.SetRed( 199 );
  eosin.SetGreen( 21 );
  eosin.SetBlue( 133 );

  image->FillBuffer( white );

  ImageType::IndexValueType coordinates[] { 0, 0 };
  ImageType::IndexType index;
  for( ImageType::IndexValueType i = ImageType::IndexValueType( 0 ); i < ImageType::IndexValueType( 128 ); ++i )
    {
    coordinates[1] = i;
    coordinates[0] = 0;
    index.SetIndex( coordinates );
    image->SetPixel( index, hematoxylin );
    coordinates[0] = 1;
    index.SetIndex( coordinates );
    image->SetPixel( index, eosin );
    coordinates[0] = 2;
    index.SetIndex( coordinates );
    image->SetPixel( index, ( hematoxylin + eosin ) / 2 - white / 100 );
    }

  ShowProgress::Pointer showProgress = ShowProgress::New();
  filter->AddObserver( itk::ProgressEvent(), showProgress );
  filter->SetInput( 0, image );   // image to be normalized using ...
  filter->SetInput( 1, image );   // reference image
  filter->GetOutput()->SetRequestedRegion( filter->GetInput( 0 )->GetRequestedRegion().GetSize() );

  using WriterType = itk::ImageFileWriter< ImageType >;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputImageFileName );
  writer->SetInput( filter->GetOutput() );
  writer->SetUseCompression( true );

  TRY_EXPECT_NO_EXCEPTION( writer->Update() );


  std::cout << "Test finished." << std::endl;
  return EXIT_SUCCESS;
}
