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

#ifndef itkStructurePreservingColorNormalizationFilter_h
#define itkStructurePreservingColorNormalizationFilter_h

#include "itkImageToImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkSmartPointer.h"

namespace itk
{

/** \class StructurePreservingColorNormalizationFilter
 *
 * \brief Filters a image by iterating over its pixels.
 *
 * Filters a image by iterating over its pixels in a multi-threaded way
 * and {to be completed by the developer}.
 *
 * \ingroup StructurePreservingColorNormalization
 *
 */
template< typename TInputImage, typename TOutputImage = TInputImage >
class StructurePreservingColorNormalizationFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN( StructurePreservingColorNormalizationFilter );

  using InputImageType = TInputImage;
  using InputRegionType = typename InputImageType::RegionType;
  using InputRegionConstIterator = typename itk::ImageRegionConstIterator< InputImageType >;
  using InputSizeType = itk::Size< InputImageType::ImageDimension >;
  using InputSizeValueType = typename InputSizeType::SizeValueType;
  using InputPixelType = typename InputImageType::PixelType;

  using OutputImageType = TOutputImage;
  using OutputRegionType = typename OutputImageType::RegionType;
  using OutputRegionIterator = typename itk::ImageRegionIterator< OutputImageType >;
  using OutputPixelType = typename OutputImageType::PixelType;

  using CalcElementType = double;
  using CalcMatrixType = vnl_matrix< CalcElementType >; // Use Eigen instead of VNL?!!!
  using CalcVectorType = vnl_vector< CalcElementType >;
  using CalcDiagMatrixType = vnl_diag_matrix< CalcElementType >;

  /** Standard class typedefs. */
  using Self = StructurePreservingColorNormalizationFilter< InputImageType, OutputImageType >;
  using Superclass = ImageToImageFilter< InputImageType, OutputImageType >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Run-time type information. */
  itkTypeMacro( StructurePreservingColorNormalizationFilter, ImageToImageFilter );

  /** Standard New macro. */
  itkNewMacro( Self );

  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int OutputImageDimension = TOutputImage::ImageDimension;

  // Pixel length (aka number of colors) must be at least 3.  Note
  // that if the pixel value is, e.g., float or unsigned char, then it
  // is a single color (e.g., grey) and it will not have "Length"
  // defined.  In such a case, the compiler will properly give an
  // error; though it may be hard for the user to understand that the
  // issue is too few colors.
  static constexpr unsigned int InputImageLength = InputPixelType::Length;
  static constexpr unsigned int OutputImageLength = OutputPixelType::Length;
  static_assert( InputImageLength >= 3,
    "itkStructurePreservingColorNormalizationFilter input image needs length (#colors) >= 3." );
  static_assert( OutputImageLength == InputImageLength,
    "StructurePreservingColorNormalizationFilter output image needs length (#colors) exactly the same as input image." );

  // This algorithm is defined for H&E (Hematoxylin (blue) and Eosin
  // (pink)), which is a total of 2 stains.  However, this approach
  // could in theory work in other circumstances.  In that case it
  // might be better to have NumberOfStains be a template parameter or
  // a setable class member.
  static constexpr unsigned int NumberOfStains = 2;

protected:
  StructurePreservingColorNormalizationFilter();
  ~StructurePreservingColorNormalizationFilter() override = default;

  void PrintSelf( std::ostream & os, Indent indent ) const override;

  void GenerateInputRequestedRegion() override;

  void BeforeThreadedGenerateData() override;

  void DynamicThreadedGenerateData( const OutputRegionType & outputRegion ) override;

  void ImageToNMF( InputRegionConstIterator &iter, CalcMatrixType &matrixV, CalcMatrixType &matrixW, CalcMatrixType &matrixH, InputPixelType &pixelUnstained ) const;

  void ImageToMatrix( InputRegionConstIterator &iter, CalcMatrixType &matrixV ) const;

  void MatrixToDistinguishers( const CalcMatrixType &matrixV, CalcMatrixType &distinguishers ) const;

  CalcMatrixType MatrixToBrightPartOfMatrix( const CalcMatrixType &matrixV ) const;

  void FirstPassDistinguishers( const CalcMatrixType &normVStart, std::array< int, NumberOfStains+1 > &firstPassDistinguisherIndices, unsigned int &numberOfDistinguishers ) const;

  void SecondPassDistinguishers( const CalcMatrixType &normVStart, const std::array< int, NumberOfStains+1 > &firstPassDistinguisherIndices, const int numberOfDistinguishers,
    const CalcMatrixType &brightV, CalcMatrixType &secondPassDistinguisherColors ) const;

  int MatrixToOneDistinguisher( const CalcMatrixType &normV, const CalcVectorType &lastOnes ) const;

  CalcMatrixType RecenterMatrix( const CalcMatrixType &normV, const CalcVectorType &firstOnes, const int row ) const;

  CalcMatrixType ProjectMatrix( const CalcMatrixType &normV, const int row ) const;

  void DistinguishersToNMFSeeds( const CalcMatrixType &distinguishers, InputPixelType &pixelUnstained, CalcMatrixType &matrixV, CalcMatrixType &matrixW,
    CalcMatrixType &matrixH ) const;

  void DistinguishersToColors( const CalcMatrixType &distinguishers, int &unstainedIndex, int &hematoxylinIndex, int &eosinIndex ) const;

  void VirtanenEuclidean( const CalcMatrixType &matrixV, CalcMatrixType &matrixW, CalcMatrixType &matrixH ) const;

  void VirtanenKLDivergence( const CalcMatrixType &matrixV, CalcMatrixType &matrixW, CalcMatrixType &matrixH ) const;

  void NMFsToImage( const CalcMatrixType &inputW, const CalcMatrixType &inputH, const CalcMatrixType &referH, const InputPixelType &referUnstained, OutputRegionIterator &out ) const;

private:
  static constexpr CalcElementType epsilon {1e-6};   // a very small matrix element
  static constexpr CalcElementType epsilon2 {epsilon * epsilon}; // a very small squared magnitude for a vector.
  static constexpr unsigned int    maxNumberOfIterations {10000u}; // For Virtanen's non-negative matrix factorization algorithm.
  static constexpr CalcElementType lambda {0.02}; // For Lasso penalty.

  const InputImageType *m_inputPtr;
  TimeStamp m_inputTimeStamp;
  CalcMatrixType m_inputW;
  CalcMatrixType m_inputH;
  const InputImageType *m_referPtr;
  TimeStamp m_referTimeStamp;
  CalcMatrixType m_referH;
  InputPixelType m_referUnstainedPixel;


#ifdef ITK_USE_CONCEPT_CHECKING
  // Add concept checking such as
  // itkConceptMacro( FloatingPointPixel, ( itk::Concept::IsFloatingPoint< typename InputImageType::PixelType > ) );
#endif
};

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkStructurePreservingColorNormalizationFilter.hxx"
#endif

#endif // itkStructurePreservingColorNormalizationFilter
