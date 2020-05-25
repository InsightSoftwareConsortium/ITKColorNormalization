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

#ifndef itkStructurePreservingColorNormalizationFilter_h
#define itkStructurePreservingColorNormalizationFilter_h

#include "itkImageToImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkSmartPointer.h"
#include "itkeigen/Eigen/Core"

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
  using OutputSizeType = itk::Size< OutputImageType::ImageDimension >;
  using OutputSizeValueType = typename OutputSizeType::SizeValueType;
  using OutputPixelType = typename OutputImageType::PixelType;

  using CalcElementType = float;
  using CalcMatrixType = Eigen::Matrix< CalcElementType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >;
  using CalcColVectorType = Eigen::Matrix< CalcElementType, Eigen::Dynamic, 1 >;
  using CalcRowVectorType = Eigen::Matrix< CalcElementType, 1, Eigen::Dynamic >;
  using CalcDiagMatrixType = Eigen::DiagonalMatrix< CalcElementType, Eigen::Dynamic >;
  using CalcUnaryFunctionPointer = CalcElementType ( * ) ( CalcElementType );

  /** Standard class typedefs. */
  using Self = StructurePreservingColorNormalizationFilter< InputImageType, OutputImageType >;
  using Superclass = ImageToImageFilter< InputImageType, OutputImageType >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Run-time type information. */
  itkTypeMacro( StructurePreservingColorNormalizationFilter, ImageToImageFilter );

  /** Standard New macro. */
  itkNewMacro( Self );

  itkGetMacro( ColorIndexSuppressedByHematoxylin, int )
  itkSetMacro( ColorIndexSuppressedByHematoxylin, int )

  itkGetMacro( ColorIndexSuppressedByEosin, int )
  itkSetMacro( ColorIndexSuppressedByEosin, int )

  static constexpr InputSizeValueType InputImageDimension = TInputImage::ImageDimension;
  static constexpr InputSizeValueType OutputImageDimension = TOutputImage::ImageDimension;

  // Pixel length (aka number of colors) must be at least 3.  Note
  // that if the pixel value is, e.g., float or unsigned char, then it
  // is a single color (e.g., grey) and it will not have "Length"
  // defined.  In such a case, the compiler will properly give an
  // error; though it may be hard for the user to understand that the
  // issue is too few colors.
  static constexpr InputSizeValueType InputImageLength = InputPixelType::Length;
  static constexpr InputSizeValueType OutputImageLength = OutputPixelType::Length;
  static_assert( InputImageLength >= 3,
    "itkStructurePreservingColorNormalizationFilter input image needs length (#colors) >= 3." );
  static_assert( OutputImageLength == InputImageLength,
    "StructurePreservingColorNormalizationFilter output image needs length (#colors) exactly the same as input image." );

  // This algorithm is defined for H&E (Hematoxylin (blue) and Eosin
  // (pink)), which is a total of 2 stains.  However, this approach
  // could in theory work in other circumstances.  In that case it
  // might be better to have NumberOfStains be a template parameter or
  // a setable class member.
  static constexpr InputSizeValueType NumberOfStains = 2;

protected:
  StructurePreservingColorNormalizationFilter();
  ~StructurePreservingColorNormalizationFilter() override = default;

  void PrintSelf( std::ostream & os, Indent indent ) const override;

  void GenerateInputRequestedRegion() override;

  void BeforeThreadedGenerateData() override;

  void DynamicThreadedGenerateData( const OutputRegionType & outputRegion ) override;

  int ImageToNMF( InputRegionConstIterator &iter, CalcMatrixType &matrixH, InputPixelType &pixelUnstained ) const;

  void ImageToMatrix( InputRegionConstIterator &iter, CalcMatrixType &matrixV ) const;

  void MatrixToDistinguishers( const CalcMatrixType &matrixV, CalcMatrixType &distinguishers ) const;

  void MatrixToBrightPartOfMatrix( CalcMatrixType &matrixV ) const;

  void FirstPassDistinguishers( const CalcMatrixType &normVStart, std::array< int, NumberOfStains+1 > &firstPassDistinguisherIndices, InputSizeValueType &numberOfDistinguishers ) const;

  void SecondPassDistinguishers( const CalcMatrixType &normVStart, const std::array< int, NumberOfStains+1 > &firstPassDistinguisherIndices, const InputSizeValueType numberOfDistinguishers,
    CalcMatrixType &secondPassDistinguisherColors ) const;

  static int MatrixToOneDistinguisher( const CalcMatrixType &normV, const CalcColVectorType &lastOnes );

  static CalcMatrixType RecenterMatrix( const CalcMatrixType &normV, const CalcColVectorType &firstOnes, const InputSizeValueType row );

  static CalcMatrixType ProjectMatrix( const CalcMatrixType &normV, const InputSizeValueType row );

  int DistinguishersToNMFSeeds( const CalcMatrixType &distinguishers, InputPixelType &pixelUnstained, CalcMatrixType &matrixH ) const;

  void DistinguishersToColors( const CalcMatrixType &distinguishers, InputSizeValueType &unstainedIndex, InputSizeValueType &hematoxylinIndex, InputSizeValueType &eosinIndex ) const;

  void VirtanenEuclidean( const CalcMatrixType &matrixV, CalcMatrixType &matrixW, const CalcMatrixType &matrixH ) const;

  void VirtanenKLDivergence( const CalcMatrixType &matrixV, CalcMatrixType &matrixW, const CalcMatrixType &matrixH ) const;

  void NMFsToImage( const CalcMatrixType &inputH, const InputPixelType &inputUnstained, const CalcMatrixType &referH, const InputPixelType &referUnstained, OutputRegionIterator &out ) const;

  int m_ColorIndexSuppressedByHematoxylin;
  int m_ColorIndexSuppressedByEosin;

  // Our installation of Eigen3 doesn't have iterators.  (They arrive
  // with Eigen 3.4.)  We define begin, cbegin, end, and cend
  // functions here.  A compiler sometimes gets segmentation fault if
  // we use the more restrictive Eigen::Matrix< ... > declarations, so
  // we have the more lax TMatrix declarations available as
  // alternates.
#define STRUCTUREPRESERVINGCOLORNORMALIZATIONFILTER_STRICT_EIGEN3_ITERATORS 0
#if STRUCTUREPRESERVINGCOLORNORMALIZATIONFILTER_STRICT_EIGEN3_ITERATORS
  template< typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols >
  static _Scalar *begin(typename Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols > &matrix);

  template< typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols >
  static const _Scalar *cbegin(const typename Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols > &matrix);

  template< typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols >
  static _Scalar *end(typename Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols > &matrix);

  template< typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols >
  static const _Scalar *cend(const typename Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols > &matrix);

#else
  template< typename TMatrix >
  static CalcElementType *begin(TMatrix &matrix);

  template< typename TMatrix >
  static const CalcElementType *cbegin(const TMatrix &matrix);

  template< typename TMatrix >
  static CalcElementType *end(TMatrix &matrix);

  template< typename TMatrix >
  static const CalcElementType *cend(const TMatrix &matrix);
#endif

  // These members are for the purpose of caching results for use the
  // next time the pipeline is run.
  const InputImageType *m_inputPtr;
  TimeStamp m_inputTimeStamp;
  CalcMatrixType m_inputH;
  InputPixelType m_inputUnstainedPixel;
  const InputImageType *m_referPtr;
  TimeStamp m_referTimeStamp;
  CalcMatrixType m_referH;
  InputPixelType m_referUnstainedPixel;

private:
  static constexpr CalcElementType biggerEpsilon {1e-3}; // a small matrix.array_inf_norm() value
  static constexpr CalcElementType epsilon {1e-6}; // a very small matrix element
  static constexpr CalcElementType epsilon2 {epsilon * epsilon}; // a very small squared magnitude for a vector.
  static constexpr InputSizeValueType maxNumberOfIterations {0}; // For Virtanen's non-negative matrix factorization algorithm.
  static constexpr CalcElementType lambda {0.02}; // For Lasso penalty.


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
