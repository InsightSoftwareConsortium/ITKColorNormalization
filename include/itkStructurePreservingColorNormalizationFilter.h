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

#include <type_traits>
#include "itkRGBPixel.h"
#include "itkRGBAPixel.h"
#include "itkVector.h"
#include "itkCovariantVector.h"
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
template< typename TImage >
class StructurePreservingColorNormalizationFilter : public ImageToImageFilter< TImage, TImage >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN( StructurePreservingColorNormalizationFilter );

  using ImageType = TImage;
  using RegionType = typename ImageType::RegionType;
  using RegionConstIterator = typename itk::ImageRegionConstIterator< ImageType >;
  using RegionIterator = typename itk::ImageRegionIterator< ImageType >;
  using SizeType = itk::Size< ImageType::ImageDimension >;
  using SizeValueType = typename SizeType::SizeValueType;
  using PixelType = typename ImageType::PixelType;

  using CalcElementType = double;
  using CalcMatrixType = Eigen::Matrix< CalcElementType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >;
  using CalcColVectorType = Eigen::Matrix< CalcElementType, Eigen::Dynamic, 1 >;
  using CalcRowVectorType = Eigen::Matrix< CalcElementType, 1, Eigen::Dynamic >;
  using CalcUnaryFunctionPointer = CalcElementType ( * ) ( CalcElementType );

  /** Standard class typedefs. */
  using Self = StructurePreservingColorNormalizationFilter< ImageType >;
  using Superclass = ImageToImageFilter< ImageType, ImageType >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Run-time type information. */
  itkTypeMacro( StructurePreservingColorNormalizationFilter, ImageToImageFilter );

  /** Standard New macro. */
  itkNewMacro( Self );

  itkGetMacro( ColorIndexSuppressedByHematoxylin, Eigen::Index )
  itkSetMacro( ColorIndexSuppressedByHematoxylin, Eigen::Index )

  itkGetMacro( ColorIndexSuppressedByEosin, Eigen::Index )
  itkSetMacro( ColorIndexSuppressedByEosin, Eigen::Index )

  // This algorithm is defined for H&E (Hematoxylin (blue) and
  // Eosin (pink)), which is a total of 2 stains.  However, this
  // approach could in theory work in other circumstances.  In that
  // case it might be better to have NumberOfStains be a template
  // parameter or a setable class member.
  static constexpr SizeValueType NumberOfStains {2};
  static constexpr CalcElementType BrightPercentileLevel {0.80};
  static constexpr CalcElementType BrightPercentageLevel {0.70};
  static constexpr CalcElementType VeryDarkPercentileLevel {0.01};

  // We have special cases for different pixel types, including: (1)
  // we refuse to process fewer than 3 colors (at compile time where
  // possible) and (2) RBGA pixels have 4 dimensions but only 3
  // colors.  We accomplish the special cases via the PixelHelper
  // class, similar in concept to the itk::PixelTraits class, but with
  // different capabilities.
  //
  // The default: for the case that the number of colors is not
  // determined at compile time, such as with a VectorImage:
  template< typename TPixelType, typename = void >
  struct PixelHelper
    {
    using PixelType = TPixelType;
    using ValueType = typename PixelType::ValueType;
    static constexpr SizeValueType NumberOfDimensions = -1;
    static constexpr SizeValueType NumberOfColors = -1;
    static constexpr typename Eigen::Index ColorIndexSuppressedByHematoxylin = -1;
    static constexpr typename Eigen::Index ColorIndexSuppressedByEosin = -1;
    static PixelType pixelInstance( unsigned numberOfDimensions ) { return PixelType {numberOfDimensions}; }
    };
  // For the case that the number of colors is implicitly set to 1 at
  // compile time.  We will refuse to compile this case via a
  // static_assert in the constructor for
  // StructurePreservingColorNormalizationFilter.
  template< typename TPixelType >
  struct PixelHelper< TPixelType, typename std::enable_if< std::is_arithmetic< TPixelType >::value >::type >
    {
    using PixelType = TPixelType;
    using ValueType = PixelType;
    static constexpr SizeValueType NumberOfDimensions = 1;
    static constexpr SizeValueType NumberOfColors = 1;
    static constexpr typename Eigen::Index ColorIndexSuppressedByHematoxylin = -1;
    static constexpr typename Eigen::Index ColorIndexSuppressedByEosin = -1;
    static PixelType pixelInstance( unsigned numberOfDimensions ) { return PixelType {}; }
    };
  // For the case that the pixel type is itk::RGBPixel:
  template< typename TScalar>
  struct PixelHelper< itk::RGBPixel< TScalar >, void >
    {
    using PixelType = itk::RGBPixel< TScalar >;
    using ValueType = typename PixelType::ValueType;
    static constexpr SizeValueType NumberOfDimensions = 3;
    static constexpr SizeValueType NumberOfColors = 3;
    static constexpr typename Eigen::Index ColorIndexSuppressedByHematoxylin = 0;
    static constexpr typename Eigen::Index ColorIndexSuppressedByEosin = 1;
    static PixelType pixelInstance( unsigned numberOfDimensions ) { return PixelType {}; }
    };
  // For the case that the pixel type is itk::RGBAPixel:
  template< typename TScalar>
  struct PixelHelper< itk::RGBAPixel< TScalar >, void >
    {
    using PixelType = itk::RGBAPixel< TScalar >;
    using ValueType = typename PixelType::ValueType;
    static constexpr SizeValueType NumberOfDimensions = 4;
    static constexpr SizeValueType NumberOfColors = 3;
    static constexpr typename Eigen::Index ColorIndexSuppressedByHematoxylin = 0;
    static constexpr typename Eigen::Index ColorIndexSuppressedByEosin = 1;
    static PixelType pixelInstance( unsigned numberOfDimensions ) { return PixelType {}; }
    };
  // For the cases that the pixel type is itk::Vector or
  // itk::CovariantVector.  If NVectorDimension is not at least 3 we
  // will refuse to compile this case via a static_assert in the
  // constructor of StructurePreservingColorNormalizationFilter.
  template< typename TScalar, unsigned int NVectorDimension >
  struct PixelHelper< itk::Vector< TScalar, NVectorDimension >, void >
    {
    using PixelType = itk::Vector< TScalar, NVectorDimension >;
    using ValueType = typename PixelType::ValueType;
    static constexpr SizeValueType NumberOfDimensions = NVectorDimension;
    static constexpr SizeValueType NumberOfColors = NVectorDimension;
    static constexpr typename Eigen::Index ColorIndexSuppressedByHematoxylin = -1;
    static constexpr typename Eigen::Index ColorIndexSuppressedByEosin = -1;
    static PixelType pixelInstance( unsigned numberOfDimensions ) { return PixelType {}; }
    };
  template< typename TScalar, unsigned int NVectorDimension >
  struct PixelHelper< itk::CovariantVector< TScalar, NVectorDimension >, void >
    {
    using PixelType = itk::CovariantVector< TScalar, NVectorDimension >;
    using ValueType = typename PixelType::ValueType;
    static constexpr SizeValueType NumberOfDimensions = NVectorDimension;
    static constexpr SizeValueType NumberOfColors = NVectorDimension;
    static constexpr typename Eigen::Index ColorIndexSuppressedByHematoxylin = -1;
    static constexpr typename Eigen::Index ColorIndexSuppressedByEosin = -1;
    static PixelType pixelInstance( unsigned numberOfDimensions ) { return PixelType {}; }
    };

  using PixelValueType = typename PixelHelper< PixelType >::ValueType;

protected:

  StructurePreservingColorNormalizationFilter();
  ~StructurePreservingColorNormalizationFilter() override = default;

  void PrintSelf( std::ostream & os, Indent indent ) const override;

  void VerifyInputInformation() const override;

  void GenerateInputRequestedRegion() override;

  void BeforeThreadedGenerateData() override;

  void DynamicThreadedGenerateData( const RegionType & outputRegion ) override;

  int ImageToNMF( RegionConstIterator &iter, CalcMatrixType &matrixH, CalcRowVectorType &unstainedPixel ) const;

  void ImageToMatrix( RegionConstIterator &inIter, SizeValueType numberOfPixels, CalcMatrixType &matrixBrightV, CalcMatrixType &matrixDarkV ) const;

  static void MatrixToDistinguishers( const CalcMatrixType &matrixV, CalcMatrixType &distinguishers );

  static void MatrixToMatrixExtremes( const CalcMatrixType &matrixV, CalcMatrixType &matrixBrightV, CalcMatrixType &matrixDarkV );

  static void FirstPassDistinguishers( const CalcMatrixType &normVStart, std::array< int, NumberOfStains+1 > &firstPassDistinguisherIndices, SizeValueType &numberOfDistinguishers );

  static void SecondPassDistinguishers( const CalcMatrixType &normVStart, const std::array< int, NumberOfStains+1 > &firstPassDistinguisherIndices, const SizeValueType numberOfDistinguishers,
    CalcMatrixType &secondPassDistinguisherColors );

  static int MatrixToOneDistinguisher( const CalcMatrixType &normV );

  static CalcMatrixType RecenterMatrix( const CalcMatrixType &normV, const SizeValueType row );

  static CalcMatrixType ProjectMatrix( const CalcMatrixType &normV, const SizeValueType row );

  int DistinguishersToNMFSeeds( const CalcMatrixType &distinguishers, CalcRowVectorType &unstainedPixel, CalcMatrixType &matrixH ) const;

  void DistinguishersToColors( const CalcMatrixType &distinguishers, SizeValueType &unstainedIndex, SizeValueType &hematoxylinIndex, SizeValueType &eosinIndex ) const;

  void NormalizeMatrixH( const CalcMatrixType &matrixDarkV, const CalcRowVectorType &unstainedPixel, CalcMatrixType &matrixH ) const;

  static void VirtanenEuclidean( const CalcMatrixType &matrixV, CalcMatrixType &matrixW, CalcMatrixType &matrixH );

  static void VirtanenKLDivergence( const CalcMatrixType &matrixV, CalcMatrixType &matrixW, CalcMatrixType &matrixH );

  void NMFsToImage( const CalcMatrixType &inputH, const CalcRowVectorType &inputUnstained, const CalcMatrixType &referH, const CalcRowVectorType &referUnstained, RegionIterator &out ) const;

  // Our installation of Eigen3 does not have iterators.  (They
  // arrive with Eigen 3.4.)  We define begin, cbegin, end, and cend
  // functions here.  A compiler sometimes gets segmentation fault if
  // we use the more restrictive Eigen::Matrix< ... > declarations, so
  // we have the more lax TMatrix declarations available as
  // alternates.
#define STRUCTUREPRESERVINGCOLORNORMALIZATIONFILTER_STRICT_EIGEN3_ITERATORS 0
#if STRUCTUREPRESERVINGCOLORNORMALIZATIONFILTER_STRICT_EIGEN3_ITERATORS
  template< typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols >
  static _Scalar *begin( Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols > &matrix );

  template< typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols >
  static const _Scalar *cbegin( const Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols > &matrix );

  template< typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols >
  static _Scalar *end( Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols > &matrix );

  template< typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols >
  static const _Scalar *cend( const Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols > &matrix );

#else
  template< typename TMatrix >
  static CalcElementType *begin( TMatrix &matrix );

  template< typename TMatrix >
  static const CalcElementType *cbegin( const TMatrix &matrix );

  template< typename TMatrix >
  static CalcElementType *end( TMatrix &matrix );

  template< typename TMatrix >
  static const CalcElementType *cend( const TMatrix &matrix );
#endif

  // These members are for the purpose of caching results for use the
  // next time the pipeline is run.
  itk::ModifiedTimeType m_ParametersMTime;
  const ImageType *m_inputPtr;
  TimeStamp m_inputTimeStamp;
  CalcMatrixType m_inputH;
  CalcRowVectorType m_inputUnstainedPixel;
  const ImageType *m_referPtr;
  TimeStamp m_referTimeStamp;
  CalcMatrixType m_referH;
  CalcRowVectorType m_referUnstainedPixel;

  Eigen::Index m_NumberOfDimensions;
  Eigen::Index m_NumberOfColors;
  Eigen::Index m_ColorIndexSuppressedByHematoxylin;
  Eigen::Index m_ColorIndexSuppressedByEosin;

private:
  static constexpr CalcElementType epsilon0 {1e-3}; // a small matrix.array_inf_norm() value
  static constexpr CalcElementType epsilon1 {1e-6}; // a very small matrix element
  static constexpr CalcElementType epsilon2 {1e-12}; // a very small squared magnitude for a vector.
  static constexpr CalcElementType lambda {0.01}; // For Lasso penalty.
  static constexpr SizeValueType maxNumberOfIterations {0}; // For Virtanen's non-negative matrix factorization algorithm.
  static constexpr SizeValueType maxNumberOfRows {100000}; // Select a subset of the pixels if the image has more than this


#ifdef ITK_USE_CONCEPT_CHECKING
  // Add concept checking such as
  // itkConceptMacro( FloatingPointPixel, ( itk::Concept::IsFloatingPoint< typename ImageType::PixelType > ) );
#endif
};

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkStructurePreservingColorNormalizationFilter.hxx"
#endif

#endif // itkStructurePreservingColorNormalizationFilter
