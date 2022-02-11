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
 * \brief This filter performs "Structure Preserving Color
 * Normalization" on an H&E image using a reference image.
 *
 * H&E (hematoxylin and eosin) are stains used to color parts of cells
 * in a histological image, often for medical diagnosis. Hematoxylin
 * is a compound that stains cell nuclei a purple-blue color. Eosin is
 * a compound that stains extracellular matrix and cytoplasm
 * pink. However, the exact color of purple-blue or pink can vary from
 * image to image, and this can make comparison of images
 * difficult. This routine addresses the issue by re-coloring one
 * image (the first image supplied to the routine) using the color
 * scheme of a reference image (the second image supplied to the
 * routine).
 *
 * Structure Preserving Color Normalization is a technique described
 * in [VPSAWBSSEN2016] and modified in [RAS2019]. The idea is to model
 * the color of an image pixel as something close to pure white, which
 * is reduced in intensity in a color-specific way via an optical
 * absorption model that depends upon the amounts of hematoxylin and
 * eosin that are present. Non-negative matrix factorization is used
 * on each analyzed image to simultaneously derive the amount of
 * hematoxylin and eosin stain at each pixel and the effective colors
 * of each stain.
 *
 * The implementation here accelerates the non-negative matrix
 * factorization by choosing the initial estimate for the color
 * absorption characteristics using a technique mimicking that
 * presented in [AGHMMSWZ2013] and [NCKZ2018]. This approach finds a
 * good solution for a non-negative matrix factorization by first
 * transforming it to the problem of finding a convex hull for a set
 * of points in a cloud.
 *
 * \ingroup StructurePreservingColorNormalization
 *
 */
template <typename TImage>
class StructurePreservingColorNormalizationFilter : public ImageToImageFilter<TImage, TImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(StructurePreservingColorNormalizationFilter);

  /** Specific class typedefs */
  using ImageType = TImage;
  using RegionType = typename ImageType::RegionType;
  using RegionConstIterator = ImageRegionConstIterator<ImageType>;
  using RegionIterator = ImageRegionIterator<ImageType>;
  using SizeType = Size<ImageType::ImageDimension>;
  using SizeValueType = typename SizeType::SizeValueType;
  using PixelType = typename ImageType::PixelType;

  using CalcElementType = double;
  using CalcMatrixType = Eigen::Matrix<CalcElementType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using CalcColVectorType = Eigen::Matrix<CalcElementType, Eigen::Dynamic, 1>;
  using CalcRowVectorType = Eigen::Matrix<CalcElementType, 1, Eigen::Dynamic>;
  using CalcUnaryFunctionPointer = CalcElementType (*)(CalcElementType);

  /** Standard class typedefs. */
  using Self = StructurePreservingColorNormalizationFilter<ImageType>;
  using Superclass = ImageToImageFilter<ImageType, ImageType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information. */
  itkTypeMacro(StructurePreservingColorNormalizationFilter, ImageToImageFilter);

  /** Standard New macro. */
  itkNewMacro(Self);

  /** If the pixel type is RGB or RGBA then
   * ColorIndexSuppressedByHematoxylin defaults to 0, indicating red.
   * Otherwise, set ColorIndexSuppressedByHematoxylin to the index in
   * the array of colors that indicates the color most suppressed by
   * hematoxylin. */
  itkGetMacro(ColorIndexSuppressedByHematoxylin, Eigen::Index);
  itkSetMacro(ColorIndexSuppressedByHematoxylin, Eigen::Index);

  /** If the pixel type is RGB or RGBA then
   * ColorIndexSuppressedByEosin defaults to 1, indicating green.
   * Otherwise, set ColorIndexSuppressedByEosin to the index in the
   * array of colors that indicates the color most suppressed by
   * eosin. */
  itkGetMacro(ColorIndexSuppressedByEosin, Eigen::Index);
  itkSetMacro(ColorIndexSuppressedByEosin, Eigen::Index);

  // This algorithm is defined for H&E (Hematoxylin (blue) and
  // Eosin (pink)), which is a total of 2 stains.  However, this
  // approach could in theory work in other circumstances.  In that
  // case it might be better to have NumberOfStains be a template
  // parameter or a setable class member.

  /** Hematoxylin and eosin; there are two stains supported. */
  static constexpr SizeValueType NumberOfStains{ 2 };
  /** For Virtanen's non-negative matrix factorization algorithm. */
  static constexpr SizeValueType maxNumberOfIterations{ 0 };
  /** Select a subset of the pixels if the image has more than this */
  static constexpr SizeValueType maxNumberOfRows{ 100000 };
  /** Colors at least this fraction as distance as first pass distinguisher are good substitutes for it. */
  static constexpr CalcElementType SecondPassDistinguishersThreshold{ 0.90 };
  /** Colors that are at least this percentile in brightness are considered bright. */
  static constexpr CalcElementType BrightPercentileLevel{ 0.80 };
  /** Colors that are at least fraction of white's brightness are considered bright. */
  static constexpr CalcElementType BrightPercentageLevel{ 0.50 };
  /** Intensity normalization of an image is based upon dark pixels of this percentil brightness. */
  static constexpr CalcElementType VeryDarkPercentileLevel{ 0.01 };
  /** A small matrix.array_inf_norm() value for checking convergence of a matrix */
  static constexpr CalcElementType epsilon0{ 1e-2 };
  /** A very small squared magnitude for a vector, to prevent division by zero. */
  static constexpr CalcElementType epsilon2{ 1e-12 };
  /** The Lasso optimization penalty. */
  static constexpr CalcElementType lambda{ 0.00 };

protected:
  // We have special cases for different pixel types, including: (1)
  // we refuse to process fewer than 3 colors (at compile time where
  // possible) and (2) RBGA pixels have 4 dimensions but only 3
  // colors.  We accomplish the special cases via the PixelHelper
  // class, similar in concept to the PixelTraits class, but with
  // different capabilities.
  //
  // The default: for the case that the number of colors is not
  // determined at compile time, such as with a VectorImage:
  template <typename TPixelType, typename = void>
  struct PixelHelper
  {
    using PixelType = TPixelType;
    using ValueType = typename PixelType::ValueType;
    static constexpr SizeValueType         NumberOfDimensions = -1;
    static constexpr SizeValueType         NumberOfColors = -1;
    static constexpr typename Eigen::Index ColorIndexSuppressedByHematoxylin = -1;
    static constexpr typename Eigen::Index ColorIndexSuppressedByEosin = -1;
    static PixelType
    pixelInstance(unsigned numberOfDimensions)
    {
      return PixelType{ numberOfDimensions };
    }
  };
  // For the case that the number of colors is implicitly set to 1 at
  // compile time.  We will refuse to compile this case via a
  // static_assert in the constructor for
  // StructurePreservingColorNormalizationFilter.
  template <typename TPixelType>
  struct PixelHelper<TPixelType, typename std::enable_if<std::is_arithmetic<TPixelType>::value>::type>
  {
    using PixelType = TPixelType;
    using ValueType = PixelType;
    static constexpr SizeValueType         NumberOfDimensions = 1;
    static constexpr SizeValueType         NumberOfColors = 1;
    static constexpr typename Eigen::Index ColorIndexSuppressedByHematoxylin = -1;
    static constexpr typename Eigen::Index ColorIndexSuppressedByEosin = -1;
    static PixelType
    pixelInstance(unsigned numberOfDimensions)
    {
      return PixelType{};
    }
  };
  // For the case that the pixel type is RGBPixel:
  template <typename TScalar>
  struct PixelHelper<RGBPixel<TScalar>, void>
  {
    using PixelType = RGBPixel<TScalar>;
    using ValueType = typename PixelType::ValueType;
    static constexpr SizeValueType         NumberOfDimensions = 3;
    static constexpr SizeValueType         NumberOfColors = 3;
    static constexpr typename Eigen::Index ColorIndexSuppressedByHematoxylin = 0;
    static constexpr typename Eigen::Index ColorIndexSuppressedByEosin = 1;
    static PixelType
    pixelInstance(unsigned numberOfDimensions)
    {
      return PixelType{};
    }
  };
  // For the case that the pixel type is RGBAPixel:
  template <typename TScalar>
  struct PixelHelper<RGBAPixel<TScalar>, void>
  {
    using PixelType = RGBAPixel<TScalar>;
    using ValueType = typename PixelType::ValueType;
    static constexpr SizeValueType         NumberOfDimensions = 4;
    static constexpr SizeValueType         NumberOfColors = 3;
    static constexpr typename Eigen::Index ColorIndexSuppressedByHematoxylin = 0;
    static constexpr typename Eigen::Index ColorIndexSuppressedByEosin = 1;
    static PixelType
    pixelInstance(unsigned numberOfDimensions)
    {
      return PixelType{};
    }
  };
  // For the cases that the pixel type is Vector or
  // CovariantVector.  If NVectorDimension is not at least 3 we
  // will refuse to compile this case via a static_assert in the
  // constructor of StructurePreservingColorNormalizationFilter.
  template <typename TScalar, unsigned int NVectorDimension>
  struct PixelHelper<Vector<TScalar, NVectorDimension>, void>
  {
    using PixelType = Vector<TScalar, NVectorDimension>;
    using ValueType = typename PixelType::ValueType;
    static constexpr SizeValueType         NumberOfDimensions = NVectorDimension;
    static constexpr SizeValueType         NumberOfColors = NVectorDimension;
    static constexpr typename Eigen::Index ColorIndexSuppressedByHematoxylin = -1;
    static constexpr typename Eigen::Index ColorIndexSuppressedByEosin = -1;
    static PixelType
    pixelInstance(unsigned numberOfDimensions)
    {
      return PixelType{};
    }
  };
  template <typename TScalar, unsigned int NVectorDimension>
  struct PixelHelper<CovariantVector<TScalar, NVectorDimension>, void>
  {
    using PixelType = CovariantVector<TScalar, NVectorDimension>;
    using ValueType = typename PixelType::ValueType;
    static constexpr SizeValueType         NumberOfDimensions = NVectorDimension;
    static constexpr SizeValueType         NumberOfColors = NVectorDimension;
    static constexpr typename Eigen::Index ColorIndexSuppressedByHematoxylin = -1;
    static constexpr typename Eigen::Index ColorIndexSuppressedByEosin = -1;
    static PixelType
    pixelInstance(unsigned numberOfDimensions)
    {
      return PixelType{};
    }
  };

public:
  // This public type depends upon the protected PixelHelper, so we
  // must go back to "public:" here.
  /** Additional specific class typedefs */
  using PixelValueType = typename PixelHelper<PixelType>::ValueType;

protected:
  StructurePreservingColorNormalizationFilter();
  ~StructurePreservingColorNormalizationFilter() override = default;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  VerifyInputInformation() const override;

  void
  GenerateInputRequestedRegion() override;

  void
  BeforeThreadedGenerateData() override;

  void
  DynamicThreadedGenerateData(const RegionType & outputRegion) override;

  int
  ImageToNMF(RegionConstIterator & iter, CalcMatrixType & matrixH, CalcRowVectorType & unstainedPixel) const;

  void
  ImageToMatrix(RegionConstIterator & iter,
                SizeValueType         numberOfPixels,
                CalcMatrixType &      matrixBrightV,
                CalcMatrixType &      matrixDarkV) const;

  static void
  MatrixToDistinguishers(const CalcMatrixType & matrixV, CalcMatrixType & distinguishers);

  static void
  MatrixToMatrixExtremes(const CalcMatrixType & matrixV, CalcMatrixType & matrixBrightV, CalcMatrixType & matrixDarkV);

  static void
  FirstPassDistinguishers(const CalcMatrixType &                normVStart,
                          std::array<int, NumberOfStains + 1> & firstPassDistinguisherIndices,
                          SizeValueType &                       numberOfDistinguishers);

  static void
  SecondPassDistinguishers(const CalcMatrixType &                      normVStart,
                           const std::array<int, NumberOfStains + 1> & firstPassDistinguisherIndices,
                           const SizeValueType                         numberOfDistinguishers,
                           CalcMatrixType &                            secondPassDistinguisherColors);

  static int
  MatrixToOneDistinguisher(const CalcMatrixType & normV);

  static CalcMatrixType
  RecenterMatrix(const CalcMatrixType & normV, const SizeValueType row);

  static CalcMatrixType
  ProjectMatrix(const CalcMatrixType & normV, const SizeValueType row);

  int
  DistinguishersToNMFSeeds(const CalcMatrixType & distinguishers,
                           CalcRowVectorType &    unstainedPixel,
                           CalcMatrixType &       matrixH) const;

  void
  DistinguishersToColors(const CalcMatrixType & distinguishers,
                         SizeValueType &        unstainedIndex,
                         SizeValueType &        hematoxylinIndex,
                         SizeValueType &        eosinIndex) const;

  void
  NormalizeMatrixH(const CalcMatrixType &    matrixDarkV,
                   const CalcRowVectorType & unstainedPixel,
                   CalcMatrixType &          matrixH) const;

  static void
  VirtanenNMFEuclidean(const CalcMatrixType & matrixV, CalcMatrixType & matrixW, CalcMatrixType & matrixH);

  static void
  VirtanenNMFKLDivergence(const CalcMatrixType & matrixV, CalcMatrixType & matrixW, CalcMatrixType & matrixH);

  void
  NMFsToImage(const CalcMatrixType &    inputH,
              const CalcRowVectorType & inputUnstained,
              const CalcMatrixType &    referH,
              const CalcRowVectorType & referUnstained,
              RegionIterator &          outIt) const;

  // Our installation of Eigen3 does not have iterators.  (They
  // arrive with Eigen 3.4.)  We define begin, cbegin, end, and cend
  // functions here.  A compiler sometimes gets segmentation fault if
  // we use the more restrictive Eigen::Matrix< ... > declarations, so
  // we have the more lax TMatrix declarations available as
  // alternates.
#define STRUCTUREPRESERVINGCOLORNORMALIZATIONFILTER_STRICT_EIGEN3_ITERATORS 0
#if STRUCTUREPRESERVINGCOLORNORMALIZATIONFILTER_STRICT_EIGEN3_ITERATORS
  template <typename TScalar, int VRows, int VCols, int VOptions, int VMaxRows, int VMaxCols>
  static TScalar *
  begin(Eigen::Matrix<TScalar, VRows, VCols, VOptions, VMaxRows, VMaxCols> & matrix);

  template <typename TScalar, int VRows, int VCols, int VOptions, int VMaxRows, int VMaxCols>
  static const TScalar *
  cbegin(const Eigen::Matrix<TScalar, VRows, VCols, VOptions, VMaxRows, VMaxCols> & matrix);

  template <typename TScalar, int VRows, int VCols, int VOptions, int VMaxRows, int VMaxCols>
  static TScalar *
  end(Eigen::Matrix<TScalar, VRows, VCols, VOptions, VMaxRows, VMaxCols> & matrix);

  template <typename TScalar, int VRows, int VCols, int VOptions, int VMaxRows, int VMaxCols>
  static const TScalar *
  cend(const Eigen::Matrix<TScalar, VRows, VCols, VOptions, VMaxRows, VMaxCols> & matrix);

#else
  template <typename TMatrix>
  static CalcElementType *
  begin(TMatrix & matrix);

  template <typename TMatrix>
  static const CalcElementType *
  cbegin(const TMatrix & matrix);

  template <typename TMatrix>
  static CalcElementType *
  end(TMatrix & matrix);

  template <typename TMatrix>
  static const CalcElementType *
  cend(const TMatrix & matrix);
#endif

  // These members are for the purpose of caching results for use the
  // next time the pipeline is run.
  ModifiedTimeType  m_ParametersMTime;
  const ImageType * m_Input;
  CalcMatrixType    m_InputH;
  CalcRowVectorType m_InputUnstainedPixel;
  const ImageType * m_Reference;
  TimeStamp         m_ReferenceTimeStamp;
  CalcMatrixType    m_ReferenceH;
  CalcRowVectorType m_ReferenceUnstainedPixel;

  Eigen::Index m_NumberOfDimensions;
  Eigen::Index m_NumberOfColors;
  Eigen::Index m_ColorIndexSuppressedByHematoxylin;
  Eigen::Index m_ColorIndexSuppressedByEosin;

private:
#ifdef ITK_USE_CONCEPT_CHECKING
  // Add concept checking such as
  // itkConceptMacro(FloatingPointPixel, (Concept::IsFloatingPoint<typename ImageType::PixelType>));
#endif
};

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkStructurePreservingColorNormalizationFilter.hxx"
#endif

#endif // itkStructurePreservingColorNormalizationFilter
