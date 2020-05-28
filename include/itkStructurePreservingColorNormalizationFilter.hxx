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

#ifndef itkStructurePreservingColorNormalizationFilter_hxx
#define itkStructurePreservingColorNormalizationFilter_hxx

#include "itkStructurePreservingColorNormalizationFilter.h"
#include "itkRGBPixel.h"
#include <numeric>

namespace itk
{

template< typename TInputImage, typename TOutputImage >
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::StructurePreservingColorNormalizationFilter()
  : m_ColorIndexSuppressedByHematoxylin( 0 ), m_ColorIndexSuppressedByEosin( 1 )
{}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "ColorIndexSuppressedByHematoxylin: " << m_ColorIndexSuppressedByHematoxylin << std::endl
     << indent << "ColorIndexSuppressedByEosin: " << m_ColorIndexSuppressedByEosin << std::endl;
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::GenerateInputRequestedRegion()
{
  // Call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // Get pointers to the input ( to be normalized ) and reference
  // images
  InputImageType *inputPtr = const_cast< InputImageType * >( this->GetInput( 0 ) );
  InputImageType *referPtr = const_cast< InputImageType * >( this->GetInput( 1 ) );

  if( inputPtr != nullptr )
    {
    inputPtr->SetRequestedRegionToLargestPossibleRegion();
    }

  if( referPtr != nullptr )
    {
    referPtr->SetRequestedRegionToLargestPossibleRegion();
    }
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::BeforeThreadedGenerateData()
{
  itkAssertOrThrowMacro( InputImageLength >= 3,
    "itkStructurePreservingColorNormalizationFilter input images need length ( #colors ) >= 3." );
  itkAssertOrThrowMacro( OutputImageLength == InputImageLength,
    "StructurePreservingColorNormalizationFilter output image needs length ( #colors ) exactly the same as the input images." );

  // Call the superclass' implementation of this method
  Superclass::BeforeThreadedGenerateData();

  // Find input and refer and make iterators for them.
  const InputImageType * const inputPtr = this->GetInput( 0 ); // image to be normalized
  const InputImageType * const referPtr = this->GetInput( 1 ); // reference image
  // For each of the two images, make sure that it was supplied, or
  // that we have it cached already.
  itkAssertOrThrowMacro( inputPtr != nullptr || m_inputPtr != nullptr, "An image to be normalized needs to be supplied as input image #0" );
  itkAssertOrThrowMacro( referPtr != nullptr || m_referPtr != nullptr, "An reference image needs to be supplied as input image #1" );

  // For each input, if there is a supplied image and it is different
  // from what we have cached then compute stuff and cache the
  // results.  These two calls to ImageToNMF could be done
  // simultaneously.
  if( inputPtr != nullptr && ( inputPtr != m_inputPtr || inputPtr->GetTimeStamp() != m_inputTimeStamp ) )
    {
    InputRegionConstIterator inputIter {inputPtr, inputPtr->GetRequestedRegion()};
    if( this->ImageToNMF( inputIter, m_inputH, m_inputUnstainedPixel ) == 0 )
      {
      m_inputPtr = inputPtr;
      m_inputTimeStamp = inputPtr->GetTimeStamp();
      { std::ostringstream mesg; mesg << "m_inputH = " << m_inputH << std::endl; std::cout << mesg.str() << std::flush; }
      }
    else
      {
      // we failed
      m_inputPtr = nullptr;
      itkAssertOrThrowMacro( m_inputPtr != nullptr, "The image to be normalized could not be processed" )
      }
    }

  if( referPtr != nullptr && ( referPtr != m_referPtr || referPtr->GetTimeStamp() != m_referTimeStamp ) )
    {
    InputRegionConstIterator referIter {referPtr, referPtr->GetRequestedRegion()};
    if( this->ImageToNMF( referIter, m_referH, m_referUnstainedPixel ) == 0 )
      {
      m_referPtr = referPtr;
      m_referTimeStamp = referPtr->GetTimeStamp();
      { std::ostringstream mesg; mesg << "m_referH = " << m_referH << std::endl; std::cout << mesg.str() << std::flush; }
      }
    else
      {
      // we failed
      m_referPtr = nullptr;
      itkAssertOrThrowMacro( m_referPtr != nullptr, "The reference image could not be processed" )
      }
    }

  if( ( m_inputH * m_referH.transpose() ).determinant() < CalcElementType( 0 ) )
    {
    // Somehow the hematoxylin and eosin rows got swapped in one of
    // the input image or reference image.  Flip them in referH to get
    // them in synch.
    static_assert( NumberOfStains == 2, "StructurePreservingColorNormalizationFilter current implementation assumes exactly two stains" );
    const CalcMatrixType referHOriginal {m_referH};
    m_referH.row( 0 ) = referHOriginal.row( 1 );
    m_referH.row( 1 ) = referHOriginal.row( 0 );
    }

  // Correct for any scaling difference between m_referH and m_inputH.
  const CalcColVectorType lastOnes {CalcColVectorType::Constant( m_referH.cols(), 1, 1.0 )};
  m_referH = ( ( ( ( m_inputH.array() * m_inputH.array() ).matrix() * lastOnes ).array() + epsilon2 )
             / ( ( ( m_referH.array() * m_referH.array() ).matrix() * lastOnes ).array() + epsilon2 ) )
    .matrix().unaryExpr( CalcUnaryFunctionPointer( std::sqrt ) ).asDiagonal() * m_referH;
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::DynamicThreadedGenerateData( const OutputRegionType & outputRegion )
{
  OutputImageType * const outputPtr = this->GetOutput();
  itkAssertOrThrowMacro( outputPtr != nullptr, "An output image needs to be supplied" )
  OutputRegionIterator outputIter {outputPtr, outputRegion};

  this->NMFsToImage( m_inputH, m_inputUnstainedPixel, m_referH, m_referUnstainedPixel, outputIter );
}


template< typename TInputImage, typename TOutputImage >
int
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::ImageToNMF( InputRegionConstIterator &inIter, CalcMatrixType &matrixH, InputPixelType &unstainedPixel ) const
{
  // To maintain locality of memory references, we are using
  // numberOfPixels as the number of rows rather than as the number of
  // columns.  With V=WH, as is standard in non-negative matrix
  // factorization, our matrices switch names and are transposed with
  // respect to the Vahadane article.  In particular, our W is a tall
  // matrix and our H is a fairly compact matrix.

  const InputSizeType size = inIter.GetRegion().GetSize();
  const InputSizeValueType numberOfPixels = std::accumulate( size.begin(), size.end(), 1, std::multiplies< InputSizeValueType >() );

  // Find distinguishers.  These are essentially the rows of matrixH.
  CalcMatrixType distinguishers;
  CalcMatrixType matrixV {numberOfPixels, InputImageLength};

  this->ImageToMatrix( inIter, matrixV );
  this->MatrixToDistinguishers( matrixV, distinguishers );

  // Use the distinguishers as seeds to the non-negative matrix
  // factorization.
  if( this->DistinguishersToNMFSeeds( distinguishers, unstainedPixel, matrixH ) != 0 )
    {
    return 1;                   // we failed.
    }

  return 0;
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::ImageToMatrix( InputRegionConstIterator &inIter, CalcMatrixType &matrixV ) const
{
  InputSizeValueType pixelIndex {0};
  for( inIter.GoToBegin(); !inIter.IsAtEnd(); ++inIter, ++pixelIndex )
    {
    InputPixelType pixelValue = inIter.Get();
    for( Eigen::Index color = 0; color < InputImageLength; ++color )
      {
      matrixV( pixelIndex, color ) = InputPixelHelper::value( pixelValue, color );
      }
    }

  // We do not want trouble with a value near zero ( when we take its
  // logarithm ) so we add a little to each value now.
  const CalcElementType nearZero {matrixV.lpNorm< Eigen::Infinity >() * epsilon1};
  matrixV = ( matrixV.array() + nearZero ).matrix();

  // Keep only pixels that are bright enough.
  this->MatrixToBrightPartOfMatrix( matrixV );
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::MatrixToBrightPartOfMatrix( CalcMatrixType &matrixV ) const
{
  // A useful vector that has a 1 for each column of matrixV.
  const CalcColVectorType lastOnes {CalcColVectorType::Constant( matrixV.cols(), 1, 1.0 )};

  // We want only the brightest pixels.  Find the 80th percentile threshold.
  const CalcColVectorType brightnessOriginal {matrixV * lastOnes};
  CalcColVectorType brightnessOrdered {brightnessOriginal};
  const CalcElementType percentileLevel {0.80};
  InputSizeValueType const quantilePosition {static_cast< InputSizeValueType >( ( brightnessOrdered.size() - 1 ) * percentileLevel )};
  std::nth_element( Self::begin( brightnessOrdered ), Self::begin( brightnessOrdered ) + quantilePosition, Self::end( brightnessOrdered ) );
  const CalcElementType percentileThreshold {brightnessOrdered( quantilePosition )};
  // Find 70% of maximum brightness
  const CalcElementType percentageLevel {0.70};
  const CalcElementType percentageThreshold {percentageLevel * *std::max_element( Self::cbegin( brightnessOriginal ), Self::cend( brightnessOriginal ) )};

  // We will keep those pixels that pass at least one of the above
  // thresholds.
  const CalcElementType brightnessThreshold {std::min( percentileThreshold, percentageThreshold )};
  InputSizeValueType numberOfRowsToKeep {0};
  for( Eigen::Index i = 0 ; i < matrixV.rows(); ++i )
    {
    if( brightnessOriginal( i ) >= brightnessThreshold )
      {
      ++numberOfRowsToKeep;
      }
    }
  CalcMatrixType brightV {numberOfRowsToKeep, matrixV.cols()};
  numberOfRowsToKeep = 0;
  for( Eigen::Index i = 0 ; i < matrixV.rows(); ++i )
    {
    if( brightnessOriginal( i ) >= brightnessThreshold )
      {
      brightV.row( numberOfRowsToKeep++ ) = matrixV.row( i );
      }
    }
  matrixV = brightV;
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::MatrixToDistinguishers( const CalcMatrixType &matrixV, CalcMatrixType &distinguishers ) const
{
  const CalcMatrixType normVStart {matrixV};

  // We will store the row ( pixel ) index of each distinguishing
  // pixel in firstPassDistinguisherIndices.
  std::array< int, NumberOfStains+1 > firstPassDistinguisherIndices {-1};
  InputSizeValueType numberOfDistinguishers {0};
  this->FirstPassDistinguishers( normVStart, firstPassDistinguisherIndices, numberOfDistinguishers );

  // Each row of secondPassDistinguisherColors is the vector of color
  // values for a distinguisher.
  CalcMatrixType secondPassDistinguisherColors {numberOfDistinguishers, matrixV.cols()};
  this->SecondPassDistinguishers( normVStart, firstPassDistinguisherIndices, numberOfDistinguishers, secondPassDistinguisherColors );

  distinguishers = secondPassDistinguisherColors;
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::FirstPassDistinguishers( const CalcMatrixType &normVStart, std::array< int, NumberOfStains+1 > &firstPassDistinguisherIndices, InputSizeValueType &numberOfDistinguishers ) const
{
  // A useful vector that has a 1 for each column of normVStart.
  const CalcColVectorType lastOnes {CalcColVectorType::Constant( normVStart.cols(), 1, 1.0 )};
  // A useful vector that has a 1 for each row of normVStart.
  const CalcColVectorType firstOnes {CalcColVectorType::Constant( normVStart.rows(), 1, 1.0 )};

  CalcMatrixType normV {normVStart};
  numberOfDistinguishers = 0;
  bool needToRecenterMatrix = true;
  while( numberOfDistinguishers <= NumberOfStains )
    {
    // Find the next distinguishing row ( pixel )
    firstPassDistinguisherIndices[numberOfDistinguishers] = this->MatrixToOneDistinguisher( normV, lastOnes );
    // If we found a distinguisher and we have not yet found
    // NumberOfStains+1 of them, then look for the next distinguisher.
    if( firstPassDistinguisherIndices[numberOfDistinguishers] >= 0 )
      {
      // We just found a distinguisher
      ++numberOfDistinguishers;
      if( numberOfDistinguishers <= NumberOfStains )
        {
        // Prepare to look for the next distinguisher
        if( needToRecenterMatrix )
          {
          normV = this->RecenterMatrix( normV, firstOnes, firstPassDistinguisherIndices[numberOfDistinguishers - 1] );
          needToRecenterMatrix = false;
          }
        else
          {
          normV = this->ProjectMatrix( normV, firstPassDistinguisherIndices[numberOfDistinguishers - 1] );
          }
        }
      }
    else
      {
      // We did not find another distinguisher.  There are no more.
      break;
      }
    }
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::SecondPassDistinguishers( const CalcMatrixType &normVStart, const std::array< int, NumberOfStains+1 > &firstPassDistinguisherIndices, const InputSizeValueType numberOfDistinguishers,
  CalcMatrixType &secondPassDistinguisherColors ) const
{
  // A useful vector that has a 1 for each column of normVStart.
  const CalcColVectorType lastOnes {CalcColVectorType::Constant( normVStart.cols(), 1, 1.0 )};
  // A useful vector that has a 1 for each row of normVStart.
  const CalcColVectorType firstOnes {CalcColVectorType::Constant( normVStart.rows(), 1, 1.0 )};

  for( int distinguisher {0}; distinguisher < numberOfDistinguishers; ++distinguisher )
    {
    CalcMatrixType normV {normVStart};
    bool needToRecenterMatrix = true;
    for( int otherDistinguisher {0}; otherDistinguisher < numberOfDistinguishers; ++otherDistinguisher )
      {
      // skip if self
      if( otherDistinguisher != distinguisher )
        {
        if( needToRecenterMatrix )
          {
          normV = this->RecenterMatrix( normV, firstOnes, firstPassDistinguisherIndices[otherDistinguisher] );
          needToRecenterMatrix = false;
          }
        else
          {
          normV = this->ProjectMatrix( normV, firstPassDistinguisherIndices[otherDistinguisher] );
          }
        }
      }
    // We have sent all distinguishers except self to the origin.
    // Whatever is far from the origin in the same direction as self
    // is a good replacement for self.  We will take an average among
    // those that are at least 80% as far as the best.  ( Note that
    // self could still be best, but not always. )
    const CalcColVectorType dotProducts {normV * normV.row( firstPassDistinguisherIndices[distinguisher] ).transpose()};
    const CalcElementType threshold {*std::max_element( Self::cbegin( dotProducts ), Self::cend( dotProducts ) ) * 999 / 1000};
    CalcRowVectorType cumulative {CalcRowVectorType::Constant( 1, normVStart.cols(), 0.0 )};
    InputSizeValueType numberOfContributions {0};
    for( InputSizeValueType row {0}; row < dotProducts.size(); ++row )
      {
      if( dotProducts( row ) >= threshold )
        {
        cumulative += normVStart.row( row );
        ++numberOfContributions;
        }
      }
    secondPassDistinguisherColors.row( distinguisher ) = cumulative / numberOfContributions;
    }
}


// static method
template< typename TInputImage, typename TOutputImage >
int
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::MatrixToOneDistinguisher( const CalcMatrixType &normV, const CalcColVectorType &lastOnes )
{
  const CalcColVectorType lengths2 = ( normV.array() * normV.array() ).matrix() * lastOnes;
  const CalcElementType * const result {std::max_element( Self::cbegin( lengths2 ), Self::cend( lengths2 ) )};
  if( *result > epsilon2 )
    {
    return std::distance( Self::cbegin( lengths2 ), result );
    }
  else
    {
    return -1;                // Nothing left to find
    }
}


// static method
template< typename TInputImage, typename TOutputImage >
typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcMatrixType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::RecenterMatrix( const CalcMatrixType &normV, const CalcColVectorType &firstOnes, const InputSizeValueType row )
{
  return normV - firstOnes * normV.row( row );
}


// static method
template< typename TInputImage, typename TOutputImage >
typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcMatrixType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::ProjectMatrix( const CalcMatrixType &normV, const InputSizeValueType row )
{
  const CalcRowVectorType rowValue {normV.row( row )};
  return normV - ( normV * rowValue.transpose() ) * ( rowValue / rowValue.squaredNorm() );
}


template< typename TInputImage, typename TOutputImage >
int
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::DistinguishersToNMFSeeds( const CalcMatrixType &distinguishers, InputPixelType &unstainedPixel, CalcMatrixType &matrixH ) const
{
  matrixH = CalcMatrixType {NumberOfStains, InputImageLength};

  const CalcRowVectorType midOnes {CalcRowVectorType::Constant( 1, matrixH.rows(), 1.0 )};
  const CalcColVectorType lastOnes {CalcColVectorType::Constant( matrixH.cols(), 1, 1.0 )};

  InputSizeValueType unstainedIndex;
  InputSizeValueType hematoxylinIndex;
  InputSizeValueType eosinIndex;
  this->DistinguishersToColors( distinguishers, unstainedIndex, hematoxylinIndex, eosinIndex );

  // If the indices unstainedIndex, hematoxylinIndex, and eosinIndex
  // are distinct then we choose a smart starting place for the
  // generic NMF algorithm.  Otherwise, we go with a guess that is
  // reasonable.
  if( unstainedIndex == hematoxylinIndex || unstainedIndex == eosinIndex || hematoxylinIndex == eosinIndex )
    {
    return 1;                   // we failed
    }

  const CalcRowVectorType unstainedCalcPixel {distinguishers.row( unstainedIndex )};
  for( Eigen::Index color = 0; color < InputImageLength; ++ color )
    {
    InputPixelHelper::value( unstainedPixel, color ) = unstainedCalcPixel( color ); // return value
    }
  const CalcRowVectorType logUnstained {unstainedCalcPixel.unaryExpr( CalcUnaryFunctionPointer( std::log ) )};
  const CalcRowVectorType logHematoxylin {logUnstained - distinguishers.row( hematoxylinIndex ).unaryExpr( CalcUnaryFunctionPointer( std::log ) )};
  const CalcRowVectorType logEosin {logUnstained - distinguishers.row( eosinIndex ).unaryExpr( CalcUnaryFunctionPointer( std::log ) )};
  // Set rows of matrixH to reflect hematoxylin and eosin.
  matrixH.row( 0 ) = logHematoxylin;
  matrixH.row( 1 ) = logEosin;

  // Make sure that each row of matrixH has unit magnitude, and each
  // element of matrixH is sufficiently non-negative.
  const auto clip = [] ( const CalcElementType &x )
    {
    return std::max( CalcElementType( 0.0 ), x );
    };
  matrixH = ( ( ( matrixH.array() * matrixH.array() ).matrix() * lastOnes ).unaryExpr( CalcUnaryFunctionPointer( std::sqrt ) ) ).asDiagonal().inverse() * matrixH;
  matrixH = matrixH.unaryExpr( clip );

  return 0;
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::DistinguishersToColors( CalcMatrixType const &distinguishers, InputSizeValueType &unstainedIndex, InputSizeValueType &hematoxylinIndex, InputSizeValueType &eosinIndex ) const
{
  // Figure out which, distinguishers are unstained ( highest
  // brightness ), hematoxylin ( suppresses red ), and eosin (
  // suppresses green ).
  const CalcColVectorType lastOnes {CalcColVectorType::Constant( distinguishers.cols(), 1, 1.0 )};
  const CalcColVectorType lengths2 {( distinguishers.array() * distinguishers.array() ).matrix() * lastOnes};
  const CalcElementType * const unstainedIterator {std::max_element( Self::cbegin( lengths2 ), Self::cend( lengths2 ) )};
  unstainedIndex =  std::distance( Self::cbegin( lengths2 ), unstainedIterator );
  // For typename RGBPixel, red is suppressed by hematoxylin and green
  // is suppressed by eosin.  More generally, the index of the color
  // suppressed by hematoxylin is indicated by
  // m_ColorIndexSuppressedByHematoxylin, and the index of the color
  // suppressed by eosin is indicated by
  // m_ColorIndexSuppressedByEosin.
  const CalcColVectorType redValues {distinguishers.col( m_ColorIndexSuppressedByHematoxylin )};
  const CalcElementType * const hematoxylinputIterator {std::min_element( Self::cbegin( redValues ), Self::cend( redValues ) )};
  hematoxylinIndex = std::distance( Self::cbegin( redValues ), hematoxylinputIterator );
  const CalcColVectorType greenValues {distinguishers.col( m_ColorIndexSuppressedByEosin )};
  const CalcElementType * const eosinputIterator {std::min_element( Self::cbegin( greenValues ), Self::cend( greenValues ) )};
  eosinIndex = std::distance( Self::cbegin( greenValues ), eosinputIterator );
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::VirtanenEuclidean( const CalcMatrixType &matrixV, CalcMatrixType &matrixW, const CalcMatrixType &matrixH ) const
{
  const auto clip = [] ( const CalcElementType &x )
    {
    return std::max( CalcElementType( 0.0 ), x );
    };
  matrixW = ( ( ( matrixV * matrixH.transpose() ).array() - lambda ).unaryExpr( clip ) + epsilon2 ).matrix() * ( matrixH * matrixH.transpose() ).inverse();

  const CalcColVectorType lastOnes {CalcColVectorType::Constant( matrixV.cols(), 1, 1.0 )};
  // Apply Virtanen's algorithm to iteratively improve matrixW and
  // matrixH.  Note that parentheses optimize the order of matrix
  // chain multiplications and affect the speed of this method.
  CalcMatrixType previousMatrixW {matrixW};
  InputSizeValueType loopIter {0};
  for( ; loopIter < maxNumberOfIterations; ++loopIter )
    {
    // Lasso term "lambda" insertion is possibly in a novel way.
    matrixW = ( matrixW.array()
              * ( ( ( ( matrixV * matrixH.transpose() ).array() - lambda ).unaryExpr( clip ) + epsilon2 )
                / ( ( matrixW * ( matrixH * matrixH.transpose() ) ).array() + epsilon2 ) ) ).matrix();
    // matrixH = ( matrixH.array()
    //           * ( ( ( matrixW.transpose() * matrixV ).array() + epsilon2 )
    //             / ( ( ( matrixW * matrixW.transpose() ) * matrixH ).array() + epsilon2 ) ) ).matrix();
    // In lieu of rigorous Lagrange multipliers, renormalize rows of
    // matrixH to have unit magnitude.
    // matrixH = ( ( ( matrixH.array() * matrixH.array() ).matrix() * lastOnes ).unaryExpr( CalcUnaryFunctionPointer( std::sqrt ) ) ).asDiagonal().inverse() * matrixH;
    if( ( loopIter & 15 ) == 15 )
      {
      if( ( matrixW - previousMatrixW ).lpNorm< Eigen::Infinity >() < epsilon0 )
        {
        break;
        }
      previousMatrixW = matrixW;
      }
    }

  // Round off values in the response, so that numbers are quite small
  // are set to zero.
  const CalcElementType maxW = matrixW.lpNorm< Eigen::Infinity >() * 15;
  matrixW = ( ( matrixW.array() + maxW ) - maxW ).matrix();
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::VirtanenKLDivergence( const CalcMatrixType &matrixV, CalcMatrixType &matrixW, const CalcMatrixType &matrixH ) const
{
  // If this method is going to get used, we may need to incorporate
  // the Lasso penalty lambda for matrixW and incorporate the Lagrange
  // multipliers to make each row of matrixH have magnitude 1.0.

  // Apply Virtanen's algorithm to iteratively improve matrixW and
  // matrixH.
  const CalcRowVectorType firstOnes {CalcRowVectorType::Constant( 1, matrixV.rows(), 1.0 )};
  const CalcColVectorType lastOnes {CalcColVectorType::Constant( matrixV.cols(), 1, 1.0 )};
  CalcMatrixType previousMatrixW {matrixW};
  for( InputSizeValueType loopIter {0}; loopIter < maxNumberOfIterations; ++loopIter )
    {
    matrixW = ( matrixW.array()
              * ( ( ( ( matrixV.array() + epsilon2 ) / ( ( matrixW * matrixH ).array() + epsilon2 ) + epsilon2 ).matrix() * matrixH.transpose() ).array()
                / ( ( ( matrixH * lastOnes ) * firstOnes ).transpose().array() + epsilon2 ) ) ).matrix();
    // matrixH = ( matrixW.array()
    //           * ( ( ( matrixW.transpose() * ( ( matrixV.array() + epsilon2 ) / ( ( matrixW * matrixH ).array() + epsilon2 ) ).matrix() ).array() + epsilon2 )
    //             / ( ( lastOnes * ( firstOnes * matrixW ) ).transpose().array() + epsilon2 ) ) ).matrix();
    if( ( loopIter & 15 ) == 15 )
      {
      if( ( matrixW - previousMatrixW ).lpNorm< Eigen::Infinity >() < epsilon0 )
        break;
      previousMatrixW = matrixW;
      }
    }

  // Round off values in the response, so that numbers are quite small
  // are set to zero.
  const CalcElementType maxW = matrixW.lpNorm< Eigen::Infinity >() * 15;
  matrixW = ( ( matrixW.array() + maxW ) - maxW ).matrix();
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::NMFsToImage( const CalcMatrixType &inputH, const InputPixelType &inputUnstained, const CalcMatrixType &referH, const InputPixelType &referUnstained, OutputRegionIterator &outputIter ) const
{
  // Read in corresponding part of the input region.
  const OutputSizeType size = outputIter.GetRegion().GetSize();
  const OutputSizeValueType numberOfPixels = std::accumulate( size.begin(), size.end(), 1, std::multiplies< OutputSizeValueType >() );
  CalcMatrixType matrixV {numberOfPixels, OutputImageLength};
  InputRegionConstIterator inputIter {m_inputPtr, m_inputPtr->GetRequestedRegion()};
  outputIter.GoToBegin();
  inputIter.GoToBegin();
  for( OutputSizeValueType pixelIndex {0}; !outputIter.IsAtEnd(); ++outputIter, ++inputIter, ++pixelIndex )
    {
    // Find input index that matches this output index.
    while ( inputIter.GetIndex() != outputIter.GetIndex() )
      {
      ++inputIter;
      }
    InputPixelType pixelValue = inputIter.Get();
    for( Eigen::Index color = 0; color < InputImageLength; ++color )
      {
      matrixV( pixelIndex, color ) = InputPixelHelper::value( pixelValue, color );
      }
    }

  // Convert matrixV using the inputUnstained pixel and a call to
  // logarithm.
  CalcRowVectorType logInputUnstained {1, InputImageLength};
  CalcRowVectorType logReferUnstained {1, InputImageLength};
  for( Eigen::Index color = 0; color < InputImageLength; ++color )
    {
    logInputUnstained[color] = std::log( CalcElementType( InputPixelHelper::value( inputUnstained, color ) ) );
    logReferUnstained[color] = std::log( CalcElementType( InputPixelHelper::value( referUnstained, color ) ) );
    }

    {
    const CalcElementType nearZero {matrixV.lpNorm< Eigen::Infinity >() * epsilon1};
    matrixV = ( matrixV.array() + nearZero ).matrix();
    }
  const CalcColVectorType firstOnes {CalcColVectorType::Constant( numberOfPixels, 1, 1.0 )};
  matrixV = ( firstOnes * logInputUnstained ) - matrixV.unaryExpr( CalcUnaryFunctionPointer( std::log ) );
  const auto clip = [] ( const CalcElementType &x )
    {
    return std::max( CalcElementType( 0.0 ), x );
    };
  matrixV = matrixV.unaryExpr( clip );
    {
    const CalcElementType nearZero {matrixV.lpNorm< Eigen::Infinity >() * epsilon1};
    matrixV = ( matrixV.array() + nearZero ).matrix();
    }

  // Find the associated matrixW
  CalcMatrixType matrixW;
  this->VirtanenEuclidean( matrixV, matrixW, inputH );

  // Use the matrixW with referH to compute updated values for
  // matrixV.
  matrixV = matrixW * referH;

  // Convert matrixV using exponentiation and the referUnstained pixel.
  matrixV = ( ( firstOnes * logReferUnstained ) - matrixV ).unaryExpr( CalcUnaryFunctionPointer( std::exp ) );

  OutputPixelType pixelValue;
  outputIter.GoToBegin();
  for( OutputSizeValueType pixelIndex {0}; !outputIter.IsAtEnd(); ++outputIter, ++pixelIndex )
    {
    for( Eigen::Index color = 0; color < InputImageLength; ++color )
      {
      OutputPixelHelper::value( pixelValue, color ) = matrixV( pixelIndex, color );
      }
    outputIter.Set( pixelValue );
    }
}

#if STRUCTUREPRESERVINGCOLORNORMALIZATIONFILTER_STRICT_EIGEN3_ITERATORS
template< typename TInputImage, typename TOutputImage >
template< typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols >
_Scalar *
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::begin( typename Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols > &matrix )
{
  return matrix.data();
}

template< typename TInputImage, typename TOutputImage >
template< typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols >
const _Scalar *
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::cbegin( const typename Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols > &matrix )
{
  return matrix.data();
}

template< typename TInputImage, typename TOutputImage >
template< typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols >
_Scalar *
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::end( typename Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols > &matrix )
{
  return matrix.data() + matrix.size();
}

template< typename TInputImage, typename TOutputImage >
template< typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols >
const _Scalar *
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::cend( const typename Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols > &matrix )
{
  return matrix.data() + matrix.size();
}

#else
template< typename TInputImage, typename TOutputImage >
template< typename TMatrix >
typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcElementType *
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::begin( TMatrix &matrix )
{
  return matrix.data();
}

template< typename TInputImage, typename TOutputImage >
template< typename TMatrix >
const typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcElementType *
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::cbegin( const TMatrix &matrix )
{
  return matrix.data();
}

template< typename TInputImage, typename TOutputImage >
template< typename TMatrix >
typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcElementType *
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::end( TMatrix &matrix )
{
  return matrix.data() + matrix.size();
}

template< typename TInputImage, typename TOutputImage >
template< typename TMatrix >
const typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcElementType *
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::cend( const TMatrix &matrix )
{
  return matrix.data() + matrix.size();
}
#endif

// Several members that are declared static constexpr are used by
// reference, and some compilers will thus demand that they be defined
// too.  We do that here.

template< typename TInputImage, typename TOutputImage >
constexpr typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::InputSizeValueType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::InputImageDimension;

template< typename TInputImage, typename TOutputImage >
constexpr typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::OutputSizeValueType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::OutputImageDimension;

template< typename TInputImage, typename TOutputImage >
template < typename TSizeValueType, typename TPixelType, typename TEnable >
constexpr TSizeValueType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::PixelHelper<TSizeValueType, TPixelType, TEnable>
::Length;

template< typename TInputImage, typename TOutputImage >
constexpr typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::InputSizeValueType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::InputImageLength;

template< typename TInputImage, typename TOutputImage >
constexpr typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::OutputSizeValueType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::OutputImageLength;

template< typename TInputImage, typename TOutputImage >
constexpr typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::InputSizeValueType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::NumberOfStains;

template< typename TInputImage, typename TOutputImage >
constexpr typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcElementType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::epsilon0;

template< typename TInputImage, typename TOutputImage >
constexpr typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcElementType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::epsilon1;

template< typename TInputImage, typename TOutputImage >
constexpr typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcElementType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::epsilon2;

template< typename TInputImage, typename TOutputImage >
constexpr typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::InputSizeValueType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::maxNumberOfIterations;

template< typename TInputImage, typename TOutputImage >
constexpr typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcElementType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::lambda;

} // end namespace itk

#endif // itkStructurePreservingColorNormalizationFilter_hxx
