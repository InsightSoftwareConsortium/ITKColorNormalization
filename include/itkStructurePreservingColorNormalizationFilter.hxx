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
{}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::GenerateInputRequestedRegion()
{
  // { std::ostringstream mesg; mesg << "Entering GenerateInputRequestedRegion" << std::endl; std::cout << mesg.str() << std::flush; }

  // Call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // Get pointers to the input and output
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
  // { std::ostringstream mesg; mesg << "Entering BeforeThreadedGenerateData" << std::endl; std::cout << mesg.str() << std::flush; }

  // Call the superclass' implementation of this method
  Superclass::BeforeThreadedGenerateData();

  // Find input, refer, output, and make iterators for them.
  const InputImageType * const inputPtr = this->GetInput( 0 ); // image to be normalized
  const InputImageType * const referPtr = this->GetInput( 1 ); // reference image
  // For each input, make sure that it was supplied, or that we have
  // it cached already.
  itkAssertOrThrowMacro( inputPtr != nullptr || m_inputPtr != nullptr, "An image to be normalized needs to be supplied as input image #0" );
  itkAssertOrThrowMacro( referPtr != nullptr || m_referPtr != nullptr, "An reference image needs to be supplied as input image #1" );

  // For each input, if there is a supplied image and it is different
  // from what we have cached then compute stuff and cache the
  // results.
  if( inputPtr != nullptr && ( inputPtr != m_inputPtr || inputPtr->GetTimeStamp() != m_inputTimeStamp ) )
    {
    InputRegionConstIterator inIter {inputPtr, inputPtr->GetRequestedRegion()};
    CalcMatrixType inputV;
    InputPixelType inputUnstainedPixel;
    if( this->ImageToNMF( inIter, inputV, m_inputW, m_inputH, inputUnstainedPixel ) == 0 )
      {
      m_inputPtr = inputPtr;
      m_inputTimeStamp = inputPtr->GetTimeStamp();
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
    InputRegionConstIterator refIter {referPtr, referPtr->GetRequestedRegion()};
    CalcMatrixType referV;
    CalcMatrixType referW;
    if( this->ImageToNMF( refIter, referV, referW, m_referH, m_referUnstainedPixel ) == 0)
      {
      m_referPtr = referPtr;
      m_referTimeStamp = referPtr->GetTimeStamp();
      }
    else
      {
      // we failed
      m_referPtr = nullptr;
      itkAssertOrThrowMacro( m_referPtr != nullptr, "The reference image could not be processed" )
      }
    }
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::DynamicThreadedGenerateData( const OutputRegionType & outputRegion )
{
  // { std::ostringstream mesg; mesg << "Entering DynamicThreadedGenerateData" << std::endl; std::cout << mesg.str() << std::flush; }

  OutputImageType * const outputPtr = this->GetOutput();
  itkAssertOrThrowMacro( outputPtr != nullptr, "An output image needs to be supplied" )
  OutputRegionIterator outIter {outputPtr, outputRegion};

  this->NMFsToImage( m_inputW, m_inputH, m_referH, m_referUnstainedPixel, outIter );
}


template< typename TInputImage, typename TOutputImage >
int
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::ImageToNMF( InputRegionConstIterator &iter, CalcMatrixType &matrixV, CalcMatrixType &matrixW, CalcMatrixType &matrixH, InputPixelType &unstainedPixel ) const
{
  // These return values will be set later, but we can free any memory
  // associated with them up front.
  matrixV = CalcMatrixType {1, 1};
  matrixW = CalcMatrixType {1, 1};
  matrixH = CalcMatrixType {1, 1};
  // { std::ostringstream mesg; mesg << "Entering ImageToNMF" << std::endl; std::cout << mesg.str() << std::flush; }

  // To maintain locality of memory references, we are using
  // numberOfPixels as the number of rows rather than as the number of
  // columns.  With V=WH, as is standard in non-negative matrix
  // factorization, our matrices switch names and are transposed with
  // respect to the Vahadane article.  In particular, our W is a tall
  // matrix and our H is a fairly compact matrix.

  const InputSizeType size {iter.GetRegion().GetSize()};
  const unsigned int numberOfPixels = std::accumulate( size.begin(), size.end(), 1, std::multiplies< InputSizeValueType >() );

  // Find distinguishers to get a very good starting point for the subsequent
  // generic NMF algorithm.
  CalcMatrixType distinguishers;
  matrixV = CalcMatrixType {numberOfPixels, InputImageLength};
  this->ImageToMatrix( iter, matrixV );
  this->MatrixToDistinguishers( matrixV, distinguishers );

  matrixW = CalcMatrixType {numberOfPixels, NumberOfStains};
  matrixH = CalcMatrixType {NumberOfStains, InputImageLength};

  // Use the distinguishers as seeds to the non-negative matrix
  // factorization.  The published SPCN algorithm uses a Euclidean
  // penalty function, so we will hard code its use here.
  if( this->DistinguishersToNMFSeeds( distinguishers, unstainedPixel, matrixV, matrixW, matrixH ) != 0 )
    {
    return 1;                   // we failed.
    }
  // { std::ostringstream mesg; mesg << "Before VirtanenEuclidean: (log) matrixH = " << matrixH << std::endl; std::cout << mesg.str() << std::flush; }
  // { std::ostringstream mesg; mesg << "Before VirtanenEuclidean: (log) matrixW = " << matrixW << std::endl; std::cout << mesg.str() << std::flush; }
  // this->VirtanenEuclidean( matrixV, matrixW, matrixH );

  // Round off values in the response, so that numbers are quite small
  // are set to zero.
  const CalcElementType maxW = matrixW.lpNorm< Eigen::Infinity >() * 15;
  matrixW = ( ( matrixW.array() + maxW ) - maxW ).matrix();
  const CalcElementType maxH = matrixH.lpNorm< Eigen::Infinity >() * 15;
  matrixH = ( ( matrixH.array() + maxW ) - maxW ).matrix();

  // { std::ostringstream mesg; mesg << "ImageToNMF: (log) matrixH = " << matrixH << std::endl; std::cout << mesg.str() << std::flush; }
  // { std::ostringstream mesg; mesg << "ImageToNMF: (log) matrixW = " << matrixW << std::endl; std::cout << mesg.str() << std::flush; }
  return 0;
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::ImageToMatrix( InputRegionConstIterator &iter, CalcMatrixType &matrixV ) const
{
  // { std::ostringstream mesg; mesg << "Entering ImageToMatrix" << std::endl; std::cout << mesg.str() << std::flush; }

  long int pixelIndex {0};
  for( iter.GoToBegin(); !iter.IsAtEnd(); ++iter, ++pixelIndex )
    {
    InputPixelType pixelValue = iter.Get();
    for( int color {0}; color < InputImageLength; ++color )
      {
      matrixV( pixelIndex, color ) = pixelValue[color];
      }
    }
  // We do not want trouble with a value near zero (when we take its
  // logarithm) so we add a little to each value now.
  const CalcElementType nearZero {matrixV.lpNorm< Eigen::Infinity >() * epsilon};
  matrixV = ( matrixV.array() + nearZero ).matrix();
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::MatrixToDistinguishers( const CalcMatrixType &matrixV, CalcMatrixType &distinguishers ) const
{
  // { std::ostringstream mesg; mesg << "Entering MatrixToDistinguishers" << std::endl; std::cout << mesg.str() << std::flush; }

  // Keep only pixels that are bright enough.
  const CalcMatrixType brightV {this->MatrixToBrightPartOfMatrix( matrixV )};

  const CalcMatrixType normVStart {brightV};

  // We will store the row (pixel) index of each distinguishing pixel
  // in firstPassDistinguisherIndices.
  std::array< int, NumberOfStains+1 > firstPassDistinguisherIndices {-1};
  unsigned int numberOfDistinguishers {0};
  this->FirstPassDistinguishers( normVStart, firstPassDistinguisherIndices, numberOfDistinguishers );

  // Each row of secondPassDistinguisherColors is the vector of color
  // values for a distinguisher.
  CalcMatrixType secondPassDistinguisherColors {numberOfDistinguishers, brightV.cols()};
  this->SecondPassDistinguishers( normVStart, firstPassDistinguisherIndices, numberOfDistinguishers, secondPassDistinguisherColors );

  distinguishers = secondPassDistinguisherColors;
}


template< typename TInputImage, typename TOutputImage >
typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcMatrixType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::MatrixToBrightPartOfMatrix( const CalcMatrixType &matrixV ) const
{
  // { std::ostringstream mesg; mesg << "Entering MatrixToBrightPartOfMatrix" << std::endl; std::cout << mesg.str() << std::flush; }

  // A useful vector that has a 1 for each column of matrixV.
  const CalcColVectorType lastOnes {CalcColVectorType::Constant( matrixV.cols(), 1, 1.0 )};

  // We want only the brightest pixels.  Find the 80th percentile threshold.
  const CalcColVectorType brightnessOriginal {matrixV * lastOnes};
  CalcColVectorType brightnessOrdered {brightnessOriginal};
  const CalcElementType percentileLevel {0.80};
  int const quantilePosition {static_cast< int >( ( brightnessOrdered.size() - 1 ) * percentileLevel )};
  std::nth_element( Self::begin( brightnessOrdered ), Self::begin( brightnessOrdered ) + quantilePosition, Self::end( brightnessOrdered ) );
  const CalcElementType percentileThreshold {brightnessOrdered( quantilePosition )};
  // Find 70% of maximum brightness
  const CalcElementType percentageLevel {0.70};
  const CalcElementType percentageThreshold {percentageLevel * *std::max_element( Self::cbegin( brightnessOriginal ), Self::cend( brightnessOriginal ) )};

  // We will keep those pixels that pass at least one of the above
  // thresholds.
  const CalcElementType brightnessThreshold {std::min( percentileThreshold, percentageThreshold )};
  unsigned int numberOfRowsToKeep {0};
  for( int i {0} ; i < matrixV.rows(); ++i )
    {
    if( brightnessOriginal( i ) >= brightnessThreshold )
      {
      ++numberOfRowsToKeep;
      }
    }
  CalcMatrixType brightV {numberOfRowsToKeep, matrixV.cols()};
  numberOfRowsToKeep = 0;
  for( int i {0} ; i < matrixV.rows(); ++i )
    {
    if( brightnessOriginal( i ) >= brightnessThreshold )
      {
      brightV.row( numberOfRowsToKeep++ ) = matrixV.row( i );
      }
    }
  return brightV;
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::FirstPassDistinguishers( const CalcMatrixType &normVStart, std::array< int, NumberOfStains+1 > &firstPassDistinguisherIndices, unsigned int &numberOfDistinguishers ) const
{
  // { std::ostringstream mesg; mesg << "Entering FirstPassDistinguishers" << std::endl; std::cout << mesg.str() << std::flush; }

  // A useful vector that has a 1 for each column of normVStart.
  const CalcColVectorType lastOnes {CalcColVectorType::Constant( normVStart.cols(), 1, 1.0 )};
  // A useful vector that has a 1 for each row of normVStart.
  const CalcColVectorType firstOnes {CalcColVectorType::Constant( normVStart.rows(), 1, 1.0 )};

  CalcMatrixType normV {normVStart};
  numberOfDistinguishers = 0;
  bool needToRecenterMatrix = true;
  while( numberOfDistinguishers <= NumberOfStains )
    {
    // Find the next distinguishing row (pixel)
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
  { std::ostringstream mesg; mesg << "firstPassDistinguisherIndices = " << firstPassDistinguisherIndices[0] << " " << firstPassDistinguisherIndices[1] << " " << firstPassDistinguisherIndices[2] << std::endl; std::cout << mesg.str() << std::flush; }
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::SecondPassDistinguishers( const CalcMatrixType &normVStart, const std::array< int, NumberOfStains+1 > &firstPassDistinguisherIndices, const int numberOfDistinguishers,
  CalcMatrixType &secondPassDistinguisherColors ) const
{
  // { std::ostringstream mesg; mesg << "Entering SecondPassDistinguishers" << std::endl; std::cout << mesg.str() << std::flush; }

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
    // those that are at least 80% as far as the best.  (Note that
    // self could still be best, but not always.)
    const CalcColVectorType dotProducts {normV * normV.row( firstPassDistinguisherIndices[distinguisher] ).transpose()};
    const CalcElementType threshold {*std::max_element( Self::cbegin( dotProducts ), Self::cend( dotProducts ) ) * 999 / 1000};
    CalcRowVectorType cumulative {CalcRowVectorType::Constant( 1, normVStart.cols(), 0.0 )};
    int numberOfContributions {0};
    for( int row {0}; row < dotProducts.size(); ++row )
      {
      if( dotProducts( row ) >= threshold )
        {
        cumulative += normVStart.row( row );
        ++numberOfContributions;
        }
      }
    // { std::ostringstream mesg; mesg << "SecondPassDistinguishers::numberOfContributions = " << numberOfContributions << std::endl; std::cout << mesg.str() << std::flush; }
    secondPassDistinguisherColors.row( distinguisher ) = cumulative / numberOfContributions;
    }
    { std::ostringstream mesg; mesg << "secondPassDistinguisherColors = " << secondPassDistinguisherColors << std::endl; std::cout << mesg.str() << std::flush; }
}


template< typename TInputImage, typename TOutputImage >
int
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::MatrixToOneDistinguisher( const CalcMatrixType &normV, const CalcColVectorType &lastOnes ) const
{
  const CalcColVectorType lengths2 = ( normV.array() * normV.array() ).matrix() * lastOnes;
  const CalcElementType * const lengths2Begin {Self::cbegin( lengths2 )};
  const CalcElementType * const result {std::max_element( lengths2Begin, Self::cend( lengths2 ) )};
  if( *result > epsilon2 )
    {
    return std::distance( lengths2Begin, result );
    }
  else
    {
    return -1;                // Nothing left to find
    }
}


template< typename TInputImage, typename TOutputImage >
typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcMatrixType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::RecenterMatrix( const CalcMatrixType &normV, const CalcColVectorType &firstOnes, const int row ) const
{
  return normV - firstOnes * normV.row( row );
}


template< typename TInputImage, typename TOutputImage >
typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcMatrixType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::ProjectMatrix( const CalcMatrixType &normV, const int row ) const
{
  const CalcRowVectorType rowValue {normV.row( row )};
  return normV - ( normV * rowValue.transpose() ) * ( rowValue / rowValue.squaredNorm() );
}


template< typename TInputImage, typename TOutputImage >
int
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::DistinguishersToNMFSeeds( const CalcMatrixType &distinguishers, InputPixelType &unstainedPixel, CalcMatrixType &matrixV, CalcMatrixType &matrixW, CalcMatrixType &matrixH ) const
{
  // { std::ostringstream mesg; mesg << "Entering DistinguishersToNMFSeeds" << std::endl; std::cout << mesg.str() << std::flush; }

  const CalcColVectorType firstOnes {CalcColVectorType::Constant( matrixW.rows(), 1, 1.0 )};
  const CalcRowVectorType midOnes {CalcRowVectorType::Constant( 1, matrixH.rows(), 1.0 )};
  const CalcColVectorType lastOnes {CalcColVectorType::Constant( matrixV.cols(), 1, 1.0 )};

  long int unstainedIndex;
  long int hematoxylinIndex;
  long int eosinIndex;
  this->DistinguishersToColors( distinguishers, unstainedIndex, hematoxylinIndex, eosinIndex );

  // If the indices unstainedIndex, hematoxylinIndex, and eosinIndex
  // are distinct then we choose a smart starting place for the
  // generic NMF algorithm.  Otherwise, we go with a guess that is
  // reasonable.
  if( unstainedIndex != hematoxylinIndex && unstainedIndex != eosinIndex && hematoxylinIndex != eosinIndex )
    {
    const CalcRowVectorType unstainedCalcPixel {distinguishers.row( unstainedIndex )};
    for( int color {0}; color < InputImageLength; ++ color )
      {
      unstainedPixel[color] = unstainedCalcPixel( color ); // return value
      }
    const CalcRowVectorType logUnstained {unstainedCalcPixel.unaryExpr( CalcUnaryFunctionPointer( std::log ) )};
    const CalcRowVectorType logHematoxylin {logUnstained - distinguishers.row( hematoxylinIndex ).unaryExpr( CalcUnaryFunctionPointer( std::log ) )};
    const CalcRowVectorType logEosin {logUnstained - distinguishers.row( eosinIndex ).unaryExpr( CalcUnaryFunctionPointer( std::log ) )};
    // Set rows of matrixH to reflect hematoxylin and eosin.
    matrixH.row( 0 ) = logHematoxylin;
    matrixH.row( 1 ) = logEosin;
    // Convert matrixV to be the exponents of decay from the unstained
    // pixel.
    matrixV = firstOnes * logUnstained - matrixV.unaryExpr( CalcUnaryFunctionPointer( std::log ) );
    }
  else
    {
    return 1;                   // we failed
    }

  // Make sure that matrixV is non-negative.
  const auto clip = [] ( const CalcElementType &x )
    {
    return std::max( CalcElementType( 0.0 ), x );
    };
  matrixV = matrixV.unaryExpr( clip );

  // Make sure that each row of matrixH has unit magnitude, and each
  // element of matrixH is sufficiently non-negative.
  matrixH = ( ( ( matrixH.array() * matrixH.array() ).matrix() * lastOnes ).unaryExpr( CalcUnaryFunctionPointer( std::sqrt ) ) ).asDiagonal().inverse() * matrixH;
  matrixH = matrixH.unaryExpr( clip );

  // Use an approximate inverse to matrixH to get an intial value of
  // matrixW, and make sure that matrixW is sufficiently non-negative.
  const CalcMatrixType kernel {( matrixH * matrixH.transpose() ).inverse()};
  // Do we really want lambda here?!!!
  matrixW = matrixV * ( matrixH.transpose() * kernel ) - ( lambda * firstOnes ) * ( midOnes * kernel );
  matrixW = matrixW.unaryExpr( clip );
  return 0;
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::DistinguishersToColors( CalcMatrixType const &distinguishers, long int &unstainedIndex, long int &hematoxylinIndex, long int &eosinIndex ) const
{
  // { std::ostringstream mesg; mesg << "Entering DistinguishersToColors" << std::endl; std::cout << mesg.str() << std::flush; }

  // Figure out which, distinguishers are unstained (highest
  // brightness), hematoxylin (suppresses red), and eosin (suppresses
  // green).
  const CalcColVectorType lastOnes {CalcColVectorType::Constant( distinguishers.cols(), 1, 1.0 )};
  const CalcColVectorType lengths2 {( distinguishers.array() * distinguishers.array() ).matrix() * lastOnes};
  const CalcElementType * const lengths2Begin {Self::cbegin( lengths2 )};
  const CalcElementType * const unstainedIterator {std::max_element( lengths2Begin, Self::cend( lengths2 ) )};
  unstainedIndex =  std::distance( lengths2Begin, unstainedIterator );
  // For typename RGBPixel, red is suppressed by hematoxylin and is
  // color 0; green is suppressed by eosin and is color 1.  What if
  // InputPixelType is some other multi-color type ... how would we
  // find a color number that is expected to be suppressed by
  // hematoxylin and a color number that is expected to be suppressed
  // by eosin?!!!
  const CalcColVectorType redValues {distinguishers.col( 0 )};
  const CalcElementType * const redValuesBegin {Self::cbegin( redValues )};
  const CalcElementType * const hematoxylinIterator {std::min_element( redValuesBegin, Self::cend( redValues ) )};
  hematoxylinIndex = std::distance( redValuesBegin, hematoxylinIterator );
  const CalcColVectorType greenValues {distinguishers.col( 1 )};
  const CalcElementType * const greenValuesBegin {Self::cbegin( greenValues )};
  const CalcElementType * const eosinIterator {std::min_element( greenValuesBegin, Self::cend( greenValues ) )};
  eosinIndex = std::distance( greenValuesBegin, eosinIterator );
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::VirtanenEuclidean( const CalcMatrixType &matrixV, CalcMatrixType &matrixW, CalcMatrixType &matrixH ) const
{
  // { std::ostringstream mesg; mesg << "Entering VirtanenEuclidean" << std::endl; std::cout << mesg.str() << std::flush; }

  const auto clip = [] ( const CalcElementType &x )
    {
    return std::max( CalcElementType( 0.0 ), x );
    };
  const CalcColVectorType lastOnes {CalcColVectorType::Constant( matrixV.cols(), 1, 1.0 )};
  // Apply Virtanen's algorithm to iteratively improve matrixW and
  // matrixH.  Note that parentheses optimize the order of matrix
  // chain multiplications and affect the speed of this method.
  CalcMatrixType previousMatrixW {matrixW};
  unsigned int iter = 0;
  for( ; iter < maxNumberOfIterations; ++iter )
    {
    // Lasso term "lambda" insertion is possibly in a novel way.
    matrixW = ( matrixW.array()
              * ( ( ( ( matrixV * matrixH.transpose() ).array() - lambda ).unaryExpr( clip ) + epsilon2 )
                / ( ( matrixW * ( matrixH * matrixH.transpose() ) ).array() + epsilon2 ) ) ).matrix();
    matrixH = ( matrixH.array()
              * ( ( ( matrixW.transpose() * matrixV ).array() + epsilon2 )
                / ( ( ( matrixW * matrixW.transpose() ) * matrixH ).array() + epsilon2 ) ) ).matrix();
    // In lieu of rigorous Lagrange multipliers, renormalize rows of
    // matrixH to have unit magnitude.
    matrixH = ( ( ( matrixH.array() * matrixH.array() ).matrix() * lastOnes ).unaryExpr( CalcUnaryFunctionPointer( std::sqrt ) ) ).asDiagonal().inverse() * matrixH;
    if( ( iter & 15 ) == 15 )
      {
      if( ( matrixW - previousMatrixW ).lpNorm< Eigen::Infinity >() < biggerEpsilon )
        {
        break;
        }
      previousMatrixW = matrixW;
      }
    }
  // { std::ostringstream mesg; mesg << "Number of iterations = " << iter << std::endl; std::cout << mesg.str() << std::flush; }
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::VirtanenKLDivergence( const CalcMatrixType &matrixV, CalcMatrixType &matrixW, CalcMatrixType &matrixH ) const
{
  // If this method is going to get used, we may need to incorporate
  // the Lasso penalty lambda for matrixW and incorporate the Lagrange
  // multipliers to make each row of matrixH have magnitude 1.0.

  // { std::ostringstream mesg; mesg << "Entering VirtanenKLDivergence" << std::endl; std::cout << mesg.str() << std::flush; }

  // Apply Virtanen's algorithm to iteratively improve matrixW and
  // matrixH.
  const CalcRowVectorType firstOnes {CalcRowVectorType::Constant( 1, matrixV.rows(), 1.0 )};
  const CalcColVectorType lastOnes {CalcColVectorType::Constant( matrixV.cols(), 1, 1.0 )};
  CalcMatrixType previousMatrixW {matrixW};
  for( unsigned int iter = 0; iter < maxNumberOfIterations; ++iter )
    {
    matrixW = ( matrixW.array()
              * ( ( ( ( matrixV.array() + epsilon2 ) / ( ( matrixW * matrixH ).array() + epsilon2 ) + epsilon2 ).matrix() * matrixH.transpose() ).array()
                / ( ( ( matrixH * lastOnes ) * firstOnes ).transpose().array() + epsilon2 ) ) ).matrix();
    matrixH = ( matrixW.array()
              * ( ( ( matrixW.transpose() * ( ( matrixV.array() + epsilon2 ) / ( ( matrixW * matrixH ).array() + epsilon2 ) ).matrix() ).array() + epsilon2 )
                / ( ( lastOnes * ( firstOnes * matrixW ) ).transpose().array() + epsilon2 ) ) ).matrix();
    if( iter & 15 == 15 )
      {
      if( ( matrixW - previousMatrixW ).lpNorm< Eigen::Infinity >() < biggerEpsilon )
        break;
      previousMatrixW = matrixW;
      }
    }
  // { std::ostringstream mesg; mesg << "final matrixH = " << matrixH << std::endl << "final matrixW = " << matrixW << std::endl; std::cout << mesg.str() << std::flush; }
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::NMFsToImage( const CalcMatrixType &inputW, const CalcMatrixType &inputH, const CalcMatrixType &referH, const InputPixelType &referUnstained, OutputRegionIterator &out ) const
{
  // { std::ostringstream mesg; mesg << "Entering NMFsToImage" << std::endl; std::cout << mesg.str() << std::flush; }

  // We will set normalizedH to referH and then manipulate the former.
  CalcMatrixType normalizedH {referH};

  if( ( inputH * normalizedH.transpose() ).determinant() < CalcElementType( 0 ) )
    {
    // Somehow the hematoxylin and eosin rows got swaped in one of the
    // input image or reference image.  Flip them back in normalizedH to
    // get them in synch.
    static_assert( NumberOfStains == 2, "StructurePreservingColorNormalizationFilter current implementation assumes exactly two stains" );
    normalizedH.row( 0 ) = referH.row( 1 );
    normalizedH.row( 1 ) = referH.row( 0 );
    }

  // Correct for any scaling difference between normalizedH and inputH.
  const CalcColVectorType lastOnes {CalcColVectorType::Constant( normalizedH.cols(), 1, 1.0 )};
  normalizedH = {( ( ( ( inputH.array() * inputH.array() ).matrix() * lastOnes ).array() + epsilon2 )
                 / ( ( ( normalizedH.array() * normalizedH.array() ).matrix() * lastOnes ).array() + epsilon2 ) ).matrix()
    .unaryExpr( CalcUnaryFunctionPointer( std::sqrt ) ).asDiagonal() * normalizedH};

  // Use the reference image's stain colors and input image's stain
  // levels to compute what the input image would look like with the
  // reference images colors.
  CalcMatrixType newV = inputW * normalizedH;

  // Exponentiate and use as a divisor of the color of an unstained pixel.
  const long numberOfRows {newV.rows()};
  const long numberOfCols {newV.cols()};
  for( unsigned int row {0}; row < numberOfRows; ++row )
    {
    for( unsigned int col {0}; col < numberOfCols; ++col )
      {
      newV( row, col ) = referUnstained[col] / exp( newV( row, col ) );
      }
    }

  // Write the pixel values into the output region.
  InputRegionConstIterator in {m_inputPtr, m_inputPtr->GetRequestedRegion()};

  in.GoToBegin();               // for indexing of input image
  long int pixelIndex {0};      // for indexing newV.
  OutputPixelType pixelValue;
  for( out.GoToBegin(); !out.IsAtEnd(); ++out, ++in, ++pixelIndex )
    {
    // Find input index that matches this output index.
    while ( in.GetIndex() != out.GetIndex() )
      {
      ++in;
      ++pixelIndex;
      }
    for( int color {0}; color < InputImageLength; ++color )
      {
      pixelValue[color] = newV( pixelIndex, color );
      }
    out.Set( pixelValue );
    }
}


template< typename TInputImage, typename TOutputImage >
template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcElementType *
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::begin(Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols > &matrix)
{
  return matrix.data();
}

template< typename TInputImage, typename TOutputImage >
template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
const typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcElementType *
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::cbegin(const Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols > &matrix)
{
  return matrix.data();
}

template< typename TInputImage, typename TOutputImage >
template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcElementType *
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::end(Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols > &matrix)
{
  return matrix.data() + matrix.size();
}

template< typename TInputImage, typename TOutputImage >
template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
const typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcElementType *
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::cend(const Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols > &matrix)
{
  return matrix.data() + matrix.size();
}


// Several members that are declared static constexpr are used by
// reference, and some compilers will thus demand that they be defined
// too.  We do that here.
template< typename TInputImage, typename TOutputImage >
const unsigned int
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::InputImageDimension;

template< typename TInputImage, typename TOutputImage >
const unsigned int
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::OutputImageDimension;

template< typename TInputImage, typename TOutputImage >
const unsigned int
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::InputImageLength;

template< typename TInputImage, typename TOutputImage >
const unsigned int
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::OutputImageLength;

template< typename TInputImage, typename TOutputImage >
const unsigned int
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::NumberOfStains;

template< typename TInputImage, typename TOutputImage >
const typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcElementType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::biggerEpsilon;

template< typename TInputImage, typename TOutputImage >
const typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcElementType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::epsilon;

template< typename TInputImage, typename TOutputImage >
const typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcElementType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::epsilon2;

template< typename TInputImage, typename TOutputImage >
const unsigned int
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::maxNumberOfIterations;

template< typename TInputImage, typename TOutputImage >
const typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcElementType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::lambda;

} // end namespace itk

#endif // itkStructurePreservingColorNormalizationFilter_hxx
