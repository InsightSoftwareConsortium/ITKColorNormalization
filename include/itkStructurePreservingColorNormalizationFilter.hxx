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
::DynamicThreadedGenerateData( const OutputRegionType & outputRegion )
{
  // Handle the case that the outputRegion size does not match the
  // input region size!!!

  // Do we have to worry that input, refer, and output might not have the same
  // physical measurements for a pixel?!!!

  // Find input, refer, output, and make iterators for them.
  const InputImageType * const inputPtr = this->GetInput( 0 ); // image to be normalized
  const InputImageType * const referPtr = this->GetInput( 1 ); // reference image
  OutputImageType * const outputPtr = this->GetOutput();
  itkAssertOrThrowMacro( inputPtr != nullptr, "An image to be normalized needs to be supplied as input image #0" );
  itkAssertOrThrowMacro( referPtr != nullptr, "An reference image needs to be supplied as input image #1" );
  itkAssertOrThrowMacro( outputPtr != nullptr, "An output image needs to be supplied" )

  InputRegionConstIterator in {inputPtr, inputPtr->GetRequestedRegion()};
  InputRegionConstIterator ref {referPtr, referPtr->GetRequestedRegion()};
  OutputRegionIterator out {outputPtr, outputRegion};
  {std::ostringstream mesg; mesg << "input region size = " << inputPtr->GetRequestedRegion().GetSize() << ", refer region size = " << referPtr->GetRequestedRegion().GetSize() << ", output region size = " << outputPtr->GetRequestedRegion().GetSize() << std::endl; std::cout << mesg.str(); } // Remove me!!!
  // {std::ostringstream mesg; mesg << "input region index = " << inputPtr->GetRequestedRegion().GetIndex() << ", refer region index = " << referPtr->GetRequestedRegion().GetIndex() << ", output region index = " << outputPtr->GetRequestedRegion().GetIndex() << std::endl; std::cout << mesg.str(); } // Remove me!!!

  // Non-negative matrix factorization produces matrices CalcMatrixType values
  // satisfying V=WH and a InputPixelType color for an unstained pixel.
  CalcMatrixType inputV;
  CalcMatrixType inputW;
  CalcMatrixType inputH;
  InputPixelType inputUnstainedPixel;
  this->ImageToNMF( in, inputV, inputW, inputH, inputUnstainedPixel );

  CalcMatrixType referV;
  CalcMatrixType referW;
  CalcMatrixType referH;
  InputPixelType referUnstainedPixel;
  this->ImageToNMF( ref, referV, referW, referH, referUnstainedPixel );

  this->NMFsToImage( inputW, inputH, referH, referUnstainedPixel, out );
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::GenerateInputRequestedRegion()
{
  // Is this implemented correctly?!!!.  Does
  // DynamicThreadedGenerateData do the right thing with respect to
  // using the largest possible region?!!!  Should GenerateData be
  // implemented too/instead?

  // Call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // Get pointers to the input and output
  InputImageType * inputPtr = const_cast< InputImageType * >( this->GetInput( 0 ) );
  InputImageType * referPtr = const_cast< InputImageType * >( this->GetInput( 1 ) );

  itkAssertOrThrowMacro( inputPtr != nullptr, "An image to be normalized needs to be supplied as input image #0" );
  itkAssertOrThrowMacro( referPtr != nullptr, "An reference image needs to be supplied as input image #1" );

  inputPtr->SetRequestedRegionToLargestPossibleRegion();
  referPtr->SetRequestedRegionToLargestPossibleRegion();
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::ImageToNMF( InputRegionConstIterator &iter, CalcMatrixType &matrixV, CalcMatrixType &matrixW, CalcMatrixType &matrixH, InputPixelType &unstainedPixel ) const
{
  // { std::ostringstream mesg; mesg << "Entering ImageToNMF" << std::endl; std::cout << mesg.str(); }

  const InputSizeType size {iter.GetRegion().GetSize()};
  const unsigned int numberOfPixels = std::accumulate( size.begin(), size.end(), 1, std::multiplies< InputSizeValueType >() );
  // To maintain locality of memory references, we are using numberOfPixels as
  // the number of rows rather than as the number of columns.  With V=WH, as is
  // standard in non-negative matrix factorization, our matrices switch names
  // and are transposed with respect to Vahadane.  In particular, our W is a
  // tall matrix and our H is a fairly compact matrix.
  matrixV = CalcMatrixType {numberOfPixels, InputImageLength};
  matrixW = CalcMatrixType {numberOfPixels, NumberOfStains};
  matrixH = CalcMatrixType {NumberOfStains, InputImageLength};

  // A vector that has a 1 for each row of matrixV.
  const CalcVectorType longOnes {matrixV.rows(), 1.0};

  // Find distinguishers to get a very good starting point for the subsequent
  // generic NMF algorithm.
  CalcMatrixType distinguishers;
  this->ImageToMatrix( iter, matrixV );
  this->MatrixToDistinguishers( matrixV, distinguishers );

  // Use the distinguishers as seeds to the non-negative matrix
  // factorization.  The published SPCN algorithm uses a Euclidean
  // penalty function, so we will hard code its use here.
  this->DistinguishersToNMFSeeds( distinguishers, longOnes, unstainedPixel, matrixV, matrixW, matrixH );
  this->VirtanenEuclid( matrixV, matrixW, matrixH );
  // this->VirtanenKLDivergence( matrixV, matrixW, matrixH );

  // Round off values in the response, so that numbers of order 1e-16
  // are set to zero.
  const CalcElementType maxW = matrixW.array_inf_norm() * 15;
  matrixW += maxW;
  matrixW -= maxW;
  const CalcElementType maxH = matrixH.array_inf_norm() * 15;
  matrixH += maxH;
  matrixH -= maxH;

  { std::ostringstream mesg; mesg << "matrixH = " << matrixH << std::endl; std::cout << mesg.str(); } // Remove me!!!
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::ImageToMatrix( InputRegionConstIterator &iter, CalcMatrixType &matrixV ) const
{
  // { std::ostringstream mesg; mesg << "Entering ImageToMatrix" << std::endl; std::cout << mesg.str(); }
  int pixelIndex {0};
  for( iter.GoToBegin(); !iter.IsAtEnd(); ++iter, ++pixelIndex )
    {
    InputPixelType pixelValue = iter.Get();
    for( int color {0}; color < InputImageLength; ++color )
      {
      matrixV.put( pixelIndex, color, pixelValue[color] );
      }
    }
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::MatrixToDistinguishers( const CalcMatrixType &matrixV, CalcMatrixType &distinguishers ) const
{
  // { std::ostringstream mesg; mesg << "Entering MatrixToDistinguishers" << std::endl; std::cout << mesg.str(); }

  // Note that in the below, parentheses optimize the order of matrix
  // chain multiplications and are very important for speed.

  // A useful value that we will use several times.
  const CalcMatrixType kernel {matrixV.transpose() * matrixV};
  // A vector that has a 1 for each column of matrixV = each column of matrixH.
  const CalcVectorType shortOnes {matrixV.cols(), 1.0};
  // A vector that has a 1 for each row of matrixV = each row of matrixW.
  const CalcVectorType longOnes {matrixV.rows(), 1.0};

  // Build a diagonal matrix from the row sums of
  // matrixV*matrixV.transpose().  Take its matrix inverse.
  const CalcDiagMatrixType rowNorms = static_cast< CalcDiagMatrixType >( matrixV * ( matrixV.transpose() * longOnes ) ).invert_in_place();

  // normV*matrixV.transpose() is a matrix that we will implicitly
  // re-center and then project out vectors from, all in the search
  // for the next distinguisher.  To start, we want each row of
  // normV*matrixV.transpose() to sum to 1.0.  This computation is
  // faster than a product of generic matrices because one factor is a
  // CalcDiagMatrixType.
  const CalcMatrixType normVStart {rowNorms * matrixV};

  // Pass #1
  //
  // We seek a row of matrixV for unstained and for each stain.
  std::array< int, NumberOfStains+1 > firstPassDistinguishers {-1};
  unsigned int numberOfDistinguishers {0};
  {
    // Find the first distinguishing row (pixel).
    CalcMatrixType normV {normVStart}; // want and have a copy constructor, not a reference
    firstPassDistinguishers[0] = this->MatrixToOneDistinguisher( kernel, shortOnes, normV );
    // If we found a distinguisher then prepare to look for the next
    // distinguisher.
    if( firstPassDistinguishers[0] >= 0 )
      {
      ++numberOfDistinguishers;
      normV = this->RecenterMatrix( longOnes, normV, firstPassDistinguishers[0] ); // (need to execute only if NumberOfStains > 0.)
      for( int nextDistinguisher {1}; nextDistinguisher <= NumberOfStains; ++nextDistinguisher )
        {
        // Find the next distinguishing row (pixel)
        firstPassDistinguishers[nextDistinguisher] = this->MatrixToOneDistinguisher( kernel, shortOnes, normV );
        // If we found a distinguisher then prepare to look for the next
        // distinguisher.
        if( firstPassDistinguishers[nextDistinguisher] < 0 )
          {
          break;
          }
        ++numberOfDistinguishers;
        if( nextDistinguisher < NumberOfStains )
          {
          normV = this->ProjectMatrix( kernel, normV, firstPassDistinguishers[nextDistinguisher] );
          }
        }
      }
  }

  // Pass #2
  //
  CalcMatrixType secondPassDistinguishers {numberOfDistinguishers, matrixV.cols()};
  {
    // Note that if ever numberOfDistinguishers is large, this can
    // instead be done in N lg N time rather than N^2 time, though it
    // is more complicated.
    for( int distinguisher {0}; distinguisher < numberOfDistinguishers; ++distinguisher )
      {
      CalcMatrixType normV {normVStart};
      bool needToRecenterMatrix = true;
      for( int nextDistinguisher {0}; nextDistinguisher < numberOfDistinguishers; ++nextDistinguisher )
        {
        // skip if self
        if( nextDistinguisher == distinguisher )
          {
          continue;
          }
        if( needToRecenterMatrix )
          {
          normV = this->RecenterMatrix( longOnes, normV, firstPassDistinguishers[nextDistinguisher] );
          needToRecenterMatrix = false;
          }
        else
          {
          normV = this->ProjectMatrix( kernel, normV, firstPassDistinguishers[nextDistinguisher] );
          }
        }
      // We have sent all distinguishers except self to the origin.
      // Whatever is far from the origin in the same direction as self
      // is a good replacement for self.  We will take an average
      // among those that are at least 80% as far as the best.  (Note
      // that often self will be best, but not always.)
      const CalcVectorType dotProducts {normV * ( kernel * normV.get_row( firstPassDistinguishers[distinguisher] ) )};
      const CalcElementType threshold {*std::max_element( dotProducts.begin(), dotProducts.end() ) * 8 / 10};
      CalcVectorType cumulative {matrixV.cols(), 0};
      int NumberOfContributions {0};
      for( int row {0}; row < dotProducts.size(); ++row )
        {
        if( dotProducts[row] >= threshold )
          {
          cumulative += matrixV.get_row( row );
          ++NumberOfContributions;
          }
        }
      secondPassDistinguishers.set_row( distinguisher, cumulative / NumberOfContributions );
      }
  }

  distinguishers = secondPassDistinguishers;
}


template< typename TInputImage, typename TOutputImage >
int
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::MatrixToOneDistinguisher( const CalcMatrixType &kernel, const CalcVectorType &shortOnes, const CalcMatrixType &normV ) const
{
  const CalcVectorType lengths2 = element_product( normV * kernel, normV ) * shortOnes;
  const CalcElementType * const result {std::max_element( lengths2.begin(), lengths2.end() )};
  if( *result > epsilon2 )
    {
    return std::distance( lengths2.begin(), result );
    }
  else
    {
    return -1;                // Nothing left to find
    }
}


template< typename TInputImage, typename TOutputImage >
typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcMatrixType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::RecenterMatrix( const CalcVectorType &longOnes, const CalcMatrixType &normV, const int row ) const
{
  return normV - outer_product( longOnes, normV.get_row( row ) );
}


template< typename TInputImage, typename TOutputImage >
typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcMatrixType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::ProjectMatrix( const CalcMatrixType &kernel, const CalcMatrixType &normV, const int row ) const
{
  const CalcVectorType rowValue {normV.get_row( row )};
  const CalcElementType squared_magnitude = dot_product( kernel * rowValue, rowValue );
  return normV - outer_product( ( normV * ( kernel * rowValue ) ), rowValue / squared_magnitude );
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::DistinguishersToNMFSeeds( const CalcMatrixType &distinguishers, const CalcVectorType &longOnes,
  InputPixelType &unstainedPixel, CalcMatrixType &matrixV, CalcMatrixType &matrixW, CalcMatrixType &matrixH ) const
{
  // { std::ostringstream mesg; mesg << "Entering DistinguishersToNMFSeeds" << std::endl; std::cout << mesg.str(); }
  int unstainedIndex;
  int hematoxylinIndex;
  int eosinIndex;
  this->DistinguishersToColors( distinguishers, unstainedIndex, hematoxylinIndex, eosinIndex );

  // If the indices unstainedIndex, hematoxylinIndex, and eosinIndex
  // are distinct then we choose a smart starting place for the
  // generic NMF algorithm.  Otherwise, we go with a guess that is
  // reasonable.
  if( unstainedIndex != hematoxylinIndex && unstainedIndex != eosinIndex && hematoxylinIndex != eosinIndex )
    {
    const CalcVectorType unstainedCalcPixel {distinguishers.get_row( unstainedIndex )};
    for( int color {0}; color < InputImageLength; ++ color )
      {
      unstainedPixel[color] = unstainedCalcPixel[color]; // return value
      }
    const CalcVectorType logUnstained {unstainedCalcPixel.apply( std::log )};
    const CalcVectorType logHematoxylin {logUnstained - distinguishers.get_row( hematoxylinIndex ).apply( std::log )};
    const CalcVectorType logEosin {logUnstained - distinguishers.get_row( eosinIndex ).apply( std::log )};
    matrixH.set_row( 0, logHematoxylin );
    matrixH.set_row( 1, logEosin );
    matrixV = outer_product( longOnes, logUnstained ) - matrixV.apply( std::log );
    }
  else
    {
    // Choose a reasonable unstained point for matrixV.  Make matrixW
    // and matrixH reasonable, or at least within the feasible region.
    CalcVectorType unstainedCalcPixel;
    const auto max = [] ( CalcElementType a, CalcElementType b ) -> CalcElementType
      {
      return std::max( a, b );
      };
    for( int color {0}; color < InputImageLength; ++ color )
      {
      CalcVectorType column {matrixV.get_column( color )};
      unstainedCalcPixel[color] = std::accumulate( column.begin(), column.end(), CalcElementType( 0 ), max );
      unstainedPixel[color] = unstainedCalcPixel[color]; // return value
      }
    const CalcVectorType logUnstained {unstainedCalcPixel.apply( std::log )};
    matrixV = outer_product( longOnes, logUnstained ) - matrixV.apply( std::log );
    std::array< CalcElementType, InputImageLength > logHematoxylin {2.772589, 1.340485, 0.7744928};
    std::array< CalcElementType, InputImageLength > logEosin {0.2479587, 2.496741, 0.6509144};
    for( int color {0}; color < InputImageLength; ++ color )
      {
      matrixH.put( 0, color, logHematoxylin[color] );
      matrixH.put( 1, color, logEosin[color] );
      }
    }
  // Note that we are not initializing matrixW.  Its will be chosen as
  // the first step in this->NMFSeedsToNMFSolution()
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::DistinguishersToColors( CalcMatrixType const &distinguishers, int &unstainedIndex, int &hematoxylinIndex, int &eosinIndex ) const
{
  // { std::ostringstream mesg; mesg << "Entering DistinguishersToColors" << std::endl; std::cout << mesg.str(); }

  // Figure out which, distinguishers are unstained, hematoxylin (suppresses
  // red), and eosin (suppresses green).
  const CalcVectorType ones {distinguishers.cols(), 1.0};
  const CalcVectorType lengths2 {element_product( distinguishers, distinguishers ) * ones};
  const typename CalcVectorType::const_iterator unstainedIterator {std::max_element( lengths2.begin(), lengths2.end() )};
  unstainedIndex = std::distance( lengths2.begin(), unstainedIterator );
  // For typename RGBPixel, red is suppressed by hematoxylin and is
  // color 0; green is suppressed by eosin and is color 1.  What if
  // InputPixelType is some other multi-color type ... how would we
  // find a color number that is expected to be suppressed by
  // hematoxylin and a color number that is expected to be suppressed
  // by eosin?
  const CalcVectorType redValues {distinguishers.get_column( 0 )};
  const typename CalcVectorType::const_iterator hematoxylinIterator {std::min_element( redValues.begin(), redValues.end() )};
  hematoxylinIndex = std::distance( redValues.begin(), hematoxylinIterator );
  const CalcVectorType greenValues {distinguishers.get_column( 1 )};
  const typename CalcVectorType::const_iterator eosinIterator {std::min_element( greenValues.begin(), greenValues.end() )};
  eosinIndex = std::distance( greenValues.begin(), eosinIterator );
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::VirtanenEuclid( const CalcMatrixType &matrixV, CalcMatrixType &matrixW, CalcMatrixType &matrixH ) const
{
  // { std::ostringstream mesg; mesg << "Entering VirtanenEuclid with initial matrixH = " << matrixH << std::endl; std::cout << mesg.str(); }
  matrixW = 1;
  for( unsigned int iter = 0; iter < numberOfIterations; ++iter )
    {
    matrixW = element_product( matrixW, element_quotient( matrixV * matrixH.transpose() + epsilon2, matrixW * ( matrixH * matrixH.transpose() ) + epsilon2 ) );
    }
  for( unsigned int iter = 0; iter < numberOfIterations; ++iter )
    {
    matrixW = element_product( matrixW, element_quotient( matrixV * matrixH.transpose() + epsilon2, matrixW * ( matrixH * matrixH.transpose() ) + epsilon2 ) );
    matrixH = element_product( matrixH, element_quotient( matrixW.transpose() * matrixV + epsilon2, ( matrixW.transpose() * matrixW ) * matrixH + epsilon2 ) );
    }
  // { std::ostringstream mesg; mesg << "final matrixH = " << matrixH << std::endl; std::cout << mesg.str(); }
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::VirtanenKLDivergence( const CalcMatrixType &matrixV, CalcMatrixType &matrixW, CalcMatrixType &matrixH ) const
{
  // { std::ostringstream mesg; mesg << "Entering VirtanenKLDivergence with initial matrixH = " << matrixH << std::endl; std::cout << mesg.str(); }
  const CalcVectorType longOnes {matrixV.rows(), 1.0};
  const CalcVectorType shortOnes {matrixV.cols(), 1.0};
  matrixW = 1;
  for( unsigned int iter = 0; iter < numberOfIterations; ++iter )
    {
    matrixW = element_product( matrixW, element_quotient( element_quotient( matrixV + epsilon2, matrixW * matrixH + epsilon2 ) * matrixH.transpose() + epsilon2,
        outer_product( longOnes, shortOnes * matrixH.transpose() ) + epsilon2 ) );
    }
  for( unsigned int iter = 0; iter < numberOfIterations; ++iter )
    {
    matrixW = element_product( matrixW, element_quotient( element_quotient( matrixV + epsilon2, matrixW * matrixH + epsilon2 ) * matrixH.transpose() + epsilon2,
        outer_product( longOnes, shortOnes * matrixH.transpose() ) + epsilon2 ) );
    matrixH = element_product( matrixH, element_quotient( matrixW.transpose() * element_quotient( matrixV + epsilon2, matrixW * matrixH + epsilon2 ) + epsilon2,
        outer_product( matrixW.transpose() * longOnes, shortOnes ) + epsilon2 ) );
    }
  // { std::ostringstream mesg; mesg << "final matrixH = " << matrixH << std::endl << "final matrixW = " << matrixW << std::endl; std::cout << mesg.str(); }
}


template< typename TInputImage, typename TOutputImage >
void
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::NMFsToImage( const CalcMatrixType &inputW, const CalcMatrixType &inputH, const CalcMatrixType &referH, const InputPixelType &referUnstained, OutputRegionIterator &out ) const
{
  // { std::ostringstream mesg; mesg << "Entering NMFsToImage" << std::endl; std::cout << mesg.str(); }

  // We will set normalizedH to referH and then manipulate the former.
  CalcMatrixType normalizedH {referH};

  if( vnl_determinant( inputH * normalizedH.transpose() ) < CalcElementType( 0 ) )
    {
    // Somehow the hematoxylin and eosin rows got swaped in one of the
    // input image or reference image.  Flip them back in normalizedH to
    // get them in synch.
    static_assert( NumberOfStains == 2, "StructurePreservingColorNormalizationFilter current implementation assumes exactly two stains" );
    normalizedH.set_row( 0, referH.get_row( 1 ) );
    normalizedH.set_row( 1, referH.get_row( 0 ) );
    }

  // Correct for any scaling difference between normalizedH and inputH.
  const CalcVectorType shortOnes {normalizedH.cols(), CalcElementType( 1 )};
  normalizedH = {static_cast< CalcDiagMatrixType >( element_quotient( element_product( inputH, inputH ) * shortOnes + epsilon2,
        element_product( normalizedH, normalizedH ) * shortOnes +epsilon2 ).apply( std::sqrt ) ) * normalizedH};

  // Use the reference image's stain colors and input image's stain
  // levels to compute what the input image would look like with the
  // reference images colors.
  CalcMatrixType newV = inputW * normalizedH;

  // Exponentiate and subtract from maximum intensity of the reference
  // image.
  const unsigned int numberOfRows {newV.rows()};
  const unsigned int numberOfCols {newV.cols()};
  for( unsigned int row {0}; row < numberOfRows; ++row )
    {
    for( unsigned int col {0}; col < numberOfCols; ++col )
      {
      newV.put( row, col, referUnstained[col] - exp( newV.get( row, col ) ) );
      }
    }

  // Write the pixel values into the output region.
  int pixelIndex {0};
  for( out.GoToBegin(); !out.IsAtEnd(); ++out, ++pixelIndex )
    {
    OutputPixelType pixelValue;
    for( int color {0}; color < InputImageLength; ++color )
      {
      pixelValue[color] = newV.get( pixelIndex, color );
      }
    out.Set( pixelValue );
    }
}


// epsilon, epsilon2, and lambda are explicitly defined here, even
// though they are declared and initialized as static constexpr
// members, because they are passed by reference in some versions of
// the implementation.
template< typename TInputImage, typename TOutputImage >
const typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcElementType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::epsilon;

template< typename TInputImage, typename TOutputImage >
const typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcElementType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::epsilon2;

template< typename TInputImage, typename TOutputImage >
const typename StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >::CalcElementType
StructurePreservingColorNormalizationFilter< TInputImage, TOutputImage >
::lambda;

} // end namespace itk

#endif // itkStructurePreservingColorNormalizationFilter_hxx
