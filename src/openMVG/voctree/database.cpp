#include "database.hpp"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/tail.hpp>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <boost/format.hpp>

namespace openMVG{
namespace voctree{

Database::Database( uint32_t num_words )
: word_files_( num_words ),
word_weights_( num_words, 1.0f ) { }

DocId Database::insert( const std::vector<Word>& document )
{
	/// @todo Evaluate whether sorting words makes much difference in speed
	DocId doc_id = database_vectors_.size( );

	// For each word, retrieve its inverted file and increment the count for doc_id.
	for ( std::vector<Word>::const_iterator it = document.begin( ), end = document.end( ); it != end; ++it )
	{
		Word word = *it;
		InvertedFile& file = word_files_[word];
		if ( file.empty( ) || file.back( ).id != doc_id )
			file.push_back( WordFrequency( doc_id, 1 ) );
		else
			file.back( ).count++;
	}

	// Precompute the document vector to compare queries against.
	database_vectors_.resize( doc_id + 1 );
	computeVector( document, database_vectors_.back( ) );

	return doc_id;
}
/**
 * \brief Find the top N matches in the database for the query document.
 *
 * \param[in]  document The query document, a set of quantized words.
 * \param[in]  N        The number of matches to return.
 * \param[out] matches  IDs and scores for the top N matching database documents.
 */
void Database::find( const std::vector<Word>& document, size_t N, std::vector<Match>& matches ) const
{
	DocumentVector query;
	// from the list of visual words associated with each feature in the document/image
	// generate the (sparse) histogram of the visual words 
	computeVector( document, query );

	// Accumulate the best N matches
	using namespace boost::accumulators;
	typedef tag::tail<left> bestN_tag;
	accumulator_set<Match, features<bestN_tag> > acc( bestN_tag::cache_size = N );

	/// @todo Try only computing distances against documents sharing at least one word
	for ( DocId i = 0; i < ( DocId ) database_vectors_.size( ); ++i )
	{
		// for each document/image in the database compute the distance between the 
		// histograms of the query image and the others
		float distance = sparseDistance( query, database_vectors_[i] );
		acc( Match( i, distance ) );
	}

	// extract the best N
	extractor<bestN_tag> bestN;
	matches.resize( std::min( N, database_vectors_.size( ) ) );
	std::copy( bestN( acc ).begin( ), bestN( acc ).end( ), matches.begin( ) );
}

DocId Database::findAndInsert( const std::vector<Word>& document, size_t N, std::vector<Match>& matches )
{
	/// @todo Can this be accelerated? Could iterate over words only once?
	find( document, N, matches );
	return insert( document );
}

/**
 * \brief Compute the TF-IDF weights of all the words. To be called after inserting a corpus of
 * training examples into the database.
 *
 * \param default_weight The default weight of a word that appears in none of the training documents.
 */
void Database::computeTfIdfWeights( float default_weight )
{
	float N = ( float ) database_vectors_.size( );
	size_t num_words = word_files_.size( );
	for ( size_t i = 0; i < num_words; ++i )
	{
		size_t Ni = word_files_[i].size( );
		if ( Ni != 0 )
			word_weights_[i] = std::log( N / Ni );
		else
			word_weights_[i] = default_weight;
	}
}

void Database::saveWeights( const std::string& file ) const
{
	std::ofstream out( file.c_str( ), std::ios_base::binary );
	uint32_t num_words = word_weights_.size( );
	out.write( ( char* ) ( &num_words ), sizeof (uint32_t ) );
	out.write( ( char* ) ( &word_weights_[0] ), num_words * sizeof (float ) );
}

void Database::loadWeights( const std::string& file )
{
	std::ifstream in;
	in.exceptions( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );

	try
	{
		in.open( file.c_str( ), std::ios_base::binary );
		uint32_t num_words = 0;
		in.read( ( char* ) ( &num_words ), sizeof (uint32_t ) );
		word_files_.resize( num_words ); // Inverted files start out empty
		word_weights_.resize( num_words );
		in.read( ( char* ) ( &word_weights_[0] ), num_words * sizeof (float ) );
	}
	catch ( std::ifstream::failure& e )
	{
		throw std::runtime_error( ( boost::format( "Failed to load vocabulary weights file '%s'" ) % file ).str( ) );
	}
}
/**
 * Given a list of visual words associated to the features of a document it computes the 
 * vector of unique weighted visual words
 * 
 * @param[in] document a list of (possibly repeated) visual words
 * @param[out] v the vector of visual words (ie the visual word histogram of the image)
 */
void Database::computeVector( const std::vector<Word>& document, DocumentVector& v ) const
{
//	for each visual word in the list
	for ( std::vector<Word>::const_iterator it = document.begin( ), end = document.end( ); it != end; ++it )
	{
		// update its weighted count inside the map
		// the map v contains only the visual words that are associated to some features
		// the visual words in v are unique unlikely the document
		Word word = *it;
		v[word] += word_weights_[word];
	}
	normalize( v );
}

/**
 * Normalize a document vector representing the histogram of visual words for a given image
 * 
 * @param[in/out] v the unnormalized histogram of visual words
 */
void Database::normalize( DocumentVector& v )
{
	float sum = 0.0f;
	for ( DocumentVector::iterator i = v.begin( ), ie = v.end( ); i != ie; ++i )
		sum += i->second;
	float inv_sum = 1.0f / sum;
	for ( DocumentVector::iterator i = v.begin( ), ie = v.end( ); i != ie; ++i )
		i->second *= inv_sum;
}

float Database::sparseDistance( const DocumentVector& v1, const DocumentVector& v2 )
{
	float distance = 0.0f;
	DocumentVector::const_iterator i1 = v1.begin( ), i1e = v1.end( );
	DocumentVector::const_iterator i2 = v2.begin( ), i2e = v2.end( );

	while ( i1 != i1e && i2 != i2e )
	{
		if ( i2->first < i1->first )
		{
			distance += i2->second;
			++i2;
		}
		else if ( i1->first < i2->first )
		{
			distance += i1->second;
			++i1;
		}
		else
		{
			distance += fabs( i1->second - i2->second );
			++i1;
			++i2;
		}
	}

	while ( i1 != i1e )
	{
		distance += i1->second;
		++i1;
	}

	while ( i2 != i2e )
	{
		distance += i2->second;
		++i2;
	}

	return distance;
}

/**
* \brief Return the size of the database in terms of number of documents
* @return the number of documents
*/
size_t Database::size() const
{
	return database_vectors_.size( );
}

} //namespace voctree
} //namespace openMVG
