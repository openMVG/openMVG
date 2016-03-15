#ifndef OPENMVG_VOCABULARY_TREE_DATABASE_HPP
#define OPENMVG_VOCABULARY_TREE_DATABASE_HPP

#include "vocabulary_tree.hpp"

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/map.hpp>

#include <openMVG/types.hpp>

#include <map>

namespace openMVG{
namespace voctree{

/**
 * @brief Struct representing a single database match.
 *
 * \c score is in the range [0,2], where 0 is best and 2 is worst.
 */
struct DocMatch
{
  DocId id;
  float score;

  DocMatch() { }

  DocMatch(DocId _id, float _score) : id(_id), score(_score) { }

  /// Allows sorting DocMatches in best-to-worst order with std::sort.

  bool operator<(const DocMatch& other) const
  {
    return score < other.score;
  }
};

typedef std::vector<DocMatch> DocMatches;

/**
 * @brief Class for efficiently matching a bag-of-words representation of a document (image) against
 * a database of known documents.
 */
class Database
{
public:
  /**
   * @brief Constructor
   *
   * If computing weights for a new vocabulary, \c num_words should be the size of the vocabulary.
   * If calling loadWeights(), it can be left zero.
   */
  Database(uint32_t num_words = 0);

  /**
   * @brief Insert a new document.
   *
   * @param doc_id Unique ID of the new document to insert
   * @param document The set of quantized words in a document/image.
   * \return An ID representing the inserted document.
   */
  DocId insert(DocId doc_id, const SparseHistogram& document);

  /**
   * @brief Perform a sanity check of the database by querying each document
   * of the database and finding its top N matches
   * 
   * @param[in] N The number of matches to return.
   * @param[out] matches IDs and scores for the top N matching database documents.
   */
   void sanityCheck(size_t N, std::map<size_t, DocMatches>& matches) const;

  /**
   * @brief Find the top N matches in the database for the query document.
   *
   * @param      document The query document, a set of quantized words.
   * @param      N        The number of matches to return.
   * @param[out] matches  IDs and scores for the top N matching database documents.
   */
  void find(const std::vector<Word>& document, size_t N, std::vector<DocMatch>& matches, const std::string &distanceMethod = "classic") const;
  
    /**
   * @brief Find the top N matches in the database for the query document.
   *
   * @param      query The query document, a normalized set of quantized words.
   * @param      N        The number of matches to return.
   * @param[out] matches  IDs and scores for the top N matching database documents.
   */
  void find(const SparseHistogram& query, size_t N, std::vector<DocMatch>& matches, const std::string &distanceMethod = "classic") const;

  /**
   * @brief Compute the TF-IDF weights of all the words. To be called after inserting a corpus of
   * training examples into the database.
   *
   * @param default_weight The default weight of a word that appears in none of the training documents.
   */
  void computeTfIdfWeights(float default_weight = 1.0f);

  /**
   * @brief Return the size of the database in terms of number of documents
   * @return the number of documents
   */
  size_t size() const;

  /// Save the vocabulary word weights to a file.
  void saveWeights(const std::string& file) const;
  /// Load the vocabulary word weights from a file.
  void loadWeights(const std::string& file);

  // Save weights and documents
  //void save(const std::string& file) const;
  //void load(const std::string& file);

  // Cereal serialize method

  template<class Archive>
  void serialize(Archive & archive)
  {
    archive(word_files_, word_weights_, database_);
  }

  const SparseHistogramPerImage& getSparseHistogramPerImage() const
  {
    return database_;
  }
  
private:

  struct WordFrequency
  {
    DocId id;
    uint32_t count;

    WordFrequency()
    {
    }

    WordFrequency(DocId _id, uint32_t _count) : id(_id), count(_count)
    {
    }

    // Cereal serialize methode

    template<class Archive>
    void serialize(Archive & archive)
    {
      archive(id, count);
    }
  };

  // Stored in increasing order by DocId
  typedef std::vector<WordFrequency> InvertedFile;

  /// @todo Use sorted vector?
  // typedef std::vector< std::pair<Word, float> > DocumentVector;
  
  friend std::ostream& operator<<(std::ostream& os, const SparseHistogram &dv);	

  std::vector<InvertedFile> word_files_;
  std::vector<float> word_weights_;
  SparseHistogramPerImage database_; // Precomputed for inserted documents


  /**
   * Normalize a document vector representing the histogram of visual words for a given image
   * @param[in/out] v the unnormalized histogram of visual words
   */
  void normalize(SparseHistogram& v) const;

  /**
   * @brief compute the sparse distance L1 between two histograms
   * 
   * @param v1 The first sparse histogram
   * @param v2 The second sparse histogram
   * @return the distance of the two histograms in norm L1
   */
  float sparseDistance(const SparseHistogram& v1, const SparseHistogram& v2, const std::string &distanceMethod = "classic") const;
};

}//namespace voctree
}//namespace openMVG

#endif
