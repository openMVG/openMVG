#include <openMVG/voctree/database.hpp>

#include <testing/testing.h>

#include <cereal/archives/binary.hpp>

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
using namespace openMVG::voctree;

int card_documents = 10;
int card_words = 12;

TEST(database, databaseIO) {

  // Create a documents vector
  vector< vector<Word> > documents_to_insert;
  documents_to_insert.resize(card_documents);
  for(int i = 0; i < documents_to_insert.size(); ++i)
  {
    documents_to_insert[i].resize(card_words);
    for(int j = 0; j < card_words; ++j)
    {
      documents_to_insert[i][j] = card_words * i + j;
    }
  }

  // Create the databases
  Database source_db( documents_to_insert.size() * documents_to_insert[0].size() ) ;
  for(int i = 0; i < documents_to_insert.size(); ++i)
    source_db.insert(i, documents_to_insert[i]);

  // Compute weights
  source_db.computeTfIdfWeights( );

  // Save the database on disk
  ofstream os("test_database.db");
  cereal::BinaryOutputArchive oarchive(os);
  oarchive(source_db);
  os.close();

  // Load the database saved on the disk
  Database reload_db;
  ifstream is("test_database.db");
  cereal::BinaryInputArchive iarchive(is);
  iarchive(reload_db);
  is.close();

  // Check databases size
  EXPECT_EQ(source_db.size(), reload_db.size());

  // Check returned matches for a given document
  for(int i = 0; i < documents_to_insert.size(); i++)
  {
    // Create match vectors
    vector<DocMatch> source_match(1), reload_match(1);
    // Query both databases with the same document
    source_db.find(documents_to_insert[i], 1, source_match);
    reload_db.find(documents_to_insert[i], 1, reload_match);
    // Check we have the same matche
    EXPECT_EQ(source_match[0].id, reload_match[0].id);
    // Check the matches scores are 0 (or near)
    EXPECT_NEAR(0, source_match[0].score, 0.001);
    EXPECT_NEAR(0, reload_match[0].score, 0.001);
  }
}

int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
