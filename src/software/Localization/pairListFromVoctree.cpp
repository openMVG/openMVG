#include <openMVG/voctree/database.hpp>
#include <openMVG/voctree/vocabulary_tree.hpp>
#include <openMVG/voctree/databaseIO.hpp>

#include <Eigen/Core>

#include <boost/program_options.hpp>

#include <iostream>
#include <fstream>
#include <ostream>
#include <string>
#include <set>
#include <chrono>

#include <openMVG/logger.hpp>

static const int DIMENSION = 128;

using namespace std;
namespace po = boost::program_options;


typedef openMVG::features::Descriptor<float, DIMENSION> DescriptorFloat;
typedef openMVG::features::Descriptor<unsigned char, DIMENSION> DescriptorUChar;

typedef size_t ImageID;

// just a list of doc id
typedef std::vector< ImageID > ListOfImageID;

// An ordered and unique list of doc id
typedef std::set< ImageID > OrderedListOfImageID;

// For each image ID it contains the  list of matching images
typedef std::map< ImageID, ListOfImageID> PairList;

// For each image ID it contains the ordered list of matching images
typedef std::map< ImageID, OrderedListOfImageID> OrderedPairList;

/**
 * Function that prints a PairList
 *
 * @param os The stream on which to print
 * @param pl The pair list
 * @return the stream
 */
std::ostream& operator<<(std::ostream& os, const PairList & pl)
{
  for(PairList::const_iterator plIter = pl.begin(); plIter != pl.end(); ++plIter)
  {
    os << plIter->first << " ";
    for(ImageID id : plIter->second)
    {
      os << id << " ";
    }
    os << "\n";
  }
  return os;
}

/**
 * Function that prints a OrderedPairList
 *
 * @param os The stream on which to print
 * @param pl The pair list
 * @return the stream
 */
std::ostream& operator<<(std::ostream& os, const OrderedPairList & pl)
{
  for(OrderedPairList::const_iterator plIter = pl.begin(); plIter != pl.end(); ++plIter)
  {
    os << plIter->first << " ";
    for(ImageID id : plIter->second)
    {
      os << id << " ";
    }
    os << "\n";
  }
  return os;
}

/**
 * It processes a pairlist containing all the matching images for each image ID and return
 * a similar list limited to a numMatches number of matching images and such that
 * there is no repetitions: eg if the image 1 matches with 2 in the list of image 2
 * there won't be the image 1
 *
 * @param[in] allMatches A pairlist containing all the matching images for each image of the dataset
 * @param[in] numMatches The maximum number of matching images to consider for each image
 * @param[out] matches A processed version of allMatches that consider only the first numMatches without repetitions
 */
void convertAllMatchesToPairList(const PairList &allMatches, const size_t numMatches, OrderedPairList &outPairList)
{
  outPairList.clear();

  PairList::const_iterator allIter = allMatches.begin();
  if(allIter == allMatches.end())
    return;

  // For the first image, just copy the numMatches first numMatches
  {
    ImageID currImageId = allIter->first;
    OrderedListOfImageID& bestMatches = outPairList[ currImageId ] = OrderedListOfImageID();
    for(size_t i = 0; i < std::min(allIter->second.size(), numMatches); ++i)
    {
      // avoid self-matching
      if(allIter->second[i] == currImageId)
        continue;

      bestMatches.insert(allIter->second[i]);
    }
    ++allIter;
  }

  // All other images
  for(; allIter != allMatches.end(); ++allIter)
  {
    ImageID currImageId = allIter->first;

    OrderedListOfImageID bestMatches;

    size_t numFound = 0;

    // otherwise check each element
    for(size_t i = 0; i < allIter->second.size(); ++i)
    {
      ImageID currMatchId = allIter->second[i];

      // avoid self-matching
      if(currMatchId == currImageId)
        continue;

      // if the currMatchId ID is lower than the current image ID and
      // the current image ID is not already in the list of currMatchId
      //BOOST_ASSERT( ( currMatchId < currImageId ) && ( outPairList.find( currMatchId ) != outPairList.end() ) );
      if(currMatchId < currImageId)
      {
        OrderedPairList::const_iterator currMatches = outPairList.find(currMatchId);
        if(currMatches != outPairList.end() &&
                currMatches->second.find(currImageId) == currMatches->second.end())
        {
          // then add it to the list
          bestMatches.insert(currMatchId);
          ++numFound;
        }
      }
      else if(currMatchId > currImageId)
      {
        // then add it to the list
        bestMatches.insert(currMatchId);
        ++numFound;
      }

      // if we are done stop
      if(numFound == numMatches)
        break;
    }

    // fill the output if we have matches
    if(!bestMatches.empty())
      outPairList[ currImageId ] = bestMatches;
  }
}


static const std::string programDescription =
        "This program is used to generate the pair list file to be passed to OpenMVG\n"
        "The pair list is a file containing, for each image of a dataset, the possible\n"
        "matching images in the dataset.\n"
        "The possible matching images are obtained in this case by querying a database\n"
        "and retrieving the list of the most similar images for each image in the dataset.\n"
        "Normally this works using a generic, pre-trained vocabulary tree, then a database\n"
        "is built with it using all the images of the dataset. Finally, each image of the \n"
        "dataset is queried in the database and the results are kept to form the pair list\n"
        "file.\n\n"
        "It takes as input either a list.txt file containing the a simple list of images (bundler\n"
        "format and older OpenMVG version format) or a sfm_data file (JSON) containing the\n"
        "list of images. In both cases it is assumed that the .desc to load are in the same\n"
        "directory as the input file\n"
        "For the vocabulary tree, it takes as input the input.tree (and the input.weight)\n"
        "file generated by createVoctree\n";

int main(int argc, char** argv)
{
  int verbosity = 1; ///< verbosity level
  string weightsName; ///< the filename for the voctree weights
  bool withWeights = false; ///< flag for the optional weights file
  string treeName; ///< the filename of the voctree
  string keylist; ///< the file containing the list of features
  std::size_t nbMaxDescriptors; ///< the file containing the list of features
  string outfile = "pairList.txt"; ///< the file in which to save the results
  bool withOutput = false; ///< flag for the optional output file
  size_t numImageQuery; ///< the number of matches to retrieve for each image

  po::options_description desc(programDescription);
  desc.add_options()
          ("help,h", "Print this message")
          ("verbose,v", po::value<int>(&verbosity)->default_value(1), "Verbosity level, 0 to mute")
          ("weights,w", po::value<string>(&weightsName), "Input name for the weight file, if not provided the weights will be computed on the database built with the provided set")
          ("tree,t", po::value<string>(&treeName)->required(), "Input name for the tree file")
          ("keylist,l", po::value<string>(&keylist)->required(), "Path to the list file (list.txt or sfm_data) generated by OpenMVG")
          ("nbMaxDescriptors,m", po::value<std::size_t>(&nbMaxDescriptors)->default_value(0), "Limit the number of descriptors you load per image. Zero means no limit.")
          (",r", po::value<size_t>(&numImageQuery)->default_value(10), "The number of matches to retrieve for each image (If 0 it will retrieve all the matches)")
          ("outfile,o", po::value<string>(&outfile)->default_value(outfile), "Name of the output file (pairList.txt by default)");


  po::variables_map vm;

  try
  {
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if(vm.count("help") || (argc == 1))
    {
      std::cout << desc << std::endl;
      return EXIT_SUCCESS;
    }

    po::notify(vm);
  }
  catch(boost::program_options::required_option& e)
  {
    std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
    std::cout << "Usage:\n\n" << desc << std::endl;
    return EXIT_FAILURE;
  }
  catch(boost::program_options::error& e)
  {
    std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
    std::cout << "Usage:\n\n" << desc << std::endl;
    return EXIT_FAILURE;
  }

  if(vm.count("weights"))
  {
    withWeights = true;
  }
  if(vm.count("outfile"))
  {
    withOutput = true;
  }

  //**********************************************************
  // Load the voctree
  //**********************************************************

  // Load vocabulary tree

  printf("Loading vocabulary tree\n");
  auto loadVoctree_start = std::chrono::steady_clock::now();
  openMVG::voctree::VocabularyTree<DescriptorFloat> tree(treeName);
  auto loadVoctree_elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - loadVoctree_start);
  std::cout << "tree loaded with" << endl << "\t" << tree.levels() << " levels" << std::endl
          << "\t" << tree.splits() << " branching factor" << std::endl
          << "\t in " << loadVoctree_elapsed.count() << " seconds" << std::endl;

  //**********************************************************
  // Create the database
  //**********************************************************

  POPART_COUT("Creating the database...");
  // Add each object (document) to the database
  openMVG::voctree::Database db(tree.words());

  if(withWeights)
  {
    POPART_COUT("Loading weights...");
    db.loadWeights(weightsName);
  }
  else
  {
    POPART_COUT("No weights specified, skipping...");
  }


  //*********************************************************
  // Read the descriptors and populate the database
  //*********************************************************

  POPART_COUT("Reading descriptors from " << keylist);
  auto detect_start = std::chrono::steady_clock::now();
  size_t numTotFeatures = openMVG::voctree::populateDatabase<DescriptorUChar>(keylist, tree, db, nbMaxDescriptors);
  auto detect_elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - detect_start);

  if(numTotFeatures == 0)
  {
    POPART_CERR("No descriptors loaded!!");
    return EXIT_FAILURE;
  }

  POPART_COUT("Done! " << db.getSparseHistogramPerImage().size() << " sets of descriptors read for a total of " << numTotFeatures << " features");
  POPART_COUT("Reading took " << detect_elapsed.count() << " sec");

  if(!withWeights)
  {
    // Compute and save the word weights
    POPART_COUT("Computing weights...");
    db.computeTfIdfWeights();
  }

  //**********************************************************
  // Query the database to get all the pair list
  //**********************************************************


  if(numImageQuery == 0)
  {
    // if 0 retrieve the score for all the documents of the database
    numImageQuery = db.size();
  }

  PairList allMatches;

  POPART_COUT("Query all documents");
  detect_start = std::chrono::steady_clock::now();
  // Now query each document
  #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for
  #endif
  for(std::size_t i = 0; i < db.getSparseHistogramPerImage().size(); ++i)
  {
    openMVG::voctree::SparseHistogramPerImage::const_iterator docIt = db.getSparseHistogramPerImage().cbegin();
    std::advance(docIt, i);
    std::vector<openMVG::voctree::DocMatch> matches;
    
    db.find(docIt->second, numImageQuery, matches);
//    POPART_COUT("query document " << docIt->first
//			<< " took " << detect_elapsed.count() 
//			<< " ms and has " << matches.size() 
//			<< " matches\tBest " << matches[0].id 
//			<< " with score " << matches[0].score);

    ListOfImageID idMatches;
    idMatches.reserve(matches.size());
    for(const openMVG::voctree::DocMatch& m : matches)
    {
      idMatches.push_back(m.id);
    }
    #ifdef OPENMVG_USE_OPENMP
      #pragma omp critical
    #endif
    {
      allMatches[ docIt->first ] = idMatches;
    }
  }
  detect_elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - detect_start);
  POPART_COUT("Query of all documents took " << detect_elapsed.count() << " sec.");

  //**********************************************************
  // process pair list
  //**********************************************************

  detect_start = std::chrono::steady_clock::now();
  OrderedPairList selectedPairs;

  POPART_COUT("Convert all matches to pairList");
  convertAllMatchesToPairList(allMatches, numImageQuery, selectedPairs);
  detect_elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - detect_start);
  POPART_COUT("Convert all matches to pairList took " << detect_elapsed.count() << " sec.");

  // write it to file
  std::ofstream fileout;
  fileout.open(outfile, ofstream::out);
  fileout << selectedPairs;
  fileout.close();

  POPART_COUT("pairList exported in: " << outfile);
  return EXIT_SUCCESS;
}
