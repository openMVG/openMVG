#pragma once

#include "openMVG/robust_estimation/rand_sampling.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_ransac_tools.hpp"
#include <limits>
#include <numeric>
#include <iostream>
#include <vector>
#include <iterator>

namespace openMVG {
namespace robust{

//*****************************************************
template<typename Kernel, typename Scorer>
double iterativeReweightedLeastSquares(const Kernel &kernel,
                                     const Scorer &scorer,
                                     typename Kernel::Model &best_model,
                                     std::vector<std::size_t> &inliers,
                                     double mtheta = std::sqrt(2),
                                     std::size_t numIter = 4)
{
  const std::size_t total_samples = kernel.NumSamples();
  const std::size_t min_samples = Kernel::MINIMUM_LSSAMPLES;
  double theta = scorer.getThreshold();  
  // used in the iterations to update (reduce) the threshold value
  const double deltaTetha = (mtheta*theta - theta) / (numIter-1);
  
  std::vector<std::size_t> all_samples(total_samples);
  std::iota(all_samples.begin(), all_samples.end(), 0);
  
  // find inliers from best model with threshold theta
  inliers.clear();
  scorer.Score(kernel, best_model, all_samples, &inliers, theta);
  
  if(inliers.size() < min_samples)
  {
    inliers.clear();
    std::cerr << "[IRLS] returning cause inliers.size() < min_samples" << std::endl;
    return std::numeric_limits<double>::infinity();
  }
  
  // LS model from the above inliers
  std::vector<typename Kernel::Model> models;
  kernel.FitLS(inliers, &models);
  assert(models.size()==1);   // LS fitting must always return 1 model
  
  // change threshold for refinement
  theta *= mtheta;
  
  // iterative refinement
  for(std::size_t i = 0; i < numIter; ++i)
  {
    // find inliers on the best-so-far model
    // @todo maybe inliers instead of all samples to save some computation
    inliers.clear();
    scorer.Score(kernel, models[0], all_samples, &inliers, theta);

    if(inliers.size() < min_samples)
    {
      inliers.clear();
      std::cerr << "[IRLS] returning cause inliers.size() < min_samples" << std::endl;
      return std::numeric_limits<double>::infinity();
    }
//    std::cout << "[IRLS] #" << i 
//            << " theta: " << theta
//            << " num inliers: " << inliers.size() << std::endl;
    
    // compute the weights for the inliers
    std::vector<double> weights;
    kernel.computeWeights(models[0], inliers, weights);
    
    // LS with weights on inliers
    models.clear();
    kernel.FitLS(inliers, &models, &weights);
    if(models.size() != 1)   // LS fitting must always return 1 model
    {
      std::cerr << "[IRLS] found "<< models.size() << " models, aborting..." << std::endl;
      return std::numeric_limits<double>::infinity();
    }
    
    // update the threshold
    theta -= deltaTetha;
  }
  
  assert(models.size()==1);
  best_model = models[0];
  inliers.clear();
  const double score = scorer.Score(kernel, best_model, all_samples, &inliers, theta);
  std::cout << "[IRLS] returning with num inliers: " << inliers.size() 
          << " and score " << score << std::endl;
  return score;
}



template<typename Kernel, typename Scorer>
double localOptimization(const Kernel &kernel, 
                       const Scorer &scorer, 
                       typename Kernel::Model &bestModel, 
                       std::vector<std::size_t> &bestInliers,
                       double mtheta = std::sqrt(2),
                       std::size_t numRep = 10,
                       std::size_t minSampleSize = 10)
{
  const std::size_t total_samples = kernel.NumSamples();
  const std::size_t min_samples = Kernel::MINIMUM_LSSAMPLES;
  assert((total_samples > min_samples) && 
          "[localOptimization] not enough data to estimate the model!");
  
  const double theta = scorer.getThreshold();
  
  std::vector<std::size_t> all_samples(total_samples);
  std::iota(all_samples.begin(), all_samples.end(), 0);
  
  std::size_t debugInit = 0;
  if(!bestInliers.empty())
  {
    debugInit = bestInliers.size(); 
    bestInliers.clear();
  }
  double bestScore = scorer.Score(kernel, bestModel, all_samples, &bestInliers, theta);
  if(debugInit != 0) assert(debugInit == bestInliers.size());
  
  // so far this is the best model
  std::size_t bestNumInliers = bestInliers.size();
  std::cout << "[localOptim] so far best num inliers: " << bestNumInliers << std::endl;
  std::cout << "[localOptim] so far best model:\n" << bestModel << std::endl;
  std::cout << "[localOptim] so far best score: " << bestScore << std::endl;
     
  // find inliers from best model with larger threshold t*m over all the samples
  std::vector<std::size_t> inliersBase;
  scorer.Score(kernel, bestModel, all_samples, &inliersBase, theta*mtheta);
  assert((inliersBase.size() > min_samples) && 
          "[localOptimization] not enough data in inliersBase to estimate the model!");
  
  // LS model from the above inliers
  std::vector<typename Kernel::Model> models;
  kernel.FitLS(inliersBase, &models);
  assert(models.size()==1);   // LS fitting must always return 1 model
  
  // find inliers with t again over all the samples
  inliersBase.clear();
  scorer.Score(kernel, models[0], all_samples, &inliersBase, theta);
    
  // sample of size sampleSize from the last best inliers
  const std::size_t sampleSize = std::min(minSampleSize, inliersBase.size()/2);
  if(sampleSize <= Kernel::MINIMUM_LSSAMPLES)
  {
    std::cout << "breaking cause sampleSize is " << sampleSize << std::endl;
    return bestScore;
  }
  
  // do numRep resampling + iterative LS
  for(std::size_t i = 0; i < numRep; ++i)
  {
    std::vector<std::size_t> sample;
    UniformSample(sampleSize, inliersBase, &sample);
    assert(sampleSize > Kernel::MINIMUM_LSSAMPLES);
    assert(sample.size() > Kernel::MINIMUM_LSSAMPLES);
  
    // LS estimation from the sample
    models.clear();
    kernel.FitLS(sample, &models);
    assert(models.size()==1);   // LS fitting must always return 1 model
  
    // IRLS 
    std::vector<std::size_t> inliers;
    const double score = iterativeReweightedLeastSquares(kernel, scorer, models[0], inliers);
    
    // store new best model if it is the case
    if((inliers.size() > bestNumInliers) || 
       ((inliers.size() == bestNumInliers) && (score < bestScore)))
    {
      bestNumInliers = inliers.size();
      bestScore = score;
      bestModel = models[0];
      bestInliers.swap(inliers);
      std::cout << "[localOptim] new best num inliers: " << bestNumInliers << std::endl;
    }
  }
  return bestScore;
}

//@todo make visible parameters for the optimization step
template<typename Kernel, typename Scorer>
typename Kernel::Model LO_RANSAC(
  const Kernel &kernel,
  const Scorer &scorer,
  std::vector<std::size_t> *best_inliers = NULL,
  double *best_score = NULL,
  bool bVerbose = true,
  double outliers_probability = 1e-2)
{
  assert(outliers_probability < 1.0);
  assert(outliers_probability > 0.0);
  std::size_t iteration = 0;
  const std::size_t min_samples = Kernel::MINIMUM_SAMPLES;
  const std::size_t total_samples = kernel.NumSamples();

  std::size_t max_iterations = 100;
  const std::size_t really_max_iterations = 4096;

  std::size_t bestNumInliers = 0;
  double bestInlierRatio = 0.0;
  typename Kernel::Model bestModel;

  // Test if we have sufficient points for the kernel.
  if (total_samples < min_samples) 
  {
    if (best_inliers) {
      best_inliers->clear();
    }
    return bestModel;
  }

  // In this robust estimator, the scorer always works on all the data points
  // at once. So precompute the list ahead of time [0,..,total_samples].
  std::vector<std::size_t> all_samples(total_samples);
  std::iota(all_samples.begin(), all_samples.end(), 0);

  std::vector<std::size_t> sample;
  for(iteration = 0; iteration < max_iterations; ++iteration) 
  {
    UniformSample(min_samples, total_samples, &sample);

    std::vector<typename Kernel::Model> models;
    kernel.Fit(sample, &models);

    // Compute the inlier list for each fit.
    for(std::size_t i = 0; i < models.size(); ++i) 
    {
      std::vector<std::size_t> inliers;
      double score = scorer.Score(kernel, models[i], all_samples, &inliers);
      std::cout << "sample=";
      std::sort(sample.begin(), sample.end());
      std::copy(sample.begin(), sample.end(),
                    std::ostream_iterator<std::size_t>(std::cout, ","));
      std::cout << "\nmodel " << i 
              << " e: " << score << std::endl;

      if (bestNumInliers <= inliers.size()) 
      {
        bestModel = models[i];
        //** LOCAL OPTIMIZATION
        std::cout << "Before Optim: num inliers: " << inliers.size() 
                << " score: " << score
                << " Kernel::MINIMUM_LSSAMPLES: " << Kernel::MINIMUM_LSSAMPLES 
                << std::endl;
        
        std::cout << "Model:\n" << bestModel << std::endl;
        
        if(inliers.size() > Kernel::MINIMUM_LSSAMPLES)
        {
          score = localOptimization(kernel, scorer, bestModel, inliers);
        }
        std::cout << "After Optim: num inliers: " << inliers.size()
                << " score: " << score << std::endl;
        std::cout << "Model:\n" << bestModel << std::endl;
        
        bestNumInliers = inliers.size();
        bestInlierRatio = inliers.size() / double(total_samples);

        if (best_inliers) 
        {
          best_inliers->swap(inliers);
        }

        if(bVerbose)
        {
          std::cout << " inliers=" << bestNumInliers << "/" << total_samples
                    << " (iter=" << iteration
                    << " ,i=" << i;
          std::cout << ",sample=";
          std::copy(sample.begin(), sample.end(),
                    std::ostream_iterator<std::size_t>(std::cout, ","));
          std::cout << ")";
        }
        if (bestInlierRatio) 
        {
          max_iterations = IterationsRequired(min_samples,
                                              outliers_probability,
                                              bestInlierRatio);
          // safeguard to not get stuck in a big number of iterations
          max_iterations = std::min(max_iterations, really_max_iterations);
          if(bVerbose)
            std::cout << " New max_iteration: " << max_iterations << std::endl;
        }
      }
    }
  }
  if (best_score)
    *best_score = bestNumInliers;
  
  if(bestNumInliers)
    kernel.Unnormalize(&bestModel);
  
  return bestModel;  
}





} // namespace robust
} // namespace openMVG