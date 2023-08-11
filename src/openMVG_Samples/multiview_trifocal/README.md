## Working document

https://docs.google.com/document/d/1Ozvt64zAWNRYE16nICrJobDnV_ByWbcVx5lGjPHqSiQ/edit


## Main branch
original trifocal branch  created by pierre: develop_keypoint_orientation_sfm @ official
working branch: develop_keypoint_orientation_sfm_tmp
  - this working branch merges the develop branch periodically
  - when everything is working, we just make a nice commit from scratch to
    OpenMVG with Gabriel as author, with pierre as coauthors


## Misc. Notes

src/software 
 - Complete software(s) build on openMVG libraries

## Inclusion into OpenMVG

### Next thing to work on

MakeInitialPair3D
sequential_SfM.cpp:188

Make trifocal one


## main() ##

main_SfM.cpp is our main

src/software/SfM/SfM_SequentialPipeline.py.in
  - calls all the executables from python.
  - main: main_SfM.cpp 
    - main funciton: ReconstructionEngine::Process() 
        SfMEngine::INCREMENTALV:
          - we are mainly targeting this one
          
        SfMEngine::INCREMENTALV2:
          - this one allows to set the triples and initializer (see\ref{pierre:comment})
  - also see main_GeometricFilter.cpp


Compute matching pairs -->  ComputeMatches               --> GeometricFilter -->
I: sfm_data.json            I: sfm_data.json, pairs.bin      I: sfm_data.json, matches.putative.bin
O: pairs.bin                O: matches.putative.bin          O: matches.f.bin

-->  main_SfM
     I: sfm_data.json, match_dir (matches.f.bin)
     O: reconstruction_sequential/

GeometricFilter 
- The toplevel .py shows it uses a fundamental matrix filtering model
- Why not essential matrix when intrinsics available?
    - TODO: investigate where intrinsics are provided, if we can fix it, and
      whether intrinsics are actually optimized / estimated
        - sfm_data has intrinsics and views (see \ref{sec:views} below)

### Sequential reconstruction engine v1

SequentialSfMReconstructionEngine::Constructor()
- initializes list of remaining images to reconstruct as all images
- TODO: Idea: remaining images as priority queue based on number of matches

SequentialSfMReconstructionEngine::Process()
- InitLandmarkTracks
  - map_tracks_: where tracks are stored
  - tracks::TracksBuilder::Build
  - tracks::TracksBuilder::Filter
  - shared_track_visibility_helper_()
    - Helper to compute if some image have some track in common
  - Analogous to TrifocalSampleApp::ComputeTracks()
    - Build
    - Filter(3)

### Useful structures:

#### Landmark 
  - a associacao de um ponto em 3D e suas observacoes em imagens: `sfm_landmark.hpp`
  - TODO: improve with tangents
  - TODO: in trifocal RANSAC project tangents
  - TODO: Um problema seria que a Landmark guarda apenas coordenadas de pontos e
    a gente precisa da tangente/orientacao. Temos que pensar como seria isso.

#### View: a struct with \label{sec:views}
  - image (caminho no disco)
  - id dos parametros intrinsecos da camera            (sometimes this is not set)
  - id dos parametros extrinsecos da camera (pose)     (sometimes this is not set)

#### Tracks

TracksBuilder
- Store tracks / correspondences
  - Track / Tracks in TracksBuilder
    - A Track 
      - is a corresponding image points with their imageId and FeatureId.
      - is called a submapTrack std::map<uint32_t, uint32_t>
    - Tracks
      - a collection of tracks is STLMAPTracks = std::map<uint32_t, submapTrack>;
    - Each element of Tracks is thus:
    (track_id, (image1_id, feature1_id), (image2_id, feature2_id),...)
    or
    {TrackIndex => {(imageIndex, featureIndex), ... ,(imageIndex, featureIndex)}
    - featureindex indexes into
      regions_per_image_.at(imageIndex).get()->Features()[featureIndex]
  - UnionFind tree
    - data structure to implement

STLMAPTracks
- tracks as a simpler map (each entry is a sequence of imageId and featureIndex):
- map_tracks[track_id].insert(feat.first);
  
    
- Build
  - Fuse pairwise corrrespondences
  - Input: PairwiseMatches
- Filter
  - Remove tracks that have conflict


#### PairwiseMatches

PairWiseMatches
- map<Pair, IndMatches>
  - pwm[image pair] = all maches between the pair of images
  
IndMatch
- i,j
- feature i matches feature j

IndMatches
- vector<IndMatch>
- all features matching between two images



### TODO 
  - search agan for all uses of the 5 pt algorithm and P3P and adapt for trifocal+p2pt
 

### How to add orientation to OpenMVG?

Pierre: 

The simplest would be to add here orientation information to store it once validated.
https://github.com/openMVG/openMVG/blob/develop/src/openMVG/sfm/sfm_landmark.hpp#L24

You would still need to modify the input feature point provider too by adding
perhaps a map of orientation here
https://github.com/openMVG/openMVG/blob/develop/src/openMVG/sfm/pipelines/sfm_features_provider.hpp#L33
Or we would always load ScaleInvariantPoint ()

The feature provider interface is the interface providing the data above the
robust estimation in the SfM pipeline.
The SfM_data container stores all the data once it is validated.
 
### Path for Trifocal

Adapted from Pierre email 29oct19

Here is the path I would follow for integration of new minimal solver.
Here I target towards contributing for the community and for large scale experiments.

1. Add the minimal solver in src/openMVG/multiview/XX.hpp XX.cpp or `src/openMVG/multiview/trifocal/*` 
2. Add corresponding unit test solver_essential_five_point_test.cpp#L257 
=> Make a PR here
3. Add a sample to demonstrate the solver to the community on real data and let people experiment with it
At this stage we will need to introduce the solver to the Robust Stage. OpenMVG
is using the Kernel concept to embed the solver, the metric and the point
selection to ACRansac in a generic way
https://github.com/openMVG/openMVG/tree/master/src/openMVG_Samples
4. Bootstrapping for Incremental SfM
OpenMVG is having two SfM pipeline
(the v1 close to my Accv2012 paper and v2 closer to what people do today to make
the pipeline a bit faster)

SequentialSfMReconstructionEngine:
As we have discussed yesterday, we could change the initial pair to a tuple, so
this way we could specify by hand the first triplet of camera to use for debugging.
https://github.com/openMVG/openMVG/blob/master/src/openMVG/sfm/pipelines/sequential/sequential_SfM.hpp#L103
Then we have to udpate AutomaticInitialPairChoice, ChooseInitialPair and MakeInitialPair3D

SequentialSfMReconstructionEngine2 already have the abstraction of the mechanism to bootstrap the reconstruction -> SfMSceneInitializer \label{pierre:comment}
Adding code in this engine could be simpler since you have to inherit and write you new SfMSceneInitializerTrifocal class

5. Bootstrap GlobalSfM
Since your solver is using less points than some other solver, we could perhaps have better results (more stable, and less outliers triplets)

