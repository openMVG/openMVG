original trifocal branch  created by pierre: develop_trifocal
working branch: trifocal_tmp
  - this working branch merges the develop branch periodically
  - when everything is working, we just make a nice commit from scratch to
    OpenMVG with Gabriel as author, with pierre as coauthors


Scrap notes


src/software 
 - Complete software(s) build on openMVG libraries


== Inclusion into OpenMVG ==

Main:

src/software/SfM/SfM_SequentialPipeline.py.in
  - calls all the executables from python.
  - main_SfM.cpp is our main file
    - ReconstructionEngine::Process() is the main function
        SfMEngine::INCREMENTALV:
          - we are mainly targeting this one
          
        SfMEngine::INCREMENTALV2:
          - this one allows to set the triples
  - also see main_GeometricFilter.cpp


Compute matching pairs -->  ComputeMatches               --> GeometricFilter -->
I: sfm_data.json            I: sfm_data.json, pairs.bin      I: sfm_data.json, matches.putative.bin
O: pairs.bin                O: matches.putative.bin          O: matches.f.bin

-->  main_SfM
     I: sfm_data.json, maatch_dir (matches.f.bin)
     O: reconstruction_sequential/



GeometricFilter 
- From the topolevel .py file, it uses a fundamental matrix filtercing model
- Why not use essential matrix when intrinsics available?
    - TODO: investigate where intrinsics are provided, if we can fix it, and
      whether intrinsics are actually optimized / estimated
        - sfm_data has intrinsics and views (see \ref{sec:views} below)


Useful structures:

Landmark 
  - a associacao de um ponto em 3D e suas observacoes em imagens: sfm_landmark.hpp

View: a struct with \label{sec:views}
  - imagem (caminho no disco)
  - id dos parametros intrinsecos da camera            (sometimes this is not set)
  - id dos parametros extrinsecos da camera (pose)     (sometimes this is not set)

Um problema seria que a Landmark guarda apenas coordenadas de pontos e a gente
precisa da tangente/orientacao. Temos que pensar como seria isso.


TODO
  - search agan for all uses of the 5 pt algorithm and P3P and adapt for trifocal+p2pt
 
