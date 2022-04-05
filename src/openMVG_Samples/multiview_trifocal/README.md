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

Useful structures:

OpenMVG chama de Landmark a associacao de um ponto em 3D e suas observacoes em imagens:

https://github.com/openMVG/openMVG/blob/develop/src/openMVG/sfm/sfm_landmark.hpp

Chama de View uma 
struct com imagem (caminho no disco), id dos parametros intrinsecos da camera e id dos parametros extrinsecos da camera (pose)


Um problema seria que a Landmark guarda apenas coordenadas de pontos e a gente precisa da tangente/orientacao. Temos que pensar como seria isso.


TODO
  - search agan for all uses of the 5 pt algorithm and adapt for trifocal
 
