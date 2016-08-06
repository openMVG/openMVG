let
     pkgs = import <nixpkgs> {};
     stdenv = pkgs.stdenv;
in rec {
          cudaEnv = stdenv.mkDerivation rec {
          name = "cuda-env";
          version = "1.1.1.1";
          src = ./.;
          buildInputs = [ pkgs.openblas pkgs.opencv3 pkgs.cudnn5_cudatoolkit80 pkgs.torchPackages.torch pkgs.torchPackages.nn pkgs.cudatoolkit pkgs.cmake pkgs.gnumake pkgs.pkgconfig pkgs.stdenv ];
     };
}
