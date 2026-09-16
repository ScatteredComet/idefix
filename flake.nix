{
  description = "A basic flake with a shell";
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  inputs.systems.url = "github:nix-systems/default";
  inputs.flake-utils = {
    url = "github:numtide/flake-utils";
    inputs.systems.follows = "systems";
  };

  outputs =
    {
      nixpkgs,
      flake-utils,
      self,
      ...
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = nixpkgs.legacyPackages.${system};

        enableMHD = false;
        enableMPI = false;
        highOrderFargo = false;
        debug = false;
        runtimeChecks = false;
        warningsAsErrors = false;
        customEquationOfState = false;
        hdf5 = false;
        evolveVectorPotential = false;

        cppApp = pkgs.stdenv.mkDerivation {
          pname = "idefix";
          version = "";
          src = self;

          nativeBuildInputs = with pkgs; [
            cmake
            ninja
          ];

          cmakeFlags = [
            "Idefix_MHD=${if enableMHD then "ON" else "OFF"}"
            "Idefix_MPI=${if enableMPI then "ON" else "OFF"}"
            "Idefix_HIGH_ORDER_FARGO=${if highOrderFargo then "ON" else "OFF"}"
            "Idefix_DEBUG=${if debug then "ON" else "OFF"}"
            "Idefix_RUNTIME_CHECKS=${if runtimeChecks then "ON" else "OFF"}"
            "Idefix_WERROR=${if warningsAsErrors then "ON" else "OFF"}"
            "Idefix_CUSTOM_EOS=${if customEquationOfState then "ON" else "OFF"}"
            "Idefix_HDF5=${if hdf5 then "ON" else "OFF"}"
            "Idefix_EVOLVE_VECTOR_POTENTIAL=${if evolveVectorPotential then "ON" else "OFF"}"
          ];

          buildInputs =
            (if enableMPI then [ pkgs.openmpi ] else [ ])
            ++ (
              if hdf5 && enableMPI then # use parallelised version of hdf5
                [ pkgs.hdf5-mpi ]
              else if hdf5 then
                [ pkgs.hdf5 ]
              else
                [ ]
            );
        };
      in
      {
        packages.default = cppApp;

        apps.default = {
          type = "app";
          program = "${cppApp}/bin/idefix";
        };

        devShells.default = pkgs.mkShell {
          inputsFrom = [ cppApp ];

          packages = with pkgs; [
            clang-tools
            gdb
            bear
          ];
        };
      }
    );
}
