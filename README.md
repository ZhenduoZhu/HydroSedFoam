# HydroSedFoam
A parallelized two-dimensional hydrodynamic, sediment transport, and bed morphology model based on OpenFOAM.

Author: Zhenduo Zhu (zhenduoz@buffalo.edu)

Department of Civil, Structural and Environmental Engineering, University at Buffalo




## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

OpenFOAM 2.3.0 need to be installed first. This model was developed and tested with OpenFOAM 2.3.0. Please note that other versions of OpenFOAM may not have the same organization of folders and files so you may not find the folders metioned below with another version.

Information on OpenFOAM 2.3.0 can be found here: https://openfoam.org/release/2-3-0/

More information on how to install OpenFOAM can be found here: https://openfoam.com/ and https://openfoam.org/

Official OpenFOAM Repository is at https://github.com/OpenFOAM

### Compilation

1. Compile 2D k-epsilon turbulent model

 * Copy the folder, kEpsilon2D, into `/src/TurbulenceModels/turbulenceModels/incompressible/RAS`
 * Replace files in `/src/TurbulenceModels/turbulenceModels/incompressible/RAS/Make` with the ones in the kEpsilon2D folder
 * Use wmake to compile it.

2. Compile HydroSedFoam

 * Copy the folder, HydroSedFoamkE into `/applications/solvers/incompressible`
 * Use wmake to compile it.

If you are using the OpenFOAM docker image, the above directories can be found in `/opt/OpenFOAM/OpenFOAM-v1712/`. You may need to login to the docker image as root to compile. To do so, edit `startOpenFoam` to read:

    docker exec -u 0 -it of_v1712 /bin/bash -rcfile /opt/OpenFOAM/setImage_v1712.sh
    
instead of

    docker exec -it of_v1712 /bin/bash -rcfile /opt/OpenFOAM/setImage_v1712.sh


### Running the Program


Case studies are provided in the folder, Cases. There are three cases:

 * Case 1: Meandering Channel Laboratory Experiments
 * Case 2: Sediment Transport and Bed Morphology
 * Case 3: Kalamazoo River, Michigan, USA

More information can be found in the reference below.



## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.


## Reference

**Title of Manuscript**:
HydroSedFoam: A new parallelized two-dimensional hydrodynamic, sediment transport, and bed morphology model

**Authors**: Zhenduo Zhu, Jessica Z. LeRoy, Bruce L. Rhoads, and Marcelo H. Garcia

**Corresponding Author**: Zhenduo Zhu (zhenduoz@buffalo.edu)

**Journal**: Computers & Geosciences

**DOI Number of Manuscript**: 

**Code Repositories**:
 * [Author's GitHub Repository](https://github.com/ZhenduoZhu/HydroSedFoam)

**Abstract**: Depth-averaged two-dimensional (2D) models are useful tools for understanding river morphodynamics through the computation of hydrodynamics, sediment transport, and an evolving river bed morphology. This paper presents a new parallelized 2D hydrodynamic, sediment transport, and bed morphology model, HydroSedFoam. The model uses the Message Passing Interface (MPI) for code parallelization and adopts a depth-averaged k-epsilon turbulence model. Three different case studies, including a laboratory experiment, an analytical solution, and a field-scale river reach, show good agreement with HydroSedFoam simulations. Further development and modification of the model are relatively straightforward to accomplish within the OpenFOAM framework.
