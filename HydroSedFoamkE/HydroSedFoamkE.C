/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    HydroSedFoamkE

Description
    Transient solver for 2D shallow-water equations with rotation, bed shear stress.

    Developed by Zhenduo Zhu at UIUC, 2015/5/28

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
//#include "turbulenceModel.H"
#include "RASModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"
	#include "readGravitationalAcceleration.H"
	#include "createFields.H"

	pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    //Z: while (runTime.loop())
	while (runTime.run())	//Z
    {	
		#include "readTimeControls.H"	//Z: used to read the control parameters used by setDeltaT
		#include "CourantNo.H"
		#include "setDeltaT.H"	//Z
		
		runTime++;	//Z
		
		Info<< "\n Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            nut = viscosity;
			
			if(kepsilonmodel)	//Z: k-epsilon turbulence model
			{
				Info<< "Running k-epsilon model\n" << endl;
			
				turbulence->correct();
				nut = turbulence->nut() + viscosity;
			}
			
			surfaceScalarField phiv("phiv", phi/fvc::interpolate(h));

            fvVectorMatrix hUEqn
            (
                fvm::ddt(hU)
              + fvm::div(phiv, hU)
			  - fvm::laplacian(nut, hU)
            );

            hUEqn.relax();

            if (pimple.momentumPredictor())
            {
                if (rotating)
                {
					if (roughness)	//Z: ADD BED SHEAR STRESS TERM
					{
						solve(hUEqn - fvc::laplacian(nut*h, U) + fvc::laplacian(nut, hU) + (F ^ hU) + (magg/Chezy2)*mag(U)*U == -magg*h*fvc::grad(h + h0));
					}
					else
					{
						solve(hUEqn - fvc::laplacian(nut*h, U) + fvc::laplacian(nut, hU) + (F ^ hU) == -magg*h*fvc::grad(h + h0));
					}	
                }
                else
                {
                    if(roughness)
					{
						solve(hUEqn - fvc::laplacian(nut*h, U) + fvc::laplacian(nut, hU) + (magg/Chezy2)*mag(U)*U == -magg*h*fvc::grad(h + h0));
					}
					else
					{
						solve(hUEqn - fvc::laplacian(nut*h, U) + fvc::laplacian(nut, hU) == -magg*h*fvc::grad(h + h0));
					}
                }

                // Constrain the momentum to be in the geometry if 3D geometry
                if (mesh.nGeometricD() == 3)
                {
                    hU -= (gHat & hU)*gHat;
                    hU.correctBoundaryConditions();
                }
            }

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                volScalarField rAU(1.0/hUEqn.A());
                surfaceScalarField ghrAUf(magg*fvc::interpolate(h*rAU));

                surfaceScalarField phih0(ghrAUf*mesh.magSf()*fvc::snGrad(h0));

                volVectorField HbyA("HbyA", hU);
                if (rotating)
                {
					if (roughness)
					{
						HbyA = rAU*(hUEqn.H() + fvc::laplacian(nut*h, U) - fvc::laplacian(nut, hU) - (F ^ hU) - (magg/Chezy2)*mag(U)*U);
					}
					else
					{
						HbyA = rAU*(hUEqn.H() + fvc::laplacian(nut*h, U) - fvc::laplacian(nut, hU) - (F ^ hU));
					}
                }
                else
                {
					if (roughness)
					{
						HbyA = rAU*(hUEqn.H() + fvc::laplacian(nut*h, U) - fvc::laplacian(nut, hU) - (magg/Chezy2)*mag(U)*U);
					}
					else
					{
						HbyA = rAU*(hUEqn.H() + fvc::laplacian(nut*h, U) - fvc::laplacian(nut, hU));
					}
                }

                surfaceScalarField phiHbyA
                (
                    "phiHbyA",
                    (fvc::interpolate(HbyA) & mesh.Sf())
                  + fvc::interpolate(rAU)*fvc::ddtCorr(h, hU, phi)
                  - phih0
                );

                while (pimple.correctNonOrthogonal())
                {
                    fvScalarMatrix hEqn
                    (
                        fvm::ddt(h)
                      + fvc::div(phiHbyA)
                      - fvm::laplacian(ghrAUf, h)
                    );

                    hEqn.solve(mesh.solver(h.select(pimple.finalInnerIter())));

                    if (pimple.finalNonOrthogonalIter())
                    {
                        phi = phiHbyA + hEqn.flux();
                    }
                }

                hU = HbyA - rAU*h*magg*fvc::grad(h + h0);

                // Constrain the momentum to be in the geometry if 3D geometry
                if (mesh.nGeometricD() == 3)
                {
                    hU -= (gHat & hU)*gHat;
                }

                hU.correctBoundaryConditions();
            }
        }
		
		dimensionedScalar hMin("hMin", dimLength, 0.01);
		bound(h, hMin);

        U == hU/h;
        hTotal == h + h0;

		//Z: SEDIMENT TRANSPORT AND GEOMORPHOLOGY
		if(sediment)
		{
			//volScalarField qstar(magSqr(U)/Chezy2/subgrav/D50 - taucb);
			//volVectorField qsbed(qsconst*8.0*pow(qstar,1.5)*U/mag(U));
			
			if(bedload)	//Meyer-Peter and Muller formula
			{
				Info<< "Calculating bedload sediment transport\n" << endl;
				qstar = magSqr(U)/Chezy2/subgrav/D50 - taucb;
				
				bound(qstar, 0);	//Z: /src/finiteVolume/cfdTools/general/bound/bound.C;
				qsbed = qsconst*8.0*pow(qstar,1.5)*U/mag(U);
				//Info<< qsbed << endl;
				if(h0 < zerotemp)
				{
					qsbed = qsbed * 0;
				}
//				if((qstar > zerotemp) && (h0 > zerotemp))	//Z: zero is the bed rock elevation (no more sediment available under zero)
//				{
//					Info<< "qstat>zero" << endl;
//					qsbed = qsconst*8.0*pow(qstar,1.5)*U/mag(U);
//					Info<< qsbed << endl;
//				}
//				else
//				{
//					qsbed = qsbed*0;
//				}
			}
			if(suspend)	//
			{
				Info<< "Calculating suspended sediment transport\n" << endl;
				surfaceScalarField phivC("phivC", phi/fvc::interpolate(h));
				volScalarField entrain
				(
					"entrain",
					1.0/(1.0/0.3+1.0/(1.3E-7*pow(0.708*(magg/Chezy2)*mag(U)/wdeposit*pow(Rep,0.6),5)))
				);
				if(h0 < zerotemp)
				{
					entrain = entrain*0;
				}
				volScalarField deposit
				(
					"deposit",
					C*2.0
				);
				suschange = deposit - entrain;

				solve
				(
					fvm::ddt(C)
					+ fvm::div(phivC, C)
					- fvm::laplacian(nut, C)
					==
					-suschange*wdeposit/h
				);

				bound(C, 0);
//				if(C < zerotemp)
//				{
//					C = C *0;
//					Info<< "Negative C\n" << endl;
//				}
			}
			if(geomorph) //bed change
			{
				Info<< "Calculating bed change\n" << endl;
				
//				if(qstar > zerotemp)
//				{
					h0 = h0 - fvc::div(qsbed) / (1.0-porobed) * runTime.deltaT() * accelbed + suschange*wdeposit*runTime.deltaT()*accelbed;
//					solve
//					(
//						fvm::ddt(h0)
//						+ fvc::div(qsbed)/ (1.0-porobed) * accelbed
//						==
//						suschange*wdeposit * accelbed
//					);
//				}
//				else
//				{
//					h0 = h0 + suschange*wdeposit*runTime.deltaT()*accelbed;
//					solve
//					(
//						fvm::ddt(h0)
//						==
//						suschange*wdeposit * accelbed
//					);
//				}
				
				if(h0 < zerotemp)
				{
					h0 = h0*0;
				}
				hTotal == h + h0;
			}
		}

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
