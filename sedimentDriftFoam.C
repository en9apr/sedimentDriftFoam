/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    scalarTransportFoam

Description
    Solves the steady or transient transport equation for a passive scalar.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

///////////////
// APR added // To access nEff in zeroSedimentationFlux BC: from simpleFoam.C
///////////////
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
//#include "kEpsilon.H"
///////////////
// APR added // Also options from simpleFoam.C
///////////////

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

		///////////////
		// APR added //
		///////////////
		
        // a sink on RHS for negative fall velocities
		if(implicitOn)
		{
		    // surface boundary condition for particles
			const fvPatch& patchSurf = mesh.boundary()[surfPatchID];
			// surface: implicit term: 
			forAll(patchSurf,surf) 
			{
				label cellnumber = patchSurf.faceCells()[surf];
				sImplicit[cellnumber] = -fallVelocityScalar.value() * patchSurf.magSf()[surf] / mesh.V()[cellnumber];  // units 1/sec				
			}
			Info<< "Implicit source term applied" << nl << endl;
		}
		

        
        // a source on RHS for negative fall velocities
		if(explicitOn)
		{
			const fvPatch& patchBed = mesh.boundary()[bedPatchID];
		    // loop over all bed cells
			forAll(patchBed,bed)
			{
				label cellbed = patchBed.faceCells()[bed];
				scalar dist_to_wall = 0.5 * mesh.V()[cellbed] / patchBed.magSf()[bed]; // in lack of another working option 
				scalar bedvelocity = mag(U[cellbed]);  // no compile errors
				scalar bedshear = 0.4 * bedvelocity / Foam::log(30.0 * dist_to_wall / roughness.value());  // log-law
				bedshear = 1000.0 * bedshear * bedshear;  // shearvelocity to bed shear stress
				scalar critshear = shields.value() * 1650.0 * 9.81 * particleSize.value(); // Shields
				scalar Tstar = (bedshear - critshear) / critshear;  // parameter in van Rijns formula
				if(Tstar < 0.0) Tstar = 0.0; // in case of negative Tstar
				scalar Dstar = particleSize.value() * 25296.0;  // parameter in van Rijns formula
				scalar vanRijn = 0.015 * particleSize.value() / dist_to_wall * Foam::pow(Tstar, 1.5) / Foam::pow(Dstar,0.3);  // concentration according to van Rijn
		
		        // pick-up rate for concrete bed, where sediment erosion is limited
				if(bedshear > critshear) 
				{
					// settling particles: 
					if(vanRijn > C[cellbed]) 
					{
						sExplicit[cellbed] = -fallVelocityScalar.value() / (2.0 * dist_to_wall) * C[cellbed];  // 1/second: 
					} 
					else 
					{
						sExplicit[cellbed] = -fallVelocityScalar.value() / (2.0 * dist_to_wall) * vanRijn;  // 1/second: 
					}
				}
			} // end of the bed
			Info<< "Explicit source term applied" << nl << endl;
		}
		///////////////
		// APR added //
		///////////////		
		

    	// solve the equation
        while (simple.correctNonOrthogonal())
        {
            // APR removed TEqn

			///////////////
			// APR added //
			///////////////
			COld = C.weightedAverage(mesh.V());

            fvScalarMatrix CEqn
            (
              fvm::ddt(C)
              + fvm::div(phiC, C)
              + fvm::Sp(sImplicit,C)
              ==
              sExplicit
			  + fvm::laplacian(turbulence->nuEff()/schmidt, C)
            );
			///////////////
			// APR added //
			///////////////


            ////////////////////////
            // APR changed T to C //
            ////////////////////////
            CEqn.relax();
            fvOptions.constrain(CEqn);
            CEqn.solve();
            fvOptions.correct(C);
            ////////////////////////
            // APR changed T to C //
            ////////////////////////


			///////////////
			// APR added //
			///////////////
			CNew = C.weightedAverage(mesh.V()); 
			Info<< "COld = " << COld.value() << endl;
			Info<< "CNew = " << CNew.value() << endl;
			Info<< "Ratio = " << (mag(CNew.value() - COld.value())/mag(CNew.value())) << endl;
			///////////////
			// APR added //
			///////////////

        }

		///////////////
		// APR added //
		///////////////
		if((mag(CNew - COld)/mag(CNew)) < concentrationRelativeTolerance)
		{
	    	U.write();
	    	p.write();
			Vs.write();
			sExplicit.write();
			sImplicit.write();
			phi.write();
			C.write();
            const volScalarField& k = mesh.lookupObject<volScalarField>("k");
            k.write();
            const volScalarField& nut = mesh.lookupObject<volScalarField>("nut");
            nut.write();
            const volScalarField& epsilon = mesh.lookupObject<volScalarField>("epsilon");
            epsilon.write();            
			Info<< "Stopping because concentration tolerance reached" << endl;
        	Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
			break;
		}
		///////////////
		// APR added //
		///////////////

        //runTime.write();

		///////////////
		// APR added //
		///////////////
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
		///////////////
		// APR added //
		///////////////

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
