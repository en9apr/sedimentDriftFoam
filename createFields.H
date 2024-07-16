/////////////////////////
// APR modified T to C //
/////////////////////////

Info<< "Reading field C\n" << endl;

volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

/////////////////////////
// APR modified T to C //
/////////////////////////



Info<< "Reading field U\n" << endl;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

////////////
// APR added: needed only for writing!
////////////
Info<< "Reading field p\n" << endl;

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
////////////
// APR added needed only for writing!
////////////


// APR removed transportProperties dictionary


///////////////
// APR added //
///////////////
volScalarField sImplicit
(
	IOobject
	(
		"sImplicit",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);

volScalarField sExplicit
(
	IOobject
	(
		"sExplicit",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);

IOdictionary sedimentProperties
(
    IOobject
    (
		"sedimentProperties",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
    )
);

dimensionedVector fallVelocity("fallVelocity", dimVelocity, sedimentProperties);

volVectorField Vs
(
    IOobject
    (
        "Vs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("Vs", U.dimensions(), vector(0.0, 0.0, 0.0))
);

Vs = fallVelocity;

dimensionedVector up("up", dimless, vector(0.0, 0.0, 0.0) );

up = fallVelocity / (vector(1.0, 1.0, 1.0) & fallVelocity);

dimensionedScalar fallVelocityScalar = fallVelocity & up.value();

// read parameters from sedimentProperties file
dimensionedScalar particleSize("particleSize", dimLength, sedimentProperties);
dimensionedScalar shields("shields", dimensionSet(0,2,-2,0,0,0,0), sedimentProperties);
dimensionedScalar roughness("roughness", dimLength, sedimentProperties);
dimensionedScalar schmidt("schmidt", dimless, sedimentProperties);
dimensionedScalar concentrationRelativeTolerance("concentrationRelativeTolerance", dimless, sedimentProperties);

dimensionedScalar COld("COld",dimless,0.0);
dimensionedScalar CNew("CNew",dimless,0.0);

const Switch implicitOn = sedimentProperties.lookupOrDefault("implicitOn", true);
const Switch explicitOn = sedimentProperties.lookupOrDefault("explicitOn", true);

// location of free surface
const word freeSurfaceName=sedimentProperties.lookup("freeSurfaceName");

label surfPatchID = mesh.boundaryMesh().findPatchID(freeSurfaceName);

if(surfPatchID < 0) 
{
	Info<< "Have not found any 'freeSurface' patches in the grid" << endl;  // error message
}

// location of bed
const word bedName=sedimentProperties.lookup("bedName");

label bedPatchID = mesh.boundaryMesh().findPatchID(bedName);

if(bedPatchID < 0)
{
	Info<< "Have not found any 'bed' patches in the grid" << endl;  // error message
}
///////////////
// APR added //
///////////////


#include "createPhi.H"


///////////////
// APR added //
///////////////
Info<< "Reading/calculating face flux field phiC \n" << endl;

surfaceScalarField phiC
(
    IOobject
    (
        "phiC",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::flux(U + Vs)
);
///////////////
// APR added //
///////////////



///////////////
// APR added // from simpleFoam.C
///////////////
singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);
///////////////
// APR added // from simpleFoam.C
///////////////



#include "createFvOptions.H"
