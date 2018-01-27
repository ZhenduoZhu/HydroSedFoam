Info<< "Reading field h\n" << endl;
volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field h0 if present\n" << endl;
volScalarField h0
(
    IOobject
    (
        "h0",
        runTime.timeName(),
		mesh,
        IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
    ),
    mesh
);
//if(bathymetry)
//{
////	Info<< "Interpolating bathymetry" << endl;
////	const int MAXN = 10000;
////	double data[2][MAXN];
////	freopen("bathymetry.txt","rb",stdin);
////	for (int counti=0;counti<MAXN;counti++) for (int countj=0;countj<3;countj++) scanf("%d",&data[countj][counti]);
//	IOdictionary bathymetry
//    (
//        IOobject
//        (
//            "bathymetry",
//            runTime.constant(),
//            mesh,
//            IOobject::READ_IF_PRESENT,
//            IOobject::NO_WRITE
//        )
//    );
//	vectorList xyz(bathymetry.lookup("xyz"));
//	//dictionary interpolationDict = mesh.solutionDict().subDict("interpolationSchemes");
//    //autoPtr<interpolation<vector> > Uinterp = interpolation<vector>::New(interpolationDict, U);
//    label cellI = mesh.findCell(xyz);
//	h0 = xyz.z();
//    //vector Ui = Uinterp->interpolate(xyz, cellI);
//}

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

Info<< "Creating field hU\n" << endl;
volVectorField hU
(
    IOobject
    (
        "hU",
        runTime.timeName(),
        mesh
    ),
    h*U,
    U.boundaryField().types()
);

Info<< "Creating field hTotal for post processing\n" << endl;
volScalarField hTotal
(
    IOobject
    (
        "hTotal",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    h+h0
);

volScalarField nut
(
    IOobject
    (
        "nut",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);
volScalarField zerotemp
(
    IOobject
    (
        "zerotemp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
	dimensionedScalar("zerotemp", dimless, 0.0)
);
volScalarField suschange
(
    IOobject
    (
        "suschange",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
	dimensionedScalar("suschange", dimless, 0.0)
);
volVectorField qsbed
(
    IOobject
    (
        "qsbed",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField qstar
(
    IOobject
    (
        "qstar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
	dimensionedScalar("qstar", dimless, 0.0)	//Z: has to initialize?
);
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
hTotal.write();

#   include "createPhi.H"

Info<< "Creating Coriolis Force" << endl;
const dimensionedVector F("F", ((2.0*Omega) & gHat)*gHat);

//Z: for k-epsilon model
singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);
