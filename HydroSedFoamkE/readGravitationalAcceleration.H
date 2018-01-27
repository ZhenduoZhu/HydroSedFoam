Info<< "\nReading gravitationalProperties" << endl;

IOdictionary gravitationalProperties
(
    IOobject
    (
        "gravitationalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const dimensionedVector g(gravitationalProperties.lookup("g"));
const Switch rotating(gravitationalProperties.lookup("rotating"));
const dimensionedVector Omega =
    rotating ? gravitationalProperties.lookup("Omega")
             : dimensionedVector("Omega", -dimTime, vector(0,0,0));
const dimensionedScalar magg = mag(g);
const dimensionedVector gHat = g/magg;
const Switch roughness(gravitationalProperties.lookup("roughness"));
const dimensionedScalar Chezy2 = gravitationalProperties.lookup("Chezy2");
const dimensionedScalar viscosity = gravitationalProperties.lookup("viscosity");
const dimensionedScalar mviscosity = gravitationalProperties.lookup("mviscosity");
const Switch sediment(gravitationalProperties.lookup("sediment"));
const Switch bedload(gravitationalProperties.lookup("bedload"));
const Switch suspend(gravitationalProperties.lookup("suspend"));
const dimensionedScalar D50 = gravitationalProperties.lookup("D50");
const dimensionedScalar taucb = gravitationalProperties.lookup("taucb");
const dimensionedScalar subgrav = gravitationalProperties.lookup("subgrav");
const dimensionedScalar qsconst = D50*sqrt(magg*subgrav*D50);
const dimensionedScalar Rep = qsconst/mviscosity;
const Switch geomorph(gravitationalProperties.lookup("geomorph"));
const dimensionedScalar porobed = gravitationalProperties.lookup("porobed");
const dimensionedScalar accelbed = gravitationalProperties.lookup("accelbed");
const dimensionedScalar wdeposit = gravitationalProperties.lookup("wdeposit");
const Switch kepsilonmodel(gravitationalProperties.lookup("kepsilonmodel"));
const Switch zerowallcell(gravitationalProperties.lookup("zerowallcell"));
