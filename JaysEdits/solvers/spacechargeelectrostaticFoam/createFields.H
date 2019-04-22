    Info<< "Reading physicalProperties\n" << endl;

    IOdictionary physicalProperties
    (
        IOobject
        (
            "physicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar epsilon0
    (
        physicalProperties.lookup("epsilon0")
    );

    dimensionedScalar k
    (
        physicalProperties.lookup("k")
    );


    Info<< "Reading field phi\n" << endl;
    volScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field rho\n" << endl;
    volScalarField rho1
    (
        IOobject
        (
            "rho1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
				dimensionedScalar("rho1",dimensionSet(0,-3,1,0,0,1,0),scalar(0))
    );

forAll(mesh.C(), celli)
{
	
  //vector Ci = mesh.C()[celli];
  // For a general case
  // S[celli] = yourFunction(Ci.x(), Ci.y(), Ci.z());

  // For your case
  //if(Ci.x() > x1 && Ci.x() < x2) S[celli] = S1;
  //else S[celli] = S2; 
	double rhovalue = 0.001259; //0.0016021766;	//1.6021766e-12;
	rho1[celli] = rhovalue;
	
}


/*
    Info<< "Calculating field rhoFlux\n" << endl;
    surfaceScalarField rhoFlux
    (
        IOobject
        (
            "rhoFlux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        -k*mesh.magSf()*fvc::snGrad(phi)
    );
*/
