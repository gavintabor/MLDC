    dimensionedScalar boxSide
    (
	"boxSide",
	dimLength,
	mesh().lookupObject<IOdictionary>("transportProperties")
    );

    dimensionedScalar boxVel
    (
	"boxVel",
	dimVelocity,
	mesh().lookupObject<IOdictionary>("transportProperties")
    );

    dimensionedScalar nu
    (
        "nu",
        dimViscosity,
        mesh().lookupObject<IOdictionary>("transportProperties")
    );

    scalar Re = (boxSide*boxVel/nu).value();

    Info << "Re = " << Re << endl;

    volScalarField x = mesh().C().component(0)/boxSide;
    volScalarField y = mesh().C().component(1)/boxSide;

    volScalarField f = pow(x,4.0)-2.0*pow(x,3.0)+pow(x,2.0);
    volScalarField g = pow(y,4.0) - pow(y,2.0);

    volScalarField F=0.2*pow(x,5.0)-0.5*pow(x,4.0)+pow(x,3.0)/3.0;
    volScalarField F1=-4.0*pow(x,6.0)+12.0*pow(x,5.0)-14.0*pow(x,4.0)
	                   +8.0*pow(x,3.0)-2.0*x*x;


    volScalarField F2=0.5*f*f;
    volScalarField G1=-24.0*pow(y,5.0)+8.0*pow(y,3.0)-4.0*y;

    volScalarField fp=4.0*pow(x,3.0)-6.0*x*x + 2.0*x;
    volScalarField fppp = 24.0*x-12.0;
    volScalarField gp=4.0*y*y*y-2.0*y;
    volScalarField gpp=12.0*y*y-2.0;
    volScalarField gppp=24.0*y;

    volScalarField b_y
    (
        IOobject
        (
            "b_y",
	    mesh().time().timeName(),
	    mesh(),
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
        ),
	((-8.0/Re)*(24.0*F+2.0*fp*gpp+fppp*g)-64.0*(F2*G1-g*gp*F1))
	*boxVel*boxVel/boxSide
    );

// ************************************************************************* //
