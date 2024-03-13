{
    Info<< "k(" << runTime.timeName() << ") = "
        << 1.0/2.0*average(magSqr(U)).value() << endl;

    Info<< "epsilon(" << runTime.timeName() << ") = "
        << (0.5*nuGlobal*average(magSqr(twoSymm(fvc::grad(U))))).value() << endl;

    Info<< "epsilon_vreman(" << runTime.timeName() << ") = "
	<< (nuGlobal*average(magSqr(fvc::grad(U)))).value() << endl;
}
