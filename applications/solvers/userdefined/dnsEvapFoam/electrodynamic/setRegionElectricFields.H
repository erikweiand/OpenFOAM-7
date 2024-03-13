// map first fluid field charge density onto electric field mesh
interpolObjEl.mapSrcToTgt(rhoQFluid, eqOp<scalar>(), rhoQEl);

volScalarField& rhoQ = rhoQEl;
volScalarField& phiQ = phiQEl;
volVectorField& ElectricFieldStrength = ElectricFieldStrengthEl;
