// map electric field properties onto fluid field mesh
interpolObjFluid.mapSrcToTgt(phiQEl, eqOp<scalar>(), phiQFluid);
interpolObjFluid.mapSrcToTgt(ElectricFieldStrengthEl, eqOp<vector>(), ElectricFieldStrengthFluid);
