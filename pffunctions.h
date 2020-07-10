#ifndef MY_PFFUNCTIONS_H
#define MY_PFFUNCTIONS_H

#include "power.h"

// by "get" they mean just print (if its enabled)
PetscErrorCode GetListOfEdges_Power(PFDATA*, PetscInt*);
PetscErrorCode FormJacobian_Power(SNES, Vec, Mat, Mat, void*);
PetscErrorCode FormJacobian_Power_private(DM, Vec, Mat, PetscInt, PetscInt, const PetscInt*, const PetscInt*, void*);
PetscErrorCode FormFunction_Power(DM, Vec, Vec, PetscInt, PetscInt, const PetscInt*, const PetscInt*, void*);
PetscErrorCode SetInitialGuess_Power(DM,Vec,PetscInt,PetscInt,const PetscInt *,const PetscInt *,void*);

#endif//MY_PFFUNCTIONS_H
