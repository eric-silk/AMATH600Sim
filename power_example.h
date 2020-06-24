#ifndef POWER_SIM_H
#define POWER_SIM_H

#include "petscsnes.h"
#include <petscdmnetwork.h>

constexpr size_t NGEN_AT_BUS_MAX = 15
constexpr size_t NLOAD_AT_BUS_MAX = 1;

struct _p_UserCtx_Power
{
  _p_UserCtx_Power() = default;
  PetscScalar Sbase;
  PetscBool jac_error; // Introduce error in the jacobian...for testing I imagine?
  PetscInt compkey_branch;
  PetscInt compkey_bus;
  PetscInt compkey_gen;
  PetscInt compkey_load;
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

#endif // POWER_SIM_H
