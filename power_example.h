#ifndef POWER_SIM_H
#define POWER_SIM_H

#include "petscsnes.h"
#include <petscdmnetwork.h>

constexpr size_t NGEN_AT_BUS_MAX = 15;
constexpr size_t NLOAD_AT_BUS_MAX = 1;

struct UserCtx_Power
{
  UserCtx_Power() = default;
  PetscScalar Sbase;
  PetscBool jac_error; // Introduce error in the jacobian...for testing I imagine?
  PetscInt compkey_branch;
  PetscInt compkey_bus;
  PetscInt compkey_gen;
  PetscInt compkey_load;
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

// Bus data, 11 columns
// TODO Change arrays to std::arrays
struct VertexPower
{
  VertexPower() = default;
  PetscInt bus_i;
  char i[20];
  char name[20];
  PetscScalar baskV;
  PetscInt ide; // Bus type code
  PetscScalar gl;  // Active shunt admittance
  PetscScalar bl;  // reactive shunt admittance
  PetscInt area; // area number
  PetscInt zone; // zone number
  PetscScalar vm; // bus voltage mag in PU
  PetscScalar va; // v phase angle
  PetsInt owner;
  PetsInt internal_i;
  PetsInt ngen;
  PetsInt gidx[NGEN_AT_BUS_MAX];
  PetsInt nload;
  PetsInt lidx[NLOAD_AT_BUS_MAX];
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));



#endif // POWER_SIM_H
