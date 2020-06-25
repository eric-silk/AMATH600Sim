#ifndef POWER_SIM_H
#define POWER_SIM_H

#include "petscsnes.h"
#include <petscdmnetwork.h>

constexpr size_t NGEN_AT_BUS_MAX = 15;
constexpr size_t NLOAD_AT_BUS_MAX = 1;

// I may take the time to do a better job of considering encapsulation. But,
// for now, it appears these really are just data structures (not fully fledged
// objects)

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
  PetscInt owner;
  PetscInt internal_i;
  PetscInt ngen;
  PetscInt gidx[NGEN_AT_BUS_MAX];
  PetscInt nload;
  PetscInt lidx[NLOAD_AT_BUS_MAX];
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

struct Load
{
  Load() = default;
  PetscInt bus_i;
  char i[20];
  char id[20];
  PetscInt status;
  PetscInt area;
  PetscInt zone;
  PetscScalar pl;
  PetscScalar ql;
  PetscScalar il;
  PetscScalar iq;
  PetscScalar yp;
  PetscScalar yq;
  PetscScalar scale_load;
  PetscInt owner;
  PetscInt internal_i;
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

struct Gen
{
  Gen() = default;
  PetscInt bus_i;
  char i[20];
  char id[20];
  PetscScalar pg;
  PetscScalar qg;
  PetscScalar qt;
  PetscScalar qb;
  PetscScalar vs;
  PetscInt ireg;
  PetscScalar mbase;
  PetscScalar zr;
  PetscScalar zx;
  PetscScalar rt;
  PetscScalar xt;
  PetscScalar gtap;
  PetscInt status;
  PetscScalar rmpct;
  PetscScalar pt;
  PetscScalar pb;
  PetscInt o1;
  PetscScalar f1;
  PetscScalar scale_gen;
  PetscInt internal_i;
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

struct EdgePower
{
  EdgePower() = default;
  PetscInt      fbus;                                                                                  
  PetscInt      tbus;                                                                                  
  char      i[20]; /* Bus Number or extended bus name*/                                                
  char      j[20]; /* Bus Number or extended bus name*/                                                
  char      ckt[20]; /* Circuit identifier. 1 by default */                                            
  PetscScalar   r; /* Branch resistance: pu */                                                         
  PetscScalar   x; /* Branch reactance: pu */                                                          
  PetscScalar   b; /* Branch charging susceptance: pu */                                               
  PetscScalar   rateA; /* rate A in MVA */                                                             
  PetscScalar   rateB; /* rate B in MVA */                                                             
  PetscScalar   rateC; /* rate C in MVA */                                                             
  PetscScalar   tapratio;                                                                              
  PetscScalar   phaseshift;                                                                            
  PetscScalar   gi; /* Complex admittance at 'i' end: pu */                                            
  PetscScalar   bi; /* Complex admittance at 'i' end: pu */                                            
  PetscScalar   gj; /* Complex admittance at 'j' end: pu */                                            
  PetscScalar   bj; /* Complex admittance at 'j' end: pu */                                            
  PetscInt  status; /* Service status */                                                               
  PetscScalar   length; /* Line length */                                                              
  PetscInt  o1; /* Owner number */                                                                     
  PetscScalar   f1; /* Fraction of ownership */                                                        
  PetscScalar   yff[2],yft[2],ytf[2],ytt[2]; /* [G,B] */                                               
  PetscInt  internal_i; /* Internal From Bus Number */                                                 
  PetscInt  internal_j; /* Internal To Bus Number */                                                   
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

struct PFData
{
    PFData() = default;
    PetscScalar sbase;
    PetscInt nbus, ngen, nbranch, nload;

    VertexPower bus;
    Load load;
    Gen gen;
    EdgePower branch;
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

extern PetscErrorCode PFReadMatPowerData(PFData*, char*);
extern PetscErrorCode GetListofEdges_Power(PFData*, PetscInt*);
extern PetscErrorCode FormJacobian_Power(SNES,Vec, Mat, Mat, void*);
extern PetscErrorCode FormJacobian_Power_private(DM, Vec, Mat, PetscInt, PetscInt, const PetscInt*, const PetscInt*, void*);
extern PetscErrorCode FormFunction_Power(DM, Vec, Vec, PetscInt, PetscInt, const PetscInt*, const PetscInt*, void*);
extern PetscErrorCode SetInitialGuess_Power(DM, Vec, PetscInt, PetscInt, const PetscInt *, const PetscInt *, void*);

#endif // POWER_SIM_H
