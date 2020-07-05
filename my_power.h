#ifndef MY_POWER_H
#define MY_POWER_H

#include <petscsnes.h>
#include <petscdmnetwork.h>

#define MAXLINE 1000
#define REF_BUS 3
#define PV_BUS 2
#define PQ_BUS 1
#define ISOLATED_BUS 4
#define NGEN_AT_BUS_MAX 15
#define NLOAD_AT_BUS_MAX 1

struct _p_UserCtx_Power
{
  // Base power in VA (S^2 = P^2 + Q^2)
  PetscScalar Sbase;
  // Do we want to introduce error into the jacobian?
  PetscBool jac_error;
  // These are the values used to indicate the type interally to DM
  // assigned by the output of DMRegisterNetworkComponent
  PetscInt compkey_branch;
  PetscInt compkey_bus;
  PetscInt compkey_gen;
  PetscInt compkey_load;
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar)); // align for arraylike access

// Why is this the only one that has the _p_ prefix
// (which I assume is "ptr") that doens't typedef a ptr?
typedef struct _p_UserCtx_Power UserCtx_Power;

// Struct representing a bus
struct _p_VERTEX_Power
{
  PetscInt bus_i; // integer bus #, used by Matpower
  char i[20]; // Bus num as a string (why a string??)
  char name[20]; // duh, but if the num can be a string, why this?
  PetscScalar basekV; // For conversion to PU, I'm sure
  PetscInt ide; // the bus type, see the macros above
  PetscScalar gl; // active shunt admittance to gnd
  PetscScalar bl; // reactive shunt admittance to gnd
  PetscInt area; // "area" number...for??
  PetscInt zone; // "zone" number...for??
  PetscScalar vm; // Voltage magnitude. in PU!
  PetscScalar va; // Voltage phase magnitude. in PU!
  PetscInt owner; // "owner" of the bus. Who can an owner be?
  PetscInt internal_i; // "internal" bus num. What's different about this than the ID?
  PetscInt ngen; // num gens at this bus (should be capped by NGEN_AT_BUS_MAX
  PetscInt gidx[NGEN_AT_BUS_MAX]; // list of the generator indices at this bus, for use with the GEN struct
  PetscInt nload;  // number of loads, duh.
  PetscInt lidx[NLOAD_AT_BUS_MAX]; // the list of load ids. Why support multiple and cap it at 1?
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar)); // align for arraylike access

typedef struct _p_VERTEX_Power *VERTEX_Power;

struct _p_LOAD
{
  PetscInt bus_i; // bus number
  char i[20]; // bus number "extended" name
  char id[20]; // load identifier, in case of multiple loads (set to 1 currently)
  PetscInt status; // What values can this take?
  PetscInt area; // area to which the load is assigned
  PetscInt zone; // zone ""...but what are the areas/zones?
  // How are these determined/enforce to be mutually exclusive?
  // e.g. what happens if both the constant MVA and constant current are given values?
  PetscScalar pl; // active power of a constant MVA load
  PetscScalar ql; // reactive component of contant MVA load
  PetscScalar ip; // active power of constant current load
  PetscScalar iq; // reactive power of constant current load
  PetscScalar yp; // "" constant admittance load
  PetscScalar yq; // ""
  PetscScalar scale_load; // No clue what this does yet.
  PetscInt owner;
  PetscInt internal_i;
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar)); // align for arraylike access


typedef struct _p_LOAD *LOAD;

// So this struct apparently can support multiple owners?
// Asking a power engineer they said this MAY be related to power schedules
// but I'm not sure. Need to ask about it. At any rate, its not REALLY a core
// concept of the sim.
struct _p_GEN
{
  PetscInt bus_i;
  char i[20];
  char id[20];
  PetscScalar pg; // active power out
  PetscScalar qg; // reactive power out
  PetscScalar qt; // max reactive power out
  PetscScalar qb; // MIN reactive power out
  PetscScalar vs; // regulated voltage set point
  PetscInt ireg; // remote bus #/ID
  PetscScalar mbase; // base MVA of machine, for PU?
  PetscScalar   zr; // Complex machine impedance ZSOURCE in pu on mbase
  PetscScalar   zx; // ----------------------"-------------------------
  PetscScalar   rt; // Step-up transformer impedance XTRAN in pu on mbase
  PetscScalar   xt; // -----------------------"--------------------------
  PetscScalar   gtap; // Step-up transformer turns ratio
  PetscInt      status; // Machine status
  PetscScalar   rmpct; // Mvar % required to hold voltage at remote bus
  PetscScalar   pt; // Gen max active power output: MW
  PetscScalar   pb; // Gen min active power output: MW
  PetscInt      o1; // Owner number
  PetscScalar   f1; // Fraction of ownership
  PetscScalar   scale_gen;
  PetscInt      internal_i; // Internal Bus Number
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _p_GEN *GEN;

struct _p_EDGE_Power
{
  PetscInt fbus; // no clue yet what this is
  PetscInt tbus; // ""
  // Bus numbers
  // Requires 2 because this is an edge, not a vertex!
  char i[20];
  char j[20];
  char ckt[20]; // circuit identifier
  PetscScalar r; // branch resistance, PU
  PetscScalar x; // branch reactance, PU
  PetscScalar b; // branch "charging suseptance", pu
  PetscScalar rateA; // rate a in MVA (rating?)
  PetscScalar rateB; // rate b in MVA (rating?)
  PetscScalar rateC; // rate c in MVA (rating?)
  PetscScalar tapratio; // transformer info?
  PetscScalar phaseshift; // for phase shifting transformers. What's the units, rads or degrees?
  PetscScalar gi; // complex admittance at "i" bus, pu
  PetscScalar bi; // ""
  PetscScalar gj; // complex admittance at "j" bus, pu
  PetscScalar bj; // ""
  PetscInt status; // service status (probably what the "status" bit in each means?
  PetscScalar length; // line length
  PetscInt o1; // owner number
  PetscScalar yff[2], yft[2], ytf[2], ytt[2]; // complex, of the form [G, B]
  PetscInt internal_i; // internal "from" bus number
  PetscInt internal_j; // internal "to" bus number
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _p_EDGE_Power *EDGE_Power;

// "PTI" format data structure
// Looks like PTI is a Siemens subsidiary
// Description from UW here:
// https://labs.ece.uw.edu/pstca/formats/pti.txt
struct _PFDATA
{
  PetscScalar sbase; // Base MVA for the system
  // These values include elements that are out of service
  PetscInt nbus, ngen, nbranch, nload;

  VERTEX_Power bus;
  LOAD load;
  GEN gen;
  EDGE_Power branch;
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

// Make this a bit more readable, I almost missed the actual type lol
typedef struct _PFDATA PFDATA;

#endif//MY_POWER_H
