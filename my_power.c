static char help[] = 
"This is a re-write of the SNES example from Petsc\n
An attempt was made to C++-ize it, but their use of array-like access\n
of structs is incompatible with C++.\n\n
Instead, I'll do my best to explain each line of non-obvious code.\n
See the source for deets.\n"

#include "power.h"
#include <petscdmnetwork.h>

// Forward declare some stuff
PetscErrorCode FormFunction(SNES snes, Vec X, Vec F, void *appctx);
PetscErrorCode SetInitialValues(DM networkdm, Vec X, void *appctx);

int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  // TODO make this an option
  char pfdata_file[PETSC_MAX_PATH_LEN] = "case9.m";
  PFDATA *pfdata;
  PetscInt numEdges = 0, numVertices = 0;
  PetscInt *edges = NULL;
  PetscInt i;
  DM networkdm;
  UserCtx_Power User;
  // Stages in the sense "we're in phase 1" not like a "stage" for a play
  // Why multiple stages though? why not just a single logger?
  PetscLogStage stage1, stage2;
}

// snes is the SNES object (duh)
// X is the current state to evaluate from
// F is where vector where the result is stored
// ctx is optional "user defined function context"; likely a way to
//   hand in arbitrary values for the evaluation. yay void ptr!
//   short for "application context" I imagine?
PetscErrorCode FormFunction(SNES snes, Vec X, Vec F, void *appctx);
{
  // TODO
  PetscErrorCode ierr;
  DM networkdm;
  UserCtx_Power *User = (UserCtx_Power*) appctx;
  Vec localX, localF;
  PetscInt nv, ne; // TODO -> num_edges, num_vertices
  const PetscInt *vtx, *edges; // TODO *vertices

  // Required at the start of each Petsc function
  // https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscFunctionBegin.html
  PetscFunctionBegin; 
  // Why don't they just do CHKERRQ(func_call())??
  // I imagine the decision to not evaluate automatically is for flexibility
  // C++ would allow exceptions, just saying...
  ierr = SNESGetDM(snes, &networkdm); CHKERRQ(ierr);
  // The use of local/global, etc. is due to the support for MPI
  // SHOULD handle "ghost" values automagically.
  ierr = DMGetLocalVector(networkdm, &localX); CHKERRQ(ierr);
  ierr = DMGetLocalVector(networkdm, &localF); CHKERRQ(ierr);
  // Set all elements of F (output vec) and localF to 0
  // acts like memset (e.g. sets all to the same scalar value)
  ierr = VecSet(F, 0.0); CHKERRQ(ierr);
  ierr = VecSet(localF, 0.0); CHKERRQ(ierr);

  // So it looks like the "get" operation gets the location of the vector, but this
  // is required to actually fill it?
  ierr = DMGlobalToLocalBegin(networkdm, X, INSERT_VALUES, localX); CHKERRQ(ierr);
  // And, it appears that it can be done async. (start, do some stuff, then come back). Neat.
  ierr = DMGlobalToLocalEnd(networkdm, X, INSERT_VALUES, localX); CHKERRQ(ierr);

  // Get the "sub" network (it IS an MPI app after all)
  // The literal "0" is a PetscInt for the subnetwork ID. Probably should be the MPI rank in the future?
  ierr = DMNetworkGetSubnetworkInfo(networkdm, 0, &nv, &ne, &vtx, &edges); CHKERRQ(ierr);
  // And this is probably the magick itself. Will annotate the implementation itself
  ierr = FormFunction_Power(networkdm, localX, localF, nv, ne, vtx, edges, User); CHKERRQ(ierr);

  // Hm. This probably returns the new state vector calculated in FormFunction_Power?
  ierr = DMRestoreLocalVector(networkdm, &localX); CHKERRQ(ierr);

  // Insert the current residual to the global
  // mode: "if ADD_VALUES ghost points from the same base point accumulate into the base point"
  //    So an additive merge of ghost points? the alternative is no parallel comms
  ierr = DMLocalToGlobalBegin(networkdm, localF, ADD_VALUES, F); CHKERRQ(ierr);
  // Again, its an async thing!
  ierr = DMLocalToGlobalEnd(networkdm, localF, ADD_VALUES, F); CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm, &localF); CHKERRQ(ierr);

  // Required at the end of PETSc functions
  // https://i.kym-cdn.com/entries/icons/facebook/000/022/978/yNlQWRM.jpg
  PetscFunctionReturn(0);
}

// Sets the initial state of the vector X
PetscErrorCode SetInitialValues(DM networkdm, Vec X, void *appctx)
{
  // So, this is used CONSTANTLY...why the local scope? Besides "avoiding globals"
  PetscErrorCode ierr;
  PetscInt vStart, vEnd, nv, ne;
  const PetscInt *vtx, *edges;
  Vec localX;
  UserCtx_Power *user_power = (UserCtx_Power*) appctx;

  PetscFunctionBegin;
  // So...I'm not sure what this is doing. Obviously getting the "vertex range" but
  // how are the vertices actually being represented here? Linked List?
  // And the results of this aren't even used as far as I can tell!
  ierr = DMNetworkGetVertexRange(networkdm, &vStart, &vEnd); CHKERRQ(ierr);

  // Self explanatory
  ierr = DMGetLocalVector(networkdm, &localX); CHKERRQ(ierr);

  ierr = VecSet(X, 0.0); CHKERRQ(ierr);
  // Seems a convenience function might be helpful (or maybe bloaty) for this style of calling
  // e.g. DMGlobalToLocal(...) w/o *Begin or *End
  ierr DMGlobalToLocalBegin(networkdm, X, INSERT_VALUES, localX); 
  ierr DMGlobalToLocalEnd(networkdm, X, INSERT_VALUES, localX); 

  ierr = DMNetworkGetSubnetworkInfo(networkdm, 0, &nv, &ne, &vtx, &edges); CHKERRQ(ierr);
  // This is from the pfffunctions, will be dissected there
  ierr = SetInitialGuess_Power(networkdm, localX, nv, ne, vtx, edges, user_power); CHKERRQ(ierr);

  // Same as before
  ierr = DMLocalToGlobalBegin(networkdm, localX, ADD_VALUES, X); CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(networkdm, localX, ADD_VALUES, X); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

