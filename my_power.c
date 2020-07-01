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
  PetscMPIInt rank;
  PetscInt eStart, eEnd, vStart, vEnd, j;
  PetscInt genj, loadj;
  Vec X, F;
  Mat J;
  SNES snes;

  // Pretty much required at the beginning. inits everything, including MPI
  // the database file, "poweroptions" here, seems to contain command line arguments commonly used
  // More investigation needed for specifics
  ierr = PetscInitialize(&argc, &argv, "poweroptions", help); if (ierr) return ierr;
  // Get the rank. Looks like there's a specific comm world for Petsc
  // Can multiple petsc apps run simultaneously on a cluster automagically? (probably not)
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  {
    // scope this section for...reasons? The comments mention something about the clang static analyzer. Dunno.
    const PetscMPIInt crank = rank; //that's c(onst) rank, not crank. I think.
    // Set up the Network object. This grid is going to be our...admittance matrix, I believe.
    ierr = DMNetworkCreate(PETSC_COMM_WORLD, &networkdm); CHKERRQ(ierr);
    // Register components. This is where we define what types can be represented (as well as what their
    // representations are) in the network. Here, its the various electrical components -- branches,
    // busses, generators, and loads. In general, though, this could be things like, idk, pipes and turbines
    // and filters in a water pipe system. Or blood vessels. (Finite elements in a heat flow problem, maybe?)
    // The output (last arg) is a (unique? universally unique?) ID to indicate the types.
    ierr = DMNetworkRegisterComponent(networkdm,
                                      "branchstruct",
                                      sizeof(struct _p_EDGE_Power),
                                      &User.compkey_branch); CHKERRQ(ierr);
    // I always thought of busses as vertexes...but, I guess (given PQ, PV, busses) that's not correct
    // They have unique behaviors rather than being a pure connection
    ierr = DMNetworkRegisterComponent(networkdm,
                                      "busstruct",
                                      sizeof(struct _p_VERTEX_POWER),
                                      &User.compkey_branch); CHKERRQ(ierr);
    // Why are neither of these "_p_VERTEX_*" or "_p_EDGE_*"?
    // Probably detailed in the header.
    ierr = DMNetworkRegisterComponent(networkdm,
                                      "genstruct",
                                      sizeof(struct _p_GEN),
                                      &User.compkey_branch); CHKERRQ(ierr);
    ierr = DMNetworkRegisterComponent(networkdm,
                                      "loadstruct",
                                      sizeof(struct _p_LOAD),
                                      &User.compkey_branch); CHKERRQ(ierr);

    // This is an interesting pattern. Why the stage and push? This ain't Git!
    ierr = PetscLogStageRegister("Read Data", &stage1); CHKERRQ(ierr);
    PetscLogStagePush(stage1);
    if (!crank)
    {
      // Only the zeroeth rank readeth the data
      // options database, NULL for default global db
      // string to prepend to name, or NULL if none
      // name of the option sought
      // maximum length of that string (ah, C strings...)
      // output, location to copy string. So...why's it NULL? Just discards it? why even have this then???
      ierr = PetscOptionsGetString(NULL, NULL, "-pfdata", pfdata_file, sizeof(pfdata_file), NULL); CHKERRQ(ierr);
      // petscnew is a memory aligned new. Kinda like mkl_alloc and such.
      ierr = PetscNew(&pfdata); CHKERRQ(ierr);
      // This is in the pffunctions I think
      ierr = PFReadMatPowerData(pfdata, pfdata_file); CHKERRQ(ierr);
      User.Sbase = pfdata->sbase;

      numEdges = pfdata->nbranch;
      numVertices = pfdata->nbus;

      // Yep, looks like edges are stored as packed tuples (e.g. [(,),(,)...(,)]
      ierr = PetscMalloc(2*numEdges, &edges); CHKERRQ(ierr);
      ierr = GetListOfEdges_Power(pfdata, edges); CHKERRQ(ierr);
    }

    // Similar to above. options db (null for default), prepended string (NULL for none)
    // the option name, and then the output (true or false)
    ierr = PetscOptionsHasName(NULL, NULL, "-jac_error", &User.jac_error); CHKERRQ(ierr);

    // Again, not really sure why the logging works this way, but this pops the stage off the stack
    PetscLogStagePop();
    // Block until all procs are synced up.
    ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);
    // Ah, maybe the log stages are to provide a formal method of delinieation?
    ierr = PetscLogStageRegister("Create network", &stage2); CHKERRQ(ierr);
    PectscLogStagePush(stage2);
    // Set number of nodes and edges
    // dm obj, global number of subnets, num of local verts, num of local edges, global
    //   number of "coupling" (probably how many this attaches to?), number of local edges for each
    //   coupling network (array? value by ref? not sure)
    ierr = DMNetworkSetSizes(networkdm, 1, &numVertices, &numEdges, 0, NULL); CHKERRQ(ierr);
    // Set up the connections
    // Given a list of edges, set 'em up. The last arg is for the coupling subnets
    ierr = DMNetworkSetEdgeList(networkdm, &edges, NULL); CHKERRQ(ierr);
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

