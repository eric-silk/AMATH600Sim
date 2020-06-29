#include "power_example.h"

// TODO change raw ptr to std::array
// This doesn't "get" anything...it just prints the netlist if desired.
// TODO rename appropriately
PetscErrorCode GetListOfEdges_Power(const PFDATA& pfdata, PetscInt *edgelist)
{
  PetscErrorCode ierr;
  PetscInt i, fbus, tbus, nbranches = pfdata.nbranch;
  EDGE_Power branch = pfdata.branch;
  PetscBool netview = PETSC_FALSE;

  PetscFunctionBegin();
  ierr = PetscOptionsHasName(NULL, NULL, "-powernet_view", &netview); CKERRQ(ierr);
  for (i = 0; i < nbranches; ++i)
  {
    fbus = branch[i].internal_i;
    tbus = branch[i].internal_j;
    // 2i and 2i+1 cuz they're stored as [(,),(,)...(,)]
    edge_list[2*i] = fbus;
    edge_list[2*i+1] = tbus;
    if (netview)
    {
      ierr = PetscPrintf(PETSC_COMM_SELF, "branch %d, bus[%d] -> bus[%d]\n", i , fbus, tbus); CHKERRQ(ierr);
    }
    if (netview)
    {
      for (i = 0; i < pfdata.nbus; ++i)
      {
        if (pfdata.bus[i].ngen)
        {
          ierr = PetscPrintf(PETSC_COMM_SELF, " bus: %D: gen\n", i); CHKERRQ(ierr);
        }
        else if (pfdata.bus[i].nload)
        {
          ierr = PetscPrintf(PETSC_COMM_SELF, " bus %D: load\n", i); CHKERRQ(ierr);
        }
        else
        {
          ierr = PetscPrintf(PETSC_COMM_SELF, " bus %D: Unknown!\n", i); CHKERRQ(ierr);
        }
      }
    }
    PetscFunctionReturn(0);
}


