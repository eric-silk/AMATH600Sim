#include "pffunctions.h"

PetscErrorCode GetListOfEdges_Power(PFDATA* pfdata, PetscInt* edgelist)
{
  PetscErrorCode ierr;
  PetscInt i, fbus, tbus, nbranches = pfdata->nbranch;
  EDGE_Power branch = pfdata->branch;
  PetscBool netview = PETSC_FALSE;

  PetscFunctionBegin;
  // This option will print out the network as a list of connected nodes
  // For instance:
  //  node1 -> node2
  //  node2 -> node3
  //  node3 -> node1
  ierr = PetscOptionsHasName(NULL, NULL, "-powernet_view", &netview); CHKERRQ(ierr);
  for (i = 0; i < nbranches; ++i)
  {
    fbus = branch[i].internal_i;
    tbus = branch[i].internal_j;
    edgelist[2*i] = fbus;
    edgelist[2*i+1] = tbus;
    if (netview)
    {
      ierr = PetscPrintf(PETSC_COMM_SELF, "branch %d, bus[%d] -> bus[%d]\n", i, fbus, tbus); CHKERRQ(ierr);
    }
  }
  if (netview)
  {
    for (i = 0; i < pfdata->nbus; i++)
    {
      if (pfdata->bus[i].ngen)
      {
        ierr = PetscPrintf(PETSC_COMM_SELF, " bus %D: gen\n", i); CHKERRQ(ierr);
      }
      else if (pfdata->bus[i].nload)
      {
        ierr = PetscPrintf(PETSC_COMM_SELF, " bus %D: load\n", i); CHKERRQ(ierr);
      }
    }
  }
  PetscFunctionReturn(0);
}

// Gee, Encapsulation would be awesome right now!
PetscErrorCode _Private_FormJacobian_Power(DM networkdm,
                                           Vec localX,
                                           Mat J,
                                           PetscInt nv,
                                           PetscInt ne,
                                           const PetscInt* vtx,
                                           const PetscInt* edges,
                                           void* appctx)
{
  PetscErrorCode ierr;
  const PetscScalar *xarr;
  PetscInt i, v, row[2], col[8], e, vfrom, vto;
  PetscInt offsetfrom, offsetto, goffsetfrom, goffsetto, numComps;
  PetscScalar values[8];
  PetscInt j, key, offset, goffset;
  PetscScalar Vm;
  UserCtx_Power *user_power = (UserCtx_Power*)appctx;
  PetscScalar Sbase = user_power->Sbase;
  VERTEX_Power bus;
  PetscBool ghostvtex;
  void* component;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(localX, &xarr); CHKERRQ(ierr);

  for (v = 0; v < nv; ++v)
  {
    // Checks if the vertex of interest is on another node (a ghost point)
    ierr = DMNetworkIsGhostVertex(networkdm, vtx[v], &ghostvtex); CHKERRQ(ierr);
    // The number of components at a vertex or edge
    // I assume "components" here are loads, generators, etc.
    ierr = DMNetworkGetNumComponents(networkdm, vtx[v], &numComps); CHKERRQ(ierr);
    for (j = 0; j < numComps; ++j)
    {
      // each local vector contains several components; get the offset
      // require to access the component of interest
      // This is why they forced alignment of the structs, I imagine
      ierr = DMNetworkGetVariableOffset(networkdm, vtx[v], &offset); CHKERRQ(ierr);
      // Same, but globally (why's this useful?)
      ierr = DMNetworkGetVariableGlobalOffset(networkdm, vtx[v], &goffset); CHKERRQ(ierr);
      ierr = DMNetworkGetComponent(networkdm, vtx[v], j, &key, &component); CHKERRQ(ierr);

      if (key == user_power->compkey_bus)
      {
        PetscInt nconnedges;
        const PetscInt *connedges;

        bus = (VERTEX_Power)(component);
        if (!ghostvtex)
        {
          // I THINK this is the way they're handling P/Theta - QU decoupling?
          // That, or more likely, they're setting the PU Voltage to 1 and angle to 0
          if (bus->ide == REF_BUS || bus->ide == ISOLATED_BUS)
          {
            row[0] = goffset;
            row[1] = goffset + 1;
            col[0] = col[2] = goffset;
            col[1] = col[3] = goffset + 1;
            values[0] = values[3] = 1.0;
            values[1] = values[2] = 0.0;
            ierr = MatSetValues(J, 2, row, 2, col, values, ADD_VALUES); CHKERRQ(ierr);
            break;
          }

          Vm = xarr[offset + 1];

          // Shunt injections
          row[0] = col[0] = goffset;
          row[1] = col[1] = goffset + 1;
          values[0] = values[1] = values[2] = values[3] = 0.0;
          if (bus->ide != PV_BUS)
          {
            // d/dVm(Vm^2*gl)
            // Original note: "Shunt injections"
            values[1] = 2.0 * Vm * bus->gl / Sbase;
            values[3] = -2.0 * Vm * bus->bl / Sbase;
          }
          ierr = MatSetValues(J, 2, row, 2, col, values, ADD_VALUES); CHKERRQ(ierr);
        }

        ierr = DMNetworkGetSupportingEdges(networkdm, vtx[v], &nconnedges, &connedges); CHKERRQ(ierr);

        for (i = 0; i < nconnedges; ++i)
        {
          EDGE_Power branch;
          VERTEX_Power busf, bust;
          PetscInt keyf, keyt;
          PetscScalar Gff, Bff, Gft, Bft, Gtf, Btf, Gtt, Btt;
          const PetscInt *cone;
          PetscScalar Vmf, Vmt, thetaf, thetat, thetaft, thetatf;

          e = connedges[i];
          ierr = DMNetworkGetComponent(networkdm, e, 0, &key, (void**)&branch); CHKERRQ(ierr);
          // Only add the branch if its in service
          if (!branch->status) continue;

          Gff = branch->yff[0];
          Bff = branch->yff[1];
          Gft = branch->yft[0];
          Bft = branch->yft[1];
          Gtf = branch->ytf[0];
          Btf = branch->ytf[1];
          Gtt = branch->ytt[0];
          Btt = branch->ytt[1];

          ierr = DMNetworkGetConnectedVertices(networkdm, edges[e], &cone); CHKERRQ(ierr);
          vfrom = cone[0];
          vto = cone[1];

          ierr = DMNetworkGetVariableOffset(networkdm, vfrom, &offsetfrom); CHKERRQ(ierr);
          ierr = DMNetworkGetVariableOffset(networkdm, vto, &offsetto); CHKERRQ(ierr);
          ierr = DMNetworkGetVariableGlobalOffset(networkdm, vfrom, &goffsetfrom); CHKERRQ(ierr);
          ierr = DMNetworkGetVariableGlobalOffset(networkdm, vto, &goffsetto); CHKERRQ(ierr);

          if (goffsetto < 0 )
          {
            goffsetto = -goffsetto - 1;
          }

          thetaf = xarr[offsetfrom];
          Vmf = xarr[offsetfrom+1];
          thetat = xarr[offsetto];
          Vmt = xarr[offsetto+1];
          thetaft = thetaf - thetat;
          thetatf = thetat - thetaf;

          ierr = DMNetworkGetComponent(networkdm, vfrom, 0, &keyf, (void**)&busf); CHKERRQ(ierr);
          ierr = DMNetworkGetComponent(networkdm, vto, 0, &keyt, (void**)&bust); CHKERRQ(ierr);

          if (vfrom == vtx[v])
          {
            if (busf->ide != REF_BUS)
            {
              row[0] = col[0] = goffsetfrom;
              col[1] = goffsetfrom + 1;
              col[2] = goffsetto;
              col[3] = goffsetto + 1;
              // df/dthetaf
              values[0] = Vmf * Vmt * (Gft*-PetscSinScalar(thetaft) + Bft*PetscCosScalar(thetaft));
              // df/dVmf
              values[1] = -2.0 * Gff * Vmf + Vmt*(Gft * PetscCosScalar(thetaft) + Bft * PetscSinScalar(thetaft));
              // df/dthetat
              values[2] = Vmf*Vmt*(Gft*PetscSinScalar(thetaft) + Bft*-PetscCosScalar(thetaft));
              // df/dVmt
              values[3] = Vmf*(Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
              ierr = MatSetValues(J, 1, row, 4, col, values, ADD_VALUES); CHKERRQ(ierr);
            }

            if (busf->ide != PV_BUS && busf->ide != REF_BUS)
            {
              row[0] = goffsetfrom + 1;
              col[0] = goffsetfrom;
              col[1] = goffsetfrom + 1;
              col[2] = goffsetto;
              col[3] = goffsetto + 1;
              values[0] = Vmf * Vmt * (Bft * PetscSinScalar(thetaft) + Gft*PetscCosScalar(thetaft));
              values[1] = -2.0 * Bff * Vmf + Vmt * (-Bft*PetscCosScalar(thetaft) + Gft * PetscSinScalar(thetaft));
              values[2] = Vmf * Vmt * (-Bft * PetscSinScalar(thetaft) + Gft * - PetscCosScalar(thetaft));
              values[3] = Vmf * (-Bft *PetscCosScalar(thetaft) + Gft * PetscSinScalar(thetaft));

              ierr = MatSetValues(J, 1, row, 4, col, values, ADD_VALUES); CHKERRQ(ierr);
            }
          }
          else // vfrom != vtx[v]
          {
            if (bust->ide != REF_BUS)
            {
              row[0] = goffsetto;
              col[0] = goffsetto;
              col[1] = goffsetto + 1;
              col[2] = goffsetfrom;
              col[3] = goffsetfrom + 1;

              // df/dthetat
              values[0] = Vmt * Vmf * (Gtf * -PetscSinScalar(thetatf) + Btf*PetscCosScalar(thetaft));
              // df/dVmt
              values[1] = 2.0 * Gtt * Vmt + Vmf * (Gtf * PetscCosScalar(thetatf) + Btf * PetscSinScalar(thetatf));
              // df/dthetaf
              values[2] = Vmt * Vmf * (Gtf*PetscSinScalar(thetatf) + Btf*-PetscCosScalar(thetatf));
              // df/dVmf
              values[3] = Vmt * (Gtf*PetscCosScalar(thetatf) + Btf * PetscSinScalar(thetatf));
            }
            if (bust->ide != PV_BUS && bust->ide != REF_BUS)
            {
              row[0] = goffsetto + 1;
              col[0] = goffsetto;
              col[1] = goffsetto + 1;
              col[2] = goffsetfrom;
              col[3] = goffsetfrom + 1;
	      values[0] =  Vmt * Vmf * (Btf * PetscSinScalar(thetatf) + Gtf * PetscCosScalar(thetatf));
	      values[1] =  -2.0 * Btt * Vmt + Vmf * (-Btf * PetscCosScalar(thetatf) + Gtf * PetscSinScalar(thetatf));
	      values[2] =  Vmt * Vmf * (-Btf * PetscSinScalar(thetatf) + Gtf * -PetscCosScalar(thetatf));
	      values[3] =  Vmt * (-Btf * PetscCosScalar(thetatf) + Gtf * PetscSinScalar(thetatf));

              ierr = MatSetValues(J, 1, row, 4, col, values, ADD_VALUES); CHKERRQ(ierr);
            }
          }
        }

        if (!ghostvtex && bus->ide == PV_BUS)
        {
          row[0] = goffset + 1;
          col[0] = goffset + 1;
          values[0] = 1.0;
          // Here's where the (optional) error is set
          if (user_power->jac_error)
          {
            values[0] = 50.0;
          }

          ierr = MatSetValues(J, 1, row, 1, col, values, ADD_VALUES); CHKERRQ(ierr);
        }
      }
    }
  }

  ierr = VecRestoreArrayRead(localX, &xarr); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode FormJacobian_Power(SNES snes, Vec X, Mat J, Mat Jpre, void* appctx)
{
  PetscErrorCode ierr;
  DM networkdm;
  Vec localX;
  PetscInt nv, ne;
  const PetscInt *vtx, *edges;

  PetscFunctionBegin;
  ierr = MatZeroEntries(J); CHKERRQ(ierr);

  ierr = SNESGetDM(snes, &networkdm); CHKERRQ(ierr);
  ierr = DMGetLocalVector(networkdm, &localX); CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(networkdm, X, INSERT_VALUES, localX); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(networkdm, X, INSERT_VALUES, localX); CHKERRQ(ierr);

  ierr = DMNetworkGetSubnetworkInfo(networkdm, 0, &nv, &ne, &vtx, &edges); CHKERRQ(ierr);
  ierr = _Private_FormJacobian_Power(networkdm, localX, J, nv, ne, vtx, edges, appctx); CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(networkdm, &localX); CHKERRQ(ierr);

  // So this actually assembles the matrix. Looks like the other functions to add/insert elements
  // just caches them, and the matrix itself isn't ready to go until these functions are called.
  // There are some special considerations if you intersperse insert vs. add directives
  // https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatAssemblyBegin.html
  ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode FormFunction_Power(DM networkdm, Vec localX, Vec localF, PetscInt nv, PetscInt ne, const PetscInt* vtx, const PetscInt* edges, void* appctx)
{
  PetscErrorCode ierr;
  UserCtx_Power *User = (UserCtx_Power*) appctx;
  PetscInt  e, v, vfrom, vto;
  const PetscScalar *xarr;
  PetscScalar *farr;
  PetscInt offsetfrom, offsetto, offset, i, j, key, numComps;
  PetscScalar Vm;
  PetscScalar Sbase=User->Sbase;
  VERTEX_Power bus=NULL;
  GEN gen;
  LOAD load;
  PetscBool ghostvtex;
  void* component;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(localX, &xarr); CHKERRQ(ierr);
  ierr = VecGetArray(localF, &farr); CHKERRQ(ierr);

  for (v=0; v<nv; v++) 
  {
    ierr = DMNetworkIsGhostVertex(networkdm, vtx[v], &ghostvtex); CHKERRQ(ierr);
    ierr = DMNetworkGetNumComponents(networkdm, vtx[v], &numComps); CHKERRQ(ierr);
    ierr = DMNetworkGetVariableOffset(networkdm, vtx[v], &offset); CHKERRQ(ierr);

    for (j = 0; j < numComps; j++)
    {
      ierr = DMNetworkGetComponent(networkdm, vtx[v], j, &key, &component); CHKERRQ(ierr);
      if (key == User->compkey_bus)
      {
        PetscInt       nconnedges;
	const PetscInt *connedges;

	bus = (VERTEX_Power)(component);
	/* Handle reference bus constrained dofs */
	if (bus->ide == REF_BUS || bus->ide == ISOLATED_BUS)
        {
	  farr[offset] = xarr[offset] - bus->va*PETSC_PI/180.0;
	  farr[offset+1] = xarr[offset+1] - bus->vm;
	  break;
	}

	if (!ghostvtex)
        {
	  Vm = xarr[offset+1];

	  /* Shunt injections */
          // Vm^2*gl, normalized to PU
          // TODO
          // Curious note:
          // both are added in the Jacobian but ONLY if its not a PV bus
          // Here, the first is ALWAYS added, irrespective of whether its a PV bus or not
	  farr[offset] += Vm*Vm*bus->gl/Sbase;
	  if(bus->ide != PV_BUS)
          {
            farr[offset+1] += -Vm*Vm*bus->bl/Sbase;
          }
	}

	ierr = DMNetworkGetSupportingEdges(networkdm, vtx[v], &nconnedges, &connedges); CHKERRQ(ierr);

	for (i=0; i < nconnedges; i++)
        {
	  EDGE_Power       branch;
	  PetscInt       keye;
          PetscScalar    Gff, Bff, Gft, Bft, Gtf, Btf, Gtt, Btt;
          const PetscInt *cone;
          PetscScalar    Vmf, Vmt, thetaf, thetat, thetaft, thetatf;

	  e = connedges[i];
	  ierr = DMNetworkGetComponent(networkdm, e, 0, &keye, (void**)&branch); CHKERRQ(ierr);
	  if (!branch->status)
          {
            continue;
          }

	  Gff = branch->yff[0];
	  Bff = branch->yff[1];
	  Gft = branch->yft[0];
	  Bft = branch->yft[1];
	  Gtf = branch->ytf[0];
	  Btf = branch->ytf[1];
	  Gtt = branch->ytt[0];
	  Btt = branch->ytt[1];

	  ierr = DMNetworkGetConnectedVertices(networkdm, e, &cone); CHKERRQ(ierr);
	  vfrom = cone[0];
	  vto = cone[1];

	  ierr = DMNetworkGetVariableOffset(networkdm, vfrom, &offsetfrom); CHKERRQ(ierr);
	  ierr = DMNetworkGetVariableOffset(networkdm, vto, &offsetto); CHKERRQ(ierr);

	  thetaf = xarr[offsetfrom];
	  Vmf = xarr[offsetfrom+1];
	  thetat = xarr[offsetto];
	  Vmt = xarr[offsetto+1];
	  thetaft = thetaf - thetat;
	  thetatf = thetat - thetaf;

          // TODO
          // These LOOK similar to equations 6.101 on pg 263 of Kundur's book?
	  if (vfrom == vtx[v])
          {
	    farr[offsetfrom]   += Gff*Vmf*Vmf + Vmf*Vmt*(Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
	    farr[offsetfrom+1] += -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
	  }
          else
          {
	    farr[offsetto]   += Gtt*Vmt*Vmt + Vmt*Vmf*(Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));
	    farr[offsetto+1] += -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
	  }
	}
      }
      else if (key == User->compkey_gen)
      {
	if (!ghostvtex)
        {
	  gen = (GEN)(component);
	  if (!gen->status)
          {
            continue;
          }
	  farr[offset] += -gen->pg/Sbase;
	  farr[offset+1] += -gen->qg/Sbase;
	}
      }
      else if (key == User->compkey_load)
      {
	if (!ghostvtex)
        {
	  load = (LOAD)(component);
	  farr[offset] += load->pl/Sbase;
	  farr[offset+1] += load->ql/Sbase;
	}
      }
    }
    if (bus && bus->ide == PV_BUS)
    {
      farr[offset+1] = xarr[offset+1] - bus->vm;
    }
  }
  ierr = VecRestoreArrayRead(localX, &xarr); CHKERRQ(ierr);
  ierr = VecRestoreArray(localF, &farr); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SetInitialGuess_Power(DM networkdm,
                                     Vec localX,
                                     PetscInt nv,
                                     PetscInt ne,
                                     const PetscInt *vtx,
                                     const PetscInt *edges,
                                     void* appctx)
{
  PetscErrorCode ierr;
  VERTEX_Power bus;
  PetscInt i;
  GEN gen;
  PetscBool ghostvtex;
  PetscScalar *xarr;
  PetscInt key, numComps, j, offset;
  void* component;
  PetscMPIInt rank;
  MPI_Comm comm;
  UserCtx_Power *User = (UserCtx_Power*) appctx;
  
  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject) networkdm, &comm); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);
  ierr = VecGetArray(localX, &xarr); CHKERRQ(ierr);
  for (i = 0; i < nv; ++i)
  {
    ierr = DMNetworkIsGhostVertex(networkdm, vtx[i], &ghostvtex); CHKERRQ(ierr);
    if (ghostvtex)
    {
      continue;
    }

    ierr = DMNetworkGetVariableOffset(networkdm, vtx[i], &offset); CHKERRQ(ierr);
    ierr = DMNetworkGetNumComponents(networkdm, vtx[i], &numComps); CHKERRQ(ierr);

    for (j = 0; j < numComps; ++j)
    {
      ierr = DMNetworkGetComponent(networkdm,vtx[i],j,&key,&component);CHKERRQ(ierr);
      if (key == User->compkey_bus)
      {
        bus = (VERTEX_Power)(component);
        xarr[offset] = bus->va * PETSC_PI/180.0;
        xarr[offset+1] = bus->vm;
      }
      else if (key == User->compkey_gen)
      {
        gen = (GEN)(component);
        if (!gen->status)
        {
          continue;
        }
        xarr[offset+1] = gen->vs;
        break;
      }
    }
  }

  ierr = VecRestoreArray(localX, &xarr); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

