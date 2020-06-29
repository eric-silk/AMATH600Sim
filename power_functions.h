#ifndef POWER_FUNCTIONS_H
#define POWER_FUNCTIONS_H

// This is meant to contain all the helper/driver functions to be used in the sim

#include "power_example.h"
#include <petscdmnetwork>

// Why is it called "Form Function" is it "forming" a function or
// is it THE form function? a snes-form function? idk.
PetscErrorCode FormFunction(SNES snes, Vec X, Vec F, void *appctx)
{
  // TODO
}

#endif//POWER_FUNCTIONS_H
