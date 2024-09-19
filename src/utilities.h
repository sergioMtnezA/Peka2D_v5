#include <stdio.h>
#include <stdlib.h>
#include "define.h"
#include "structs.h"

// ENDIAN ##########################################################################
#ifndef _SWAP_ENDIAN
#define _SWAP_ENDIAN

EXPORT_DLL int IsLittleEndian(void);

/******************************************************************************
  FUNCTION: SwapEndian
  PURPOSE: Swap the byte order of a structure
  EXAMPLE: float F=123.456;; SWAP_FLOAT(F);
******************************************************************************/

#define SWAP_SHORT(Var)  Var = *(short*)         SwapEndian((void*)&Var, sizeof(short))
#define SWAP_USHORT(Var) Var = *(unsigned short*)SwapEndian((void*)&Var, sizeof(short))
#define SWAP_LONG(Var)   Var = *(long*)          SwapEndian((void*)&Var, sizeof(long))
#define SWAP_ULONG(Var)  Var = *(unsigned long*) SwapEndian((void*)&Var, sizeof(long))
#define SWAP_INT(Var)    Var = *(int*)           SwapEndian((void*)&Var, sizeof(int))
#define SWAP_FLOAT(Var)  Var = *(float*)         SwapEndian((void*)&Var, sizeof(float))
#define SWAP_DOUBLE(Var) Var = *(double*)        SwapEndian((void*)&Var, sizeof(double))

EXPORT_DLL void *SwapEndian(void* Addr, const int Nb);

#endif
// ENDIAN ##########################################################################
