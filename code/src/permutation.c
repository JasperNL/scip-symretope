/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   symmetry.c
 * @ingroup OTHER_CFILES
 * @brief  methods for permutations
 * @author Jasper van Doornmalen
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "permutation.h"
#include "scip/scip.h"


/** Compute the greatest common divisor of two nonnegative integers
 * @param a First integer, nonnegative
 * @param b Second integer, nonnegative
 * @return The greatest common divisor of a and b.
 */
SCIP_Longint gcd(
   SCIP_Longint a,
   SCIP_Longint b
)
{
   assert( a >= 0 );
   assert( b >= 0 );
   do
   {
      if ( a > b )
         a = a % b;
      else
         b = b % a;
   } while( a > 0 && b > 0 );

   assert( a == 0 || b == 0 );
   if ( a > b )
      return a;
   else
      return b;
}

/** Compute the least common multiple of two nonnegative integers
 * @param a First integer, nonnegative
 * @param b Second integer, nonnegative
 * @return The least common multiple of a and b.
 */
SCIP_Longint lcm(
   SCIP_Longint a,
   SCIP_Longint b
)
{
   return a / gcd(a, b) * b;
}

/** Get a permutation object from a permutation array.
 * @param scip The SCIP instance.
 * @param perm The permutation on 0..nvars-1, where entry i of this vector is the value that entry i permutes to.
 * @param nvars The number of variables in the array.
 * @return A SCIP_PERMUTATION object for this permutation.
 */
SCIP_RETCODE SCIPgetPermutation(
   SCIP* scip,
   int* perm,
   int nvars,
   SCIP_PERMUTATION* permutation
)
{
   int i;
   int j;
   SCIP_Longint order;
   int ncycles;
   int prevcyclemaxindex;
   int cyclemaxindex;
   SCIP_Bool ismonotone;
   SCIP_Bool isordered;
   int* checked;
   int thiscyclesize;
   int maxcyclesize;
   int ndescendpoints;
   int** cycles;
   int* cycleblock;
   int* varcycle;
   int* varcyclepos;
   int* cyclelengths;
   // int* cyclelengthsind;
   // int* diffcyclelengths;
   // int ndiffcyclelengths;
   // int maxndiffcyclelengths;
   int cycleid;
   int cycleblockpos;

   /* Careful: This is not a copy! */
   permutation->perm = perm;

   /* Get the group order */
   assert( nvars > 0 );
   order = 1;
   ncycles = 0;
   prevcyclemaxindex = -1;
   cyclemaxindex = -1;
   ismonotone = TRUE;
   isordered = TRUE;

   SCIP_CALL( SCIPallocClearBufferArray(scip, &checked, nvars) );
   maxcyclesize = 0;
   for (i = 0; i < nvars; ++i)
   {
      /* If this index is already processed, skip this iteration. */
      if ( checked[i] )
         continue;

      ++ncycles;
      j = i;
      thiscyclesize = 0;
      ndescendpoints = 0;
      cyclemaxindex = j;
      do
      {
         if ( cyclemaxindex < j )
            cyclemaxindex = j;
         if ( j < prevcyclemaxindex )
            isordered = FALSE;
         if ( perm[j] < j )
            ++ndescendpoints;

         checked[j] = TRUE;
         j = perm[j];
         ++thiscyclesize;
      } while (j != i);
      assert( j == i );

      /* COMMENT: Yes, but then you need to adapt the theory */
      if ( ndescendpoints > 1 )  /* Investigate: Could we also count "a single ascend point" as monotone? */
         ismonotone = FALSE;
      prevcyclemaxindex = cyclemaxindex;

      order = lcm(order, thiscyclesize);
      if ( maxcyclesize < thiscyclesize )
         maxcyclesize = thiscyclesize;
   }

   #ifndef NDEBUG
   /* Debugging: Check if all entries are checked */
   for (i = 0; i < nvars; ++i)
   {
      assert( checked[i] );
   }
   #endif

   SCIPfreeBufferArray(scip, &checked);

   /* If this is an orbitope/orbisack, then do something. */
   if ( order <= 2 )
   {
      /* This is only possible in the manual case, as for symmetry breaking we check if it's an orbisack before */
      SCIPdebugMessage("Permutation defines an orbitope action.");
   }

   permutation->order = order;
   permutation->nvars = nvars;
   permutation->ncycles = ncycles;
   permutation->ismonotone = ismonotone;
   permutation->isordered = isordered;

   /* Compute the cycle decomposition */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &cycles, ncycles ) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &cycleblock, nvars ) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &cyclelengths, ncycles ) );
   // SCIP_CALL( SCIPallocBlockMemoryArray(scip, &cyclelengthsind, ncycles) );
   // /* For the different cycle lengths: How many different cycle lengths could there exist?
   //  * In the worst case, the maximal k, where 1 + 2 + ... + k <= n
   //  * That is a quadratic formula with k = (sqrt(1 + 8n) - 1) / 2 as answer.
   //  */
   // maxndiffcyclelengths = (((int) SQRT((SCIP_Real) (1 + 8 * nvars))) / 2) + 1;
   // SCIP_CALL( SCIPallocBlockMemoryArray(scip, &diffcyclelengths, maxndiffcyclelengths) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &varcycle, nvars ) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &varcyclepos, nvars ) );

   for (i = 0; i < nvars; ++i)
   {
      varcycle[i] = -1;
      varcyclepos[i] = -1;
   }

   cycleid = 0;
   cycleblockpos = 0;
   // ndiffcyclelengths = 0;
   for (i = 0; i < nvars; ++i)
   {
      /* If this index is already processed, don't process. */
      if ( varcycle[i] >= 0 )
         continue;

      /* Store where the cycle starts */
      cycles[cycleid] = &cycleblock[cycleblockpos];

      /* Walk a circle in the cycle */
      j = i;
      thiscyclesize = 0;
      do
      {
         varcycle[j] = cycleid;
         varcyclepos[j] = thiscyclesize;
         cycleblock[cycleblockpos] = j;
         j = perm[j];
         ++thiscyclesize;
         ++cycleblockpos;
      } while (j != i);
      assert( j == i );

      cyclelengths[cycleid] = thiscyclesize;

      /* Check if there is a cycle of this length. */
      // if ( thiscyclesize <= 1 )
      // {
      //    cyclelengthsind[cycleid] = -1;
      // }
      // else
      // {
      //    for (j = 0; j < ndiffcyclelengths; ++j)
      //    {
      //       if ( diffcyclelengths[j] == thiscyclesize )
      //       {
      //          cyclelengthsind[cycleid] = j;
      //          break;
      //       }
      //    }
      //    /* If the loop did not break */
      //    if (j == ndiffcyclelengths)
      //    {
      //       assert( ndiffcyclelengths < maxndiffcyclelengths );
      //       cyclelengthsind[cycleid] = ndiffcyclelengths;
      //       diffcyclelengths[ndiffcyclelengths++] = thiscyclesize;
      //    }
      // }

      ++cycleid;
   }
   assert( cycleid == ncycles );
   assert( cycleblockpos == nvars );

   permutation->cycles = cycles;
   permutation->cycleblock = cycleblock;
   permutation->ncycles = ncycles;
   permutation->cyclelengths = cyclelengths;
   // permutation->cyclelengthsind = cyclelengthsind;
   // permutation->diffcyclelengths = diffcyclelengths;
   // permutation->ndiffcyclelengths = ndiffcyclelengths;
   // permutation->maxndiffcyclelengths = maxndiffcyclelengths;
   permutation->maxcyclesize = maxcyclesize;
   permutation->varcycle = varcycle;
   permutation->varcyclepos = varcyclepos;

   return SCIP_OKAY;
}

/** Given a SCIP_PERMUTATION object, give the entry on which we apply the permutation pow times.
 * @param perm The SCIP_PERMUTATION object
 * @param index The index that we want to permute.
 * @param pow The power for which we want to permute.
 * @return The entry.
 */
int permGet(
   SCIP_PERMUTATION* perm,                   /**< the permutation */
   int index,                                /**< entry to permute */
   int pow                                   /**< power to permute */
)
{
   int cycleid;
   int cyclelen;
   int pos;

   assert( perm != NULL );
   assert( index >= 0 );
   assert( index <= perm->nvars );
   assert( perm->varcycle != NULL );
   assert( perm->cyclelengths != NULL );
   assert( perm->varcyclepos != NULL );
   assert( perm->cycles != NULL );

   cycleid = perm->varcycle[index];
   assert( cycleid < perm->ncycles );

   cyclelen = perm->cyclelengths[cycleid];

   /* Get position of index in cycle */
   pos = perm->varcyclepos[index];
   assert( pos >= 0 );
   assert( pos < cyclelen );

   /* Update this position by "pow" places. */
   pos = (pos + pow) % cyclelen;
   if ( pos < 0 )
      pos += cyclelen;
   assert( pos >= 0 );
   assert( pos < cyclelen );

   /* Return value of new position */
   return perm->cycles[cycleid][pos];
}

/** Given a SCIP_PERMUTATION object, give the permutation array that maps 0..nvars to the permutation raised to a power.
 * @param perm The SCIP_PERMUTATION object
 * @param pow The power for which we want to permute.
 * @param arr The allocated array to store the permutation in.
 * @param nvars The number of variables of the permutation. The array size should match with this.
 * @return SCIP_OKAY if successful.
 */
SCIP_RETCODE getPermArray(
   SCIP_PERMUTATION* perm,
   SCIP_Longint pow,
   int* arr,
   int nvars
)
{
   int* cycle;
   int cyclelen;
   int c;
   int i;

   assert( perm != NULL );
   assert( perm->nvars == nvars );
   assert( perm->order >= 0 );
   assert( perm->cycles != NULL );
   assert( perm->cyclelengths != NULL );
   assert( arr != NULL );
   assert( nvars >= 0 );

   /* Make power positive. */
   pow = pow % perm->order;
   if ( pow < 0 )
      pow += perm->order;

   for (c = 0; c < perm->ncycles; ++c)
   {
      cycle = perm->cycles[c];
      cyclelen = perm->cyclelengths[c];
      for (i = 0; i < cyclelen; ++i)
         arr[cycle[i]] = cycle[(i + pow) % cyclelen];
   }

   return SCIP_OKAY;
}

/** Free a SCIP_PERMUTATION object.
 * @param scip the SCIP instance
 * @param permutation The SCIP_PERMUTATION object
 * @param freePerm Whether the field "perm" of SCIP_PERMUTATION must be freed, as well. \
 *    This is a reference to the original permutation array that was used for this object.
 * @return SCIP_OKAY if successful.
 */
SCIP_RETCODE SCIPfreePermutationContents(
   SCIP* scip,
   SCIP_PERMUTATION* permutation,
   SCIP_Bool freePerm
)
{
   if ( freePerm )
   {
      SCIPfreeBlockMemoryArray(scip, &permutation->perm, permutation->nvars);
   }
   SCIPfreeBlockMemoryArray(scip, &permutation->cycles, permutation->ncycles);
   SCIPfreeBlockMemoryArray(scip, &permutation->cycleblock, permutation->nvars);
   SCIPfreeBlockMemoryArray(scip, &permutation->cyclelengths, permutation->ncycles);
   // SCIPfreeBlockMemoryArray(scip, &permutation->cyclelengthsind, permutation->ncycles);
   // SCIPfreeBlockMemoryArray(scip, &permutation->diffcyclelengths, permutation->maxndiffcyclelengths);
   SCIPfreeBlockMemoryArray(scip, &permutation->varcycle, permutation->nvars);
   SCIPfreeBlockMemoryArray(scip, &permutation->varcyclepos, permutation->nvars);
   return SCIP_OKAY;
}
