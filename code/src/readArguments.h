/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*    This file is part of "mipsymmetries" a collection of                   */
/*    symmetry handling methods.                                             */
/*                                                                           */
/*    Copyright (C) 2010-2017  Marc Pfetsch, Thomas Rehn, Christopher Hojny  */
/*                                                                           */
/*                                                                           */
/*    Based on SCIP  --- Solving Constraint Integer Programs                 */
/*       SCIP is distributed under the terms of the SCIP Academic Licence,   */
/*       see file COPYING in the SCIP distribution.                          */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"

/**@file   readArguments.h
 * @brief  read comand line arguments
 * @author Marc Pfetsch
 */

#ifndef __READARGUMENTS_H__
#define __READARGUMENTS_H__

#include <scip/scip.h>

#ifdef __cplusplus
extern "C" {
#endif


/** get problem name
 *
 *  Returns 0 on error
 */
extern
int getProblemName(
   const char*           filename,           /**< filename */
   char*                 probName,           /**< problem name (on exit) */
   int                   maxSize             /**< maximal size of probName string */
   );

/** read comand line arguments */
extern
SCIP_RETCODE readArguments(
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array with shell parameters */
   const char**          filename,           /**< file name from arguments */
   const char**          solutionfile,       /**< solution file */
   const char**          writesolfilename,   /**< write solution file name */
   const char**          settingsname,       /**< name of settings file */
   SCIP_Real*            timeLimit,          /**< time limit read from arguments */
   SCIP_Real*            memLimit,           /**< memory limit read from arguments */
   SCIP_Longint*         nodeLimit,          /**< node limit read from arguments */
   int*                  dispFreq,           /**< display frequency */
   SCIP_Bool*            onlypre,            /**< Only run preprocessing? */
   int*                  permseed,           /**< seed for permutations */
   int*                  randseed,           /**< seed for randomization */
   SCIP_Real*            cutoffvalue         /**< value of cutoff (SCIP_INVALID otherwise) */
   );

#ifdef __cplusplus
}
#endif

#endif
