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

/**@file   readArguments.c
 * @brief  read comand line arguments
 * @author Marc Pfetsch
 */

#include "readArguments.h"
#include <string.h>
#include <assert.h>

/** get problem name
 *
 *  Returns 0 on error
 */
int getProblemName(
   const char*           filename,           /**< filename */
   char*                 probName,           /**< problem name (on exit) */
   int                   maxSize             /**< maximal size of probName string */
   )
{
   int i = 0;
   int result = 1;
   int l;
   int j;

   /* first find end of string */
   while ( filename[i] != 0) /* find end of string */
      ++i;
   l = i;                    /* end of string */

   /* if we found ".gz" */
   if ( l > 3 && filename[l-3] == '.' && filename[l-2] == 'g' && filename[l-1] == 'z' )
   {
      l = l - 4;
      i = l;
   }

   /* go back until '.' or '/' appears */
   while ((i > 0) && (filename[i] != '.') && (filename[i] != '/'))
      --i;
   assert(i > 0);

   /* if we found '.', search for '/' */
   if (filename[i] == '.')
   {
      l = i;
      while ((i > 0) && (filename[i] != '/'))
	 --i;
   }

   /* correct counter */
   if (filename[i] == '/')
      ++i;

   /* copy name */
   j = 0;
   while ( (i < l) && (filename[i] != 0) )
   {
      probName[j++] = filename[i++];
      if ( j >= maxSize-1 )
      {
	 result = 0;
	 break;
      }
   }
   probName[j] = 0;
   return result;
}


/** read comand line arguments */
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
   SCIP_Real*            cutoffvalue         /**< value of cutoff if setcutoff is true */
   )
{  /*lint --e{818}*/
   int i;
   char usage[SCIP_MAXSTRLEN];
   int status;

   assert( argc > 0 );
   assert( argv != NULL );
   assert( filename != NULL );
   assert( solutionfile != NULL );
   assert( writesolfilename != NULL );
   assert( settingsname != NULL );
   assert( timeLimit != NULL );
   assert( memLimit != NULL );
   assert( nodeLimit != NULL );
   assert( dispFreq != NULL );
   assert( onlypre != NULL );
   assert( permseed != NULL );
   assert( cutoffvalue != NULL );

   /* init usage text */
   status = snprintf(usage, SCIP_MAXSTRLEN, "usage: %s <file> [-l <solution file>] [-w <write solution file>] [-s <setting file>] [-t <time limit>] [-m <mem limit>] [-n <node limit>] [-d <display frequency>] [-p <seed>] [-setcutoff <value>] [-O]", argv[0]);
   if ( status < 0 || status > SCIP_MAXSTRLEN )
   {
      SCIPerrorMessage("string not long enough to hold usage message.\n");
      return SCIP_OKAY;
   }

   /* init arguments */
   *timeLimit = 1e20;
   *memLimit = 1e20;
   *nodeLimit = SCIP_LONGINT_MAX;
   *filename = NULL;
   *solutionfile = NULL;
   *writesolfilename = NULL;
   *settingsname = NULL;
   *dispFreq = -1;
   *onlypre = FALSE;
   *permseed = -1;
   *randseed = -1;
   *cutoffvalue = SCIP_INVALID;

   /* check all arguments */
   for (i = 1; i < argc; ++i) /*lint -e850*/
   {
      /* check for solution name */
      if ( ! strcmp(argv[i], "-l") )
      {
         if ( *solutionfile != NULL )
         {
            fprintf(stderr, "%s\n", usage);
            return SCIP_ERROR;
         }
         if ( i == argc-1 )
         {
            fprintf(stderr, "No solution file name supplied.\n");
            fprintf(stderr, "%s\n", usage);
            return SCIP_ERROR;
         }
         ++i;
         *solutionfile = argv[i];
         assert( i < argc );
      }
      /* check for solution name */
      else if ( ! strcmp(argv[i], "-w") )
      {
         if ( *writesolfilename != NULL )
         {
            fprintf(stderr, "%s\n", usage);
            return SCIP_ERROR;
         }
         if ( i == argc-1 )
         {
            fprintf(stderr, "No solution file name to write to supplied.\n");
            fprintf(stderr, "%s\n", usage);
            return SCIP_ERROR;
         }
         ++i;
         *writesolfilename = argv[i];
         assert( i < argc );
      }
      /* check for settings name */
      else if ( ! strcmp(argv[i], "-s") )
      {
         if ( *settingsname != NULL )
         {
            fprintf(stderr, "%s\n", usage);
            return SCIP_ERROR;
         }
         if ( i == argc-1 )
         {
            fprintf(stderr, "No setting file name supplied.\n");
            fprintf(stderr, "%s\n", usage);
            return SCIP_ERROR;
         }
         ++i;
         *settingsname = argv[i];
         assert( i < argc );
      }
      /* check for time limit */
      else if ( ! strcmp(argv[i], "-t") )
      {
         if ( i == argc-1 )
         {
            fprintf(stderr, "No time limit supplied.\n");
            fprintf(stderr, "%s\n", usage);
            return SCIP_ERROR;
         }
         ++i;
         *timeLimit = atof(argv[i]);
         assert( i < argc );
      }
      /* check for memory limit */
      else if ( ! strcmp(argv[i], "-m") )
      {
         if ( i == argc-1 )
         {
            fprintf(stderr, "No memory limit supplied.\n");
            fprintf(stderr, "%s\n", usage);
            return SCIP_ERROR;
         }
         ++i;
         *memLimit = atof(argv[i]);
         assert( i < argc );
      }
      /* check for node limit */
      else if ( ! strcmp(argv[i], "-n") )
      {
         if ( i == argc-1 )
         {
            fprintf(stderr, "No node limit supplied.\n");
            fprintf(stderr, "%s\n", usage);
            return SCIP_ERROR;
         }
         ++i;
         *nodeLimit = atol(argv[i]);
         assert( i < argc );
      }
      /* check for display frequency */
      else if ( ! strcmp(argv[i], "-d") )
      {
         if ( i == argc-1 )
         {
            fprintf(stderr, "No display frequency supplied.\n");
            fprintf(stderr, "%s\n", usage);
            return SCIP_ERROR;
         }
         ++i;
         *dispFreq = atoi(argv[i]);
         assert( i < argc );
      }
      /* check for permutation seed */
      else if ( ! strcmp(argv[i], "-p") )
      {
         if ( i == argc-1 )
         {
            fprintf(stderr, "No permutation seed given.\n");
            fprintf(stderr, "%s\n", usage);
            return SCIP_ERROR;
         }
         ++i;
         *permseed = atoi(argv[i]);
         assert( i < argc );
      }
      /* check for permutation seed */
      else if ( ! strcmp(argv[i], "-seed") )
      {
         if ( i == argc-1 )
         {
            fprintf(stderr, "No random seed given.\n");
            fprintf(stderr, "%s\n", usage);
            return SCIP_ERROR;
         }
         ++i;
         *randseed = atoi(argv[i]);
         assert( i < argc );
      }
      else if ( ! strcmp(argv[i], "-setcutoff") )
      {
         if ( i == argc-1 )
         {
            fprintf(stderr, "No cutoff value supplied.\n");
            fprintf(stderr, "%s\n", usage);
            return SCIP_ERROR;
         }
         ++i;
         *cutoffvalue = atof(argv[i]);
         assert( i < argc );
      }
      else if ( ! strcmp(argv[i], "-O") )
      {
         ++i;
         *onlypre = TRUE;
      }
      else
      {
         /* if filename is already specified */
         if ( *filename != NULL )
         {
            fprintf(stderr, "Filename alread specified.\n");
            return SCIP_ERROR;
         }
         else
            *filename = argv[i];
      }
   }

   if ( *filename == NULL )
   {
      fprintf(stderr, "No filename supplied.\n");
      fprintf(stderr, "%s\n", usage);
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}
