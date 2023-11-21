/**@file   main.cpp
 * @brief  main file for symretope propagation code
 * @author Jasper van Doornmalen, Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <scip/scip.h>

#include <scip/scipdefplugins.h>
#include "readArguments.h"
#include "cons_symretope.h"

#include <iostream>
#include <fstream>

/** run scip with commandline arguments */
static
SCIP_RETCODE runSCIP(
   int                   argc,
   char**                argv
   )
{
   SCIP* scip = 0;
   const char* filename;
   const char* settingsname;
   const char* solutionfile;
   const char* writesolfilename;
   SCIP_Real timeLimit;
   SCIP_Real memLimit;
   SCIP_Longint nodeLimit;
   SCIP_Bool onlypre = FALSE;
   SCIP_RETCODE retcode;
   SCIP_Real cutoffvalue = SCIP_INVALID;
   int dispFreq;
   int permseed;
   int randseed;

   /* parse command line arguments */
   retcode = readArguments(argc, argv, &filename, &solutionfile, &writesolfilename, &settingsname,
      &timeLimit, &memLimit, &nodeLimit, &dispFreq, &onlypre, &permseed, &randseed, &cutoffvalue);
   if ( retcode != SCIP_OKAY )
      exit(1);
   assert( filename != 0 );

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   SCIPprintVersion(scip, 0);

   SCIPinfoMessage(scip, 0, "\n");
   SCIPinfoMessage(scip, 0, "Symretope propagation methods - (c) Jasper van Doornmalen, Christopher Hojny.\n");
   SCIPinfoMessage(scip, 0, "[GitHash: %s]\n", SYMGITHASH);
   SCIPinfoMessage(scip, 0, "\n");

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL( SCIPincludeConshdlrSymretope(scip) );


   /* --------------------------------------------------------------------- */
   /* handle permutations */
   if ( permseed >= 0 )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "randomization/permutationseed", permseed) );
   }

   if ( randseed >= 0 )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "randomization/randomseedshift", randseed) );
   }

   /* --------------------------------------------------------------------- */

   /* set time, node, and memory limit */
   if ( ! SCIPisInfinity(scip, timeLimit) )
      SCIP_CALL( SCIPsetRealParam(scip, "limits/time", timeLimit) );
   if ( ! SCIPisInfinity(scip, memLimit) )
      SCIP_CALL( SCIPsetRealParam(scip, "limits/memory", memLimit) );
   if ( nodeLimit < SCIP_LONGINT_MAX )
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", nodeLimit) );
   if ( dispFreq >= 0 )
      SCIP_CALL( SCIPsetIntParam(scip, "display/freq", dispFreq) );

   /* check for parameters */
   if ( settingsname != 0 )
   {
      if ( SCIPfileExists(settingsname) )
      {
         SCIPinfoMessage(scip, 0, "reading parameter file <%s> ...\n\n", settingsname);
         SCIP_CALL( SCIPreadParams(scip, settingsname) );
      }
      else
      {
         SCIPerrorMessage("parameter file <%s> not found - using default parameters.\n", settingsname);
      }
   }

   /* output changed parameters */
   SCIPinfoMessage(scip, 0, "Changed settings:\n");
   SCIP_CALL( SCIPwriteParams(scip, 0, FALSE, TRUE) );
   SCIPinfoMessage(scip, 0, "\n");

   /* solve problem */
   if ( ! onlypre )
      SCIPinfoMessage(scip, 0, "\nsolving problem ...\n\n");
   else
      SCIPinfoMessage(scip, 0, "\nrunning preprocessing ...\n\n");

   /* read problem */
   SCIP_CALL( SCIPreadProb(scip, filename, 0) );

   /* possibly read solution */
   if ( solutionfile != 0 )
   {
      SCIP_CALL( SCIPreadSol(scip, solutionfile) );
   }

   /* handle cutoff */
   if ( cutoffvalue != SCIP_INVALID )
   {
      SCIPinfoMessage(scip, 0, "\nSetting cutoff value to %g.\n\n", cutoffvalue);
      SCIP_CALL( SCIPsetObjlimit(scip, cutoffvalue) );
      SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   }

   SCIP_CALL( SCIPsetIntParam(scip, "constraints/symresack/sepafreq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/symretope/sepafreq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/orbisack/sepafreq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/orbitope/sepafreq", -1));

   if ( onlypre )
   {
      SCIP_CALL( SCIPpresolve(scip) );
   }
   else
   {
      SCIP_CALL( SCIPsolve(scip) );
   }

   SCIP_CALL( SCIPprintStatistics(scip, 0) );
#if 0
   if ( SCIPgetBestSol(scip) != 0 )
   {
      SCIP_CALL( SCIPprintSol(scip, SCIPgetBestSol(scip), 0, FALSE) );
   }
#endif

   if ( writesolfilename != NULL )
   {
      FILE* file;

      file = fopen(writesolfilename, "w");
      if( file == NULL )
      {
         SCIPwarningMessage(scip, "error creating file <%s>\n", writesolfilename);
      }
      else
      {
         SCIP_Bool printzeros;

         SCIPinfoMessage(scip, file, "solution status: ");
         SCIP_CALL_FINALLY( SCIPprintStatus(scip, file), fclose(file) );

         SCIP_CALL_FINALLY( SCIPgetBoolParam(scip, "write/printzeros", &printzeros), fclose(file) );

         SCIPinfoMessage(scip, file, "\n");
         SCIP_CALL_FINALLY( SCIPprintBestSol(scip, file, printzeros), fclose(file) );

         SCIPinfoMessage(scip, NULL, "written solution information to file <%s>\n", writesolfilename);
         fclose(file);
      }
   }

   SCIP_CALL( SCIPfreeProb(scip) );

   // SCIPprintMemoryDiagnostic(scip);

   /* free */
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}



/** main function */
int
main(
   int                   argc,
   char**                argv
   )
{
   SCIP_RETCODE retcode;

   retcode = runSCIP(argc, argv);
   if ( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
