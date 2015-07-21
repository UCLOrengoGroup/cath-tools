/************************************************************************/
/*>WHOLEPDB *doReadWholePDB(FILE *fp, int *natom, BOOL AllAtoms, 
                            int OccRank, int ModelNum, BOOL doWhole)
   -----------------------------------------------------------------
   Input:   FILE     *fp      A pointer to type FILE in which the
                              .PDB file is stored.
            BOOL     AllAtoms TRUE:  ATOM & HETATM records
                              FALSE: ATOM records only
            int      OccRank  Occupancy ranking
            int      ModelNum NMR Model number (0 = all)
            BOOL     doWhole  Read headers and footers
   Output:  int      *natom   Number of atoms read. -1 if error.
   Returns: WHOLEPDB *wpdb    A pointer to the WHOLE PDB structure
                              containing the PDB linked list

   Reads a PDB file into a PDB linked list. The OccRank value indicates
   occupancy ranking to read for partial occupancy atoms.
   If any partial occupancy atoms are read the global flag 
   gPDBPartialOcc is set to TRUE.

   04.11.88 V1.0  Original
   07.02.89 V1.1  Ignore records which aren't ATOM or HETATM
   28.03.90 V1.2  Altered field widths to match PDB standard better
                  See notes above for deviations
   28.06.90 V1.2a Buffer size increased to 85 chars.
   15.02.91 V1.2b Changed comment header to match new standard.
   07.01.92 V1.3  Ignores blank lines properly
   11.05.92 V1.4  Check on EOF in while() loop, memset() buffer. 
                  ANSIed.
   01.06.92 V1.5  Documented for autodoc
   19.06.92 V1.6  Corrected use of stdlib
   01.10.92 V1.7  Changed to use fgets()
   10.06.93 V1.9  Returns 0 on failure rather than exiting
                  Replaced SIZE with sizeof(PDB) directly
   17.06.93 V2.0  Rewritten to use fsscanf()
   08.07.93 V2.1  Split from ReadPDB()
   09.07.93 V2.2  Modified to return pointer to PDB. Rewrote allocation
                  scheme.
   17.03.94 V2.3  Handles partial occupancies
                  Sets natom to -1 if there was an error to distinguish 
                  from no atoms.
                  Handles atom names which start in column 13 rather
                  than column 14. This is allowed in the standard, but
                  very rare.
                  Sets flag for partials.
   06.04.94 V2.4  Atom names starting in column 13 have their first
                  character moved to the end if it is a digit.
   03.10.94 V2.5  Check residue number as well as atom name when running
                  through alternative atoms for partial occupancy
                  Moved increment of NPartial, so only done if there
                  is space in the array. If OccRank is 0, all atoms are
                  read regardless of occupancy.
   06.03.95 V2.7  Added value for NMR model to read (0 = all)
                  No longer static. Sets gPDBMultiNMR if ENDMDL records
                  found.
   13.01.97 V2.8  Added check on return from fsscanf. Blank lines used
                  to result in duplication of the previous line since
                  fsscanf() does not reset the variables on receiving
                  a blank line. Also fixed in fsscanf().
   25.02.98 V2.9  Added code to read gzipped PDB files transparently
                  when GUNZIP_SUPPORT is defined
   17.08.98 V2.10 Added case to popen() for SunOS
   08.10.99 V2.11 Initialise CurIns and CurRes
   15.02.01 V2.12 Added atnam_raw
   27.04.05 V2.14 Added another atnam_raw for multiple occupancies
   03.06.05 V2.15 Added altpos
   14.10.05 V2.16 Modified detection of partial occupancy. handles
                  residues like 1zeh/B16 where a lower partial is
                  erroneously set to zero
   12.03.07 V3.0  Modified from doReadPDB() to extract headers and footers
                  doReadPDB() is now a wrapper to this
   29.04.08 V-.-. Modified doReadWholePDB() to set pointers wpdb->header
                  and wpdb->trailer to NULL and to set wpdb->natoms.
   02.05.08 V-.-  Added SIGATM, ANISOU and SIGUIJ recordtypes for PDB
                  Format Version. 3.0.1 files. 

*/
WHOLEPDB *doReadWholePDB(FILE *fpin,
                         int  *natom,
                         BOOL AllAtoms,
                         int  OccRank,
                         int  ModelNum,
                         BOOL doWhole)
{
   char     record_type[8],
            atnambuff[8],
            *atnam,
            atnam_raw[8],
            resnam[8],
            chain[4],
            insert[4],
            buffer[160],
            CurAtom[8],
            cmd[80],
            CurIns = ' ',
            altpos;
   int      atnum,
            resnum,
            CurRes = 0,
            NPartial,
            ModelCount = 1;
   FILE     *fp = fpin;
   double   x,y,z,
            occ,
            bval;
   WHOLEPDB *wpdb = NULL;
   PDB      *pdb  = NULL,
            *p,
            multi[MAXPARTIAL];   /* Temporary storage for partial occ   */

#ifdef GUNZIP_SUPPORT
   int      signature[3],
            i,
            ch;
#endif

   *natom         = 0;
   CurAtom[0]     = '\0';
   NPartial       = 0;
   gPDBPartialOcc = FALSE;
   gPDBMultiNMR   = FALSE;
   cmd[0]         = '\0';

#ifdef GUNZIP_SUPPORT
   /* See whether this is a gzipped file                                */
   for(i=0; i<3; i++)
      signature[i] = fgetc(fpin);
   for(i=2; i>=0; i--)
      ungetc(signature[i], fpin);
   if((signature[0] == (int)0x1F) &&
      (signature[1] == (int)0x8B) &&
      (signature[2] == (int)0x08))
   {
      /* It is gzipped so we'll open gunzip as a pipe and send the data
         through that into a temporary file
      */
      sprintf(cmd,"gunzip >/tmp/readpdb_%d",(int)getpid());
      if((fp = (FILE *)popen(cmd,"w"))==NULL)
      {
         *natom = (-1);
         return(NULL);
      }
      while((ch=fgetc(fpin))!=EOF)
         fputc(ch, fp);
      pclose(fp);

      /* We now reopen the temporary file as our PDB input file         */
      sprintf(cmd,"/tmp/readpdb_%d",(int)getpid());
      if((fp = fopen(cmd,"r"))==NULL)
      {
         *natom = (-1);
         return(NULL);
      }
   }
#endif   



   /* Parse PDB File */
   while(fgets(buffer,159,fp))
   {


     /* Parse Coordinate Records */
     
     /* Deal with Model Numbers  */
      if(ModelNum != 0)          /* We are interested in model numbers  */
      {
         if(!strncmp(buffer,"ENDMDL",6))
         {
            ModelCount++;
         }

         if(ModelCount < ModelNum)   /* Haven't reached the right model */
            continue;
         else if(ModelCount > ModelNum)    /* Gone past the right model */
            break;
      }

      if(!strncmp(buffer,"ENDMDL",6))
         gPDBMultiNMR   = TRUE;
      

      /* Parse ATOM and HETATM Records */
      if(!strncmp(buffer, "ATOM  ", 6) ||
         !strncmp(buffer, "HETATM", 6))
      {
         if(fsscanf(buffer,
                    "%6s%5d%1x%5s%4s%1s%4d%1s%3x%8lf%8lf%8lf%6lf%6lf",
                    record_type,&atnum,atnambuff,resnam,chain,&resnum,
                    insert,&x,&y,&z,&occ,&bval) != EOF)
         {
            if((!strncmp(record_type,"ATOM  ",6)) || 
               (!strncmp(record_type,"HETATM",6) && AllAtoms))
            {
               /* Copy the raw atom name                                */
               /* 03.06.05 Note: this reads the alternate atom position 
                  as well as the atom name - changes in FixAtomName() now
                  strip that. We now copy only the first 4 characters 
                  into atnam_raw and put the 5th character into altpos
               */
               strncpy(atnam_raw, atnambuff, 4);
               atnam_raw[4] = '\0';
               altpos = atnambuff[4];
               
               /* Fix the atom name accounting for start in column 
                  13 or 14
               */
               atnam = FixAtomName(atnambuff, occ);
               
               /* Check for full occupancy. If occupancy is 0.0 assume 
                  that it is actually fully occupied; the column just 
                  hasn't been filled in correctly
                  
                  04.10.94 Read all atoms if OccRank is 0
                  
                  14.10.05 Now takes an atom as full occupancy:
                           if occ==1.0
                           if occ==0.0 and altpos==' '
                           if OccRank==0
                           This fixes problems where a lower (partial)
                           occupancy has erroneously been set to zero
               */
               if(((altpos == ' ') && (occ < (double)SMALL)) ||
                  (occ > (double)0.999) || 
                  (OccRank == 0))
               {
                  /* Trim the atom name to 4 characters                 */
                  atnam[4] = '\0';
                  
                  if(NPartial != 0)
                  {
                     if(!StoreOccRankAtom(OccRank,multi,NPartial,&pdb,&p,
                                          natom))
                     {
                        if(pdb != NULL) FREELIST(pdb, PDB);
                        *natom = (-1);
                        if(cmd[0]) unlink(cmd);
                        return(NULL);
                     }
                     
                     /* Set partial occupancy counter to 0              */
                     NPartial = 0;
                  }
                  
                  /* Allocate space in the linked list                  */
                  if(pdb == NULL)
                  {
                     INIT(pdb, PDB);
                     p = pdb;
                  }
                  else
                  {
                     ALLOCNEXT(p, PDB);
                  }
                  
                  /* Failed to allocate space; free up list so far & 
                     return
                  */
                  if(p==NULL)
                  {
                     if(pdb != NULL) FREELIST(pdb, PDB);
                     *natom = (-1);
                     if(cmd[0]) unlink(cmd);
                     return(NULL);
                  }
                  
                  /* Increment the number of atoms                      */
                  (*natom)++;
                  
                  /* Store the information read                         */
                  p->atnum  = atnum;
                  p->resnum = resnum;
                  p->x      = (REAL)x;
                  p->y      = (REAL)y;
                  p->z      = (REAL)z;
                  p->occ    = (REAL)occ;
                  p->bval   = (REAL)bval;
                  p->altpos = altpos;    /* 03.06.05 Added this one     */
                  p->next   = NULL;
                  strcpy(p->record_type, record_type);
                  strcpy(p->atnam,       atnam);
                  strcpy(p->atnam_raw,   atnam_raw);
                  strcpy(p->resnam,      resnam);
                  strcpy(p->chain,       chain);
                  strcpy(p->insert,      insert);
               }
               else   /* Partial occupancy                              */
               {
                  /* Set flag to say we've got a partial occupancy atom */
                  gPDBPartialOcc = TRUE;
                  
                  /* First in a group, store atom name                  */
                  if(NPartial == 0)
                  {
                     CurIns = insert[0];
                     CurRes = resnum;
                     strncpy(CurAtom,atnam,8);
                  }
                  
                  if(strncmp(CurAtom,atnam,strlen(CurAtom)-1) || 
                     resnum != CurRes || 
                     CurIns != insert[0])
                  {
                     /* Atom name has changed 
                        Select and store the OccRank highest occupancy 
                        atom
                     */
                     if(!StoreOccRankAtom(OccRank,multi,NPartial,&pdb,&p,
                                          natom))
                     {
                        if(pdb != NULL) FREELIST(pdb, PDB);
                        *natom = (-1);
                        if(cmd[0]) unlink(cmd);
                        return(NULL);
                     }
                     
                     /* Reset the partial atom counter                  */
                     NPartial = 0;
                     strncpy(CurAtom,atnam,8);
                     CurRes = resnum;
                     CurIns = insert[0];
                  }
                  
                  if(NPartial < MAXPARTIAL)
                  {
                     /* Store the partial atom data                     */
                     multi[NPartial].atnum  = atnum;
                     multi[NPartial].resnum = resnum;
                     multi[NPartial].x      = (REAL)x;
                     multi[NPartial].y      = (REAL)y;
                     multi[NPartial].z      = (REAL)z;
                     multi[NPartial].occ    = (REAL)occ;
                     multi[NPartial].bval   = (REAL)bval;
                     multi[NPartial].next   = NULL;
                     strcpy(multi[NPartial].record_type, record_type);
                     strcpy(multi[NPartial].atnam,       atnam);
                     /* 27.04.05 - added this line                      */
                     strcpy(multi[NPartial].atnam_raw,   atnam_raw);
                     strcpy(multi[NPartial].resnam,      resnam);
                     strcpy(multi[NPartial].chain,       chain);
                     strcpy(multi[NPartial].insert,      insert);
                     /* 03.06.05 - added this line                      */
                     multi[NPartial].altpos = altpos;
                     
                     NPartial++;
                  }
               }
               
            }
         }
	 continue;
      }


      /* Parse SIGATM, ANISOU and SIGUIJ records */
      if(!strncmp(buffer, "SIGATM", 6) ||
         !strncmp(buffer, "ANISOU", 6) ||
         !strncmp(buffer, "SIGUIJ", 6))
	{
	  continue;
	}


      /* Parse TER records */
      if(!strncmp(buffer, "TER", 3))
	{
	  continue;
	}


      /* Store Header and Trailer */
      if(doWhole)
	{
	  if(wpdb == NULL)
	    {
	      /* Create the wholepdb structure                         */
	      if((wpdb = (WHOLEPDB *)malloc(sizeof(WHOLEPDB)))==NULL)
		{
		  FREELIST(pdb, PDB);
		  return(NULL);
		}
	      /* Added 29.04.08 */
	      wpdb->header  = NULL;
	      wpdb->trailer = NULL;
	    }

	  if(*natom)    /* Footer                                     */
	    {
	      if((wpdb->trailer = StoreString(wpdb->trailer, buffer))
		 == NULL)
		return(NULL);
	    }
	  else          /* Header                                     */
	    {
	      if((wpdb->header = StoreString(wpdb->header, buffer))
		 == NULL)
		return(NULL);
	    }
	}
       

   } /* End of Parse PDB File loop */


   if(NPartial != 0)
   {
      if(!StoreOccRankAtom(OccRank,multi,NPartial,&pdb,&p,natom))
      {
         if(pdb != NULL) FREELIST(pdb, PDB);
         *natom = (-1);
         if(cmd[0]) unlink(cmd);
         return(NULL);
      }
   }

   if(cmd[0]) unlink(cmd);


   if(wpdb == NULL)
   {
      /* Create the wholepdb structure                         */
      if((wpdb = (WHOLEPDB *)malloc(sizeof(WHOLEPDB)))==NULL)
      {
         FREELIST(pdb, PDB);
         return(NULL);
      }
      /* Added 29.04.08 */
      wpdb->header  = NULL;
      wpdb->trailer = NULL;
   }
   wpdb->natoms = *natom; /* Added 29.04.08 */
   wpdb->pdb = pdb;

   /* Return pointer to start of linked list                            */
   return(wpdb);
}

