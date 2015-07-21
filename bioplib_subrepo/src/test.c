#include <stdio.h>
#include "bioplib/pdb.h"

int main(int argc, char **argv)
{
   WHOLEPDB *wpdb;
   PDB *pdb;
   FILE *in;
   int natoms;
   
   in=fopen("test.pdb", "r");
   
   pdb = ReadPDB(in, &natoms);
   WritePDB(stdout, pdb);
   
/*
   wpdb = ReadWholePDB(in);
   
   WriteWholePDB(stdout, wpdb);
*/

   return(0);
}
