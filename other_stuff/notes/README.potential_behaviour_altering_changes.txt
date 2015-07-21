Investigate:
Changes involved in moving to reading dssp+sec
Changes involved in dropping secs
Changes arising from investigating if it's a problem that sec files give values of 'H' or 'S' but the SSAP source code compares against 'H' and 'E'
Changes arising from fixing the asymmetric results spotted by Benoit


Fairly trivial
Change the superposition to be based on residue pairs with a vaguely decent score (but keep still calculate the RMSD the same way)
