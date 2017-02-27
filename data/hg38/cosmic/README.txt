# The COSMIC database now requires users to agree to a license (and pay a fee for non-academic use).
# If you want RADIA to annotate mutations that are present in COSMIC, please download the hg38 
# COSMIC database and create .bed files where the first column is the chromosome, the second 
# column is the startCoordinate, the third column is the stopCoordinate, and the fourth column
# is the COSMIC Id.  Be aware of the off-by-one problem.  The RADIA VCFs are 1-based and most
# database dumps are 0-based.  The RADIA VCF coordinate will match the stopCoordinate column
# in the COSMIC .bed file.  If you need assistance, please contact me.
