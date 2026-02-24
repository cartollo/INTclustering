for k in {2..10}; do
  Rscript compDisMatrix.R -c  results/maximum_ward.D2_NOWHITE_CZM_isclr/out_create_dismatrix_maximum_ward.D2_NOWHITE_CZM_iscore_isclr.rds -f results/maximum_ward.D2_NOWHITE_CZM_isclr/out_create_dismatrix_maximum_ward.D2_NOWHITE_CZM_isclr.rds -n $k
done