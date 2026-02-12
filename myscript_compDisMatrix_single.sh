for k in {2..10}; do
  Rscript compDisMatrix.R -c results/bray_curtis_ward.D_NOWHITE_skip_ispca/out_create_dismatrix_bray_curtis_ward.D_NOWHITE_skip_iscore_ispca.rds -f results/bray_curtis_ward.D_NOWHITE_skip_ispca/out_create_dismatrix_bray_curtis_ward.D_NOWHITE_skip_ispca.rds -n $k
done