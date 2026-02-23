for k in {2..10}; do
  Rscript compDisMatrix.R -c results/mahalanobis_ward.D2_NOWHITE_CZM_isclr/out_create_dismatrix_mahalanobis_ward.D2_NOWHITE_CZM_iscore_isclr.rds -f results/mahalanobis_ward.D2_NOWHITE_CZM_isclr/out_create_dismatrix_mahalanobis_ward.D2_NOWHITE_CZM_isclr.rds -n $k
done