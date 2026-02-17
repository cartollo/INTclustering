for k in {5..5}; do
  Rscript compDisMatrix.R -c out_create_dismatrix_euclidean_ward.D2_NOWHITE_CZM_iscore_isclr.rds -f results/mahalanobis_ward.D_NOWHITE_CZM_iscore_isclr_isrelab_sil/out_create_dismatrix_mahalanobis_ward.D_NOWHITE_CZM_isclr_isrelab_0.0.05.rds -n $k
done