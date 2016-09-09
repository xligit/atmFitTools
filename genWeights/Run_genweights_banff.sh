#!/bin/tcsh

if ($#argv == 0) then
	./genweights_BANFF.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/numu_x_numu/root/'*'.root /disk2/usr5/xiaoyue/skmc/test_fc/missert_rootfiles/postfit_data_1p1h_biascorrection_20160310.root /disk2/usr5/xiaoyue/t2kmc14a/nu_mode/spline_weight
	./genweights_BANFF.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/numubar_x_numubar/root/'*'.root /disk2/usr5/xiaoyue/skmc/test_fc/missert_rootfiles/postfit_data_1p1h_biascorrection_20160310.root /disk2/usr5/xiaoyue/t2kmc14a/nu_mode/spline_weight
	./genweights_BANFF.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/nue_x_nue/root/'*'.root /disk2/usr5/xiaoyue/skmc/test_fc/missert_rootfiles/postfit_data_1p1h_biascorrection_20160310.root /disk2/usr5/xiaoyue/t2kmc14a/nu_mode/spline_weight
	./genweights_BANFF.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/nuebar_x_nuebar/root/'*'.root /disk2/usr5/xiaoyue/skmc/test_fc/missert_rootfiles/postfit_data_1p1h_biascorrection_20160310.root /disk2/usr5/xiaoyue/t2kmc14a/nu_mode/spline_weight
	./genweights_BANFF.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/numu_x_nue/root/'*'.root /disk2/usr5/xiaoyue/skmc/test_fc/missert_rootfiles/postfit_data_1p1h_biascorrection_20160310.root /disk2/usr5/xiaoyue/t2kmc14a/nu_mode/spline_weight
	./genweights_BANFF.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/numubar_x_nuebar/root/'*'.root /disk2/usr5/xiaoyue/skmc/test_fc/missert_rootfiles/postfit_data_1p1h_biascorrection_20160310.root /disk2/usr5/xiaoyue/t2kmc14a/nu_mode/spline_weight

	./genweights_BANFF.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/numu_x_numu/root/'*'.root /disk2/usr5/xiaoyue/skmc/test_fc/missert_rootfiles/postfit_data_1p1h_biascorrection_20160310.root /disk2/usr5/xiaoyue/t2kmc14a/antinu_mode/spline_weight
	./genweights_BANFF.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/numubar_x_numubar/root/'*'.root /disk2/usr5/xiaoyue/skmc/test_fc/missert_rootfiles/postfit_data_1p1h_biascorrection_20160310.root /disk2/usr5/xiaoyue/t2kmc14a/antinu_mode/spline_weight
	./genweights_BANFF.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/nue_x_nue/root/'*'.root /disk2/usr5/xiaoyue/skmc/test_fc/missert_rootfiles/postfit_data_1p1h_biascorrection_20160310.root /disk2/usr5/xiaoyue/t2kmc14a/antinu_mode/spline_weight
	./genweights_BANFF.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/nuebar_x_nuebar/root/'*'.root /disk2/usr5/xiaoyue/skmc/test_fc/missert_rootfiles/postfit_data_1p1h_biascorrection_20160310.root /disk2/usr5/xiaoyue/t2kmc14a/antinu_mode/spline_weight
	./genweights_BANFF.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/numu_x_nue/root/'*'.root /disk2/usr5/xiaoyue/skmc/test_fc/missert_rootfiles/postfit_data_1p1h_biascorrection_20160310.root /disk2/usr5/xiaoyue/t2kmc14a/antinu_mode/spline_weight
	./genweights_BANFF.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/numubar_x_nuebar/root/'*'.root /disk2/usr5/xiaoyue/skmc/test_fc/missert_rootfiles/postfit_data_1p1h_biascorrection_20160310.root /disk2/usr5/xiaoyue/t2kmc14a/antinu_mode/spline_weight


else

    ./genweights_BANFF.pl /disk2/usr5/xiaoyue/skmc/test_fc/root/'*'.root /disk2/usr5/xiaoyue/skmc/test_fc/missert_rootfiles/postfit_data_1p1h_biascorrection_20160310.root /disk2/usr5/xiaoyue/skmc/test_fc/weights

endif
