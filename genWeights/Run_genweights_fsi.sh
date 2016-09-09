#!/bin/tcsh

if ($#argv == 0) then

    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/numu_x_numu/root/'*'.root 0 t2k +1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/numubar_x_numubar/root/'*'.root 0 t2k +1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/nue_x_nue/root/'*'.root 0 t2k +1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/nuebar_x_nuebar/root/'*'.root 0 t2k +1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/numu_x_nue/root/'*'.root 0 t2k +1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/numubar_x_nuebar/root/'*'.root 0 t2k +1

    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/numu_x_numu/root/'*'.root 1 t2k +1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/numubar_x_numubar/root/'*'.root 1 t2k +1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/nue_x_nue/root/'*'.root 1 t2k +1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/nuebar_x_nuebar/root/'*'.root 1 t2k +1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/numu_x_nue/root/'*'.root 1 t2k +1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/numubar_x_nuebar/root/'*'.root 1 t2k +1

    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/numu_x_numu/root/'*'.root 2 t2k +1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/numubar_x_numubar/root/'*'.root 2 t2k +1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/nue_x_nue/root/'*'.root 2 t2k +1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/nuebar_x_nuebar/root/'*'.root 2 t2k +1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/numu_x_nue/root/'*'.root 2 t2k +1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/numubar_x_nuebar/root/'*'.root 2 t2k +1

    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/numu_x_numu/root/'*'.root 0 t2k -1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/numubar_x_numubar/root/'*'.root 0 t2k -1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/nue_x_nue/root/'*'.root 0 t2k -1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/nuebar_x_nuebar/root/'*'.root 0 t2k -1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/numu_x_nue/root/'*'.root 0 t2k -1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/numubar_x_nuebar/root/'*'.root 0 t2k -1

    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/numu_x_numu/root/'*'.root 1 t2k -1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/numubar_x_numubar/root/'*'.root 1 t2k -1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/nue_x_nue/root/'*'.root 1 t2k -1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/nuebar_x_nuebar/root/'*'.root 1 t2k -1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/numu_x_nue/root/'*'.root 1 t2k -1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/numubar_x_nuebar/root/'*'.root 1 t2k -1

    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/numu_x_numu/root/'*'.root 2 t2k -1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/numubar_x_numubar/root/'*'.root 2 t2k -1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/nue_x_nue/root/'*'.root 2 t2k -1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/nuebar_x_nuebar/root/'*'.root 2 t2k -1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/numu_x_nue/root/'*'.root 2 t2k -1
    ./genweights_FSI.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/numubar_x_nuebar/root/'*'.root 2 t2k -1

else

./genweights_FSI.pl /disk2/usr5/xiaoyue/skmc/test_fc/root/jan14sk4_skdetsim13p90_neut532.reduc'*'.root 0 sk

./genweights_FSI.pl /disk2/usr5/xiaoyue/skmc/test_fc/root/jan14sk4_skdetsim13p90_neut532.reduc'*'.root 1 sk

./genweights_FSI.pl /disk2/usr5/xiaoyue/skmc/test_fc/root/jan14sk4_skdetsim13p90_neut532.reduc'*'.root 2 sk

endif
