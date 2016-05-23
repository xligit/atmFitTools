#!/bin/tcsh

if ($#argv == 0) then
    ./genweights_throws.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/numu_x_numu/root/ +1 0 numu_x_numu t2k
    ./genweights_throws.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/numubar_x_numubar/root/ +1 0 numubar_x_numubar t2k
    ./genweights_throws.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/nue_x_nue/root/ +1 0 nue_x_nue t2k
    ./genweights_throws.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/nuebar_x_nuebar/root/ +1 0 nuebar_x_nuebar t2k
    ./genweights_throws.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/numu_x_nue/root/ +1 1 numu_x_nue t2k
    ./genweights_throws.pl /disk/usr4/okumura/t2k_14a_root/nu-mode/numubar_x_nuebar/root/ +1 1 numubar_x_nuebar t2k


    ./genweights_throws.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/numu_x_numu/root/ -1 0 numu_x_numu t2k
    ./genweights_throws.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/numubar_x_numubar/root/ -1 0 numubar_x_numubar t2k
    ./genweights_throws.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/nue_x_nue/root/ -1 0 nue_x_nue t2k
    ./genweights_throws.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/nuebar_x_nuebar/root/ -1 0 nuebar_x_nuebar t2k
    ./genweights_throws.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/numu_x_nue/root/ -1 1 numu_x_nue t2k
    ./genweights_throws.pl /disk/usr4/okumura/t2k_14a_root/antinu-mode/numubar_x_nuebar/root/ -1 1 numubar_x_nuebar t2k

else
    ./genweights_throws.pl /disk2/usr5/xiaoyue/skmc/test_fc/root/ +1 0 sk sk
endif
