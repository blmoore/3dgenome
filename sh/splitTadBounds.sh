## h1 split TAD boundaries
python ../py/binAroundBed.py ../data/bedfiles/h1_both_tads.bed 500000 > ../data/nonrepo/inbed/hbb.bins
python ../py/binAroundBed.py ../data/bedfiles/h1_ctcf_tads.bed 500000 > ../data/nonrepo/inbed/hcb.bins
python ../py/binAroundBed.py ../data/bedfiles/h1_yy1_tads.bed 500000 > ../data/nonrepo/inbed/hyb.bins
python ../py/binAroundBed.py ../data/bedfiles/h1_none_tads.bed 500000 > ../data/nonrepo/inbed/hnb.bins

# gm12878
python ../py/binAroundBed.py ../data/bedfiles/gm_both_tads.bed 500000 > ../data/nonrepo/inbed/gbb.bins
python ../py/binAroundBed.py ../data/bedfiles/gm_ctcf_tads.bed 500000 > ../data/nonrepo/inbed/gcb.bins
python ../py/binAroundBed.py ../data/bedfiles/gm_yy1_tads.bed 500000 > ../data/nonrepo/inbed/gyb.bins
python ../py/binAroundBed.py ../data/bedfiles/gm_none_tads.bed 500000 > ../data/nonrepo/inbed/gnb.bins

## k562
python ../py/binAroundBed.py ../data/bedfiles/k5_both_tads.bed 500000 > ../data/nonrepo/inbed/kbb.bins
python ../py/binAroundBed.py ../data/bedfiles/k5_ctcf_tads.bed 500000 > ../data/nonrepo/inbed/kcb.bins
python ../py/binAroundBed.py ../data/bedfiles/k5_yy1_tads.bed 500000 > ../data/nonrepo/inbed/kyb.bins
python ../py/binAroundBed.py ../data/bedfiles/k5_none_tads.bed 500000 > ../data/nonrepo/inbed/knb.bins
