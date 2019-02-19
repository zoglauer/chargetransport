.com define_common_blks
.com create_save_1strip
.com qtrans_new
.com field3d.pro
.com vectgrid.pro

Det_ID='10'

define_common_blks,Det_ID,/cgs
batdum,Det_ID
make_vfield,Det_ID
create_save_1strip,Det_ID,'standard',1.0,1.0,/ps
print,'DONE!!'

