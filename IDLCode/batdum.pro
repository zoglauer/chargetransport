pro batdum,num

MAXI=10000

define_common_blks,num,/cgs
print,"Fields round 1"
field3d,num,'deleteme1.idlsave',maxits=MAXI
print,"Fields round 2"
field3d,num,'deleteme2.idlsave',maxits=MAXI
print,"Fields round 3"
field3d,num,'idlsave/NCT14_D'+num+'_efield.idlsave',/conduct,maxits=MAXI,tol=1e-12
print,"Fields round 4"
field3d,num,'idlsave/NCT14_D'+num+'_wfield.idlsave',/wfield,/conduct,maxits=MAXI,tol=1e-12

end
