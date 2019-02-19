pro batdum,num
define_common_blks,num,/cgs
field3d,num,maxits=10000
field3d,num,maxits=10000
field3d,num,'idlsave/NCT14_D'+num+'_efield.idlsave',/conduct,maxits=10000,tol=1e-12
field3d,num,'idlsave/NCT14_D'+num+'_wfield.idlsave',/wfield,/conduct,maxits=10000,tol=1e-12
end
