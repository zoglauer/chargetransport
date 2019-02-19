;.com define_common_blks.pro
;define_common_blks,/cgs
;.com field3d.pro
;.com vectgrid.pro
pro make_vfield,num
restore,'idlsave/NCT14_D'+num+'_efield.idlsave'
eugrid=temporary(eugrid[*,*,5:107-5-1]);kair:kmax-kair-1
egrid=vectgrid(eugrid)
save,filename='idlsave/NCT14_D'+num+'_v_efield.idlsave',egrid

;restore,'zk_wfieldx.idlsave'
;eugrid=temporary(eugrid[*,*,5:107-5-1])
;wgrid=vectgrid(eugrid)
;save,filename='vwfield.idlsave',wgrid
end
