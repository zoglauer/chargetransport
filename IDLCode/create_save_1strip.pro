pro create_save_1strip,num,f_idx,cal_fe,cal_fh,ps=ps;,ztot=ztot,xtot=xtot,ytot=ytot,emax=emax
nfit=3
no_idlsave=1
if no_idlsave eq 1 then begin
  s_filename='NCT14_D'+num+'_1Pixel_z38.idlsave'
  s1_filename='idlsave/NCT14_D'+num+'_gridZ_'+f_idx+'.idlsave'
  ntot=lonarr(3)
  ;if not keyword_set(emax) then energy=59.5 else energy=emax
  
  ntot[0]=10l & ntot[1]=10l & ntot[2]=38l
  ;ntot[0]=10l & ntot[1]=10l & ntot[2]=300l
  
  t_evt=ntot[0]*ntot[1]*ntot[2]
  x=fltarr(t_evt)
  y=fltarr(t_evt)
  z=fltarr(t_evt)
  e=fltarr(t_evt)
  ;print,'Total event:',t_evt
  if ntot[0] eq 1 then dx=0.0 else $ 
  dx=(0.49-0.315)/float(ntot[0]-1)
  if ntot[1] eq 1 then dy=0.0 else $
  dy=(0.49-0.315)/float(ntot[1]-1)
  if ntot[2] eq 1 then dz=0.0 else $
  dz=(1.49-0.01)/float(ntot[2]-1)
  ;dz=(1.49-0.01)/(38-1)
  print,'dx,dy,dz',dx,dy,dz
  
  cont=0
  ;energy=60.
  energy=100.
  for k=0,ntot[2]-1l do begin
    for j=0,ntot[1]-1l do begin
      for i=0,ntot[0]-1l do begin 
      ;for k=0,37 do begin
        if ntot[0] eq 1 then x[cont]=0.4 else x[cont] =0.315+i*dx
        if ntot[1] eq 1 then y[cont]=0.4 else y[cont] =0.315+j*dy
	;x[cont] =0.4
	;y[cont] =0.4
        if ntot[2] eq 1 then z[cont]=0.75 else z[cont] =0.01+k*dz
        e[cont] =energy 
        cont=cont+1 
      endfor
    endfor
  endfor
  qtrans_new,z,x,y,e,q,num,fac_e1=cal_fe,fac_h1=cal_fh
  ;s_filename='strip1_signal.idlsave'
  ;save,q,ntot,filename=s_filename
  save,q,filename=s_filename
  print,'Qtran ... DONE! Now caluclate the time and CTD'
  t=findgen(800)*0.5
  t_info=fltarr(3,t_evt)
  for i=0,t_evt-1 do begin
    p_time=cal_time(t,q.ano[i,*],0,0.)
    n_time=cal_time(t,q.cat[i,*],0,0.)
    t_info[0,i]=p_time-n_time
    t_info[2,i]=n_time
    t_info[1,i]=q.z0[i]
    if i mod 500 eq 0 then print,i
  endfor
  
  pos=where(abs(t_info[0,*]) lt 500)
  ctd0=t_info[0,pos]
  depth0=t_info[1,pos]
  fit_vt=poly_fit(ctd0,depth0,nfit,yfit,err)
  fit_vd=poly_fit(depth0,ctd0,nfit,yfit,err)
  print,'Save idlsave file, filename: NCT14_D1_gridZ_'+f_idx+'.idlsave'
  ;printf,10,f_idx+'  '+strcompress(cal_fe)+'  '+strcompress(cal_fh)+'  '+$
  ;                     strcompress(fit_v[0])+'  '+strcompress(fit_v[1])+'  '+$
;		       strcompress(fit_v[2])+'  '+strcompress(fit_v[3])
  save,t_evt,t_info,fit_vd,fit_vt,filename=s1_filename 
    ;filename='gridz_idlsave/NCT09_D2_gridZ_'+f_idx+'.idlsave'
    ; filename='~/Desktop/Berkeley/QTransSim/idlsave/NCT14_D'+num+'_gridZ_'+f_idx+'.idlsave'
  ;print,'========'
  ;save,ctd0,depth0,filename='dcal_try_ctdz_'+f_idx+'_1.idlsave'
endif else begin
  restore,'dcal_try_ctdz_'+f_idx+'.idlsave'
endelse

  fit_vt=poly_fit(ctd0,depth0,nfit,yfit,err)
  fit_vd=poly_fit(depth0,ctd0,nfit,yfit,err)
  print,'fitting parameter: (ctd2z)'
  for i=0,nfit do print,fit_vt[i]
  print,'fitting parameter: (z2ctd)'
  ;for i=0,nfit do print,fit_vd[i]
  FORMAT2='(%"%11.5f %11.5f %11.5f %11.5f")'
  print,format=FORMAT2,fit_vd[0],fit_vd[1],fit_vd[2],fit_vd[3]
  ;printer=1
  ;if(printer eq 0) then begin
  ;set_plot,'X'
  ;DEVICE,PSEUDO_COLOR=8
  if keyword_set(ps) then begin
  !P.Multi=[0,1,1,0,0]
  !P.CharSize=1.5
  set_plot,'ps'
  device,filename='z_ctd_curve.ps'
  endif
  ;Window,1,XSize=900,YSize=300,XPos=0,YPos=200
  ;setcolors, NAMES=cnames, VALUES=cindx
  ;endif else begin
  ;keywords = PSConfig()
  ;thisDevice = !D.Name
  ;Set_Plot, 'PS'
  ;Device, _EXTRA=keywords
  ;endelse 
  ;!P.Multi=[0,2,1,0,0]
  ;!P.CharSize=1.5
  ;Window,1,XSize=900,YSize=300,XPos=0,YPos=200
  if keyword_set(ps) then begin
  plot,ctd0,depth0,psym=1,xtitle='CTD (ns) D7',ytitle='depth (cm)'
  ctdx=findgen(100)*4.-200
  zy=fltarr(100)
  for i=0,nfit do begin
    zy=zy+fit_vt[i]*ctdx^i
  endfor
  oplot,ctdx,zy,color=100,thick=2
  ;oplot,[-165.245,111.822],[1.5,0]
   
  ;plot,depth0,ctd0,psym=1,ytitle='CTD (ns)',xtitle='depth (cm)'
  ;depx=findgen(100)*0.15
  ;zt=fltarr(100)
  ;for i=0,nfit do zt=zt+fit_vd[i]*depx^i
  ;oplot,depx,zt,color=100,thick=2
  endif

  if keyword_set(ps) then begin
  device,/close
  set_plot,'X'
  endif
  ;oplot,[1.5,0],[-165.245,111.822]
  ;if(printer eq 1) then begin
  ;Device, /Close_File
  ;Set_Plot, thisDevice
  ;endif   
end
