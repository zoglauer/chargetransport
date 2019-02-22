;*****************************************************************************
; FILE_NAME: QTRANS
;
; PURPOSE: Simulates the charge transport and signal induction on several 
;          strips of a Ge detector
;          
; INCLUDED PROCEDURES/FUNCTIONS:
;         INITIAL_POS
;         QMOVE()
;         DRIFT_DIR()
;         QTRANS
;
; FULL DOCUMENTATION: Search for ';+****' below          
;*****************************************************************************

;ntostr added by SEB, can't find Susan's program
function ntostr,ndum
cdum=strtrim(ndum,2)
return,cdum
end

;*****************************************************************************
; PROCEDURE NAME: INITIAL_POS
;
; PURPOSE: To sort out the many ways in which simulated event 
;           parameters can be chosen
;
; SYNTAX: initial_pos, zin, xin, yin, ein, qin, num, atten, random=random, 
;         event_atten=event_atten, expdist=expdist
; 
; WAYS TO DEFINE INITIAL PARAMS:
;     Every input event in the simulation must be defined by a 3D initial 
;   interaction position (Xin, Yin, Zin) and an energy (ein or atten). There
;   are several different ways to define these using input params.
;   INPUTS are described in the QTRANS main documentation below.
;    
; OUTPUTS: q - charge associated with ein
;          num - number of total events
;          atten - chosen attenuation
;*****************************************************************************
PRO initial_pos, zin, xin, yin, ein, qin, num, atten, random=random, event_atten=event_atten, expdist=expdist

COMMON FLDGRID
COMMON PIXEL
COMMON DETECTOR
QE=1.60217733d-19 ; elementary charge

;Determine number of objects
expdist=keyword_set(expdist)? expdist:0
num = ( expdist GT 1 )? expdist: $
   (n_elements(zin) > (n_elements(xin)*(1-keyword_set(random)))) $
  > (n_elements(yin)*(1-keyword_set(random)))

IF (num EQ 0) THEN $
  IF (n_elements(ein) GT 1) THEN num = n_elements(ein) $
  ELSE num=2000
                  ;default 2000 objects, 122 keV exp dist

;---------------------------------
;DETERMINE X,Y,Z INITIAL POSITIONS
;---------------------------------

;choose random grid of x,y from midgap_left to midgap_rt of center strip 
xstr='' & ystr='' & xset=0 & yset=0
IF keyword_set(random) OR (n_elements(xin) EQ 0) THEN BEGIN 
  xin=(randomu(seed,num)*istrip +(imax/2-istrip/2))*DX 
  xstr=' X' & xset=1
ENDIF ELSE IF n_elements(xin) NE num THEN message, $
  'X0 must have same number of elements as Z0'

IF keyword_set(random) OR (n_elements(yin) EQ 0) then BEGIN
  yin=(randomu(seed,num)*istrip +(jmax/2-istrip/2))*DY
  ystr='Y ' & yset=1
ENDIF ELSE IF n_elements(yin) NE num THEN message, $
  'Y0 must have same number of elements as Z0'

both=(xset+yset EQ 2)? ',':' '
IF (xset+yset NE 0) THEN print,'Using random'+xstr+both+ystr+'start positions'

;--------------------------
;DETERMINE INITIAL CHARGE
;--------------------------

ne0=n_elements(ein)
CASE 1 OF 
  
  (ne0 EQ 1): BEGIN
    energy_type = ntostr(ein[0])+' keV'
    qin=replicate(QE*ein[0]*1e3/2.98d, num)
    atten=152
    ;atten=get_atten(ein[0])
  END

  (ne0 EQ 0): BEGIN
    energy_type = 'by default 122 keV'
    qin=replicate(QE*122.0e3/2.98d, num)
    ein=[122.0]
    atten=152d
    ;atten=get_atten(122)
  END

  ELSE: BEGIN
    IF (ne0 NE num) THEN message,'E0 must have same number of elements as Z0'
    energy_type = 'Input Energies'
    qin=QE*ein*1e3/2.98d
    atten=0
  ENDELSE

ENDCASE

atten=keyword_set(event_atten)? event_atten: atten
print,'Energies are '+energy_type

;---------------------------
;FINALLY GET Z DISTRIBUTION
;---------------------------
IF keyword_set(expdist) OR (n_elements(zin) EQ 0) THEN BEGIN
  MN=0.0 ;cm
  MX=1.5 ;cm
  zin=randomu_exp(num,-1.*atten,MN,MX,/error)  ; -1.* added MEB, 10/18/05 to fix AC/DC illumination
  IF zin[0] EQ 0 THEN BEGIN
    zin=randomu(seed,num)*(MX-MN) + MN
    print,'Attenuation undefined, random Z Dist being used' 
  ENDIF ELSE print,'Exponential Z0, attenuation: ',ntostr(atten)+' cm^-1'
endif
;* note that expdist is used if z0 is NOT defined even if expdist is not set.

return
END  ; initial_pos

; Note the electric and weighting fields assume that the bias voltage is
; applied to the electrode at z=HEIGHT, and that the signal is taken from
; the ground electrode at z=0.

;*****************************************************************************
; FUNCTION NAME: DRIFT_DIR
;
; PURPOSE: Find the drift velocity and direction of holes and
;          electrons given the vector electric field.
;
; INPUTS: e - vector electric field
;
; OUTPUT: Drift Velocity (mag and dir) for each hole and electron
;
; OPTIONAL INPUT: /drift111 - set to use drift velocities for the 
;                 <111> crystal orientation. Default is <100>
;
; NOTES: Drift velocities come from Ottaviani,Canali & Quaranta, 
;         IEEE Trans. Sci., Vol NS 22, pp. 192-204, Feb 1975.
;*****************************************************************************  
FUNCTION drift_dir,e,st=st,drift111=drift111	; q is currently unused.
; Calculate the direction (and magnitude) of the charge displacement
; e[] is assumed to have 4 dimensions, with the last one representing
;  the 3 spatial directions.
  COMMON FLDGRID
  COMMON STEP	; DTSTEP,NUMITS, dtstep_tolerance,efield_tolerance
		; Only DTSTEP and efield_tolerance are used here.
  COMMON FACTOR_CAL
tmp=size(e,/dimensions) & n0=tmp[0] & n1=tmp[1] & n2=tmp[2]

;Experimental values from which to interpolate drift velocity
IF keyword_set(drift111) THEN BEGIN ;<111>
  e_e=[0,30., 50,100,200,500,1000,2000,5000]
  vd_e=[0,1.10000E+06,1.70000E+06,2.70000E+06,4.10000E+06,6.20000E+06,7.50000E+06,8.90000E+06,9.20000E+06]
  e_h=[0,10,20.,50.,100.,200.,500,1000,2000,5000,10000]
  vd_h=[0,420000.,860000.,2.10000E+06,3.20000E+06,4.30000E+06,5.70000E+06,6.90000E+06,7.90000E+06,8.90000E+06,9.10000E+06]
ENDIF ELSE BEGIN ;<100>
  e_e=[0,30.,50.,100.,200.,500.,1000.,2000.,5000.] ;[V/cm]
  vd_e=[0,1.10000E+06,1.70000E+06,2.70000E+06,4.10000E+06,9.30000E+06,1.12500E+07,1.33500E+07,1.38000E+07] ;[cm/s]
 ; changed first entry from 0 to 0.1 MEB 2006 (for log interpolation)
  e_h=[0.1,10,20.,50,100,200,500.,1000,2000,5000,10000.] ;[V/cm]  
  vd_h=[0,420000.,860000.,2.10000E+06,3.20000E+06,4.83700E+06,6.41250E+06,7.76250E+06,8.88750E+06,1.00130E+07,1.02380E+07] ;[cm/s]
endelse

;low efield approximation:
;mu0_e=36666.7
;;mu0_e = 36000.0 ;electron mobility (77k) cm^2/Vs, KNOLL p357 (3rd)
;;mu0_h=42000.0 ;hole mobility (77k) cm^2/Vs, KNOLL p357 (3rd)
;mu0_h=43000.0

;;Emin_e=30.  ;'low field' is E < Emin_e 
;;Emin_h=10.  ;'low field' is E < Emin_h  ;THIS IS OFF FOR NOW

;Emin_e=0.  ;'low field' is E < Emin_e 
;Emin_h=0.  ;'low field' is E < Emin_h


eabs=dblarr(n0,n1,n2,/nozero) ; magnitude of e[]
eabs[*,*,0]=sqrt(total(e^2,3)) ;[nx*ny,2,3]
eabs[*,*,1]=eabs[*,*,0]
eabs[*,*,2]=eabs[*,*,0]
vd=dblarr(n0,n1,n2)

;get electron velocities
;vd[*,0,*]= [(eabs[*,0,*] lt Emin_e) * (mu0_e*eabs[*,0,*]) + $
;              (eabs[*,0,*] ge Emin_e) * $
;              interpol(vd_e,e_e,eabs[*,0,*])] * $
;  (-1)*e[*,0,*]/(eabs[*,0,*] > efield_tolerance) * $
;  (eabs[*,0,*] gt efield_tolerance)

;vd[*,0,*]= interpol(vd_e,e_e,eabs[*,0,*]) * $  ; spline added MEB 2006
vd[*,0,*]= interpol(vd_e,e_e,eabs[*,0,*], /spline) * $
  (-1)*e[*,0,*]/(eabs[*,0,*] > efield_tolerance) * $
  (eabs[*,0,*] gt efield_tolerance)

electron_scale_factor=fac_e;130./85. ;to match the experimental range in tdif
vd[*,0,*]=electron_scale_factor * temporary(vd[*,0,*])

;get hole velocities
;vd[*,1,*]= [(eabs[*,1,*] lt Emin_h) * (mu0_h*eabs[*,1,*]) + $
;              (eabs[*,1,*] ge Emin_h) * $
;              interpol(vd_h,e_h,eabs[*,1,*])] * $
;  e[*,1,*]/(eabs[*,1,*] > efield_tolerance) * $
;  (eabs[*,1,*] gt efield_tolerance)

;vd[*,1,*]=  interpol(vd_h,alog(e_h),alog(eabs[*,1,*])) * $
vd[*,1,*]=  interpol(vd_h,e_h,eabs[*,1,*], /spline) * $  ; alog added here MEB 2006
  e[*,1,*]/(eabs[*,1,*] > efield_tolerance) * $
  (eabs[*,1,*] gt efield_tolerance)
hole_scale_factor=fac_h;135./110. ;to match the experimental range in tdif
vd[*,1,*]=hole_scale_factor * temporary(vd[*,1,*])

return, vd*DTSTEP 
END ; drift_dir()


;*****************************************************************************
; FUNCTION NAME: QMOVE
;
; PURPOSE: This is the main function for charge transport
; 
; Inputs and keywords described in QTRANS documentation below
;*****************************************************************************
FUNCTION qmove,x0,y0,z0,q0,egrid,wgrid,$
               tau_surf=tau_surf,notrap=notrap,save_path=save_path,$
               multich=multich,rstep=rstep,numsites=numsites,$
               zipfilename=zipfilename,sim=sim, $
               drift111=drift111

COMMON PIXEL
COMMON FLDGRID
COMMON DETECTOR
COMMON STEP
COMMON DIR
;HERE IS WHERE I RESTORE CATHODE FIELD - OFF
;print,'Restoring cathode field'
;restore,'~/current_fields/wfield_cath_c16.sav'
;print,'done'
tmp=size(x0,/dimensions) & nx=tmp[0]>1 ;ny has been removed because there is now only 1 dimension
en2q = 1.60217733d-19 * 1d3 / 2.98d   ;(C/keV) for a 1MeV event
;-------------
;ARRAY SET UP
;-------------
qano=dblarr(nx,numits+3) ; charge induced at the anode in double precision
qcat=qano	; ditto for the cathode

x=dblarr(nx,2,3) ;this array will store the current absolute position at each iteration

;-----------------
;DEFINE CONSTANTS
;-----------------
;add random step variables
;Xrms = sqrt( (8.61d-8 * 72.)/511. ) * 2.99d10 * DTSTEP
;this is the random thermal step per dimension

kT_e = 0.0253 * (85./293.) ;V, value of kT/e for T=72K
;diffusion coefficient for holes:
d_h = 4.2e4 *kT_e ;cm^2/s
;diffusion coefficient for e-:
d_e = 3.6e4 *kT_e ;cm^2/s
;sigma = sqrt(2*D*DTSTEP)

SEP = 0.200 ;cm. This is the separation between strips if 2 are 'active'.
 
IGAP=(ISTRIP-ICNTCT)+1  ;info about the grid, ISTRIP, ICNTCT etc stored in DEFINE_COMMON_BLKS.PRO 
NSTRIP=fix((IMAX-1)/ISTRIP)
D111=keyword_set(drift111) ;not set is default=drift100 set

;------------------------------
;SET INITIAL X, Y, Z POSITIONS
;------------------------------
x[*,0,0]=x0 & x[*,1,0]=x0
x[*,0,1]=y0 & x[*,1,1]=y0
x[*,0,2]=z0 & x[*,1,2]=z0

;-----------------------
;RECORD INFO FOR TESTING
;-----------------------

if keyword_set(save_path) then begin ;save full path for testing
  follow=where(x[*,0,2] gt 0,nfollow) ;whose path is needed? 
  paths=dblarr(nfollow,numits+1,2,3) ;store x,y,z pos at each iter
  paths[*,0,*,*]=x[follow,*,*] 
  ef=dblarr(nfollow,numits+1,2,3) ;store electric field Ex, Ey, Ez at each iter
  wf=ef ;anode weighting field
  wfc=ef ;cathode weighting field
  tx1=ef
  
  IF keyword_set(multich) THEN wf2=wf ;second anode weighting field
    print,'Following ',ntostr(nfollow),' charges...'
endif else begin
  print,'no charges followed...'
  ngs=0
endelse
IF D111 THEN print,'Using <111> direction'
;------------------
;DEFINE FIELD GRIDS
;------------------
;These will hold the interpolated field values at each iteration 
efield=dblarr(nx,2,3,/nozero)
wfield=dblarr(nx,2,/nozero)
wcfield=dblarr(nx,2,/nozero) ;weighting field at cathode

if keyword_set(multich) then begin
  wfield2=wfield ;store a second wfield array if multiple channels are active
  wcfield2=wcfield
  qano2=qano ;charge array for second anode channel
  qcat2=qcat
endif

;-------------------
;DEFINE GAP GEOMETRY
;-------------------
;find the position of the right side of each strip on the anode, etc... 
ano_right=((ICNTCT/2+indgen(NSTRIP+1)*ISTRIP) < IMAX)*DX
ano_left=[0,ano_right(0:NSTRIP-1)+IGAP*DX]
cat_right=((ICNTCT/2+indgen(NSTRIP+1)*ISTRIP) < JMAX)*DY
cat_left=[0,ano_right(0:NSTRIP-1)+IGAP*DY]

gap_grid=intarr(IMAX) ;a 1D grid marking the gap at an electrode plane, along y at the cathode, along x at the anode. 
for i=0,NSTRIP-1 do gap_grid[ano_right(i)/DX : ano_left(i+1)/DX] = 1.0
scalexy=replicate(1000.0,nx)  ;same as X, larger than txstep will ever be (see use below)

;*if dx is ever not equal to dy, add a second array gap_gridy

;-------------------------------------
;BEGIN MAIN CHARGE TRANSPORT ITERATION
;-------------------------------------
;openw,mcheck,'check_step_081021.out',/get_lun
k=0L	; Iteration counter, retained for `safety', to avoid infinite looping.
cath=dblarr(nx,2) ;use because two transforms must be made to get to wcfield. 
repeat begin
  k=k+1L	; Increment loop counter.
  dxabs=0b	; Free unused memory from previous iteration.

  ;-----------------------------------------------------------------------
  ;Find the electric and weighting field vectors at the current positions.
  ;-----------------------------------------------------------------------

  for dir=0,2 do begin
    efield[*,*,dir] = interpolate(egrid[*,*,*,dir],	x[*,*,0]/DX,$
							x[*,*,1]/DY,$
							x[*,*,2]/DZA )
  endfor;dir=0,2						
  wfield[*,*] = interpolate(wgrid[*,*,*],x[*,*,0]/DX,$
					 x[*,*,1]/DY,$
					 x[*,*,2]/DZA )
  cath[*,*] = interpolate(wgrid[*,*,*],	x[*,*,1]/DY,$
                                        (JMAX-1) - x[*,*,0]/DX,$
			 		(KMAX-2*KAIR-1) - x[*,*,2]/DZA )
;    wcfield[*,*,dir] = interpolate(wgrid_cath[*,*,*,dir],x[*,*,0]/DX,$
;                                                         x[*,*,1]/DY,$
;                                                         x[*,*,2]/DZA )
;;  endfor ; dir=0,2
  ;transform vector components of anode weighting field to cathode WF
  wcfield=cath
;;  wcfield[*,*]=-1*cath[*,*,1]
;;  wcfield[*,*]=cath[*,*,0]
;;  wcfield[*,*]=-1*cath[*,*,2]

;  if k eq 1 then save,wcfield,filename='~/qtrans/wcfield0.sav'
;Consider cross-talk effect
  if keyword_set(multich) then begin ;store for anode only
      wfield2[*,*] = interpolate(wgrid[*,*,*],((x[*,*,0] - sep)/DX),$
						x[*,*,1]/DY,$
						x[*,*,2]/DZA )
      wcfield2[*,*] = interpolate(wgrid[*,*,*],	(x[*,*,1] - sep)/DY,$
                                             (JMAX-1) - x[*,*,0]/DX,$
			              (KMAX-2*KAIR-1) - x[*,*,2]/DZA )
    ;our second channel is to the right in X (on the anode side), 
    ;hence ITs weighting field is the same as wgrid for a position 
    ;SEP microns to the left in X.
    ;;within_bound = ((x[*,*,0]-SEP) ge 0) and ( (x[*,*,0]-SEP) le DX*(IMAX-1))
    ;;wfield2[*,*,0]=wfield2[*,*,0]*within_bound
    ;;within_bound=0b             ; Free memory.
  endif

  ;if k eq 1 then wcfield0[*,*,*]=wcfield[*,*,0,*]

  ;-----------------------
  ;CHECK FOR OUT OF BOUNDS
  ;-----------------------
  ; Null the fields in the x direction where the position is out of bound.
  ; Note that the keyword `missing' for interpolate() cannot be used,
  ;  because we want to null each field only in one direction.
  within_bound = (x[*,*,0] ge 0) and (x[*,*,0] le DX*(IMAX-1))
  efield[*,*,0]=efield[*,*,0]*within_bound
  ;;wfield[*,*,0]=wfield[*,*,0]*within_bound
  ;;wcfield[*,*,0]=wcfield[*,*,0]*within_bound
  within_bound=0b	; Free memory.
  ; Ditto, in the y direction.
  within_bound = (x[*,*,1] ge 0) and (x[*,*,1] le DY*(JMAX-1))
  efield[*,*,1]=efield[*,*,1]*within_bound
  ;;wfield[*,*,1]=wfield[*,*,1]*within_bound
  ;;wcfield[*,*,1]=wcfield[*,*,1]*within_bound

  ;;if keyword_set(multich) then wfield2[*,*,1]=wfield2[*,*,1]*within_bound
  ;;within_bound=0b	; Free memory.
  ; In the z direction, null the electric field where the position is either
  ; out of bound or on a surface; the weighting field is unchanged.
  efield[*,*,2]=efield[*,*,2]*(	 (x[*,*,2] gt 0  )$
				   and	(x[*,*,2] lt HEIGHT) )

  ;-------------------------------------------------
  ; Find the movement induced by the electric field
  ;-------------------------------------------------
  txstep=drift_dir(efield,drift111=d111)

  ;---------------------------
  ;ADD THERMAL MOTION (IF SET)
  ;---------------------------
  ; add random thermal step if keyword /rstep is set
  if keyword_set(rstep) then for dir = 0, 2 do begin
    txstep[*,0,dir] = txstep[*,0,dir] + randomu(seed, nx, /normal)* $
      sqrt(2*D_e*DTSTEP)
    txstep[*,1,dir] = txstep[*,1,dir] + randomu(seed, nx, /normal)* $
      sqrt(2*D_h*DTSTEP)
  ENDFOR

  ;----------------------------
  ;CHECK FOR BOUNDARY OVERSHO0T
  ;----------------------------
  ; Check that the boundary wasn't overshot, and correct if it is.
  ; testing change
  tx=x+txstep
  too_high = (tx[*,*,2] ge HEIGHT) ;e- hit anode
  too_low  = (tx[*,*,2] le 0.0d  ) ;hole hit cathode
 
  ;tx=0b

  ;;reset scalexy to arbitrary high value
  ;scalexy[*]=1000  
  ;gapx = gap_grid(floor(x[*,0,0]/DX) < IMAX-1)*too_high[*,0] 
  ;for i=0,NSTRIP do $
  ;  scalexy[*] = ((scalexy[*] < abs(ano_left(i)-x[*,0,0])) < $
  ;               abs(x[*,0,0]-ano_right(i)))
  ;scalexy[*] = scalexy[*] < abs(txstep[*,0,0])
  ;sign=((txstep[*,0,0] lt 0)*(-2) + 1)

  ;txstep[*,0,0] = (1b-too_high[*,0]) * txstep[*,0,0] + $
  ;  gapx* sign * scalexy

  ;reset scalexy to arbitrary high value
  ;scalexy[*]=1000
  ;gapy = gap_grid(floor(x[*,0,1]/DY) < JMAX-1)*too_low[*,0] 
  ;for i=0,NSTRIP do $
  ;  scalexy[*] = ((scalexy[*] < abs(cat_left(i)-x[*,0,1])) < $
  ;               abs(x[*,0,1]-cat_right(i)))
  ;scalexy[*] = scalexy[*] < abs(txstep[*,0,1])
  ;sign=((txstep[*,0,1] lt 0)*(-2) + 1)

  ;txstep[*,0,1] = (1b-too_low[*,0]) * txstep[*,0,1] + $
  ;  gapy* sign * scalexy 

  ;sign=0b
  ;scale42high=(HEIGHT-x[*,*,2])/(abs(txstep[*,*,2])>dxstep_tolerance)
  ;scale42low =(x[*,*,2]-0.0d  )/(abs(txstep[*,*,2])>dxstep_tolerance)
  
  ;txstep[*,*,2]$
  ;  = (1b-too_low-too_high)	* txstep[*,*,2]$
  ;  +too_high			* txstep[*,*,2] * scale42high $
  ;  +too_low			* txstep[*,*,2] * scale42low
  ;scale42high=0 & scale42low=0	; Free unused memory.


  ;qtrap=q0                      ;initially no charge is trapped

  ;---------------------------------
  ;CALCULATE INDUCED CHARGE FOR STEP
  ;---------------------------------
  qano(*,k)=q0*(wfield(*,0)-wfield(*,1))
  qcat(*,k)=q0*(wcfield(*,0)-wcfield(*,1))
  ;printf,mcheck,strcompress(x[0,0,2])+' '+strcompress(x[0,1,2])+' '+$
  ;              strcompress(qano[0,k]/en2q)+' '+$
  ;              strcompress(qcat[0,k]/en2q)+' '+strcompress(wfield[0,0])+$
  ;		strcompress(wfield[0,1])
  ;get the charge induced on the second anode channel if /multich is set
  
  ;Consider cross-talk effect
  if keyword_set(multich) then begin
    qano2(*,k)=q0*(wfield2(*,0)-wfield2(*,1))
    qcat2(*,k)=q0*(wcfield2(*,0)-wcfield2(*,1))
  endif
  ;save each iteration (i.e. the path taken by each charge) if keyword is set
  ;;if keyword_set(save_path) then begin
  ;;  paths[*,k,*,*]=x[follow,*,*]
  ;;  ef[*,k,*,*]=efield[follow,*,*]
  ;;  wfc[*,k,*,*]=wcfield[follow,*,*]
  ;;  wf[*,k,*,*]=wfield[follow,*,*]
  ;;  tx1[*,k,*,*]=txstep[follow,*,*]
  ;;  if keyword_set(multich) then wf2[*,k,*,*]=wfield2[follow,*,*]
    
  ;;endif

  ;-------------------------------
  ;MOVE CHARGES ALONG EFIELD LINES
  ;-------------------------------
; Now move the charges and ...
  blas_axpy,x,1.0d,txstep
  x[*,*,2]=x[*,*,2]*(1b-too_low-too_high) + HEIGHT*too_high ; + 0.0*too_low

 
  ;----------------------------------------
  ;TEST FOR TOTAL COLLECTION OF ALL CHARGES
  ;----------------------------------------
  dxabs=sqrt(total(txstep^2,3))
  txstep=0b                    
;endrep until ((total(dxabs ge dxstep_tolerance) eq 0) or (k ge NUMITS))
endrep until (k ge NUMITS) ;SEB 5/25/05
;free_lun,mcheck
;------------------------------------
;END MAIN CHARGE TRANSPORT ITERATION
;------------------------------------

if (k gt NUMITS)$
 then message,'Number of iterations exceeds expected value!  Forced exit.',$
	/continue
 
;=save paths for testing=
;;if keyword_set(save_path) then BEGIN

;;  path1=create_struct('x',fltarr(k,2), $
;;                      'y',fltarr(k,2), $
;;                      'z',fltarr(k,2), $
;;                      'ef',fltarr(k,2,3), $
;;                      'wf',fltarr(k,2,3), $
;;                      'wfc',fltarr(k,2,3), $
;;                      'q',fltarr(k,2),$
;;                      'tx',fltarr(k,2,3))
;;  path=replicate(path1, nfollow)
;;  FOR m=0,nfollow-1 DO BEGIN
;;    path[m].x=paths(follow[m],indgen(k),*,0)
;;    path[m].y=paths(follow[m],indgen(k),*,1)
;;    path[m].z=paths(follow[m],indgen(k),*,2)
;;    path[m].ef=temporary(ef(follow[m],indgen(k),*,*))
;;    path[m].wf=temporary(wf(follow[m],indgen(k),*,*))
;;    path[m].wfc=temporary(wfc(follow[m],indgen(k),*,*))
;;    path[m].tx=temporary(tx1(follow[m],indgen(k),*,*))
;;    path[m].q[*,0]=temporary(qano[follow[m],0:k-1])/en2q
;;    path[m].q[*,1]=temporary(qcat[follow[m],0:k-1])/en2q
;;  ENDFOR

;;  IF (strpos(save_path,'.fits'))[0] NE -1 THEN BEGIN
;;    mwrfits,path,outfile_dir+save_path,/create 
;;    pathstyle='fits format'
;;  ENDIF ELSE BEGIN
;;    save,path,filename=outfile_dir+save_path
;;    pathstyle='IDL Save format'
;;  ENDELSE
  
;;  print,'Followed charges saved in '+pathstyle+' to: ',outfile_dir+save_path
  
;;endif
;==over==
;save,qano,qcat,filename='~/qtrans/outfiles/qanoqcat_testing.sav'

;-------------------------------------
;COMBINE WAVEFORMS INTO SINGLE EVENTS
;-------------------------------------
;combine waveforms if an evt occurs at more than 1 site (i.e. is compton scattered)
if keyword_set(numsites) then begin
  combine_waveforms,qano,qcat,numsites,index_keep
  nx0=nx
  nx=n_elements(qano[*,0])
  print,'Combined waveforms ... ',nx0,' reduced to ',nx
endif

;---------------------------
;CALCULATE TIME DIFFERENCES
;---------------------------
;---For comparison purposes only...
;do time positioning stuff
;if (not keyword_set(sim)) then begin 
  ;half_timeano=fltarr(nx)
  ;half_timecat=fltarr(nx)
  ;xint=indgen(k)*DTSTEP*1e9  ;times in ns
  ;for i=0L,nx-1 do begin
  ;  qano_mx=max(qano[i,*])      ;assume min is 0.0
  ;  qcat_mn=min(qcat[i,*])      ;assume max is 0.0
;    qcat_mn=-1.*qano_mx         ;force max charge to be the same.

  ;  half_timeano[i]= interpol(xint,qano(i,0:k-1)/qano_mx,0.5)
; DTSTEP*(where(abs(qano[i,*]/qano_mx - 0.5) eq $
;                                   min(abs(qano[i,*]/qano_mx - 0.5))))[0]
    
  ;  half_timecat[i]= interpol(xint,qcat(i,0:k-1)/qcat_mn,0.5)
;DTSTEP*(where(abs(qcat[i,*]/qcat_mn - 0.5) eq $
;                                   min(abs(qcat[i,*]/qcat_mn - 0.5))))[0]
  ;endfor
  ;*******NOTE changed scales to ns********

  ;tdif = half_timeano - half_timecat
;  save,tdif,filename = '~/qtrans/geant_tdif_comparison.sav'
;  print,'saved tdif in: ','~/qtrans/geant_tdif_comparison.sav'
;endif
;---


;------------------------
;CHANGE Q UNITS TO keV!!!
;------------------------
qano = float (temporary(qano[*,0:k-1])/en2q)   ;to save storage space
qcat = float (temporary(qcat[*,0:k-1])/en2q)
IF keyword_set(multich) THEN BEGIN
  qano2 = float (temporary(qano2[*,0:k-1])/en2q)
  qcat2 = float (temporary(qcat2[*,0:k-1])/en2q)
  print,"qano is setting"
ENDIF
;=test x move by zk=
;print,'now plot!'
;Window,1,XSize=900,YSize=600,XPos=0,YPos=200
;plot,my_movex(0,*),my_movex(1,*),psym=7
;=test over=
;zipfilename=keyword_set(zipfilename)? zipfilename: 'sim_zip.zip'

;-----------------------------------------
;RETURN STRUCTURE OF TOTAL INDUCED CHARGE
;-----------------------------------------
case 1 of
  
  keyword_set(sim): BEGIN
    ;convert induced charge into simulated oscilloscope output
    sim_to_oscilloscope, qano, qcat, flags, noise_level = 0.0005
    return, flags
  end
  Keyword_set(multich): return,{ano:qano[*,0:k-1], cat:qcat[*,0:k-1],$
                                ano2:qano2[*,0:k-1],$
				cat2:qcat2[*,0:k-1],$
                                x0:x0,  $
                                y0:y0,$
                                z0:z0, $
                                e0:fltarr(nx)} 
  else: return,{ano:qano[*,0:k-1], cat:qcat[*,0:k-1],$
                x0:x0,  $
                y0:y0,$
                z0:z0, $
                e0:fltarr(nx)} 
endcase
END                             ; qmove


;+*****************************************************************************
; NAME: QTRANS
;
; PURPOSE: To simulate the charge transport and signal induction on an 
;          active anode (DC) and cathode (AC) strip in a 
;          Ge detector given input event parameters. 
;
; SYNTAX: qtrans,z0,x0,y0,e0,q,outfile=outfile,
;           /notrap, /anode_irradiate, egrid=egrid,wgrid=wgrid,
;          (through DRIFT_DIR): /drift111=drift111
;          (through INITIAL_POS): EXPDIST=NUM OR 1, /RANDOM,  
;               EVENT_ATTEN= # in cm^-1 
;
; INPUTS: 
;   z0 - input array of NUM initial Z (depth) positions (0-1.1 cm)
;         *if undefined, /EXPDIST is automatically set 	
;   x0 - input array of NUM initial X positions (0-0.8 cm)
;         *if undefined, /RANDOM is set for x0 ONLY.
;   y0 - input array of NUM initial Y positions (0-0.8 cm)
;         *if undefined, /RANDOM is set for y0 ONLY.
;   e0 - input array of NUM initial energies  (keV) 
;         *set e0 to a single energy if the same energy is 
;          desired for all positions
;        ** if undefined,  Default = 122 keV (NOTE: this does NOT set the
;           attenuation for 122 keV! It only affects the final signal
;           amplitude.)
;
;  NOTE: The total number of events NUM is defined by EXPDIST (if set and > 1),
;  the array size of z0 (if defined), x0 & y0 (if defined and /RANDOM is not
;  set), e0 (if defined and > 1 element), or the default = 2000 in that order. 
;    
; OPTIONAL INPUTS:    *KEYWORDS ALWAYS OVERRIDE INPUT ARRAYS*
;
;    /NOTRAP - Set to exclude trapping - CURRENTLY THIS IS DEFAULT!
;    /ANODE_IRRADIATE - set this to transform z0 => 1.1 - z0. If z0 is an
;                      exponential distribution, this places the majority of
;                      events near z0=1.1 cm, as if the ANODE (DC) side of the
;                      detector were irradiated. The default is Cathode (AC)
;                      irradiation which places the majority of events near
;                      z0=0.0 cm.
;     
;   (from DRIFT_DIR function)
;    /DRIFT111 - set to use drift velocities for <111> crystal orientation.
; 
;   (from INITIAL_POS procedure)
;    EXPDIST=NUM - Set to choose NUM random z0 positions from an exponential
;                  distribution. Default is Cathode (AC) irradiation of the
;                  detector which places the majority of events near z0=0.0 cm. 
;
;                  NOTE:
;                  -EXPDIST overwrites values (if present) in z0
;                  -Set EXPDIST=1 (same as /EXPDIST) if NUM is already 
;                   defined by e0, x0, or y0. (NUM set by EXPDIST will 
;                   dominate if [SA1]not equal to 1)
;                  -Attenuation of exponential is taken from a single e0
;                   value of 60, 122, or 662 keV (hard-coded) OR keyword
;                   EVENT_ATTEN (see below).
;                  -**if no attenuation value is available from above 
;                   sources, then NUM random z0 positions are chosen from a FLAT
;                   distribution (Default setting of e0 does NOT set the
;                   attenuation!! - to use 122 keV atten, set e0=[122] ) 
;
;    /RANDOM  -  Set to choose NUM random values for BOTH x0 and y0 from 
;                flat distributions spanning the mid-gap to mid-gap region
;                around the central strip (0.3 - 0.5 cm).
;
;                NOTE: 
;                -Random values overwrite elements (if present) in x0 and y0.
;                -The number of elements initially in x0 & y0 is nulled when
;                 /RANDOM is set. Thus, the input number of events (NUM) will
;                 be taken from EXPDIST (if set), z0 (if defined), e0 (if
;                 defined), or default (2000) in that order. 
;                -To randomly set only ONE of x0 or y0, pass in an undefined
;                 array for that parameter rather than using /RANDOM manually.
;
;    EVENT_ATTEN = ATTEN - Set to ATTEN in cm^-1 to choose Z-positions from
;                   an exponential with this attenuation. Applies ONLY if
;                   /EXPDIST is set. 
;
;                   NOTE:
;                   -Overrides value determined by e0 (if present). 
;
; OUTPUT: q - a structure with TAGS described below. 
;      Let: 
;      n_evt = total number of input events
;      n_time_steps = total number of time steps before the last charge 
;                     hit a detector edge
;
;      The tags of q are:     
;      {
;       ANO:  Array[n_evt, n_time_steps] containing the time resolved signal
;             induced on the anode in units of keV. 
;       CAT:   Array[n_evt, n_time_steps] containing the time resolved signal
;              induced on the cathode in units of keV.
;       HTIME_ANO:  Array[n_evt] containing the half-rise time of the anode
;                   signal in ns.
;       TDIF:   Array[n_evt] containing the collection time difference between
;               the anode half-rise time and cathode half-rise time (Anode -
;               Cathode) in ns. 
;       X0: Array[n_evt] containing the initial X position of the event (cm).
;       Y0: Array[n_evt] containing the in:qitial Y position of the event (cm).
;       Z0: Array[n_evt] containing the initial Z position of the event (cm).
;       E0: Array[n_evt] containing the initial energy of the event in keV.
;       }
; 
;    - q can be output as a named variable 'q' in the QTRANS call, 
;      or saved in a file. Use OUTFILE = 'outfilename.sav' to save as an 
;      IDL save file or OUTFILE = 'outfilename.fits' to save as a FITS 
;      file (recommended). Files are saved to the hard-coded directory
;      OUTFILE_DIRECTORY = '~samrose/qtrans/outfiles/'.
;
; HISTORY: Based on QTRANS written by S. Boggs (see QTRANS.PRO for more)
;          Mod by S. Amrose UCB/SSL   2002
;***************************************************************************** 
PRO qtrans_new,z0,x0,y0,e0,q,num,outfile=outfile,$
               notrap=notrap,$
               anode_irradiate=anode_irradiate,$
               egrid=egrid,wgrid=wgrid,_ref_extra=ex,fac_e1=fac_e1,fac_h1=fac_h1
               ;through DRIFT_DIR: drift111=drift111
               ;through INITIAL_POS: EXPDIST=NUM OR 1, /RANDOM,  
               ;EVENT_ATTEN= # in cm^-1 

IF n_params() EQ 0 THEN BEGIN
  print,'syntax- qtrans,z0,x0,y0,e0 (keV),q,outfile=outfile,notrap=notrap,z0_expdist=z0_expdist,event_energy=event_energy,anode_irradiate=anode_irradiate,egrid=egrid,wgrid=wgrid,_ref_extra=ex'
  return
ENDIF
;+==
COMMON FACTOR_CAL
if not keyword_set(fac_e1) then fac_e=1.0 else fac_e=fac_e1
if not keyword_set(fac_h1) then fac_h=1.0 else fac_h=fac_h1
print,'fac_e/fac_h',fac_e,fac_h
;+==
!except=1
start_time = systime(1)
printf,-2,'Start time: ',systime(0)

;------------------------
;Define the COMMON BLOCKS
;------------------------
;resolve_routine,'qmove',/no_recompile,/is_function
define_common_blks,num,/cgs	; This defines all numeric parameters.
;*Note that you must run define_common_blks before successfully compiling qtrans!
COMMON FLDGRID
COMMON PIXEL
COMMON DETECTOR
COMMON DIR ;,  outfile_dir
field_dir='idlsave/'
wfieldfile='NCT14_D'+num+'_wfield.idlsave' ;Use wfiled to calculate dQ
efieldfile='NCT14_D'+num+'_v_efield.idlsave'

;outfile_dir='~/nct/sigsims/'
;outfile_dir='~/Desktop/Am_dcal/idlsave/'
outfile = 'outfilename.sav'

;VOLT= keyword_set(voltage)? voltage : -1500.d
VOLT=keyword_set(voltage)?voltage : 1200.d

;no trapping yet
notrap=keyword_set(notrap)? 1:0
notrap=1 & print,'No Trapping'  ;hard coded disable of trapping

QE=1.60217733d-19 ; in double precision, as opposed to QELEC in COMMON DETECTOR
;----------------------
;RESTORE VECTOR FIELDS
;----------------------
printf,-2,format='("Restoring vector fields...",$)'
if not keyword_set(egrid) then begin
  restore,field_dir+efieldfile
  print,'restored ',field_dir+efieldfile
endif else egrid=egrid
if not keyword_set(wgrid) then begin
  restore,field_dir+wfieldfile
  print,'                          restored ',field_dir+wfieldfile
  wgrid=eugrid[*,*,5:107-5-1]
  eugrid=0
endif else wgrid=wgrid
print,''

;----------------------------
;SORT OUT INITIAL POSITIONS
;----------------------------

initial_pos,z0,x0,y0,e0,q0,num,atten,_extra=ex ;KEYWORDS: EXPDIST=NUM OR 1, /RANDOM,  
                                               ;EVENT_ATTEN= # in cm^-1 
print,'Total Objects: ',ntostr(num)

;-----------------------------
;ANODE or CATHODE IRRADIATION?
;-----------------------------
if keyword_set(anode_irradiate) then BEGIN
  z0=1.5-z0                     ;zmax =1.5(cm)
  print,'ANODE (DC) IRR'
ENDIF ELSE IF keyword_set(expdist) THEN print,'CATHODE (AC) IRR'

printf,-2,format=$
 '("Transporting charges/tracing field lines (will take forever)...",$)'

;save_path=keyword_set(save_path)? save_path: 0
;multich=keyword_set(multich)? multich:0
;rstep=keyword_set(rstep)? rstep:0
;numsites=keyword_set(numsites)? numsites:0
;zipfilename=zipfilename

;------------------
;TRANSPORT CHARGES
;------------------
q=qmove(x0,y0,z0,q0,egrid,wgrid,notrap=notrap,multich=1,_extra=ex)
			; This is the main procedure for charge transport.
if not keyword_set(egrid) then egrid=0b
if not keyword_set(wgrid) then  wgrid=0b
q.e0=e0

printf,-2,'DONE!!'
;------------------------------------------
;RECORD OUTPUT INTO A FITSFILE OR .SAV FILE
;------------------------------------------
if keyword_set(outfile) then BEGIN
  IF (strpos(outfile,'.fits'))[0] NE -1 THEN BEGIN
    mwrfits,q,outfile_dir+outfile,/create 
    outstyle='fits format'
  ENDIF ELSE BEGIN
    save,q,filename=outfile
    outstyle='IDL Save format'
  ENDELSE

  print,'Saved outfile in '+outstyle+' to: ',outfile_dir+outfile
endif

end_time = systime(1)
printf,-2,'Time taken=',(end_time-start_time)/60.,' minutes.'

;RETURN
END ; qtrans

