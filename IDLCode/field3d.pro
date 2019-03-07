;*****************************************************************************
; FILE NAME: FIELD3D
;
; PURPOSE: Calculate a scalar electric or weighting field using relaxation
; INCLUDES: Geometry Functions
;           Field3d
; FULL DOCUMENTATION: Search for ';+****' or scroll below
;***************************************************************************

;*****************************************************************************
;simple_test_geom defines the contact/detector geometry. It sets up a 
;simple system with one air/Ge bndry in the middle:
;XXXXXXXXXXX GND
;*********** Air
;........... Ge
;@@@@@@@@@@@ HV 
;This one designed for ELECTRIC FIELD         
;***************************************************************************** 
FUNCTION simple_test_geom,IMAX,JMAX,KMAX,KJUMP              
COMMON TYPES
pnt=intarr(imax,jmax,kmax,/nozero)
pnt(*,*,*)=semiconductor  ;define all as Semiconductor 

pnt(*,*,KJUMP:KMAX-2)=air ;redefine upper section as air        

pnt(*,*,KMAX-1)=gnd  ;define top plane as ground
pnt(*,*,0)=hv ;define bottom plane as HV       
      
RETURN,pnt    
END   

;***************************************************************************
;efgeom_plus  - just like efgeom, but with two 1cm air buffers at either z end
; See efgeom in FIELD3D_GEOMFUNCT.PRO
;*****************************************************************************
FUNCTION efgeom_plus,IMAX,JMAX,KMAX,KAIR          
COMMON PIXEL
COMMON TYPES

pnt=intarr(IMAX,JMAX,KMAX,/nozero)      
pnt(*,*,*)=semiconductor  ;define all as Semiconductor 
  
;Define a ground plane on each side.
pnt(*,*,0)=gnd
pnt(*,*,kmax-1)=gnd

if kair ne 0 then begin
  pnt(*,*,1:kair-1)=air         ;Air between ground plane and HV strips
  pnt(*,*,kmax-kair:kmax-2)=air ;air between ground plane and ground strips
endif

for i=0,IMAX-1 do begin
 ip=(i+ICNTCT/2) MOD ISTRIP
 if (ip lt ICNTCT) then begin
    pnt(i,0:JMAX-1,kmax-kair-1)=gnd   ;redefine X-strips as ground       
 endif
endfor       

for j=0,JMAX-1 do begin
 jp=(j+ICNTCT/2) MOD ISTRIP
 if (jp lt ICNTCT) then begin
    pnt(0:IMAX-1,j,kair)=hv   ;redefine Y-strips as HV       
 endif
endfor       
RETURN, pnt
END

;*****************************************************************
;WFGEOM defines the contact/detector geometry
;This one designed for WEIGHTING FIELD
;*****************************************************************
FUNCTION wfgeom,IMAX,JMAX,KMAX,KAIR
COMMON TYPES
COMMON PIXEL
pnt=intarr(IMAX,JMAX,KMAX,/nozero)      
pnt(*,*,*)=semiconductor  ;define all as Semiconductor 
  
;Define a ground plane on each side.
pnt(*,*,0)=gnd
pnt(*,*,kmax-1)=gnd

pnt(*,*,1:kair-1)=air ;Air between ground plane and HV strips
pnt(*,*,kmax-kair:kmax-2)=air ;air between ground plane and ground strips

for i=0,IMAX-1 do begin
 ip=(i+ICNTCT/2) MOD ISTRIP
 if (ip lt ICNTCT) then begin
    pnt(i,0:JMAX-1,kmax-kair-1)=gnd   ;redefine X-strips as ground       
 endif
endfor       

;set central strip in upper plane to HV =1.0 for wfield
pnt((imax/2)-(ICNTCT/2):(imax/2)+(ICNTCT/2),*,kmax-kair-1)=HV

for j=0,JMAX-1 do begin
 jp=(j+ICNTCT/2) MOD ISTRIP
 if (jp lt ICNTCT) then begin
    pnt(0:IMAX-1,j,kair)=gnd   ;redefine Y-strips as ground       
 endif
endfor       

return,pnt
END ; wfgeom



;*****************************************************************
; PROCEDURE NAME: Planar_boundary
; PURPOSE: Apply planar boundary conditions
;*****************************************************************
pro planar_bndry,a,b,c,d,g,dz1
COMMON FLDGRID

;Right edge of detector
;blas_axpy,g,1.,replicate(1.0/DX^2,JMAX),2,[IMAX-1,0,0],3,KIND
R = replicate(1.0/DX^2, JMAX)
FOR K = 0, KMAX-1 DO G[IMAX-1, *, K] = G[IMAX-1, *, K] + R

;Left edge of detector
;blas_axpy,g,1.,replicate(1.0/DX^2,JMAX),2,[0,0,0],3,KIND
R = replicate(1.0/DX^2, JMAX)
FOR K = 0, KMAX-1 DO G[0, *, K] = G[0, *, K] + R

;Back edge of detector
;blas_axpy,g,1.,replicate(1.0/DY^2,IMAX),1,[0,JMAX-1,0],3,KIND
R = replicate(1.0/DY^2, IMAX)
FOR K = 0, KMAX-1 DO G[*, JMAX-1, K] = G[*, JMAX-1, K] + R

;Front edge of detector
;blas_axpy,g,1.,replicate(1.0/DY^2,IMAX),1,[0,0,0],3,KIND
R = replicate(1.0/DY^2, IMAX)
FOR K = 0, KMAX-1 DO G[*, 0, K] = G[*, 0, K] + R

;Top edge of detector
;blas_axpy,g,1.,replicate(1.0/(dz1[0,0,KMAX-1])^2,IMAX),1,[0,0,KMAX-1],2,JIND	;Only true when DZ1=DZ2
R = replicate(1.0/(dz1[0,0,KMAX-1])^2,IMAX)
FOR J = 0, JMAX-1 DO G[*, J, KMAX-1] = G[*, J, KMAX-1] + R

;Bottom edge of detector
;blas_axpy,g,1.,replicate(1.0/(dz1[0,0,0])^2,IMAX),1,[0,0,0],2,JIND	;Only true when DZ1=DZ2
R = replicate(1.0/(dz1[0,0,0])^2,IMAX)
FOR J = 0, JMAX-1 DO G[*, J, 0] = G[*, J, 0] + R

return
end ;planar_bndry



;*****************************************************************
; PROCEDURE NAME: Symmetric_boundary
; PURPOSE: Apply symmetric boundary conditions
;*****************************************************************
pro symmetric_bndry,a,b,c,d,g,dz1
COMMON FLDGRID

;Boundaries values


;Right edge of detector
; KIND1=indgen(KMAX)
; replicate_inplace,b,2.0/dx^2,2,[IMAX-1,0,0],3,KIND1
b[IMAX-1, *, *] = 2.0/dx^2

;note that it is not necessary to set the i+1 value (a) to zero, this is done by resizing the array below. Same goes for the remaining arrays. 
;Left edge of detector
;KIND2=indgen(KMAX)
;replicate_inplace,a,2.0/dx^2,2,[0,0,0],3,KIND2
a[0, *, *] = 2.0/dx^2

;Back edge of detector
;KIND3=indgen(KMAX)
;replicate_inplace,d,2.0/dy^2,1,[0,JMAX-1,0],3,KIND3
d[*, JMAX-1, *] = 2.0/dy^2

;Front edge of detector
;KIND4=indgen(KMAX)
;replicate_inplace,c,2.0/dy^2,1,[0,0,0],3,KIND4
c[*, 0, *] = 2.0/dy^2

;Top edge of detector
; Orig: blas_axpy,g,1.,replicate(1.0/(dz1[0,0,KMAX-1])^2,IMAX),1,[0,0,KMAX-1],2,JIND	;Only true when DZ1=DZ2
; AZ:
H = replicate(1.0/(dz1[0,0,KMAX-1])^2,IMAX)
for J = 0, JMAX-1 DO G[*, J, KMAX-1] = G[*, J, KMAX-1] + H
 
;Bottom edge of detector
; Orig: blas_axpy,g,1.,replicate(1.0/(dz1[0,0,0])^2,IMAX),1,[0,0,0],2,JIND	;Only true when DZ1=DZ2
; AZ:
H = replicate(1.0/(dz1[0,0,0])^2,IMAX)
for J = 0, JMAX-1 DO G[*, J, 0] = G[*, J, 0] + H

     
return
end ;symmetric_bndry



;+****************************************************************
; NAME: FIELD3d
;	
; PURPOSE: Calculate a scalar electric or weighting field using 
;          relaxation techniques 
;	
; CALLING SEQUENCE: field3d,outfile,/wfield,/conduct,
;                   maxits=maxits,tol=tol,/condsurf100
;	
; OPTIONAL INPUTS:
;       /wfield - make a weighting field instead of electric 
;                 field (default)
;       /conduct - Use surface conductivity boundary conditions
;       /condsurf100 - Multiply surface conductivity by 100
;      * always set ONE of the next two keywords *
;       maxits = Max number of iterations in relaxation to use
;       tol = Stop relaxation when the absolute difference 
;             between iterations is less than TOL  
;
; OUTPUTS:  outfile - Name of IDL save file to take output 
;                     field array EUGRID. 
;		
; NOTES: -This program will NOT compile unless you first call the 
;        procedure DEFINE_COMMON_BLKS.PRO, which defines the 
;        common blocks used in FIELD3D. Do this at least once
;        per IDL session: IDL> define_common_blks,/cgs. 
;        Define_common_blks is also called IN field3d, so you 
;        don't need to manually run it again even if changes are made.
;
;        -The output array is saved as EUGRID for both electric
;        and weighting fields.  
;	
; EXAMPLE:
;	- make a weighting field
;    IDL> define_common_blks,/cgs   ;if not already called 
;    IDL> field3d, 'wfield_out.sav', /conduct, /wfield, tol=1e-6
;
; PROCEDURES CALLED:
;	DEFINE_COMMON_BLKS
;
; REVISION HISTORY: Based on FIELD3D by Steve Boggs, Mark Peng 
;                   Mod by Susan Amrose    UCB/SSL   2002
;-****************************************************************
PRO field3d,num,outfile,wfield=wfield,conduct=conduct,$
             maxits=maxits,tol=tol,condsurf100=condsurf100


  start_time=systime(1)
                       
; LAPLACE'S EQUATION AND BOUNDARY CONDITIONS                  
;      1 = Ground                 
;      2 = High Voltage (at volt, below)                 
;      4 = Air                 
;      0 = Germanium          
;   The variable pnt2(i,j) is different.  If pnt2(i,j) is 0, then the pot- 
;   ential of the point is fixed.  If pnt2(i,j) is 1, then the potential can 
;   be varied during the relaxation. */                 
; GND=1 & HV=2 & Air=4 & Semiconductor=0 & Steering=5
; These type numbers are no longer used explicitly in the programme;
; they are retained here for backward compatibility.
; Fixed=0b & Variable=1b
; Note: It used to be (fixed,variable)=(1,0).  Hubert changed it on 1- 5-2000
; in order to facilitate simple matrix operation.  See the main loop below.

  ;--------------------
  ; DEFINE COMMON BLKS
  ;--------------------
  ;notice that this is different from the previous versions which use /si. I am using /cgs to comply with qtrans.pro
 ; define_common_blks,num,/cgs 
  COMMON PIXEL
  COMMON fldgrid
  COMMON TYPES
  COMMON DETECTOR


  ;------------------
  ; SET CONDUCTIVITY
  ;------------------
  ;a keyword to multiply the conductivity by 100
  if keyword_set(condsurf100) then begin 
    CNDSURF=100.*1.0E-12
    print,'Conductivity of surface multiplied by 100'
  ENDIF


  ;-------------------------
  ; CHOOSE WFIELD OR EFIELD
  ;-------------------------
  if keyword_set(wfield) then begin ;Compute a weighting field
    VOLT=1.0
    NELEC=0.0
    bndry_type='planar_bndry'   ;use Planar_bndry conditions. 
    geom_funct='wfgeom'
  endif else begin
   ; VOLT=double(-1000.)         ;electrode HV [Volts]
   ; VOLT=double(1000.)           ;electrode HV [Volts] :D1,D4,D6,D7,D11,D12
   ; VOLT=double(1200.)           ;electrode HV ;D2,D9,D10 
   ; VOLT=double(1500.)            ;electrode HV ;D3,D5,D8,
           ;     D1    D2    D3    D4    D5    D6    D7    D8    D9    D10   D11   D12
     VOLT_TEST=[1000.,1200.,1500.,1000.,1500.,1000.,1000.,1500.,1200.,1200.,1000.,1000.]
     VOLT=VOLT_TEST[num-1]
    ;NOTE that NELEC is also set in define_common_blks. I don't
    ;know why it's redefined below. It causes confusion.
 ;   NELEC=double(-4.3e9)       ; changed MEB 060817 for new run
   ; NELEC=double(-6.0e9)       ;electron charge density [cm-3];
    ;NELEC=double(-2.0e9)        ;electron charge density [cm-3] :D12
    ; NELEC=double(-5.0e9)         ;D7
    ; NELEC=double(-4.0e9)         ;D6
    ; NELEC=double(-4.8e9)         ;D1 
    ; NELEC=double(-7.0e9)         ;D2
    ; NELEC=double(-8.0e9)         ;D5
    ; NELEC=double(-3.3e9)         ;D8
    ; NELEC=double(-1.0e9)         ;D11
                ;    D1    D2     D3      D4     D5     D6     D7     D8     D9      D10    D11    D12
     NELEC_TEST=[-4.8e9,-7.0e9,-1.0e10,-5.0e9,-8.0e9,-4.0e9,-5.0e9,-3.3e9,-1.0e10,-6.0e9,-1.0e9,-2.0e9]
     NELEC=NELEC_TEST[num-1]
     print,VOLT,NELEC
;  NELEC=double(-5.0e10)  ;CHANGE THIS BACK TO THE ABOVE FOR NORMAL RUNS
;  print,'YOU ARE USING A DIFFERENT NELEC!!!!!'
    bndry_type='symmetric_bndry' ;use sym bndry condition
    geom_funct='efgeom_plus'
  endelse


 pnt=call_FUNCTION(geom_funct,IMAX,JMAX,KMAX,KAIR)

  
  ;-------------------------
  ; SET REMAINING CONSTANTS
  ;-------------------------
  RHO=QELEC*NELEC               ;internal charge density [C/cm-3]
  ;Define the maximum # of iterations 
  if not keyword_set(maxits) then MAXITS=10000 
  tol=keyword_set(tol)? tol:0 ;when the change in field values between iterations for all points is less than tol, the program is exited before reaching the max number of iterations. If it reaches MAXITS defined above (intended to be very high if it is not set by user) it will exit regardless of whether tol has been set.  


  ;------------------
  ;INITIALIZE FIELDS
  ;------------------
; Here we find the coefficients of the differential equation for the
; potential, using the geometry and electric field parameters
; read in previously.
  dz1=replicate(DZA,IMAX,JMAX,KMAX)
  ;the delta z value is different in the 'air' portion.
  dz1[*,*,kmax-kair:KMAX-1]=DZB 
  dz1[*,*,0:kair-1]=DZB
  dz2=dz1
  
  ;replicate_inplace,dz2,DZB,1,[0,0,kmax-kair-1],2,JIND
  ;replicate_inplace,dz2,DZB,1,[0,0,kair],2,JIND
  dz2[*,*,kmax-kair-1] = DZB
  dz2[*,*,kair] = DZB
  
  DZ1KJUMP=DZA
  DZ2KJUMP=DZB
  
  u= VOLT   * (pnt eq HV) 	;Start at zero potential

  ;--------------------------------
  ; DEFINE RELAXATION COEFFICIENTS
  ;--------------------------------
; Note that the following coefficients are for semiconductor and air only.
; HV and ground points are fixed, and pnt2[] takes care of this.
; See the main loop below.

; Coefficients for the bulk
  a=replicate(1.0/DX^2,IMAX,JMAX,KMAX)
  b=a
  c=replicate(1.0/DY^2,IMAX,JMAX,KMAX)
  d=c
  e=2.0/(dz2*(dz1+dz2))
  f=2.0/(dz1*(dz1+dz2))
  g=-2.0/DX^2-2.0/DY^2-2.0/(dz1*dz2) 
  
  h=-RHO/EPS * (pnt eq Semiconductor) ;or GND or HV, but these won't matter.
  h(*,*,kair)=0.0 
  h(*,*,kmax-kair-1)=0.0   ;set h = 0 at air/Ge bndry. This is very important! 
  pnt2=	(pnt eq Semiconductor) or $
    (pnt eq Air)                ;But allow h to change


  ;-------------------------
  ; CHECK FOR INITIAL GUESS
  ;-------------------------
  if keyword_set(u0)$
    then blas_axpy,u,1.0,u0*pnt2 ;Use initial guess of u[] if provided.
  

  ;----------------------------------------------------
  ; SET EDGE BNDRY VALUES, USE SYM OR PLANAR CONDITIONS
  ;----------------------------------------------------
  print,bndry_type
  call_procedure,bndry_type,a,b,c,d,g,dz1 



  ;-----------------------------------
  ; SET CONDUCTIVITY BNDRY CONDITIONS
  ;-----------------------------------
;if 0 then begin ;removing for very_simple_geo test
  If keyword_set(conduct) then BEGIN
    print,"Creating conductivity boundary conditions"
;   Boundary between air and detector

;   lower bndry, air to Ge
    ;replicate_inplace,f,0.0,1,[0,0,kair],2,JIND
    ;replicate_inplace,g,-2./dx^2-2./dy^2-(CNDBULK/CNDSURF)*(1./DZ1KJUMP),1,[0,0,kair],2,JIND
    ;replicate_inplace,e,(CNDBULK/CNDSURF)*(1./DZ1KJUMP),1,[0,0,kair],2,JIND
    
    f[*, *, KAIR] = 0.0
    g[*, *, KAIR] = -2./dx^2-2./dy^2-(CNDBULK/CNDSURF)*(1./DZ1KJUMP)
    e[*, *, KAIR] = (CNDBULK/CNDSURF)*(1./DZ1KJUMP)
    
;   upper bndry, Ge to air 
    ;replicate_inplace,e,0.0,1,[0,0,kmax-kair-1],2,JIND
    ;replicate_inplace,g,-2./dx^2-2./dy^2-(CNDBULK/CNDSURF)*(1./DZ1KJUMP),1,[0,0,kmax-kair-1],2,JIND
    ;replicate_inplace,f,(CNDBULK/CNDSURF)*(1./DZ1KJUMP),1,[0,0,kmax-kair-1],2,JIND
    
    e[*, *, KMAX-KAIR-1] = 0.0
    g[*, *, KMAX-KAIR-1] = -2./dx^2-2./dy^2-(CNDBULK/CNDSURF)*(1./DZ1KJUMP)
    f[*, *, KMAX-KAIR-1] = (CNDBULK/CNDSURF)*(1./DZ1KJUMP)
    
  endif else begin
  
    print,"Creating no-conductivity boundary conditions"


  ;----------------------------------------------
  ; OR ELSE SET NO-CONDUCTIVITY BNDRY CONDITIONS
  ;---------------------------------------------- 
;   Boundary between air and detector
;     lower bndry, air to Ge
    ;replicate_inplace,a,0.0,1,[0,0,kair],2,JIND 
    ;replicate_inplace,b,0.0,1,[0,0,kair],2,JIND
    ;replicate_inplace,c,0.0,1,[0,0,kair],2,JIND
    ;replicate_inplace,d,0.0,1,[0,0,kair],2,JIND
    ;replicate_inplace,f,(EPS0/EPS)*1./DZ2KJUMP,1,[0,0,kair],2,JIND 
    ;replicate_inplace,g,(-1.0)*((EPS0/EPS/DZ2KJUMP)+1.0/DZ1KJUMP),1,[0,0,kair],2,JIND
    ;replicate_inplace,e,1.0/DZ1KJUMP,1,[0,0,kair],2,JIND
    
    a[*, *, kair] = 0.0
    b[*, *, kair] = 0.0
    c[*, *, kair] = 0.0
    d[*, *, kair] = 0.0
    e[*, *, kair] = 1.0/DZ1KJUMP
    f[*, *, kair] = (EPS0/EPS)*1./DZ2KJUMP
    g[*, *, kair] = (-1.0)*((EPS0/EPS/DZ2KJUMP)+1.0/DZ1KJUMP)
    
    
;     upper bndry, Ge to air
;     notice that e,f switch when going air/Ge or Ge/air.
    ;replicate_inplace,a,0.0,1,[0,0,kmax-kair-1],2,JIND 
    ;replicate_inplace,b,0.0,1,[0,0,kmax-kair-1],2,JIND
    ;replicate_inplace,c,0.0,1,[0,0,kmax-kair-1],2,JIND
    ;replicate_inplace,d,0.0,1,[0,0,kmax-kair-1],2,JIND
    ;replicate_inplace,e,(EPS0/EPS)*1./DZ2KJUMP,1,[0,0,kmax-kair-1],2,JIND 
    ;replicate_inplace,g,(-1.0)*((EPS0/EPS/DZ2KJUMP)+1.0/DZ1KJUMP),1,[0,0,kmax-kair-1],2,JIND
    ;replicate_inplace,f,1.0/DZ1KJUMP,1,[0,0,kmax-kair-1],2,JIND
    
    a[*, *, kmax-kair-1] = 0.0
    b[*, *, kmax-kair-1] = 0.0
    c[*, *, kmax-kair-1] = 0.0
    d[*, *, kmax-kair-1] = 0.0
    e[*, *, kmax-kair-1] = (EPS0/EPS)*1./DZ2KJUMP
    f[*, *, kmax-kair-1] = 1.0/DZ1KJUMP
    g[*, *, kmax-kair-1] = (-1.0)*((EPS0/EPS/DZ2KJUMP)+1.0/DZ1KJUMP)
    
    print, 'AZ: Are we sure that there is no copy-and-paste error here (switch index e & f)' 
    
  endelse
;endif


  ;------------------
  ; BEGIN RELAXATION
  ;------------------
;  Now we are going to solve the 3-D boundary value problem for the  
;  potential everywhere in the detector.  This is accomplished via 
;  a relaxation method from Numerical Recipies (p. 659).                 
            
  print,"Running relaxation..." 
; Free up unused array memory for better performance during main iteration.
  pnt=0b
  dz1=0b
  dz2=0b
; Resize coefficient arrays for efficiency in the main loop.
; I know this is rather stupid, but it is a hysteresis,
; and it may be easier to declare the coefficients above this way, too.
  a = a[0:IMAX-2,*,*]
  b = b[1:IMAX-1,*,*]
  c = c[*,0:JMAX-2,*]
  d = d[*,1:JMAX-1,*]
  e = e[*,*,0:KMAX-2]
  f = f[*,*,1:KMAX-1]
; save,a,b,c,d,e,f,g,h,filename='field3dabcdefgh_010629.sav'  

  u_odd = u                     ; for odd-even ordering
  rjac = (  cos(2*!pi/float(IMAX))/float(DX)^2 $
            + cos(2*!pi/float(JMAX))/float(DY)^2 $
            + cos(!pi/float(KMAX))/float(DZA)^2 ) / (  1/float(DX)^2 $
                                                       + 1/float(DY)^2 $
                                                       + 1/float(DZA)^2 ) ;(17.5.24)
  
; This formula complies with (17.5.24) in Numerical Recipes,
; but in fact, Boggs' is good enough for MAXITS=1000.
; Actually, the optimal r_jacobi is dependent on the z-coordinate, as dz 
; varies. So, it may be better to set rjac and omega as arrays.



  ;------------------
  ; FIRST INTERATION
  ;------------------
; We first do half an iteration `by hand', for efficiency, as the 
; value of omega in the first iteration follows a different rule.


  omega = 1.0                   ;(17.5.30, n=0) Chebyshev acceleration
  
; Update u[].
  resid = g*u-h                 ;(17.5.28)
  resid[0:IMAX-2,*,*] = resid[0:IMAX-2,*,*] + a * u_odd[1:IMAX-1,*,*]
  resid[1:IMAX-1,*,*] = resid[1:IMAX-1,*,*] + b * u_odd[0:IMAX-2,*,*]
  resid[*,0:JMAX-2,*] = resid[*,0:JMAX-2,*] + c * u_odd[*,1:JMAX-1,*]
  resid[*,1:JMAX-1,*] = resid[*,1:JMAX-1,*] + d * u_odd[*,0:JMAX-2,*]
  resid[*,*,0:KMAX-2] = resid[*,*,0:KMAX-2] + e * u_odd[*,*,1:KMAX-1]
  resid[*,*,1:KMAX-1] = resid[*,*,1:KMAX-1] + f * u_odd[*,*,0:KMAX-2]
  ; blas_axpy,u,(-1.)*omega,resid/g*pnt2 ;(17.5.29)
  u = u + (-1.)*omega * resid/g*pnt2
  ; Note that pnt2[] keeps the fixed points unchanged.

  omega = 1.0/(1.0-0.5*rjac^2)  ;(17.5.30, n=1/2) Chebyshev acceleration
;tol=fltarr(maxits+1)
  
;num0pnt2=n_elements(where(pnt2 eq 0)) ;the number of 0's which are averaged in to the mean each time. 
;npnt2=n_elements(pnt2)
  


  
  ;----------------------------------
  ; MAIN LOOP - REMAINING ITERATIONS
  ;----------------------------------
; Iterate until MAXITS is reached OR the TOL condition is satisfied 
  n=0L
  while n lt long(MAXITS) DO BEGIN
    print,"Iteration ",n,"of ",MAXITS
; Update u_odd[].
    resid = g*u_odd-h           ;(17.5.28)
    resid[0:IMAX-2,*,*] = resid[0:IMAX-2,*,*] + a * u[1:IMAX-1,*,*]
    resid[1:IMAX-1,*,*] = resid[1:IMAX-1,*,*] + b * u[0:IMAX-2,*,*]
    resid[*,0:JMAX-2,*] = resid[*,0:JMAX-2,*] + c * u[*,1:JMAX-1,*]
    resid[*,1:JMAX-1,*] = resid[*,1:JMAX-1,*] + d * u[*,0:JMAX-2,*]
    resid[*,*,0:KMAX-2] = resid[*,*,0:KMAX-2] + e * u[*,*,1:KMAX-1]
    resid[*,*,1:KMAX-1] = resid[*,*,1:KMAX-1] + f * u[*,*,0:KMAX-2]
    ; blas_axpy,u_odd,(-1.)*omega,resid/g*pnt2 ;(17.5.29)
    u_odd = u_odd + (-1.)*omega * resid/g*pnt2
    
    omega = 1.0/(1.0-0.25*omega*rjac^2) ;(17.5.30, n>1/2) Chebyshev acceleration
    
; Update u[].
    resid = g*u-h               ;(17.5.28)
    resid[0:IMAX-2,*,*] = resid[0:IMAX-2,*,*] + a * u_odd[1:IMAX-1,*,*]
    resid[1:IMAX-1,*,*] = resid[1:IMAX-1,*,*] + b * u_odd[0:IMAX-2,*,*]
    resid[*,0:JMAX-2,*] = resid[*,0:JMAX-2,*] + c * u_odd[*,1:JMAX-1,*]
    resid[*,1:JMAX-1,*] = resid[*,1:JMAX-1,*] + d * u_odd[*,0:JMAX-2,*]
    resid[*,*,0:KMAX-2] = resid[*,*,0:KMAX-2] + e * u_odd[*,*,1:KMAX-1]
    resid[*,*,1:KMAX-1] = resid[*,*,1:KMAX-1] + f * u_odd[*,*,0:KMAX-2]
    ; blas_axpy,u,(-1.)*omega,resid/g*pnt2 ;(17.5.29)
    u = u + (-1.)*omega * resid/g*pnt2
    
;  tol(n)=mean(resid*pnt2/g)*npnt2/float(npnt2-num0pnt2)
;print,max(abs(u)),' max abs(u)',min(abs(u(where(u*pnt2 ne 0)))),' min ne 0'

    n=n+1+(abs(mean(resid*pnt2/g)) lt tol)*maxits
    omega = 1.0/(1.0-0.25*omega*rjac^2) ;(17.5.30, n>1/2) Chebyshev acceleration
  ENDWHILE



  ;---------------------------
  ; DONE! SAVE RESULT TO FILE
  ;--------------------------- 
  eugrid=FLOAT(temporary(u)) ; FLOAT required for GDL compatibility (IDL automaticaaly creates floats)
  SAVE,eugrid,tol,FILENAME=outfile ;Save potential grid to file 

  
  ;---------------------------
  ; REPORT STATS FOR THIS RUN
  ;---------------------------
  end_time = systime(1)
  printf,-2,'Time taken=',(end_time-start_time)/60.,' minutes.'
  printf,-2,'Average last residue=',mean(resid*pnt2/g),' Volts.'
  printf,-2,'Maximum last residue=',max(resid/g*pnt2),' Volts.'
  
RETURN                 
END 
