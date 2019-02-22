; Current revision: 1.5?
; Not sure if we want to make this a permanent one though.
; $Id: define_common_blks.pro,v 1.4 2001/03/29 21:19:02 hubert Exp hubert $
;
; $Log: define_common_blks.pro,v $
; Revision 1.4  2001/03/29 21:19:02  hubert
; Initial revision; goes with r1.4 of other programmes.
;

pro define_common_blks,num,SI=SI,cgs=cgs

; In qtrans.pro, lengths are defined in [cm], as cgs formulae are used.
; In field3d.pro, SI formulae are used, and lengths are defined in [m].
; This has been done so that all formulae are directly out of textbooks
; (thus, error in converting between SI and cgs is eliminated).
; Someday, we may want to do something to unify this situation.

; If neither /SI or /cgs is set, we return an error
;  instead of assuming either, lest we assume it wrong.
  if not (keyword_set(SI) xor keyword_set(cgs))$
    then message,'Either /SI or /cgs, but not both, must be set!'

; Guidelines to defining the grid dimensions
; -To take advantage of the reflective symmetry of the pixel plane, one should
;  always model only one quadrant of it.
;
; -As it is implemented currently, the true boundaries of the grid are i=-0.5,
;  i=IMAX-0.5, j=-0.5 and j=JMAX-0.5.  The modelled geometry is assumed to
;  be reflectively symmetric about these lines.  This implies that features
;  positioned at grid boundaries must each be an odd number of grid cells
;  wide, with the central grid cell right outside the outermost grid point;
;  eg, a 5x5 pixel contact placed at the origin will occupy the cells within
;  i < 2 and j < 2.
;
; -Metal-semiconductor boundaries must be placed along grid points, and their
;  potentials should be set fixed to that of the metal;
;  consequently, metallic features must be at least 1 cell wide.
;
; -Due to the way field gradients are calculated in procedure vectgrid in
;  qtrans.pro, semiconductor features must be bound by at least 5 rows/columns
;  of grid points (ie, at least 4 cells wide if it is not at a grid boundary,
;  and 9 cells wide if it is, with the 5th cell at the boundary).
;  These numbers can be reduced to 4, 3 and 7 respectively, but the
;  computation will be quite a bit more complicated.

;(This is the xy-plane at k=KJUMP)
;
;   |<--      IPIX       -->|		"x" grid points
; -x-x-x-x-x-x-x-x-x-x-x-x-x-x-		Pixel boundaries should be at centres
;  |_|   |_____________|   |_|		 of grid cells, while other boundaries
;    |<->|<--ICNTCT -->| ->| |<-	 should be at grid points.
;    IGAP                 ISTEER

common pixel,ICNTCT,ISTRIP
  ; used in all user-level procedures
  ;  (field3d, capcalc, qtrans, qshare, plotfield, and their subroutines)
 
ICNTCT=36;76;13 ;number of grid points in a contact strip (including both endpoints)
ISTRIP=40;80;16 ;number of grid points in a strip/gap pair

common fldgrid,IMAX,IIND,JMAX,JIND,KMAX,KIND,KAIR,DX,DY,DZA,DZB,HEIGHT,NSTRIP
  ; used in field3d, capcalc, qtrans and plotfield (and their subroutines)

;Define the size of the GRID                 

;imax=10
;jmax=10
;kmax=21
;kair=0
IMAX=161;321;161;65                         ;(4*ISTRIP+1)# of x gridpoints        
JMAX=161;321;161;65                         ;# of y gridpoints          
KMAX=107;171;107                         ;# of z gridpoints          
;kmax=30
;KMAX=77
kair=5;10;5                          ;size of air buffer top and bottom

IIND=indgen(IMAX)
JIND=indgen(JMAX)
KIND=indgen(KMAX)

NSTRIP=fix((IMAX-1)/ISTRIP)        ;istrip defined above         
DX=double(8.e-3)/float(imax-1)
DY=double(8.e-3)/float(jmax-1)
;DX=double(2.E-3)/float(ISTRIP)  ;grid spacing in [meters]          
;DY=double(2.E-3)/float(ISTRIP) ;grid spacing in [meters]          
DZA=double(1.5E-2)/float(kmax-2*kair-1)  ;grid spacing in [meters]          
DZB=double(1.0E-2)/float(kair)  ;grid spacing in [meters]   
HEIGHT=double(1.5e-2)   ;depth of Ge in [meters]
  ; First define all lengths in metres, and if keyword_set(cgs), scale them
  ;  to centimetres afterwards.  This saves us from writing each length twice.
if keyword_set(cgs) then begin
  DX =DX *100.                  ;grid spacing in [cm]
  DY =DY *100.                  ;grid spacing in [cm]
  HEIGHT=HEIGHT*100.            ;thickness of detector in [cm]
  DZA=DZA*100.                  ;grid spacing in [cm]	(semiconductor)
  DZB=DZB*100.                  ;grid spacing in [cm]	(air)
endif                           ; keyword_set(cgs)

; LAPLACE'S EQUATION AND BOUNDARY CONDITIONS
;   The variable pnt(i,j,k) defines the nature of the point: either contact
;   or type of material.
;      1 = Ground
;      2 = High Voltage (at volt, above)
;      4 = Air
;      0 = Germanium, CdZnTe or other semiconductor
;      5 = Steering electrode voltage
;   The variable pnt2(i,j,k) is different.  If pnt2(i,j,k) is 0, then the pot-
;   ential of the point is fixed.  If pnt2(i,j,k) is 1, then the potential can
;   be varied during the relaxation. */

common types,GND,HV,Air,Semiconductor,Fixed,Variable
  ; used in field3d and qtrans (and their subroutines)
  GND=1 & HV=2 & Air=4 & Semiconductor=0
  ; These type numbers are no longer used explicitly in the programme;
  ; they are retained here for backward compatibility.
  Fixed=0b & Variable=1b
  ; Note: It used to be (fixed,variable)=(1,0).  Hubert changed it
  ; on 1- 5-2000 in order to facilitate simple matrix operation.
  ; See the main loop in field3d.

common detector,EPS0,EPS, NELEC,QELEC, CNDBULK,CNDSURF, trap_tau_surf
  ; All but the last three variables are used in field3d,
  ;  only EPS0 and EPS are used in capcalc, while
  ;  the last three variables are used in qmove() in qtrans.pro.
  ; This common block is not really necessary in terms of functionality,
  ;  but these variables are all properties of the detector that may change.
;trap_tau_surf= [145.d-9,145.d-9] ;[s] tau values for surface trapping [0] = electron, [1]=hole

trap_tau_surf= [1032.75d-9,1e9] 
if keyword_set(SI) then begin
  EPS0=8.854188e-12             ;permittivity of free space [F/m]
  EPS=16.0*EPS0                 ;permittivity of germanium [F/m]
;NOTE, NELEC really redefined in feild3d.pro -- need to fix!!!
; NELEC=1.0e16                  ;electron charge density [m-3] in germanium
   NELEC=-3.9e15                ;/* charge impurity density of 1e16 */
                                ;/* per cubic meter.  This is the same */
                                ;/* as 1e10 per cubic centimeter */
  QELEC=1.602177e-19            ;electron charge [Coulomb]

  CNDBULK=5.0E-10               ;Bulk conductivity [ohm-1 m-1] (!!NOT RESISTIVITY!!)
;  CNDSURF=1.0E-12               ;Surface conductivity [ohm-1]
CNDSURF=1.0E-16
endif else begin
   ; XTRAP[] is only used for germanium detectors, not CdZnTe.
    ;XTRAP=[2.5d,$	; Electron trapping length [cm]
    ;       0.1d ]	; Hole     trapping length [cm], ~0.025
    ; mu    == slope of velocity vs electric field curve;
    ; tau   == trapping time; so that
    ; mutau == trapping length per unit field strength.
 mutau=[ 1.5d-3,$              ;[cm/V] for electron
	    1.0d-5 ]	;[cm/V] for hole
 mutau0=mutau*1.0
 
; mutau0[] is used to enable a different set of (smaller) values
	; of mutau at the surfaces of the detector (the cathode and anode
	; planes).  Increased charge trapping at the surfaces is observed
	; experimentally and reported in the literature.
	; It is currently unused, but it may be needed in the future.
 EPS0=8.854188e-14              ;permittivity of free space [F/cm]
 EPS=16.0*EPS0                  ;permittivity of germanium [F/cm]
; NELEC=1.0e10                  ;electron charge density [cm-3] in germanium
 ;NELEC=-3.5e9 ;NCT'05          ;/* charge impurity density of 1e16 */
 ;NELEC=-2.0e9   ;NCT'14_D12
; NELEC=-5.0e9   ;NCT'14_D7
 ; NELEC=-4.0e9  ;NCT'14_D6
;  NELEC=-4.8e9  ;NCT'14_D1
 ; NELEC=-7.0e9  ;NCT'14_D2
      ;    D1    D2    D3      D4     D5     D6     D7     D8     D9      D10    D11    D12
   TEST=[-4.8e9,-7.0e9,-1.0e10,-5.0e9,-8.0e9,-4.0e9,-5.0e9,-3.3e9,-1.0e10,-6.0e9,-1.0e9,-2.0e9]
 ; NELEC=-8.0e9  ;NCT'14_D5
  NELEC = TEST[num-1]
  ;print,NELEC
 ; NELEC=-3.3e9   ;NCT'14_D8
 ; NELEC=-1.0e9   ;NCT'14_D11
 ;NELEC=-6.0e9   ;NCT'09_D5       ;/* per cubic meter.  This is the same */
                                ;/* as 1e10 per cubic centimeter */
 QELEC=1.602177e-19             ;electron charge [Coulomb] (no change for cgs)
 
 CNDBULK=5.0E-12                ;Bulk conductivity [ohm-1 cm-1] (!!NOT RESISTIVITY!!)
; CNDSURF=1.0E-12                ;Surface conductivity [ohm-1] (no change for cgs)
CNDSURF=1.0E-16
endelse                         ; not keyword_set(SI)

common step, DTSTEP,NUMITS, dxstep_tolerance,efield_tolerance
common FACTOR_CAL, fac_e, fac_h
  ; NUMITS and dxstep_tolerance are used in qmove(), and
  ; DXSTEP and efield_tolerance in drift_dir(); both defined in qtrans.pro.
  ; Neither is this common block necessary in terms of functionality,
  ;  but these variables are all user-defined parameters.

 DTSTEP=0.5e-9                ; [s], the constant time step for tracing charges =2 ns
; DTSTEP=0.125e-9
 ;MAXTIME=300.e-9 ;[s] 300ns
 MAXTIME=400.e-9
 NUMITS=long(MAXTIME/DTSTEP)

  ; The `tolerance' level for errors accumulated from numerical calculations;
  ; they are needed because numerical calculations seldom yield an exact zero
  ; when the true answer is zero.
 dxstep_tolerance=1.e-5	; should be less than DXSTEP, obviously.
 efield_tolerance=1.0e-20

COMMON DIR,  outfile_dir
outfile_dir='.'

end
