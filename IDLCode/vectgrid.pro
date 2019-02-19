FUNCTION vectgrid,ugrid
;PRO vectgrid
;Derive the electric and weighting fields at the grid point from the scalar
;potentials at those points.

;------------------
;GET COMMON BLOCKS
;------------------
resolve_routine,'field3d',/no_recompile,/compile_full_file ;for geometry function
COMMON FLDGRID ; IMAX,JMAX,KMAX,KAIR,DX,DY,DZA,DZB,HEIGHT,NSTRIP
COMMON PIXEL   ; ICNTCT,ISTRIP
COMMON TYPES   ; GND,HV,Air,Semiconductor,Fixed,Variable

;-----------------
;SET UP CONSTANTS
;-----------------
zmax=KMAX-2*KAIR
wucgrid=ugrid  ;copy input array
;This field (wucgrid) provides only the scalar potential field.

;---------------------------
;ESTABLISH DETECTOR GEOMETRY
;---------------------------

; Complication to this field calculation business
; -At the boundaries of the grid, we only know the potential field on one side
;  (within the grid), but not the other.  Thus, a different formula is used.
; -On the pixel plane, metal electrodes (contact and steering) create
;  additional boundaries, which also call for a different formula.  Thus,
;  the geometry of the pixel plane becomes relevant in this calculation.


;--------------------------------------
;LOCATE METAL-SEMICONDUCTOR BOUNDARIES
;--------------------------------------
; index small				large
;	-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-
;	  metal  |  semiconductor  |  metal
;	      l_edge		r_edge
;
; l_edge == left edge of a semiconductor region, r_edge == right edge

; The following grids are (IMAX-5) by JMAX and IMAX by (JMAX-5),
; excluding two and three rows/columns on each side of the pixelgeom[] grid.
;
; 0 1 2 3		IMAX-1
; x-x-x-x- . . . -x-x-x-x-x
;    |<- xl_edge[] ->|
;      |<- xr_edge[] ->|
;
; Note that i=2 and i=IMAX-3 are guaranteed not to be a right edge and a left
; edge respectively, by the definition of the geometry in field3d.pro;
; same in the y-direction.

;lower contact plane:
; First, reconstruct the pixel-plane geometry.
pixelgeom=(efgeom_plus(imax,jmax,kmax,kair))[*,*,kair]
	; efgeom() is defined in field3d.pro.
; Define all metal areas as GND here.
;  Note that there is no type `air' on the pixel plane.
pixelgeom = GND * (temporary(pixelgeom) ne Semiconductor)

yl_edge_l=    (pixelgeom[*,3:JMAX-3] eq Semiconductor)$
  and (pixelgeom[*,1:JMAX-5] eq GND)
yr_edge_l=    (pixelgeom[*,2:JMAX-4] eq Semiconductor)$
  and (pixelgeom[*,4:JMAX-2] eq GND)	
; Grids to mask out metal interiors on the pixel plane
;  (the name should actually be x_not_metal_interior)
x_nonmetal=intarr(IMAX,JMAX,2,/nozero) ;last index, 0=lower, 1=upper
x_nonmetal[1:IMAX-2,*,0]=	   (pixelgeom[0:IMAX-3,*] eq Semiconductor)$
			or (pixelgeom[2:IMAX-1,*] eq Semiconductor)
 ; The following 2 statements are simplified by the restrictions on
 ; the grid definition in field3d.pro.
x_nonmetal[0,*,0]=	   (pixelgeom[1,*]	  eq Semiconductor)
x_nonmetal[IMAX-1,*,0]=	   (pixelgeom[IMAX-2,*]	  eq Semiconductor)

y_nonmetal=intarr(IMAX,JMAX,2,/nozero)
y_nonmetal[*,1:JMAX-2,0]=    (pixelgeom[*,0:JMAX-3] eq Semiconductor)$
			or (pixelgeom[*,2:JMAX-1] eq Semiconductor)
y_nonmetal[*,0,0]=	   (pixelgeom[*,1]	  eq Semiconductor)
y_nonmetal[*,JMAX-1,0]=	   (pixelgeom[*,JMAX-2]	  eq Semiconductor)

;upper contact plane:
pixelgeom=(efgeom_plus(imax,jmax,kmax,kair))[*,*,kmax-kair-1]
	; efgeom() is defined in field3d.pro.
; Define all metal areas as GND here.
;  Note that there is no type `air' on the pixel plane.
pixelgeom = GND * (temporary(pixelgeom) ne Semiconductor)
xl_edge_u=    (pixelgeom[3:IMAX-3,*] eq Semiconductor)$
	and (pixelgeom[1:IMAX-5,*] eq GND	   )
xr_edge_u=    (pixelgeom[2:IMAX-4,*] eq Semiconductor)$
	and (pixelgeom[4:IMAX-2,*] eq GND	   )
; Grids to mask out metal interiors on the pixel plane
;  (the name should actually be x_not_metal_interior)
;these arrays are defined above
x_nonmetal[1:IMAX-2,*,1]=  (pixelgeom[0:IMAX-3,*] eq Semiconductor)$
			or (pixelgeom[2:IMAX-1,*] eq Semiconductor)
 ; The following 2 statements are simplified by the restrictions on
 ; the grid definition in field3d.pro.
x_nonmetal[0,*,1]=	   (pixelgeom[1,*]	  eq Semiconductor)
x_nonmetal[IMAX-1,*,1]=	   (pixelgeom[IMAX-2,*]	  eq Semiconductor)

y_nonmetal[*,1:JMAX-2,1]=  (pixelgeom[*,0:JMAX-3] eq Semiconductor)$
			or (pixelgeom[*,2:JMAX-1] eq Semiconductor)
y_nonmetal[*,0,1]=	   (pixelgeom[*,1]	  eq Semiconductor)
y_nonmetal[*,JMAX-1,1]=	   (pixelgeom[*,JMAX-2]	  eq Semiconductor)

pixelgeom=0b	; save memory

; Now, we are ready to calculate the electric fields.
dwucdx=fltarr(IMAX,JMAX,zmax)
dwucdy=fltarr(IMAX,JMAX,zmax)
dwucdz=fltarr(IMAX,JMAX,zmax)

;Find the x-derivatives.
  ;Simple 4th-order derivative at ordinary points
dwucdx[2:IMAX-3,*,*]$
 =(-1./DX)*( (2./ 3.)*(wucgrid[3:IMAX-2,*,*]-wucgrid[1:IMAX-4,*,*])$
	    -(1./12.)*(wucgrid[4:IMAX-1,*,*]-wucgrid[0:IMAX-5,*,*]) )

  ;3rd-order fit to field at left boundaries
dwucdx[0:1,*,*]=(-1./DX)*(  3.    *(wucgrid[1:2,*,*]-wucgrid[0:1,*,*])$
			  - 1.5   *(wucgrid[2:3,*,*]-wucgrid[0:1,*,*])$
			  +(1./3.)*(wucgrid[3:4,*,*]-wucgrid[0:1,*,*]) )

  ;3rd-order fit to field at right boundaries
dwucdx[IMAX-2:IMAX-1,*,*]$
 =(-1./DX)*(  3.    *(wucgrid[IMAX-2:IMAX-1,*,*]-wucgrid[IMAX-3:IMAX-2,*,*])$
	    - 1.5   *(wucgrid[IMAX-2:IMAX-1,*,*]-wucgrid[IMAX-4:IMAX-3,*,*])$
	    +(1./3.)*(wucgrid[IMAX-2:IMAX-1,*,*]-wucgrid[IMAX-5:IMAX-4,*,*]) )

  ; Metal-semiconductor boundaries on pixel plane
    ; upper plane Semiconductor left edges 
dwucdx[2:IMAX-4,*,ZMAX-1]$	; First, reset the fields to null.
 = (xl_edge_u eq 0) * dwucdx[2:IMAX-4,*,ZMAX-1]$
  + xl_edge_u       * (-1./DX)$	; Next, calculate the boundaries properly.
  *(  3.    *(wucgrid[3:IMAX-3,*,ZMAX-1]-$
              wucgrid[2:IMAX-4,*,ZMAX-1])$
      - 1.5   *(wucgrid[4:IMAX-2,*,ZMAX-1]-$
                wucgrid[2:IMAX-4,*,ZMAX-1])$
      +(1./3.)*(wucgrid[5:IMAX-1,*,ZMAX-1]-$
                wucgrid[2:IMAX-4,*,ZMAX-1]) )
;stop=here
;recall that wucgrid has been resized to be zmax elements BEFORE entering program
xl_edge_u=0b
    ; upper plane Semiconductor right edges 
dwucdx[3:IMAX-3,*,ZMAX-1]$	; First, reset the fields to null.
  = (xr_edge_u eq 0) * dwucdx[3:IMAX-3,*,ZMAX-1]$
  + xr_edge_u       * (-1./DX)$	; Next, calculate the boundaries properly.
*(  3.    *(wucgrid[3:IMAX-3,*,ZMAX-1]-$
              wucgrid[2:IMAX-4,*,ZMAX-1])$
      - 1.5   *(wucgrid[3:IMAX-3,*,ZMAX-1]-$
                wucgrid[1:IMAX-5,*,ZMAX-1])$
      +(1./3.)*(wucgrid[3:IMAX-3,*,ZMAX-1]-$
                wucgrid[0:IMAX-6,*,ZMAX-1]) )
;old way, this uses U(i-1)-U(i) when it should be U(i)-U(i-1)
;  *(  3.    *(wucgrid[2:IMAX-4,*,ZMAX-1]-$
;              wucgrid[3:IMAX-3,*,ZMAX-1])$
;      - 1.5   *(wucgrid[1:IMAX-5,*,ZMAX-1]-$
;                wucgrid[3:IMAX-3,*,ZMAX-1])$
;      +(1./3.)*(wucgrid[0:IMAX-6,*,ZMAX-1]-$
;                wucgrid[3:IMAX-3,*,ZMAX-1]) )
xr_edge_u=0b

    ; Finally, mask out all the metal interior on both the upper and lower plane.
dwucdx[*,*,ZMAX-1]= dwucdx[*,*,ZMAX-1] * x_nonmetal[*,*,1]
dwucdx[*,*,0]= dwucdx[*,*,0] * x_nonmetal[*,*,0] ;lower plane
x_nonmetal=0b

;save,dwucdx,filename='~/temp_wgrid_dwucdx.sav'
;print,'saved: ~/temp_wgrid_dwucdx.sav'

;Find the y-derivatives.
  ;Simple 4th-order derivative at ordinary points
dwucdy[*,2:JMAX-3,*]$
 =(-1./DY)*( (2./ 3.)*(wucgrid[*,3:JMAX-2,*]-wucgrid[*,1:JMAX-4,*])$
	    -(1./12.)*(wucgrid[*,4:JMAX-1,*]-wucgrid[*,0:JMAX-5,*]) )

  ;3rd-order fit to field at left boundaries
dwucdy[*,0:1,*]=(-1./DY)*(  3.    *(wucgrid[*,1:2,*]-wucgrid[*,0:1,*])$
			  - 1.5   *(wucgrid[*,2:3,*]-wucgrid[*,0:1,*])$
			  +(1./3.)*(wucgrid[*,3:4,*]-wucgrid[*,0:1,*]) )

  ;3rd-order fit to field at right boundaries
dwucdy[*,JMAX-2:JMAX-1,*]$
 =(-1./DY)*(  3.    *(wucgrid[*,JMAX-2:JMAX-1,*]-wucgrid[*,JMAX-3:JMAX-2,*])$
	    - 1.5   *(wucgrid[*,JMAX-2:JMAX-1,*]-wucgrid[*,JMAX-4:JMAX-3,*])$
	    +(1./3.)*(wucgrid[*,JMAX-2:JMAX-1,*]-wucgrid[*,JMAX-5:JMAX-4,*]) )

  ; Metal-semiconductor boundaries on pixel plane
    ; lower plane Semiconductor left edges 

dwucdy[*,2:JMAX-4,0]$	; First, reset the fields to null.
 = (yl_edge_l eq 0) * dwucdy[*,2:JMAX-4,0]$
  + yl_edge_l       * (-1./DY)$	; Next, calculate the boundaries properly.
  *(  3.    *(wucgrid[*,3:JMAX-3,0]-$
              wucgrid[*,2:JMAX-4,0])$
      - 1.5   *(wucgrid[*,4:JMAX-2,0]-$
                wucgrid[*,2:JMAX-4,0])$
      +(1./3.)*(wucgrid[*,5:JMAX-1,0]-$
                wucgrid[*,2:JMAX-4,0]) )
;recall that wucgrid has been resized to be zmax elements
yl_edge_l=0b
    ; lower plane Semiconductor right edges 
dwucdy[*,3:JMAX-3,0]$	; First, reset the fields to null.
  = (yr_edge_l eq 0) * dwucdy[*,3:JMAX-3,0]$
  + yr_edge_l       * (-1./DY)$	; Next, calculate the boundaries properly.
*(  3.    *(wucgrid[*,3:JMAX-3,0]-$
              wucgrid[*,2:JMAX-4,0])$
      - 1.5   *(wucgrid[*,3:JMAX-3,0]-$
                wucgrid[*,1:JMAX-5,0])$
      +(1./3.)*(wucgrid[*,3:JMAX-3,0]-$
                wucgrid[*,0:JMAX-6,0]) )
;old way, this uses U(i-1)-U(i) when it should be U(i)-U(i-1)
;  *(  3.    *(wucgrid[2:JMAX-4,*,ZMAX-1]-$
;              wucgrid[3:JMAX-3,*,ZMAX-1])$
;      - 1.5   *(wucgrid[1:JMAX-5,*,ZMAX-1]-$
;                wucgrid[3:JMAX-3,*,ZMAX-1])$
;      +(1./3.)*(wucgrid[0:JMAX-6,*,ZMAX-1]-$
;                wucgrid[3:JMAX-3,*,ZMAX-1]) )
yr_edge_l=0b

    ; Finally, mask out all the metal interior on both the upper and lower plane.
dwucdy[*,*,ZMAX-1]= dwucdy[*,*,ZMAX-1] * y_nonmetal[*,*,1]
dwucdy[*,*,0]= dwucdy[*,*,0] * y_nonmetal[*,*,0] ;lower plane
y_nonmetal=0b

;Find the z-derivatives.
  ; ordinary points
dwucdz[*,*,2:zmax-3]$
  =(-1./DZA)*( (2./ 3.)*(wucgrid[*,*,3:zmax-2]-$
                         wucgrid[*,*,1:ZMAX-4])$
               -(1./12.)*(wucgrid[*,*,4:ZMAX-1]-$
                          wucgrid[*,*,0:ZMAX-5]) )

  ; bottom boundaries
dwucdz[*,*,0:1]=(-1./DZA)*(  3.    *(wucgrid[*,*,1:2]-wucgrid[*,*,0:1])$
			   - 1.5   *(wucgrid[*,*,2:3]-wucgrid[*,*,0:1])$
			   +(1./3.)*(wucgrid[*,*,3:4]-wucgrid[*,*,0:1]) )
;check for sign consistency between three-point fit at edges and simple slope using nearest neighbor
Check=(-1./DZA)*(wucgrid[*,*,1:2]-wucgrid[*,*,0:1])
signchange= (abs(check)/check) ne (abs(dwucdz[*,*,0:1])/dwucdz[*,*,0:1])
;if there is a sign change, use the simple slope. 
dwucdz[*,*,0:1] = signchange*Check + (1b-signchange)*dwucdz[*,*,0:1]
print,total(signchange),' sign switches for bottom bndry' 
  ; top boundaries
dwucdz[*,*,zmax-2:zmax-1]$
  =(-1./DZA)*(  3.    *(wucgrid[*,*,zmax-2:zmax-1]-$
                        wucgrid[*,*,zmax-3:zmax-2])$
                - 1.5   *(wucgrid[*,*,zmax-2:zmax-1]-$
                          wucgrid[*,*,zmax-4:zmax-3])$
                +(1./3.)*(wucgrid[*,*,zmax-2:zmax-1]-$
                          wucgrid[*,*,zmax-5:zmax-4]))

;check for sign consistency between three-point fit at edges and simple slope using nearest neighbor
Check=(-1./DZA)*(wucgrid[*,*,zmax-2:zmax-1]-$
                 wucgrid[*,*,zmax-3:zmax-2])
signchange= (abs(check)/check) ne $
  (abs(dwucdz[*,*,zmax-2:zmax-1])/dwucdz[*,*,zmax-2:zmax-1])
;if there is a sign change, use the simple slope. 
dwucdz[*,*,zmax-2:zmax-1] = signchange*Check + (1b-signchange)*dwucdz[*,*,zmax-2:zmax-1]
print,total(signchange),' sign switches for top bndry' 

wucgrid=0b

vgrid=fltarr(IMAX,JMAX,zmax,3,/nozero)
vgrid[*,*,*,0] = dwucdx
vgrid[*,*,*,1] = dwucdy
vgrid[*,*,*,2] = dwucdz


RETURN,vgrid
END ; vectgrid()


