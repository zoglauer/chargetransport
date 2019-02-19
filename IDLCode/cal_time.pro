; BIPOLAR_TRANSFER_FAST
; Transfer function in time domain for bipolar shaper from NCT fast
; channel electronics.  The processing chain consists of a (CR)^2 and
; RC.  Time 't' is in nanoseconds.
;
; MEB 071211

function bipolar_transfer_fast, t, s_flag;, flag2
  ; s_flag is for the sign of the input signal , 0 is positive signal,1 is negative signal
  if s_flag eq 0 then begin
    t5    = 335800.0d
    gain0 = 7.388661d
  endif else begin
    t5    = 374000.0d
    gain0 = -6.657754d
  endelse
  ; RC time constants (ns) for the modified fast channel
  ;==for error test
  ;idx=[0,0,0]
  ;tt1=[1.0908d2,1.1127d2,1.0691d2]
  ;tt2=[7.7880d1,7.9445d1,7.6330d1]
  ;tt3=[4.870d1, 4.9679d1,4.7731d1]
  ;t1=tt1[idx[0]]
  ;t2=tt2[idx[1]]
  ;t3=tt3[idx[2]]
  ;================
  t1 = 1.0908d2 
  t2 = 7.788d1
  t3 = 4.870d1
  t4 = 100000.0d
  
  ; five poles, s1 and s2 are complex number
  a1= -1./t1 & b1=sqrt(t1/t2-1.)/t1
  
  s1=dcomplex(a1,-1.*b1)
  s2=dcomplex(a1,b1)
  s3=-1./t3
  s4=-1./t4
  s5=-1./t5
  
  ; coefficients of partial fraction expansion:
  c1=s1^4/( (s1-s2)*(s1-s3)*(s1-s4)*(s1-s5) )
  c2=s2^4/( (s2-s1)*(s2-s3)*(s2-s4)*(s2-s5) )
  c3=s3^4/( (s3-s1)*(s3-s2)*(s3-s4)*(s3-s5) )
  c4=s4^4/( (s4-s1)*(s4-s2)*(s4-s3)*(s4-s5) )
  c5=s5^4/( (s5-s1)*(s5-s2)*(s5-s3)*(s5-s4) )
  
  h_t=c1*exp(s1*t)+c2*exp(s2*t)+c3*exp(s3*t)+c4*exp(s4*t)+c5*exp(s5*t)
  ; gain
  ;gain = 5.977376d0
  gain = 9.029197d0
  return, (t ge 0.)*(gain/t3)*gain0*real_part(h_t)
end

; CAL_TIME
; return the time which crossing the 0 point
; Input : time , signal data and timing flag.
; Output: crossing time
; 
; writen by MEB 071211
; modified by ZKLiu 071217

function cal_time, t, data,flag,trigger;,idx;,flag2
  check_trig=0
  check=0
  ;flag =0 , find the "+" -> "-" zero-point
  ;flag =1 , find the "-" -> "+" zero-point
  ;flag =2 , fing the trigger time
  ;flag =3 , return bipolar response
  ;trigger is the threshold level..
  n = n_elements(t)
  ;check the signal is AC or DC
  if data[n-1] lt 0. then s_flag=1 else s_flag=0
  
  dt = t[1]-t[0]
  
  if flag ne 3 then begin
    i=0l
    if flag eq 1 then old_val=99.0 else old_val=-99.0
    this_val=0.0
    time0=0.0
    time1=0.0
    while check eq 0 do begin 
      this_val = total(bipolar_transfer_fast(t[i]-t[0:i],s_flag) * data[0:i])*dt
      ;find the trigger time
      if flag eq 2 then begin
        if check_trig eq 0 then begin
          if (this_val gt trigger) and (old_val lt trigger) then begin
	    sn=this_val
            so=old_val
            tn=t[i]
            to=t[i-1]
            time1=(tn*so-to*sn-trigger*dt)/(so-sn)
            check_trig=1
            check=1
          endif
        endif
      endif
      
      ;find '+'->'-' zero-crossing time
      if flag eq 0 then begin
        if (old_val gt 0.) and (this_val lt 0.) then begin
          sn=this_val
          so=old_val
          tn=t[i]
          to=t[i-1]
          time0=(tn*so-to*sn)/(so-sn)
          check=1
        endif
      endif else if flag eq 1 then begin
        if (old_val lt 0.) and (this_val gt 0.) then begin
          sn=this_val
  	  so=old_val
	  tn=t[i]
	  to=t[i-1]
	  time0=(tn*so-to*sn)/(so-sn)
	  check=1
        endif
      endif
      old_val=this_val
      i=i+1l
      if i gt n-1l then begin 
        ;no cross-zero point
        check=2
        time0=-999.0
      endif  
    endwhile

    if flag eq 2 then return,time1 else return,time0 
    
  endif else begin
    signal2 = dblarr(n)
    for i = 0L, n-1 do begin
      signal2[i] = total(bipolar_transfer_fast(t[i]-t[0:i],s_flag) * data[0:i])*dt
    endfor
    return, signal2
  endelse
end

function bipolar_response_fast_test,t,data
  n = n_elements(t)
  signal2 = dblarr(n)
  dt = t[1]-t[0]
  if data[n-1] lt 0 then s_flag=1 else s_flag=0

  for i = 0L, n-1 do begin
    signal2[i] = total(bipolar_transfer_fast(t[i]-t[0:i],s_flag) * data[0:i])*dt
  endfor
  return, signal2
end
