PRO IRX_beta_plot

;   .comp /home/hdominguez/idl/pro/new/IRX_beta_plot.pro

dir='/home/hdominguez/SHARDS/DR2_gamma/redanddead_gamma'

readcol, '/home/hdominguez/SHARDS/DR2_gamma/redanddead_gamma/ID_ra_dec_SFR_beta_1.0_1.5_3889_sql.lst', num, id, ra, dec, zspec, zflag, z,  LIR, SFR_IR_all, SFR_w, SFR_1600, SFR_2800, beta, av_beta,  Format='I, A,F,F,F,F,F,F,f,f,f,f,', skipline=1


;z es combinación de zphot con zspec (cuando existe)
;renombramos para no cambiar código
zphot=z


wok=where(sfr_w gt 0)
wnok=where(sfr_w lt 0)

;readcol, '/home/hdominguez/SHARDS/DR2_gamma/redanddead_gamma/out_HDFN_SHARDS_iDR2gammaP01.20150304_z_1.00_1.50_long_info',  a,a,a,a,a,a,a,a,id_lowsf,format='a,a,a,a,a,a,a,a,a',skipline=44, numline=23

readcol, '/home/hdominguez/SHARDS/DR2_gamma/redanddead_gamma/out_HDFN_SHARDS_iDR2_P01.20150630_z_1.00_1.50_plot',  a,a,a,a,a,a,a,a,id_plot,beta,irx,sfr_ir, sfr_uv,usage, format='a,a,a,a,a,a,a,a,a,f,f,f',skipline=609, numline=1203-610

low_ir=where(usage eq 17.)



readcol, '/home/hdominguez/SHARDS/DR2_gamma/redanddead_gamma/out_HDFN_SHARDS_iDR2_P01.20150630_z_1.00_1.50_plot', N ,object ,RA  ,DEC ,redshift ,q,VJ ,UV ,mass , sSFR ,SFR_UV_obs ,SFR_IR2 ,SFR_tot , Av, beta2,  mag_H, mag_j,format='i,a,f,f,f',skipline=1221



low_ir_all=fltarr(n_elements(low_ir))
for i=0, n_elements(low_ir)-1 do low_ir_all(i)=where(id_plot(low_ir(i)) eq object)


low_ir_all=fltarr(n_elements(low_ir))
for i=0, n_elements(low_ir)-1 do low_ir_all(i)=where(id_plot(low_ir(i)) eq id)

stop


;--------------File selección Pablo-------------------------------

; # Considering galaxies at   1.000 <= z <=   1.500
; # Considering galaxies with log(M/M_sun) >=  7.00
; # Considering galaxies with sSFR (M_sun/yr) <=  0.20
; # Considering galaxies detected in the MIR/FIR (in sSFR criterion)
; # Plotting UVJ diagrams
; # Number of galaxies after cut in redshift and mass:  3657
; # Number of galaxies selected from UVJ diagram:   112
; # Number of galaxies with Wuyts SFR(IR)  593 out of 3657 (016.2%)
; # Number of galaxies with TIR   SFR(IR)  593 out of 3657 (016.2%)
; # Fit to IRX (1600) vs. beta: a= 30.55 b= 14.15
; # Cross point of Meurer's IRX-beta curve and mine (1600):  -0.607
; # Fit to IRX (2800) vs. beta: a=  8.09 b=  3.02
; # Cross point of Meurer's IRX-beta curve and mine (2800):  -0.973
; # SFR^c(2800)/SFR^c(1600): av=   1.2372 st=   1.2313 me=   1.0866 mi=   1.0006 ma=  25.6503   3332
; # Lost galaxies: 20016 ; final_number= 3657
; # Number of galaxies selected from sSFR (2800) vs. mass diagram:   237
; # sSFR stats for UVJ selected galaxies: av=   0.0870 st=   0.0710 me=   0.0764 mi=   0.0098 ma=   0.3647     84
; # Number of galaxies selected from sSFR (1600) vs. mass diagram:   243
; # sSFR stats for UVJ selected galaxies: av=   0.0470 st=   0.0547 me=   0.0298 mi=   0.0030 ma=   0.3853     84


; # Total number of UVJ-based quiescent galaxies: 112 . With V-J>=1: 84
; # Total number of UVJ-based massive [log(M/M_sun) >= 10.0] quiescent galaxies: 85 

; # Total number of sSFR-based quiescent galaxies: 237
;# Total number of sSFR-based massive [log(M/M_sun) >= 10.0] quiescent galaxies: 152 
;# Total number of quiescent galaxies (any selection): 254
;# Total number of [log(M/M_sun) >= 10.0] quiescent galaxies (any selection): 157 
;# Total number of [log(M/M_sun) >= 10.0] <=12.0 quiescent galaxies (any selection): 154 


;======================SFR====================================

;-------------------calculamos Auv-----------------------------

;transformamos IRX-Beta en valores de Auv

;----------------1600------------------------
  
;----------------1600------------------------
           
fit_a_a1600=30.55
fit_b_a1600=14.15

irx_fit=fit_a_a1600+fit_b_a1600*beta
beta1=where(beta lt -0.607)
beta11=where(beta ge -0.607)

A_1600=2.5*alog10(irx_fit+1) 
a_1600(beta1)=1.99*beta(beta1)+4.43



;----------------2800------------------------

fit_a_a2800=8.09
fit_b_a2800=3.02

irx_fit_2800=fit_a_a2800+fit_b_a2800*beta
beta2=where(beta lt -0.973)

A_2800=2.5*alog10(irx_fit_2800+1)                                 ; fit Pablo, Ya es A_2800
a_1600_nuv=1.99*beta+4.43                                         ; Meurer 1999, es A_1600 


a_2800_meurer=a_1600_nuv*(1.7957/2.4665)                          ; de A_1600 a A_2800 
a_2800(beta2)=a_2800_meurer(beta2)                                ; Meurer 1999 para beta < -1.9: 


; Mail de pablo, Calzetti ext. law
;A_2800=1.7957*A_V
;A_1600=2.4665*A_V
 

A_V1=a_1600/2.4665             ; extinción en V
A_V2=a_2800/1.7957
 

;============FITS PLOTS=============================

;irx=sfr_w/sfr_1600
;irx2=sfr_w/sfr_2800

; Pablo PG
 x=findgen(200)/20-4    ; todos los valores de beta

y1=fit_a_a1600+fit_b_a1600*x
y2=fit_a_a2800+fit_b_a2800*x


;Líneas plot beta

  ;Meurer 1999
  a3=1.99*x+4.43
  ;Takeuchi 2012
  a2=3.06+1.58*x 
  ;this work
  a4=8.09+3.02*x


logirx=alog10(10^(0.4*a3)-1)+0.076                       ; Meurer 1999
logirx1=alog10(10^(0.4*a3*(1.7957/2.4665))-1)+0.076      ; de 1600 a 2800
ym=10^(logirx)
ym1=10^(logirx1)


logirx=alog10(10^(0.4*a2)-1)+0.076      ; Takeuchi 2012
yt=10^(logirx)

;logirx=alog10(10^(0.4*a4)-1)+0.076      ; This work 2015
;yp=10^(logirx)


stop






;------IR por debajo de lim IR a su z-----------------
;fitexy,  zphot(wnok), abs(sfr_w(wnok)), wnoka, wnokb, x_sig=0.1, y_sig=0.1
;low_ir=where(sfr_w lt (zphot*wnokb+wnoka) and sfr_w gt 0) ; 167 con SFR_wuyts menor que su límite IR ; 174



;errors
;------------

err_beta=0.18 ;rel

err_irx=0.18   ;rel
;err_irx=0.91 ;med



print, 'EMPIEZA PLOT'
stop


;====================PLOT=====================================

  set_plot,'PS'
  Device, /Helvetica, /Color, Bits_per_Pixel=8, File=dir+'/synthesizer_offset/files/IRX_Beta_paper_plot.eps', XSize=15, YSize=13


plot, beta, 10^irx, psym=sym(1), symsize=0.2, /ylog, yrange=[0.2,2.e2],xrange=[-2.4,2.4], /ystyle , /xstyle, xtitle='UV slope (!4b!3)', ytitle='IRX = SFR!dIR!N/SFR!dUV!N', thick=4., xthick=4., ythick=4., charsize=1.5, charthick=4.

oploterror, [2.,2.],[0.5,0.5], [err_beta, err_beta], [err_irx*0.5, err_irx*0.5], thick=3.
;oploterror, [2.,2.],[2.,2.], [err_beta, err_beta], [err_irx, err_irx], thick=3.

;original lines
;oplot, x, y2, thick=6., line=3
oplot, x, ym1, col=cgcolor('grn6'), thick=6., line=3.

;used fit
change=closest(x, -0.97)
oplot, x(0:change), ym1(0:change), thick=8., col=2
oplot, x(change:n_elements(x)-1), y2(change:n_elements(x)-1), thick=6., col=2


;oplot, x, yp, thick=6.
;oplot, x, yt, col=4, thick=4., line=2.
;xyouts, .1,.4, 'Takeuchi et al. (2012)', col=4 , charsize=1.2, charthick=3.

;xyouts, 0., .3, 'Meurer et al. (1999)', col=3 , charsize=1.2, charthick=3.
;xyouts, 0.,.5, 'Fit to IR-faint sample', charsize=1.2, charthick=3., col=2

legend, ['This work'], line=0, thick=6., color=2, box=0, pos=[0.45,0.3], charsize=1., charthick=3, /normal
legend, ['Meurer+1999'], line=3, thick=6., color=cgcolor('grn6'), box=0, pos=[0.45,0.25], charsize=1., charthick=3, /normal



oplot, [-0.97,-0.97 ],[0.1, 2000], line=2, thick=2

oplot, beta(low_ir), 10^irx(low_ir), psym=sym(1), col=2,symsize=0.8


 Device, /Close_File
 set_plot,'x'

print, 'FIN plot beta'



stop


END
