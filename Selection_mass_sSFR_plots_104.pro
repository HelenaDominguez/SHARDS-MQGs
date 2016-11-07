pro Selection_mass_sSFR_plots_104

;     .comp /home/hdominguez/idl/pro/new/Selection_mass_sSFR_plots_104.pro
;     .comp /home/hdominguez/idl/sym.pro



dir='/home/hdominguez/SHARDS/DR2_gamma/redanddead_gamma'


;M >= 7
readcol, dir+'/out_HDFN_SHARDS_iDR2gammaP01.20150511_z_1.00_1.50_whole_sample_norep',num, id,RA,DEC,z,q,VJ,UV,mass,sSFR,SFR_UV_obs,SFR_IR, SFR_tot , AV,mag_H, mag_j,a,a,a,sel_string, Format='F, A,F,F,F,F,F,F,f,f,f,f,f,f,f,f,a,a,a,a', skipline=1
; same order than previous files

;Mk=Msalp/1.7

mass_salp=mass
mass_krou=mass_salp+alog10(1/1.7)
mass=mass_krou

m1_all=where(mass lt 10.)  ;3180
m2_all=where(mass ge 10.)  ;477


;selected galaxies
;------------------
readcol, dir+'/out_HDFN_SHARDS_iDR2gammaP01.20150305_z_1.00_1.50_104',num_rd, id_rd,RA_rd,DEC_rd,z_rd,q_rd,VJ_rd,UV_rd,mass_rd,sSFR_rd,SFR_UV_obs_rd,SFR_IR_rd, SFR_tot_rd , AV_rd,mag_H, mag_j,a,a,a,sel_string, Format='F, A,F,F,F,F,F,F,f,f,f,f,f,f,f,f,a,a,a,a', skipline=1
; same order than previous files


mass_rd_salp=mass_rd
mass_rd_krou=mass_rd_salp+alog10(1/1.7)
mass_rd=mass_rd_krou







;----------------selection-----------

UVJ=where(sel_string eq 'UVJ')         ;65
low_SF=where(sel_string eq 'sSFR')     ;39
IR_det_rd=where(sfr_ir_rd gt 0)        ;22   
noIR_det_rd=where(sfr_ir_rd le 0)      ;82
                                
;-----------------------------------
stop

rd_all=fltarr(n_elements(id_rd))   ; son los indices de los subsamples en la muestra entera
uvj_all=fltarr(n_elements(uvj))
lowsf_all=fltarr(n_elements(low_sf))
both=fltarr(n_elements(low_sf))
ir_rd_all=fltarr(n_elements(ir_det_rd))
id_104_all=fltarr(n_elements(id_rd))

ll=0
for ll=0, n_elements(id_rd)-1 do begin &$
rd_all(ll)=where(id eq id_rd(ll))  &$
endfor

for ll=0, n_elements(uvj)-1 do begin &$
uvj_all(ll)=where(id eq id_rd(uvj(ll)))  &$
endfor

ll=0
for ll=0, n_elements(low_sf)-1 do begin &$
lowsf_all(ll)=where(id eq id_rd(low_sf(ll)))  &$
both(ll)=where(id_rd(uvj) eq id_rd(low_sf(ll)))  &$
   endfor

ll=0
for ll=0, n_elements(ir_det_rd)-1 do begin &$
ir_rd_all(ll)=where(id eq id_rd(ir_det_rd(ll)))  &$
endfor


ll=0
for ll=0, n_elements(id_rd)-1 do begin &$
id_104_all(ll)=where(id eq id_rd(ll))  &$
endfor



stop


xxuvj=vj
yyuvj=uv


;=========================SELECCION===================================

;UVJ selected =139 (sin cortes en masa, etc)
;---------------------
rd_uvj=where((yyuvj ge 1.3) and (xxuvj le 1.6) and (yyuvj ge (0.875*(xxuvj)+0.6)))


; UVJ mass  gt 10. = 109
;-----------------
rd_uvj_mass=where((yyuvj ge 1.3) and (xxuvj le 1.6) and (yyuvj ge (0.875*(xxuvj)+0.6)) and mass gt 10)


; UVJ no IR = 111
;-----------------
rd_uvj_no_IR=where((yyuvj ge 1.3) and (xxuvj le 1.6) and (yyuvj ge (0.875*(xxuvj)+0.6))  and sfr_ir le 0)


; UVJ  IR = 28
;-----------------
rd_uvj_IR=where((yyuvj ge 1.3) and (xxuvj le 1.6) and (yyuvj ge (0.875*(xxuvj)+0.6)) and sfr_ir gt 0)


; UVJ IR mass gt 10. = 27
;-----------------
rd_uvj_IR_mass=where((yyuvj ge 1.3) and (xxuvj le 1.6) and (yyuvj ge (0.875*(xxuvj)+0.6)) and sfr_ir gt 0 and mass gt 10.)



; UVJ no IR & mass RAINBOW gt 10. = 82
;-----------------
rd_uvj_no_IR_mass=where((yyuvj ge 1.3) and (xxuvj le 1.6) and (yyuvj ge (0.875*(xxuvj)+0.6)) and sfr_ir le 0 and mass gt 10)

;UVJ mass gt 10., ssfr gt 0.2 =4
;------------------------------
rd_uvj_no_IR_mass_high_ssfr=where((yyuvj ge 1.3) and (xxuvj le 1.6) and (yyuvj ge (0.875*(xxuvj)+0.6)) and sfr_ir le 0 and mass gt 10 and ssfr gt 0.2)


stop
;======================================================

;sSFR lt 0.2 =237
;---------------------
rd_ssfr=where(ssfr lt 0.2)


;sSFR lt 0.2 & IR = 25
;---------------------
rd_ssfr_ir=where(ssfr lt 0.2 and sfr_IR gt 0)


;sSFR lt 0.2 & mass RAINBOW gt 10. =139
;------------------------------
rd_ssfr_mass=where(ssfr lt 0.2 and mass gt 10)


;sSFR lt 0.2 & IR mass gt 10. = 25
;---------------------
rd_ssfr_ir_mass=where(ssfr lt 0.2 and sfr_IR gt 0 and mass gt 0)


;sSFR lt 0.2 & IR mass gt 10. & UVJ = 11  (we recover 11 IR detected in UVJ region)
;---------------------
rd_ssfr_ir_mass_UVJ=where((yyuvj ge 1.3) and (xxuvj le 1.6) and (yyuvj ge (0.875*(xxuvj)+0.6)) and ssfr lt 0.2 and sfr_IR gt 0 and mass gt 0)



;=============== total sample =============
;Las ssfr < 0.2 + las 4 UVJ con ssfr > 0.2 (el resto de UVJ son low SSFR)
total_rd=[rd_ssfr_mass, rd_uvj_no_IR_mass_high_ssfr]    ;143



;flag_uvj=fltarr(n_elements(z))
;flag_uvj(rd_uvj)=1       
;flag_uvj(uvj_all)=1 siempre ok :)


xuvj=fltarr(n_elements(xxuvj))
yuvj=fltarr(n_elements(yyuvj))


for i=0, n_elements(z) -1 do begin    &$

   sigma = 0.03   &$
   meanx = xxuvj(i)   &$
   meany=yyuvj(i)   &$
   arrayx = RANDOMN(seed, 1000, /normal)   &$
   arrayy = RANDOMN(seed, 1000, /normal)   &$
   xuvj(i)= arrayx(0) * sigma + meanx   &$
   yuvj(i)=arrayy(0) * sigma + meany   &$


endfor


stop


;---------------------UVJ----------------------------------

a=[0.8,1.6]
b=[1.3,2.0]

a1=[0., 0.8]
b1=[1.3, 1.3]

a2=[1.6, 1.6]
b2=[2.0, 2.5]

a3=[1., 1.]
b3=[2.5, 1.5]


err_uv=0.09
err_vj=0.08


stop

print, 'EMPIEZA PLOT UVJ'
;----------------------------------------------------------------------------------------------
stop

set_plot,'ps'
device, filename=dir+'/synthesizer_offset/files/selection_UVJ_104_Krou_err.eps', xsize=9, ysize=7, /inch,/col

!p.multi=0
plot, xuvj, yuvj, psym=3, xrange=[0.2,2.2], yrange=[1.0, 2.3], xtitle='V - J ', ytitle='U - V',  /xstyle, /ystyle, thick=3., charthick=4., xthick=5., ythick=5., charsize=2., position=[0.15,0.15,0.92,0.92]



;-----all-------------
oplot, xuvj, yuvj, psym=sym(1),symsize=0.5, thick=3, col=13

;-----------UVJ---------
oplot, xuvj(uvj_all), yuvj(uvj_all), psym=sym(1), col=3,  thick=2. , symsize=1.3  ; corte en masa 

;----IR------
oplot,  xuvj(ir_rd_all), yuvj(ir_rd_all), psym=sym(1), col=4,symsize=2.0, thick=4.
;oplot,  xuvj(ir_rd_all), yuvj(ir_rd_all), psym=sym(2), col=8,symsize=1.2, thick=4.
;oplot,  xuvj(ir_rd_all), yuvj(ir_rd_all), psym=sym(7),symsize=1.2, thick=4.

;----low sSFR--------
oplot, xuvj(lowsf_all), yuvj(lowsf_all), psym=sym(1), col=2,  thick=4., symsize=1.3

;----selection------
oplot, xuvj(id_104_all), yuvj(id_104_all), psym=sym(6), symsize=1.3, thick=2


;----errors----
oploterror, [2.,2.],[1.15,1.15], [err_vj, err_vj], [err_uv, err_uv], thick=3.


oplot, a, b, thick=3.
oplot, a1, b1, thick=3.
oplot, a2, b2, thick=3.
oplot, a3, b3, line=2, thick=3.

legend, ['UVJ & no IR'], psym=sym(1), color=3, box=0, pos=[0.16,0.85], charsize=1.5, charthick=3, /normal
legend, ['sSFR [Gyr!u-1!N] < 0.2'], psym=sym(1), color=2, box=0, pos=[0.16,0.8], charsize=1.5, charthick=3, /normal
legend, ['IR detected MQGs'], psym=sym(1),symsize=2.0, color=4, box=0, pos=[0.16,0.75], charsize=1.5, charthick=3, /normal
legend, [''], psym=sym(1), color=2,symsize=1.3, box=0, pos=[0.16,0.75], charsize=1.5, charthick=3, /normal



xyouts,1.62 ,2.2, '1.0 < z < 1.5', charthick=4, charsize=1.5

plot, xuvj, yuvj, psym=3, xrange=[0.2,2.2], yrange=[1.0, 2.3], xtitle='V - J ', ytitle='U - V',  /xstyle, /ystyle, thick=3., charthick=4., xthick=5., ythick=5., charsize=2., position=[0.15,0.15,0.92,0.92], /nodata, /noerase


 Device, /Close_File
 set_plot,'x'

print, 'FIN plot'


err_mass_med=0.1   ;dex
;err_sSFR_med=0.35  ; relative error
err_sSFR_med=0.03  ; median err 


print, 'EMPIEZA PLOT sSFR'
stop
;-------------------------sSFR vs Mass---------------------------------

set_plot,'ps'
device, filename=dir+'/synthesizer_offset/files/selection_sSFR_104_Krou_err_med.eps', xsize=9, ysize=7, /inch,/col

!p.multi=0

; total 2800
;------------

plot,  10^(mass), ssfr, psym=1, yrange=[0.008, 3.],xrange=[10^(9.7), 1.e12], ytitle='sSFR [Gyr!u-1!N]', xtitle='M [M!9!dn!3!n]', thick=2., xthick=5., ythick=5. , /xstyle, /ystyle, /ylog, /xlog,charsize=2.0,  charthick=4,/nodata, position=[0.15,0.15,0.92,0.92]


ok_plot=where(mass lt 10. or ssfr ge 0.2)

oplot,  10^(mass(ok_plot)), ssfr(ok_plot), col=13,  psym=sym(1),symsize=0.5

;-----UVJ------
oplot,  10^(mass(uvj_all)), ssfr(uvj_all), psym=sym(1), col=3,  thick=2., symsize=1.3

;----IR------
oplot,  10^(mass(ir_rd_all)), ssfr(ir_rd_all), psym=sym(1), col=4,symsize=2.,  thick=4.
;oplot,  10^(mass(ir_rd_all)), ssfr(ir_rd_all), psym=sym(2), col=8,symsize=1.2,  thick=4.
;oplot,  10^(mass(ir_rd_all)), ssfr(ir_rd_all), psym=sym(7), symsize=1.2,  thick=4.


;----low sSFR--------
oplot,  10^(mass(lowsf_all)), ssfr(lowsf_all), psym=sym(1), col=2,  thick=2. , symsize=1.3


;----Final sample (117)
oplot,10^mass(id_104_all), ssfr(id_104_all), psym=sym(6), symsize=1.3, thick=2

;----errors----
;cuentas el 6/11/15
;oploterror, [10^11.7,10^11.7],[0.02,0.02], [err_mass_med*10^11.7*2.3, err_mass_med*10^11.7*2.3], [err_sSFR_med*0.02, err_sSFR_med*0.02], thick=3.
oploterror, [10^11.7,10^11.7],[.06,.06], [err_mass_med*10^11.7*2.3, err_mass_med*10^11.7*2.3], [0.03, 0.03], thick=3.

xx=[0.01,1.e13]
yy=[0.2,0.2]

oplot, xx, yy, line=2, thick=3.

xyouts,0.7, 0.85,/norm, '1.0 < z < 1.5', charthick=4, charsize=1.5

 Device, /Close_File
 set_plot,'x'

print, 'FIN plot'



stop
stop




;MAL!! ;!!COMPRUEBA SSFR!!
;=========================================================
; total 2800
;------------

plot,  10^(mass), sfr_total/10^(mass)*1.e9, psym=1, yrange=[0.01, 100.],xrange=[1.e8, 1.e13], ytitle='sSFR(total, 2800) Gyr-1', xtitle='Mass (Msun)', thick=2., xthick=2., ythick=2. , /xstyle, /ystyle, /ylog, /xlog,charsize=1., /nodata

;oplot,  10^(mass), sfr_total/10^(mass)*1.e9, psym=3,col=13
oplot,  10^(mass), ssfr_total, psym=3,col=13

;!!COMPRUEBA SSFR!!

oplot,  10^(mass_rd(uvj)), ssfr_rd(uvj), psym=sym(1), col=3
oplot,  10^(mass_rd(low_sf)), ssfr_rd(low_sf), psym=sym(1), col=2
oplot,  10^(mass_rd(ir_det_rd)), ssfr_rd(ir_det_rd), psym=sym(5), col=5,symsize=0.5
oplot,  10^(mass_rd), ssfr_rd, psym=sym(6), col=8, symsize=1.3


xx=[0.01,1.e13]
yy=[0.2,0.2]

oplot, xx, yy, line=2




stop



END
