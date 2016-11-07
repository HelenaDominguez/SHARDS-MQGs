PRO SFH_tracks

; .comp /home/hdominguez/idl/pro/new/SFH_tracks.pro

code='BC'
model='txp'
imf='krou'
met_string='_met'
;met_string=''
nsim=300

Msun='M'+sunsymbol()



dir='/home/hdominguez/SHARDS/DR2_gamma/redanddead_gamma'

filelist=dir+'/ID_SHARDS_candels_zphotoiris_spectra_files_104.lst'
readcol, filelist, source, zcomb, G141_file, G102_file, format='a, f, a, a'

id_list=source


readcol, dir+'/synthesizer_offset/files/analyze_parameters_BC_txp_krou_flags_104.lst',   shards_id , nclus , z ,  log_tau,  log_age, log_age_MW,  av , met, mass,ratio, flag_d4, flag_mg, flag_deg, flag_sol, format='a,i, f, f, f'

age_gy=10^log_age*1.e-9
age_mw_gy=10^log_age_mw*1.e-9
tau_my=10^log_tau*1.e-6



restore, dir+'/synthesizer_offset/files/index_BC03_txp_krou_104_d4_err.sav', /verbose


;selected galaxies
;------------------
readcol, dir+'/out_HDFN_SHARDS_iDR2gammaP01.20150305_z_1.00_1.50_104',num_rd, id_rd,RA_rd,DEC_rd,z_rd,q_rd,VJ_rd,UV_rd,mass_rd,sSFR_rd,SFR_UV_obs_rd,SFR_IR_rd, SFR_tot_rd , AV_rd,mag_H, mag_j,a,a,a,sel_string, Format='F, A,F,F,F,F,F,F,f,f,f,f,f,f,f,f,a,a,a,a', skipline=1
; same order than previous files


age_model=fltarr(n_elements(id_list), 20000)
sfr_model=fltarr(n_elements(id_list), 20000)
mass_model=fltarr(n_elements(id_list), 20000)
lbt_model=fltarr(n_elements(id_list), 20000)
ssfr_model=fltarr(n_elements(id_list), 20000)


age_bf=fltarr(n_elements(id_list))  ;Gyr
tau_bf=fltarr(n_elements(id_list)) ;Myr
mass_bf=fltarr(n_elements(id_list)) ;Msun
log_tau_bf=fltarr(n_elements(id_list))
sfr_model_max=fltarr(n_elements(id_list))
tq=fltarr(n_elements(id_list))

now=fltarr(n_elements(id_list))
p1=fltarr(n_elements(id_list))
p2=fltarr(n_elements(id_list))
p3=fltarr(n_elements(id_list))

a11=fltarr(n_elements(id_list))
a1=fltarr(n_elements(id_list))
a2=fltarr(n_elements(id_list))
a3=fltarr(n_elements(id_list))
a4=fltarr(n_elements(id_list))
a5=fltarr(n_elements(id_list))
a6=fltarr(n_elements(id_list))


for j=0, n_elements(id_list) -1 do begin    &$
;for j=123, 134 do begin 

id=id_list(j)   &$
print, j, id

w=where(shards_id eq id and flag_sol eq 1 )   ;choose good solution

if w(0) ne -1 then begin
num=nclus(w)+1
if code eq 'BC' then  file_sfh=dir+'/synthesizer_offset/res/analyze_rejoin5/clusters/'+id+'.BCall.bc2003.stelib.'+model+'.'+imf+'.sfh.clust'+strtrim(num,1)+'.res'   &$


readcol,file_sfh,age_sfh ,SFR_old , lmassold, time_you ,SFR_you, lmassyou, SFR_sfh, mass_sfh, lbt
readcol,file_sfh,x,x,x,x,x,x,x,x,x, age,x,tau, format='a,a,a,a,a,a,a,a,a,a,a,a,a', numline=1

age_bf(j)=float(age)
tau_bf(j)=float(tau)
log_tau_bf(j)=alog10(tau_bf(j)*1.e6)

age_model(j,*)=age_sfh
sfr_model(j,*)=sfr_sfh
mass_model(j,*)=mass_sfh
lbt_model(j,*)=lbt
sfr_model_max(j)=max(sfr_sfh)
ssfr_sfh=(sfr_sfh/10^mass_sfh)*1.e9
ssfr_model(j,*)=ssfr_sfh




now(j)=closest(lbt, 0)  
if age_bf(j) ge 2.0 then p1(j)=closest(lbt, -2.0) else p1(j)=-1 ; 2Gyr before obs
if age_bf(j) ge 1.0 then p2(j)=closest(lbt, -1.0) else p2(j)=-1 ; 1 Gyr before obs
p3(j)=closest(lbt, 1.) ; 1 Gyr after obs


mass_bf(j)=mass_sfh(now(j))
;ssfr_bf(j)=ssfr_model(j,now(j))

q=closest(ssfr_sfh, 0.2)    ; sSFR , 0.2
tq(j)=age_sfh(now(j))-age_sfh(q)  ; time since quenching en Gyr



a1(j)=closest(age_sfh, 0.1)  ;100
a2(j)=closest(age_sfh, 0.5)   ;500
a3(j)=closest(age_sfh, 1.0)   ;1000
a4(j)=closest(age_sfh, 2.0)   ;2000
a5(j)=closest(age_sfh, 3.0)   ;3000
a6(j)=closest(age_sfh, 5.0)   ;5000

endif

endfor


stop

;-----------------Different epochs (always same position in sfh files)----------------
a1=100  ; 0.1 Gyr
a2=500  ; 0.5
a3=1000 ; 1.
a4=2000 ; 2.
a5=3000 ; 3.
a6=4000 ; 4.
a7=5000 ; 5.
;------------------------------------------------------------


;==============Definitions to remove model index=====

sfr_01=fltarr(n_elements(id_list))
sfr_now=fltarr(n_elements(id_list))

age_cl=fltarr(n_elements(id_list))
log_tau_cl=fltarr(n_elements(id_list))
mass_cl=fltarr(n_elements(id_list))
ngal_cl=fltarr(n_elements(id_list))

for l=0, n_elements(id_list)-1 do begin &$

sfr_01(l)=sfr_model(l,100) &$ 
sfr_now(l)=sfr_model(l,now(l)) &$



endfor



;---limites en SFR < 10⁻5----
;-------------------------------------------

      

limit=where(alog10(sfr_now) lt -3. )   

arrayx = RANDOMN(seed, n_elements(limit), /normal)
sfr_now(limit)=10^(-2.5+arrayx*0.1)




; MS Rodigiero 2011 (Salpeter IMF)
;-------------------
;logSFR=a*logM+b
a_ms=0.79
b_ms=-6.42




;MS Speagle 2014 (Kroupa IMF)
;-----------------

mass_ms=findgen(100)/10.+9

tz00=6    ;z=0.9 ;+1 Gyr  ;p3
tz0=5     ;z=1.2 ;0 Gyr   ;now
tz2=4     ;z=1.6 ;-1 Gyr  ;p2
tz1=3     ;z=2.1 ;-2 Gyr  ;p1


sfr_ms00=(0.84-0.026*tz00)*mass_ms-(6.51-0.11*tz00)
sfr_ms0=(0.84-0.026*tz0)*mass_ms-(6.51-0.11*tz0)
sfr_ms1=(0.84-0.026*tz1)*mass_ms-(6.51-0.11*tz1)
sfr_ms2=(0.84-0.026*tz2)*mass_ms-(6.51-0.11*tz2)



;MS Elbaz (Salpeter)
;---------------
;SFR=7.2*[M/(10^10)]^0.9
;SFR_kr=0.7*SFR_salp

mass_mse=10^(mass_ms-10)
SFR_mse=0.7*(7.2*mass_mse^0.9)



;colorin=[cgcolor('Sky Blue'),cgcolor('Royal Blue'), cgcolor('teal'), cgcolor('Gold'),cgcolor('Orange'),cgcolor('Crimson')]

colorin=[cgcolor('Sky Blue'),cgcolor('Blue'), cgcolor('Grn5'),cgcolor('Orange'),cgcolor('Crimson')]



;===================PLOT AGES==================================

w=where(age_bf gt 0.)
age_vect=[0.1, 0.5, 1., 2., 4.]
tmp=fltarr(n_elements(id_list))



print, 'Empieza Plot tracks'
stop

set_plot,'ps'
device, filename=dir+'/synthesizer_offset/files/SFH_tracks_'+code+'_'+model+'_'+imf+''+met_string+'_only.eps',/col, /portrait, /times


;===============tracks===================================================

!p.multi=0

plot, mass_model, sfr_model,  xrange=[9.0, 11.5],yrange=[1.e-3, 1.e4], /ylog, xtitle='log (M/'+Msun+')', ytitle='SFR!dSED!N ['+Msun+'/yr]',thick=2., xthick=4., ythick=4., charthick=3., charsize=1.2, psym=1, /nodata, /xstyle, /ystyle, position=[0.12,0.12,0.95,0.95]

;-----whole tracks-----
for i =0, n_elements(w)-1 do begin &$
oplot, mass_model(w(i),0:now(w(i))), sfr_model(w(i),0:now(w(i))), col=13, thick=1., line=0  &$
;oplot, mass_model(w(i),now(w(i)):19999), sfr_model(w(i),now(w(i)):19999), col=13, line=2, thick=1.  &$
endfor




;-----MS-------


cgColorFill, [9.,9.,11.5, 11.5, 9.], [10^(sfr_ms1(0)-0.2),10^(sfr_ms1(0)+0.2) ,10^(sfr_ms1(25)+0.2) ,10^(sfr_ms1(25)-0.2),10^(sfr_ms1(0)-0.2)],  COL=cgcolor('YGB2')


oplot, mass_ms, 10^sfr_ms0,  thick=9., col=cgcolor('Dark Gray'),line=3           ;z=1.2
oplot, mass_ms, 10^sfr_ms1,  thick=9., col=cgcolor('Dark Gray')  ;z=2.0
oplot, mass_ms, SFR_mse,  thick=9., col=cgcolor('Slate Gray'), line=2    ;z=1.0


;-----points diferrent ages------
for i =0, n_elements(w)-1 do begin &$
if now(w(i)) gt a1 then oplot, [mass_model(w(i),a1),mass_model(w(i),a1)],  [sfr_model(w(i),a1),sfr_model(w(i),a1)], psym=sym(1), symsize=0.6, col=colorin(0), thick=4. &$
if now(w(i)) gt a2 then oplot, [mass_model(w(i),a2),mass_model(w(i),a2)],  [sfr_model(w(i),a2),sfr_model(w(i),a2)], psym=sym(1), symsize=0.6, col=colorin(1), thick=4. &$
if now(w(i)) gt a3 then oplot, [mass_model(w(i),a3),mass_model(w(i),a3)],  [sfr_model(w(i),a3),sfr_model(w(i),a3)], psym=sym(1), symsize=0.6, col=colorin(2), thick=4. &$
if now(w(i)) gt a4 then oplot, [mass_model(w(i),a4),mass_model(w(i),a4)],  [sfr_model(w(i),a4),sfr_model(w(i),a4)], psym=sym(1), symsize=0.6, col=colorin(3), thick=4. &$
;if now(w(i)) gt a5 then oplot, [mass_model(w(i),a5),mass_model(w(i),a5)],  [sfr_model(w(i),a5),sfr_model(w(i),a5)], psym=sym(6), symsize=1., col=colorin(4), thick=2. &$
if now(w(i)) gt a6 then oplot, [mass_model(w(i),a6),mass_model(w(i),a6)],  [sfr_model(w(i),a6),sfr_model(w(i),a6)], psym=sym(1), symsize=0.6, col=colorin(4), thick=4. &$

  endfor





xyouts,10.3,10^3.6 , 'MS @ z=2.0 Speagle+14', charthick=2.5, col=cgcolor('dark Grey'), charsize=1.
xyouts,10.3, 10^3.3 ,'MS @ z=1.2 Speagle+14', charthick=2.5, col=cgcolor('dark Grey'), charsize=1.
xyouts,10.3, 10^3.0 ,'MS @ z=1.0 Elbaz+07', charthick=2.5, col=cgcolor('Dark Grey'), charsize=1.

oplot, [11.15,11.4],[10^3.6,10^3.6], thick=5., col=cgcolor('Dark Grey')
oplot, [11.15,11.4],[10^3.3,10^3.3], thick=5.,col=cgcolor('Dark Grey'), line=3
oplot, [11.15,11.4],[10^3.0,10^3.0], thick=5., col=cgcolor('Slate Grey') ,line=2


xyouts,0.25 ,0.45,/norm ,'t!d0!N=0.1 Gyr', col=colorin(0), charthick=3.5
xyouts,0.25 , 0.40,/norm,'t!d0!N=0.5 Gyr', col=colorin(1), charthick=3.5
xyouts,0.25 , 0.35,/norm,'t!d0!N=1.0 Gyr', col=colorin(2), charthick=3.5
xyouts,0.25 , 0.30,/norm,'t!d0!N=2.0 Gyr', col=colorin(3), charthick=3.5
xyouts,0.25 , 0.25,/norm,'t!d0!N=4.0 Gyr', col=colorin(4), charthick=3.5


plot, mass_model, sfr_model,  xrange=[9.0, 11.5],yrange=[1.e-3, 1.e4], /ylog, xtitle='log (M/'+Msun+')', ytitle='SFR!dSED!N ['+Msun+'/yr]',thick=2., xthick=4., ythick=4., charthick=3., charsize=1.2, psym=1, /nodata, /xstyle, /ystyle, position=[0.12,0.12,0.95,0.95], /noerase




device,/close & set_plot,'x'

stop



;--------------colors by age bin--------------------

;for i=0, n_elements(id_list)-1 do tmp(i)=closest(age_vect, age_gy(i))

;w0=where(tmp eq 0)
;w1=where(tmp eq 1)
;w2=where(tmp eq 2)
;w3=where(tmp eq 3)
;w4=where(tmp eq 4)


;w0=where(age_mw_gy le 0.5)
;w1=where(age_mw_gy gt  0.5 and age_mw_gy le 1.0)
;w2=where(age_mw_gy gt  1.0 and age_mw_gy le 2.0)
;w3=where(age_mw_gy gt  2.0 and age_mw_gy le 4.0)
;w4=where(age_mw_gy gt  4.0)


w0=where(age_gy le 0.5)
w1=where(age_gy gt  0.5 and age_gy le 1.0)
w2=where(age_gy gt  1.0 and age_gy le 2.0)
w3=where(age_gy gt  2.0 and age_gy le 4.0)
w4=where(age_gy gt  4.0)

restore, dir+'/synthesizer_offset/files/BestFit_SFH_BC_txp_krou_met_104_l100.sav', /verbose
;% RESTORE: Restored variable: ID_LIST.
;% RESTORE: Restored variable: AGE_BF.
;% RESTORE: Restored variable: TAU_BF.
;% RESTORE: Restored variable: MASS_BF.
;% RESTORE: Restored variable: SFR_BF.
;% RESTORE: Restored variable: SSFR_BF.
;% RESTORE: Restored variable: TQ.
;% RESTORE: Restored variable: T_LIRG.
;% RESTORE: Restored variable: T_ULIRG.
;% RESTORE: Restored variable: F_LIRG.
;% RESTORE: Restored variable: F_ULIRG.
;% RESTORE: Restored variable: MASS_LGYR.
;% RESTORE: Restored variable: SFR_L100.
;% RESTORE: Restored variable: SSFR_L100.

stop

sfr_now_original=sfr_now
;sfr_now=sfr_tot_rd
sfr_now=sfr_l100

limit=where(alog10(sfr_now) lt -3. )   

arrayx = RANDOMN(seed, n_elements(limit), /normal)
sfr_now(limit)=10^(-2.5+arrayx*0.1)


print, 'Empieza Plot observed'
stop

set_plot,'ps'
device, filename=dir+'/synthesizer_offset/files/SFH_tracks_'+code+'_'+model+'_'+imf+''+met_string+'_obs_SFR_SED_l100.eps',/col, /portrait, /times


colours
;===============tracks===================================================

!p.multi=0

plot, mass_model, sfr_model,  xrange=[9.0, 11.5],yrange=[1.e-3, 1.e4], /ylog, xtitle='log (M/'+Msun+')', ytitle='SFR!dSED!N ['+Msun+'/yr]',thick=2., xthick=4., ythick=4., charthick=3., charsize=1.2, psym=1, /nodata, /xstyle, /ystyle, position=[0.12,0.12,0.95,0.95]

;-----whole tracks-----
for i =0, n_elements(w)-1 do begin &$
oplot, mass_model(w(i),0:now(w(i))), sfr_model(w(i),0:now(w(i))), col=13, thick=1., line=0  &$
;oplot, mass_model(w(i),now(w(i)):19999), sfr_model(w(i),now(w(i)):19999), col=13, line=2, thick=3.  &$
endfor


;-----MS-------


cgColorFill, [9.,9.,11.5, 11.5, 9.], [10^(sfr_ms1(0)-0.2),10^(sfr_ms1(0)+0.2) ,10^(sfr_ms1(25)+0.2) ,10^(sfr_ms1(25)-0.2),10^(sfr_ms1(0)-0.2)],  COL=cgcolor('YGB2')

oplot, mass_ms, 10^sfr_ms0,  thick=9., col=cgcolor('Dark Gray'),line=3           ;z=1.2
oplot, mass_ms, 10^sfr_ms1,  thick=9., col=cgcolor('Dark Gray')  ;z=2.0
oplot, mass_ms, SFR_mse,  thick=9., col=cgcolor('Slate Gray'), line=2    ;z=1.0


;----observed------

;for i =0, n_elements(w0)-1 do begin &$
oplot, mass_bf(w0), sfr_now(w0), psym=sym(1),thick=4., col=colorin(0), symsize=1.2  &$
oplot, mass_bf(w0), sfr_now(w0), psym=sym(6),thick=4, symsize=1.2  &$
;endfor


;for i =0, n_elements(w1)-1 do begin &$
oplot, mass_bf(w1), sfr_now(w1), psym=sym(1),thick=4., col=colorin(1), symsize=1.2  &$
oplot, mass_bf(w1), sfr_now(w1), psym=sym(6),thick=4, symsize=1.2  &$
;endfor


;for i =0, n_elements(w2)-1 do begin &$
oplot, mass_bf(w2), sfr_now(w2), psym=sym(1),thick=4., col=colorin(2), symsize=1.2  &$
oplot, mass_bf(w2), sfr_now(w2), psym=sym(6),thick=4, symsize=1.2  &$
;endfor


;for i =0, n_elements(w3)-1 do begin &$
oplot, mass_bf(w3), sfr_now(w3), psym=sym(1),thick=4., col=colorin(3), symsize=1.2  &$
oplot, mass_bf(w3), sfr_now(w3), psym=sym(6),thick=4, symsize=1.2  &$
;endfor


;for i =0, n_elements(w4)-1 do begin &$
oplot, mass_bf(w4), sfr_now(w4), psym=sym(1),thick=4., col=colorin(4), symsize=1.2  &$
oplot, mass_bf(w4), sfr_now(w4), psym=sym(6),thick=4, symsize=1.2  &$
;endfor




oplot, [9.1, 9.1], [1.e3, 1.e3] ,psym=sym(1), thick=4
xyouts,9.22 ,10^(2.9) ,'Obs. Galaxies, SFR!dSED!N', charthick=2.5
;xyouts,9.22 ,10^(2.9) ,'Obs. Galaxies, SFR!dUV!N', charthick=2.5

xyouts,0.2 ,0.45,/norm ,'t!d0!N [Gyr] < 0.5', col=colorin(0), charthick=3.5
xyouts,0.2 , 0.40,/norm,'0.5 < t!d0!N [Gyr] < 1.0', col=colorin(1), charthick=3.5
xyouts,0.2 , 0.35,/norm,'1.0 < t!d0!N [Gyr] < 2.0', col=colorin(2), charthick=3.5
xyouts,0.2 , 0.30,/norm,'2.0 < t!d0!N [Gyr] < 4.0', col=colorin(3), charthick=3.5
xyouts,0.2 , 0.25,/norm,'t!d0!N [Gyr] > 4.0', col=colorin(4), charthick=3.5


plot, mass_model, sfr_model,  xrange=[9.0, 11.5],yrange=[1.e-3, 1.e4], /ylog, xtitle='log (M/'+Msun+')', ytitle='SFR!dSED!N ['+Msun+'/yr]',thick=2., xthick=4., ythick=4., charthick=3., charsize=1.2, psym=1, /nodata, /xstyle, /ystyle, position=[0.12,0.12,0.95,0.95], /noerase


device,/close & set_plot,'x'

stop





;===================PLOT ALL==================================


stop
set_plot,'ps'
device, filename=dir+'/synthesizer/plots_mine/SFH_tracks_'+code+'_'+model+'_'+imf+''+met_string+'_colors1.eps',/col, /portrait, /times


;===============tracks===================================================

!p.multi=0

plot, mass_model, sfr_model,  xrange=[9.0, 11.7],yrange=[1.e-3, 5000], /ylog, xtitle='log (M/'+Msun+')', ytitle='SFR ('+Msun+'/yr)',thick=2., xthick=4., ythick=4., charthick=3., charsize=1.5, psym=1, /nodata, /xstyle


;for i =0, n_elements(id_list)-1 do begin &$
;oplot, mass_model(i,0:now(i)), sfr_model(i,0:now(i)), col=13, thick=1., line=0  &$
;oplot, mass_model(i,now(i):19999), sfr_model(i,now(i):19999), col=13, line=1, thick=1.  &$
;endfor


for i =0, n_elements(id_list)-1 do begin &$
if now(i) gt a1 then oplot, [mass_model(i,a1),mass_model(i,a1)],  [sfr_model(i,a1),sfr_model(i,a1)], psym=sym(1), symsize=1., col=colorin(0), thick=2. &$
if now(i) gt a2 then oplot, [mass_model(i,a2),mass_model(i,a2)],  [sfr_model(i,a2),sfr_model(i,a2)], psym=sym(1), symsize=1., col=colorin(1), thick=2. &$
if now(i) gt a3 then oplot, [mass_model(i,a3),mass_model(i,a3)],  [sfr_model(i,a3),sfr_model(i,a3)], psym=sym(1), symsize=1., col=colorin(2), thick=2. &$
if now(i) gt a4 then oplot, [mass_model(i,a4),mass_model(i,a4)],  [sfr_model(i,a4),sfr_model(i,a4)], psym=sym(1), symsize=1., col=colorin(3), thick=2. &$
;if now(i) gt a5 then oplot, [mass_model(i,a5),mass_model(i,a5)],  [sfr_model(i,a5),sfr_model(i,a5)], psym=sym(1), symsize=1., col=colorin(4), thick=2. &$
if now(i) gt a6 then oplot, [mass_model(i,a6),mass_model(i,a6)],  [sfr_model(i,a6),sfr_model(i,a6)], psym=sym(1), symsize=1., col=colorin(5), thick=2. &$

endfor


oplot, mass_model(*,a1), sfr_model(*,a1), psym=8, symsize=0.6, col=colorin(0), thick=2. &$
oplot, mass_model(*,a2), sfr_model(*,a2), psym=8, symsize=0.6, col=colorin(1), thick=2. &$
oplot, mass_model(*,a3), sfr_model(*,a3), psym=8, symsize=0.6, col=colorin(2), thick=2. &$
oplot, mass_model(*,a4), sfr_model(*,a4), psym=8, symsize=0.6, col=colorin(3), thick=2. &$
;oplot, mass_model(*,a5), sfr_model(*,a5), psym=8, symsize=0.6, col=colorin(4), thick=2. &$
oplot, mass_model(*,a6), sfr_model(*,a6), psym=8, symsize=0.6, col=colorin(5), thick=2. &$

oplot, alog10(mass_bf), sfr_now, psym=4,thick=4.


oplot, mass_ms, 10^sfr_ms0,  thick=4.
oplot, mass_ms, 10^sfr_ms1,  thick=1.

xyouts,0.8 ,0.8,/norm ,'MS @ z=2.0', charthick=2.5
xyouts,0.8 ,0.7,/norm ,'MS @ z=1.2', charthick=2.5

oplot, [9.1, 9.1], [1., 1.] ,psym=4, thick=4
;xyouts,9.22 ,.99 ,'Obs. Galaxies', charthick=2.5


xyouts,0.2 ,0.45,/norm ,'0.1 Gyr', col=colorin(0), charthick=3.5
xyouts,0.2 , 0.40,/norm,'0.5 Gyr', col=colorin(1), charthick=3.5
xyouts,0.2 , 0.35,/norm,'1.0 Gyr', col=colorin(2), charthick=3.5
xyouts,0.2 , 0.30,/norm,'2.0 Gyr', col=colorin(3), charthick=3.5
xyouts,0.2 , 0.25,/norm,'4.0 Gyr', col=colorin(4), charthick=3.5
;xyouts,0.2, 0.2,/norm,'4.0 Gyr', col=colorin(5), charthick=3.5

device,/close & set_plot,'x'



print, 'FIN PLOT'
stop


print, 'PLOT age vs Age-MW'
stop


set_plot,'ps'
device, filename=dir+'/synthesizer_offset/files/Age_vs_Age_MW_104.eps',/col, /portrait, /times

plot, age_gy, age_mw_gy, psym=1, xtitle='Age (Gyr)', ytitle='Age!dM!N (Gyr)',thick=2., xthick=4., ythick=4., charthick=3., charsize=1.5
oplot, findgen(10), findgen(10)

device,/close & set_plot,'x'



plot, alog10(mass_bf), sfr_now, psym=1 ,  xrange=[9.8, 11.7],yrange=[1.e-8, 5000], /ylog

nosfr=where(sfr_now lt 1.e-6)
help, nosfr
print, log_tau_bf(nosfr)
print, alog10(age_bf(nosfr)*1.e9)

stop





END
