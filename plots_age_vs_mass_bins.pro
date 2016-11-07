PRO plots_age_vs_mass_bins

;    .comp /home/hdominguez/idl/pro/new/plots_age_vs_mass_bins.pro

code='BC'
model='txp'
imf='krou'
;met_string='_met'
met_string=''
nsim=300

Msun='M'+sunsymbol()
Zsun='Z'+sunsymbol()

dir='/home/hdominguez/SHARDS/DR2_gamma/redanddead_gamma'

filelist=dir+'/ID_SHARDS_candels_zphotoiris_spectra_files_104.lst'
readcol, filelist, source, zcomb, G141_file, G102_file, format='a, f, a, a'


;Sol  indices
readcol, dir+'/synthesizer_offset/files/analyze_parameters_BC_txp_krou_flags_104.lst',   shards_id , nclus , z ,  log_tau,  log_age, log_age_MW,  av , met, mass,ratio, flag_d4, flag_mg, flag_deg, flag_sol, format='a, f, f, f'

;Most sign.
;readcol, dir+'/synthesizer_offset/files/analyze_parameters_BC_txp_krou_most_significant_sol_104.lst',   shards_id , nclus , z ,  log_tau,  log_age, log_age_MW,  av , met, mass,ratio, format='a, f, f, f'

;SSP solutions
;readcol, dir+'/synthesizer_offset/files/analyze_parameters_SSP_txp_krou_101_MostSig.lst',   shards_id , nclus , z ,  log_tau,  log_age, log_age_MW,  av , met, mass,ratio, format='a, f, f, f'


age_gy=10^log_age*1.e-9
age_mw_gy=10^log_age_mw*1.e-9
tau_my=10^log_tau*1.e-6



restore, dir+'/synthesizer_offset/files/index_BC03_txp_krou_104_d4_err.sav', /verbose


; Form analyze_errors
;-----------------
;err_tau=40.1550
;err_age=0.375800
;err_av=0.122500
;err_mass=0.106000


; Form analyze_errors
;-----------------
err_tau=0.16     ;relative
;err_age=0.12     ;relative

;err_tau=21       ;med
err_age=0.21     ;med

err_av=0.1
err_mass=0.05


;Primary solutions
;===================================
; By the moment select only the most probable

;ps=where(flag_sol ge 1)
ok=where(mass ge 10.)


onesol=where(ratio eq 100.)  ;48
degsol=where(ratio lt 100.)  ;56






stop

;=============
;quartiles
;=============

qage=quartile(age_mw_gy)
qmass=quartile(mass)
qtau=quartile(tau_my)
qav=quartile(av)
qmet=quartile(met)


qage_ps=quartile(age_mw_gy(ok))
qmass_ps=quartile(mass(ok))
qtau_ps=quartile(tau_my(ok))
qav_ps=quartile(av(ok))
qmet_ps=quartile(met(ok))


;OLD vs YOUNG
;==================

;old=where(age_mw_gy(ps(ok)) ge 2.0)    ;17 ~ 16%
;young=where(age_mw_gy(ps(ok)) lt 2.0)  ;87

old=where(age_mw_gy ge 2.0)    ;16 ~ 16%
young=where(age_mw_gy lt 2.0)  ;88

qmass1=quartile(mass(young))
qmass2=quartile(mass(old))

qtau1=quartile(tau_my(young))
qtau2=quartile(tau_my(old))


qav1=quartile(av(young))
qav2=quartile(av(old))


qmet1=quartile(met(young))
qmet2=quartile(met(old))


stop


;==========================================
;                PLOT
;==========================================


x=mass
y=age_mw_gy
;y=age_gy


xmin=9.9
xmax=11.5
ymin=0
ymax=4.5
;ymax=5.5


wtau1=where(tau_my lt 80.)
wtau2=where(tau_my ge 80. and tau_my lt 400.)
wtau3=where(tau_my ge 400. )


wtau1_ps=where(tau_my(onesol)  lt 80.)
wtau2_ps=where(tau_my(onesol)  ge 80. and tau_my(onesol)  lt 400.)
wtau3_ps=where(tau_my(onesol)  ge 400. )

wtau1_deg=where(tau_my(degsol)  lt 80.)
wtau2_deg=where(tau_my(degsol)  ge 80. and tau_my(degsol)  lt 400.)
wtau3_deg=where(tau_my(degsol)  ge 400. )


;-----median and scatter in each bin-----

med_mass1=median(mass(wtau1))
med_age1=median(age_mw_gy(wtau1))

med_mass2=median(mass(wtau2))
med_age2=median(age_mw_gy(wtau2))

med_mass3=median(mass(wtau3))
med_age3=median(age_mw_gy(wtau3))

stv_mass1=stdev(mass(wtau1))
stv_age1=stdev(age_mw_gy(wtau1))

stv_mass2=stdev(mass(wtau2))
stv_age2=stdev(age_mw_gy(wtau2))

stv_mass3=stdev(mass(wtau3))
stv_age3=stdev(age_mw_gy(wtau3))

;-------------

med_mass=[med_mass1,med_mass2,med_mass3]
stv_mass=[stv_mass1,stv_mass2,stv_mass3]

med_age=[med_age1,med_age2,med_age3]
stv_age=[stv_age1,stv_age2,stv_age3]


; error bars including 68% of data
;-----------------------------

;=====age========
num1=fix(n_elements(wtau1)*0.136)  ;1sigma lower lim
num2=fix(n_elements(wtau1)*0.842)  ;1 sigma upper lim

ord1=sort(log_age_mw(wtau1))
age_ord1=log_age_mw(wtau1(ord1))

age_inf1=10^age_ord1(num1)*1.e-9
age_sup1=10^age_ord1(num2)*1.e-9

;---
num1=fix(n_elements(wtau2)*0.136)  ;1sigma lower lim
num2=fix(n_elements(wtau2)*0.842)  ;1 sigma upper lim

ord2=sort(log_age_mw(wtau2))
age_ord2=log_age_mw(wtau2(ord2))

age_inf2=10^age_ord2(num1)*1.e-9
age_sup2=10^age_ord2(num2)*1.e-9

;---
num1=fix(n_elements(wtau3)*0.136)  ;1sigma lower lim
num2=fix(n_elements(wtau3)*0.842)  ;1 sigma upper lim

ord3=sort(log_age_mw(wtau3))
age_ord3=log_age_mw(wtau3(ord3))

age_inf3=10^age_ord3(num1)*1.e-9
age_sup3=10^age_ord3(num2)*1.e-9



;=====Mass========
num1=fix(n_elements(wtau1)*0.136)  ;1sigma lower lim
num2=fix(n_elements(wtau1)*0.842)  ;1 sigma upper lim

ord1=sort(mass(wtau1))
mass_ord1=mass(wtau1(ord1))

mass_inf1=mass_ord1(num1)
mass_sup1=mass_ord1(num2)

;----
num1=fix(n_elements(wtau2)*0.136)  ;1sigma lower lim
num2=fix(n_elements(wtau2)*0.842)  ;1 sigma upper lim

ord2=sort(mass(wtau2))
mass_ord2=mass(wtau2(ord2))

mass_inf2=mass_ord2(num1)
mass_sup2=mass_ord2(num2)

;----
num1=fix(n_elements(wtau3)*0.136)  ;1sigma lower lim
num2=fix(n_elements(wtau3)*0.842)  ;1 sigma upper lim

ord3=sort(mass(wtau3))
mass_ord3=mass(wtau3(ord3))

mass_inf3=mass_ord3(num1)
mass_sup3=mass_ord3(num2)


print, 'empieza PLOT Tau colors'

stop
;=============================================


psfile=strtrim(dir+'/synthesizer_offset/files/AgeMW_vs_Mass_tau_bins_'+code+'_'+model+'_'+imf+''+met_string+'_PrimarySolutions.eps',2)

;psfile=strtrim(dir+'/synthesizer_offset/files/AgeMW_vs_Mass_tau_bins_'+code+'_'+model+'_'+imf+''+met_string+'_SSP.eps',2)

;psfile=strtrim(dir+'/synthesizer_offset/files/AgeMW_vs_Mass_tau_bins_'+code+'_'+model+'_'+imf+''+met_string+'_MostSign_Solutions.eps',2)

;psfile=strtrim(dir+'/synthesizer_offset/files/Age_vs_Mass_tau_bins_'+code+'_'+model+'_'+imf+''+met_string+'_PrimarySolutions.eps',2)

set_plot, 'ps'
device, filename=psfile,/portrait,/times, /color

;================Tau=========================

;wp1=wtau1
;wp2=wtau2
;wp3=wtau3


wp1=wtau1_ps  ;46
wp2=wtau2_ps  ;45
wp3=wtau3_ps  ;13


colorin=fltarr(n_elements(x))

colorin(onesol(wtau1_ps))=cgcolor('Blu6')
colorin(onesol(wtau2_ps))=cgcolor('GRN5')
colorin(onesol(wtau3_ps))=cgcolor('RYB1')

colorin(degsol(wtau1_deg))=cgcolor('Blu6')
colorin(degsol(wtau2_deg))=cgcolor('GRN5')
colorin(degsol(wtau3_deg))=cgcolor('RYB1')

;-------------------------------------

!p.multi=0
plot, x,y , psym=3, xrange=[xmin, xmax], yrange=[ymin, ymax], xtitle='log (M/'+Msun+')', ytitle=' !S!A-!R!Nt!dM!N [Gyr]',$
    charsize=1.2, thick=2., charthick=4., xthick=4., ythick=4.,/xstyle, /ystyle, /nodata, position=[0.12,0.12,0.95,0.95]


cgcolorfill, [qmass_ps(1), qmass_ps(1), qmass_ps(3), qmass_ps(3), qmass_ps(1)],[qage_ps(1), qage_ps(3), qage_ps(3), qage_ps(1), qage_ps(1)],color=cgcolor('RED2')

for i=0, n_elements(onesol)-1 do begin  
oplot,  [x(onesol(i)),x(onesol(i))] ,[y(onesol(i)),y(onesol(i))], psym=sym(1), symsize=1.3,col=colorin(onesol(i))
oplot,  [x(onesol(i)),x(onesol(i))] ,[y(onesol(i)),y(onesol(i))], psym=sym(6), symsize=1.3
endfor


for i=0, n_elements(degsol)-1 do begin  
oplot,  [x(degsol(i)),x(degsol(i))] ,[y(degsol(i)),y(degsol(i))], psym=sym(2), symsize=1.3,col=colorin(degsol(i))
oplot,  [x(degsol(i)),x(degsol(i))] ,[y(degsol(i)),y(degsol(i))], psym=sym(7), symsize=1.3
endfor

xyouts,0.18, 0.87,/norm,'!4s!3 [Myr] < 80     ',  charthick=3.,  charsize=1., col=cgcolor('blu6')
xyouts,0.18, 0.82,/norm,'80 Myr < !4s!3 [Myr] < 400 ',  charthick=3.,  charsize=1., col=cgcolor('grn5')
xyouts,0.18, 0.77,/norm,'!4s!3 [Myr] > 400      ',  charthick=3.,  charsize=1., col=cgcolor('Ryb1')


;-----median values and errors------
;=====================================

oplot, [median(mass), median(mass)],[median(age_mw_gy),median(age_mw_gy)], psym=sym(4), col=8, symsize=4.0
oplot, [median(mass), median(mass)],[median(age_mw_gy),median(age_mw_gy)], psym=sym(9), symsize=4.0


;-----

oplot, [mass_inf1, mass_sup1],[med_age1,med_age1],  thick=5., col=cgcolor('blu6')          ;68% in mass
oplot, [mass_inf1, mass_inf1], [med_age1+0.05,med_age1-0.05],  thick=5., col=cgcolor('blu6') ;err hat
oplot, [mass_sup1, mass_sup1], [med_age1+0.05,med_age1-0.05],  thick=5., col=cgcolor('blu6')

oplot, [med_mass1, med_mass1],[age_inf1,age_sup1],  thick=5., col=cgcolor('blu6')          ;68% in age
oplot, [med_mass1+0.01, med_mass1-0.01],[age_inf1,age_inf1],  thick=5., col=cgcolor('blu6') ;hat
oplot, [med_mass1+0.01, med_mass1-0.01],[age_sup1,age_sup1],  thick=5., col=cgcolor('blu6') 

oplot, [med_mass1, med_mass1],[med_age1,med_age1],  psym=sym(4),symsize=2.0, thick=3., col=cgcolor('blu6')
oplot, [med_mass1, med_mass1],[med_age1,med_age1],  psym=sym(9),symsize=2.0, thick=3.


;----

oplot, [mass_inf2, mass_sup2],[med_age2,med_age2],  thick=5., col=cgcolor('grn5')          ;68% in mass
oplot, [mass_inf2, mass_inf2], [med_age2+0.05,med_age2-0.05],  thick=5., col=cgcolor('grn5') ;err hat
oplot, [mass_sup2, mass_sup2], [med_age2+0.05,med_age2-0.05],  thick=5., col=cgcolor('grn5')

oplot, [med_mass2, med_mass2],[age_inf2,age_sup2],  thick=5., col=cgcolor('grn5')          ;68% in age
oplot, [med_mass2+0.01, med_mass2-0.01],[age_inf2,age_inf2],  thick=5., col=cgcolor('grn5') ;hat
oplot, [med_mass2+0.01, med_mass2-0.01],[age_sup2,age_sup2],  thick=5., col=cgcolor('grn5') 

oplot, [med_mass2, med_mass2],[med_age2,med_age2],  psym=sym(4),symsize=2.0, thick=3., col=cgcolor('grn5')
oplot, [med_mass2, med_mass2],[med_age2,med_age2],  psym=sym(9),symsize=2.0, thick=3.

;---------------------


oplot, [mass_inf3, mass_sup3],[med_age3,med_age3],  thick=5., col=cgcolor('Ryb1')          ;68% in mass
oplot, [mass_inf3, mass_inf3], [med_age3+0.05,med_age3-0.05],  thick=5., col=cgcolor('Ryb1') ;err hat
oplot, [mass_sup3, mass_sup3], [med_age3+0.05,med_age3-0.05],  thick=5., col=cgcolor('Ryb1')

oplot, [med_mass3, med_mass3],[age_inf3,age_sup3],  thick=5., col=cgcolor('Ryb1')          ;68% in age
oplot, [med_mass3+0.01, med_mass3-0.01],[age_inf3,age_inf3],  thick=5., col=cgcolor('Ryb1') ;hat
oplot, [med_mass3+0.01, med_mass3-0.01],[age_sup3,age_sup3],  thick=5., col=cgcolor('Ryb1') 

oplot, [med_mass3, med_mass3],[med_age3,med_age3],  psym=sym(4),symsize=2.0, thick=3., col=cgcolor('Ryb1')
oplot, [med_mass3, med_mass3],[med_age3,med_age3],  psym=sym(9),symsize=2.0, thick=3.



;oplot, [11.1, 11.1],[0.7,0.7]+0.02, psym=sym(4), col=8, symsize=1.5
oplot, [11.1, 11.1],[0.7,0.7]+0.02, psym=sym(9), symsize=1.5

xyouts, 11.15, 0.7, 'Median values',  charthick=3.,  charsize=1.2

oplot, [11.1, 11.1],[0.5,0.5]+0.02, psym=sym(1)
xyouts, 11.14, 0.5, '1 solution',  charthick=3.,  charsize=1.2

oplot, [11.1, 11.1],[0.3,0.3]+0.02, psym=sym(2)
xyouts, 11.14, 0.3, '> 1 solutions',  charthick=3.,  charsize=1.2


;errors
;-------

;oploterror, [11.3, 11.3],[1.2,1.2] ,err_mass, err_age*1.2, thick=4, errcol=cgcolor("Black"), psym=3, symsize=0.5;, /nohat
oploterror, [11.3, 11.3],[1.2,1.2] ,err_mass, err_age, thick=4, errcol=cgcolor("Black"), psym=3, symsize=0.5;, /nohat


device,/close 
set_plot,'x'

print, 'FIN PLOT'
STOP
STOP



x=mass(ps)
y=age_mw_gy(ps)

xmin=9.7
xmax=11.5
ymin=0
ymax=5.



print, 'empieza PLOT CLASS'

stop
;=============================================


psfile=strtrim(dir+'/synthesizer/plots_mine/AgeMW_vs_Mass_class_'+code+'_'+model+'_'+imf+''+met_string+'_PrimarySolutions.eps',2)
set_plot, 'ps'
device, filename=psfile,/portrait,/times, /color

;================Tau=========================

wp1=bc
wp2=gv
wp3=psb
wp4=oq


colorin1=4
colorin2=3
colorin3=8
colorin4=2

; cte size
ss(0:116)=1.

w=ps

;-------------------------------------

!p.multi=0
plot, x(w),y(w) , psym=3, xrange=[xmin, xmax], yrange=[ymin, ymax], xtitle='log Mass ('+Msun+')', ytitle='Mass-Weighted Age (Gyr)',$
    charsize=1.5, thick=2., charthick=4., xthick=4., ythick=4.,/xstyle, /ystyle, /nodata


cgcolorfill, [qmass_ps(1), qmass_ps(1), qmass_ps(3), qmass_ps(3), qmass_ps(1)],[qage_ps(1), qage_ps(3), qage_ps(3), qage_ps(1), qage_ps(1)],color=cgcolor('Light grey')


oplot,  x(wp1),y(wp1), psym=sym(1), symsize=1.,col=cgcolor('Blu8')
oplot,  x(wp2),y(wp2), psym=sym(1), symsize=1.,col=cgcolor('GRN6')
oplot,  x(wp3),y(wp3), psym=sym(1), symsize=1.,col=cgcolor('RYB1')

oplot,  x(wp1),y(wp1), psym=sym(6), symsize=1.
oplot,  x(wp2),y(wp2), psym=sym(6), symsize=1.
oplot,  x(wp3),y(wp3), psym=sym(6), symsize=1.

;different sizes
;---------------


for i=0, n_elements(wp1)-1 do begin &$
oplot,  [x(wp1(i)),x(wp1(i))], [y(wp1(i)),y(wp1(i))], psym=sym(1), symsize=ss(wp1(i)),col=colorin1 &$
oplot,  [x(wp1(i)),x(wp1(i))], [y(wp1(i)),y(wp1(i))], psym=sym(6), symsize=ss(wp1(i)),col=0, thick=2. &$
endfor


for i=0, n_elements(wp2)-1 do begin &$
oplot,  [x(wp2(i)),x(wp2(i))], [y(wp2(i)),y(wp2(i))], psym=sym(1), symsize=ss(wp2(i)),col=colorin2 &$
oplot,  [x(wp2(i)),x(wp2(i))], [y(wp2(i)),y(wp2(i))], psym=sym(6), symsize=ss(wp2(i)),col=0, thick=2. &$
endfor


for i=0, n_elements(wp3)-1 do begin &$
oplot,  [x(wp3(i)),x(wp3(i))], [y(wp3(i)),y(wp3(i))], psym=sym(1), symsize=ss(wp3(i)),col=colorin3 &$
oplot,  [x(wp3(i)),x(wp3(i))], [y(wp3(i)),y(wp3(i))], psym=sym(6), symsize=ss(wp3(i)),col=0, thick=2. &$
endfor




for i=0, n_elements(wp4)-1 do begin &$
oplot,  [x(wp4(i)),x(wp4(i))], [y(wp4(i)),y(wp4(i))], psym=sym(1), symsize=ss(wp4(i)),col=colorin4 &$
oplot,  [x(wp4(i)),x(wp4(i))], [y(wp4(i)),y(wp4(i))], psym=sym(6), symsize=ss(wp4(i)),col=0, thick=2. &$
endfor



xyouts,10., 4.4,'Blue Cloud',  charthick=2.,  charsize=1., col=colorin1
xyouts,10., 4.2,'Green Valley',  charthick=2.,  charsize=1., col=colorin2
xyouts,10., 4.,'Post-SB',  charthick=2.,  charsize=1., col=colorin3
xyouts,10., 3.8,'Old Quiescent',  charthick=2.,  charsize=1., col=colorin4



device,/close 
set_plot,'x'

print, 'FIN PLOT'
STOP
STOP

print, 'empieza PLOT AV'

stop
;=============================================


psfile=strtrim(dir+'/synthesizer/plots_mine/AgeMW_vs_Mass_av_bins_'+code+'_'+model+'_'+imf+''+met_string+'.eps',2)
set_plot, 'ps'
device, filename=psfile,/portrait,/times, /color

;================AV=========================

wp1=wav1
wp2=wav2
wp3=wav3

colorin=fltarr(n_elements(tau_my))
colorin(wp1)=cgcolor('Steel Blue')
colorin(wp2)=cgcolor('Forest Green')
colorin(wp3)=cgcolor('Red')



;-------------------------------------

!p.multi=0
plot, x(w),y(w) , psym=3, xrange=[xmin, xmax], yrange=[ymin, ymax], xtitle='log Mass ('+Msun+')', ytitle='Mass-Weighted Age (Gyr)',$
    charsize=1., thick=2., charthick=2., xthick=2., ythick=2.,/xstyle, /ystyle

for i=0, n_elements(wp1)-1 do begin &$
oplot,  [x(wp1(i)),x(wp1(i))], [y(wp1(i)),y(wp1(i))], psym=sym(1), symsize=ss(wp1(i)),col=colorin(wp1(i)) &$
oplot,  [x(wp1(i)),x(wp1(i))], [y(wp1(i)),y(wp1(i))], psym=sym(6), symsize=ss(wp1(i)),col=0, thick=2. &$
endfor


for i=0, n_elements(wp2)-1 do begin &$
oplot,  [x(wp2(i)),x(wp2(i))], [y(wp2(i)),y(wp2(i))], psym=sym(1), symsize=ss(wp2(i)),col=colorin(wp2(i)) &$
oplot,  [x(wp2(i)),x(wp2(i))], [y(wp2(i)),y(wp2(i))], psym=sym(6), symsize=ss(wp2(i)),col=0, thick=2. &$
endfor


for i=0, n_elements(wp3)-1 do begin &$
oplot,  [x(wp3(i)),x(wp3(i))], [y(wp3(i)),y(wp3(i))], psym=sym(1), symsize=ss(wp3(i)),col=colorin(wp3(i)) &$
oplot,  [x(wp3(i)),x(wp3(i))], [y(wp3(i)),y(wp3(i))], psym=sym(6), symsize=ss(wp3(i)),col=0, thick=2. &$
endfor




oplot, x(ww), y(ww), psym=1, thick=4. , col=200

;xyouts, ,'OK solutions ('+strcompress(string(n_elements(w)), /remove_all)+')'

xyouts,0.20, 0.9,/norm,  'Good solutions ('+strcompress(string(n_elements(w)), /remove_all)+')',  charthick=2.,  charsize=.8

xyouts,0.20, 0.85,/norm,'       AV < 0.5',  charthick=2.,  charsize=.8, col=cgcolor('Royal Blue')
xyouts,0.20, 0.80,/norm,' 0.5 < AV < 1.0',  charthick=2.,  charsize=.8, col=cgcolor('Forest Green')
xyouts,0.20, 0.75,/norm,'       AV > 1.0',  charthick=2.,  charsize=.8, col=cgcolor('Red')



device,/close 
set_plot,'x'

print, 'FIN PLOT'
STOP
STOP



print, 'empieza PLOT Met'

stop
;=============================================


psfile=strtrim(dir+'/synthesizer/plots_mine/AgeMW_vs_Mass_met_bins_'+code+'_'+model+'_'+imf+''+met_string+'.eps',2)
set_plot, 'ps'
device, filename=psfile,/portrait,/times, /color

;================Met========================

wp1=wm1
wp2=wm2
wp3=wm3

colorin=fltarr(n_elements(tau_my))
colorin(wp1)=cgcolor('Steel Blue')
colorin(wp2)=cgcolor('Forest Green')
colorin(wp3)=cgcolor('Red')



;-------------------------------------

!p.multi=0
plot, x(w),y(w) , psym=3, xrange=[xmin, xmax], yrange=[ymin, ymax], xtitle='log Mass ('+Msun+')', ytitle='Mass-Weighted Age (Gyr)',$
    charsize=1., thick=2., charthick=2., xthick=2., ythick=2.,/xstyle, /ystyle

for i=0, n_elements(wp1)-1 do begin &$
oplot,  [x(wp1(i)),x(wp1(i))], [y(wp1(i)),y(wp1(i))], psym=sym(1), symsize=ss(wp1(i)),col=colorin(wp1(i)) &$
oplot,  [x(wp1(i)),x(wp1(i))], [y(wp1(i)),y(wp1(i))], psym=sym(6), symsize=ss(wp1(i)),col=0, thick=2. &$
endfor


for i=0, n_elements(wp2)-1 do begin &$
oplot,  [x(wp2(i)),x(wp2(i))], [y(wp2(i)),y(wp2(i))], psym=sym(1), symsize=ss(wp2(i)),col=colorin(wp2(i)) &$
oplot,  [x(wp2(i)),x(wp2(i))], [y(wp2(i)),y(wp2(i))], psym=sym(6), symsize=ss(wp2(i)),col=0, thick=2. &$
endfor


for i=0, n_elements(wp3)-1 do begin &$
oplot,  [x(wp3(i)),x(wp3(i))], [y(wp3(i)),y(wp3(i))], psym=sym(1), symsize=ss(wp3(i)),col=colorin(wp3(i)) &$
oplot,  [x(wp3(i)),x(wp3(i))], [y(wp3(i)),y(wp3(i))], psym=sym(6), symsize=ss(wp3(i)),col=0, thick=2. &$
endfor




oplot, x(ww), y(ww), psym=1, thick=4. , col=200

;xyouts, ,'OK solutions ('+strcompress(string(n_elements(w)), /remove_all)+')'

xyouts,0.20, 0.9,/norm,  'Good solutions ('+strcompress(string(n_elements(w)), /remove_all)+')',  charthick=2.,  charsize=.8

xyouts,0.20, 0.85,/norm,'Z = 0.4x'+Zsun+'  ',  charthick=2.,  charsize=.8, col=cgcolor('Royal Blue')
xyouts,0.20, 0.80,/norm,'Z = '+Zsun+'',  charthick=2.,  charsize=.8, col=cgcolor('Forest Green')
xyouts,0.20, 0.75,/norm,'Z = 2.5x'+Zsun+'',  charthick=2.,  charsize=.8, col=cgcolor('Red')



device,/close 
set_plot,'x'

print, 'FIN PLOT'
STOP
STOP


END
