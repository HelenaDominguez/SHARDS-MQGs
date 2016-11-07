PRO plots_age_vs_tau_bins

;    .comp /home/hdominguez/idl/pro/new/plots_age_vs_tau_bins.pro

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

;Primary Solutions from index
readcol, dir+'/synthesizer_offset/files/analyze_parameters_BC_txp_krou_flags_104.lst',   shards_id , nclus , z ,  log_tau,  log_age, log_age_MW,  av , met, mass,ratio, flag_d4, flag_mg, flag_deg, flag_sol, format='a, f, f, f'


;Most significant solutions
;readcol, dir+'/synthesizer_offset/files/analyze_parameters_BC_txp_krou_most_significant_sol_104.lst',   shards_id , nclus , z ,  log_tau,  log_age, log_age_MW,  av , met, mass,ratio, format='a, f, f, f'


;SSP solutions
;readcol, dir+'/synthesizer_offset/files/analyze_parameters_SSP_txp_krou_101_MostSig.lst',   shards_id , nclus , z ,  log_tau,  log_age, log_age_MW,  av , met, mass,ratio, format='a, f, f, f'



age_gy=10^log_age*1.e-9
age_mw_gy=10^log_age_mw*1.e-9
tau_my=10^log_tau*1.e-6



restore, dir+'/synthesizer_offset/files/index_BC03_txp_krou_104_d4_err.sav', /verbose



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


onesol=where(ratio eq 100.)
degsol=where(ratio lt 100.)


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


young=where(age_mw_gy(ok) lt 2.0)
;int=where(age_mw_gy(ps(ok)) ge 1.0 and age_mw_gy(ps(ok)) lt 2.0)
old=where(age_mw_gy(ok) ge 2.0)


qmass1=quartile(mass(young))
qmass2=quartile(mass(old))

qtau1=quartile(tau_my(young))
qtau2=quartile(tau_my(old))


qav1=quartile(av(young))
qav2=quartile(av(old))


qmet1=quartile(met(young))
qmet2=quartile(met(old))


wm1=where(mass lt 10.5)
wm2=where(mass ge 10.5  and  mass lt 10.8)
wm3=where(mass ge 10.8)


wm1_ps=where(mass(onesol) lt 10.5)
wm2_ps=where(mass(onesol) ge 10.5  and  mass(onesol) lt 10.8)
wm3_ps=where(mass(onesol) ge 10.8)


wm1_deg=where(mass(degsol) lt 10.5 )
wm2_deg=where(mass(degsol) ge 10.5  and  mass(degsol) lt 10.8)
wm3_deg=where(mass(degsol) ge 10.8)



;-----median and scatter in each bin-----

med_tau1=median(tau_my(wm1))
med_age1=median(age_mw_gy(wm1))

med_tau2=median(tau_my(wm2))
med_age2=median(age_mw_gy(wm2))

med_tau3=median(tau_my(wm3))
med_age3=median(age_mw_gy(wm3))

stv_tau1=stdev(tau_my(wm1))
stv_age1=stdev(age_mw_gy(wm1))

stv_tau2=stdev(tau_my(wm2))
stv_age2=stdev(age_mw_gy(wm2))

stv_tau3=stdev(tau_my(wm3))
stv_age3=stdev(age_mw_gy(wm3))

stop

;------log scale------

;=====tau========
num1=fix(n_elements(wm1)*0.136)  ;1sigma lower lim
num2=fix(n_elements(wm1)*0.842)  ;1 sigma upper lim

ord1=sort(log_tau(wm1))
tau_ord1=log_tau(wm1(ord1))

tau_inf1=10^tau_ord1(num1)*1.e-6
tau_sup1=10^tau_ord1(num2)*1.e-6
delta_tau1=tau_sup1-tau_inf1

;------

num1=fix(n_elements(wm2)*0.136)  ;1sigma lower lim
num2=fix(n_elements(wm2)*0.842)  ;1 sigma upper lim

ord2=sort(log_tau(wm2))
tau_ord2=log_tau(wm2(ord2))

tau_inf2=10^tau_ord2(num1)*1.e-6
tau_sup2=10^tau_ord2(num2)*1.e-6
delta_tau2=tau_sup2-tau_inf2

;------

num1=fix(n_elements(wm3)*0.136)  ;1sigma lower lim
num2=fix(n_elements(wm3)*0.842)  ;1 sigma upper lim

ord3=sort(log_tau(wm3))
tau_ord3=log_tau(wm3(ord3))

tau_inf3=10^tau_ord3(num1)*1.e-6
tau_sup3=10^tau_ord3(num2)*1.e-6
delta_tau3=tau_sup3-tau_inf3

;------


;=====age========
num1=fix(n_elements(wm1)*0.136)  ;1sigma lower lim
num2=fix(n_elements(wm1)*0.842)  ;1 sigma upper lim

ord1=sort(log_age_mw(wm1))
age_ord1=log_age_mw(wm1(ord1))

age_inf1=10^age_ord1(num1)*1.e-9
age_sup1=10^age_ord1(num2)*1.e-9

;---
num1=fix(n_elements(wm2)*0.136)  ;1sigma lower lim
num2=fix(n_elements(wm2)*0.842)  ;1 sigma upper lim

ord2=sort(log_age_mw(wm2))
age_ord2=log_age_mw(wm2(ord2))

age_inf2=10^age_ord2(num1)*1.e-9
age_sup2=10^age_ord2(num2)*1.e-9

;---
num1=fix(n_elements(wm3)*0.136)  ;1sigma lower lim
num2=fix(n_elements(wm3)*0.842)  ;1 sigma upper lim

ord3=sort(log_age_mw(wm3))
age_ord3=log_age_mw(wm3(ord3))

age_inf3=10^age_ord3(num1)*1.e-9
age_sup3=10^age_ord3(num2)*1.e-9





stop






x=tau_my
y=age_mw_gy
;y=age_gy


xmin=2.
xmax=1000.
ymin=0
ymax=4.5
;ymax=5.5


print, 'Empieza PLOT'
stop

;psfile=strtrim(dir+'/synthesizer_offset/files/AgeMW_vs_Tau_mass_bins_'+code+'_'+model+'_'+imf+''+met_string+'_SSP.eps',2)

psfile=strtrim(dir+'/synthesizer_offset/files/AgeMW_vs_Tau_mass_bins_'+code+'_'+model+'_'+imf+''+met_string+'_PrimarySolutions.eps',2)


;psfile=strtrim(dir+'/synthesizer_offset/files/AgeMW_vs_Tau_mass_bins_'+code+'_'+model+'_'+imf+''+met_string+'_MostSign_Solutions.eps',2)


;psfile=strtrim(dir+'/synthesizer_offset/files/Age_vs_Tau_mass_bins_'+code+'_'+model+'_'+imf+''+met_string+'_PrimarySolutions.eps',2)

set_plot, 'ps'
device, filename=psfile,/portrait,/times, /color

;================class=========================


colorin=fltarr(n_elements(x))

colorin(onesol(wm1_ps))=cgcolor('Blu6')
colorin(onesol(wm2_ps))=cgcolor('GRN5')
colorin(onesol(wm3_ps))=cgcolor('RYB1')


;colorin(onesol(wm1_ps))=4
;colorin(onesol(wm2_ps))=3
;colorin(onesol(wm3_ps))=2


;colorin(degsol(wm1_deg))=4
;colorin(degsol(wm2_deg))=3
;colorin(degsol(wm3_deg))=2


colorin(degsol(wm1_deg))=cgcolor('Blu6')
colorin(degsol(wm2_deg))=cgcolor('GRN5')
colorin(degsol(wm3_deg))=cgcolor('RYB1')



; cte size
;ss(0:116)=1.

;-------------------------------------

!p.multi=0
plot, x,y , psym=3, xrange=[xmin, xmax], yrange=[ymin, ymax], xtitle='!4s!3  [Myr] ', ytitle=' !S!A-!R!Nt!dM!N  [Gyr]',$
    charsize=1.2, thick=2., charthick=4., xthick=4., ythick=4.,/xstyle, /ystyle, /xlog, /nodata, position=[0.12,0.12,0.95,0.95]

cgcolorfill, [qtau_ps(1), qtau_ps(1), qtau_ps(3), qtau_ps(3), qtau_ps(1)],[qage_ps(1), qage_ps(3), qage_ps(3), qage_ps(1), qage_ps(1)],color=cgcolor('Red2')


for i=0, n_elements(onesol)-1 do begin  &$
oplot,  [x(onesol(i)),x(onesol(i))] ,[y(onesol(i)),y(onesol(i))], psym=sym(1), symsize=1.3,col=colorin(onesol(i))  &$
oplot,  [x(onesol(i)),x(onesol(i))] ,[y(onesol(i)),y(onesol(i))], psym=sym(6), symsize=1.3  &$
endfor


for i=0, n_elements(degsol)-1 do begin   &$
oplot,  [x(degsol(i)),x(degsol(i))] ,[y(degsol(i)),y(degsol(i))], psym=sym(2), symsize=1.3,col=colorin(degsol(i)) &$
oplot,  [x(degsol(i)),x(degsol(i))] ,[y(degsol(i)),y(degsol(i))], psym=sym(7), symsize=1.3 &$
endfor




;-----median values------
;==============================

;-----------------
oplot, [tau_inf1, tau_sup1],[med_age1,med_age1],  thick=5., col=cgcolor('blu6')          ;68% in mass
oplot, [tau_inf1, tau_inf1], [med_age1+0.05,med_age1-0.05],  thick=5., col=cgcolor('blu6') ;err hat
oplot, [tau_sup1, tau_sup1], [med_age1+0.05,med_age1-0.05],  thick=5., col=cgcolor('blu6')

oplot, [med_tau1, med_tau1],[age_inf1,age_sup1],  thick=5., col=cgcolor('blu6')          ;68% in age
oplot, [med_tau1+4, med_tau1-2],[age_inf1,age_inf1],  thick=5., col=cgcolor('blu6') ;hat
oplot, [med_tau1+4, med_tau1-2],[age_sup1,age_sup1],  thick=5., col=cgcolor('blu6') 

oplot, [med_tau1, med_tau1],[med_age1,med_age1],  psym=sym(4),symsize=2.0, thick=3., col=cgcolor('blu6')
oplot, [med_tau1, med_tau1],[med_age1,med_age1],  psym=sym(9),symsize=2.0, thick=3.


;-----------------
oplot, [tau_inf2, tau_sup2],[med_age2,med_age2],  thick=5., col=cgcolor('grn5')          ;68% in mass
oplot, [tau_inf2, tau_inf2], [med_age2+0.05,med_age2-0.05],  thick=5., col=cgcolor('grn5') ;err hat
oplot, [tau_sup2, tau_sup2], [med_age2+0.05,med_age2-0.05],  thick=5., col=cgcolor('grn5')

oplot, [med_tau2, med_tau2],[age_inf2,age_sup2],  thick=5., col=cgcolor('grn5')          ;68% in age
oplot, [med_tau2+10, med_tau2-7],[age_inf2,age_inf2],  thick=5., col=cgcolor('grn5') ;hat
oplot, [med_tau2+10, med_tau2-7],[age_sup2,age_sup2],  thick=5., col=cgcolor('grn5') 

oplot, [med_tau2, med_tau2],[med_age2,med_age2],  psym=sym(4),symsize=2.0, thick=3., col=cgcolor('grn5')
oplot, [med_tau2, med_tau2],[med_age2,med_age2],  psym=sym(9),symsize=2.0, thick=3.

;-----------------
oplot, [tau_inf3, tau_sup3],[med_age3,med_age3],  thick=5., col=cgcolor('Ryb1')          ;68% in mass
oplot, [tau_inf3, tau_inf3], [med_age3+0.05,med_age3-0.05],  thick=5., col=cgcolor('Ryb1') ;err hat
oplot, [tau_sup3, tau_sup3], [med_age3+0.05,med_age3-0.05],  thick=5., col=cgcolor('Ryb1')

oplot, [med_tau3, med_tau3],[age_inf3,age_sup3],  thick=5., col=cgcolor('Ryb1')          ;68% in age
oplot, [med_tau3+15, med_tau3-15],[age_inf3,age_inf3],  thick=5., col=cgcolor('Ryb1') ;hat
oplot, [med_tau3+15, med_tau3-15],[age_sup3,age_sup3],  thick=5., col=cgcolor('Ryb1') 

oplot, [med_tau3, med_tau3],[med_age3,med_age3],  psym=sym(4),symsize=2.0, thick=3., col=cgcolor('Ryb1')
oplot, [med_tau3, med_tau3],[med_age3,med_age3],  psym=sym(9),symsize=2.0, thick=3.


oplot, [median(tau_my), median(tau_my)],[median(age_mw_gy),median(age_mw_gy)], psym=sym(4), col=8, symsize=4.0
oplot, [median(tau_my), median(tau_my)],[median(age_mw_gy),median(age_mw_gy)], psym=sym(9), symsize=4.0

xyouts,0.18, 0.87,/norm,'log (M/'+Msun+')  < 10.5',  charthick=3.,  charsize=1., col=cgcolor('Blu6')
xyouts,0.18, 0.82,/norm,'10.5 < log (M/'+Msun+') < 10.8',  charthick=3.,  charsize=1., col=cgcolor('GRN5')
xyouts,0.18, 0.77,/norm,'log (M/'+Msun+') > 10.8',  charthick=3.,  charsize=1., col=cgcolor('RYB1')




;errors
;-------

oploterror, [500, 500],[0.7,0.7] ,0.4*err_tau*500., err_age, thick=4, errcol=cgcolor("Black"), psym=3, symsize=0.5;, /nohat
;oploterror, [500, 500],[0.7,0.7] ,err_tau, err_age, thick=4, errcol=cgcolor("Black"), psym=3, symsize=0.5;, /nohat

device,/close 
set_plot,'x'

print, 'FIN PLOT'
STOP
STOP



oploterror, [med_tau1, med_tau1],[med_age1,med_age1],   stv_age1, psym=sym(4), col=cgcolor('blu6'), errcol=cgcolor('blu6'),symsize=2.0, errthick=5.
oplot, [med_tau1, med_tau1],[med_age1,med_age1],  psym=sym(4),symsize=2.0, thick=3., col=cgcolor('blu6')
oplot, [med_tau1, med_tau1],[med_age1,med_age1],  psym=sym(9),symsize=2.0, thick=3.
oplot, [tau_inf1, tau_sup1],[med_age1,med_age1],  thick=5., col=cgcolor('blu6')          ;68% in tau
oplot, [tau_inf1, tau_inf1], [med_age1+0.05,med_age1-0.05],  thick=5., col=cgcolor('blu6') ;err hat
oplot, [tau_sup1, tau_sup1], [med_age1+0.05,med_age1-0.05],  thick=5., col=cgcolor('blu6')



oploterror, [med_tau2, med_tau2],[med_age2,med_age2],  stv_age2, psym=sym(4), col=cgcolor('grn5'), errcol=cgcolor('grn5'),symsize=2.0, errthick=5.
oplot, [med_tau2, med_tau2],[med_age2,med_age2],  psym=sym(4),symsize=2.0, thick=3., col=cgcolor('grn5')
oplot, [med_tau2, med_tau2],[med_age2,med_age2],  psym=sym(9),symsize=2.0, thick=3.
oplot, [tau_inf2, tau_sup2],[med_age2,med_age2],  thick=5., col=cgcolor('grn5')          ;68% in tau
oplot, [tau_inf2, tau_inf2], [med_age2+0.05,med_age2-0.05],  thick=5., col=cgcolor('grn5') ;err hat
oplot, [tau_sup2, tau_sup2], [med_age2+0.05,med_age2-0.05],  thick=5., col=cgcolor('grn5')


oploterror, [med_tau3, med_tau3],[med_age3,med_age3],   stv_age3, psym=sym(4), col=cgcolor('Ryb1'), errcol=cgcolor('Ryb1'), symsize=2.0, errthick=5.
oplot, [med_tau3, med_tau3],[med_age3,med_age3],  psym=sym(4),symsize=2.0, thick=3., col=cgcolor('ryb1')
oplot, [med_tau3, med_tau3],[med_age3,med_age3],  psym=sym(9),symsize=2.0, thick=3.
oplot, [tau_inf3, tau_sup3],[med_age3,med_age3],  thick=5., col=cgcolor('ryb1')          ;68% in tau
oplot, [tau_inf3, tau_inf3], [med_age3+0.05,med_age3-0.05],  thick=5., col=cgcolor('ryb1') ;err hat
oplot, [tau_sup3, tau_sup3], [med_age3+0.05,med_age3-0.05],  thick=5., col=cgcolor('ryb1')


oplot, [median(tau_my), median(tau_my)],[median(age_mw_gy),median(age_mw_gy)], psym=sym(4), col=8, symsize=4.0
oplot, [median(tau_my), median(tau_my)],[median(age_mw_gy),median(age_mw_gy)], psym=sym(9), symsize=4.0

xyouts,0.18, 0.87,/norm,'log (M/'+Msun+')  < 10.5',  charthick=3.,  charsize=1., col=cgcolor('Blu6')
xyouts,0.18, 0.82,/norm,'10.5 < log (M/'+Msun+') < 10.8',  charthick=3.,  charsize=1., col=cgcolor('GRN5')
xyouts,0.18, 0.77,/norm,'log (M/'+Msun+') > 10.8',  charthick=3.,  charsize=1., col=cgcolor('RYB1')




END

