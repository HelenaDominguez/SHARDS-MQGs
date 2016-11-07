PRO plots_age_vs_tau_ssfr

;  .comp /home/hdominguez/idl/pro/new/plots_age_vs_tau_ssfr.pro

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



;selected galaxies
;------------------
readcol, dir+'/out_HDFN_SHARDS_iDR2gammaP01.20150305_z_1.00_1.50_104',num_rd, id_rd_104,RA_rd,DEC_rd,z_rd,q_rd,VJ_rd,UV_rd,mass_rd,sSFR_rd,SFR_UV_obs_rd,SFR_IR_rd, SFR_tot_rd , AV_rd,mag_H, mag_j,a,a,a,sel_string, Format='F, A,F,F,F,F,F,F,f,f,f,f,f,f,f,f,a,a,a,a', skipline=1
; same order than previous files



readcol, dir+'/synthesizer_offset/files/analyze_parameters_BC_txp_krou_flags_104.lst',   shards_id , nclus , z ,  log_tau,  log_age, log_age_MW,  av , met, mass,ratio, flag_d4, flag_mg, flag_deg, flag_sol, format='a, f, f, f'

age_gy=10^log_age*1.e-9
age_mw_gy=10^log_age_mw*1.e-9
tau_my=10^log_tau*1.e-6



restore, dir+'/synthesizer_offset/files/index_BC03_txp_krou_104_d4_err.sav', /verbose


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


;rename
sfr_now=sfr_bf
sfr_bf=sfr_l100
ssfr_now=ssfr_bf
ssfr_bf=ssfr_l100

stop


pc_mass=mass_lgyr*100./10^mass_bf  ;% mass formed in the last Gyr
pc_mass(where(age_gy lt 1.0))=100.  ; galaxies younger than 1 Gyr have formed all their mass in the last Gyr

pc1=where(pc_mass le 10.)    ; 48
pc2=where(pc_mass gt 10. and pc_mass le 50.) ;18
pc3=where(pc_mass gt 50.)  ;38



size_vect=fltarr(n_elements(source))

size_vect(pc1)=1.0
size_vect(pc2)=2.0
size_vect(pc3)=3.0

symbol_vect=intarr(n_elements(source))
symbol_vect(where(sel_string eq 'UVJ'))=1
symbol_vect(where(sel_string eq 'sSFR'))=4


symbol_vect2=intarr(n_elements(source))
symbol_vect2(where(sel_string eq 'UVJ'))=6
symbol_vect2(where(sel_string eq 'sSFR'))=9



stop



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

ps=where(flag_sol ge 1)
ok=where(mass(ps) ge 10.)


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



;OLD vs YOUNG
;==================


young=where(age_mw_gy lt 2.0)
;int=where(age_mw_gy ge 1.0 and age_mw_gy lt 2.0)
old=where(age_mw_gy ge 2.0)


qmass1=quartile(mass(young))
qmass2=quartile(mass(old))

qtau1=quartile(tau_my(young))
qtau2=quartile(tau_my(old))


qav1=quartile(av(young))
qav2=quartile(av(old))


qmet1=quartile(met(young))
qmet2=quartile(met(old))



;========================
ssfr_lim=where(ssfr_bf lt 1.e-4)
ssfr_lim2=where(ssfr_bf gt 0.1)

ssfr_bf_original=ssfr_bf
ssfr_bf(ssfr_lim)=1.e-4
;ssfr_bf(ssfr_lim2)=0.3

stop


ssfr_red=-4.
;ssfr_blue=-0.2
ssfr_blue=-0.1


;ssfr_red=1.e-4
;ssfr_blue=0.3

a=255./(ssfr_red-ssfr_blue)
b=(-ssfr_blue)*a

z=alog10(ssfr_bf)
;z=ssfr_bf
colorin=a*z+b


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

psfile=strtrim(dir+'/synthesizer_offset/files/AgeMW_vs_Tau_sSFR_colors_'+code+'_'+model+'_'+imf+''+met_string+'_PrimarySolutions_size_mass_l100.eps',2)

;psfile=strtrim(dir+'/synthesizer_offset/files/Age_vs_Tau_mass_bins_'+code+'_'+model+'_'+imf+''+met_string+'_PrimarySolutions.eps',2)

set_plot, 'ps'
device, filename=psfile,/portrait,/times, /color

;================class=========================


colours
; cte size
;ss(0:116)=1.

;-------------------------------------


!p.multi=0
plot, x,y , psym=3, xrange=[xmin, xmax], yrange=[ymin, ymax], xtitle='!4s!3  [Myr] ', ytitle='!S!A-!R!Nt!dM!N [Gyr]',$
    charsize=1.2, thick=2., charthick=4., xthick=4., ythick=4.,/xstyle, /ystyle, /xlog, /nodata, position=[0.12,0.12,0.95,0.95]
;oplot, [1.5,1.5],[1.9, 2.1], col=4, thick=3.

cgcolorfill, [qtau(1), qtau(1), qtau(3), qtau(3), qtau(1)],[qage(1), qage(3), qage(3), qage(1), qage(1)],color=cgcolor('Light grey')


   cgLoadCT, 25, NColors=10, Bottom=1
   ;cgColorbar, NColors=10, Bottom=1, /Discrete, Range=[-2.,-4.], charsize=0.8, charthick=2.5, position=[0.2,0.85,0.9,0.88]
   cgColorbar, NColors=10, Bottom=1, /Discrete, Range=[-1.,-4.], charsize=1.0, charthick=2.5, position=[0.2,0.85,0.9,0.88]


cgLoadCT, 25

for i=0, n_elements(x)-1 do begin  &$
oplot,  [x(i),x(i)] ,[y(i),y(i)], symsize=size_vect(i),col=colorin(i), psym=sym(symbol_vect(i))  &$
oplot,  [x(i),x(i)] ,[y(i),y(i)], symsize=size_vect(i), col=cgcolor('black'), psym=sym(symbol_vect2(i))  &$
endfor




;errors
;-------

;oploterror, [500, 500],[0.7,0.7] ,err_tau*500, err_age*0.7, thick=4, errcol=cgcolor("Black"), psym=3, symsize=0.5;, /nohat
;oploterror, [500, 500],[0.7,0.7] ,err_tau, err_age, thick=4, errcol=cgcolor("Black"), psym=3, symsize=0.5;, /nohat
oploterror, [500, 500],[0.7,0.7] ,0.4*err_tau*500., err_age, thick=4, errcol=cgcolor("Black"), psym=3, symsize=0.5;, /nohat


;sizes legend
;-------------

oplot, [3.,3.],[3.5,3.5], psym=sym(6), symsize=1., col=cgcolor('black'), thick=1.3
xyouts, 4., 3.5-0.05, '< 10% mass in last Gyr', charsize=1.2, charthick=2.,col=cgcolor('black')

oplot, [3.,3.],[3.2,3.2], psym=sym(6), symsize=2., col=cgcolor('black'), thick=1.3
xyouts, 4., 3.2-0.05, '10-50% mass in last Gyr', charsize=1.2, charthick=2.,col=cgcolor('black')


oplot, [3.,3.],[2.9,2.9], psym=sym(6), symsize=3.5, col=cgcolor('black'), thick=1.3
xyouts, 4., 2.9-0.05, '> 50% mass in last Gyr', charsize=1.2, charthick=2.,col=cgcolor('black')


oplot, [100.,100.],[3.5,3.5], psym=sym(6), symsize=1.2, col=cgcolor('black'), thick=1.3
xyouts, 120., 3.5-0.05, 'UVJ selected', charsize=1.2, charthick=2.,col=cgcolor('black')


oplot, [100.,100.],[3.2,3.2], psym=sym(9), symsize=1.2, col=cgcolor('black'), thick=1.3
xyouts, 120., 3.2-0.05, 'sSFR selected', charsize=1.2, charthick=2.,col=cgcolor('black')



xyouts, 20, 4.23,'log sSFR!dSED!N [Gyr!u-1!N]', charsize=1.2, charthick=3.,col=cgcolor('black')

device,/close 
set_plot,'x'

print, 'FIN PLOT'
STOP
STOP



END

