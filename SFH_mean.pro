PRO SFH_mean

; .comp /home/hdominguez/idl/pro/new/SFH_mean.pro

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


readcol, dir+'/synthesizer_offset/files/analyze_parameters_BC_txp_krou_flags_104.lst',   shards_id , nclus , z_list ,  log_tau,  log_age, log_age_MW,  av , met, mass,ratio, flag_d4, flag_mg, flag_deg, flag_sol, format='a,i, f, f, f'

age_gy=10^log_age*1.e-9
age_mw_gy=10^log_age_mw*1.e-9
tau_my=10^log_tau*1.e-6



;readcol, dir+'/out_HDFN_SHARDS_iDR2gammaP01.20150305_z_1.00_1.50_135',num_rd, id_rd,RA_rd,DEC_rd,z_rd,q_rd,VJ_rd,UV_rd,mass_rd,sSFR_rd,SFR_UV_obs_rd,SFR_IR_rd, SFR_tot_rd , AV_rd,mag_H, mag_j, Format='F, A,F,F,F,F,F,F,f,f,f,f,f,f,f,f,a,f', skipline=1
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
now=fltarr(n_elements(id_list))


age_peak=fltarr(n_elements(id_list))
sfr_peak=fltarr(n_elements(id_list))
deltat=fltarr(n_elements(id_list))
age_uz=fltarr(n_elements(id_list))
age_upeak=fltarr(n_elements(id_list))
z_peak=fltarr(n_elements(id_list))
age_uform=fltarr(n_elements(id_list))
z_form=fltarr(n_elements(id_list))


age_peak(0:n_elements(id_list)-1)=-1.
sfr_peak(0:n_elements(id_list)-1)=-1.
deltat(0:n_elements(id_list)-1)=-1.
age_uz(0:n_elements(id_list)-1)=-1.
age_upeak(0:n_elements(id_list)-1)=-1.
z_peak(0:n_elements(id_list)-1)=-1.
age_uform(0:n_elements(id_list)-1)=-1.
z_form(0:n_elements(id_list)-1)=-1.


print, age_peak

for j=0, n_elements(id_list) -1 do begin    &$
;for j=26, 28 do begin 

id=id_list(j)   &$
z=zcomb(j)

print, j,'  ',  id,'  ', z

w=where(shards_id eq id and flag_sol ge 1 )   ;choose good solution

print, w

if w(0) ne -1 then begin

	num=nclus(w)+1
	if code eq 'BC' then  file_sfh=dir+'/synthesizer_offset/res/analyze_rejoin5/clusters/'+id+'.BCall.bc2003.stelib.'+model+'.'+imf+'.sfh.clust'+strtrim(num,1)+'.res'   &$


	readcol,file_sfh,age_sfh ,SFR_old , lmassold, time_you ,SFR_you, lmassyou, SFR_sfh, mass_sfh, lbt
	readcol,file_sfh,x,x,x,x,x,x,x,x,x, x, age,x,tau, format='a,a,a,a,a,a,a,a,a,a,a,a,a', numline=1

	age_bf(j)=float(age)
	tau_bf(j)=float(tau)
	log_tau_bf(j)=alog10(tau_bf(j)*1.e6)

	age_model(j,*)=age_sfh
	sfr_model(j,*)=sfr_sfh
	mass_model(j,*)=mass_sfh
	lbt_model(j,*)=lbt
	;ssfr_model(j,*)=ssfr_sfh

	now(j)=closest(lbt, 0)  
	mass_bf(j)=mass_sfh(now(j))


	 peak=where(sfr_sfh eq max(sfr_sfh))    ; where SFR is maximal
	 age_peak(j)=age_sfh(peak(0)) 
	 sfr_peak(j)=sfr_sfh(peak(0))


	deltat(j)=age_bf(j)-age_peak(j)         ; time since SFR peak
	age_uz(j)=getage(z)                     ; age universe at obs z

	age_upeak(j)=age_uz(j)-deltat(j)        ;age univ at SFR peak
	age_uform(j)=age_uz(j)-age_bf(j)        ;age univ at gal formation
	z_peak(j)=getredshift(age_upeak(j))     ;z univ at SFR peak
	z_form(j)=getredshift(age_uform(j))     ;z univ at SFR peak

endif  else begin


print, id, ' has no good solutions'

age_peak(j)=-1.
sfr_peak(j)=-1.
deltat(j)=-1.
age_uz(j)=-1.
age_upeak(j)=-1.
z_peak(j)=-1.
age_uform(j)=-1.
z_form(j)=-1.

print,  age_peak(j)

endelse

endfor





odd=where(age_uz - age_bf lt 0.1 or z_peak gt 10. and z_peak ne -1.)
w=where(age_uz - age_bf gt 0.1 and  z_peak lt 10. and z_peak ne -1.)



;--------mass---------------
m1=where((age_uz-age_bf) gt 0.1 and  z_peak lt 10. and mass_bf le 10.5 and z_peak ne -1.)                       ;46
m2=where((age_uz-age_bf) gt 0.1  and  z_peak lt 10. and mass_bf gt 10.5 and mass_bf le 10.8 and z_peak ne -1.)  ;36
m3=where((age_uz-age_bf) gt 0.1  and  z_peak lt 10. and mass_bf gt 10.8 and z_peak ne -1.)                      ;22



stop

;============mean===========================

val_z1=mean(z_peak(m1))   ;z at SFR peak
val_z2=mean(z_peak(m2))
val_z3=mean(z_peak(m3))

val_au1=mean(age_upeak(m1)) ; Age univ. at z peak
val_au2=mean(age_upeak(m2))
val_au3=mean(age_upeak(m3))


val_zf1=mean(z_form(m1))    ;z formation galaxy
val_zf2=mean(z_form(m2))
val_zf3=mean(z_form(m3))

val_auf1=mean(age_uform(m1)) ; age univ. formation galaxy
val_auf2=mean(age_uform(m2))
val_auf3=mean(age_uform(m3))


val_a1=mean(age_peak(m1)) ; age galaxy SFR peak
val_a2=mean(age_peak(m2))
val_a3=mean(age_peak(m3))

val_s1=mean(sfr_peak(m1)) ; SFR max
val_s2=mean(sfr_peak(m2))
val_s3=mean(sfr_peak(m3))

val_t1=mean(tau_bf(m1)) ; tau
val_t2=mean(tau_bf(m2))
val_t3=mean(tau_bf(m3))


;============median===========================

;val_z1=median(z_peak(m1))   ;z at SFR peak
;val_z2=median(z_peak(m2))
;val_z3=median(z_peak(m3))

;val_au1=median(age_upeak(m1)) ; Age univ. at z peak
;val_au2=median(age_upeak(m2))
;val_au3=median(age_upeak(m3))


;val_zf1=median(z_form(m1))    ;z formation galaxy
;val_zf2=median(z_form(m2))
;val_zf3=median(z_form(m3))

;val_auf1=median(age_uform(m1)) ; age univ. formation galaxy
;val_auf2=median(age_uform(m2))
;val_auf3=median(age_uform(m3))


;val_a1=median(age_peak(m1)) ; age galaxy SFR peak
;val_a2=median(age_peak(m2))
;val_a3=median(age_peak(m3))

val_s1=median(sfr_peak(m1)) ; SFR max
val_s2=median(sfr_peak(m2))
val_s3=median(sfr_peak(m3))

val_t1=median(tau_bf(m1)) ; tau
val_t2=median(tau_bf(m2))
val_t3=median(tau_bf(m3))


;=========SFH mean=========


tt=findgen(1000)/100
t=10^tt                                 ; tiempo sin log, eje de las x en la integral

zz=getredshift(tt)                      ;redshift at each t



tau1=val_t1*1.e6
tau2=val_t2*1.e6
tau3=val_t3*1.e6


sfr1_aux=(2.71828)^(-t/tau1)*t
sfr2_aux=(2.71828)^(-t/tau2)*t
sfr3_aux=(2.71828)^(-t/tau3)*t

sfr1=(sfr1_aux*val_s1)/max(sfr1_aux)  ; normalizamos al maximo de SFR en cada bin de masa
sfr2=(sfr2_aux*val_s2)/max(sfr2_aux) 
sfr3=(sfr3_aux*val_s3)/max(sfr3_aux)  

;-----------------------------
sfr1_max=max(sfr1)            ;comprobamos con datos medios
peak1=where(sfr1 eq sfr1_max)
tmax1=t(peak1)

sfr2_max=max(sfr2)
peak2=where(sfr2 eq sfr2_max)
tmax2=t(peak2)

sfr3_max=max(sfr3)
peak3=where(sfr3 eq sfr3_max)
tmax3=t(peak3)
;------------------------------------

;z range 1.0-1.5

x1=5.
x2=getage(1.5)

y1=5.
y2=280


str_tau1=strcompress(string(val_t1,'(i3)'),/remove_all)
str_tau2=strcompress(string(val_t2,'(i3)'),/remove_all)
str_tau3=strcompress(string(val_t3,'(i3)'),/remove_all)

str_tau1='80'
str_tau2='200'
str_tau3='250'


;===================PLOT==================================

print, 'EMPIEZA PLOT SFH medias mass'
stop


set_plot,'ps'
device, filename=dir+'/synthesizer_offset/files/SFH_mean_'+code+'_'+model+'_'+imf+''+met_string+'_mass_bins_median.eps',/col, /portrait, /times

trange=[2., 5.]

;z range


plot, t*1.e-9+val_auf1, sfr1, yrange=[0,y2],xrange=trange, xtitle='Age of Universe [Gyr]', ytitle='SFR!dSED!N ['+Msun+'/yr]', thick=4., xthick=4., ythick=4., charsize=1.5, charthick=4., /nodata, xstyle=9,/ystyle,  position=[0.17, 0.15, 0.9, 0.87], ytickinterval=100


cgColorFill, [x1, x2, x2, x1], [y1, y1, y2, y2],  COLOR=cgcolor('ryb4')

oplot, t*1.e-9+val_auf1, sfr1, col=4, thick=4.
;oplot, [tmax1*1.e-9+val_auf1, tmax1*1.e-9+val_auf1],[0,sfr1_max], col=4, line=2, thick=3
;oplot, [ tmax1*1.e-9+val_auf1,  tmax1*1.e-9+val_auf1+val_t1*1.e-3],[sfr1_max/2., sfr1_max/2.], col=4, thick=5

oplot, t*1.e-9+val_auf2, sfr2, col=3, thick=4.
;oplot, [tmax2*1.e-9+val_auf2, tmax2*1.e-9+val_auf2],[0,sfr2_max], col=3, line=2, thick=3
;oplot, [ tmax2*1.e-9+val_auf2,  tmax2*1.e-9+val_auf2+val_t2*1.e-3],[sfr2_max/2., sfr2_max/2.], col=3, thick=5


oplot, t*1.e-9+val_auf3, sfr3, col=2, thick=4.
;oplot, [tmax3*1.e-9+val_auf3, tmax3*1.e-9+val_auf3],[0,sfr3_max], col=2, line=2, thick=3
;oplot, [ tmax3*1.e-9+val_auf3,  tmax3*1.e-9+val_auf3+val_t3*1.e-3],[sfr3_max/2., sfr3_max/2.], col=2, thick=5


xyouts, 0.2,0.83,/norm,  'log(M/'+Msun+')=10.0-10.5, <!4s!3>='+str_tau1+' Myr, <!S!A-!R!Nt!dM!N>=1.1 Gyr', charsize=1., charthick=3., col=4
xyouts, 0.2,0.78,/norm, 'log(M/'+Msun+')=10.5-10.8, <!4s!3>='+str_tau2+' Myr, <!S!A-!R!Nt!dM!N>=1.4 Gyr', charsize=1., charthick=3., col=3
xyouts, 0.2,0.73,/norm,  'log(M/'+Msun+')=10.8-11.5, <!4s!3>='+str_tau3+' Myr, <!S!A-!R!Nt!dM!N>=1.6 Gyr', charsize=1., charthick=3., col=2


xval=[getage(5.), getage(4.),getage(3.),getage(2.),getage(1.5), getage(1.)]
lb=[5, 4, 3, 2, 1.5, 1]


axis,xaxis=1,xr=trange,xst=1,xtickv=xval,xtickname=string(lb,f='(F4.1)'), xticks=n_elements(xval)-1,chars=1., chart=4, xthick=4.


xyouts, [0.43], [0.95],' redshift ', /normal, charthick=4, charsize=1.5, color=cgcolor("Black")
	

plot, t*1.e-9+val_auf1, sfr1, yrange=[0,y2],xrange=trange, thick=4., xthick=4., ythick=4., charsize=1.5, charthick=4., /nodata, xstyle=9,/ystyle, position=[0.17, 0.15, 0.9, 0.87], /noerase, ytickinterval=100


device,/close & set_plot,'x'



print, 'FIN PLOT'
stop




;-------age------
;m1=where((age_uz-age_bf) gt 0.1 and  z_peak lt 10. and age_bf le 1.)                   ;26
;m2=where((age_uz-age_bf) gt 0.1  and  z_peak lt 10. and age_bf gt 1. and age_bf le 2.) ;45
;m3=where((age_uz-age_bf) gt 0.1  and  z_peak lt 10. and age_bf gt 2.)                  ;33

;-------age MW------
m1=where((age_uz-age_bf) gt 0.1 and  z_peak lt 10. and age_mw_gy le 1.)                       ;39
m2=where((age_uz-age_bf) gt 0.1  and  z_peak lt 10. and age_mw_gy gt 1. and age_mw_gy le 2.)  ;50
m3=where((age_uz-age_bf) gt 0.1  and  z_peak lt 10. and age_mw_gy gt 2.)                      ;15



stop

;============mean===========================

val_z1=mean(z_peak(m1))   ;z at SFR peak
val_z2=mean(z_peak(m2))
val_z3=mean(z_peak(m3))

val_au1=mean(age_upeak(m1)) ; Age univ. at z peak
val_au2=mean(age_upeak(m2))
val_au3=mean(age_upeak(m3))


val_zf1=mean(z_form(m1))    ;z formation galaxy
val_zf2=mean(z_form(m2))
val_zf3=mean(z_form(m3))

val_auf1=mean(age_uform(m1)) ; age univ. formation galaxy
val_auf2=mean(age_uform(m2))
val_auf3=mean(age_uform(m3))


val_a1=mean(age_peak(m1)) ; age galaxy SFR peak
val_a2=mean(age_peak(m2))
val_a3=mean(age_peak(m3))

val_s1=mean(sfr_peak(m1)) ; SFR max
val_s2=mean(sfr_peak(m2))
val_s3=mean(sfr_peak(m3))

val_t1=mean(tau_bf(m1)) ; tau
val_t2=mean(tau_bf(m2))
val_t3=mean(tau_bf(m3))


;============median===========================

;val_z1=median(z_peak(m1))   ;z at SFR peak
;val_z2=median(z_peak(m2))
;val_z3=median(z_peak(m3))

;val_au1=median(age_upeak(m1)) ; Age univ. at z peak
;val_au2=median(age_upeak(m2))
;val_au3=median(age_upeak(m3))


;val_zf1=median(z_form(m1))    ;z formation galaxy
;val_zf2=median(z_form(m2))
;val_zf3=median(z_form(m3))

;val_auf1=median(age_uform(m1)) ; age univ. formation galaxy
;val_auf2=median(age_uform(m2))
;val_auf3=median(age_uform(m3))


;val_a1=median(age_peak(m1)) ; age galaxy SFR peak
;val_a2=median(age_peak(m2))
;val_a3=median(age_peak(m3))

val_s1=median(sfr_peak(m1)) ; SFR max
val_s2=median(sfr_peak(m2))
val_s3=median(sfr_peak(m3))

val_t1=median(tau_bf(m1)) ; tau
val_t2=median(tau_bf(m2))
val_t3=median(tau_bf(m3))


;=========SFH mean=========


tt=findgen(1000)/100
t=10^tt                                 ; tiempo sin log, eje de las x en la integral

zz=getredshift(tt)                      ;redshift at each t



tau1=val_t1*1.e6
tau2=val_t2*1.e6
tau3=val_t3*1.e6


sfr1_aux=(2.71828)^(-t/tau1)*t
sfr2_aux=(2.71828)^(-t/tau2)*t
sfr3_aux=(2.71828)^(-t/tau3)*t

sfr1=(sfr1_aux*val_s1)/max(sfr1_aux)  ; normalizamos al maximo de SFR en cada bin de masa
sfr2=(sfr2_aux*val_s2)/max(sfr2_aux) 
sfr3=(sfr3_aux*val_s3)/max(sfr3_aux)  

;-----------------------------
sfr1_max=max(sfr1)            ;comprobamos con datos medios
peak1=where(sfr1 eq sfr1_max)
tmax1=t(peak1)

sfr2_max=max(sfr2)
peak2=where(sfr2 eq sfr2_max)
tmax2=t(peak2)

sfr3_max=max(sfr3)
peak3=where(sfr3 eq sfr3_max)
tmax3=t(peak3)
;------------------------------------

;z range 1.0-1.5

x1=5.
x2=getage(1.5)

y1=0.
y2=250




str_tau1=strcompress(string(val_t1,'(i3)'),/remove_all)
str_tau2=strcompress(string(val_t2,'(i3)'),/remove_all)
str_tau3=strcompress(string(val_t3,'(i3)'),/remove_all)


str_tau1='60'
str_tau2='200'
str_tau3='400'

;===================PLOT==================================

print, 'EMPIEZA PLOT SFH medias ages'
stop


set_plot,'ps'
device, filename=dir+'/synthesizer_offset/files/SFH_mean_'+code+'_'+model+'_'+imf+''+met_string+'_ageMW_bins_mean_tau.eps',/col, /portrait, /times

trange=[1.5, 5.]

;z range


plot, t*1.e-9+val_auf1, sfr1, yrange=[0.,y2],xrange=trange, xtitle='', ytitle='SFR!dSED!N ['+Msun+'/yr]', thick=4., xthick=4., ythick=4., charsize=1.5, charthick=4., /nodata, xstyle=9, position=[0.17, 0.15, 0.9, 0.87]


cgColorFill, [x1, x2, x2, x1], [y1, y1, y2, y2],  COLOR=cgcolor('ryb4')

oplot, t*1.e-9+val_auf1, sfr1, col=4, thick=4.
oplot, t*1.e-9+val_auf2, sfr2, col=3, thick=4.
oplot, t*1.e-9+val_auf3, sfr3, col=2, thick=4.


xyouts, 0.2,0.8,/norm,  '!S!A-!R!Nt!dM!N < 1.0 Gyr, <!4s!3>='+str_tau1+' Myr, <log(M/'+Msun+')>=10.4', charsize=1., charthick=3., col=4
xyouts, 0.2,0.75,/norm, '!S!A-!R!Nt!dM!N = 1.0-2.0 Gyr, <!4s!3>='+str_tau2+' Myr, <log(M/'+Msun+')>=10.5', charsize=1., charthick=3., col=3
xyouts, 0.2,0.7,/norm,  '!S!A-!R!Nt!dM!N > 2.0 Gyr, <!4s!3>='+str_tau3+' Myr, <log(M/'+Msun+')>=10.7', charsize=1., charthick=3., col=2



xval=[getage(5.), getage(4.),getage(3.),getage(2.),getage(1.5), getage(1.)]
lb=[5, 4, 3, 2, 1.5, 1]


axis,xaxis=1,xr=trange,xst=1,xtickv=xval,xtickname=string(lb,f='(F4.1)'), xticks=n_elements(xval)-1,chars=1., chart=4, xthick=4.


xyouts, [0.43], [0.95],' redshift ', /normal, charthick=4, charsize=1.5, color=cgcolor("Black")
	

plot, t*1.e-9+val_auf1, sfr1, yrange=[0,y2],xrange=trange, xtitle='Age of Universe [Gyr]', ytitle='', thick=4., xthick=4., ythick=4., charsize=1.5, charthick=4., /nodata, xstyle=9, position=[0.17, 0.15, 0.9, 0.87], /noerase


device,/close & set_plot,'x'



print, 'FIN PLOT'
stop





END

