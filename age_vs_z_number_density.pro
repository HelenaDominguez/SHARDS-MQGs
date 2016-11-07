PRO age_vs_z_number_density


;  .comp /home/hdominguez/idl/pro/new/age_vs_z_number_density.pro

code='BC'
model='txp'
imf='krou'
;met_string='_met'
met_string=''
nsim=300

Msun='M'+sunsymbol()
Zsun='Z'+sunsymbol()



; ----
; Data

; ----

dir='/home/hdominguez/SHARDS/DR2_gamma/redanddead_gamma'
Msun='M'+sunsymbol()

filelist=dir+'/ID_SHARDS_candels_zphotoiris_spectra_files_104.lst'
readcol, filelist, source, zcomb, G141_file, G102_file, format='a, f, a, a'

readcol, dir+'/synthesizer_offset/files/analyze_parameters_BC_txp_krou_flags_104.lst',   shards_id , nclus , z ,  tau,  age, age_MW,  av , met, mass, format='a, f, f, f'

age_mw_gy=10^age_mw*1.e-9
age_gy=10^age*1.e-9


;====== Choi values=======
readcol, dir+'/synthesizer_offset/files/choi_ages.txt', age_U_choi, age_choi

age_u_choi1=age_u_choi[0:3]
age_u_choi2=age_u_choi[4:7]
age_u_choi3=age_u_choi[8:11]

age_choi1=age_choi[0:3]
age_choi2=age_choi[4:7]
age_choi3=age_choi[8:11]


;====== mendel plot=======
readcol, dir+'/synthesizer_offset/files/Mendel_ages.txt', age_U_mendel, age_mendel

;Mendel
age_u_m=[age_u_mendel(0), age_u_mendel(0)]
age_m=[age_mendel(0), age_mendel(0)]

;Whitaker
age_u_w=[age_u_mendel(1), age_u_mendel(1)]
age_w=[age_mendel(1), age_mendel(1)]

age_u_w=[getage(1.64),getage(1.64)]  ;we use median z from W13 paper insetad of Mendel plot

;Ondera
age_u_o=[age_u_mendel(2), age_u_mendel(2)]
age_o=[age_mendel(2), age_mendel(2)]

;Belli
age_u_be=[age_u_mendel(3), age_u_mendel(4)]
age_be=[age_mendel(3), age_mendel(4)]

;Schiavom
age_u_s=[age_u_mendel(5), age_u_mendel(5)]
age_s=[age_mendel(5), age_mendel(5)]

stop

;=========Barro 2015==============

age_b=[0.75,1.07,1.14,1.31]
age_u_b=[getage(1.6), getage(1.6),getage(1.6),getage(1.6)]

; very similar ages, we use mean
age_b=[0.75,1.2]
age_u_b=[getage(1.6), getage(1.6)]

;=======Toft 2012 =========

z_toft=mean([2.03,1.82,1.80,2.19])
log_age_toft=mean([8.9,8.6,9.18,9.24])

age_toft=10^log_age_toft*1.e-9
age_u_toft=getage(z_toft)

;======Van de Sande==========

z_vds=mean([1.71,1.60,2.02,1.94,1.44])
log_age_vds=mean([8.74,8.96,9.41,8.74,8.82])

age_vds=10^log_age_vds*1.e-9
age_u_vds=getage(z_vds)


;========Krogager=====

age_u_k=getage(2.35)
age_k=0.98

;=======Beddregal========

age_u_bed=getage(1.5)
;age_bed=[3.1,1.5] ;ages on-off red sequence
age_bed=[2.0,2.0] ;median age


;=====my sample 3 bins=============


age1=where(age_mw_gy lt 1.)                       ;38
age2=where(age_mw_gy ge 1. and age_mw_gy lt 2.  ) ;50
age3=where(age_mw_gy ge 2.) 


stop


x=[getage(median(zcomb(age1))),getage(median(zcomb(age2))),getage(median(zcomb(age3)))]
;y=[mean(age_mw_gy(age1)),mean(age_mw_gy(age2)),mean(age_mw_gy(age3))]
y=[median(age_mw_gy(age1)),median(age_mw_gy(age2)),median(age_mw_gy(age3))]
yerr=[stdev(age_mw_gy(age1)),stdev(age_mw_gy(age2)),stdev(age_mw_gy(age3))]

stop

;=================tau 400 =====================

tau_my=10^tau*1.e-6

w=where(tau_my gt 400 and  tau_my lt 401)
;w=where(tau_my eq max(tau_my))


id=shards_id(w)
z_obs=z(w)

num=fix(nclus(w)+1)

file_sfh=dir+'/synthesizer_offset/res/analyze_rejoin5/clusters/'+id+'.BCall.bc2003.stelib.'+model+'.'+imf+'.sfh.clust'+strtrim(num,1)+'.res'  

readcol,file_sfh,age_sfh ,SFR_old , lmassold, time_you ,SFR_you, lmassyou, SFR_sfh, mass_sfh, lbt
readcol,file_sfh,a,a,a,a,a,a,a,a,a, a, age,a,tau, format='a,a,a,a,a,a,a,a,a,a,a,a,a', numline=1

age_bf=float(age)
tau_bf=float(tau)
log_age=alog10(age_sfh*1.e9)
log_tau_bf=alog10(tau_bf*1.e6)



age_obs=getage(z_obs)
age_univ=getage(z_obs)+lbt

;stop

ww=where(log_age gt 0 and log_age lt 9.9)
log_age_obs_mw=fltarr(n_elements(ww))
for i=0, n_elements(ww)-1 do log_age_obs_mw(i)=age_mw_txp(log_age(ww(i)),log_tau_bf(0))

age_obs_mw=10^log_age_obs_mw*1.e-9


;stop

now=closest(lbt, 0) 

age_form=getage(z_obs)-age_sfh(now)
z_form=getredshift(age_form)         ; zform=2.0

;=========================================


print, 'empieza Plot Number density vs z'
stop



set_plot,'ps'
device, filename=dir+'/synthesizer_offset/files/Age_vs_age_universe.eps',/col, /portrait, /times


plot, x, y, psym=sym(1), xtitle='Age of Universe [Gyr]', ytitle='Stellar Age [Gyr]', yrange=[0., 7.], xrange=[2., 12.],  symsize=2,/ystyle, thick=6, charthick=3,	xthick=4, ythick=4, chars=1.5, /nodata, xstyle=9,position=[0.17, 0.15, 0.9, 0.87], ytickinterval=1., yminor=1

cgColorFill, [getage(1.5),getage(1.5),getage(1.),getage(1.),getage(1.5)], [0.,7.,7.,0.,0.], col=cgcolor('Light Yellow')



;----lines------
oplot, [4.5,11.5],[0, 7.],line=1, thick=4, col=13
oplot, [3.2,10.2],[0, 7.],line=3, thick=4, col=13
oplot, [2.2,9.2],[0, 7.],line=4, thick=4, col=13
oplot, age_univ, age_obs_mw, line=2,col=13, thick=4.


;----my data------
plotsym, 3, thick=5, /fill
oplot, [x(0),x(0)], [y(0),y(0)], psym=8, symsize=2.0/1.3, col=cgcolor('crimson')
oplot, [x(1),x(1)], [y(1),y(1)], psym=8, symsize=2.0, col=cgcolor('Crimson')
oplot, [x(2),x(2)], [y(2),y(2)], psym=8, symsize=2.0/3.1, col=cgcolor('Crimson')



plotsym, 3, thick=5
oplot, [x(0),x(0)], [y(0),y(0)], psym=8, symsize=1.8
oplot, [x(1),x(1)], [y(1),y(1)], psym=8, symsize=2.0
oplot, [x(2),x(2)], [y(2),y(2)], psym=8, symsize=1.2


;choi
oplot, age_u_choi1, age_choi1, psym=sym(5), col=cgcolor('BLU5'), thick=2., symsize=1.2
oplot, age_u_choi1, age_choi1, psym=6,  thick=1., symsize=1.2

oplot, age_u_choi2, age_choi2, psym=sym(5), col=cgcolor('BLU3'), thick=3., symsize=1.2
oplot, age_u_choi2, age_choi2, psym=6,  thick=1., symsize=1.2

oplot, age_u_choi3, age_choi3, psym=sym(5), col=cgcolor('BLU7'), thick=4., symsize=1.2
oplot, age_u_choi3, age_choi3, psym=6,  thick=1., symsize=1.2



;others
oplot, age_u_m, age_m, psym=sym(1),thick=2., symsize=1.5, col=cgcolor('Maroon')
oplot, age_u_m, age_m, psym=sym(6),thick=2., symsize=1.5


oplot, age_u_o, age_o, psym=sym(3), thick=2., symsize=1.5, col=cgcolor('dark khaki')
oplot, age_u_o, age_o, psym=sym(8), thick=2., symsize=1.5

oplot, age_u_be, age_be, psym=sym(4),thick=2., symsize=1.5, col=cgcolor('Sea Green')
oplot, age_u_be, age_be, psym=sym(9),thick=2., symsize=1.5


oplot, age_u_s, age_s, psym=sym(4), thick=2., symsize=1.5, col=cgcolor('TG5')
oplot, age_u_s, age_s, psym=sym(9), thick=2., symsize=1.5

oplot, age_u_b, age_b, psym=sym(5), thick=2., symsize=1.5, col=cgcolor('orange')
oplot, [age_u_b(0),age_u_b(0)] , [age_b(0),age_b(0)] , psym=sym(5), thick=2., symsize=1.5, col=cgcolor('ORG6')
oplot, age_u_b, age_b, psym=sym(10), thick=2., symsize=1.5


oplot, [age_u_toft,age_u_toft], [age_toft,age_toft], psym=sym(3), thick=2., symsize=1.5, col=cgcolor('Pink')
oplot, [age_u_toft,age_u_toft], [age_toft,age_toft], psym=sym(8), thick=2., symsize=1.5


oplot, [age_u_vds,age_u_vds], [age_vds,age_vds], psym=sym(2), thick=2., symsize=1.5, col=cgcolor('Indian Red')
oplot,[age_u_vds,age_u_vds], [age_vds,age_vds] , psym=sym(7), thick=2., symsize=1.5



oplot, [age_u_k,age_u_k], [age_k,age_k], psym=sym(1), thick=2., symsize=1.5, col=cgcolor('Pur6')
oplot,[age_u_k,age_u_k], [age_k,age_k] , psym=sym(6), thick=2., symsize=1.5



oplot, [age_u_bed,age_u_bed] , age_bed, psym=sym(1), thick=2., symsize=1.5, col=cgcolor('Green Yellow')
oplot, [age_u_bed,age_u_bed], age_bed, psym=sym(6), thick=2., symsize=1.5


oplot, age_u_w, age_w, psym=sym(4),thick=2., symsize=1.5, col=cgcolor('orchid')
oplot, age_u_w, age_w, psym=sym(9),thick=2., symsize=1.5



;====lines
oplot, [8.,8.8],[1.5,1.5],line=1,col=13, thick=4.
xyouts, 9.,1.5, '!4s!3=100 Myr, z!df!N=1.4', charthick=1.5, charsize=0.8

oplot, [8.,8.8],[1.1,1.1],line=2,col=13, thick=4.
xyouts, 9.,1.1, '!4s!3=400 Myr, z!df!N=2.0', charthick=1.5, charsize=0.8

oplot, [8.,8.8],[0.7,0.7],line=3,col=13, thick=4.
xyouts, 9.,0.7, '!4s!3=100 Myr, z!df!N=2.0', charthick=1.5, charsize=0.8

oplot, [8.,8.8],[0.3,0.3],line=4,col=13, thick=4.
xyouts, 9.,0.3, '!4s!3=100 Myr, z!df!N=3.0', charthick=1.5, charsize=0.8




;legend
;=========================

plotsym, 3, thick=5, /fill
legend, ['This work'], psym=8, color=cgcolor("crimson"), box=0, pos=[0.2,0.85], charsize=1., charthick=2, /normal
plotsym, 3, thick=2
legend, ['This work'], psym=8, box=0, pos=[0.2,0.85], charsize=1., charthick=1.5, /normal


legend, ['Krogager+2013'], psym=sym(1), color=cgcolor("PUR6"), box=0, pos=[0.2,0.80], charsize=1., charthick=2, /normal
legend, ['Toft+2012'], psym=sym(3), color=cgcolor("pink"), box=0, pos=[0.2,0.77], charsize=1., charthick=2, /normal
legend, ['Whitaker+2013'], psym=sym(4), color=cgcolor("orchid"), box=0, pos=[0.2,0.74], charsize=1., charthick=2, /normal
legend, ['Mendel+2015'], psym=sym(1), color=cgcolor("Maroon"), box=0, pos=[0.2,0.71], charsize=1., charthick=2, /normal
legend, ['Van de Sande+2013'], psym=sym(2), color=cgcolor("Indian red"), box=0, pos=[0.2,0.68], charsize=1., charthick=2, /normal
legend, ['Barro+2015'], psym=sym(5), color=cgcolor("Orange"), box=0, pos=[0.2,0.65], charsize=1., charthick=2, /normal
legend, ['Barro+2015, quenching'], psym=sym(5), color=cgcolor("Org6"), box=0, pos=[0.2,0.62], charsize=1., charthick=2, /normal
legend, ['Onodera+2015'], psym=sym(3), color=cgcolor("dark khaki"), box=0, pos=[0.2,0.59], charsize=1., charthick=2, /normal
legend, ['Bedregal+2013'], psym=sym(1), color=cgcolor("Green yellow"), box=0, pos=[0.2,0.56], charsize=1., charthick=2, /normal
legend, ['Belli+2015'], psym=sym(4), color=cgcolor("Sea Green"), box=0, pos=[0.2,0.53], charsize=1., charthick=2, /normal
legend, ['Schiavon+2006'], psym=sym(4), color=cgcolor("TG5"), box=0, pos=[0.2,0.50], charsize=1., charthick=2, /normal
legend, [''], psym=sym(5), color=cgcolor("BLU3"), box=0, pos=[0.2,0.47], charsize=1., charthick=2, /normal
legend, [''], psym=sym(5), color=cgcolor("BLU5"), box=0, pos=[0.23,0.47], charsize=1., charthick=2, /normal
legend, ['Choi+2014'], psym=sym(5), color=cgcolor("BLU7"), box=0, pos=[0.26,0.47], charsize=1., charthick=2, /normal





plot, x, y, psym=sym(1), xtitle='', ytitle='', yrange=[0., 7.], xrange=[2., 12.],  symsize=2,/ystyle, thick=6, charthick=3,	xthick=4, ythick=4, chars=1.5, /nodata, xstyle=9,position=[0.17, 0.15, 0.9, 0.87], ytickinterval=1., yminor=1, /noerase


xval=[getage(5.), getage(4.),getage(3.),getage(2.),getage(1.5), getage(1.), getage(.5),  getage(.2)]
lb=[5, 4, 3, 2, 1.5, 1, 0.5, 0.2]

trange=[2., 12.]

axis,xaxis=1,xr=trange,xst=1,xtickv=xval,xtickname=string(lb,f='(F4.1)'), xticks=n_elements(xval)-1,chars=1.2, chart=3, xthick=4.


xyouts, [0.43], [0.95],' redshift ', /normal, charthick=3, charsize=1.5, color=cgcolor("Black")
	


device,/close 
set_plot,'x'
stop



;area
area=111.60                   ; arcmin2 (RB)
area_deg=area*(1/3600.)       ;deg2
area_sky=41253.               ;deg2
portion=area_deg/area_sky


;vol @ diff. z
vol1=vcomoving(1.5, /Mpc)-vcomoving(1., /Mpc)
v1=vol1*portion  ;Mpc3



;----dens by ages------



age1=where(age_mw_gy lt 1.)                       ;38
age2=where(age_mw_gy ge 1. and age_mw_gy lt 2.  ) ;50
age3=where(age_mw_gy ge 2.) 


n1=n_elements(age1) ;26
n2=n_elements(age2) ;61
n3=n_elements(age3) ;23
;n4=n_elements(age4) ;10

dens1=n1/v1
dens2=n2/v1
dens3=n3/v1
;dens4=n4/v1
dens_tot=104./v1



;my sample 3 bins

x=[getage(1.2),getage(1.2),getage(1.2)]
y=[mean(age_mw_gy(age1)),mean(age_mw_gy(age2)),mean(age_mw_gy(age3))]
z=[dens1,dens2, dens3]
yerr=[stdev(age_mw_gy(age1)),stdev(age_mw_gy(age2)),stdev(age_mw_gy(age3))]




;Whitaker
;------------------

nw=177
volw=vcomoving(2.2, /Mpc)-vcomoving(1.4, /Mpc)

areaw=900.
areaw_deg=areaw*(1/3600.)   ;deg2
area_sky=41253.             ;deg2

portionw=areaw_deg/area_sky

vw=volw*portionw
densw=nw/vw
zw=[densw, densw]


;Whitaker
;------------------

nc=37000
volc=vcomoving(0.09, /Mpc)-vcomoving(0.07, /Mpc)
areac_deg=9400.
area_sky=41253.             ;deg2
portionc=areac_deg/area_sky

vc=volw*portionc
densc=nc/vc

zc=[densc, densc]




stop



END


