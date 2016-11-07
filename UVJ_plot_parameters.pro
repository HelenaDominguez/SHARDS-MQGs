PRO  UVJ_plot_parameters

;  .comp /home/hdominguez/idl/pro/new/UVJ_plot_parameters.pro

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
analyze_dir='analyze_rejoin4'

Msun='M'+sunsymbol()

filelist=dir+'/ID_SHARDS_candels_zphotoiris_spectra_files_104.lst'
readcol, filelist, source, zcomb, G141_file, G102_file, format='a, f, a, a'


;Most significant solutions
readcol, dir+'/synthesizer_offset/files/analyze_parameters_BC_txp_krou_most_significant_sol_104.lst',   shards_id , nclus , z ,  log_tau,  log_age, log_age_MW,  av , met, mass,ratio, format='a, f, f, f'


;Primary Solutions
;readcol, dir+'/synthesizer_offset/files/analyze_parameters_BC_txp_krou_flags_104.lst',   shards_id , nclus , z ,  log_tau,  log_age, log_age_MW,  av , met, mass, ratio, format='a, f, f, f'




age_gy=10^log_age*1.e-9
age_mw_gy=10^log_age_mw*1.e-9
tau_my=10^log_tau*1.e-6
;ps=where(nclus eq 0)

stop

deg=where(ratio lt 100.)
onesol=where(ratio eq 100.)

vect_symbol=intarr(n_elements(shards_id))
vect_symbol2=intarr(n_elements(shards_id))

vect_symbol(deg)=2
vect_symbol2(deg)=7

vect_symbol(onesol)=1
vect_symbol2(onesol)=6

;selected galaxies
;------------------
readcol, dir+'/out_HDFN_SHARDS_iDR2gammaP01.20150305_z_1.00_1.50_104',num_rd, id_rd,RA_rd,DEC_rd,z_rd,q_rd,VJ_rd,UV_rd,mass_rd,sSFR_rd,SFR_UV_obs_rd,SFR_IR_rd, SFR_tot_rd , AV_rd,mag_H, mag_j,a,a,a,sel_string, Format='F, A,F,F,F,F,F,F,f,f,f,f,f,f,f,f,a,a,a,a', skipline=1
; same order than previous files


restore, dir+'/synthesizer_offset/files/BestFit_SFH_BC_txp_krou_met_104_l100.sav', /verbose
;% RESTORE: Restored variable: ID_LIST.
;% RESTORE: Restored variable: AGE_BF.
;% RESTORE: Restored variable: TAU_BF.
;% RESTORE: Restored variable: MASS_BF.
;% RESTORE: Restored variable: SFR_BF.
;% RESTORE: Restored variable: SSFR_BF.
;% RESTORE: Restored variable: TQ.


;rename
sfr_now=sfr_bf
sfr_bf=sfr_l100
ssfr_now=ssfr_bf
ssfr_bf=ssfr_l100


stop

;==============quartiles===================



;OLD vs YOUNG
;==================


old=where(age_mw_gy ge 2.0)    ;15 ~ 15%
int=where(age_mw_gy ge 1.0 and age_mw_gy lt 2.0)  ;50
young=where(age_mw_gy lt 1.0)  ;39

qmass1=quartile(mass(young))
qmass2=quartile(mass(int))
qmass3=quartile(mass(old))


qage1=quartile(age_mw_gy(young))
qage2=quartile(age_mw_gy(int))
qage3=quartile(age_mw_gy(old))

qage01=quartile(age_gy(young))
qage02=quartile(age_gy(int))
qage03=quartile(age_gy(old))

qtau1=quartile(tau_my(young))
qtau2=quartile(tau_my(int))
qtau3=quartile(tau_my(old))


qav1=quartile(av(young))
qav2=quartile(av(int))
qav3=quartile(av(old))


qmet1=quartile(met(young))
qmet2=quartile(met(int))
qmet3=quartile(met(old))


qssfr1=quartile(ssfr_bf(young))
qssfr2=quartile(ssfr_bf(int))
qssfr3=quartile(ssfr_bf(old))


qtq1=quartile(tq(young))
qtq2=quartile(tq(int))
qtq3=quartile(tq(old))

qz1=quartile(zcomb(young))
qz2=quartile(zcomb(int))
qz3=quartile(zcomb(old))


stop


;UVJ vs sSFR
;==================


uvj=where(sel_string eq 'UVJ')    ;65
lowsf=where(sel_string eq 'sSFR')  ;39



qmass4=quartile(mass(uvj))
qmass5=quartile(mass(lowsf))

qage4=quartile(age_mw_gy(uvj))
qage5=quartile(age_mw_gy(lowsf))

qage04=quartile(age_gy(uvj))
qage05=quartile(age_gy(lowsf))


qtau4=quartile(tau_my(uvj))
qtau5=quartile(tau_my(lowsf))


qav4=quartile(av(uvj))
qav5=quartile(av(lowsf))


qmet4=quartile(met(uvj))
qmet5=quartile(met(lowsf))

qssfr4=quartile(ssfr_bf(uvj))
qssfr5=quartile(ssfr_bf(lowsf))

qtq4=quartile(tq(uvj))
qtq5=quartile(tq(lowsf))

qz4=quartile(zcomb(uvj))
qz5=quartile(zcomb(lowsf))


stop


;mass
;==================


m1=where(mass le 10.5)    ;45
m2=where(mass gt 10.5 and mass le 10.8) ;36
m3=where(mass gt 10.8) ;23


qmass6=quartile(mass(m1))
qmass7=quartile(mass(m2))
qmass8=quartile(mass(m3))

qage06=quartile(age_gy(m1))
qage07=quartile(age_gy(m2))
qage08=quartile(age_gy(m3))

qage6=quartile(age_mw_gy(m1))
qage7=quartile(age_mw_gy(m2))
qage8=quartile(age_mw_gy(m3))

qtau6=quartile(tau_my(m1))
qtau7=quartile(tau_my(m2))
qtau8=quartile(tau_my(m3))

qav6=quartile(av(m1))
qav7=quartile(av(m2))
qav8=quartile(av(m3))

qmet6=quartile(met(m1))
qmet7=quartile(met(m2))
qmet8=quartile(met(m3))

qssfr6=quartile(ssfr_bf(m1))
qssfr7=quartile(ssfr_bf(m2))
qssfr8=quartile(ssfr_bf(m3))

qtq6=quartile(tq(m1))
qtq7=quartile(tq(m2))
qtq8=quartile(tq(m3))


qz6=quartile(zcomb(m1))
qz7=quartile(zcomb(m2))
qz8=quartile(zcomb(m3))


stop

; ==========Gaussian Noise===============

   sigma = 0.025
  
   arrayx = RANDOMN(seed, 1000, /normal)
   arrayy = RANDOMN(seed, 1000, /normal)
   xuvj= arrayx * sigma + vj_rd
   yuvj=arrayy * sigma + uv_rd


;---------------------UVJ region----------------------------------

a=[0.8,1.6]
b=[1.3,2.0]

a1=[0., 0.8]
b1=[1.3, 1.3]

a2=[1.6, 1.6]
b2=[2.0, 2.5]

a3=[1., 1.]
b3=[2.5, 1.5]


;-----------young-old separation--------------

       C0 = 0.59
        xjoin = (C0 - 2.85) / (-1.25-0.88)
        yjoin = -1.25 * xjoin + 2.85

        xjoin1 = !x.crange(0)
        yjoin1 =  -1.25 * xjoin1 + 2.85





;-------------positions---------------------

; 2x2
p1=[0.1,0.5, 0.5, 0.9]
p11=[0.15,0.53, 0.45, 0.55]

p2=[0.5,0.5, 0.9, 0.9]
p22=[0.55,0.53, 0.85, 0.55]

p3=[0.1,0.1, 0.5, 0.5]
p33=[0.15,0.13, 0.45, 0.15]

p4=[0.5,0.1, 0.9, 0.5]
p44=[0.55,0.13, 0.85, 0.15]


; 2x3
;----------------------
m=0.08  ; margen
p1=[m,m+2*(1-2*m)/3, 1/2., 1-m]
p2=[1/2.,m+2*(1-2*m)/3, 1-m, 1-m]

p3=[m,m+(1-2*m)/3, 1/2., m+2*(1-2*m)/3]
p4=[1/2.,m+(1-2*m)/3, 1-m, m+2*(1-2*m)/3]

p5=[m,m, 1/2., m+(1-2*m)/3]
p6=[1/2.,m, 1-m, m+(1-2*m)/3]


;3 bins in age
;--------------


wage1=where(age_mw_gy lt 1.)                          ; 39
wage2=where(age_mw_gy ge 1. and age_mw_gy lt 2.)      ; 48
wage3=where(age_mw_gy ge 2.)                          ; 17



;3 bins in tau
;--------------


wtau1=where(tau_my lt 80.0)                       ;47
wtau2=where(tau_my ge 80. and tau_my lt 400.)     ;43
wtau3=where(tau_my ge 400.)                       ;14



;3 bins in av
;--------------


wav1=where(av lt .5)                  ;42
wav2=where(av ge .5 and av lt 1.)     ;42
wav3=where(av ge 1.)                  ;20



;3 bins in met
;--------------

wm1=where(met lt -0.1) ;24
wm2=where(met eq 0.)   ;43
wm3=where(met ge .3)   ;37


;3 bins in tq
;--------------

wtq1=where(tq lt 0.5)                  ;24
wtq2=where(tq ge 0.5 and tq lt 1.0 )   ;45
wtq3=where(tq ge 1.0)                  ;35


;3 bins in ssfr
;--------------

limit=where(alog10(ssfr_bf) lt -5. )   

ssfr=ssfr_bf
arrayx = RANDOMN(seed, n_elements(limit), /normal)
ssfr(limit)=10^(-5.+arrayx*0.1)

ws1=where(alog10(ssfr) lt -4.)                             ;45
ws2=where(alog10(ssfr)  ge -4. and alog10(ssfr) lt -2. )   ;20
ws3=where(alog10(ssfr)  ge -2.0)                           ;39



;sizes for plot

ss=fltarr(N_elements(mass))
ss[0:N_elements(mass)-1]=.8

  




print, 'empieza UVJ PLOT para '+code+'_'+model+'_'+imf+''+met_string+'
stop



psfile=strtrim(dir+'/synthesizer_offset/files/UVJ_3bins_'+code+'_'+model+'_'+imf+''+met_string+'_PrimarySolutions_6panels_l100_new.eps',2)
set_plot, 'ps'
device, filename=psfile,/portrait,/times, /color, xsize=9,ysize=12, /inch

;============age  1 x 1 =====================

colorin=fltarr(n_elements(age_mw_gy))



colorin1=cgcolor('BLUe')
colorin2=cgcolor('GRN5')
colorin3=cgcolor('RED6')



wp1=(wage1)  ;39
wp2=(wage2)  ;48
wp3=(wage3)  ;17


colorin=strarr(n_elements(xuvj))
colorin(wp1)=cgcolor('BLUe')
colorin(wp2)=cgcolor('GRN5')
colorin(wp3)=cgcolor('RED6')


plot, xuvj, yuvj, psym=3, xrange=[0.80,1.9], yrange=[1.1, 2.3], ytitle='U - V ',  /xstyle, /ystyle, thick=2., charthick=3.,  charsize=1.2, xthick=4., ythick=4., position=p1, xtickformat='(A1)'




for i=0, n_elements(xuvj)-1 do begin &$
oplot,  [xuvj(i),xuvj(i)], [yuvj(i),yuvj(i)], psym=sym(vect_symbol(i)), symsize=ss(i),col=colorin(i) &$
oplot,  [xuvj(i),xuvj(i)], [yuvj(i),yuvj(i)], psym=sym(vect_symbol2(i)), symsize=ss(i),col=0, thick=2. &$
endfor


;for i=0, n_elements(wp1)-1 do begin &$
;oplot,  [xuvj(wp1(i)),xuvj(wp1(i))], [yuvj(wp1(i)),yuvj(wp1(i))], psym=sym(1), symsize=ss(wp1(i)),col=colorin1 &$
;oplot,  [xuvj(wp1(i)),xuvj(wp1(i))], [yuvj(wp1(i)),yuvj(wp1(i))], psym=sym(6), symsize=ss(wp1(i)),col=0, thick=2. &$
;endfor



;for i=0, n_elements(wp2)-1 do begin &$
;oplot,  [xuvj(wp2(i)),xuvj(wp2(i))], [yuvj(wp2(i)),yuvj(wp2(i))], psym=sym(1), symsize=ss(wp2(i)),col=colorin2 &$
;oplot,  [xuvj(wp2(i)),xuvj(wp2(i))], [yuvj(wp2(i)),yuvj(wp2(i))], psym=sym(6), symsize=ss(wp2(i)),col=0, thick=2. &$
;endfor



;for i=0, n_elements(wp3)-1 do begin &$
;oplot,  [xuvj(wp3(i)),xuvj(wp3(i))], [yuvj(wp3(i)),yuvj(wp3(i))], psym=sym(1), symsize=ss(wp3(i)),col=colorin3 &$
;oplot,  [xuvj(wp3(i)),xuvj(wp3(i))], [yuvj(wp3(i)),yuvj(wp3(i))], psym=sym(6), symsize=ss(wp3(i)),col=0, thick=2. &$
;endfor


colours ;to be in black!
oplot, a, b, thick=3., col=13
oplot, a1, b1, thick=3., col=13
oplot, a2, b2, thick=3., col=13
;oplot, a3, b3, line=2, thick=3.
;oplot,[xjoin,xjoin1],[yjoin,yjoin1],line=2,color=0,thick=2

xyouts,1.4, 1.4,'!S!A-!R!Nt!dM!N [Gyr] < 1.0',  charthick=3.,  charsize=.8, col=colorin1
xyouts,1.4, 1.3,'1.0 < !S!A-!R!Nt!dM!N [Gyr] < 2.0 ',  charthick=3.,  charsize=.8, col=colorin2
xyouts,1.4, 1.2,'!S!A-!R!Nt!dM!N [Gyr] > 2.0 ',  charthick=3.,  charsize=.8, col=colorin3

legend, ['1 solution'], psym=sym(1), color=cgcolor("black"), box=0, pos=[0.1,0.92], charsize=1., charthick=2, /normal
legend, ['> 1 solutions'], psym=sym(2), color=cgcolor("black"), box=0, pos=[0.1,0.9], charsize=1., charthick=2, /normal

;============Tau  1 X 2 =====================

;colorin1=cgcolor('GRN3')
;colorin2=cgcolor('GRN5')
;colorin3=cgcolor(' GRN8')


wp1=(wtau3)  ;54
wp2=(wtau2)  ;37
wp3=(wtau1)  ;22


colorin=strarr(n_elements(xuvj))
colorin(wp1)=cgcolor('BLUe')
colorin(wp2)=cgcolor('GRN5')
colorin(wp3)=cgcolor('RED6')



plot, xuvj, yuvj, psym=3, xrange=[0.80,1.9], yrange=[1.1, 2.3], /xstyle, /ystyle, thick=2., charthick=3., charsize=1.2,  xthick=4., ythick=4., position=p2 , xtickformat='(A1)',  ytickformat='(A1)',/noerase  




for i=0, n_elements(xuvj)-1 do begin &$
oplot,  [xuvj(i),xuvj(i)], [yuvj(i),yuvj(i)], psym=sym(vect_symbol(i)), symsize=ss(i),col=colorin(i) &$
oplot,  [xuvj(i),xuvj(i)], [yuvj(i),yuvj(i)], psym=sym(vect_symbol2(i)), symsize=ss(i),col=0, thick=2. &$
endfor




colours ;to be in black!
oplot, a, b, thick=3., col=13
oplot, a1, b1, thick=3., col=13
oplot, a2, b2, thick=3., col=13
;oplot, a3, b3, line=2, thick=3.

; oplot,[xjoin,xjoin1],[yjoin,yjoin1],line=2,color=0,thick=2

xyouts,1.4, 1.4,'!4s!3 [Myr] < 80     ',  charthick=3.,  charsize=.8, col=colorin3
xyouts,1.4, 1.3,'80 < !4s!3 [Myr] < 400 ',  charthick=3.,  charsize=.8, col=colorin2
xyouts,1.4, 1.2,'!4s!3 [Myr] > 400     ',  charthick=3.,  charsize=.8, col=colorin1


;============AV  2 x 1=====================

;colorin1=cgcolor('ORG3')
;colorin2=cgcolor('ORG6')
;colorin3=cgcolor(' ORG8')


wp1=wav1
wp2=wav2  
wp3=wav3  


colorin=strarr(n_elements(xuvj))
colorin(wp1)=cgcolor('BLUe')
colorin(wp2)=cgcolor('GRN5')
colorin(wp3)=cgcolor('RED6')

plot, xuvj, yuvj, psym=3, xrange=[0.80,1.9], yrange=[1.1, 2.3],   xtitle='' ,ytitle='U - V ', /xstyle, /ystyle, thick=2., charthick=3.,charsize=1.2, xthick=4., ythick=4., position=p3,  xtickformat='(A1)', /noerase  



for i=0, n_elements(xuvj)-1 do begin &$
oplot,  [xuvj(i),xuvj(i)], [yuvj(i),yuvj(i)], psym=sym(vect_symbol(i)), symsize=ss(i),col=colorin(i) &$
oplot,  [xuvj(i),xuvj(i)], [yuvj(i),yuvj(i)], psym=sym(vect_symbol2(i)), symsize=ss(i),col=0, thick=2. &$
endfor



colours ;to be in black!
oplot, a, b, thick=3., col=13
oplot, a1, b1, thick=3., col=13
oplot, a2, b2, thick=3., col=13
;oplot, a3, b3, line=2, thick=3.
; oplot,[xjoin,xjoin1],[yjoin,yjoin1],line=2,color=0,thick=2


xyouts,1.4, 1.4,' A!dV!N [mag] < 0.5       ',  charthick=3.,  charsize=.8, col=colorin1
xyouts,1.4, 1.3,' 0.5 < A!dV!N [mag]< 1.0',  charthick=3.,  charsize=.8, col=colorin2
xyouts,1.4, 1.2,' A!dV!N [mag] > 1.0        ',  charthick=3.,  charsize=.8, col=colorin3



;===========Metallicity 2 x 2 =====================


;colorin1=cgcolor('Plum')
;colorin2=cgcolor('Dark Orchid')
;colorin3=cgcolor('Maroon')

wp1=wm1
wp2=wm2   
wp3=wm3


colorin=strarr(n_elements(xuvj))
colorin(wp1)=cgcolor('BLUe')
colorin(wp2)=cgcolor('GRN5')
colorin(wp3)=cgcolor('RED6')


plot, xuvj, yuvj, psym=3, xrange=[0.80,1.9], yrange=[1.1, 2.3],   xtitle='' , /xstyle, /ystyle, thick=2., charthick=3.,charsize=1.2, xthick=4., ythick=4., position=p4,  ytickformat='(A1)',    xtickformat='(A1)', /noerase 




for i=0, n_elements(xuvj)-1 do begin &$
oplot,  [xuvj(i),xuvj(i)], [yuvj(i),yuvj(i)], psym=sym(vect_symbol(i)), symsize=ss(i),col=colorin(i) &$
oplot,  [xuvj(i),xuvj(i)], [yuvj(i),yuvj(i)], psym=sym(vect_symbol2(i)), symsize=ss(i),col=0, thick=2. &$
endfor




colours ;to be in black!
oplot, a, b, thick=3., col=13
oplot, a1, b1, thick=3., col=13
oplot, a2, b2, thick=3., col=13



xyouts,1.4, 1.4,'Z = 0.4 '+Zsun+'',  charthick=3.,  charsize=.8, col=colorin1
xyouts,1.4, 1.3,'Z = '+Zsun+'',  charthick=3.,  charsize=.8, col=colorin2
xyouts,1.4, 1.2,'Z = 2.5 '+Zsun+'',  charthick=3.,  charsize=.8, col=colorin3

;============sSFR 3 x 1=====================



;colorin1=cgcolor('BLU3')
;colorin2=cgcolor('BLU6')
;colorin3=cgcolor('BLU8')

 

wp1=ws3
wp2=ws2  
wp3=ws1 


colorin=strarr(n_elements(xuvj))
colorin(wp1)=cgcolor('BLUe')
colorin(wp2)=cgcolor('GRN5')
colorin(wp3)=cgcolor('RED6')

plot, xuvj, yuvj, psym=3, xrange=[0.80,1.9], yrange=[1.1, 2.3],   xtitle='V - J ' ,ytitle='U - V ', /xstyle, /ystyle, thick=2., charthick=3.,charsize=1.2, xthick=4., ythick=4., position=p5,  /noerase  



for i=0, n_elements(xuvj)-1 do begin &$
oplot,  [xuvj(i),xuvj(i)], [yuvj(i),yuvj(i)], psym=sym(vect_symbol(i)), symsize=ss(i),col=colorin(i) &$
oplot,  [xuvj(i),xuvj(i)], [yuvj(i),yuvj(i)], psym=sym(vect_symbol2(i)), symsize=ss(i),col=0, thick=2. &$
endfor





colours ;to be in black!
oplot, a, b, thick=3., col=13
oplot, a1, b1, thick=3., col=13
oplot, a2, b2, thick=3., col=13




xyouts,1.3, 1.4,'sSFR!dSED!N [Gyr!u-1!N] < 10!u-4!N',  charthick=3.,  charsize=.8, col=colorin3
xyouts,1.3, 1.3,'10!u-4!N < sSFR!dSED!N [Gyr!u-1!N] < 10!u-2!N',  charthick=3.,  charsize=.8, col=colorin2
xyouts,1.3, 1.2,'sSFR!dSED!N [Gyr!u-1!N] > 10!u-2!N',  charthick=3.,  charsize=.8, col=colorin1


;=========== Quenching time 3 x 2=====================

;colorin1=cgcolor('RYB3')
;colorin2=cgcolor('RYB2')
;colorin3=cgcolor('RYB1')




wp1=wtq1
wp2=wtq2   
wp3=wtq3


colorin=strarr(n_elements(xuvj))
colorin(wp1)=cgcolor('BLUe')
colorin(wp2)=cgcolor('GRN5')
colorin(wp3)=cgcolor('RED6')


plot, xuvj, yuvj, psym=3, xrange=[0.80,1.9], yrange=[1.1, 2.3],   xtitle='V - J' , /xstyle, /ystyle, thick=2., charthick=3.,charsize=1.2, xthick=4., ythick=4., position=p6,  ytickformat='(A1)',/noerase 



for i=0, n_elements(xuvj)-1 do begin &$
oplot,  [xuvj(i),xuvj(i)], [yuvj(i),yuvj(i)], psym=sym(vect_symbol(i)), symsize=ss(i),col=colorin(i) &$
oplot,  [xuvj(i),xuvj(i)], [yuvj(i),yuvj(i)], psym=sym(vect_symbol2(i)), symsize=ss(i),col=0, thick=2. &$
endfor


colours ;to be in black!
oplot, a, b, thick=3., col=13
oplot, a1, b1, thick=3., col=13
oplot, a2, b2, thick=3., col=13


xyouts,1.4, 1.4,'t!dq!N [Gyr] < 0.5',  charthick=3.,  charsize=.8, col=colorin1
xyouts,1.4, 1.3,'0.5 < t!dq!N [Gyr] < 1.0 ',  charthick=3.,  charsize=.8, col=colorin2
xyouts,1.4, 1.2,'t!dq!N [Gyr] > 1.0  ',  charthick=3.,  charsize=.8, col=colorin3



 Device, /Close_File
 set_plot,'x'





print, 'FIN plot'




stop





END
