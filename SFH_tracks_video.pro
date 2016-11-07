PRO SFH_tracks_video

; .comp /home/hdominguez/idl/pro/new/SFH_tracks_video.pro

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
restore,'/home/hdominguez/Ceverino/SFH_104.save', /verbose


;selected galaxies
;------------------
readcol, dir+'/out_HDFN_SHARDS_iDR2gammaP01.20150305_z_1.00_1.50_104',num_rd, id_rd,RA_rd,DEC_rd,z_rd,q_rd,VJ_rd,UV_rd,mass_rd,sSFR_rd,SFR_UV_obs_rd,SFR_IR_rd, SFR_tot_rd , AV_rd,mag_H, mag_j,a,a,a,sel_string, Format='F, A,F,F,F,F,F,F,f,f,f,f,f,f,f,f,a,a,a,a', skipline=1
; same order than previous files



stop

;-----------------Different epochs (always same position in sfh files)----------------
a1=100  ; 0.1 Gyr
a2=500  ; 0.5
a3=1000 ; 1.
a4=2000 ; 2.
a5=3000 ; 3.
a6=4000 ; 4.
a7=5000 ; 5.


;wage=[100,200,500,1000,2000,3000,4000,5000]
;legend=['t!d0!N=0.1 Gyr','t!d0!N=0.2 Gyr','t!d0!N=0.5 Gyr','t!d0!N=1.0 Gyr','t!d0!N=2.0 Gyr','t!d0!N=3.0 Gyr', 't!d0!N=4.0 Gyr','t!d0!N=5.0 Gyr']
;y_pos=[0.45,0.4,0.35,0.3,0.25,0.2,0.15, 0.1]+0.025


;57 steps
wage1=[1,2,3,4,5,8,10,15, 20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,160,170,180,190]


wage0=[1,2,3,4,5,8]
wage1=(findgen(103))*5+10
wage2=(findgen(28)+11)*50
wage3=(findgen(30)+20)*100


;wage3=(findgen(25)+26)*100


wage=[wage0, wage1, wage2, wage3]
wtime=wage/1000. ;time in Gyr




delta_sfr=fltarr(n_elements(source), n_elements(wage)-1)
median_delta_sfr=fltarr(n_elements(wage)-1)
for i=0, n_elements(source)-1 do begin &$
for j=0, n_elements(wage)-2 do begin  &$
delta_SFR(i,j)=alog10(sfr_model(i,wage(j)))-alog10(sfr_model(i,wage(j+1)))  &$
median_delta_sfr(j)=median(abs(delta_sfr(*,j)))
endfor
endfor



stop

;------------------------------------------------------------



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



colorin=[cgcolor('Sky Blue'),cgcolor('Royal Blue'),cgcolor('BLU7'), cgcolor('teal'), cgcolor('Lime Green'),cgcolor('Orange'),cgcolor('Crimson'), cgcolor('RED8') ]

loadct,13
colorin=findgen(n_elements(wage))*1.5+5



;colorin=[cgcolor('Sky Blue'),cgcolor('Blue'), cgcolor('Grn5'),cgcolor('Orange'),cgcolor('Crimson')]



;===================PLOT AGES==================================

w=where(log_age gt 0)

print, 'Empieza Plot tracks'
stop

for j=0, n_elements(wage)-1 do begin   ; cicle in step
;for j=0, 10 do begin

print, j

num=j

set_plot,'ps'
device, filename=dir+'/synthesizer_offset/files/plots_MS/SFH_tracks_'+code+'_'+model+'_'+imf+''+met_string+'_'+strtrim(num,1)+'_step.eps',/col, /portrait, /times

;stop

;===============tracks===================================================

!p.multi=0

plot, mass_model(0,*), sfr_model(0,*),  xrange=[9.0, 11.5],yrange=[1.e-3, 1.e4], /ylog, xtitle='log (M/'+Msun+')', ytitle='SFR!dSED!N ['+Msun+'/yr]',thick=2., xthick=4., ythick=4., charthick=3., charsize=1.2, psym=1, /nodata, /xstyle, /ystyle, position=[0.12,0.12,0.95,0.95]


;-----MS-------
oplot, mass_ms, 10^sfr_ms0,  thick=9., col=cgcolor('Grey')
oplot, mass_ms, 10^sfr_ms1,  thick=9., col=cgcolor('Grey'), line=2
oplot, mass_ms, SFR_mse,  thick=9., col=cgcolor('Dark Grey')




for i =0, n_elements(w)-1 do begin &$
;tracks
colours
oplot, mass_model(w(i),0:wage(j)), sfr_model(w(i),0:wage(j)), col=13, thick=1., line=1 &$
endfor

loadct,13
for i =0, n_elements(w)-1 do begin &$
oplot, [mass_model(w(i),wage(j)),mass_model(w(i),wage(j))],  [sfr_model(w(i),wage(j)),sfr_model(w(i),wage(j))], psym=sym(1), symsize=0.8, col=colorin(j), thick=4. &$
legend_str=strcompress(string(wtime(j),'(f5.3)'),/remove_all)

xyouts,0.25 ,0.4,/norm ,'t='+legend_str+' Gyr', col=colorin(j), charthick=3.0, charsize=0.8
endfor






xyouts,10.3,10^3.6 , 'MS @ z=2.0 Speagle+14', charthick=2.5, col=cgcolor('dark Grey'), charsize=1.
xyouts,10.3, 10^3.3 ,'MS @ z=1.2 Speagle+14', charthick=2.5, col=cgcolor('dark Grey'), charsize=1.
xyouts,10.3, 10^3.0 ,'MS @ z=1.0 Elbaz+07', charthick=2.5, col=cgcolor('Dark Grey'), charsize=1.

oplot, [11.15,11.4],[10^3.6,10^3.6], thick=5., col=cgcolor('Grey') ,line=2
oplot, [11.15,11.4],[10^3.3,10^3.3], thick=5.,col=cgcolor('Grey')
oplot, [11.15,11.4],[10^3.0,10^3.0], thick=5., col=cgcolor('dark Grey')



plot, mass_model, sfr_model,  xrange=[9.0, 11.5],yrange=[1.e-3, 1.e4], /ylog, xtitle='log (M/'+Msun+')', ytitle='SFR!dSED!N ['+Msun+'/yr]',thick=2., xthick=4., ythick=4., charthick=3., charsize=1.2, psym=1, /nodata, /xstyle, /ystyle, position=[0.12,0.12,0.95,0.95], /noerase


device,/close & set_plot,'x'


endfor    ;cicle age step

stop






print, 'FIN PLOT'
stop

stop


END
