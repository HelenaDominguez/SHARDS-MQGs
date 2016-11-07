PRO IR_z_plot

;   .comp /home/hdominguez/idl/pro/new/IR_z_plot.pro

Msun='M'+sunsymbol()

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


low_ir_all=fltarr(n_elements(low_ir))
for i=0, n_elements(low_ir)-1 do low_ir_all(i)=where(id_plot(low_ir(i)) eq id)

stop

;error SFR_IR
err_SFR=0.12  ;rel
;err_SFR=1.38  ;med
err_z=0.005

print, 'EMPIEZA PLOT'
stop


;====================PLOT=====================================

  set_plot,'PS'
  Device, /Helvetica, /Color, Bits_per_Pixel=8, File=dir+'/synthesizer_offset/files/IR_z_paper_plot.eps', XSize=15, YSize=13


plot, z, sfr_ir_all, psym=3, /ylog, xrange=[0.95,1.55],yrange=[0.7, 200], /ystyle, xtitle='redshift',ytitle='SFR!dIR!N ['+Msun+'/yr]', thick=4., xthick=4., ythick=4., charsize=1.5, charthick=4, /xstyle

oplot, z, sfr_ir_all,  psym=sym(1), symsize=0.2

oplot, z(low_ir_all), sfr_ir_all(low_ir_all), psym=sym(1), col=2,symsize=0.8
oplot, [0.9,1.6],[3.5,50], col=4, thick=6

oploterror, [1.5,1.5],[5.,5.], [err_z*1.5, err_z*1.5], [err_sfr*5., err_sfr*5.], thick=3.
;oploterror, [1.5,1.5],[5.,5.], [err_z*1.5, err_z*1.5], [err_sfr, err_sfr], thick=3.

xyouts, 0.55,0.3, 'MIPS-faint limit (30 !4l!3Jy 5!4r!3)', charsize=1., charthick=3., col=4, /norm
xyouts,0.55, 0.25 ,'MIPS-faint sample',col=2, charsize=1., charthick=3., /norm
xyouts,0.55, 0.2 ,'MIPS-detected', charsize=1., charthick=3., /norm


;legend, ['MIPS-faint limit'], line=0, thick=6., color=4, box=0, pos=[0.65,0.4], charsize=1., charthick=3, /normal
;legend, ['MIPS-faint limit (30 !4l!3Jy @ 5!4r!3)'], box=0, pos=[0.7,0.35], charsize=.8, charthick=3, /normal
;legend, ['MIPS-faint sample'], psym=sym(1), col=2,symsize=0.8, box=0, pos=[0.65,0.3], charsize=1., charthick=3, /normal
;legend, ['MIPS-detected'],psym=sym(1), symsize=0.2 , box=0, pos=[0.65,0.25], charsize=1., charthick=3, /normal


 Device, /Close_File
 set_plot,'x'

print, 'FIN plot IR-z'


stop
END
