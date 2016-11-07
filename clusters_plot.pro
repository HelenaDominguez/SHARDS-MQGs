PRO clusters_plot


; .comp /home/hdominguez/idl/pro/new/clusters_plot.pro

code='BC'
model='txp'
imf='krou'
met_string='_met'
;met_string=''
nsim=300

Msun='M'+sunsymbol()
Zsun='Z'+sunsymbol()

dir='/home/hdominguez/SHARDS/DR2_gamma/redanddead_gamma'

analyze_dir='analyze_rejoin5'
index_dir='analyze_rejoin5'
;index_dir='analyze_index_tunned'


filelist=dir+'/ID_SHARDS_candels_zphotoiris_spectra_files_104.lst'
readcol, filelist, source, zcomb, G141_file, G102_file, format='a, f, a, a'


filelist=dir+'/ID_ra_dec_zspec_104_sql.lst'
readcol, filelist, source2, ra, dec, zspec, zspec_flag, zcomb2, format='a, f, f, f'


readcol, dir+'/synthesizer_offset/files/analyze_parameters_BC_txp_krou.lst',   shards_id , nclus , z ,  log_tau,  log_age, log_age_MW,  av , met, mass, format='a, f, f, f'

age_gy=10^log_age*1.e-9
age_mw_gy=10^log_age_mw*1.e-9
tau_my=10^log_tau*1.e-6


;================No cicle, only one galaxy plot===================

;id='SHARDS10008577'

id='SHARDS20000663'


if code eq 'BC' then  file_sfh=dir+'/synthesizer_offset/res/'+id+'_1pop.BCall.bc2003.stelib.csq.'+model+'.hr.'+imf+'_M_0.1_100.calz.res'

		data=leefile(file_sfh)

	
		nparam=n_elements(data(*,0))
	
        	name=string(data(1,*))
		met=float(data(3,*))
		log_tau=float(data(4,*))
		log_age=float(data(5,*))
		av=float(data(6,*))
		mass=float(data(nparam-7,*))
		zmod=float(data(nparam-1,*))
		



		; Simlulation to cluster assignment
                ;----------------------------------
		if code eq 'BC' then file_N=dir+'/synthesizer_offset/res/'+analyze_dir+'/clusters/'+id+'.BCall.bc2003.stelib.'+model+'.'+imf+'.clusters.dat'
                readcol, file_N, nsim, ncluster, /silent
	
        	ncl=ncluster+1    
		
		cl1=where(ncl eq 1)
		cl2=where(ncl eq 2)
		cl3=where(ncl eq 3)
		cl4=where(ncl eq 4)
		cl5=where(ncl eq 5)
		cl6=where(ncl eq 6)
		cl7=where(ncl eq 7)
		
	
		Num_cl=max(ncl)



;==========median cluster params==================

 clus_pars=fltarr(num_cl,5)

       for j=1, num_cl do begin  &$
        print, j

	if code eq 'BC' then file_cl=dir+'/synthesizer_offset/res/'+analyze_dir+'/clusters/'+id+'.lala.BCall.bc2003.stelib.txp.krou.clust'+strtrim(j,1)+'.res' 

                ; Cluster parameters
                ;---------------------- 
                nlines=7

		readcol, file_cl, tau_cl, age_cl,met_cl, av_cl, format='x,x,x,x,x,x,a,x,a,x,a,x,a', skipline=nlines, numline=1 

                       age_mw_cl=age_mw_txp(age_cl(0),tau_cl(0)) &$

	               
                       clus_pars(j-1,0)=tau_cl       &$
                       clus_pars(j-1,1)=age_cl       &$
                       clus_pars(j-1,2)=age_mw_cl    &$			
                       clus_pars(j-1,3)=av_cl        &$
		       clus_pars(j-1,4)=met_cl       &$

endfor

stop


;=========gaussian noise===========

sigma_tau=0.05
sigma_age=0.05
;sigma_age_mw=0.1
sigma_av=0.05
sigma_m=0.05



   array1 = RANDOMN(seed, N_elements(log_tau), /normal)   
   array2 = RANDOMN(seed, N_elements(log_tau), /normal)   
   array3 = RANDOMN(seed, N_elements(log_tau), /normal)   
   array4 = RANDOMN(seed, N_elements(log_tau), /normal)   
   array5 = RANDOMN(seed, N_elements(log_tau), /normal)  



   log_taug=array1*sigma_tau + log_tau
   log_ageg =array2*sigma_age + log_age
   avg  =array3*sigma_av  + av
   massg=array4*sigma_m + mass
  


age_gy=10^log_age*1.e-9
tau_my=10^log_tau*1.e-6


ageg_gy=10^log_ageg*1.e-9
taug_my=10^log_taug*1.e-6

age_clus=10^clus_pars(*,1)*1.e-9
tau_clus=10^clus_pars(*,0)*1.e-6


;====================density matrix============================0

   x1 = taug_my(cl1)
   y1 = ageg_gy(cl1)

   x2 = taug_my(cl2)
   y2 = ageg_gy(cl2)

   x3 = taug_my(cl3)
   y3 = ageg_gy(cl3)


;------cl1------------------

   xrange = [100, 300]
   yrange = [1., 3.]
   
   xbinsize = 20
   ybinsize = 0.2
   
 density1 = Hist_2D(x1, y1, Min1=xrange[0], Max1=xrange[1], Bin1=xbinsize,  Min2=yrange[0], Max2=yrange[1], Bin2=ybinsize) 
 maxdensity1=total(density1)


xcont1=findgen(n_elements(density1(*, 0)))*xbinsize+xrange(0)+xbinsize
ycont1=findgen(n_elements(density1(0,*)))*ybinsize+yrange(0)+ybinsize

xcont2=findgen(n_elements(density1(*, 0)))*xbinsize+xrange(0)
ycont2=findgen(n_elements(density1(0,*)))*ybinsize+yrange(0)

xcont_cl1=(xcont1+xcont2)/2
ycont_cl1=(ycont1+ycont2)/2


order=reverse(sort(density1))
density1_ord=density1(order)
density1_ord_norm=float(cumulate(density1_ord))/N_elements(x1)  ; matriz de densidades

ind1=fltarr(n_elements(order))
ind2=fltarr(n_elements(order))
for l=0, n_elements(order)-1 do begin &$
ind2(l)=order(l)/11 &$
ind1(l)=order(l)-ind2(l)*11 &$
endfor


array_ref=fltarr(11,11)
array_dens=fltarr(11,11)
for l=0, n_elements(order)-1 do begin &$
array_ref(ind1(l),ind2(l))=density1_ord(l) &$
array_dens(ind1(l),ind2(l))=density1_ord_norm(l) &$
endfor


;-----cl2------------------------

   xrange = [200, 450]
   yrange = [2., 4.6]
   
   xbinsize = 40
   ybinsize = 0.4
 density2 = Hist_2D(x2, y2, Min1=xrange[0], Max1=xrange[1], Bin1=xbinsize,  Min2=yrange[0], Max2=yrange[1], Bin2=ybinsize) 


xcont1=findgen(n_elements(density2(*, 0)))*xbinsize+xrange(0)+xbinsize
ycont1=findgen(n_elements(density2(0,*)))*ybinsize+yrange(0)+ybinsize

xcont2=findgen(n_elements(density2(*, 0)))*xbinsize+xrange(0)
ycont2=findgen(n_elements(density2(0,*)))*ybinsize+yrange(0)

xcont_cl2=(xcont1+xcont2)/2
ycont_cl2=(ycont1+ycont2)/2


order=reverse(sort(density2))
density2_ord=density2(order)
density2_ord_norm=float(cumulate(density2_ord))/N_elements(x2)  ; matriz de densidades

ind1=fltarr(n_elements(order))
ind2=fltarr(n_elements(order))
for l=0, n_elements(order)-1 do begin &$
ind2(l)=order(l)/7 &$
ind1(l)=order(l)-ind2(l)*7 &$
endfor


array_ref2=fltarr(7,7)
array_dens2=fltarr(7,7)
for l=0, n_elements(order)-1 do begin &$
array_ref2(ind1(l),ind2(l))=density2_ord(l) &$
array_dens2(ind1(l),ind2(l))=density2_ord_norm(l) &$
endfor


;-----cl3------------------------

   xrange = [0, 70.]
   yrange = [1., 2.4]
   
   xbinsize = 10
   ybinsize = 0.2

 density3 = Hist_2D(x3, y3, Min1=xrange[0], Max1=xrange[1], Bin1=xbinsize,  Min2=yrange[0], Max2=yrange[1], Bin2=ybinsize) 


xcont1=findgen(n_elements(density3(*, 0)))*xbinsize+xrange(0)+xbinsize
ycont1=findgen(n_elements(density3(0,*)))*ybinsize+yrange(0)+ybinsize

xcont2=findgen(n_elements(density3(*, 0)))*xbinsize+xrange(0)
ycont2=findgen(n_elements(density3(0,*)))*ybinsize+yrange(0)

xcont_cl3=(xcont1+xcont2)/2
ycont_cl3=(ycont1+ycont2)/2


order=reverse(sort(density3))
density3_ord=density3(order)
density3_ord_norm=float(cumulate(density3_ord))/N_elements(x3)  ; matriz de densidades

ind1=fltarr(n_elements(order))
ind2=fltarr(n_elements(order))
for l=0, n_elements(order)-1 do begin &$
ind2(l)=order(l)/8 &$
ind1(l)=order(l)-ind2(l)*8 &$
endfor


array_ref3=fltarr(8,8)
array_dens3=fltarr(8,8)
for l=0, n_elements(order)-1 do begin &$
array_ref3(ind1(l),ind2(l))=density3_ord(l) &$
array_dens3(ind1(l),ind2(l))=density3_ord_norm(l) &$
endfor




stop


;===================PLOT==================================

	colour=[cgColor("Black"),cgColor("Red"),cgColor("Royal Blue"),cgColor("Orange"),cgColor("Sea Green"),cgColor("Turquoise"),cgColor("Blue Violet"),cgColor("Yellow"),cgColor("Violet"),cgColor("Grey")]
	
print, 'EMPIEZA PLOT'
stop

set_plot,'ps'
device, filename=dir+'/synthesizer_offset/files/Clusters_'+code+'_'+model+'_'+imf+''+met_string+'_'+id+'_new.eps',/col, /portrait, /times, xsize=7, ysize=7, /inch

!p.multi=0

xmin=0.
xmax=500.

ymin=1.
ymax=5.

plot, taug_my, ageg_gy,psym=3, xrange=[xmin, xmax], yrange=[ymin, ymax], /xstyle, /ystyle, xtitle='!4s!3 [Myr]', ytitle='t!d0!N [Gyr]', thick=4, charthick=4, xthick=4, ythick=4, chars=1.5, xtickinterval=100, position=[0.1,0.1,0.95,0.95], xminor=1, yminor=1




oplot, taug_my(cl1), ageg_gy(cl1),psym=sym(1), thick=4, symsize=0.5, col=colour(1)
oplot, taug_my(cl2), ageg_gy(cl2),psym=sym(1), thick=4, symsize=0.5,col=colour(2)
oplot, taug_my(cl3), ageg_gy(cl3),psym=sym(1), thick=4,symsize=0.5, col=colour(3)


   cgContour, array_dens, xcont_cl1, ycont_cl1, LEVELS=[0.68], /overplot, C_Colors=['black'], C_Annotation=['1!4r!3'], C_Thick=4, C_CharThick=4., c_labels=0

   cgContour, array_dens2, xcont_cl2, ycont_cl2, LEVELS=[0.70], /OnImage, C_Colors=['black'], C_Annotation=['1!4r!3'], C_Thick=4, C_CharThick=4, c_labels=0

        cgContour, array_dens3, xcont_cl3, ycont_cl3, LEVELS=[0.68], /OnImage, C_Colors=['black'], C_Annotation=['1!4r!3'], C_Thick=4, C_CharThick=4 , c_labels=0


oplot, tau_clus, age_clus,psym=sym(1), symsize=2., col=cgcolor("Sea Green")
oplot, tau_clus, age_clus,psym=sym(6), symsize=2.

i=28

ra_str=strcompress(string(ra(i),'(f7.3)'),/remove_all)
dec_str=strcompress(string(dec(i),'(f7.4)'),/remove_all)

;xyouts, [0.2], [0.85], 'SHARDS '+ra_str+'+'+dec_str+'', /normal, charthick=3, charsize=1.2, color=cgcolor("Black")
xyouts, [0.12], [0.87], 'SHARDSJ123620.3+6200844.3', /normal, charthick=3, charsize=1.5, color=cgcolor("Black") ;SHARDS20000663
xyouts, [0.12], [0.82], 'z!dsp!N='+strcompress(string(zcomb(i),'(f6.3)'),/remove_all), /normal, charthick=3, charsize=1.5, color=cgcolor("Black")



;xyouts, [0.2], [0.87], id, /normal, charthick=4, charsize=1.2, color=cgcolor("Black")
;xyouts, [0.2], [0.80], 'z='+strcompress(string(zcomb(i),'(f6.3)'),/remove_all), /normal, charthick=4, charsize=1.2, color=cgcolor("Black")
;xyouts, 0.5, 0.2, /norm, 'BC03 - txp - Krou - 3 Z', charsize=1.2, charthick=4.


;oplot, [300, 300],[1.5, 1.5], psym=sym(1), symsize=2., thick=5., col=cgcolor("Sea Green")
;oplot, [300, 300],[1.5, 1.5], psym=sym(6), symsize=2., thick=5.
xyouts, 250, 1.5, 'Cluster mean solutions', charsize=1.5, charthick=4., col=cgcolor("Sea Green")


device,/close 
set_plot,'x'


print, 'FIN PLOT'
stop


;===================PLOT histogram==================================


set_plot,'ps'
device, filename=dir+'/synthesizer/plots_mine/Clusters_'+code+'_'+model+'_'+imf+''+met_string+'_'+id+'_hist.eps',/col, /portrait, /times


;!p.multi=[0,1,2,0]

plothist, tau_my, bin=20., xtitle='!4s!3 (Myr)', ytitle='Number solutions', thick=4, charthick=4, xthick=4, ythick=4, chars=1., position=[0.1,0.58,0.95,0.95], /noerase

;plothist, taug_my, bin=2, /overplot, line=2.
plothist, taug_my(cl1), bin=20.,  col=colour(1), /overplot, thick=4
plothist, taug_my(cl2), bin=20.,  col=colour(2), /overplot, thick=4
;plothist, taug_my(cl3), bin=2.,  col=colour(3), /overplot, thick=4

oplot, [tau_clus(0), tau_clus(0)], [0,250], line=2, col=cgcolor("orange"), thick=4
oplot, [tau_clus(1), tau_clus(1)], [0,250], line=2, col=cgcolor("orange"), thick=4
;oplot, [tau_clus(2), tau_clus(2)], [0,150], line=2, col=cgcolor("orange"), thick=4



plothist, age_gy, bin=.2, xtitle='Age (Gyr)', ytitle='Number solutions', thick=4, charthick=4, xthick=4, ythick=4, chars=1., position=[0.1,0.1,0.95,0.47], /noerase
plothist, ageg_gy(cl1), bin=.2,  col=colour(1), /overplot, thick=4
plothist, ageg_gy(cl2), bin=.2,  col=colour(2), /overplot, thick=4
;plothist, ageg_gy(cl3), bin=.2,  col=colour(3), /overplot, thick=4

oplot, [age_clus(0), age_clus(0)], [0,250], line=2, col=cgcolor("orange"), thick=4
oplot, [age_clus(1), age_clus(1)], [0,250], line=2, col=cgcolor("orange"), thick=4
;oplot, [age_clus(2), age_clus(2)], [0,150], line=2, col=cgcolor("orange"), thick=4



device,/close 
set_plot,'x'

stop



stop

END
