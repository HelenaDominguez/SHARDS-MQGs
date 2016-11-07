PRO  SED_paper




;  .comp /home/hdominguez/idl/pro/new/SED_paper.pro

code='BC'
model='txp'
imf='krou'
;met_string='_met'
met_string=''
nsim=300

analyze_dir='analyze_rejoin5'

Msun='M'+sunsymbol()
Zsun='Z'+sunsymbol()



; ----
; Data

; ----

dir='/home/hdominguez/SHARDS/DR2_gamma/redanddead_gamma'
Msun='M'+sunsymbol()

filelist=dir+'/ID_SHARDS_candels_zphotoiris_spectra_files_104.lst'
readcol, filelist, source, zcomb, G141_file, G102_file, format='a, f, a, a'

filelist=dir+'/ID_ra_dec_zspec_104_sql.lst'
readcol, filelist, source2, ra, dec, zspec, zspec_flag, zcomb2, format='a, f, f, f'


;properties and flags
readcol, dir+'/synthesizer_offset/files/analyze_parameters_BC_txp_krou.lst',   shards_id , nclus , z_list ,  log_tau,  log_age, log_age_MW,  av , met, mass,ratio_file,  format='a,i, f, f, f'



;flag_d4, flag_mg, flag_deg, flag_sol

;Index
restore, dir+'/synthesizer_offset/files/index_BC03_txp_krou_104_d4_err.sav', /verbose

;F24 lim
readcol, dir+'/synthesizer_offset/files/F24_output_BC_txp_krou_met_135.lst', shards_id_24, i_24, nclus_24, z_24 ,Mass_24,AV_24,F24, IR_det, format='a, i, i, f'



print, 'Version different clusters, begining PLOT!  for '+code+' ' +model+'  '+imf+''+met_string+''



N_clus=fltarr(n_elements(source))
Nbands=fltarr(n_elements(source))

SN_g102=fltarr(n_elements(source))
SN_g141=fltarr(n_elements(source))


num=where(source eq 'SHARDS20000663')
;num=where(source eq 'SHARDS10003973')

stop


psfile=strtrim(dir+"/synthesizer_offset/files/SED_plot_clusters_"+code+"_"+model+"_"+imf+""+met_string+"_"+source(num)+".eps",2)


set_plot, 'ps'
device, filename=psfile, /portrait, /times, /color
!p.multi=0

for i=num(0), num(0) do begin
	id=source(i)
	z=zcomb(i)

	print, i
        print, i
        print, i
        print, i

	print, id, z


               
		;========solutions=======================

		w=where(shards_id eq source(i))
                ;w=5
                f24_sol=f24(w)
                ir_det_sol=ir_det(w)


               ;if flag_sol(w(0)) gt 0. then begin  ; plot good gal
               ;if flag_sol(w(0)) le 0 then begin  ; plot bad gal               
       	 	
        
       	 	
		; ----------
		; Photometry
		; ----------
		
		if code eq 'BC' then file=dir+'/synthesizer_offset/res/'+id+'_1pop.BCall.bc2003.stelib.csq.'+model+'.hr.'+imf+'_M_0.1_100.calz.res'
                if code eq 'M05' then  file=dir+'/synthesizer/res/'+id+'_1pop.M05.none.none.csq.exp.hr.krou_M_0.1_100.calz.res'
		
                strspawn='grep filter '+file+' > file_filt'
		spawn, strspawn
		readcol, 'file_filt', filtname, Eff_wl_micron, width, ind_temp, format='x,a,f,x,f,x,x,x,i', /silent

		Eff_wl=Eff_wl_micron*10000.
		Eff_wl_rf=Eff_wl_micron*10000./(1+z)
		width=width*10000./(1+z)
				
		cat_input_phot=dir+'/synthesizer_offset/data/observed/'+id+'.phot'
		dataphot=leefile(cat_input_phot)
	
		sz=size(dataphot)
		nfiltphot=dataphot(sz(1)-1)

		flux_pe=fltarr(nfiltphot)
		err_pe=fltarr(nfiltphot)
        	filt_pe=strarr(nfiltphot)

		ind1=where(dataphot eq 'flux:')
		ind2=where(dataphot eq 'err:')
		ind3=where(dataphot eq 'filt:')

		for j=0, nfiltphot-1 do begin &$		
			flux_pe(j)=dataphot(j+ind1+1) &$
			err_pe(j)=dataphot(j+ind2+1) &$
			filt_pe(j)=dataphot(j+ind3+1) &$
		endfor

		
                filt_used=where(ind_temp eq 1 and flux_pe gt 0 and err_pe gt 0)

		width_cat=width(filt_used)
		filt_cat=filtname(filt_used)		
		lambda_cat=Eff_wl(filt_used)
		lambda_cat_rf=Eff_wl_rf(filt_used)
		flux_cat=flux_pe(filt_used)
		err_cat=err_pe(filt_used)
                Nbands(i)=n_elements(flux_cat)

	err_cat_original=err_cat
        err_cat=3*err_cat
              
              print, 'Nbands synth=', n_elements(flux_cat)
              ;print, 'maximal wl difference=', min(lambda_cat_eazy(sort(lambda_cat_eazy))/lambda_cat(sort(lambda_cat)))
		




		; -------------------------------------------------
		; Flux model for each cluster
		; -------------------------------------------------
		

              

		; Simlulation to cluster assignment
                ;----------------------------------
		if code eq 'BC' then file_N=dir+'/synthesizer_offset/res/'+analyze_dir+'/clusters/'+id+'.BCall.bc2003.stelib.'+model+'.'+imf+'.clusters.dat'
		readcol, file_N, nsim, ncluster, /silent


              Nc=max(ncluster)
              ncl=nc+1     ;number clusters
              Num_cl=ncl
	        
              c1=where(ncluster eq 0)
              c2=where(ncluster eq 1)
	      c3=where(ncluster eq 2)
              c4=where(ncluster eq 3)
              c5=where(ncluster eq 4)
	      c6=where(ncluster eq 5)


              ratio=[n_elements(c1)/1000., n_elements(c2)/1000.,n_elements(c3)/1000.,n_elements(c4)/1000.,n_elements(c5)/1000.,n_elements(c6)/1000.]*100.


		; Best-fit model from synthesizer (originally in rest-frame)
                ;-------------------------------------------------------------
		if code eq 'BC' then file_bestfit_synth=dir+'/synthesizer_offset/res/plots/'+id+'_1pop.BCall.bc2003.stelib.csq.'+model+'.hr.'+imf+'_M_0.1_100.calz.'+id+'r.ANGSTROMS_SED'
                if code eq 'M05' then  file_bestfit_synth=dir+'/synthesizer_offset/res/plots/'+id+'_1pop.M05.none.none.csq.exp.hr.krou_M_0.1_100.calz.'+id+'r.ANGSTROMS_SED'
           
	        print,'reading ', file_bestfit_synth
		readcol, file_bestfit_synth, lambda_synth_rf, x,x,x,x,x,x,x,flux_temp_synth, /silent
	
		lambda_synth=lambda_synth_rf*(1+z)
		

                 ;other clusters
	         ;----------------------------------------------------
                 ;-----------------------------------------------------
		dist=lumdist(z)               					; en Mpc
		dist_cm=DOUBLE(dist*1E6*3.086E18)   			; Mpc --> pc --> cm 
	
		j=1
		if code eq 'BC' then file_cl=dir+'/synthesizer_offset/res/'+analyze_dir+'/clusters/'+id+'.lala.BCall.bc2003.stelib.txp.krou.clust'+strtrim(j,1)+'.res' 
                if code eq 'M05' then file_cl=dir+'/synthesizer_offset/res/'+analyze_dir+'/clusters/'+id+'.M05.none.none.exp.krou.clust'+strtrim(j,1)+'.res' 
		readcol, file_cl, wl_cl1_RF, Lnu_cl1, /silent
		readcol, file_cl, wl_cl_RF, Lnu_cl, /silent
              
		clus=fltarr(num_cl+1,n_elements(Lnu_cl))
		;clus(0,*) = wl_cl_rf*(1+z)                         ; observed wl       
                clus(0,*) = wl_cl_rf
		clus(1,*) = Lnu_cl
	
		clus_fnu=fltarr(num_cl+1,n_elements(Lnu_cl))	
		fnu=lnu_cl/(4*!pi*(dist_cm)^2)    				; erg s-1 Hz-1 cm-2  F=L/(4pid2)
		fnu_Jy=fnu*1.e23                   				; Jy
		fnu_muJy=fnu_Jy*1.e6              				; muJy
	        
                ;z=z(0) ;transform array in float
	        a=where(wl_cl_rf*(1+z) gt 5000. and wl_cl_rf*(1+z) lt 10000.)      ; whole range
		factor=median(flux_temp_synth(a)/fnu_muJy(a))     ; flux conversion factor from the data
		
		clus_fnu(0,*) = wl_cl_rf*(1+z)   ;observed
		clus_fnu(1,*) = fnu_muJy*factor

                ; Cluster parameters
                ;---------------------- 
		

                ;if i eq 0 then nlines=6
                ;if i gt 0 then nlines=7
                ;nlines=7
               ; if i eq 132 then begin  nlines=6 
                ;  endif else begin
                    nlines=7
                 ; endelse

                ;stop


		readcol, file_cl, tau_cl, age_cl,met_cl, av_cl, format='x,x,x,x,x,x,a,x,a,x,a,x,a', skipline=nlines, numline=1 

                       age_mw_cl=age_mw_txp(age_cl(0),tau_cl(0)) &$

	               clus_pars=fltarr(num_cl,5)
                       clus_pars(0,0)=tau_cl       &$
                       clus_pars(0,1)=age_cl       &$
                       clus_pars(0,2)=age_mw_cl    &$			
                       clus_pars(0,3)=av_cl        &$
		       clus_pars(0,4)=met_cl       &$

	        ;-----------------------------------------------------
		
		for j=2,Num_cl do begin &$
	
			if code eq 'BC' then file_cl=dir+'/synthesizer_offset/res/'+analyze_dir+'/clusters/'+id+'.lala.BCall.bc2003.stelib.txp.krou.clust'+strtrim(j,1)+'.res' 
                        if code eq 'M05' then file_cl=dir+'/synthesizer_offset/res/'+analyze_dir+'/clusters/'+id+'.M05.none.none.exp.krou.clust'+strtrim(j,1)+'.res'  &$
			readcol, file_cl, wl_cl_RF, Lnu_cl, /silent &$
			
			clus(j,*) = Lnu_cl &$
			
			fnu=lnu_cl/(4*!pi*(dist_cm)^2)     &$
			fnu_Jy=fnu*1.e23                   &$
			fnu_muJy=fnu_Jy*1.e6               &$

  
			wl_clus_order=wl_cl_rf(sort(wl_cl_rf))   &$
			fnu_muJy_order=fnu_muJy(sort(wl_cl_rf))   &$
			
			
                       a=where(wl_cl_rf*(1+z) gt 5000. and wl_cl_rf*(1+z) lt 100000.)      ; whole range
		       factor=median(flux_temp_synth(a)/fnu_muJy_order(a))   &$
		       clus_fnu(j,*) = fnu_muJy*factor               &$
			

                ; Cluster parameters
                ;----------------------
		readcol, file_cl, tau_cl, age_cl,met_cl, av_cl, format='x,x,x,x,x,x,a,x,a,x,a,x,a', skipline=7, numline=1 

                       age_mw_cl=age_mw_txp(age_cl(0),tau_cl(0)) &$
                
                       clus_pars(j-1,0)=tau_cl       &$
                       clus_pars(j-1,1)=age_cl       &$
                       clus_pars(j-1,2)=age_mw_cl    &$			
                       clus_pars(j-1,3)=av_cl        &$
		       clus_pars(j-1,4)=met_cl       &$

	
		endfor
	
               ; Masses from  analyze_parameters file
               ;---------------------------------

               w=where(shards_id eq id)
               ;w=where(shards_id eq id(0))

               mass_cl=n_elements(w)
               mass_cl=mass(w)
		
               ;============= filter names =========================
                
               filt_cat_str=strarr(n_elements(lambda_cat))
               filt_cat_str2=strarr(n_elements(lambda_cat))

               for k=0,n_elements(lambda_cat)-1 do begin &$
			
			tmp=strsplit(filt_cat(k),'_',/extract) &$
			filt_cat_str(k)=tmp(0) &$
			if n_elements(tmp) gt 1. then filt_cat_str2(k)=tmp(1) &$
				
		endfor
	
               ;============= shards filters=========================
		shards_b=where(filt_cat_str eq 'shards')
			
                width_cat_shards=width_cat(shards_b)		
		filt_cat_shards=filt_cat(shards_b)
                lambda_cat_shards=lambda_cat(shards_b)
                lambda_cat_shards_RF=lambda_cat(shards_b)/(1+z)
                flux_cat_shards=flux_cat(shards_b)
                err_cat_shards=err_cat(shards_b)
                 
		
               ;============= G102 filters=========================
                
		g102_b=where(filt_cat_str eq 'g102')

	if n_elements(g102_b) gt 1 then begin  
	
                width_cat_g102=width_cat(g102_b)		
		filt_cat_g102=filt_cat(g102_b)
                lambda_cat_g102=lambda_cat(g102_b)
                lambda_cat_g102_RF=lambda_cat(g102_b)/(1+z)
                flux_cat_g102=flux_cat(g102_b)
                err_cat_g102=err_cat(g102_b)


                 ; S/N
                 ;-----------------------------
		 g102_SN=3*flux_cat_g102/err_cat_g102  ;x3 as we have increased err to 3 sig.
		 SN_g102(i)=median(g102_SN)


		;upper y lower limits con errores
                ;-------------------------------
	        flux_g102_high=flux_cat_g102+err_cat_g102
		flux_g102_low=flux_cat_g102-err_cat_g102
                 
	endif                  		     			
		
               ;============= G141 filters=========================
                
		g141_b=where(filt_cat_str eq 'g141')

	if n_elements(g141_b) gt 1 then begin 	

                width_cat_g141=width_cat(g141_b)		
		filt_cat_g141=filt_cat(g141_b)
                lambda_cat_g141=lambda_cat(g141_b)
                lambda_cat_g141_RF=lambda_cat(g141_b)/(1+z)
                flux_cat_g141=flux_cat(g141_b)
                err_cat_g141=err_cat(g141_b)
                       
		; S/N
                ;-----------------------------
		g141_SN=3*flux_cat_g141/err_cat_g141  ;x3 as we have increased err to 3 sig.
		SN_g141(i)=median(g141_SN)


		;upper y lower limits con errores
                ;-------------------------------
	        flux_g141_high=flux_cat_g141+err_cat_g141
		flux_g141_low=flux_cat_g141-err_cat_g141
                 	
             		
        endif
		
               ;============= Broad Band=========================
                
		broad_b=where(filt_cat_str ne  'shards' and filt_cat_str ne  'g102'  and filt_cat_str ne  'g141')
			
                width_cat_bb=width_cat(broad_b)		
		filt_cat_bb=filt_cat(broad_b)
                lambda_cat_bb=lambda_cat(broad_b)
                lambda_cat_bb_RF=lambda_cat(broad_b)/(1+z)
                flux_cat_bb=flux_cat(broad_b)
                err_cat_bb=err_cat(broad_b)
                 			



			;==============================================
			;
			;               OUPUT PLOT
			; 
			;==============================================
			;==============================================
                        cc=0.07
                        bb=0.08
                        aa=0.24

			p1=[0.11, 0.10, 0.53,0.92]
			p2=[0.65, cc+2*aa+2*bb, 0.98, cc+3*aa+2*bb] 
              		p3=[0.65,cc+aa+bb , 0.98, cc+2*aa+bb]
			p4=[0.65,cc , 0.98, cc+aa]
 			
			;p1=[0.10, 0.10, 0.60,0.90]
			;p2=[0.72, 0.65, 0.95, 0.90] 
              		;p3=[0.72, 0.35, 0.95, 0.60]
			;p4=[0.72, 0.10, 0.95, 0.30]


		        xrangex=[0.3, 10.]

			ymax=max(flux_cat)*2
                        ymin=min(flux_cat)*0.4
			;ymax=100.
			;ymin=0.01
                        

			yrange=ymax-ymin

			min_mag=-2.5*alog10(ymax)+23.9
			max_mag=-2.5*alog10(ymin)+23.9


      

	colour=[cgColor("Red"),cgColor("Royal Blue"),cgColor("Orange"),cgColor("Turquoise"),cgColor("grn5"),cgColor("Blue Violet"),cgColor("Yellow"),cgColor("Violet"),cgColor("Grey")]



                        ;================= SED PLOT =======================
	
			plot, lambda_synth*1.e-4, flux_temp_synth,  /ylog,/xlog, xrange=xrangex,yrange=[ymin, ymax], $
			xtitle='!4k!3!dobs!N [!4l!3m]', ytitle='F!d!7m!5!N [!7l!3Jy]', $
			ystyle=9, thick=6, charthick=3, $
			xthick=3, ythick=3, chars=1., xstyle=9, /nodata, position=p1

			;Emission lines
			;------------------------------------
			oplot, [6562.,6562.]*1.e-4*(1+z),[ymin,ymax], line=2, col=13, thick=2.   ; Halpha
			oplot, [4861.,4861.]*1.e-4*(1+z),[ymin,ymax], line=2, col=13 , thick=2.   ; Hbeta
			oplot, [3726.,3726.]*1.e-4*(1+z),[ymin,ymax], line=2, col=13, thick=2.    ; OII



			; plotting different models
			;----------------------------
			for p=1,num_cl do oplot, clus_fnu(0,*)*1.e-4, clus_fnu(p,*), col=colour(p-1), thick=2
			oplot, clus_fnu(0,*)*1.e-4, clus_fnu(1,*), col=colour(0), thick=2
			

                        ;----------Broad Band---------
			oploterror, lambda_cat_bb*1.e-4, flux_cat_bb,err_cat_bb, thick=4, errcol=cgcolor("Black"), psym=3, symsize=0.5, /nohat
			oplot, lambda_cat_bb*1.e-4, flux_cat_bb, thick=2, col=cgcolor("White"), psym=sym(4), symsize=1.2
			oplot, lambda_cat_bb*1.e-4, flux_cat_bb, thick=4, col=cgcolor("Black"), psym=sym(9), symsize=1.2
                        

                        ;----------SHARDS---------
			oploterror, lambda_cat_shards*1.e-4, flux_cat_shards,err_cat_shards, thick=4, errcol=cgcolor("Black"), psym=3, symsize=0.5, /nohat
			oplot, lambda_cat_shards*1.e-4, flux_cat_shards, thick=2, col=cgcolor("White"), psym=sym(1), symsize=.7
			oplot, lambda_cat_shards*1.e-4, flux_cat_shards, thick=8, col=cgcolor("Black"), psym=sym(6), symsize=.7
                        
	    
                        ;----------G141 & G142 --------- 
		;	if n_elements(g102_b) gt 1 then oplot, lambda_cat_g102*1.e-4, flux_cat_g102, thick=4, col=cgcolor("Lime Green"), psym=10
		;	if n_elements(g141_b) gt 1 then oplot, lambda_cat_g141*1.e-4, flux_cat_g141, thick=4, col=cgcolor("Green"), psym=10

			if n_elements(g102_b) gt 1 then oplot, lambda_cat_g102*1.e-4, flux_cat_g102, thick=4, col=cgcolor("Black"), psym=10
			if n_elements(g141_b) gt 1 then oplot, lambda_cat_g141*1.e-4, flux_cat_g141, thick=4, col=cgcolor("Charcoal"), psym=10

                        ; labels
			;---------------------------- 
			   if ir_det_sol(0) eq -1. then str_ir='F(24) < 30 !4l!3Jy (5!4r!3) '  
                           if ir_det_sol(0) eq 1. then  str_ir='F(24) > 30 !4l!3Jy (5!4r!3) '  ; IR detection

			ra_str=strcompress(string(ra(i),'(f7.3)'),/remove_all)
			dec_str=strcompress(string(dec(i),'(f7.4)'),/remove_all)

                        ;xyouts, [0.15], [0.87], 'SHARDS '+ra_str+'+'+dec_str+'', /normal, charthick=3, charsize=0.9, color=cgcolor("Black")
                        if id  eq 'SHARDS10003973' then xyouts, [0.13], [0.87], 'SHARDSJ123727.9+622034.7', /normal, charthick=3, charsize=0.8, color=cgcolor("Black") ;SHARDS10003937
                        if id  eq 'SHARDS20000663' then xyouts, [0.13], [0.87], 'SHARDSJ123620.3+620844.3', /normal, charthick=3, charsize=0.8, color=cgcolor("Black") ;SHARDS20000663
 
			xyouts, [0.13], [0.84], 'z!dsp!N='+strcompress(string(z,'(f6.3)'),/remove_all), /normal, charthick=3, charsize=0.8, color=cgcolor("Black")
			xyouts, [0.13], [0.81], str_ir, /normal, charthick=3, charsize=0.8, color=cgcolor("Black")
	         	xyouts, [0.13], [0.78], 'N bands='+strcompress(string(Nbands(i),'(i4)'),/remove_all), /normal, charthick=3, charsize=0.8, color=cgcolor("Black")

                        youts=[0.33, 0.3, 0.27, 0.24, 0.21]-0.12
			xouts=[ 0.25, 0.29, 0.33, 0.37, 0.41, 0.44, 0.47, 0.51, 0.54, 0.59]-0.04


                        xyouts, xouts(0)+0.01, [0.27], 'Z', /normal, charthick=2, charsize=0.5, color=cgcolor("Black")
                        xyouts, xouts(1)+0.01, [0.27], '!4s!3', /normal, charthick=2, charsize=0.5, color=cgcolor("Black")
                        xyouts, xouts(2), [0.27], '  t!d0!N', /normal, charthick=2, charsize=0.5, color=cgcolor("Black")
                        xyouts, xouts(3)+0.01, [0.27], '!S!A-!R!Nt!dM!N', /normal, charthick=2, charsize=0.5, color=cgcolor("Black")
                        xyouts, xouts(4)+0.01, [0.27], 'A!dV!N', /normal, charthick=2, charsize=0.5, color=cgcolor("Black")
                        xyouts, xouts(5), [0.27], 'log(M)', /normal, charthick=2, charsize=0.5, color=cgcolor("Black")
                        xyouts, xouts(6)+0.01, [0.27], 'F!dpr!N(24)', /normal, charthick=2, charsize=0.5, color=cgcolor("Black")
                        xyouts, xouts(7)+0.01, [0.27], 'Sign.', /normal, charthick=2, charsize=0.5, color=cgcolor("Black")
                                ;xyouts, xouts(7), [0.6], 'IR', /normal, charthick=3, charsize=0.5, color=cgcolor("Black")

       
                        xyouts, xouts(0), [0.24], '['+Zsun+']', /normal, charthick=2, charsize=0.5, color=cgcolor("Black")
                        xyouts, xouts(1), [0.24], '[Myr]', /normal, charthick=2, charsize=0.5, color=cgcolor("Black")
                        xyouts, xouts(2), [0.24], '[Gyr]', /normal, charthick=2, charsize=0.5, color=cgcolor("Black")
                        xyouts, xouts(3), [0.24], '[Gyr]', /normal, charthick=2, charsize=0.5, color=cgcolor("Black")
                        xyouts, xouts(4), [0.24], '[mag]', /normal, charthick=2, charsize=0.5, color=cgcolor("Black")
                        xyouts, xouts(5)+0.01, [0.24], '['+Msun+']', /normal, charthick=2, charsize=0.5, color=cgcolor("Black")
                        xyouts, xouts(6)+0.01, [0.24], '[!4l!3Jy]', /normal, charthick=2, charsize=0.5, color=cgcolor("Black")
       			xyouts, xouts(7)+0.01, [0.24], '%', /normal, charthick=2, charsize=0.5, color=cgcolor("Black")




                        for s=0,num_cl-1 do begin
                           str_met1=strcompress(string(10^clus_pars(s, 4),'(f3.1)'),/remove_all)        ; met
                           str_tau1=strcompress(string(10^clus_pars(s,0)*1.e-6,'(i6)'),/remove_all)   ; tau
                           str_age=strcompress(string(10^clus_pars(s,1)*1.e-9,'(f5.1)'),/remove_all)  ;  age
                           str_age1=strcompress(string(10^clus_pars(s,2)*1.e-9,'(f5.1)'),/remove_all) ; mass weighted age
                           str_av1=strcompress(string(clus_pars(s,3),'(f5.1)'),/remove_all)           ; extinction
                           str_mass=strcompress(string(mass_cl(s),'(f5.2)'),/remove_all)              ; mass
                           str_ratio=strcompress(string(round(ratio(s)),'(i)'),/remove_all)                  ; ratio model
                           str_f24=strcompress(string(f24_sol(s),'(f5.1)'),/remove_all)     ; 24 um flux
                           ;if ir_det_sol(s) eq -1. then str_ir='NO'  
                           ;if ir_det_sol(s) eq 1. then  str_ir='YES'  ; IR detection

            
                           xyouts, xouts(0), youts(s), str_met1, /normal, charthick=2, charsize=0.5, color=colour(s)
                           xyouts, xouts(1), youts(s), str_tau1, /normal, charthick=2, charsize=0.5, color=colour(s)
                           xyouts, xouts(2), youts(s), str_age, /normal, charthick=2, charsize=0.5, color=colour(s)
                           xyouts, xouts(3), youts(s), str_age1, /normal, charthick=2, charsize=0.5, color=colour(s)
                           xyouts, xouts(4), youts(s), str_av1, /normal, charthick=2, charsize=0.5, color=colour(s)
                           xyouts, xouts(5), youts(s), str_mass, /normal, charthick=2, charsize=0.5, color=colour(s)
                           xyouts, xouts(6)+0.01, youts(s), str_f24, /normal, charthick=2, charsize=0.5, color=colour(s)
                           xyouts, xouts(7)+0.01, youts(s), str_ratio, /normal, charthick=2, charsize=0.5, color=colour(s)
            

endfor

                          ;======= simbols legend =============
                         
                       if id  eq 'SHARDS10003973' then begin

                        legend, ['Broad-Band data'], psym=sym(9), col=cgcolor("Black"), charthick=3, charsize=0.9,pos=[0.29,0.45], /norm, box=0
                        legend, ['SHARDS data'], psym=sym(6), col=cgcolor("Black"), charthick=3, charsize=0.9,pos=[0.29,0.42], /norm, box=0
                        legend, ['G102 data'], line=0,number=0.3, thick=3,col=cgcolor("Black"), charthick=3, charsize=0.9,pos=[0.29,0.39], /norm, box=0
                        legend, ['G141 data'], line=0,number=0.1, thick=3,col=cgcolor("Charcoal"), charthick=3, charsize=0.9,pos=[0.29,0.36], /norm, box=0
                        
                       endif
                        


			;----------axis---------------------------

			cgAxis, XAxis=1, YLog=1, XRange=xrangex/(1+z), /Save, ythick=3. ,chars=0.8, chart=3, ystyle=1, col=0, xthick=3., xstyle=9, xtitle=' '
			xyouts, [0.3], [0.97],'!4k!3!drest-frame!N [!4l!3m] ', /normal, charthick=3, charsize=1., color=cgcolor("Black")
	
			cgAxis, YAxis=1, YLog=0, YRange=[max_mag, min_mag], /Save, ytitle='AB mag', ythick=4. ,chars=1., chart=3, ystyle=1, col=0 ; lo ponemos al final para no cambiar escala


                       
                       ;==================== zoom Mg D4000 ================

				xrangex=[0.2, 0.45]     ; en RF!!
                                ymin=0.8*min(flux_cat_shards)
                                ymax=1.5*max(flux_cat_shards)
                       


                        plot, lambda_cat_shards_RF*1.e-4, flux_cat_shards,xrange=xrangex,yrange=[ymin, ymax], $
			xtitle='!4k!3!drest-frame!N [!4l!3m]', ytitle='F!d!7m!5!N [!7l!3Jy]',  thick=6, charthick=2., $
			xthick=3, ythick=3, chars=.7, /nodata,  /noerase, /xstyle, /ystyle,/ylog , position=p2

                        ;bandas indices

                        yfill1=ymin+0.002
			yfill2=ymax-0.1
                        
                         ;-----D4000--------

                         xfill1=[3750,4050]
		         xfill2=[3950,4250]


			colorbands=[cgcolor('RYB6'),cgcolor('Rose'), cgcolor('Light Cyan')]
			for f=0,n_elements(xfill1)-1 do begin &$
				cgColorFill, [xfill1(f), xfill1(f), xfill2(f), xfill2(f), xfill1(f)]*1.e-4, [yfill1, yfill2, yfill2, yfill1, yfill1],  COLOR=colorbands(f) &$
			endfor 
                         

   
                        ;------MgUV------ 

                                lm1=2525   		 
				lm2=2625
				lm3=2725
				lm4=2825
                          xfill1=[lm1,lm2,lm3]
			  xfill2=[lm2,lm3,lm4]


			for f=0,n_elements(xfill1)-1 do begin &$
				cgColorFill, [xfill1(f), xfill1(f), xfill2(f), xfill2(f), xfill1(f)]*1.e-4, [yfill1, yfill2, yfill2, yfill1, yfill1], COLOR=colorbands(f) &$
			endfor 
                      
   
			;Emission lines
			;------------------------------------
			oplot, [6562.,6562.]*1.e-4,[ymin,ymax], line=2, col=13   ; Halpha
			oplot, [4861.,4861.]*1.e-4,[ymin,ymax], line=2, col=13  ; Hbeta
			oplot, [3726.,3726.]*1.e-4,[ymin,ymax], line=2, col=13   ; OII




			; plotting different models
			;----------------------------
			for p=1,num_cl do oplot, clus_fnu(0,*)/(1+z)*1.e-4, clus_fnu(p,*), col=colour(p-1), thick=1
			oplot, clus_fnu(0,*)/(1+z)*1.e-4, clus_fnu(1,*), col=colour(0), thick=1
			

                        ;----------Broad Band---------
			oploterror, lambda_cat_bb_RF*1.e-4, flux_cat_bb,err_cat_bb, thick=4, errcol=cgcolor("Black"), psym=3, symsize=0.5, /nohat
			oplot, lambda_cat_bb_RF*1.e-4, flux_cat_bb, thick=2, col=cgcolor("White"), psym=sym(4), symsize=.9
			oplot, lambda_cat_bb_RF*1.e-4, flux_cat_bb, thick=4, col=cgcolor("Black"), psym=sym(9), symsize=.9
                        

                        ;----------SHARDS---------
			oploterror, lambda_cat_shards_RF*1.e-4, flux_cat_shards,err_cat_shards, thick=4, errcol=cgcolor("Black"), psym=3, symsize=0.5, /nohat
			oplot, lambda_cat_shards_RF*1.e-4, flux_cat_shards, thick=2, col=cgcolor("White"), psym=sym(1), symsize=.5
			oplot, lambda_cat_shards_RF*1.e-4, flux_cat_shards, thick=8, col=cgcolor("Black"), psym=sym(6), symsize=.5
                        
                        
                         ;----------Index-----------------
			
		str_d4=strcompress(string(d4_corr(i),'(f6.2)'),/remove_all) 
                str_mg=strcompress(string(mg_corr(i),'(f6.2)'),/remove_all) 
		str_d4_err=strcompress(string(d4_err(i),'(f6.2)'),/remove_all) 
		str_mg_err=strcompress(string(mg_err(i),'(f6.2)'),/remove_all) 
	
		symb=plusminus()

                xyouts, [0.67],[0.92],'D4000='+str_d4+''+symb+''+str_d4_err+'', /normal, charthick=1.5, charsize=0.6
		xyouts, [0.67],[0.90],'Mg!dUV!N='+str_mg+''+symb+' '+str_mg_err+'', /normal, charthick=1.5, charsize=0.6

		xyouts, [0.8],[0.92],'SHARDS region', /normal, charthick=2.5, charsize=0.8
               
	
       
                        ; reploteamos ejes por zonas indices

			plot, lambda_cat_shards_RF*1.e-4, flux_cat_shards,xrange=xrangex,yrange=[ymin, ymax], $
			xtitle='!4k!3!drest-frame!N [!4l!3m]', ytitle='',  thick=6, charthick=2., $
			xthick=3, ythick=3, chars=.7, /nodata,  /noerase, /xstyle, /ystyle,/ylog , position=p2

                        ;================ zoom en G102  ==========================
;stop

                 if n_elements(g102_b) gt 1 then begin

                       

			yg_max=max(flux_g102_high)
			yg_min=min(flux_g102_low)

                        ok_g102=where(lambda_cat_g102*1.e-4 gt 0.9 and lambda_cat_g102*1.e-4 lt 1.21)
			ok_g141=where(lambda_cat_g141*1.e-4 gt 0.9 and lambda_cat_g141*1.e-4 lt 1.21)

			plot, lambda_synth*1.e-4, flux_temp_synth, xrange=[0.9 ,1.21]/(1+z), yrange=[yg_min*0.80, yg_max*1.2], /xstyle, /ystyle,thick=2,xthick=3,ythick=3.,chars=0.7,chart=2 , xtitle='!4k!3!drest-frame!N [!4l!3m]', ytitle='F!d!7m!5!N [!7l!3Jy]', /nodata, position=p3, /noerase, xtickinterval=0.05
                        ;title=' '+spectra_name(w(j))+' '
			;position=[0.6,0.20,0.83,0.45]
                        ;position=[0.2,0.6,0.4,0.80]

			cgColorFill, [lambda_cat_g102_rf*1.e-4, reverse(lambda_cat_g102_rf)*1.e-4, lambda_cat_g102_rf(0)*1.e-4],[flux_g102_high, reverse(flux_g102_low), flux_g102_high(0)], color=cgcolor("Medium Grey")

	
			;Emission lines
			;------------------------------------
			oplot, [6562.,6562.]*1.e-4,[yg_min*0.80, yg_max*1.2], line=2, col=13   ; Halpha
			oplot, [4861.,4861.]*1.e-4,[yg_min*0.80, yg_max*1.2], line=2, col=13   ; Hbeta
			oplot, [3726.,3726.]*1.e-4,[yg_min*0.80, yg_max*1.2], line=2, col=13   ; OII


			

			; plotting different models
			;----------------------------
			for p=1,num_cl do oplot, clus_fnu(0,*)*1.e-4/(1+z), clus_fnu(p,*), col=colour(p-1), thick=1
			oplot, clus_fnu(0,*)*1.e-4/(1+z), clus_fnu(1,*), col=colour(0), thick=1
			

			;oplot, lambda_cat_g102_rf(ok_g102)*1.e-4, flux_cat_g102(ok_g102), thick=3, col=cgcolor("Lime Green"), psym=10
			;if n_elements(ok_g141) gt 1 then oplot, lambda_cat_g141_rf(ok_g141)*1.e-4, flux_cat_g141(ok_g141), thick=3, col=cgcolor("Lime Green"), psym=10
                        
                        oplot, lambda_cat_g102_rf(ok_g102)*1.e-4, flux_cat_g102(ok_g102), thick=4, col=cgcolor("Black"), psym=10
			if n_elements(ok_g141) gt 1 then oplot, lambda_cat_g141_rf(ok_g141)*1.e-4, flux_cat_g141(ok_g141), thick=4, col=cgcolor("Charcoal"), psym=10


			;xyouts,chars=0.7,chart=3, 1.02/(1+z), yg_max+1., 'class='+strcompress(string(fix(class(w(j)))),/remove_all)+'
			xyouts,chars=0.5,chart=1.5, 0.91/(1+z), yg_max, '<S/N> ='+strcompress(string(round(SN_g102(i)),'(i)'),/remove_all)+'
			legend, ['3!4r!3 errors'], psym=sym(5), symsize=1., col=cgcolor("Medium Grey"), charthick=1.5, charsize=0.5,pos=[0.65,0.64], /norm, box=0

			xyouts, [0.8],[0.6],'G102 region', /normal, charthick=2.5, charsize=0.8
			
	       endif   ;(end no G102)
			;-----------------------------------------------

                       

			
                        ;================ zoom en G141  ==========================

              if n_elements(g141_b) gt 1 then begin                  
                    

			yg_max=max(flux_g141_high)
			yg_min=min(flux_g141_low)

                        if g102_b(0) ne -1 then ok_g102=where(lambda_cat_g102*1.e-4 gt 1.1 and lambda_cat_g102*1.e-4 lt 1.7)
			if g141_b(0) ne -1 then ok_g141=where(lambda_cat_g141*1.e-4 gt 1.1 and lambda_cat_g141*1.e-4 lt 1.7)

			plot, lambda_synth*1.e-4, flux_temp_synth, xrange=[1.1 ,1.7]/(1+z), yrange=[yg_min*0.8, yg_max*1.2], /xstyle, /ystyle,thick=2,xthick=3,ythick=3.,chars=0.7,chart=2 , xtitle='!4k!3!drest-frame!N [!4l!3m]', ytitle='F!d!7m!5!N [!7l!3Jy]', /nodata, position=p4, /noerase
                        ;title=' '+spectra_name(w(j))+' '



			cgColorFill, [lambda_cat_g141_rf*1.e-4, reverse(lambda_cat_g141_rf)*1.e-4, lambda_cat_g141_rf(0)*1.e-4],[flux_g141_high, reverse(flux_g141_low), flux_g141_high(0)], color=cgcolor("Medium Grey")
  
			
			;Emission lines
			;------------------------------------
			oplot, [6562.,6562.]*1.e-4,[yg_min*0.80, yg_max*1.2], line=2, col=13   ; Halpha
			oplot, [4861.,4861.]*1.e-4,[yg_min*0.80, yg_max*1.2], line=2, col=13    ; Hbeta
			oplot, [3726.,3726.]*1.e-4,[yg_min*0.80, yg_max*1.2], line=2, col=13  ; OII





			; plotting different models
			;----------------------------
			for p=1,num_cl do oplot, clus_fnu(0,*)*1.e-4/(1+z), clus_fnu(p,*), col=colour(p-1), thick=1
			oplot, clus_fnu(0,*)*1.e-4/(1+z), clus_fnu(1,*), col=colour(0), thick=1
			
		        ; if n_elements(ok_g102) gt 1 then oplot, lambda_cat_g102_rf(ok_g102)*1.e-4, flux_cat_g102(ok_g102), thick=3, col=cgcolor("Lime Green"), psym=10
			;oplot, lambda_cat_g141_rf(ok_g141)*1.e-4, flux_cat_g141(ok_g141), thick=3, col=cgcolor("Lime Green"), psym=10


		         if n_elements(ok_g102) gt 1 then oplot, lambda_cat_g102_rf(ok_g102)*1.e-4, flux_cat_g102(ok_g102), thick=4, col=cgcolor("Black"), psym=10
			oplot, lambda_cat_g141_rf(ok_g141)*1.e-4, flux_cat_g141(ok_g141), thick=4, col=cgcolor("Charcoal"), psym=10


			;xyouts,chars=0.7,chart=3, 1.02/(1+z), yg_max+1., 'class='+strcompress(string(fix(class(w(j)))),/remove_all)+'
			xyouts,chars=0.5,chart=1.5, 1.12/(1+z), yg_max, '<S/N> ='+strcompress(string(round(SN_g141(i)),'(i)'),/remove_all)+'
			legend, ['3!4r!3 errors'], psym=sym(5), symsize=1., col=cgcolor("Medium Grey"), charthick=1.5, charsize=0.5,pos=[0.65,0.32], /norm, box=0
			xyouts, [0.8],[0.28],'G141 region', /normal, charthick=2.5, charsize=0.8

		endif   ;(end no G141)
			;-----------------------------------------------

                ; endif   ; selected galaxies

endfor




device,/close 
set_plot,'x'


stop


END






