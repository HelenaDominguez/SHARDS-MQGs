PRO  SED_index_tracks_paper

;  .comp /home/hdominguez/idl/pro/new/SED_index_tracks_paper.pro

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

analyze_dir='analyze_rejoin5'
index_dir='analyze_rejoin5'
;index_dir='analyze_index_tunned'


Msun='M'+sunsymbol()

filelist=dir+'/ID_SHARDS_candels_zphotoiris_spectra_files_104.lst'
readcol, filelist, source, zcomb, G141_file, G102_file, format='a, f, a, a'

filelist=dir+'/ID_ra_dec_zspec_104_sql.lst'
readcol, filelist, source2, ra, dec, zspec, zspec_flag, zcomb2, format='a, f, f, f'

readcol, dir+'/synthesizer_offset/files/analyze_parameters_BC_txp_krou.lst',   shards_id , nclus , z ,  log_tau,  log_age, log_age_MW,  av , met, mass, format='a, f, f, f'

age_gy=10^log_age*1.e-9
age_mw_gy=10^log_age_mw*1.e-9
tau_my=10^log_tau*1.e-6

;D4000 & MgUV values
restore, dir+'/synthesizer_offset/files/index_BC03_txp_krou_104_d4_err.sav', /verbose

;F24 lim
readcol, dir+'/synthesizer_offset/files/F24_output_BC_txp_krou_met_135.lst', shards_id_24, i_24, nclus_24, z_24 ,Mass_24,AV_24,F24, IR_det, format='a, i, i, f'




print, 'Version different clusters, begining PLOT!  for '+code+' ' +model+'  '+imf+''+met_string+''

 delta_wl1=fltarr(n_elements(source))
 delta_wl2=fltarr(n_elements(source))
 delta_wl=fltarr(n_elements(source))

N_clus=fltarr(n_elements(source))
Nbands=fltarr(n_elements(source))
id_list=strarr(n_elements(source))

stop


;i=92
;for i=0, n_elements(source)-1 do begin
for i=28, 28 do begin
	id=source(i)
	z=zcomb(i)

	print, i
        print, i
        print, i
        print, i

	print, id, z

       

       

		;========solutions=======================

		w=where(shards_id eq source(i))

		age_sol=age_gy(w)
		age_mw_sol=age_mw_gy(w)
		tau_sol=tau_my(w)
		mass_sol=mass(w)
		av_sol=av(w)
		met_sol=met(w)
                f24_sol=f24(w)
                ir_det_sol=ir_det(w)

		;ratio=ngal(w)*100/300

		wage=fltarr(n_elements(w))
		wtau=fltarr(n_elements(w))
		wav=fltarr(n_elements(w))


               ;if n_elements(w) gt 1 and min(mass_sol) ge 10.0 then begin   ; plot for degenerate galaxies only
               ;if n_elements(w) eq 1 and min(mass_sol(0)) ge 10.0 and age_mw_sol(0) gt 1.0 and age_mw_sol(0) lt 2.0 then begin   ; plot for degenerate galaxies only
       	 	
               id_list(i)=id

		; ----------
		; Photometry
		; ----------
		
		if code eq 'BC' then file=dir+'/synthesizer_offset/res/'+id+'_1pop.BCall.bc2003.stelib.csq.'+model+'.hr.'+imf+'_M_0.1_100.calz.res'
                if code eq 'M05' then  file=dir+'/synthesizer/res'+met_string+'/'+id+'_1pop.M05.none.none.csq.exp.hr.krou_M_0.1_100.calz.res'
		
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
                wl_pe=fltarr(nfiltphot)

		ind1=where(dataphot eq 'flux:')
		ind2=where(dataphot eq 'err:')
		ind3=where(dataphot eq 'filt:')
                ind4=where(dataphot eq 'wl:')

		for j=0, nfiltphot-1 do begin &$		
			flux_pe(j)=dataphot(j+ind1+1) &$
			err_pe(j)=dataphot(j+ind2+1) &$
			filt_pe(j)=dataphot(j+ind3+1) &$
			wl_pe(j)=dataphot(j+ind4+1) &$
		endfor

		;stop

                filt_used=where(ind_temp eq 1 and flux_pe gt 0 and err_pe gt 0)

		width_cat=width(filt_used)
		filt_cat=filtname(filt_used)		
		lambda_cat=Eff_wl(filt_used)
		lambda_cat_res=wl_pe(filt_used)
		lambda_cat_rf=Eff_wl_rf(filt_used)
		flux_cat=flux_pe(filt_used)
		err_cat=err_pe(filt_used)
                Nbands(i)=n_elements(flux_cat)

	       delta_wl(i)=max(abs(eff_wl(filt_used)*1.e-1-wl_pe(filt_used)))

              err_cat_original=err_cat
              err_cat=3*err_cat_original

              print, 'Nbands synth=', n_elements(flux_cat)
              print, 'maximal wl difference = ',max(abs(eff_wl(filt_used)*1.e-1-wl_pe(filt_used)))
              
		

;endfor


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
	
	        a=where(wl_cl_rf*(1+z) gt 5000. and wl_cl_rf*(1+z) lt 10000.)      ; whole range
		factor=median(flux_temp_synth(a)/fnu_muJy(a))     ; flux conversion factor from the data
		
		clus_fnu(0,*) = wl_cl_rf*(1+z)   ;observed
		clus_fnu(1,*) = fnu_muJy*factor

                ; Cluster parameters
                ;---------------------- 
		



                nlines=7

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
                 
		;stop

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
		 g102_SN=flux_cat_g102/err_cat_g102
		 SN_g102=median(g102_SN)


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
		g141_SN=flux_cat_g141/err_cat_g141
		SN_g141=median(g141_SN)


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
                 			

               ;============ index filters ====================

;stop

              print, mg_corr(i)
              print, d4_corr(i)

              if finite(mg_corr(i)) eq 1 then begin
                 file_mg=dir+'/synthesizer_offset/res/'+index_dir+'/'+id+'_mguvdaddi_'+id+'.res'
                 readcol, file_mg,x, band_mg, wl_mg, a,a,a,a,a,a,a, format='a, a, f, f,f,f,f,f,f', skipline=1
                 
                 blue_mg=where(band_mg eq 'BLUE01' or band_mg eq 'BLUE02' )
		 red_mg=where(band_mg eq 'RED01')

                 mg_b=fltarr(n_elements(wl_mg))
                 for k=0, n_elements(wl_mg)-1 do mg_b(k)=where(lambda_cat*1.e-1 eq wl_mg(k))
                 for k=0, n_elements(wl_mg)-1 do mg_b(k)=closest(lambda_cat*1.e-1, wl_mg(k))
                delta_wl1(i)=max(wl_mg-lambda_cat(mg_b)*1.e-1)
              endif

              if finite(d4_corr(i)) eq 1 then begin
                 file_d4=dir+'/synthesizer_offset/res/'+index_dir+'/'+id+'_d4000_'+id+'.res'
                 readcol, file_d4,x, band_d4, wl_d4, a,a,a,a,a,a,a, format='a, a, f, f,f,f,f,f,f', skipline=1

                 blue_d4=where(band_d4 eq 'BLUE01' or band_d4 eq 'BLUE02' )
		 red_d4=where(band_d4 eq 'RED01')            

                 d4_b=fltarr(n_elements(wl_d4))
                for k=0, n_elements(wl_d4)-1 do d4_b(k)=where(lambda_cat*1.e-1 eq wl_d4(k))
                for k=0, n_elements(wl_d4)-1 do d4_b(k)=closest(lambda_cat*1.e-1, wl_d4(k))
                delta_wl2(i)=max(wl_d4-lambda_cat(d4_b)*1.e-1)
             endif 
                

stop

                        ;==============================================
			;
			;               OUPUT PLOT
			; 
			;==============================================
			;==============================================



psfile=strtrim(dir+"/synthesizer_offset/files/SED_plot_index_tracks_"+code+"_"+model+"_"+imf+""+met_string+"_20000663_zoom.ps",2)


set_plot, 'ps'
device, filename=psfile, /portrait, /times, /color, xsize=9, ysize=3, /inch
!p.multi=0
		
                 ; horizontal           
                 p1=[0.1,0.15,0.95,0.93] 
                 ;p2=[0.15,0.15,0.9,0.45]


                ;vertical 
                ; p1=[0.15,0.15,0.5,0.9]
                ; p2=[0.58,0.15,0.9,0.9]

	colour=[cgColor("Red"),cgColor("Royal Blue"),cgColor("Orange"),cgColor("Turquoise"),cgColor("grn5"),cgColor("Blue Violet"),cgColor("Yellow"),cgColor("Violet"),cgColor("Grey")]


                       
                       ;==================== zoom Mg D4000 ================

				xrangex=[0.2, 0.45]     ; en RF!!
                                ymin=0.8*min(flux_cat_shards)
                                ymax=1.5*max(flux_cat_shards)
                       


                        plot, lambda_cat_shards_RF*1.e-4, flux_cat_shards,xrange=xrangex,yrange=[ymin, ymax], $
			xtitle='!4k!3!drest-frame!N [!4l!3m]', ytitle='F!d!7m!5!N [!7l!3Jy]',  thick=6, charthick=3., $
			xthick=3, ythick=3, chars=1., /nodata, /xstyle, ystyle=9,/ylog , position=p1

                        ;bandas indices

                        yfill1=ymin+0.002
			yfill2=ymax-0.2
                        
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
			oplot, [6562.,6562.]*1.e-4,[ymin,ymax], line=1, col=13   ; Halpha
			oplot, [4861.,4861.]*1.e-4,[ymin,ymax], line=1, col=13   ; Hbeta
			oplot, [3726.,3726.]*1.e-4,[ymin,ymax], line=1, col=13   ; OII




			; plotting different models
			;----------------------------
			for p=1,num_cl do oplot, clus_fnu(0,*)/(1+z)*1.e-4, clus_fnu(p,*), col=colour(p-1), thick=5
			oplot, clus_fnu(0,*)/(1+z)*1.e-4, clus_fnu(1,*), col=colour(0), thick=4
			
                        ;----------Broad Band---------
			oploterror, lambda_cat_bb_RF*1.e-4, flux_cat_bb,err_cat_bb, thick=4, errcol=cgcolor("Black"), psym=3, symsize=0.5, /nohat
			oplot, lambda_cat_bb_RF*1.e-4, flux_cat_bb, thick=2, col=cgcolor("White"), psym=sym(4), symsize=1.5
			oplot, lambda_cat_bb_RF*1.e-4, flux_cat_bb, thick=4, col=cgcolor("Black"), psym=sym(9), symsize=1.5
                        

                        ;----------SHARDS---------
			oploterror, lambda_cat_shards_RF*1.e-4, flux_cat_shards,err_cat_shards, thick=4, errcol=cgcolor("Black"), psym=3, symsize=0.5, /nohat
			oplot, lambda_cat_shards_RF*1.e-4, flux_cat_shards, thick=2, col=cgcolor("White"), psym=sym(1), symsize=1.5
			oplot, lambda_cat_shards_RF*1.e-4, flux_cat_shards, thick=8, col=cgcolor("Black"), psym=sym(6), symsize=1.5
                        

                        ;----------Index filters---------
			if finite(mg_corr(i)) eq 1 then oplot, lambda_cat_RF(mg_b(blue_mg))*1.e-4, flux_cat(mg_b(blue_mg)), thick=4, col=cgcolor("BLU6"), psym=sym(1), symsize=1.5
			if finite(mg_corr(i)) eq 1 then oplot, lambda_cat_RF(mg_b(red_mg))*1.e-4, flux_cat(mg_b(red_mg)), thick=4, col=cgcolor("RED6"), psym=sym(1), symsize=1.5
			if finite(mg_corr(i)) eq 1 then oplot, lambda_cat_RF(mg_b)*1.e-4, flux_cat(mg_b), thick=4, col=cgcolor("Black"), psym=sym(6), symsize=1.5
                        

		        if finite(d4_corr(i)) eq 1 then oplot, lambda_cat_RF(d4_b(blue_d4))*1.e-4, flux_cat(d4_b(blue_d4)), thick=4, col=cgcolor("BLU6"), psym=sym(1), symsize=1.5
			if finite(d4_corr(i)) eq 1 then oplot, lambda_cat_RF(d4_b(red_d4))*1.e-4, flux_cat(d4_b(red_d4)), thick=4, col=cgcolor("RED6"), psym=sym(1), symsize=1.5
			if finite(d4_corr(i)) eq 1 then oplot, lambda_cat_RF(d4_b)*1.e-4, flux_cat(d4_b), thick=4, col=cgcolor("Black"), psym=sym(6), symsize=1.5


                        
                         ;----------Index-----------------
			
		str_d4=strcompress(string(d4_corr(i),'(f6.2)'),/remove_all) 
		str_mg=strcompress(string(mg_corr(i),'(f6.2)'),/remove_all) 
		str_d4_err=strcompress(string(d4_err(i),'(f6.2)'),/remove_all) 
		str_mg_err=strcompress(string(mg_err(i),'(f6.2)'),/remove_all) 

                symb=plusminus()

                xyouts, [0.13],[0.65],'D4000='+str_d4+' '+symb+' '+str_d4_err+'', /normal, charsize=1., charthick=3
	        xyouts, [0.13],[0.6],'Mg!dUV!N='+str_mg+' '+symb+' '+str_mg_err+'', /normal, charsize=1., charthick=3	



                        ; labels
			;----------------------------
                if ir_det_sol(0) eq 1. then str_ir='F(24) > 30 !4l!3Jy (5!4r!3)'
                if ir_det_sol(0) eq -1. then str_ir='F(24) < 30 !4l!3Jy (5!4r!3)'
		
			ra_str=strcompress(string(ra(i),'(f7.3)'),/remove_all)
			dec_str=strcompress(string(dec(i),'(f7.4)'),/remove_all)


                        ;xyouts, [0.17], [0.87], 'SHARDS '+ra_str+'+'+dec_str+'', /normal, charthick=3, charsize=0.9, color=cgcolor("Black")
 			xyouts, [0.13], [0.87], 'SHARDSJ123620.3+6200844.3', /normal, charthick=3, charsize=0.9, color=cgcolor("Black") ;SHARDS20000663
			xyouts, [0.13], [0.82], 'z!dsp!N='+strcompress(string(z,'(f6.3)'),/remove_all), /normal, charthick=3, charsize=0.9, color=cgcolor("Black")
	         	xyouts, [0.13], [0.77], 'N bands='+strcompress(string(Nbands(i),'(i4)'),/remove_all), /normal, charthick=3, charsize=0.9, color=cgcolor("Black")
                        xyouts, [0.13], [0.72], str_ir, /normal, charthick=3, charsize=0.9, color=cgcolor("Black")

	            ;solution parameters
                    ;---------------------

		;ylegend=0.8	
		ylegend=0.5   ;horizontal
		youts=[0.40, 0.35, 0.30, 0.25, 0.2, 0.15, 0.1, 0.5]-0.05 ;horizontal plot
		;youts=[0.40, 0.35, 0.30, 0.25, 0.2, 0.15, 0.1, 0.5]+0.35
		;xouts=[ 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]
                xouts=[ 0.6, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88]

			xyouts, xouts(0)+0.01, [ylegend], 'Z', /normal, charthick=3, charsize=0.8, color=cgcolor("Black")
			xyouts, xouts(1)+0.01, [ylegend], '!4s!3', /normal, charthick=3, charsize=0.8, color=cgcolor("Black")
			xyouts, xouts(2)+0.01, [ylegend], 't!d0!N', /normal, charthick=3, charsize=0.8, color=cgcolor("Black")
			xyouts, xouts(3)+0.01, [ylegend], '!S!A-!R!Nt!dM!N', /normal, charthick=3, charsize=0.8, color=cgcolor("Black")
			xyouts, xouts(4)+0.01, [ylegend], 'A!dV!N', /normal, charthick=3, charsize=0.8, color=cgcolor("Black")
			xyouts, xouts(5), [ylegend], 'log(M)', /normal, charthick=3, charsize=0.8, color=cgcolor("Black")
			xyouts, xouts(6)+0.005, [ylegend], 'F!dpr!N(24)', /normal, charthick=3, charsize=0.7, color=cgcolor("Black")
			xyouts, xouts(7)+0.01, [ylegend], 'Sign.', /normal, charthick=3, charsize=0.8, color=cgcolor("Black")
                        
		       
			xyouts, xouts(0), [ylegend-0.05], '['+Zsun+']', /normal, charthick=3, charsize=0.8, color=cgcolor("Black")
			xyouts, xouts(1), [ylegend-0.05], '[Myr]', /normal, charthick=3, charsize=0.8, color=cgcolor("Black")
			xyouts, xouts(2), [ylegend-0.05], '[Gyr]', /normal, charthick=3, charsize=0.8, color=cgcolor("Black")
			xyouts, xouts(3), [ylegend-0.05], '[Gyr]', /normal, charthick=3, charsize=0.8, color=cgcolor("Black")
			xyouts, xouts(4), [ylegend-0.05], '[mag]', /normal, charthick=3, charsize=0.8, color=cgcolor("Black")
			xyouts, xouts(5)+0.01, [ylegend-0.05], '['+Msun+']', /normal, charthick=3, charsize=0.8, color=cgcolor("Black")
                        xyouts, xouts(6)+0.01, [ylegend-0.05], '[!4l!3Jy]', /normal, charthick=3, charsize=0.8, color=cgcolor("Black")
                        xyouts, xouts(7)+0.01, [ylegend-0.05], ' %', /normal, charthick=3, charsize=0.8, color=cgcolor("Black")

		        



        	w=where(shards_id eq source(i))

		age_sol=age_gy(w)
		age_mw_sol=age_mw_gy(w)
		tau_sol=tau_my(w)
		mass_sol=mass(w)
		av_sol=av(w)
		met_sol=met(w)
		;ratio=ngal(w)*100/300

             for l=0, n_elements(w) -1 do begin   ;cicle in sol
		   

		str_met1=strcompress(string(10^met_sol(l),'(f3.1)'),/remove_all)    ; met
		str_tau1=strcompress(string(tau_sol(l),'(i6)'),/remove_all)         ; tau
		str_age=strcompress(string(age_sol(l),'(f5.1)'),/remove_all)    ; mass weighted age
		str_age1=strcompress(string(age_mw_sol(l),'(f5.1)'),/remove_all)    ; mass weighted age
		str_av1=strcompress(string(av_sol(l),'(f6.2)'),/remove_all)         ; extinction
		str_mass=strcompress(string(mass_sol(l),'(f5.2)'),/remove_all)      ; mass
		str_ratio=strcompress(string(round(ratio(l)),'(i)'),/remove_all)           ; ratio model
                str_f24=strcompress(string(f24_sol(l),'(f6.1)'),/remove_all)            ; F24 from model

			    
			    xyouts, xouts(0), youts(l), str_met1, /normal, charthick=4, charsize=0.8, color=colour(l)
			    xyouts, xouts(1), youts(l), str_tau1, /normal, charthick=4, charsize=0.8, color=colour(l)
			    xyouts, xouts(2), youts(l), str_age, /normal, charthick=4, charsize=0.8, color=colour(l)
			    xyouts, xouts(3), youts(l), str_age1, /normal, charthick=4, charsize=0.8, color=colour(l)
			    xyouts, xouts(4), youts(l), str_av1, /normal, charthick=4, charsize=0.8, color=colour(l)
			    xyouts, xouts(5), youts(l), str_mass, /normal, charthick=4, charsize=0.8, color=colour(l)
			    xyouts, xouts(6)+0.01, youts(l), str_f24, /normal, charthick=4, charsize=0.8, color=colour(l)
			    xyouts, xouts(7)+0.01, youts(l), str_ratio, /normal, charthick=4, charsize=0.8, color=colour(l)



			 endfor   ;cicle sol


       ; lamnvda Rest frame

                        plot, lambda_cat_shards_RF*1.e-4, flux_cat_shards,xrange=xrangex,yrange=[ymin, ymax], $
			xtitle='', ytitle='',  thick=6, charthick=3., $
			xthick=3, ythick=3, chars=1., /nodata, /xstyle, ystyle=9,/ylog , position=p1, /noerase


	;-----AB mag------


			min_mag=-2.5*alog10(ymax)+23.9
			max_mag=-2.5*alog10(ymin)+23.9

	cgAxis, YAxis=1, YLog=0, YRange=[max_mag, min_mag], /Save, ytitle='AB mag', ythick=4. ,chars=1., chart=3, ystyle=1, col=0, ytickinterval=2. ; lo ponemos al final para no cambiar escala


device,/close 
set_plot,'x'


	;==========================PLOT index tracks==============================================
			


print, 'EMPIEZA plot tracks'
stop

psfile=strtrim(dir+"/synthesizer_offset/files/SED_plot_index_tracks_"+code+"_"+model+"_"+imf+""+met_string+"_20000663_tracks.ps",2)


set_plot, 'ps'
device, filename=psfile, /portrait, /times, /color, xsize=7, ysize=7, /inch
!p.multi=0
	

	
		;========solutions=======================

		w=where(shards_id eq source(i))

		age_sol=age_gy(w)
		age_mw_sol=age_mw_gy(w)
		tau_sol=tau_my(w)
		mass_sol=mass(w)
		av_sol=av(w)
		met_sol=met(w)
		;ratio=ngal(w)*100/300

		wage=fltarr(n_elements(w))
		wtau=fltarr(n_elements(w))
		wav=fltarr(n_elements(w))


		av_vect=[0., 1., 2.]
		av_str_vect=['0.00', '1.00', '2.00']


		thick_vect=[3., 2.5, 2., 1.5, 1., 1., 1.]
		size_vect=[3., 2.5, 2., 1.5, 1., 1., 1.]



		plot, d4_corr,mg_corr, xrange=[1.0, 2.6], yrange=[1.0, 2.2], /xstyle, /ystyle, xtitle='D4000', ytitle='Mg!dUV!N',xthick=4, ythick=4, charthick=3., charsize=1.5, /nodata, position=[0.11,0.1,0.95,0.95], /noerase



		   for l=0, n_elements(w) -1 do begin   ;cicle in sol
		    
		wav(l)=closest(av_vect, av_sol(l))
		av_str=av_str_vect(wav(l))

		;===================  File with index  =====================

		; index carmen
		;---------------------------------------------------------------
		index_file=dir+'/synthesizer_offset/index_lib/version_bc2003_sfh_txp_lib_stelib_res_hr_imf_krou_mlo_0.1_mup_100_met_m62_av_'+av_str+'.indices.dat'
	        data=leefile(index_file)



		; # Cols: 
		; #    1 - Time in log(yr)
		; #    2 - Measured index identifier
		; #    3 - tauold Gyr : 000.00100000
		; #    4 - tauold Gyr : 000.02000000
		; #    5 - tauold Gyr : 000.05000000
		; #    6 - tauold Gyr : 000.07500000
		; #    7 - tauold Gyr : 000.10000000
		; #    8 - tauold Gyr : 000.20000000
		; #    9 - tauold Gyr : 000.50000000
		; #   10 - tauold Gyr : 001.00000000
		; #   11 - tauold Gyr : 002.00000000
		; #   12 - tauold Gyr : 004.00000000
		; #   13 - tauold Gyr : 007.00000000
		; #   14 - tauold Gyr : 009.00000000
		; #   15 - tauold Gyr : 011.00000000
		; #   16 - tauold Gyr : 013.00000000
		; #   17 - tauold Gyr : 015.00000000
		; #   18 - tauold Gyr : 020.00000000


		index_name=data(1,*)
		ages=10^data(0,*)*1.e-9                                       ;edades modelos = filas
		taus=[-1, -1, 1.,20, 50, 75, 100., 200., 500., 1000., 2000.]  ;taus modelos   = columnas, -1. to skip name and age columns

		wmg=where(index_name eq 'MgUV')
		wd4=where(index_name eq 'D4000')
                wdn4=where(index_name eq 'Dn4000')

		wa0=closest(ages(wd4), 0.5)  ; row where ages are those (same for D4 and Mg if include wd4 or wmg)
		wa1=closest(ages(wd4), 0.8)
		wa2=closest(ages(wd4), 1.0)
		wa3=closest(ages(wd4), 1.2)
		wa4=closest(ages(wd4), 1.4)
		wa5=closest(ages(wd4), 1.8)
		wa6=closest(ages(wd4), 2.0)
		wa7=closest(ages(wd4), 2.5)
                wa8=closest(ages(wd4), 3.0)
                wa9=closest(ages(wd4), 4.0)
                wa10=closest(ages(wd4), 5.0)


		;vect_ages=[wa0, wa1, wa2, wa3, wa4, wa5,wa6, wa7, wa8, wa9, wa10]
        	vect_ages=[wa0,  wa2, wa6, wa8, wa10]
		vect_tau=[2, 3, 4, 5]
		size_vect_ages=[0.3,0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9,2.1, 2.2]*3.0




		   wage(l)=closest(ages(wd4), age_sol(l))  ; nearest model age
		   wtau(l)=closest(taus, tau_sol(l))       ;nearest tau age
		 
;stop

		;track
		oplot, data(wtau(l), wd4(1:n_elements(wd4)-1)), data(wtau(l),wmg(1:n_elements(wmg)-1)), thick=6., col=colour(l)

		;ages
		for n=0, n_elements(vect_ages)-1 do oplot, [data(wtau(l), wd4(vect_ages(n))),data(wtau(l), wd4(vect_ages(n)))], [data(wtau(l),wmg(vect_ages(n))),data(wtau(l),wmg(vect_ages(n)))], col=colour(l), psym=sym(6), symsize=size_vect_ages(n)

		;closest age
		oplot, [data(wtau(l), wd4(wage(l))),data(wtau(l), wd4(wage(l)))], [data(wtau(l),wmg(wage(l))),data(wtau(l),wmg(wage(l)))], col=colour(l), psym=sym(4), symsize=3.8

		;index form phot
		if finite(d4_corr(i)) eq 1 and  finite(mg_corr(i)) eq 1 then begin
		oplot, [d4_corr(i),d4_corr(i) ], [mg_corr(i),mg_corr(i)], psym=sym(1), symsize=2.
		oploterror, [d4_corr(i),d4_corr(i) ], [mg_corr(i),mg_corr(i)],  [d4_err(i),d4_err(i) ], [mg_err(i),mg_err(i)],psym=sym(1), symsize=2.
		endif

                  endfor   ;cicle sol

          

legend, [' Indices from SHARDS data'], psym=sym(1), symsize=size_vect_ages(4), box=0, pos=[0.15,0.9], charsize=1.5, charthick=4, /normal
legend, [' t!d0!N= 0.5 Gyr'], psym=sym(6), symsize=size_vect_ages(0), box=0, pos=[0.15,0.85], charsize=1.5, charthick=4, /normal
legend, [' t!d0!N= 1.0 Gyr'], psym=sym(6), symsize=size_vect_ages(1), box=0, pos=[0.15,0.8], charsize=1.5, charthick=4, /normal
legend, [' t!d0!N= 2.0 Gyr'], psym=sym(6), symsize=size_vect_ages(2), box=0, pos=[0.15,0.75], charsize=1.5, charthick=4, /normal
legend, [' t!d0!N= 3.0 Gyr'], psym=sym(6), symsize=size_vect_ages(3), box=0, pos=[0.15,0.7], charsize=1.5, charthick=4, /normal
legend, [' t!d0!N= 5.0 Gyr'], psym=sym(6), symsize=size_vect_ages(4), box=0, pos=[0.15,0.65], charsize=1.5, charthick=4, /normal



		endfor            ;cicle gal








;stop

device,/close 
set_plot,'x'

stop
stop



;file with id list
;--------------------
 openw,1, dir+'/synthesizer_offset/files/id_ages_12.lst', width=1000
printf, 1, 'shards_id     '
for r =0, 134 do  printf,1, id_list(r)
close, 1
                         


stop

END

