PRO combine_2pop_SED_paper_both

;   .comp /home/hdominguez/idl/pro/new/combine_2pop_SED_paper_both.pro

dir='/home/hdominguez/SHARDS/DR2_gamma/redanddead_gamma'
Msun='M'+sunsymbol()


burst='senior'
;burst='interm'


cat_input_phot=dir+'/synthesizer_offset/Stacks/synthesizer_stack_HDFN_1Gyr_or_less/data/observed/seds_stack_1Gyr_or_less.cat'


  dataphot=leefile(cat_input_phot)
	
		sz=size(dataphot)
		nfiltphot=dataphot(sz(1)-1)



		flux_y_cat=fltarr(nfiltphot)
		err_y_cat=fltarr(nfiltphot)
        	wl_y_cat=strarr(nfiltphot)

		ind1=3
		ind2=3+nfiltphot
		ind3=3+2*nfiltphot

		for j=0, nfiltphot-1 do begin &$		
			flux_y_cat(j)=float(dataphot(2*j+ind1)) &$
			err_y_cat(j)=float(dataphot(2*j+ind1+1)) &$
			wl_y_cat(j)=float(dataphot(j+ind3)) &$
		endfor

;flux is Fnu
;----------------------
ok=where(flux_y_cat ne max(flux_y_cat))  ; eliminamos punto outlier

flux_y=flux_y_cat(ok)
err_y=err_y_cat(ok)
wl_y=wl_y_cat(ok)*1.e-3  ;wl in um

large_err=where(err_y gt flux_y)

err_y(large_err)=0.98*flux_y(large_err)

err_rel=err_y/flux_y


;------Filters-------

filt_min=fltarr(n_elements(wl_y))
filt_max=fltarr(n_elements(wl_y))

for i=0, n_elements(ok)-2 do begin

bb=strcompress(string(ok(i+1)), /remove_all)

if i+1 lt 10 then filt_file=strtrim(dir+"/synthesizer_offset/Stacks/synthesizer_stack_HDFN_1Gyr_or_less/FILTERS/stack/filter_00"+bb+".res",1)
if i+1 ge 10 and i+1 lt 100 then filt_file=strtrim(dir+"/synthesizer_offset/Stacks/synthesizer_stack_HDFN_1Gyr_or_less/FILTERS/stack/filter_0"+bb+".res",1)
if i+1 ge 100 then filt_file=strtrim(dir+"/synthesizer_offset/Stacks/synthesizer_stack_HDFN_1Gyr_or_less/FILTERS/stack/filter_"+bb+".res",1)

readcol, filt_file, range, trans
aa=where(trans eq 1)

filt_min(i)=range(aa(0))*1.e-4
filt_max(i)=range(aa(n_elements(aa)-1))*1.e-4

endfor


;=====1 sigma error shadow=====

max_flux=flux_y*(1+err_rel)
min_flux=flux_y*(1-err_rel)

;----- make square vector -----

flux_up=fltarr(n_elements(max_flux)*3-3)
flux_low=fltarr(n_elements(max_flux)*3-3)
wl_err=fltarr(n_elements(max_flux)*3-3)

stop

for j=1, n_elements(wl_y)-2 do begin &$
print, j

   ;delta1=(wl_y(j)-wl_y(j-1))/2  &$
   ;delta2=(wl_y(j+1)-wl_y(j))/2  &$

;delta=max(delta1, delta2)

  ; wl_err(3*j)=wl_y(j)-delta   &$
  ; wl_err(3*j+1)=wl_y(j)   &$
  ; wl_err(3*j+2)=wl_y(j)+delta  &$
   

   wl_err(3*j)=filt_min(j)   &$
   wl_err(3*j+1)=wl_y(j)   &$
   wl_err(3*j+2)=filt_max(j)  &$
   
if filt_min(j) gt wl_y(j) then  print, 'ERRORR!!' 
if filt_max(j) lt wl_y(j) then print, 'ERRORR!!'

   flux_up(3*j+1)=max_flux(j)   &$
   flux_up(3*j)=max_flux(j)   &$
   flux_up(3*j+2)=max_flux(j)   &$

   flux_low(3*j+1)=min_flux(j)    &$
   flux_low(3*j)=min_flux(j)   &$
   flux_low(3*j+2)=min_flux(j)  &$

endfor

nn=n_elements(wl_err)

wl_err=wl_err(3:nn-1)
flux_up=flux_up(3:nn-1)
flux_low=flux_low(3:nn-1)

;-----models----------

if burst eq 'senior' then readcol, dir+'/synthesizer_offset/Stacks/template_bc2003_krou_Zsun_350Myr_3.0Gyr_Av=04.ANGSTROMS_SED', wl_mod_o_aa, L_nu_o, L_lambda 
if burst eq 'interm' then readcol, dir+'/synthesizer_offset/Stacks/template_bc2003_krou_Zsun_200Myr_1.4Gyr_Av=04.ANGSTROMS_SED', wl_mod_o_aa, L_nu_o, L_lambda 

readcol, dir+'/synthesizer_offset/Stacks/template_bc2003_krou_Zsun_50Myr_0.8Gyr_Av=08.ANGSTROMS_SED', wl_mod_y_aa, L_nu_y, L_lambda_y 



wl_mod_y=wl_mod_y_aa*1.e-4  ;um
wl_mod_o=wl_mod_o_aa*1.e-4  ;um



;-----2 pop models-----

if burst eq 'senior' then begin
L_nu_2popa=0.2*L_nu_o+0.8*L_Nu_y
L_nu_2popb=0.3*L_nu_o+0.7*L_Nu_y
L_nu_2popc=0.5*L_nu_o+0.5*L_Nu_y
L_nu_2popd=0.8*L_nu_o+0.2*L_Nu_y
endif 



if burst eq 'interm' then begin
L_nu_2popa=0.01*L_nu_o+0.99*L_Nu_y
L_nu_2popb=0.05*L_nu_o+0.95*L_Nu_y
L_nu_2popc=0.1*L_nu_o+0.9*L_Nu_y
L_nu_2popd=0.2*L_nu_o+0.8*L_Nu_y

endif 

;----Nomarlize 3.6 um with respect to young----

range_stack=where(wl_y gt 0.9 and  wl_y lt 3.)   ; 1.5*(1+z)=3.6, son lam_RF

range_mod_y=where(wl_mod_y gt 0.9 and  wl_mod_y lt 3.)
range_mod_o=where(wl_mod_o gt 0.9 and  wl_mod_o lt 3.)


Ky=flux_y(range_stack)/mean(L_nu_y(range_mod_y))
Ko=mean(L_nu_o(range_mod_o))/mean(L_nu_y(range_mod_y))
Ka=mean(L_nu_2popa(range_mod_o))/mean(L_nu_y(range_mod_y))
Kb=mean(L_nu_2popb(range_mod_o))/mean(L_nu_y(range_mod_y))
Kc=mean(L_nu_2popc(range_mod_o))/mean(L_nu_y(range_mod_y))
Kd=mean(L_nu_2popd(range_mod_o))/mean(L_nu_y(range_mod_y))



;=====R=50===============


region=where(wl_mod_y gt 0.1 and wl_mod_y lt 4.)

fres_y=fltarr(n_elements(region))
fres_2popa=fltarr(n_elements(region))
fres_2popb=fltarr(n_elements(region))
fres_2popc=fltarr(n_elements(region))
fres_2popd=fltarr(n_elements(region))

resolution='50'
if resolution eq '2500' then res_value=2500
if resolution eq '1000' then res_value=1000
if resolution eq '50' then res_value=50

print, res_value


for i=0, n_elements(region)-1 do begin &$

wl=wl_mod_y(region(i))

delta_lam=wl/res_value &$
wl_up=wl+delta_lam/2 &$
wl_d=wl-delta_lam/2 &$

wres2=closest(wl_mod_y, wl_up) &$
wres1=closest(wl_mod_y, wl_d) &$


fres_y(i)=mean(L_nu_y(wres1:wres2)) &$
fres_2popa(i)=mean(L_nu_2popa(wres1:wres2)) &$
fres_2popb(i)=mean(L_nu_2popb(wres1:wres2)) &$
fres_2popc(i)=mean(L_nu_2popc(wres1:wres2)) &$
fres_2popd(i)=mean(L_nu_2popd(wres1:wres2)) &$

endfor





;-----Grism data------------

delta_wl=fltarr(n_elements(wl_y)-1)
for i=0, n_elements(wl_y)-2 do delta_wl(i)=wl_y(i+1)-wl_y(i)

grism=where(delta_wl le 0.002)
no_grism=where(delta_wl gt 0.002)


grism=where(wl_y gt 0.5 and wl_y lt 0.9)
no_grism=where(wl_y lt 0.5 or wl_y gt 0.9)

;=========SFH mean=========

tt=findgen(1000)/100
t=10^tt                                 ; tiempo sin log, eje de las x en la integral

t2=t+3.e9


;tau
;-------
if burst eq 'senior' then val_t1=400
if burst eq 'interm' then val_t1=200
val_t2=50

tau1=val_t1*1.e6
tau2=val_t2*1.e6


;---SFR---
sfr1_aux=(2.71828)^(-t/tau1)*t   ;old
sfr2_aux=(2.71828)^(-t/tau2)*t   ;young

k1=integral(t, sfr1_aux)
k2=integral(t, sfr2_aux)


Mass_mature=10^10.5

if burst eq 'senior' then begin

sfr1a=(sfr1_aux*mass_mature*0.2)/k1   ;20% old, 80% young
sfr2a=(sfr2_aux*mass_mature*0.8)/k2

sfr1b=(sfr1_aux*mass_mature*0.3)/k1  ;30% old, 70% young
sfr2b=(sfr2_aux*mass_mature*0.7)/k2

sfr1c=(sfr1_aux*mass_mature*0.5)/k1  ;50% old, 50% young
sfr2c=(sfr2_aux*mass_mature*0.5)/k2

sfr1d=(sfr1_aux*mass_mature*0.8)/k1  ;20% old, 80% young
sfr2d=(sfr2_aux*mass_mature*0.2)/k2

endif 


if burst eq 'interm' then begin


sfr1a=(sfr1_aux*mass_mature*0.01)/k1   ;1% int, 99% young
sfr2a=(sfr2_aux*mass_mature*0.99)/k2

sfr1b=(sfr1_aux*mass_mature*0.05)/k1  ;5% int, 95% young
sfr2b=(sfr2_aux*mass_mature*0.95)/k2

sfr1c=(sfr1_aux*mass_mature*0.1)/k1  ;10% int, 90% young
sfr2c=(sfr2_aux*mass_mature*0.9)/k2

sfr1d=(sfr1_aux*mass_mature*0.2)/k1  ;20% int, 80% young
sfr2d=(sfr2_aux*mass_mature*0.8)/k2

endif


print, 'Empieza plot '+burst+''
stop
;====================PLOT=====================================


  set_plot,'PS'
  Device, /Helvetica, /Color, Bits_per_Pixel=8, File=dir+'/synthesizer_offset/files/combine_2pop_SED_paper_'+burst+'_prueba2.eps', XSize=15, YSize=22

p1=[0.11,0.45,0.95,0.95]
p2=[0.11,0.27,0.95,0.44]
p3a=[0.11,0.06,0.31,0.18]
p3b=[0.32,0.06,0.52,0.18]
p3c=[0.53,0.06,0.73,0.18]
p3d=[0.74,0.06,0.94,0.18]


xmin=0.18
xmax=3.75



if burst eq 'senior' then begin color_burst=cgcolor('RED6')
 yrange1=10^16.
 yrange2=10^19.5
endif 

if burst eq 'interm' then begin
color_burst=cgcolor('GRN8')
 yrange1=10^16.6
 yrange2=10^19.5 
endif 


print, color_burst, yrange1
stop

;------------Plot whole SED Young -----------------
;===================================================
plot, wl_y, flux_y/ky(0), /xlog, /ylog, psym=1,xtitle='',ytitle='L!d!7m!5!N(erg s!u-1!N Hz!u-1!N'+Msun+'!u-1!N)',xrange=[xmin, xmax],yrange=[yrange1,yrange2], /xstyle, /ystyle, thick=3, charthick=2.,charsize=1., xthick=3., ythick=3.,xtickformat='(A1)', /nodata, position=p1

xyouts, 0.5-0.01,0.0075,'0.5', charthick=2.,charsize=.75
xyouts, 2.0-0.01,0.0075,'2', charthick=2.,charsize=.75
xyouts, 3.0-0.01,0.0075,'3', charthick=2.,charsize=.75


out=where(flux_low/Ky(0) lt yrange1)
flux_low(out)=yrange1*Ky(0)


;stop

ind2=where(wl_err gt 0. and wl_err le 0.83)


for l=0, n_elements(wl_y)-3 do begin &$

j=3*l+1 &$

cgColorFill, [wl_err(j-1), wl_err(j-1:j+1), reverse(wl_err(j-1:j+1)), wl_err(j-1)], [flux_up(j-1) , flux_up(j-1:j+1), reverse(flux_low(j-1:j+1)),flux_low(j-1)]/ky(0),  COLOR=cgcolor('Light grey') &$

endfor

;cgColorFill, [wl_err(ind2(0)), wl_err(ind2), reverse(wl_err(ind2)), wl_err((ind2(0)))], [flux_up((ind2(0))) , flux_up(ind2), reverse(flux_low(ind2)),flux_low((ind2(0)))]/ky(0),  COLOR=cgcolor('Light grey')

;ind1=where(wl_y ge 0.83)

;cgColorFill, [wl_y(ind1(0)), wl_y(ind1), reverse(wl_y(ind1)), wl_y(ind1)], [max_flux(ind1(0)) , max_flux(ind1), reverse(min_flux(ind1)),min_flux(ind1(0))]/ky(0),  COLOR=cgcolor('light grey')


oplot, wl_mod_y, L_nu_y, col=cgcolor('Blu7'), Thick=2.  ;young SFH
oplot, wl_mod_y, L_nu_o, col=color_burst, Thick=2.  ;old SFH
if burst eq 'senior' then oplot, wl_mod_o, L_nu_2popa/Ka(0), col=cgcolor('Blu4')
if burst eq 'interm' then oplot, wl_mod_o, L_nu_2popb/Ka(0), col=cgcolor('TG5')
oplot, wl_mod_o, L_nu_2popc/kc(0), col=cgcolor('RYB2')


oplot,wl_y(no_grism), flux_y(no_grism)/ky(0) , thick=2, col=cgcolor("White"), psym=sym(1), symsize=.5
oplot, wl_y(no_grism), flux_y(no_grism)/ky(0), thick=0.5, col=cgcolor("Black"), psym=sym(6), symsize=.5
oplot,wl_y(grism), flux_y(grism)/Ky(0) , thick=2, col=cgcolor("black"), psym=10   



xyouts, 0.45,0.8, /norm, 'Mature: t!d0!N=0.8 Gyr, !4s!3=50 Myr, A!dV!N=0.8 mag', col=cgcolor('Blu7'), charthick=2.5, charsize=.8
if burst eq 'senior' then xyouts, 0.45,0.78, /norm,'Senior: t!d0!N=3.0 Gyr, !4s!3=400 Myr, A!dV!N=0.4 mag', col=cgcolor('Red6'), charthick=2.5, charsize=.8
if burst eq 'interm' then xyouts, 0.45,0.78, /norm, 'Interm.: t!d0!N=1.4 Gyr, !4s!3=200 Myr, A!dV!N=0.5 mag', col=cgcolor('GRN8'), charthick=2.5, charsize=.8

if burst eq 'senior' then begin
xyouts, 0.55,0.75, /norm,'80% Mature + 20% Senior', col=cgcolor('Blu4'), charthick=2.5, charsize=.8
xyouts, 0.55,0.73, /norm,'50% Mature + 50% Senior', col=cgcolor('RYB2'), charthick=2.5, charsize=.8
endif


if burst eq 'interm' then begin
xyouts, 0.55,0.75, /norm,'95% Mature + 5% Interm.', col=cgcolor('TG5'), charthick=2.5, charsize=.8
xyouts, 0.55,0.73, /norm,'90% Mature + 10% Interm.', col=cgcolor('RYB2'), charthick=2.5, charsize=.8
endif




plot, wl_y, flux_y/ky(0), /xlog, /ylog, psym=1,xtitle='',ytitle='',xrange=[xmin, xmax],yrange=[yrange1, yrange2], /xstyle, /ystyle, thick=3, charthick=2.,charsize=1., xthick=3., ythick=3.,xtickformat='(A1)', /nodata, position=p1, /noerase


;------------Plot model range -----------------
;=============================================

aa=fltarr(n_elements(wl_mod_y))
aa[0:n_elements(wl_mod_y)-1]=1.

x1=xmin
x2=xmax

y1=1.08
y2=0.92

yy1=1.16
yy2=0.84


if burst eq 'interm' then begin
ymin=0.5
ymax=1.7
endif 

if burst eq 'senior' then begin
ymin=0.4 
ymax=1.4
endif


print, ymax


plot, wl_mod_y, L_nu_y/L_nu_2popa, /xlog,xtitle='!4k!3!drest-frame!N [!4l!3m]',ytitle='F!d2pop!N/F!dmature!N',xrange=[xmin,xmax],yrange=[ymin,ymax], /xstyle, /ystyle, thick=3, charthick=2.,charsize=.8, xthick=3., ythick=3.,yminor=1, yticklen=0.01,  /nodata, position=p2, /noerase



cgColorFill, [x1, x2, x2, x1], [yy1, yy1, yy2, yy2],  COLOR=cgcolor('light grey')
cgColorFill, [x1, x2, x2, x1], [y1, y1, y2, y2],  COLOR=13


;oplot, wl_mod_y, L_nu_2popa/(Ka(0)*L_nu_y), col=cgcolor('Blu3')
;oplot, wl_mod_y, L_nu_2popb/(Kb(0)*L_nu_y), col=cgcolor('GRN2')
;oplot, wl_mod_y, L_nu_2popc/(Kc(0)*L_nu_y), col=cgcolor('RYB3')
;oplot, wl_mod_y, L_nu_2popd/(Kd(0)*L_nu_y), col=cgcolor('Thistle')

oplot, wl_mod_y(region),fres_2popa/(Ka(0)*fres_y), col=cgcolor('Blu4'), Thick=4.
oplot, wl_mod_y(region),fres_2popb/(Kb(0)*fres_y), col=cgcolor('TG5'), Thick=4.
oplot, wl_mod_y(region),fres_2popc/(Kc(0)*fres_y), col=cgcolor('RYB2'), Thick=4.
oplot, wl_mod_y(region),fres_2popd/(Kd(0)*fres_y), col=cgcolor('Maroon'), Thick=4.


oplot,  wl_mod_y, aa, line=1.

if burst eq 'senior' then begin
xyouts, 0.6,0.34, /norm,'a) 80% Mature + 20% Senior', col=cgcolor('Blu4'), charthick=2.5, charsize=.8
xyouts, 0.6,0.32, /norm,'b) 70% Mature + 30% Senior', col=cgcolor('TG5'), charthick=2.5, charsize=.8
xyouts, 0.6,0.30, /norm,'c) 50% Mature + 50% Senior', col=cgcolor('RYB2'), charthick=2.5, charsize=.8
xyouts, 0.6,0.28, /norm,'d) 20% Mature + 80% Senior', col=cgcolor('Maroon'), charthick=2.5, charsize=.8

legend, ['1!4r!3 stack error'], psym=sym(5), box=0, pos=[0.6,0.44], charsize=.8, charthick=2, /normal, color=13
legend, ['2!4r!3 stack error'], psym=sym(5), box=0, pos=[0.6,0.42], charsize=.8, charthick=2, /normal, color=cgcolor('light grey')

endif

if burst eq 'interm' then begin
xyouts, 0.6,0.425, /norm,'a) 99% Mature + 1% Interm.', col=cgcolor('Blu4'), charthick=2.5, charsize=.8
xyouts, 0.6,0.405, /norm,'b) 95% Mature + 5% Interm', col=cgcolor('TG5'), charthick=2.5, charsize=.8
xyouts, 0.6,0.385, /norm,'c) 90% Mature + 10% Interm.', col=cgcolor('RYB2'), charthick=2.5, charsize=.8
xyouts, 0.6,0.365, /norm,'d) 80% Mature + 20% Interm.', col=cgcolor('Maroon'), charthick=2.5, charsize=.8

legend, ['1!4r!3 stack error'], psym=sym(5), box=0, pos=[0.6,0.32], charsize=.8, charthick=2, /normal, color=13
legend, ['2!4r!3 stack error'], psym=sym(5), box=0, pos=[0.6,0.30], charsize=.8, charthick=2, /normal, color=cgcolor('light grey')

endif



plot, wl_mod_y, L_nu_y/L_nu_2popa, /xlog,xtitle='',ytitle='',xrange=[xmin,xmax],yrange=[ymin,ymax], /xstyle, /ystyle, thick=3, charthick=2.,charsize=.8, xthick=3., ythick=3.,yminor=1, yticklen=0.01,   /nodata, position=p2, /noerase

xyouts, 0.2-0.01,ymin-0.1,'0.2', charthick=2.2,charsize=.8
xyouts, 0.3-0.01,ymin-0.1,'0.3', charthick=2.2,charsize=.8
xyouts, 0.4-0.01,ymin-0.1,'0.4', charthick=2.2,charsize=.8
xyouts, 0.5-0.01,ymin-0.1,'0.5', charthick=2.2,charsize=.8
xyouts, 2.0-0.01,ymin-0.1,'2', charthick=2.2,charsize=.8
xyouts, 3.0-0.01,ymin-0.1,'3', charthick=2.2,charsize=.8



;------------Plot SFH -----------------
;=============================================

trange=[1.4,5.8]

if burst eq 'senior' then begin
deltat1=2.0
deltat2=4.
endif

if burst eq 'interm' then begin
deltat1=3.5
deltat2=4.
endif

;para no salirse del plot
t_orig=t

sfr1a_orig=sfr1a
sfr1b_orig=sfr1b
sfr1c_orig=sfr1c
sfr1d_orig=sfr1d

sfr2a_orig=sfr2a
sfr2b_orig=sfr2b
sfr2c_orig=sfr2c
sfr2d_orig=sfr2d



ww1=where(t_orig*1.e-9+deltat1 lt trange(1)) 
w1=n_elements(ww1)-1
t1=t_orig(0:w1)*1.e-9+deltat1

sfr1a=sfr1a_orig(0:w1)
sfr1b=sfr1b_orig(0:w1)
sfr1c=sfr1c_orig(0:w1)
sfr1d=sfr1d_orig(0:w1)

ww2=where(t_orig*1.e-9+deltat2 lt trange(1)) 
w2=n_elements(ww2)-1
t2=t_orig(0:w2)

t2=t_orig(0:w2)*1.e-9+deltat2
sfr2a=sfr2a_orig(0:w2)
sfr2b=sfr2b_orig(0:w2)
sfr2c=sfr2c_orig(0:w2)
sfr2d=sfr2d_orig(0:w2)


if burst eq 'senior' then yrange_sfh=205
if burst eq 'interm' then yrange_sfh=250

print, yrange_sfh
;-----
plot, t*1.e-9, sfr1a, yrange=[0,yrange_sfh],xrange=trange, xtitle='', ytitle='SFR ['+Msun+'/yr]', thick=3., xthick=3., ythick=3., charsize=.8, charthick=2.5 , xtickinterval=1., xminor=1,xstyle=9,  /ystyle, /nodata, /noerase, position=p3a

oplot, t1, sfr1a, col=cgcolor('Blu4'), thick=4.
oplot, t2, sfr2a, col=cgcolor('Blu4'), thick=4.

inf_line1=fltarr(n_elements(t1))
inf_line1[0:n_elements(t1)-1]=0.1

inf_line2=fltarr(n_elements(t2))
inf_line2[0:n_elements(t2)-1]=0.1

cgColorFill, [t1(0),t1,reverse(t1),t1(0)], [sfr1a(0), sfr1a, inf_line1, sfr1a(0)],  COLOR=cgcolor('Blu4')
cgColorFill, [t2(0),t2,reverse(t2),t2(0)], [sfr2a(0), sfr2a, inf_line2, sfr2a(0)],  COLOR=cgcolor('Blu4')

xyouts, 0.11,0.16, /norm,' a) ', col=cgcolor('Blu4'), charthick=2.5, charsize=.8


xval=[getage(5.), getage(3.),getage(2.), getage(1.5)]
lb=[5, 3, 2, 1.5]
axis,xaxis=1,xr=trange,xst=1,xtickv=xval,xtickname=string(lb,f='(F4.1)'), xticks=n_elements(xval)-1,chars=.8, chart=2.5, xthick=2.5

plot, t*1.e-9, sfr1a, yrange=[0,yrange_sfh],xrange=trange, xtitle='', ytitle='SFR ['+Msun+'/yr]', thick=3., xthick=3., ythick=3., charsize=.8, charthick=2.5 , xtickinterval=1., xminor=1,xstyle=9, /ystyle,  /nodata, /noerase, position=p3a


;-----
plot, t*1.e-9, sfr1b, yrange=[0,yrange_sfh],xrange=trange, xtitle='', ytitle='', thick=3., xthick=3., ythick=3., charsize=.8, charthick=2.5 ,ytickformat='(A1)', xtickinterval=1.,xminor=1,xstyle=9, /ystyle,/nodata,/noerase,  position=p3b

oplot, t1, sfr1b, col=cgcolor('TG5'), thick=4.
oplot, t2, sfr2b, col=cgcolor('TG5'), thick=4.


cgColorFill, [t1(0),t1,reverse(t1),t1(0)], [sfr1b(0), sfr1b, inf_line1, sfr1b(0)],  COLOR=cgcolor('TG5')
cgColorFill, [t2(0),t2,reverse(t2),t2(0)], [sfr2b(0), sfr2b, inf_line2, sfr2b(0)],  COLOR=cgcolor('TG5')

xyouts, 0.32,0.16, /norm,' b) ', col=cgcolor('TG5'), charthick=2.5, charsize=.8

xyouts, 4.2,yrange_sfh+50, 'redshift', charthick=2.5, charsize=.8
xyouts, 4.,-80, 'Age of Universe [Gyr]', charthick=2.5, charsize=.8

axis,xaxis=1,xr=trange,xst=1,xtickv=xval,xtickname=string(lb,f='(F4.1)'), xticks=n_elements(xval)-1,chars=.8, chart=2.5, xthick=2.5

plot, t*1.e-9, sfr1b, yrange=[0,yrange_sfh],xrange=trange, xtitle='', ytitle='', thick=3., xthick=3., ythick=3., charsize=.8, charthick=2.5 ,ytickformat='(A1)', xtickinterval=1.,xminor=1,xstyle=9,  /ystyle,/nodata,/noerase,  position=p3b

;----
plot, t*1.e-9, sfr1c, yrange=[0,yrange_sfh],xrange=trange, xtitle='', ytitle='', thick=3., xthick=3., ythick=3., charsize=.8, charthick=2.5 ,ytickformat='(A1)', xtickinterval=1.,xstyle=9,  /ystyle,xminor=1,/nodata, /noerase, position=p3c

oplot, t1, sfr1c, col=cgcolor('RYB2'), thick=4.
oplot, t2, sfr2c, col=cgcolor('RYB2'), thick=4.



cgColorFill, [t1(0),t1,reverse(t1),t1(0)], [sfr1c(0), sfr1c, inf_line1, sfr1c(0)],  COLOR=cgcolor('Ryb2')
cgColorFill, [t2(0),t2,reverse(t2),t2(0)], [sfr2c(0), sfr2c, inf_line2, sfr2c(0)],  COLOR=cgcolor('Ryb2')


xyouts, 0.53,0.16, /norm,' c) ', col=cgcolor('RYB2'), charthick=2.5, charsize=.8

axis,xaxis=1,xr=trange,xst=1,xtickv=xval,xtickname=string(lb,f='(F4.1)'), xticks=n_elements(xval)-1,chars=.8, chart=2.5, xthick=2.5

plot, t*1.e-9, sfr1c, yrange=[0,yrange_sfh],xrange=trange, xtitle='', ytitle='', thick=3., xthick=3., ythick=3., charsize=.8, charthick=2.5 ,ytickformat='(A1)', xtickinterval=1.,xstyle=9,  /ystyle,xminor=1,/nodata, /noerase, position=p3c

;----
plot, t*1.e-9, sfr1d, yrange=[0,yrange_sfh],xrange=trange, xtitle='', ytitle='', thick=3., xthick=3., ythick=3., charsize=.8, charthick=2.5 ,ytickformat='(A1)', xtickinterval=1.,xstyle=9,  /ystyle,xminor=1,/nodata, /noerase, position=p3d

oplot, t1, sfr1d, col=cgcolor('Maroon'), thick=4.
oplot, t2, sfr2d, col=cgcolor('Maroon'), thick=4.


cgColorFill, [t1(0),t1,reverse(t1),t1(0)], [sfr1a(0), sfr1d, inf_line1, sfr1d(0)],  COLOR=cgcolor('Maroon')
cgColorFill, [t2(0),t2,reverse(t2),t2(0)], [sfr2a(0), sfr2d, inf_line2, sfr2d(0)],  COLOR=cgcolor('Maroon')


xyouts, 0.74,0.16, /norm,' d) ', col=cgcolor('Maroon'), charthick=2.5, charsize=.8


axis,xaxis=1,xr=trange,xst=1,xtickv=xval,xtickname=string(lb,f='(F4.1)'), xticks=n_elements(xval)-1,chars=.8, chart=2.5, xthick=2.5

plot, t*1.e-9, sfr1d, yrange=[0,yrange_sfh],xrange=trange, xtitle='', ytitle='', thick=3., xthick=3., ythick=3., charsize=.8, charthick=2.5 ,ytickformat='(A1)', xtickinterval=1.,xstyle=9,  /ystyle,xminor=1,/nodata, /noerase, position=p3d

 Device, /Close_File
 set_plot,'x'

print, 'FIN plot '

stop




END

