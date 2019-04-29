;goto,jdmskip
;goto,jdmskip2
; 2004Ju16 JDM	Modified to use oifits data as input.
; ** NEED to ADD something to deal with closure phase error offset
; for broad band data!
; 2005May29 JDm Added systematic + calsize error into error bar.
; 2006Nov20 JDM Generalized for MIRC.
; 2006Dec03 JDM more and more more...
; 2006Dec21 JDM even more.
; 2006Dec26 JDM added scatter error determination
; 2007Aug07 JDM to be added:
;	      adding daqsync option for normalization.
; 2008Jul17 JDM adding better treatment of correlated errors due to cal diam errors
; 2009Sep19 JDM xchan
; 2009Oct25 JDM deadling with data with NO chopper data.
;	?     ability to plot rawvis power vis photometry within a block
; 	?     postscript output
;	?     ability to save and restore all the clicking!!
;               future calibrator studies. NOCAL_ prefix.
; 2010Feb09 JDM fixed bugs in t3amp flagging.
; 2010Feb13 JDM added better autoflaggin when photometry estimates fail (==0.0)
;		fixed some bugs where did not handle flags in some channels right
;		if there were some NANs or 0000s in the photometry files.
;;		.. lots of small changes to make things consistent when t3amps are flagged.
; 2010May04 JDM will automatically create RAW files without any calibration.
;2011Jul12  JDM working to adapt for 6 telescopes!! lots of data!
;2012Apr24  JDM mirc_average created from mirc_reduce
;		will not depend on calibrator SIZES
;		will not do an AVG output.
;		will use new prefix. MIRC_L1 (Following level1, level2, distinction of PTI/KI)
;		will use MINFILE and MAXFILE so integration times are never too long...
; 2016Apr06 GHS Added option to reject outliers using sigma clipping.
; 2016Sep15 GHS Added sigma clipping algorithm to pipeline_mirc6b
;
skipahead=0
rerun:
;goto,plotonly
;goto,skipv2
; TO FIX
;  1 Must ADD CALIBRATOR error later after the caltables... somehow.
;   can't be treated like an extra error since there is too much averaging with MIRC
;   data
;  2 when plotting data, divide by a smoothed version of the photometry signals since we can tolerate low
;    snr on flux estimate in individual files.


;restore,'testing_altair.idlvar'


noclick=0 ; 0 = interactive,  1=automatic (for testing)


          ; USER ADJUSTABLE !! could be made adaptable.
minfiles=5 ; requires this number of files for averaging [over-rides NSPLIT]

;maxfiles=15 ; could be changed by hand for crazy sources.. like very wide field sources for instance!
avg_time=2.5 ; minutes
nsplit_max=1000; just for setting up arrays.
;  smoothlength = float(dialog_input(prompt="Smoothing (1/3) Length (hrs):",initial=smoothlength))


mode2 = 'split' ;+strtrim(string(nsplit),2)


use_scatter =0 ;

hey=dialog_checklist(['Normal Data Selection (RECOMMENDED)','FAST-TRACK (quicklook only)'], tit='Choose level of data selection:',$
     /exclusive, initial=0)

q=reform(where(hey eq 1))
if q eq 0 then FASTTRACK=0
if q eq 1 then FASTTRACK=1

if q eq 2 then begin
restoreflags = 1
restore,'lastflags_all.idlvar'
endif else restoreflags=0



smooth_fun=3./60 ; 3 minute averaging for displaying data.

output_path="./"
oifits_path='./OIFITS/' ; save oifits here. end path with /
;oifits_path='./'
infile=dialog_pickfile(tit="Choose date.idlvar file from mirc_process2:",filter="*_info.idlvar")
restore,infile

in=strpos(infile,'/',/reverse_search)
if in eq -1 then path = './'
if in ne -1 then path = strmid(infile,0,in[0]+1)

mirclog = mirc_dirswap(mirclog,path,filename=filename)

read_mirclog,mirclog,instructure0,num_data=num_data ; num blocks

if (keyword_set(num_data) eq 0 )then num_data=instructure[0].num_datablocks

;--------------------------------------------
;read_mircarray,file_mircarray,mircarray

; The following information is in mircarray now.
ntel=fix(mircarray.mirc)

nbaselines=ntel*(ntel-1)/2
ntriangles=ntel*(ntel-1)*(ntel-2)/6
if ntel eq 4 then begin
    ihelp=[0,0,0,1,1,2] ; map (tel i, tel j) -> baseline #
    jhelp=[1,2,3,2,3,3]
    itrihelp=[0,0,0,1] ; tel-> tri
    jtrihelp=[1,1,2,2]
    ktrihelp=[2,3,3,3]
    itribase=[0,0,1,3]	      ; base-> tri
    jtribase=[3,4,5,5]
    ktribase=[1,2,2,4]
endif
if ntel eq 6 then begin
ntel=mircarray.ntel
nbaselines=mircarray.nbaselines
ntriangles=mircarray.ntriangles
nquads = mircarray.nquads
  ihelp=*mircarray.ihelp
    jhelp=*mircarray.jhelp
    itrihelp=*mircarray.itrihelp
    jtrihelp=*mircarray.jtrihelp
    ktrihelp=*mircarray.ktrihelp
    itribase=*mircarray.itribase
    jtribase=*mircarray.jtribase
    ktribase=*mircarray.ktribase
    ; FB 2016
    iquadhelp=*mircarray.iquadhelp
    jquadhelp=*mircarray.jquadhelp
    kquadhelp=*mircarray.kquadhelp
    lquadhelp=*mircarray.lquadhelp
    iquadbase=*mircarray.iquadbase
    jquadbase=*mircarray.jquadbase
    kquadbase=*mircarray.kquadbase
    lquadbase=*mircarray.lquadbase
endif

;---------------------------------------------


pos=strpos(infile,'_info.idlvar')
fluxfile=strmid(infile,0,pos)+'_flux.idlvar'
if file_test(fluxfile) then  restore,fluxfile else begin
 chopon=0
 daqsync=0
endelse


pos=strpos(infile,'_info.idlvar')
xchanfile=strmid(infile,0,pos)+'_xchan.idlvar'
if file_test(xchanfile) then restore,xchanfile
if keyword_set(xchan) eq 0 then xchan=0


norm_list=['SHUTTER','FIBER']
if chopon eq 1 then norm_list=[norm_list,'CHOP']
if daqsync eq 1 then norm_list=[norm_list,'DAQSYNC']
if xchan eq 1 then norm_list=[norm_list,'XCHAN']

result=dialog_checklist(norm_list,tit="Choose NORM Methods to Analyze:",/nonexclusive,$
  initial=-1 )
if total(result) eq 0 then begin
 print,'must select at least one method!!'
 stop
endif
norm_labels=norm_list(where(result eq 1,normct))
num_norms=normct
;stop ;JDM


if keyword_set(comboinfo) eq 0 then begin
combofile=strmid(infile,0,pos)+'_comboinfo.idlvar'
if nlines(combofile) ne 0 then restore,combofile else begin
 combofile=dialog_pickfile(tit="Choose combofile:",$
 filter="*_comboinfo.idlvar")
endelse
restore,combofile
endif

bestfp=comboinfo.bestfp
bestlambda=comboinfo.bestlambda

; JDM
temp_string='ls '+strmid(infile,0,pos)+'*_reduced.idlvar'
spawn,temp_string,list
if n_elements(list) eq 1 then begin
   reducefile=list[0]
 endif else  begin
result=dialog_checklist(list,tit="Choose REDUCTION:",/exclusive,$
  initial=-1 )
in=where( result eq 1)
reducefile=list[in]

endelse

;outfile=strmid(infile,0,pos)+'_reduced.idlvar'
restore,reducefile

tn=tag_names(psinfo)
test=where(tn eq 'SPLODGES_BASE_V2',ct)
if ct eq 1 then begin

hey=dialog_checklist(['Default Selection Display (default)','New cross-talk corrected Selection Display'], tit='Choose display type:',$
     /exclusive, initial=0)

q=reform(where(hey eq 1))
if q eq 0 then newsplodge=0
if q eq 1 then newsplodge=1

endif


original_vis2info=vis2info
original_tpinfo=tpinfo

dins=where(instructure.type eq 'DATA')
nblocks=num_data
blocknums=instructure(dins).block
starnames=instructure(dins).star

targets=strcompress(starnames + '(' + string(fix(blocknums))+')',/rem)
targets_file=strcompress(starnames + '.' + string(fix(blocknums))+'.',/rem)

print,'Targets: ',targets
avg_time= float(dialog_input(prompt="Maximum Averaging Time (minutes)",initial=avg_time))

if avg_time lt 0.5 then minfiles=3 ;


; make copy of vis2info and tpinfo for storing information
; on the errors due to diameter errors -- this will be inefficient since most of the
; structure's entries are not needed to be duplicated..
vis2info_cal=vis2info
tpinfo_cal=tpinfo

header0={source:'  ', calib:0, target:0, $
         calsize:0.00, calsizeerr:0.00,block:0}
headers=replicate(header0,nblocks)
headers.source = strcompress(starnames)
headers.block = blocknums

headers(*).target=1
ntargets=n_elements(targets)

; Now I will go through and strip data out of the array of pointers.
; This is a pain but probably necessary.
nfiles =0
for b=0,nblocks-1 do begin
  nfiles0=n_elements( *misc(b).num_data)
  nfiles=nfiles+nfiles0
endfor

;get some info
;    ave_bg   = *misc(0).ave_bg
;    dims=size(ave_bg,/dim)
;    dimx=dims(0)
;    dimy=dims(1)
;    dimz=dims(2)
;    u=fltarr(nbaselines,dimy)
;    v=fltarr(nbaselines,dimy)
dimy = comboinfo.nwave


; nfiles SHOULD EQUAL the number of elements in in UVINFO. except for badfiles.

vis2_template = {uthrs:0.0, vis2:fltarr(nbaselines,dimy), vis2_err:fltarr(nbaselines,dimy), diam_vis2_err:fltarr(nbaselines,dimy),diam_cal:'',$
   uvin:0, block:0, file:0, target:0, calib:0,badflag:intarr(nbaselines,dimy)}
rawvis2_all=replicate(vis2_template,nfiles)

original_rawvis2_all=rawvis2_all

calvis2_shut=replicate(vis2_template,nfiles)
calvis2_chop=replicate(vis2_template,nfiles)
calvis2_fiber=replicate(vis2_template,nfiles)
calvis2_xchan=replicate(vis2_template,nfiles)
calvis2_daqsync=replicate(vis2_template,nfiles)
splodges_base = fltarr(nfiles,dimy*nbaselines)
splodges_tel = fltarr(nfiles,dimy*ntel)

flux_tel  = fltarr(nfiles,ntel*num_norms)  ; 3= number of photometry methods.

c=0

for b=0,nblocks-1 do begin
  nfiles0=n_elements( *misc(b).num_data)
  rawvis2=(*vis2info(b).rawvis2)
  rawvis2_err = (*vis2info(b).rawvis2_err)
  original_rawvis2=(*original_vis2info(b).rawvis2)
  original_rawvis2_err = (*original_vis2info(b).rawvis2_err)

  rawvis2_err_cal= (*vis2info_cal(b).rawvis2_err)
  datafiles=*misc(b).num_data

  for f=0,nfiles0-1 do begin
   uvin=where(uvinfo.star eq headers(b).source and uvinfo.block eq headers(b).block and uvinfo.filenum eq datafiles(f)) ;crosscheck
   rawvis2_all(c).uthrs=uvinfo(uvin).time/3600. ; assumes all data in night on same utdate.
   rawvis2_all(c).vis2=rawvis2(*,*,f)
   rawvis2_all(c).vis2_err=rawvis2_err(*,*,f)
   original_rawvis2_all(c).vis2=original_rawvis2(*,*,f)
   original_rawvis2_all(c).vis2_err=original_rawvis2_err(*,*,f)

   rawvis2_all(c).diam_vis2_err=rawvis2_err_cal(*,*,f)
   rawvis2_all(c).uvin=uvin(0)
   rawvis2_all(c).block=headers(b).block
   rawvis2_all(c).file=datafiles(f)
   rawvis2_all(c).target=headers(b).target
   rawvis2_all(c).calib=headers(b).calib
   rawvis2_all(c).diam_cal=headers(b).source

   if keyword_set(xchan) then calvis2_xchan(c).vis2 = (*vis2info(b).vis2_norms_xchan)(*,*,f)
   if keyword_set(chopon) then calvis2_chop(c).vis2 = (*vis2info(b).vis2_norms_chop)(*,*,f)
   if keyword_set(daqsync) then calvis2_daqsync(c).vis2 = (*vis2info(b).vis2_norms_daqsync)(*,*,f)
   calvis2_shut(c).vis2 = (*vis2info(b).vis2_norms_shut)(*,*,f)
   calvis2_fiber(c).vis2 = (*vis2info(b).vis2_norms_fiber)(*,*,f)

   splodges_base(c,*)=(*psinfo(b).splodges_base)(f,*)
   splodges_tel(c,*)=(*psinfo(b).splodges_tel)(f,*)
   tn=tag_names(psinfo)
   test=where(tn eq 'SPLODGES_BASE_V2',ct)
   if ct eq 1 and newsplodge eq 1 then begin
      splodges_base(c,*)=(*psinfo(b).splodges_base_v2)(f,*)
      splodges_tel(c,*)=(*psinfo(b).splodges_tel_v2)(f,*)
   endif


   for t=0,ntel-1 do begin

      if t eq 0 then begin
        if chopon then  temp_chop     = (*flux_chop(b).b1)(f)
        if daqsync then temp_daqsync  = (*flux_daqsync(b).b1)(f)
        if xchan then   temp_xchan    = (*flux_xchan(b).b1)(f)
        temp_fiber = (*flux_fiber(b).b1)(f)
        temp_shut  = (*flux_shut(b).b1)(f)
      endif
     if t eq 1 then begin
       if chopon then  temp_chop     = (*flux_chop(b).b2)(f)
        if daqsync then temp_daqsync  = (*flux_daqsync(b).b2)(f)
        if xchan then   temp_xchan    = (*flux_xchan(b).b2)(f)
        temp_fiber = (*flux_fiber(b).b2)(f)
        temp_shut  = (*flux_shut(b).b2)(f)
      endif
     if t eq 2 then begin
        if chopon then  temp_chop     = (*flux_chop(b).b3)(f)
        if daqsync then temp_daqsync  = (*flux_daqsync(b).b3)(f)
        if xchan then   temp_xchan    = (*flux_xchan(b).b3)(f)
        temp_fiber = (*flux_fiber(b).b3)(f)
        temp_shut  = (*flux_shut(b).b3)(f)
      endif
     if t eq 3 then begin
        if chopon then  temp_chop     = (*flux_chop(b).b4)(f)
        if daqsync then temp_daqsync  = (*flux_daqsync(b).b4)(f)
        if xchan then   temp_xchan    = (*flux_xchan(b).b4)(f)
        temp_fiber = (*flux_fiber(b).b4)(f)
        temp_shut  = (*flux_shut(b).b4)(f)
      endif
     if t eq 4 then begin
        if chopon then  temp_chop     = (*flux_chop(b).b5)(f)
        if daqsync then temp_daqsync  = (*flux_daqsync(b).b5)(f)
        if xchan then   temp_xchan    = (*flux_xchan(b).b5)(f)
        temp_fiber = (*flux_fiber(b).b5)(f)
        temp_shut  = (*flux_shut(b).b5)(f)
      endif
     if t eq 5 then begin
        if chopon then  temp_chop     = (*flux_chop(b).b6)(f)
        if daqsync then temp_daqsync  = (*flux_daqsync(b).b6)(f)
        if xchan then   temp_xchan    = (*flux_xchan(b).b6)(f)
        temp_fiber = (*flux_fiber(b).b6)(f)
        temp_shut  = (*flux_shut(b).b6)(f)
      endif


for nn=0,num_norms-1 do begin
    if norm_labels[nn] eq 'SHUTTER'  then $
     flux_tel(c,t*num_norms+nn) = temp_shut
    if norm_labels[nn] eq 'FIBER'  then $
     flux_tel(c,t*num_norms+nn) = temp_fiber
    if norm_labels[nn] eq 'CHOP'  then $
     flux_tel(c,t*num_norms+nn) = temp_chop
    if norm_labels[nn] eq 'DAQSYNC'  then $
     flux_tel(c,t*num_norms+nn) = temp_daqsync
    if norm_labels[nn] eq 'XCHAN'  then $
     flux_tel(c,t*num_norms+nn) = temp_xchan
endfor
   endfor; tel
   c=c+1
  endfor;

 endfor;blocks

 ; Do normalization for flux_tel
; Normalize each channel by median of the flux for all measurements of the same SOURCE.
uniq_sources = (headers.source)(uniq(headers.source,sort(headers.source)))
uvgoodlist = rawvis2_all.uvin


for us=0,n_elements(uniq_sources)-1 do begin

      target_index = where(uvinfo(uvgoodlist).star eq uniq_sources(us))

  for t=0,ntel-1 do begin
     normvalue = median(flux_tel(target_index,(t*num_norms):(t*num_norms+num_norms-1)))
     ;JDM inconsistent.
     ; norm by median of all norm methods -- not individually.
   for ps=0,num_norms-1 do begin
     normvalue = median(flux_tel(target_index,t*num_norms+ps) ) ; each norm method separately
     if normvalue eq 0 then normvalue=1.0
     flux_tel(target_index,t*num_norms+ps) = flux_tel(target_index,t*num_norms+ps) / $
            normvalue
   endfor
  endfor
endfor

flux_tel_flag = flux_tel ; Setup for Flagging too: (nfiles,ntel*num_norms)
flux_tel_flag(*)=0


 ;Check for NANs.... should be removed elsewhere -- print WARNING!
 ;  Also mark bad if calvis2 = 0.00 exactly or exactly -1. indicates problem with photometry

; JDM
; only bother to check if tick mark!!

inshut=where(norm_labels eq 'SHUT')
infiber=where(norm_labels eq 'FIBER')
inchop=where(norm_labels eq 'CHOP')
indaqsync=where(norm_labels eq 'DAQSYNC')
inxchan=where(norm_labels eq 'XCHAN')

if inshut ne -1 then begin
   in=where(calvis2_shut.vis2 ne calvis2_shut.vis2 or $
            calvis2_shut.vis2 eq 0.000 or $
            calvis2_shut.vis2 eq -1,ct)
   if ct gt 0 then begin
     print,' NANs or 000000s or -1 in calvis2_shut!: ',ct
     ac=array_coords(in, (calvis2_shut.vis2) )
     calvis2_shut(ac(*,2)).vis2(ac(*,0),ac(*,1)) = -1
    ;alarmbeep,10
   endif
endif

if infiber ne -1 then begin
   in=where(calvis2_fiber.vis2 ne calvis2_fiber.vis2 or $
            calvis2_fiber.vis2 eq 0.000 or $
            calvis2_fiber.vis2 eq -1,ct)
   if ct gt 0 then begin
     print,' NANs or 0000000s  or -1 in calvis2_fiber!: ',ct
      ac=array_coords(in, (calvis2_fiber.vis2) )
     calvis2_fiber(ac(*,2)).vis2(ac(*,0),ac(*,1)) = -1
     ;alarmbeep,10
   endif
endif

if inchop ne -1 then begin
   in=where(calvis2_chop.vis2 ne calvis2_chop.vis2 or $
            calvis2_chop.vis2 eq 0.0000 or $
            calvis2_chop.vis2 eq -1.0,ct)
   if ct gt 0 then begin
     print,' NANs or 000000s or -1 in calvis2_chop!: ',ct
      ac=array_coords(in, (calvis2_chop.vis2) )
     calvis2_chop(ac(*,2)).vis2(ac(*,0),ac(*,1)) = -1
     ;alarmbeep,10
   endif
endif

if indaqsync ne -1 then begin
   in=where(calvis2_daqsync.vis2 ne calvis2_daqsync.vis2 or $
            calvis2_daqsync.vis2 eq 0.000 or $
            calvis2_daqsync.vis2 eq -1,ct)
   if ct gt 0 then begin
     print,' NANs or 00000s or -1 in calvis2_daqsync!: ',ct
      ac=array_coords(in, (calvis2_daqsync.vis2) )
     calvis2_daqsync(ac(*,2)).vis2(ac(*,0),ac(*,1)) = -1
     ;alarmbeep,10
   endif
endif

if inxchan ne -1 then begin
   in=where(calvis2_xchan.vis2 ne calvis2_xchan.vis2 or $
            calvis2_xchan.vis2 eq 0.000 or $
            calvis2_xchan.vis2 eq -1,ct)
   if ct gt 0 then begin
     print,' NANs or 00000s or 01 in calvis2_xchan!: ',ct

;     ac=array_coords(in, (calvis2_xchan.vis2) )
;print,'method1' ; new faster method!
     test = calvis2_xchan.vis2
     test[in]=-1
     calvis2_xchan.vis2=test
print,'method2'
;     calvis2_xchan(ac(*,2)).vis2(ac(*,0),ac(*,1)) = -1 ;SLOW STEP
     ;alarmbeep,10

   endif
endif



inc=where(rawvis2_all.calib eq 1)
int=where(rawvis2_all.target eq 1)
intc=where(rawvis2_all.calib eq 1 or rawvis2_all.target eq 1)


vis2_chop=rawvis2_all
vis2_shut=rawvis2_all
vis2_fiber=rawvis2_all
vis2_daqsync=rawvis2_all
vis2_xchan=rawvis2_all

spec=fix(dimy/2)-1. ; pick middle of spectral channel for doing estimates.
;spec=4 ; testing..
smoothlength=2.0 ;hrs

;----------------
; FLAG all files based on appearance of FRINGES or not
;----------------
skiptemp:
splodges_base_flag = splodges_base
splodges_tel_flag  = splodges_tel
splodges_base_flag(*)=0 ; 0=good, 1=bad
splodges_tel_flag(*)=0

splodges_base_autoflag = splodges_base
splodges_tel_autoflag  = splodges_tel
splodges_base_autoflag(*)=0 ; 0=good, 1=bad
splodges_tel_autoflag(*)=0

;restore,'temp_flag.idlvar'

;setup markers
temp=rawvis2_all.block - shift(rawvis2_all.block,-1)
markers = where(temp ne 0,ctmark)+1
;ctmark=ctmark+1
;markers=[markers,nfiles]
;if ctmark ne (nblocks) then stop ; PROBLEM!
if ctmark eq 0 then ctmark=1 ;JDM . in case of only 1 block

cleanplot
erase
!p.multi=[0,1,2,0,0]

loadct,3 ; or 10?
loadct,10
;image_cont,jhist(splodges_base),/noc,color=
xsize=!d.x_size*.9
ysize=!d.y_size*.45

imdisp,jhist(splodges_base)*(1-splodges_base_flag),  position=[.15,.55,.95,.95],outpos=basepos,/usepos,/erase,$
   $aspect = 1.2*ysize/xsize,$
     /axis,ticklen=0, yticks=5,ytickname=replicate(' ',6),$
    xticks=5, xtickname=replicate(' ',6)

for nm=0,ctmark-1 do begin
  oplot,[markers(nm),markers(nm)],[0,1e6]
  if (nm mod 4 eq 0) then xyouts,markers(nm),0,'!C'+targets(nm),align=1
  if (nm mod 4 eq 1) then xyouts,markers(nm),0,'!C!C'+targets(nm),align=1
  if (nm mod 4 eq 2) then xyouts,markers(nm),0,'!C!C!C'+targets(nm),align=1
  if (nm mod 4 eq 3) then xyouts,markers(nm),0,'!C!C!C!C'+targets(nm),align=1
endfor

for nm=0,nbaselines-1 do begin
  oplot,[0,1e6],[dimy*nm,dimy*nm],thick=3
  xyouts,0,dimy*nm+dimy/2,strcompress(get_basename(nm,mircarray=mircarray),/rem)+'  ',align=1,charsize=2
endfor

baseplot_x =!X
baseplot_y =!Y
;----------------- TEL PLOT -----------
imdisp,jhist(splodges_tel)*(1-splodges_tel_flag), position=[.15,.05,.95,.45],/usepos,$
   aspect = 1.2*ysize/xsize, $
    /axis, ticklen=0, yticks=5,ytickname=replicate(' ',6),$
    xticks=5, xtickname=replicate(' ',6)
for nm=0,ctmark-1 do begin
  oplot,[markers(nm),markers(nm)],[0,1e6]
;  if (nm mod 4 eq 0) then xyouts,markers(nm),0,'!C'+targets(nm),align=1
;  if (nm mod 4 eq 1) then xyouts,markers(nm),0,'!C!C'+targets(nm),align=1
;  if (nm mod 4 eq 2) then xyouts,markers(nm),0,'!C!C!C'+targets(nm),align=1
;  if (nm mod 4 eq 3) then xyouts,markers(nm),0,'!C!C!C!C'+targets(nm),align=1
endfor

for nm=0,ntel-1 do begin
  oplot,[1,1e6],[dimy*nm,dimy*nm],thick=3
  if nm eq 0 then xyouts,0,dimy*nm+dimy/2,mircarray.b1,align=1,charsize=4
  if nm eq 1 then xyouts,0,dimy*nm+dimy/2,mircarray.b2,align=1,charsize=4
  if nm eq 2 then xyouts,0,dimy*nm+dimy/2,mircarray.b3,align=1,charsize=4
  if nm eq 3 then xyouts,0,dimy*nm+dimy/2,mircarray.b4,align=1,charsize=4
  if nm eq 4 then xyouts,0,dimy*nm+dimy/2,mircarray.b5,align=1,charsize=4
  if nm eq 5 then xyouts,0,dimy*nm+dimy/2,mircarray.b6,align=1,charsize=4

endfor
telplot_x=!x
telplot_y=!y

splodges_base_flag=splodges_base_autoflag
splodges_tel_flag=splodges_tel_autoflag
splodges_vis2_flag=fltarr(nbaselines,nfiles) ; 0=good, 1=bad
splodges_t3_flag  =fltarr(ntriangles,nfiles) ; 0=good, 1=bad
splodges_t4_flag  =fltarr(nquads,nfiles) ; 0=good, 1=bad

;********************************************
;offer chance for automatic removal of areas of low flux on the telescopes
;********************************************
ans_cut = 'N'
;print,'Would you like to apply automatic clipping to the telescope fluxes?'
;print,'<N> or Y'
;read,ans_cut
if (ans_cut ne 'Y') then goto,splodge_replot
fluxclip:

sigrej = 5.0
ansrej = 'N'
print,'Rejection threshold (Max Flux/ Threshold): '
print,sigrej
print,'Would you like to change rejection threshold?'
print,'Y of <N>'
read,ansrej
if (ansrej eq 'Y') then begin
   print,'Enter new rejection threshold:'
   read,sigrej
endif

;replot_fluxtel
splodges_base_autoflagtemp=splodges_base_autoflag
ans_cont = ''
bounds=[-1,markers]
for i=1,ctmark do begin
  for j=1,ntel do begin
    if bounds[i-1] GT 0 then x0=bounds[i-1] else x0=bounds[i-1]+1
    ;if i LT (ctmark-1) then x1=bounds[i] else x1=bounds[i]-1
    x1=bounds[i]-1
    y0=(dimy*(j-1)) ;bottom of telescope zone
    y1=(dimy*(j)-1) ;top of telescope zone
    ;x0=bounds[i-1]
    ;x1=bounds[i]
    ;y0=(dimy*(j-1)) ;bottom of telescope zone
    ;y1=(dimy*(j)-1) ;top of telescope zone
    print,"x0 ",x0
    print,"x1 ",x1
    help,splodges_tel
    ypick=fix (mean([y0,y1])/dimy)
    yinr=[ypick*dimy,ypick*dimy+dimy-1]
    telf=splodges_tel[x0:x1,y0:y1] ;telscope obs
    coll_telf = total(telf,2) ;flux
    corr_telf_low = coll_telf LT max(coll_telf)/sigrej ;max flux comparison
    s=size(corr_telf_low)
    for k=y0,y1 do begin
      splodges_tel_autoflag[x0:x1,k]=corr_telf_low
      ;print,splodges_tel_autoflag[x0:x1,k]
    endfor
    ; now must also do all the BASELINES with this TEL.
    in=[where(ihelp eq ypick),where(jhelp eq ypick)]
    in=in(where(in ne -1,ct))
    print,"ct: ",ct
    print,"dimy ",dimy
    for tempin=0,ct-1 do begin
      print,'flagging bad baselines: ',in[tempin], "name: ",strcompress(get_basename(in[tempin],mircarray=mircarray),/rem)
      ytemp0=in(tempin)*dimy
      ytemp1=ytemp0+dimy-1
      for k=ytemp0,ytemp1 do begin
          ;print,"k ", k," corr_telf_low ", corr_telf_low
          splodges_base_autoflagtemp[x0:x1,k]=corr_telf_low
      endfor
    endfor
      print,"new telescope j ", j
        splodges_base_autoflag=splodges_base_autoflagtemp+splodges_base_autoflag
        splodges_base_autoflag[where(splodges_base_autoflag GT 1)]=1
  endfor

endfor
splodges_base_flag=splodges_base_autoflag
splodges_tel_flag=splodges_tel_autoflag
splodges_vis2_flag=fltarr(nbaselines,nfiles) ; 0=good, 1=bad
splodges_t3_flag  =fltarr(ntriangles,nfiles) ; 0=good, 1=bad
;now try to get the BASELINES from the TEL
for base=0,nbaselines-1 do begin
 badin=where(splodges_base_flag(*,base*dimy) eq 1,ctbad)
 if ctbad gt 0 then begin
   splodges_vis2_flag(base,badin)=1
   ;Next flag all Triangles that include this baseline.
   in=[where(itribase eq base),where(jtribase eq base),where(ktribase eq base)]
   in=in(where(in ne -1,ct))
   for tempin=0,ct-1 do begin ;
     splodges_t3_flag(in(tempin),badin)=1
   endfor
  endif ;bad
endfor

;-----------------------
 ;Apply splodge flags to the
 ;vis2 data.
;
rawvis2_all_original=rawvis2_all
for s=0,dimy-1 do begin
 rawvis2_all_original.badflag(*,s) = splodges_vis2_flag
endfor
;-------------

;REPLOT
splodge_replot:

cleanplot
erase
!p.multi=[0,1,2,0,0]

loadct,3 ; or 10?
loadct,10
;image_cont,jhist(splodges_base),/noc,color=
xsize=!d.x_size*.9
ysize=!d.y_size*.45

imdisp,jhist(splodges_base)*(1-splodges_base_flag),  position=[.15,.55,.95,.95],outpos=basepos,/usepos,/erase,$
   $aspect = 1.2*ysize/xsize,$
     /axis,ticklen=0, yticks=5,ytickname=replicate(' ',6),$
    xticks=5, xtickname=replicate(' ',6)

for nm=0,ctmark-1 do begin
  oplot,[markers(nm),markers(nm)],[0,1e6]
  if (nm mod 4 eq 0) then xyouts,markers(nm),0,'!C'+targets(nm),align=1
  if (nm mod 4 eq 1) then xyouts,markers(nm),0,'!C!C'+targets(nm),align=1
  if (nm mod 4 eq 2) then xyouts,markers(nm),0,'!C!C!C'+targets(nm),align=1
  if (nm mod 4 eq 3) then xyouts,markers(nm),0,'!C!C!C!C'+targets(nm),align=1
endfor

for nm=0,nbaselines-1 do begin
  oplot,[0,1e6],[dimy*nm,dimy*nm],thick=3
  xyouts,0,dimy*nm+dimy/2,strcompress(get_basename(nm,mircarray=mircarray),/rem)+'  ',align=1,charsize=2
endfor

baseplot_x =!X
baseplot_y =!Y
;----------------- TEL PLOT -----------
imdisp,jhist(splodges_tel)*(1-splodges_tel_flag), position=[.15,.05,.95,.45],/usepos,$
   aspect = 1.2*ysize/xsize, $
    /axis, ticklen=0, yticks=5,ytickname=replicate(' ',6),$
    xticks=5, xtickname=replicate(' ',6)
for nm=0,ctmark-1 do begin
  oplot,[markers(nm),markers(nm)],[0,1e6]
;  if (nm mod 4 eq 0) then xyouts,markers(nm),0,'!C'+targets(nm),align=1
;  if (nm mod 4 eq 1) then xyouts,markers(nm),0,'!C!C'+targets(nm),align=1
;  if (nm mod 4 eq 2) then xyouts,markers(nm),0,'!C!C!C'+targets(nm),align=1
;  if (nm mod 4 eq 3) then xyouts,markers(nm),0,'!C!C!C!C'+targets(nm),align=1
endfor

for nm=0,ntel-1 do begin
  oplot,[1,1e6],[dimy*nm,dimy*nm],thick=3
  if nm eq 0 then xyouts,0,dimy*nm+dimy/2,mircarray.b1,align=1,charsize=4
  if nm eq 1 then xyouts,0,dimy*nm+dimy/2,mircarray.b2,align=1,charsize=4
  if nm eq 2 then xyouts,0,dimy*nm+dimy/2,mircarray.b3,align=1,charsize=4
  if nm eq 3 then xyouts,0,dimy*nm+dimy/2,mircarray.b4,align=1,charsize=4
  if nm eq 4 then xyouts,0,dimy*nm+dimy/2,mircarray.b5,align=1,charsize=4
  if nm eq 5 then xyouts,0,dimy*nm+dimy/2,mircarray.b6,align=1,charsize=4

endfor
telplot_x=!x
telplot_y=!y

print,' Double click in a square to flag whole block as bad (<0.5 sec)'
print,' Click twice over a range to flag a specific RANGE '
print,' Click to the right of plot to REPLOT (allows to resize window)'
print,' Click below all plots to reset bad flags!!'
print,' Click to the left of plots to continue..'
print,' Click above the plots to clip based on distance from average flux'


dclick = 0.5 ; set time for double click
;dclick = 1.0 ; at Nice lag


if  restoreflags eq 1 then goto,skipplot1
cursor,x0,y0,/down
cursor,x0,y0,/up

time0=systime(1)
if y0 gt dimy*13+dimy/2 then begin
  ans_cut='Y'
  goto, fluxclip ; above plot.. go to fluxclipping
endif

if x0 gt nfiles then goto,splodge_replot
if y0 lt !y.crange(0) then begin
  splodges_base_flag(*)=0
  splodges_tel_flag(*)=0
  goto,splodge_replot
endif


if x0 ge 0 then begin ; selection!
print,' Next Click...'
cursor,x1,y1,/down
cursor,x1,y1,/up
time1=systime(1)
print,"x0: ",x0," x1 ",x1
; convert y0 to the correct  value
if y0 le (dimy*ntel) then figclick=0 $; =0 tel, =1 base
 else begin  ; top click
 figclick=1
 normxy0=convert_coord(x0,y0,/data,/to_norm)
 normxy1=convert_coord(x1,y1,/data,/to_norm)
 !x=baseplot_x
 !y=baseplot_y
 newxy0=convert_coord(normxy0(0),normxy0(1),/norm,/to_data)
 newxy1=convert_coord(normxy1(0),normxy1(1),/norm,/to_data)
 x0=newxy0(0)
 x1=newxy1(0)
 y0=newxy0(1)
 y1=newxy1(1)
 print,"x0: ",x0," x1 ",x1
 !x=telplot_x ;return
 !y=telplot_y
 endelse


; now re-order
; First Determine which block and which baseline/tel.
  bpick=rawvis2_all(nint( mean([x0,x1]))).block
  ypick=fix (mean([y0,y1])/dimy) ;just sets to mean of y clicks
  yinr=[ypick*dimy,ypick*dimy+dimy-1] ;

 if (time1-time0 le dclick) then begin;  Double Click!
  xin=where(rawvis2_all.block eq bpick)
 endif else begin ; two single clicks
  print, 'two single clicks!'

  xin0=indgen( abs(x1-x0)+1)+nint(min([x0,x1]))
  ; make sure in range.
  xin0=xin0(where(xin0 ge 0 and xin0 lt nfiles))
  ; make sure all in the same block bick
  xin=xin0(where(rawvis2_all(xin0).block eq bpick,ct))
  if ct eq 0 then goto,splodge_replot  ; error catch
 endelse

 xinr=[min(xin),max(xin)]
 if figclick eq 1 then begin
   splodges_base_flag(xinr(0):xinr(1),yinr(0):yinr(1))=1
 endif
 if figclick eq 0 then begin
   splodges_tel_flag(xinr(0):xinr(1),yinr(0):yinr(1))=1
   ; now must also do all the BASELINES with this TEL.
   in=[where(ihelp eq ypick),where(jhelp eq ypick)]
   in=in(where(in ne -1,ct))
   for tempin=0,ct-1 do begin
     print,'flagging bad baselines: ',in[tempin]
     ytemp0=in(tempin)*dimy
     ytemp1=ytemp0+dimy-1
     print,"x0: ",xinr(0)," x1 ",xinr(1)
     splodges_base_flag(xinr(0):xinr(1),ytemp0:ytemp1)=1
   endfor
 endif

endif  ; clicked in xrange


if x0 ge 0 then goto,splodge_replot ; else continue!

save,splodges_base_flag,splodges_tel_flag,file='temp_flag.idlvar'
skipplot1:

;----------------------------------------------------
; OK .. Now take results of the splodge clicks and convert to
; a more usable format.
;

; Easiest way to proceed to do this baseline by baselines.
for base=0,nbaselines-1 do begin
 badin=where(splodges_base_flag(*,base*dimy) eq 1,ctbad)
 if ctbad gt 0 then begin
   splodges_vis2_flag(base,badin)=1
   ; Next flag all Triangles that include this baseline.
   in=[where(itribase eq base),where(jtribase eq base),where(ktribase eq base)]
   in=in(where(in ne -1,ct))
   for tempin=0,ct-1 do begin ;
     splodges_t3_flag(in(tempin),badin)=1
   endfor
 endif ;bad
endfor

;-----------------------
; Apply splodge flags to the
; vis2 data.
;
rawvis2_all_original=rawvis2_all
for s=0,dimy-1 do begin
 rawvis2_all_original.badflag(*,s) = splodges_vis2_flag
endfor


;---------------------------------------------------------------
;  Next, allow easy flagging of 'bad' photometry.
for t=0,ntel-1 do begin
 for ps=0,num_norms-1 do begin
  flux_tel_flag(*,t*num_norms+ps)=splodges_tel_flag(*,t*dimy)
 endfor
endfor

; mark bad fluxes if the flux_tel is exactly equal to 0.0000
inbad = where(flux_tel eq 0.00,ctbad)
if ctbad ne 0 then begin
 print,ctbad,' values marked bad because of 0.00 0 flux'
 flux_tel_flag[inbad]=1
endif

testfluxplots:

; SETUP
;
syms=[1,2,4,5,6,7] ; Symbols to use for different norm methods.
labels=[norm_labels,'ALL']
yrecords=replicate(!Y,num_norms+1)
cleanplot



;----------
; Loop over each TELESCOPE, presenting all the norm methods at once for easy comparison....

for t=0,ntel-1 do begin

if t eq 0 then telname=mircarray.b1
if t eq 1 then telname=mircarray.b2
if t eq 2 then telname=mircarray.b3
if t eq 3 then telname=mircarray.b4
if t eq 4 then telname=mircarray.b5
if t eq 5 then telname=mircarray.b6


replot_fluxtel:

; Re-calculate medians:
; Normalize each channel by median of the flux for all measurements of the same SOURCE.
uniq_sources = (headers.source)(uniq(headers.source,sort(headers.source)))
uvgoodlist = rawvis2_all.uvin
if n_elements(uvgoodlist) ne n_elements(flux_tel(*,0)) then begin
      print,'problem here.. check this out. '
      stop
    endif


for us=0,n_elements(uniq_sources)-1 do begin

   target_index = where(uvinfo(uvgoodlist).star eq uniq_sources(us))

  for tt=0,ntel-1 do begin
     vals=flux_tel(target_index,(tt*num_norms):(tt*num_norms+num_norms-1))
     val_flags = flux_tel_flag(target_index,(tt*num_norms):(tt*num_norms+num_norms-1))
     goodin=where(val_flags eq 0,ct)
     if ct gt 0 then normvalue = median(vals(goodin)) else normvalue =1
     ; norm by median of all norm methods -- not individually.
   for ps=0,num_norms-1 do begin
     vals = flux_tel(target_index,tt*num_norms+ps)
     val_flags = flux_tel_flag(target_index,tt*num_norms+ps)
     goodin=where(val_flags eq 0,ct)
     if ct gt 0 then normvalue = median(vals(goodin)) else normvalue =1
	if normvalue eq 0 then normvalue=1.0
     flux_tel(target_index,tt*num_norms+ps) = flux_tel(target_index,tt*num_norms+ps) / $
            normvalue
   endfor
  endfor
endfor



erase


xpos = fltarr(num_norms*2,nfiles) ; norm coords of all data points for data selection.
ypos = fltarr(num_norms*2,nfiles)


!p.multi=[0,1,num_norms+1] ; plot CHOP, FIBER, SHUT, ALL combined.
multiplot

 for ps=0,num_norms-1 do begin
  if ps eq 0 then bigtitle='Telescope '+telname else bigtitle=' '
  goodin=where(flux_tel_flag(*,t*num_norms+ps) eq 0,ct)

  if ct gt 0 then begin ; at least one good data point exists
     plot,goodin,flux_tel(goodin,t*num_norms+ps),psym=syms(ps),$
         xr=[0,nfiles-1],ytit=labels(ps),xst=1,$
         xticks=1, xtickname=replicate(' ',2),yminor=1,tit=bigtitle
    temppos=convert_coord(goodin,flux_tel(goodin,t*num_norms+ps),$
         /data,/to_norm)
    xpos(ps,goodin)=reform(temppos(0,*))
    ypos(ps,goodin)=reform(temppos(1,*))

  endif else plot,[0,nfiles-1],[0,0],psym=syms(ps),$
                   xticks=1, xtickname=replicate(' ',2),$
                   xr=[0,nfiles-1],xst=1,yr=[0,1],ytit=labels(ps),yminor=1,tit=bigtitle
  oplot,[0,1e6],[0,0],line=2
  oplot,[0,1e6],[1,1],line=1
  for nm=0,ctmark-1 do begin
   oplot,[markers(nm),markers(nm)],[-1e6,1e6]
  endfor

  yrecords(ps)=!Y
  multiplot
 endfor ;ps
  yrall=[min(yrecords(0:num_norms-1).crange(0)), max(yrecords(0:num_norms-1).crange(1))]

;do the loop again and this time OVERPLOT all together for easy flagging.
 ; mark the distinction between the plots with a thick line.
 oplot,[!x.crange(0),!x.crange(1)],[!y.crange(0),!y.crange(0)],thick=5
 for ps=0,num_norms-1 do begin
   goodin=where(flux_tel_flag(*,t*num_norms+ps) eq 0,ct)
   if ct gt 0 then begin ; at least one good data point exists
    if ps eq 0 then  plot,goodin,flux_tel(goodin,t*num_norms+ps),$
                             psym=syms(ps),xr=[0,nfiles-1],xst=1,yr=yrall,ytit='ALL',$
                             yminor=1,xticks=1, xtickname=replicate(' ',2) $
    else oplot,goodin,flux_tel(goodin,t*num_norms+ps),$
                             psym=syms(ps)
    temppos=convert_coord(goodin,flux_tel(goodin,t*num_norms+ps),$
              /data,/to_norm)
    xpos(num_norms+ps,goodin)=reform(temppos(0,*))
    ypos(num_norms+ps,goodin)=reform(temppos(1,*))
   endif else begin ; no good data for the tel/norm combination
      if ps eq 0 then  plot,[0,nfiles-1],[0,0],psym=syms(ps),xr=[0,nfiles-1],xst=1,$
                            yr=yrall,ytit='ALL',$
                            yminor=1,xticks=1, xtickname=replicate(' ',2) $
                 else oplot,[0,nfiles-1],[0,0],psym=syms(ps)
   endelse

  endfor ;ps
  for nm=0,ctmark-1 do begin
    oplot,[markers(nm),markers(nm)],[-1e6,1e6]
    if (nm mod 4 eq 0) then xyouts,markers(nm),!y.crange(0),'!C'+targets(nm),align=1
    if (nm mod 4 eq 1) then xyouts,markers(nm),!y.crange(0),'!C!C'+targets(nm),align=1
    if (nm mod 4 eq 2) then xyouts,markers(nm),!y.crange(0),'!C!C!C'+targets(nm),align=1
    if (nm mod 4 eq 3) then xyouts,markers(nm),!y.crange(0),'!C!C!C!C'+targets(nm),align=1
   endfor
   oplot,[0,1e6],[0,0],line=2
   oplot,[0,1e6],[1,1],line=1
   yrecords(num_norms)=!y


; Now do the data selection part of the tool,
;   double-click will mark whole block.
;   singel click will find the closest point and remove it.
;   NO RANGE in this tool since the hope is that the fringe selection will be practical
;   for marking extended periods of last fringes/telescopes.

print,' Double click in a square to flag whole block as bad (<0.5 sec)'
print,' Click once to mark a single point as bad'
print,' Click to the right of plot to REPLOT (allows to resize window)'
print,' Click below all plots to reset bad flags!!'
print,' Click to the left of plots to continue..'


if restoreflags eq 1 then goto,skipplot2
dclick = 0.5 ; set time for double click
cursor,x0,y0,/down
cursor,x0,y0,/up
time0=systime(1)
dclickstat=0
while ( (systime(1)-time0) le dclick) do begin
 cursor,x1,y1,/nowait
 if !mouse.button eq 1 then dclickstat = 1
endwhile
wait,.2


if x0 ge 0 then begin ; do something (otherwise continue to next tel).

  if (y0 lt !y.crange(0)) then begin ; reset flags for this telescope...
   for ps=0,num_norms-1 do begin
    flux_tel_flag(*,t*num_norms+ps)=splodges_tel_flag(*,t*dimy)
   endfor
   goto,replot_fluxtel
  endif

  if (x0 gt !x.crange(1)) then goto,replot_fluxtel

  ; If the program gets to here, then one must have clicked
  ; above the  bottom axis and inside the horizontal boundaries of the plot.

  temppos = convert_coord(x0,y0,/data,/to_norm)
  xnorm=temppos(0)
  ynorm=temppos(1)

; Find closest Point
  ri2at,xpos - xnorm,ypos-ynorm,distances,angles
  close_in =(where(distances eq min(distances)))(0) ; closest point
  ac=array_coords(close_in,xpos)

  print,ac
  if ac(0) ge num_norms then norm_in=findgen(num_norms) else norm_in=ac(0)
  if dclickstat eq 1 then begin
    bpick=rawvis2_all(nint(x0)).block
    xin=where(rawvis2_all.block eq bpick)
  endif

  if dclickstat eq 0 then begin
    xin=ac(1)
  endif

  for nm=0,n_elements(norm_in)-1 do begin
    yin=norm_in(nm)
    flux_tel_flag(xin,t*num_norms+yin)=1
  endfor
 ; now replot
 goto,replot_fluxtel

endif ; x0>=0


;window,1
;image_cont,flux_tel_flag,/noc
;wset,0
skipplot2:

endfor ;tel
cleanplot ; clean up the messy multiplot settings.

;stop
; Now.. the last set of flags were potentially specific to the various norm methods
; and thus I need to insert this information into the following loops.
;
; NOTE:
;  I WILL MARK BOTH T3AMP and  T3PHI bad
;  based onthe flux tel flagging.
;  MABYE only need to do T3AMP.. not sure.
skiptemp2:
;----------------------------------------------------

; Apply flags

; Loop through flux calibration for visibility flagging
for psn=0,num_norms-1 do begin ; all types are supported

   rawvis2_all=rawvis2_all_original ; will get reset each time we change normalization
				    ; since each method may have different
				    ; bad files..... unfortunately.
   temp=rawvis2_all

   ; Now ADD the BAD FLAGS for FLUX_TEL_FLAG
   for t=0,ntel-1 do begin
      badin = where(flux_tel_flag(*,t*num_norms+psn) eq 1,ct)
      if ct gt 0 then begin
         bin=[where(ihelp eq t),where(jhelp eq t)]
         bin=bin(where(bin ne -1,ctbase))
         for b=0,ctbase-1 do begin ; should be  ntel-1
            rawvis2_all(badin).badflag(bin(b),*)=1
         endfor
      endif
   endfor

   ; Also must look for -1 values in the calvis2_XXX.vis2
   if norm_labels[psn] eq 'CHOP' then begin
      calvis2_temp = calvis2_chop
   endif
   if norm_labels[psn] eq 'FIBER' then begin
      calvis2_temp = calvis2_fiber
   endif
   if norm_labels[psn] eq 'SHUTTER' then begin
      calvis2_temp = calvis2_shut
   endif
   if norm_labels[psn] eq 'DAQSYNC' then begin
      calvis2_temp = calvis2_daqsync
   endif
   if norm_labels[psn] eq 'XCHAN' then begin
      calvis2_temp = calvis2_xchan
   endif

   ; JDM SAVE BADFLAGS.
   for base=0,nbaselines-1 do begin
      for s=0,dimy-1 do begin

         badin=where(calvis2_temp.vis2[base,s] eq -1,ctbad)
         if ctbad gt 0 then rawvis2_all[badin].badflag[base,s] = 1

         if norm_labels[psn] eq 'CHOP' then begin
            vis2_chop.badflag(base,s)=rawvis2_all.badflag(base,s)
         endif
         if norm_labels[psn] eq 'FIBER' then begin
            vis2_fiber.badflag(base,s)=rawvis2_all.badflag(base,s)
         endif
         if norm_labels[psn] eq 'SHUTTER' then begin
            vis2_shut.badflag(base,s)=rawvis2_all.badflag(base,s)
         endif
         if norm_labels[psn] eq 'DAQSYNC' then begin
            vis2_daqsync.badflag(base,s)=rawvis2_all.badflag(base,s)
         endif
         if norm_labels[psn] eq 'XCHAN' then begin
            vis2_xchan.badflag(base,s)=rawvis2_all.badflag(base,s)
         endif

      endfor  ; loop over spectral channels
   endfor   ; loop over baselines

endfor   ; loop over the NORM methods.


;----------------------------------------------------
; Sigma clipping for visibility flagging

; Sigma clipping is performed on each spectral channel
; (for each baseline and each target block).
; After bad data rejection is finished for each channel,
; the flags are applied globally to all channels (i.e.,
; if a data point is clipped from one channel, then it
; will be clipped from all channels).
; Set skip_edge to specify which channels to include in sigma clipping:
; - skip_edge = 1: exclude edge channels (0 and 7) when applying sigma clipping flags
; - skip_edge = 0: include all 8 channels when applying sigma clipping flags
skip_edge = 1 ; This is set for vis2, T3phi, and T3amp
repeatsigmaclipvis2:
if (fasttrack eq 0) then begin
   ans_clip = ''
;   print,'Would you like to apply sigma clipping to the visibilities?'
   print,'Would you like to apply sigma clipping?'
   print,'(applies to visibilities, closure phases, and closure amplitudes)'
   print,'<Y> or N'
   read,ans_clip
endif else ans_clip = 'N'

if (ans_clip eq 'N') then goto,skipclip_all

sigrej = 3.0
ansrej = 'N'
print,'Sigma clipping rejection threshold (N*stdev): '
print,sigrej
print,'Would you like to change rejection threshold?'
print,'Y of <N>'
read,ansrej
if (ansrej eq 'Y') then begin
   print,'Enter new rejection threshold:'
   read,sigrej
endif
ans_cont = ''

print,'------------------------------'
print,'Sigma clipping of visibilities'
print,'Rejection threshold (N*sigma): ',sigrej
if (skip_edge) then begin
   print,'NOTE: Edge channels will not be used to apply sigma clipping flags'
   print,'Hit enter to continue'
   read, ans_cont
endif

; Loop through flux calibration for sigma clipping of visibilities
for psn=0,num_norms-1 do begin ; all types are supported

   if norm_labels[psn] eq 'CHOP' then begin
      calvis2_temp = calvis2_chop
      rawvis2_all=vis2_chop
   endif
   if norm_labels[psn] eq 'FIBER' then begin
      calvis2_temp = calvis2_fiber
      rawvis2_all=vis2_fiber
   endif
   if norm_labels[psn] eq 'SHUTTER' then begin
      calvis2_temp = calvis2_shut
      rawvis2_all=vis2_shut
   endif
   if norm_labels[psn] eq 'DAQSYNC' then begin
      calvis2_temp = calvis2_daqsync
      rawvis2_all=vis2_daqsync
   endif
   if norm_labels[psn] eq 'XCHAN' then begin
      calvis2_temp = calvis2_xchan
      rawvis2_all=vis2_xchan
   endif

   ut=rawvis2_all.uthrs
   !p.multi=0
   titlabel='('+norm_labels(psn)+' NORM)'

   ; Default flags for sigma clipping
   ;flag_sigclip = intarr(nbaselines,dimy,nfiles)
   flag_sigclip = rawvis2_all.badflag

   ; Loop through baselines for visibility flagging
   for base=0,nbaselines-1 do begin

      inc0=where(rawvis2_all.calib eq 1 and rawvis2_all.badflag(base,spec) eq 0,ct_c)
      int0=where(rawvis2_all.target eq 1 and rawvis2_all.badflag(base,spec) eq 0,ct_t)
      intc0=where( (rawvis2_all.calib eq 1 or rawvis2_all.target eq 1) and $
                   rawvis2_all.badflag(base,spec) eq 0,ct_tc)

      if (ct_tc eq 0) then goto,skipclip0 ; Don't require there aren't fringes on this baseline

      ; Loop through observing blocks and spectral channels, and do sigma clipping
      for bp=0,nblocks-1 do begin

         ; Check if data exists for this target block
         inbp=where(rawvis2_all(intc0).block eq headers(bp).block,ctbp)
         if (ctbp eq 0) then goto,skipclip1 ; No data on this target block

         maxbad = 0.0

         for ns=0, dimy-1 do begin

            nbad_sum = 0.0
            nbad = 1 ; set to go into while loop
            while (nbad gt 0) do begin

               intc0=where( (rawvis2_all.calib eq 1 or rawvis2_all.target eq 1) and $
                            flag_sigclip(base,ns,*) eq 0,ct_tc)

               v2=rawvis2_all.vis2(base,ns) / calvis2_temp.vis2(base,ns)
               v2_err=rawvis2_all.vis2_err(base,ns) / calvis2_temp.vis2(base,ns)

               ; Look for divide by zero and set to 0
               inInf =where(v2 eq 1./0 or v2 eq -1./0 or v2_err eq 1./0 or v2_err eq -1./0 or v2 ne v2 or v2_err ne v2_err,ctInf)
               if ctInf gt 0 then begin
                  v2[inInf]=0.0
                  v2_err[inInf]=1.
               endif

               v2temp = v2(intc0)
               v2temp_err= v2_err(intc0)
               uttemp = ut(intc0)
               blocktemp = rawvis2_all(intc0).block

               inbp=where(blocktemp eq headers(bp).block,ctbp)

               medv2 = median(v2temp(inbp))
               ;stdevv2 = stdev(v2temp(inbp)) ; need to replace with appropriate statisic
               stdevv2=robust_sigma(v2temp(inbp))
               igood = where((v2temp(inbp) le medv2 + sigrej*stdevv2) and (v2temp(inbp) ge medv2 - sigrej*stdevv2) and (v2temp(inbp) ge 0.0), ngood)
               ibad = where((v2temp(inbp) gt medv2 + sigrej*stdevv2) or (v2temp(inbp) lt medv2 - sigrej*stdevv2) or ((v2temp(inbp)+v2temp_err(inbp)) lt 0.0), nbad)
               ;igood2 = where((v2temp(inbp) le medv2 + sigrej*stdevv3) and (v2temp(inbp) ge medv2 - sigrej*stdevv3), ngood2)
               ;ibad2 = where((v2temp(inbp) gt medv2 + sigrej*stdevv3) or (v2temp(inbp) lt medv2 - sigrej*stdevv3), nbad2)
               ;print,'Number of good measurements: ', ngood
               ;print,'Number of bad measurements: ', nbad
               ;print,'Hit enter to continue?'
               ;read,ans_cont

               ;if (nbad gt 0) then rawvis2_all(intc0(inbp(ibad))).badflag(base,ns) = 1
               if (nbad gt 0) then flag_sigclip(base,ns,intc0(inbp(ibad))) = 1

               nbad_sum = nbad_sum + nbad
               if (nbad_sum ge maxbad) then maxbad = nbad_sum

            endwhile   ; while data points are still rejected

            ; Plot good data and sigma clipped data
            ; Need to read flags prior to sigma clipping (use rawvis2_all.badflag)
            intc0=where( (rawvis2_all.calib eq 1 or rawvis2_all.target eq 1) and $
                         rawvis2_all.badflag(base,ns) eq 0,ct_tc)

            inbp=where(rawvis2_all(intc0).block eq headers(bp).block,ctbp)

            if (ns eq 0) then begin
               plotsym,0
               v2spec=rawvis2_all(intc0(inbp)).vis2(base,*) / calvis2_temp(intc0(inbp)).vis2(base,*)

               minvis = min(v2spec)
               maxvis = max(v2spec)
               vrange = maxvis - minvis

               mintime = min(ut(intc0(inbp)))
               maxtime = max(ut(intc0(inbp)))
               trange = maxtime - mintime

               ; Figure out number of rows and columns to plot all spectral channels
               ; in one multiplot:
               nrow = fix(sqrt(dimy))        ; 2 rows for 8 spectral channels
               ncol = fix(round(dimy/nrow))  ; 4 columns for 8 spectral channels

               erase
               multiplot,/reset
               cleanplot
               !P.multi = [0,ncol,nrow]
               multiplot,[ncol,nrow],/init,/doxaxis,/doyaxis,ygap=0.03,xgap=0.02
            endif

            igood = where(flag_sigclip(base,ns,intc0(inbp)) eq 0,ngood)
            ibad = where(flag_sigclip(base,ns,intc0(inbp)) eq 1,nbad)
            ;igood2 = where(flag_sigclip(base,ns,intc0(inbp)) eq 0,ngood2)
            ;ibad2 = where(flag_sigclip(base,ns,intc0(inbp)) eq 1,nbad2)

            multiplot,/doxaxis,/doyaxis,ygap=0.03,xgap=0.02
            ploterror,ut(intc0(inbp(igood))),v2(intc0(inbp(igood))),v2_err(intc0(inbp(igood))),psym=4,$
                      xrange=[mintime-0.1*trange,maxtime+0.1*trange], $
                      yrange=[minvis-0.1*vrange,maxvis+0.1*vrange], $
                      xtit="Time (hrs)",ytit="Visibility**2",xstyle=1,ystyle=1,$
                      title="Baseline Num: "+strcompress(string(fix(base)),/remove_all)+'('+$
                      get_basename(base,mircarray=mircarray)+')'+titlabel+$
                      ' Chan '+ strcompress(string(fix(ns)),/remove_all)
            if (nbad gt 0) then oploterror,ut(intc0(inbp(ibad))),v2(intc0(inbp(ibad))),v2_err(intc0(inbp(ibad))),psym=7,errcolor=200,color=200
            xyouts,mintime-0.05*trange,maxvis+0.05*vrange,headers(bp).source+' ('+strcompress(string(headers(bp).block),/remove_all)+')'
            ;if (nbad2 gt 0) then oploterror,ut(intc0(inbp(ibad2))),v2(intc0(inbp(ibad2))),v2_err(intc0(inbp(ibad2))),psym=7,errcolor=400,color=400
            ;xyouts,mintime-0.05*trange,maxvis+0.05*vrange,headers(bp).source+' ('+strcompress(string(headers(bp).block),/remove_all)+')'

         endfor   ; loop over spectral channels

         ; Count number of rejected measurements (across all channels)
         nflag_sum = 0
         if (skip_edge) then begin ; skip edge channels when counting flags
            for nf=0, ctbp-1 do begin
               bflag = flag_sigclip(base,1:dimy-2,intc0(inbp(nf)))
               indflag = where(bflag eq 1, nflag)
               if (nflag gt 0) then nflag_sum = nflag_sum + 1
            endfor
         endif else begin ; include edge channels when counting flags
            for nf=0, ctbp-1 do begin
               bflag = flag_sigclip(base,*,intc0(inbp(nf)))
               indflag = where(bflag eq 1, nflag)
               if (nflag gt 0) then nflag_sum = nflag_sum + 1
            endfor
         endelse

         print,'Block: ',headers(bp).block,'  Source: ',headers(bp).source
         print,'Diamonds - good data; X - sigma clipped'
         print,'Number of measurements on target (per spectral channel): '
         print, n_elements(inbp)
         print,'Maximum number of sigma clipped measurements per spectral channel: '
         print, maxbad
         print,'Number of uniquely rejected measurements (not including edge channels): '
         print, nflag_sum
         print,'Do you want to reject these measurements?'
         print,'<Y> or N:'
         read,ans_cont

         if (ans_cont eq 'N') then flag_sigclip(base,*,intc0(inbp)) = 0

         skipclip1:

      endfor                    ; loop over observing blocks

      ; Loop through all files on this baseline - set bad flags for all
      ; spectral channels

      if (skip_edge) then begin ; skip edge channels when applying sigma clipping flags
         for nf=0,nfiles-1 do begin
            bflag = flag_sigclip(base,1:dimy-2,nf)
            indflag = where(bflag eq 1, nflag)
            if (nflag gt 0) then rawvis2_all(nf).badflag(base,*) = 1
         endfor
      endif else begin ; include all channels when applying sigma clipping flags
         for nf=0,nfiles-1 do begin
            bflag = flag_sigclip(base,*,nf)
            indflag = where(bflag eq 1, nflag)
            if (nflag gt 0) then rawvis2_all(nf).badflag(base,*) = 1
         endfor
      endelse

      skipclip0:

      print,' Calculating Calibration Table for All Spectral Channels...'
      for s=0,dimy-1 do begin

      ; re-checking badflags.
         intc0=where( (rawvis2_all.calib eq 1 or rawvis2_all.target eq 1) and $
                      rawvis2_all.badflag(base,s) eq 0,ct_tc)

         if norm_labels[psn] eq 'CHOP' then begin
            vis2_chop.badflag(base,s)=rawvis2_all.badflag(base,s)
         endif
         if norm_labels[psn] eq 'FIBER' then begin
            vis2_fiber.badflag(base,s)=rawvis2_all.badflag(base,s)
         endif
         if norm_labels[psn] eq 'SHUTTER' then begin
            vis2_shut.badflag(base,s)=rawvis2_all.badflag(base,s)
         endif
         if norm_labels[psn] eq 'DAQSYNC' then begin
            vis2_daqsync.badflag(base,s)=rawvis2_all.badflag(base,s)
         endif
         if norm_labels[psn] eq 'XCHAN' then begin
            vis2_xchan.badflag(base,s)=rawvis2_all.badflag(base,s)
         endif

      endfor   ; loop over spectral channels.

   endfor    ; baselines

endfor   ; loop over the NORM methods.


goto,repeatsigmaclipvis2


skipclip_all:
print,'----------------------------------------'
print,'Finished sigma clipping the visibilities.'
print,'Begin visual inspection and flagging of visibilities.'
print,'Hit enter to continue'
read,ans_cont

multiplot,/reset
cleanplot
!P.multi=0
;----------------------------------------------------
; Visual data inspection

; Loop through flux calibration for visibility flagging
for psn=0,num_norms-1 do begin ; all types are supported

   if norm_labels[psn] eq 'CHOP' then begin
      calvis2_temp = calvis2_chop
      rawvis2_all=vis2_chop
   endif
   if norm_labels[psn] eq 'FIBER' then begin
      calvis2_temp = calvis2_fiber
      rawvis2_all=vis2_fiber
   endif
   if norm_labels[psn] eq 'SHUTTER' then begin
      calvis2_temp = calvis2_shut
      rawvis2_all=vis2_shut
   endif
   if norm_labels[psn] eq 'DAQSYNC' then begin
      calvis2_temp = calvis2_daqsync
      rawvis2_all=vis2_daqsync
   endif
   if norm_labels[psn] eq 'XCHAN' then begin
      calvis2_temp = calvis2_xchan
      rawvis2_all=vis2_xchan
   endif

   ; Loop through baselines for visibility flagging
   ;skipstart:
   for base=0,nbaselines-1 do begin

      zoomup=0
      ;skipsetup:
      ;smoothlength = float(dialog_input(prompt="Smoothing (1/3) Length (hrs):",initial=smoothlength))
      skipsetup1:

      inc0=where(rawvis2_all.calib eq 1 and rawvis2_all.badflag(base,spec) eq 0,ct_c)
      int0=where(rawvis2_all.target eq 1 and rawvis2_all.badflag(base,spec) eq 0,ct_t)
      intc0=where( (rawvis2_all.calib eq 1 or rawvis2_all.target eq 1) and $
                   rawvis2_all.badflag(base,spec) eq 0,ct_tc)

      if ct_tc eq 0 then goto,skipthis1 ; Don't require there e fringes in all baselines.

      mintimes=min(rawvis2_all(intc0).uthrs)
      maxtimes=max(rawvis2_all(intc0).uthrs)
      num_interp = fix((maxtimes-mintimes+.4)/(smoothlength/5.) ) > 20
      uthrs=pvector([mintimes-.2,maxtimes+.2,num_interp])
      numsrc=n_elements(int0)


      !p.multi=0
      titlabel='('+norm_labels(psn)+' NORM)'

      if norm_labels[psn] eq 'CHOP' then begin

      ; result_fun = caltables_psnorm(rawvis2_all(intc0).uthrs,rawvis2_all(intc0).uthrs,$
      ;    rawvis2_all(intc0).vis2(base,spec),rawvis2_all(intc0).vis2_err(base,spec),$
      ;    calvis2_chop(intc0).vis2(base,spec),smooth_fun,booterror=booterror_fun,$
      ;    nsamps = 5)

         v2=rawvis2_all.vis2(base,spec) / calvis2_chop.vis2(base,spec)
         v2_err=rawvis2_all.vis2_err(base,spec) / calvis2_chop.vis2(base,spec)
      endif

      if norm_labels[psn] eq 'FIBER' then begin

      ; result_fun = caltables_psnorm(rawvis2_all(intc0).uthrs,rawvis2_all(intc0).uthrs,$
      ;    rawvis2_all(intc0).vis2(base,spec),rawvis2_all(intc0).vis2_err(base,spec),$
      ;    calvis2_fiber(intc0).vis2(base,spec),smooth_fun,booterror=booterror_fun,nsamps=5)

         v2=rawvis2_all.vis2(base,spec) / calvis2_fiber.vis2(base,spec)
         v2_err=rawvis2_all.vis2_err(base,spec) / calvis2_fiber.vis2(base,spec)
      endif

      if norm_labels[psn] eq 'SHUTTER' then begin

      ;  result_fun = caltables_psnorm(rawvis2_all(intc0).uthrs,rawvis2_all(intc0).uthrs,$
      ;    rawvis2_all(intc0).vis2(base,spec),rawvis2_all(intc0).vis2_err(base,spec),$
      ;    calvis2_shut(intc0).vis2(base,spec),smooth_fun,booterror=booterror_fun,nsamps=5)

         v2=rawvis2_all.vis2(base,spec) / calvis2_shut.vis2(base,spec)
         v2_err=rawvis2_all.vis2_err(base,spec) / calvis2_shut.vis2(base,spec)
      endif

      if norm_labels[psn] eq 'DAQSYNC' then begin

      ;  result_fun = caltables_psnorm(rawvis2_all(intc0).uthrs,rawvis2_all(intc0).uthrs,$
      ;    rawvis2_all(intc0).vis2(base,spec),rawvis2_all(intc0).vis2_err(base,spec),$
      ;    calvis2_daqsync(intc0).vis2(base,spec),smooth_fun,booterror=booterror_fun,nsamps=5)

         v2=rawvis2_all.vis2(base,spec) / calvis2_daqsync.vis2(base,spec)
         v2_err=rawvis2_all.vis2_err(base,spec) / calvis2_daqsync.vis2(base,spec)
      endif
      if norm_labels[psn] eq 'XCHAN' then begin
      ;  result_fun = caltables_psnorm(rawvis2_all(intc0).uthrs,rawvis2_all(intc0).uthrs,$
      ;    rawvis2_all(intc0).vis2(base,spec),rawvis2_all(intc0).vis2_err(base,spec),$
      ;    calvis2_xchan(intc0).vis2(base,spec),smooth_fun,booterror=booterror_fun,nsamps=5)

         v2=rawvis2_all.vis2(base,spec) / calvis2_xchan.vis2(base,spec)
         v2_err=rawvis2_all.vis2_err(base,spec) / calvis2_xchan.vis2(base,spec)
      endif


      ; Look for divide by zero and set to 0
      inInf =where(v2 eq 1./0 or v2 eq -1./0 or v2_err eq 1./0 or v2_err eq -1./0 or v2 ne v2 or v2_err ne v2_err,ctInf)
      if ctInf gt 0 then begin
         v2[inInf]=0.0
         v2_err[inInf]=1.
      endif



      plotsym,0

      ut=rawvis2_all.uthrs

      if (zoomup eq 0) then begin ; standard scaling.
         yr=[  min([0,v2(intc0)]),max(v2(intc0))]
         xr=[mintimes-.2,maxtimes+.2]
      endif else begin
         yr=zoomyr
         xr=zoomxr              ;defined elsewheres
      endelse


      ploterror,ut(intc0),v2(intc0),v2_err(intc0), psym=4,$
                xtit="Time (hrs)",ytit="Visibility**2",xr=xr,xst=1,yr=yr,$
                title="Baseline Num:"+string(fix(base))+'('+$
                get_basename(base,mircarray=mircarray)+')'+titlabel
      ;plotsym,0,/fill
      ;oploterror,ut(intc0),v2(intc0),v2_err(intc0),psym=8

      oplot,[-1e6,1e6],[0,0],line=1

      ; Draw all boundaries between blocks and label.
      for bp=0,nblocks-1 do begin
         inbp=where(rawvis2_all.block eq headers(bp).block,ctbp)
         bptimes=rawvis2_all(inbp).uthrs
         oplot,[min(bptimes),min(bptimes)],!y.crange,line=2
         ;  oplot,[max(bptimes),max(bptimes)],!y.crange,line=2
         xyouts,min(bptimes),!y.crange(1),'!C'+targets(bp),orient=-90
      endfor


      if (zoomup eq 1) then legend,['ZOOM UP'],box=0


      ;oplot,ut(intc0),result_fun,thick=2
      ;oplot,ut(intc0),result_fun+booterror_fun,thick=.5
      ;oplot,ut(intc0),result_fun-booterror_fun,thick=.5

      print,' '
      print,' Click left of axis to continue (or to goback to unzoomed view)'
      print,' Click right of axis to change smoothing length'
      print,' Click near datapoint to remove a time sample.'
      print,' Click below axis twice to specify a range to remove.'
      print,' Click above plot twice to specify a range to ZOOM up'
      ;stop
      if noclick eq 1 or fasttrack eq 1 then x0=!x.crange(0)-1 else $
         cursor,x0,y0
      wait,.2
      if (x0 ge !x.crange(1)) then begin
         smoothlength = float(dialog_input(prompt="Smoothing (1/3) Length (hrs):",initial=smoothlength))
         goto,skipsetup1
      endif

      if (x0 le !x.crange(1) and x0 ge !x.crange(0)) then begin

      ;  identify bad thing:
         if (y0 ge !y.crange(0) and y0 lt !y.crange(1)) then begin

            temp=convert_coord(x0,y0 ,/data,/to_device)
            temp1=convert_coord(ut[intc0],v2[intc0],/data,/to_device)

      ; ri2at,(x0-ut(intc0))/(xr[1]-xr[0]), (y0-v2(intc0))/(yr(1)-yr(0)), tdiff,angle
            ri2at, (temp[0])-reform(temp1[0,*]) , (temp[1])-reform(temp1[1,*]),$
                   tdiff,angle


            in0=where(tdiff eq min(tdiff))
            rawvis2_all(intc0(in0)).badflag(base,*)=1 ; mark all waves
         endif
         if (y0 lt !y.crange(0)) then begin
            print,'Click 2nd point for range'
            cursor,x1,y1
            wait,.2
            xmin=min([x0,x1])
            xmax=max([x0,x1])

            in0=where(  ut(intc0) ge xmin and ut(intc0) le xmax,ctbad )
            for cb=0,ctbad-1 do begin
               rawvis2_all(intc0(in0(cb))).badflag(base,*)=1 ; mark all waves
            endfor

         endif                      ;kill-range
         if (y0 ge !y.crange(1)) then begin ; new zoomwindow
            print,'Click 2nd point for zoom range'
            cursor,x1,y1
            wait,.2
            xmin=min([x0,x1])
            xmax=max([x0,x1])
            zoomxr=[xmin,xmax]

            in0=where(  ut(intc0) ge xmin and ut(intc0) le xmax,ctbad )
            if ctbad gt 0 then begin
               zoomyr=[  min([0,v2(intc0(in0))]),max(v2(intc0(in0)))]
               zoomup=1
            endif

         endif                  ; zoom kill

         goto,skipsetup1
      endif
      if (zoomup eq 1) then begin
         zoomup=0
         print,' Back to normal view of all data'
         goto,skipsetup1        ; jdm should rewrite this using routines... or gui.
      endif


      ; where do I save this calibrated info?
      ; oh right -- I created some arrays called  vis2_chop, etc. use them!!
      ;----Save calibration factors and the badflag arrays but do not apply yet. This will be done
      ; in last step after correlation analysis.

      skipthis1:
      ; This is a good time to deal with the correlated errors from the calibrator....but how? -jdm
      print,' Calculating Calibration Table for All Spectral Channels...'
      for s=0,dimy-1 do begin

      ; re-checking badflags.
         intc0=where( (rawvis2_all.calib eq 1 or rawvis2_all.target eq 1) and $
                      rawvis2_all.badflag(base,s) eq 0,ct_tc)

         if norm_labels[psn] eq 'CHOP' then begin
            vis2_chop.badflag(base,s)=rawvis2_all.badflag(base,s)
         endif
         if norm_labels[psn] eq 'FIBER' then begin
            vis2_fiber.badflag(base,s)=rawvis2_all.badflag(base,s)
         endif
         if norm_labels[psn] eq 'SHUTTER' then begin
            vis2_shut.badflag(base,s)=rawvis2_all.badflag(base,s)
         endif
         if norm_labels[psn] eq 'DAQSYNC' then begin
            vis2_daqsync.badflag(base,s)=rawvis2_all.badflag(base,s)
         endif
         if norm_labels[psn] eq 'XCHAN' then begin
            vis2_xchan.badflag(base,s)=rawvis2_all.badflag(base,s)
         endif

      endfor   ; loop over spectral channels.

   endfor    ; baselines
   skipclicks1:
endfor   ; loop over the NORM methods.

;---------------------------------------------------------------------------
; **********
; Now do T3DATA
;***********
; Now I will go through and strip data out of the array of pointers.
; This is a pain but probably necessary.
t3_template = {uthrs:0.0, t3amp:fltarr(ntriangles,dimy), $
  t3amp_err:fltarr(ntriangles,dimy), diam_t3amp_err:fltarr(ntriangles,dimy),diam_cal:'', t3phi:fltarr(ntriangles,dimy), $
  t3phi_err:fltarr(ntriangles,dimy),$
   uvin:0,$
   block:0, file:0, target:0, calib:0,$
   badflag_t3amp:intarr(ntriangles,dimy), badflag_t3phi:intarr(ntriangles,dimy)}
rawt3_all=replicate(t3_template,nfiles)
original_rawt3_all=rawt3_all
calt3_shut=replicate(t3_template,nfiles)
calt3_chop=replicate(t3_template,nfiles)
calt3_fiber=replicate(t3_template,nfiles)
calt3_daqsync=replicate(t3_template,nfiles)
calt3_xchan=replicate(t3_template,nfiles)


c=0
for b=0,nblocks-1 do begin
   nfiles0=n_elements( *misc(b).num_data)
   rawt3amp=(*tpinfo(b).tp)
   rawt3amp_err = (*tpinfo(b).tp_err)
   original_rawt3amp=(*original_tpinfo(b).tp)
   original_rawt3amp_err = (*original_tpinfo(b).tp_err)

   rawt3amp_err_cal = (*tpinfo_cal(b).tp_err)
   rawt3phi=(*tpinfo(b).cp)
   rawt3phi_err = (*tpinfo(b).cp_err) ; Assume tp bias was already corrected
  				      ; in previous steps, if necessary

   datafiles=*misc(b).num_data

   for f=0,nfiles0-1 do begin
      uvin=where(uvinfo.star eq headers(b).source and uvinfo.block eq headers(b).block and uvinfo.filenum eq datafiles(f)) ;crosscheck
      rawt3_all(c).uthrs=uvinfo(uvin).time/3600. ; assumes all data in night on same utdate.
      rawt3_all(c).t3amp=rawt3amp(*,*,f)
      rawt3_all(c).t3amp_err=rawt3amp_err(*,*,f)
      original_rawt3_all(c).t3amp=original_rawt3amp(*,*,f)
      original_rawt3_all(c).t3amp_err=original_rawt3amp_err(*,*,f)

      rawt3_all(c).diam_t3amp_err = rawt3amp_err_cal(*,*,f)
      rawt3_all(c).diam_cal=headers(b).source
      rawt3_all(c).t3phi=rawt3phi(*,*,f)
      rawt3_all(c).t3phi_err=rawt3phi_err(*,*,f)

      rawt3_all(c).uvin=uvin(0)
      rawt3_all(c).block=headers(b).block
      rawt3_all(c).file=datafiles(f)
      rawt3_all(c).target=headers(b).target
      rawt3_all(c).calib=headers(b).calib

      for t=0,ntriangles-1 do begin ; setup t3amp calibration.

         ;must do this for shut/fiber/chop more ! pain!
         calt3_shut(c).t3amp(t,*)= sqrt(  $
              (*vis2info(b).vis2_norms_shut)(itribase(t),*,f) * $
              (*vis2info(b).vis2_norms_shut)(jtribase(t),*,f) * $
              (*vis2info(b).vis2_norms_shut)(ktribase(t),*,f)  )

         if keyword_set(daqsync) eq 1 then $
            calt3_daqsync(c).t3amp(t,*)= sqrt(  $
                 (*vis2info(b).vis2_norms_daqsync)(itribase(t),*,f) * $
                 (*vis2info(b).vis2_norms_daqsync)(jtribase(t),*,f) * $
                 (*vis2info(b).vis2_norms_daqsync)(ktribase(t),*,f)  )

         if keyword_set(xchan) eq 1 then $
            calt3_xchan(c).t3amp(t,*)= sqrt(  $
                 (*vis2info(b).vis2_norms_xchan)(itribase(t),*,f) * $
                 (*vis2info(b).vis2_norms_xchan)(jtribase(t),*,f) * $
                 (*vis2info(b).vis2_norms_xchan)(ktribase(t),*,f)  )

         calt3_fiber(c).t3amp(t,*)= sqrt(  $
              (*vis2info(b).vis2_norms_fiber)(itribase(t),*,f) * $
              (*vis2info(b).vis2_norms_fiber)(jtribase(t),*,f) * $
              (*vis2info(b).vis2_norms_fiber)(ktribase(t),*,f)  )

         if keyword_set(chopon) eq 1 then $
            calt3_chop(c).t3amp(t,*)= sqrt(  $
                 (*vis2info(b).vis2_norms_chop)(itribase(t),*,f) * $
                 (*vis2info(b).vis2_norms_chop)(jtribase(t),*,f) * $
                 (*vis2info(b).vis2_norms_chop)(ktribase(t),*,f)  )
      endfor                    ; Filling in all the norm factors for the triangles.

      c=c+1
   endfor                       ; f
endfor                          ;blocks

;Check for NANs.... should be removed elsewhere -- print WARNING!
if inshut ne -1 then begin
   in=where(calt3_shut.t3amp ne calt3_shut.t3amp or $
            calt3_shut.t3amp eq 0,ct)
   if ct gt 0 then begin
      print,' NANs  or 00s in calt3_shut!: ',ct
      ac=array_coords(in, (calt3_shut.t3amp) )
      calt3_shut(ac(*,2)).t3amp(ac(*,0),ac(*,1)) = -1
   endif
endif

if infiber ne -1 then begin
   in=where(calt3_fiber.t3amp ne calt3_fiber.t3amp or $
            calt3_fiber.t3amp eq 0,ct)
   if ct gt 0 then begin
      print,' NANs or 00s in calt3_fiber!: ',ct
      ac=array_coords(in, (calt3_fiber.t3amp) )
      calt3_fiber(ac(*,2)).t3amp(ac(*,0),ac(*,1)) = -1
   endif
endif

if inchop ne -1 then begin
   in=where(calt3_chop.t3amp ne calt3_chop.t3amp or $
            calt3_chop.t3amp eq 0.0 ,ct)
   if ct gt 0 then begin
      print,' NANs or 0s in calt3_chop!: ',ct
      ac=array_coords(in, (calt3_chop.t3amp) )
      calt3_chop(ac(*,2)).t3amp(ac(*,0),ac(*,1)) = -1
   endif
endif

if indaqsync ne -1 then begin
   in=where(calt3_daqsync.t3amp ne calt3_daqsync.t3amp or $
            calt3_daqsync.t3amp eq 0.0 ,ct)
   if ct gt 0 then begin
      print,' NANs  000s. in calt3_daqsync!: ',ct
      ac=array_coords(in, (calt3_daqsync.t3amp) )
      calt3_daqsync(ac(*,2)).t3amp(ac(*,0),ac(*,1)) = -1
   endif
endif

if inxchan ne -1 then begin
   in=where(calt3_xchan.t3amp ne calt3_xchan.t3amp or $
            calt3_xchan.t3amp eq 0.0  ,ct)
   if ct gt 0 then begin
      print,' NANs 000s. in calt3_xchan!: ',ct
      ; Faster way
      ;ac=array_coords(in, (calt3_xchan.t3amp) )
      test = calt3_xchan.t3amp
      test[in]=-1
      calt3_xchan.t3amp=test

;     test_ac = array_indices(calt3_xchan.t3amp,in)
;     calt3_xchan(ac(*,2)).t3amp(ac(*,0),ac(*,1)) = -1

   endif
endif

;-----------------------
; Apply splodge flags to the
; vis2 data.
;
rawt3_all_original=rawt3_all
for s=0,dimy-1 do begin
   rawt3_all_original.badflag_t3phi(*,s) = splodges_t3_flag
   rawt3_all_original.badflag_t3amp(*,s) = splodges_t3_flag
endfor


tempskip2:

;--------------------------------------------------------
; Do t3amp and t3phi separately.
;--------------------------------------------------------

;-------------------------------------------------------
; next up the closure phase.
;---------------------------------------------------------
rawt3_all=rawt3_all_original

inc=where(rawt3_all.calib eq 1)
int=where(rawt3_all.target eq 1)
intc=where(rawt3_all.calib eq 1 or rawt3_all.target eq 1)

t3_chop=rawt3_all ;init
t3_shut=rawt3_all
t3_fiber=rawt3_all
t3_daqsync=rawt3_all
t3_xchan=rawt3_all

;spec=fix(dimy/2)-1. ; pick middle of spectral channel for doing estimates.
; defined above
smoothlength=3


; Now ADD the BAD FLAGS for FLUX_TEL_FLAG
; NOTE: Only mark t3phi as BAD if ALL of the t3amp (all NORM methods)
; are also flagged.... [this is tricky]
; This whole marking could be optional.

for t=0,ntel-1 do begin
   flags = flux_tel_flag(*,(t*num_norms):(t*num_norms+num_norms-1))
   if num_norms gt 1 then flags = total(flags,2)
   badin =where(flags eq num_norms,ct) ; telflux marked bad for all norm methods!
   if ct gt 0 then begin
      tin=[where(itrihelp eq t),where(jtrihelp eq t), where(ktrihelp eq t)]
      tin=tin(where(tin ne -1,cttri))
      for t2=0,cttri-1 do begin
         rawt3_all(badin).badflag_t3amp(tin(t2),*)=1
         rawt3_all(badin).badflag_t3phi(tin(t2),*)=1
      endfor
   endif
endfor
jdmskip2:


;----------------------------------------------------
; Sigma clipping of closure phases

repeatsigmaclipcp:
if (fasttrack eq 0) then begin
   ans_clip = '
   print,'Would you like to apply sigma clipping to closure phases?'
   print,'<Y> or N'
   read,ans_clip
endif else ans_clip = 'N'

if (ans_clip eq 'N') then goto,skipclip_cp_all

sigrej = 3.0
ansrej = 'N'
print,'Sigma clipping rejection threshold: '
print,sigrej
print,'Would you like to change rejection threshold?'
print,'Y of <N>'
read,ansrej
if (ansrej eq 'Y') then begin
   print,'Enter new rejection threshold:'
   read,sigrej
endif
ans_cont = ''

print,'--------------------------------'
print,'Sigma clipping of closure phases'
print,'Rejection threshold (N*sigma): ',sigrej
if (skip_edge) then begin
   print,'NOTE: Edge channels will not be used to apply sigma clipping flags'
   print,'Hit enter to continue'
   read, ans_cont
endif

ut=rawt3_all.uthrs
flag_sigclip = rawt3_all.badflag_t3phi

for tri=0,ntriangles-1 do begin

   inc0=where(rawt3_all.calib eq 1 and rawt3_all.badflag_t3phi(tri,spec) eq 0)
   int0=where(rawt3_all.target eq 1 and rawt3_all.badflag_t3phi(tri,spec) eq 0)
   intc0=where( (rawt3_all.calib eq 1 or rawt3_all.target eq 1) and $
                rawt3_all.badflag_t3phi(tri,spec) eq 0,ct_tc)

   ; Check it there are data on this triangle
   if ct_tc eq 0 then goto,skipclip_cp0

   ; Loop through observing blocks and spectral channels, and do sigma clipping
   for bp=0,nblocks-1 do begin

      ; Check if data exists for this target block
      inbp=where(rawt3_all(intc0).block eq headers(bp).block,ctbp)
      if (ctbp eq 0) then goto,skipclip_cp1 ; No data on this target block

      maxbad = 0.0

      for ns=0, dimy-1 do begin

         nbad_sum = 0.0
         nbad = 1               ; set to go into while loop
         while (nbad gt 0) do begin

            ; Extract data for this triangle and spectral channel
            intc0=where( (rawt3_all.calib eq 1 or rawt3_all.target eq 1) and $
                         flag_sigclip(tri,ns,*) eq 0,ct_tc)

            t3phi=rawt3_all.t3phi(tri,ns)
            t3phi_err=rawt3_all.t3phi_err(tri,ns)

            t3phi_temp = t3phi(intc0)
            t3phi_err_temp = t3phi_err(intc0)
            uttemp = ut(intc0)
            blocktemp = rawt3_all(intc0).block

            ; Extract data block for this spectral channel
            inbp=where(blocktemp eq headers(bp).block,ctbp)
            ;
            ; medt3=atan(total(sin(t3phi_temp(inbp))),total(cos(t3phi_temp(inbp))))
            ;errweight=t3phi_err_temp(inbp)
            ;result = boot_angles_mod(t3phi_temp(inbp),weights=1.)
            ;medt3=result[0]
            ;medt3 = median(t3phi_temp(inbp))
            angle=mod360(t3phi_temp(inbp))
            ;print,"Angle: ",angle
            rp0=total(cos(!pi*angle/180.))
            ip0=total(sin(!pi*angle/180.))
            medt3=180.*atan(ip0,rp0)/!pi
            print,"Mean: ", medt3

            if n_elements(t3phi_temp(inbp)) gt 2 then begin
              stdevhold=0.0
              for i=0,n_elements(angle)-1 do begin
                stdevhold=stdevhold+angle_diff(medt3,angle(i))^2.0
              endfor
              stdevt3=sqrt(stdevhold/n_elements(angle))
              print,"stdevt3: ",stdevt3
            ;  stdevt3 = robust_sigma(angle)
              ;stdevt3 = mod360((angle_stats(t3phi_temp))(1))
            ;  print,stdevt3
            endif else begin
              print,"Number of data points 1 or 0! Skipping this channel."
              goto,notenoughpoints1
            endelse
            ;stdevt3=robust_sigma(t3phi_temp(inbp))
            ;igood = where((mod360(t3phi_temp(inbp)) le mod360(medt3 + sigrej*stdevt3)) and (mod360(t3phi_temp(inbp)) ge mod360(medt3 - sigrej*stdevt3)), ngood)
            ;ibad = where((mod360(t3phi_temp(inbp)) gt mod360(medt3 + sigrej*stdevt3)) or (mod360(t3phi_temp(inbp)) lt mod360(medt3 - sigrej*stdevt3)), nbad)
            diff=abs(angle_diff(t3phi_temp(inbp),medt3))
            print,"diff: ",diff
            igood = where( diff le sigrej*stdevt3, ngood)
            ibad = where(diff gt sigrej*stdevt3, nbad)
            if (nbad gt 0) then flag_sigclip(tri,ns,intc0(inbp(ibad))) = 1
            nbad_sum = nbad_sum + nbad
            if (nbad_sum ge maxbad) then maxbad = nbad_sum

         endwhile   ; while data points are still rejected

         ; Plot good data and sigma clipped data
         ; Need to read flags prior to sigma clipping (use rawt3_all.badflag)
         intc0=where( (rawt3_all.calib eq 1 or rawt3_all.target eq 1) and $
                      rawt3_all.badflag_t3phi(tri,spec) eq 0,ct_tc)

         ; Extract data block for this spectral channel
         inbp=where(rawt3_all(intc0).block eq headers(bp).block,ctbp)

         if (ns eq 0) then begin
            plotsym,0
            t3phi_spec=rawt3_all(intc0(inbp)).t3phi(tri,*)

            mint3 = min(t3phi_spec)
            maxt3 = max(t3phi_spec)
            midt3 = (maxt3 + mint3)/2.0
            t3range = maxt3 - mint3

            mintime = min(ut(intc0(inbp)))
            maxtime = max(ut(intc0(inbp)))
            utrange = maxtime - mintime

            ; Figure out number of rows and columns to plot all spectral channels
            ; in one multiplot:
            nrow = fix(sqrt(dimy))        ; 2 rows for 8 spectral channels
            ncol = fix(round(dimy/nrow))  ; 4 columns for 8 spectral channels

            erase
            multiplot,/reset
            cleanplot
            !P.multi = [0,ncol,nrow]
            multiplot,[ncol,nrow],/init,/doxaxis,/doyaxis,ygap=0.03,xgap=0.02
         endif

         igood = where(flag_sigclip(tri,ns,intc0(inbp)) eq 0,ngood)
         ibad = where(flag_sigclip(tri,ns,intc0(inbp)) eq 1,nbad)

         multiplot,/doxaxis,/doyaxis,ygap=0.03,xgap=0.02
         if (n_elements(t3phi(intc0(inbp(igood)))) gt 2) then begin
         ploterror,ut(intc0(inbp(igood))),t3phi(intc0(inbp(igood))),t3phi_err(intc0(inbp(igood))), psym=4, $
                   xrange=[mintime-0.1*trange,maxtime+0.1*utrange], $
                   yrange=[mint3-0.1*t3range,maxt3+0.1*t3range], $
                   xtit="Time (hrs)",ytit="Closure Phase (degs)",xstyle=1,ystyle=1,$
                   title="Triangle Num: "+strcompress(string(fix(tri)),/remove_all)+$
                   '('+get_triname(tri,mircarray=mircarray)+')'+$
                   ' Chan '+ strcompress(string(fix(ns)),/remove_all)

         if (nbad gt 0) then oploterror,ut(intc0(inbp(ibad))),t3phi(intc0(inbp(ibad))),t3phi_err(intc0(inbp(ibad))), psym=7,errcolor=200,color=200

         xyouts,mintime-0.05*utrange,maxt3+0.05*t3range,headers(bp).source+' ('+strcompress(string(headers(bp).block),/remove_all)+')'
         print,"Bad values: ", t3phi(intc0(inbp(ibad)))
       endif else begin
        goto,notenoughpoints1
      endelse

         notenoughpoints1:

      endfor                    ; loop over spectral channels

      ; Count number of rejected measurements (across all channels)
      nflag_sum = 0
      if (skip_edge) then begin ; skip edge channels when counting flags
         for nf=0, ctbp-1 do begin
            bflag = flag_sigclip(tri,1:dimy-2,intc0(inbp(nf)))
            indflag = where(bflag eq 1, nflag)
            if (nflag gt 0) then nflag_sum = nflag_sum + 1
         endfor
      endif else begin          ; include edge channels when counting flags
         for nf=0, ctbp-1 do begin
            bflag = flag_sigclip(tri,*,intc0(inbp(nf)))
            indflag = where(bflag eq 1, nflag)
            if (nflag gt 0) then nflag_sum = nflag_sum + 1
         endfor
      endelse

      print,'Block: ',headers(bp).block,'  Source: ',headers(bp).source
      print,'Diamonds - good data; X - sigma clipped'
      print,'Number of measurements on target (per spectral channel): '
      print, n_elements(inbp)
      print,'Maximum number of sigma clipped measurements per spectral channel: '
      print, maxbad
      print,'Number of uniquely rejected measurements (not including edge channels): '
      print, nflag_sum
      print,'Do you want to reject these measurements?'
      print,'<Y> or N:'
      read,ans_cont

      if (ans_cont eq 'N') then flag_sigclip(tri,*,intc0(inbp)) = 0

      skipclip_cp1:

   endfor                       ; loop over observing blocks

   ; Loop through all files on this triangle - set bad flags for all
   ; spectral channels

   if (skip_edge) then begin   ; skip edge channels when applying sigma clipping flags
      for nf=0,nfiles-1 do begin
         bflag = flag_sigclip(tri,1:dimy-2,nf)
         indflag = where(bflag eq 1, nflag)
         if (nflag gt 0) then rawt3_all(nf).badflag_t3phi(tri,*) = 1
      endfor
   endif else begin    ; include all channels when applying sigma clipping flags
      for nf=0,nfiles-1 do begin
         bflag = flag_sigclip(tri,*,nf)
         indflag = where(bflag eq 1, nflag)
         if (nflag gt 0) then rawt3_all(nf).badflag_t3phi(tri,*) = 1
      endfor
   endelse

   skipclip_cp0:

   ; where do I save this calibrated info?
   ;t3_chop.. etc.

   ; norm factors don't matter for this calibration. EXCEPT they have a different set of badflags..
   ; Now, because the closure phase is *required* for averaging of the
   ; triple amp, I will flag all the t3amp to be consistent with t3phi.
   print,' Calculating Calibration Table for All Spectral Channels...'
   for s=0,dimy-1 do begin

      inc0=where(rawt3_all.calib eq 1 and rawt3_all.badflag_t3phi(tri,s) eq 0)
      int0=where(rawt3_all.target eq 1 and rawt3_all.badflag_t3phi(tri,s) eq 0)
      intc0=where( (rawt3_all.calib eq 1 or rawt3_all.target eq 1) and $
                   rawt3_all.badflag_t3phi(tri,s) eq 0)


      t3_chop.badflag_t3amp(tri,s)=rawt3_all.badflag_t3phi(tri,s)
      t3_chop.badflag_t3phi(tri,s)=rawt3_all.badflag_t3phi(tri,s)
      t3_fiber.badflag_t3amp(tri,s)=rawt3_all.badflag_t3phi(tri,s)
      t3_fiber.badflag_t3phi(tri,s)=rawt3_all.badflag_t3phi(tri,s)
      t3_shut.badflag_t3amp(tri,s)=rawt3_all.badflag_t3phi(tri,s)
      t3_shut.badflag_t3phi(tri,s)=rawt3_all.badflag_t3phi(tri,s)
      t3_daqsync.badflag_t3amp(tri,s)=rawt3_all.badflag_t3phi(tri,s)
      t3_daqsync.badflag_t3phi(tri,s)=rawt3_all.badflag_t3phi(tri,s)
      t3_xchan.badflag_t3amp(tri,s)=rawt3_all.badflag_t3phi(tri,s)
      t3_xchan.badflag_t3phi(tri,s)=rawt3_all.badflag_t3phi(tri,s)

   endfor    ; spec channels

endfor    ; triangles



goto,repeatsigmaclipcp
skipclip_cp_all:
print,'------------------------------------------'
print,'Finished sigma clipping the closure phases.'
print,'Begin visual inspection and flagging of closure phases.'
print,'Hit enter to continue'
read,ans_cont

multiplot,/reset
cleanplot
!P.multi=0
;----------------------------------------------------
; Visual flagging closure phases
for tri=0,ntriangles-1 do begin

   zoomup=0
   skipsetup_t3phi:
   ;smoothlength = float(dialog_input(prompt="Smoothing (1/3) Length (hrs):",initial=smoothlength))
   skipsetup1_t3phi:

   inc0=where(rawt3_all.calib eq 1 and rawt3_all.badflag_t3phi(tri,spec) eq 0)
   int0=where(rawt3_all.target eq 1 and rawt3_all.badflag_t3phi(tri,spec) eq 0)
   intc0=where( (rawt3_all.calib eq 1 or rawt3_all.target eq 1) and $
                rawt3_all.badflag_t3phi(tri,spec) eq 0,ct_tc)

   if ct_tc eq 0 then goto,skipthis2

   mintimes=min(rawt3_all(intc0).uthrs)
   maxtimes=max(rawt3_all(intc0).uthrs)
   num_interp = fix((maxtimes-mintimes+.4)/(smoothlength/5.) ) > 20
   uthrs=pvector([mintimes-.2,maxtimes+.2,num_interp])
   numsrc=n_elements(int0)

   !p.multi=0

   ; no difference between norm methods..

   ;result_fun = caltables(uthrs,rawt3_all(intc0).uthrs,$
   ;     rawt3_all(intc0).t3phi(tri,spec),rawt3_all(intc0).t3phi_err(tri,spec),$
   ;     smooth_fun,booterror=booterror_fun,/phase,nsamps=5)


   plotsym,0,/fill
   ut=rawt3_all.uthrs
   t3phi=rawt3_all.t3phi(tri,spec)
   t3phi_err=rawt3_all.t3phi_err(tri,spec)

   if (zoomup eq 0) then begin  ; standard scaling.
      yr=[  min([0,t3phi(intc0)]),max(t3phi(intc0))]
      xr=[mintimes-.2,maxtimes+.2]
   endif else begin
      yr=zoomyr
      xr=zoomxr                 ;defined elsewheres
   endelse
   ploterror,ut(intc0),t3phi(intc0),t3phi_err(intc0), psym=4,$
             xtit="Time (hrs)",ytit="Closure Phase (degs)",xr=xr,xst=1,yr=yr,$
             title="Triangle Num:"+string(fix(tri))+'('+get_triname(tri,mircarray=mircarray)+')'

   plotsym,0,/fill
   oplot,[-1e6,1e6],[0,0],line=1

   ; Draw all boundaries between blocks and label.
   for bp=0,nblocks-1 do begin
      inbp=where(rawvis2_all.block eq headers(bp).block,ctbp)
      bptimes=rawvis2_all(inbp).uthrs
      oplot,[min(bptimes),min(bptimes)],!y.crange,line=2
      xyouts,min(bptimes),!y.crange(1),'!C'+targets(bp),orient=-90
   endfor


   if (zoomup eq 1) then legend,['ZOOM UP'],box=0


   ;oplot,uthrs,result_fun,thick=2
   ;oplot,uthrs,result_fun+booterror_fun,thick=.5
   ;oplot,uthrs,result_fun-booterror_fun,thick=.5


   print,' '
   print,' Click left of axis to continue (or to goback to unzoomed view)'
   print,' Click right of axis to change smoothing length'
   print,' Click near datapoint to remove a time sample.'
   print,' Click below axis twice to specify a range to remove.'
   print,' Click above plot twice to specify a range to ZOOM up'
   ;stop
   if noclick eq 1 or fasttrack eq 1 then x0=!x.crange(0)-1 else $
      cursor,x0,y0
   wait,.2
   if (x0 ge !x.crange(1)) then begin
      smoothlength = float(dialog_input(prompt="Smoothing (1/3) Length (hrs):",initial=smoothlength))
      goto,skipsetup1_t3phi
   endif
   if (x0 le !x.crange(1) and x0 ge !x.crange(0)) then begin
      ;  identify bad thing:
      if (y0 ge !y.crange(0) and y0 lt !y.crange(1)) then begin
         ; ri2at,(x0-ut(intc0))/(xr[1]-xr[0]), (y0-t3phi(intc0))/(yr(1)-yr(0)), tdiff,angle
         temp=convert_coord(x0,y0 ,/data,/to_device)
         temp1=convert_coord(ut[intc0],t3phi[intc0],/data,/to_device)
         ri2at, (temp[0])-reform(temp1[0,*]) , (temp[1])-reform(temp1[1,*]),$
                tdiff,angle

         in0=where(tdiff eq min(tdiff))
         rawt3_all(intc0(in0)).badflag_t3phi(tri,*)=1 ; mark all waves
      endif
      if (y0 lt !y.crange(0)) then begin
         print,'Click 2nd point for range'
         cursor,x1,y1
         wait,.2
         xmin=min([x0,x1])
         xmax=max([x0,x1])

         in0=where(  ut(intc0) ge xmin and ut(intc0) le xmax,ctbad )
         for cb=0,ctbad-1 do begin
            rawt3_all(intc0(in0(cb))).badflag_t3phi(tri,*)=1 ; mark all waves
         endfor

      endif                         ;kill-range
      if (y0 ge !y.crange(1)) then begin ; new zoomwindow
         print,'Click 2nd point for zoom range'
         cursor,x1,y1
         wait,.2
         xmin=min([x0,x1])
         xmax=max([x0,x1])
         zoomxr=[xmin,xmax]

         in0=where(  ut(intc0) ge xmin and ut(intc0) le xmax,ctbad )
         zoomyr=[  min([0,t3phi(intc0(in0))]),max(t3phi(intc0(in0)))]
         zoomup=1

      endif                     ; zoom kill

      goto,skipsetup1_t3phi
   endif
   if (zoomup eq 1) then begin
      zoomup=0
      print,' Back to normal view of all data'
      goto,skipsetup1_t3phi     ; jdm should rewrite this using routines... or gui.
   endif

   ; where do I save this calibrated info?
   ;t3_chop.. etc.

   ; norm factors don't matter for this calibration. EXCEPT they have a different set of badflags..
   ; Now, because the closure phase is *required* for averaging of the
   ; triple amp, I will flag all the t3amp to be consistent with t3phi.
   print,' Calculating Calibration Table for All Spectral Channels...'
   skipthis2:
   for s=0,dimy-1 do begin

      inc0=where(rawt3_all.calib eq 1 and rawt3_all.badflag_t3phi(tri,s) eq 0)
      int0=where(rawt3_all.target eq 1 and rawt3_all.badflag_t3phi(tri,s) eq 0)
      intc0=where( (rawt3_all.calib eq 1 or rawt3_all.target eq 1) and $
                   rawt3_all.badflag_t3phi(tri,s) eq 0)


      t3_chop.badflag_t3amp(tri,s)=rawt3_all.badflag_t3phi(tri,s)
      t3_chop.badflag_t3phi(tri,s)=rawt3_all.badflag_t3phi(tri,s)
      t3_fiber.badflag_t3amp(tri,s)=rawt3_all.badflag_t3phi(tri,s)
      t3_fiber.badflag_t3phi(tri,s)=rawt3_all.badflag_t3phi(tri,s)
      t3_shut.badflag_t3amp(tri,s)=rawt3_all.badflag_t3phi(tri,s)
      t3_shut.badflag_t3phi(tri,s)=rawt3_all.badflag_t3phi(tri,s)
      t3_daqsync.badflag_t3amp(tri,s)=rawt3_all.badflag_t3phi(tri,s)
      t3_daqsync.badflag_t3phi(tri,s)=rawt3_all.badflag_t3phi(tri,s)
      t3_xchan.badflag_t3amp(tri,s)=rawt3_all.badflag_t3phi(tri,s)
      t3_xchan.badflag_t3phi(tri,s)=rawt3_all.badflag_t3phi(tri,s)

   endfor   ; spec channels

endfor  ; triangles


;----------------------------------------------------
; Now do t3amp
; First apply flags based on telescope fluxes

skipv2:
for psn=0,num_norms-1 do begin  ;

   rawt3_all=rawt3_all_original ; will get reset each time we change normalization since each method may have different
   ;bad files..... unfortunately.

   ; The closure phase set might have marked some files as bad and so
   ; we need to reflect that here. Unfortunatley, this requires a case statement.
   if norm_labels[psn] eq 'CHOP' then begin
      rawt3_all.badflag_t3amp = t3_chop.badflag_t3amp
      rawt3_all.badflag_t3phi = t3_chop.badflag_t3phi
      calt3_temp = calt3_chop
   endif
   if norm_labels[psn] eq 'FIBER' then begin
      rawt3_all.badflag_t3amp = t3_fiber.badflag_t3amp
      rawt3_all.badflag_t3phi = t3_fiber.badflag_t3phi
      calt3_temp = calt3_fiber
   endif
   if norm_labels[psn] eq 'SHUTTER' then begin
      rawt3_all.badflag_t3amp = t3_shut.badflag_t3amp
      rawt3_all.badflag_t3phi = t3_shut.badflag_t3phi
      calt3_temp = calt3_shut
   endif
   if norm_labels[psn] eq 'DAQSYNC' then begin
      rawt3_all.badflag_t3amp = t3_daqsync.badflag_t3amp
      rawt3_all.badflag_t3phi = t3_daqsync.badflag_t3phi
      calt3_temp = calt3_daqsync
   endif
   if norm_labels[psn] eq 'XCHAN' then begin
      rawt3_all.badflag_t3amp = t3_xchan.badflag_t3amp
      rawt3_all.badflag_t3phi = t3_xchan.badflag_t3phi
      calt3_temp = calt3_xchan
   endif

   ; Now ADD the BAD FLAGS for FLUX_TEL_FLAG
   ;MARK
   for t=0,ntel-1 do begin
      badin = where(flux_tel_flag(*,t*num_norms+psn) eq 1,ct)
      if ct gt 0 then begin
         tin=[where(itrihelp eq t),where(jtrihelp eq t), where(ktrihelp eq t)]
         tin=tin(where(tin ne -1,cttri))
         for t2=0,cttri-1 do begin
            rawt3_all(badin).badflag_t3amp(tin(t2),*)=1
            rawt3_all(badin).badflag_t3phi(tin(t2),*)=1
         endfor
      endif
   endfor

   ; Now add in bad flags based on calt3_temp [nans and 00s]
   for tri =0,ntriangles -1 do begin
      for s=0,dimy-1 do begin
         badin=where(calt3_temp.t3amp[tri,s] eq -1,ctbad)
         if ctbad gt 0 then rawt3_all[badin].badflag_t3amp[tri,s]=1
      endfor
   endfor

   ; Now save badflag informations to the appropriate structure
   for tri=0,ntriangles-1 do begin
      for s=0,dimy-1 do begin

         if norm_labels[psn] eq 'CHOP' then begin
            t3_chop.badflag_t3amp(tri,s)=rawt3_all.badflag_t3amp(tri,s)
         endif
         if norm_labels[psn] eq 'FIBER' then begin
            t3_fiber.badflag_t3amp(tri,s)=rawt3_all.badflag_t3amp(tri,s)
         endif
         if norm_labels[psn] eq 'SHUT' then begin
            t3_shut.badflag_t3amp(tri,s)=rawt3_all.badflag_t3amp(tri,s)
         endif
         if norm_labels[psn] eq 'DAQSYNC' then begin
            t3_daqsync.badflag_t3amp(tri,s)=rawt3_all.badflag_t3amp(tri,s)
         endif
         if norm_labels[psn] eq 'XCHAN' then begin
            t3_xchan.badflag_t3amp(tri,s)=rawt3_all.badflag_t3amp(tri,s)
         endif

      endfor                    ; loop over spectral channels.
   endfor                       ; triangles
endfor                          ; loop over psnorm methods.

;----------------------------------------------------
; Sigma clipping of closure amplitudes
repeatsigmaclipt3a:
if (fasttrack eq 0) then begin
   ans_clip = ''
   print,'Would you like to apply sigma clipping to closure amplitudes?'
   print,'<Y> or N'
   read,ans_clip
endif else ans_clip = 'N'

if (ans_clip eq 'N') then goto,skipclip_amp_all

sigrej = 3.0
ansrej = 'N'
print,'Sigma clipping rejection threshold: '
print,sigrej
print,'Would you like to change rejection threshold?'
print,'Y of <N>'
read,ansrej
if (ansrej eq 'Y') then begin
   print,'Enter new rejection threshold:'
   read,sigrej
endif
ans_cont = ''

print,'------------------------------------'
print,'Sigma clipping of closure amplitudes'
print,'Rejection threshold (N*sigma): ',sigrej
if (skip_edge) then begin
   print,'NOTE: Edge channels will not be used to apply sigma clipping flags'
   print,'Hit enter to continue'
   read, ans_cont
endif

for psn=0,num_norms-1 do begin  ;

   ; The closure phase set might have marked some files as bad and so
   ; we need to reflect that here. Unfortunatley, this requires a case statement.
   if norm_labels[psn] eq 'CHOP' then begin
      rawt3_all = t3_chop
      calt3_temp = calt3_chop
   endif
   if norm_labels[psn] eq 'FIBER' then begin
      rawt3_all = t3_fiber
      calt3_temp = calt3_fiber
   endif
   if norm_labels[psn] eq 'SHUTTER' then begin
      rawt3_all = t3_shut
      calt3_temp = calt3_shut
   endif
   if norm_labels[psn] eq 'DAQSYNC' then begin
      rawt3_all = t3_daqsync
      calt3_temp = calt3_daqsync
   endif
   if norm_labels[psn] eq 'XCHAN' then begin
      rawt3_all = t3_xchan
      calt3_temp = calt3_xchan
   endif

   ut=rawt3_all.uthrs
   flag_sigclip = rawt3_all.badflag_t3amp

   for tri=0,ntriangles-1 do begin

      inc0=where(rawt3_all.calib eq 1 and rawt3_all.badflag_t3amp(tri,spec) eq 0)
      int0=where(rawt3_all.target eq 1 and rawt3_all.badflag_t3amp(tri,spec) eq 0)
      intc0=where( (rawt3_all.calib eq 1 or rawt3_all.target eq 1) and $
                   rawt3_all.badflag_t3amp(tri,spec) eq 0,ct_tc)

      ; Check it there are data on this triangle
      if (ct_tc eq 0) then goto,skipclip_amp0

      ; Loop through observing blocks and spectral channels, and do sigma clipping
      for bp=0,nblocks-1 do begin

         ; Check if data exist for this target block
         inbp=where(rawt3_all(intc0).block eq headers(bp).block,ctbp)
         if (ctbp eq 0) then goto,skipclip_amp1 ; No data on this target block

         maxbad = 0.0

         for ns=0, dimy-1 do begin

            nbad_sum = 0.0
            nbad = 1            ; set to go into while loop
            while (nbad gt 0) do begin

               ; Extract data for this triangle and spectral channel
               intc0=where( (rawt3_all.calib eq 1 or rawt3_all.target eq 1) and $
                            flag_sigclip(tri,ns,*) eq 0,ct_tc)

               t3amp = rawt3_all.t3amp(tri,ns) / calt3_temp.t3amp(tri,ns)
               t3amp_err = rawt3_all.t3amp_err(tri,ns) / calt3_temp.t3amp(tri,ns)

               t3amp_temp = t3amp(intc0)
               t3amp_err_temp = t3amp_err(intc0)
               uttemp = ut(intc0)
               blocktemp = rawt3_all(intc0).block

               ; Extract data block for this spectral channel
               inbp=where(blocktemp eq headers(bp).block,ctbp)

               medt3 = median(t3amp_temp(inbp))
               ;stdevt3 = stdev(t3amp_temp(inbp))
               stdevt3=robust_sigma(t3amp_temp(inbp))
               igood = where((t3amp_temp(inbp) le medt3 + sigrej*stdevt3) and (t3amp_temp(inbp) ge medt3 - sigrej*stdevt3), ngood)
               ibad = where((t3amp_temp(inbp) gt medt3 + sigrej*stdevt3) or (t3amp_temp(inbp) lt medt3 - sigrej*stdevt3), nbad)

               if (nbad gt 0) then flag_sigclip(tri,ns,intc0(inbp(ibad))) = 1

               nbad_sum = nbad_sum + nbad
               if (nbad_sum ge maxbad) then maxbad = nbad_sum

            endwhile   ; while data points are still rejected


            ; Plot good data and sigma clipped data
            ; Need to read flags prior to sigma clipping (use rawt3_all.badflag)
            intc0=where( (rawt3_all.calib eq 1 or rawt3_all.target eq 1) and $
                         rawt3_all.badflag_t3amp(tri,spec) eq 0,ct_tc)

            ; Extract data block for this spectral channel
            inbp=where(rawt3_all(intc0).block eq headers(bp).block,ctbp)

            if (ns eq 0) then begin
               plotsym,0

               t3amp_spec=rawt3_all(intc0(inbp)).t3amp(tri,*) / calt3_temp(intc0(inbp)).t3amp(tri,*)

               mint3 = min(t3amp_spec)
               maxt3 = max(t3amp_spec)
               midt3 = (maxt3 + mint3)/2.0
               t3range = maxt3 - mint3

               mintime = min(ut(intc0(inbp)))
               maxtime = max(ut(intc0(inbp)))
               utrange = maxtime - mintime

               ; Figure out number of rows and columns to plot all spectral channels
               ; in one multiplot:
               nrow = fix(sqrt(dimy))        ; 2 rows for 8 spectral channels
               ncol = fix(round(dimy/nrow))  ; 4 columns for 8 spectral channels

               erase
               multiplot,/reset
               cleanplot
               !P.multi = [0,ncol,nrow]
               multiplot,[ncol,nrow],/init,/doxaxis,/doyaxis,ygap=0.03,xgap=0.02
            endif

            igood = where(flag_sigclip(tri,ns,intc0(inbp)) eq 0,ngood)
            ibad = where(flag_sigclip(tri,ns,intc0(inbp)) eq 1,nbad)

            multiplot,/doxaxis,/doyaxis,ygap=0.03,xgap=0.02

            ploterror,ut(intc0(inbp(igood))),t3amp(intc0(inbp(igood))),t3amp_err(intc0(inbp(igood))), psym=4, $
                      xrange=[mintime-0.1*trange,maxtime+0.1*utrange], $
                      yrange=[mint3-0.1*t3range,maxt3+0.1*t3range], $
                      xtit="Time (hrs)",ytit="Triple Amplitude",xstyle=1,ystyle=1,$
                      title="Triangle Num: "+strcompress(string(fix(tri)),/remove_all)+$
                      '('+get_triname(tri,mircarray=mircarray)+')'+titlabel+$
                      ' Chan '+ strcompress(string(fix(ns)),/remove_all)


            if (nbad gt 0) then oploterror,ut(intc0(inbp(ibad))),t3amp(intc0(inbp(ibad))),t3amp_err(intc0(inbp(ibad))), psym=7,errcolor=200,color=200

            xyouts,mintime-0.05*utrange,maxt3+0.05*t3range,headers(bp).source+' ('+strcompress(string(headers(bp).block),/remove_all)+')'

         endfor                 ; loop over spectral channels

         ; Count number of rejected measurements (across all channels)
         nflag_sum = 0
         if (skip_edge) then begin ; skip edge channels when counting flags
            for nf=0, ctbp-1 do begin
               bflag = flag_sigclip(tri,1:dimy-2,intc0(inbp(nf)))
               indflag = where(bflag eq 1, nflag)
               if (nflag gt 0) then nflag_sum = nflag_sum + 1
            endfor
         endif else begin       ; include edge channels when counting flags
            for nf=0, ctbp-1 do begin
               bflag = flag_sigclip(tri,*,intc0(inbp(nf)))
               indflag = where(bflag eq 1, nflag)
               if (nflag gt 0) then nflag_sum = nflag_sum + 1
            endfor
         endelse

         print,'Block: ',headers(bp).block,'  Source: ',headers(bp).source
         print,'Diamonds - good data; X - sigma clipped'
         print,'Number of measurements on target (per spectral channel): '
         print, n_elements(inbp)
         print,'Maximum number of sigma clipped measurements per spectral channel: '
         print, maxbad
         print,'Number of uniquely rejected measurements (not including edge channels): '
         print, nflag_sum
         print,'Do you want to reject these measurements?'
         print,'<Y> or N:'
         read,ans_cont

         if (ans_cont eq 'N') then flag_sigclip(tri,*,intc0(inbp)) = 0

         skipclip_amp1:

      endfor                    ; loop over observing blocks

      ; Loop through all files on this triangle - set bad flags for all
      ; spectral channels

      if (skip_edge) then begin  ; skip edge channels when applying sigma clipping flags
         for nf=0,nfiles-1 do begin
            bflag = flag_sigclip(tri,1:dimy-2,nf)
            indflag = where(bflag eq 1, nflag)
            if (nflag gt 0) then rawt3_all(nf).badflag_t3amp(tri,*) = 1
         endfor
      endif else begin  ; include all channels when applying sigma clipping flags
         for nf=0,nfiles-1 do begin
            bflag = flag_sigclip(tri,*,nf)
            indflag = where(bflag eq 1, nflag)
            if (nflag gt 0) then rawt3_all(nf).badflag_t3amp(tri,*) = 1
         endfor
      endelse

      skipclip_amp0:

      ; where do I save this calibrated info?
      ;t3_chop.. etc.
      ;----Save calibration factors and the badflag arrays but do not apply yet. This will be done
      ; in last step after correlation analysis.
      ; This is a good time to deal with the correlated errors from the calibrator....but how? -jdm
      print,' Calculating Calibration Table for All Spectral Channels...'

      for s=0,dimy-1 do begin

         inc0=where(rawt3_all.calib eq 1 and rawt3_all.badflag_t3amp(tri,s) eq 0)
         int0=where(rawt3_all.target eq 1 and rawt3_all.badflag_t3amp(tri,s) eq 0)
         intc0=where( (rawt3_all.calib eq 1 or rawt3_all.target eq 1) and $
                      rawt3_all.badflag_t3amp(tri,s) eq 0)

         if norm_labels[psn] eq 'CHOP' then begin
            t3_chop.badflag_t3amp(tri,s)=rawt3_all.badflag_t3amp(tri,s)
         endif
         if norm_labels[psn] eq 'FIBER' then begin
            t3_fiber.badflag_t3amp(tri,s)=rawt3_all.badflag_t3amp(tri,s)
         endif
         if norm_labels[psn] eq 'SHUT' then begin
            t3_shut.badflag_t3amp(tri,s)=rawt3_all.badflag_t3amp(tri,s)
         endif
         if norm_labels[psn] eq 'DAQSYNC' then begin
            t3_daqsync.badflag_t3amp(tri,s)=rawt3_all.badflag_t3amp(tri,s)
         endif
         if norm_labels[psn] eq 'XCHAN' then begin
            t3_xchan.badflag_t3amp(tri,s)=rawt3_all.badflag_t3amp(tri,s)
         endif


      endfor                    ; loop over spectral channels.
   endfor                       ; triangles
endfor                          ; loop over psnorm methods.


goto,repeatsigmaclipt3a
skipclip_amp_all:
print,'----------------------------------------------'
print,'Finished sigma clipping the closure amplitudes'
print,'Begin visual inspection and flagging of closure amplitudes.'
print,'Hit enter to continue'
read,ans_cont

multiplot,/reset
cleanplot
!P.multi=0
;----------------------------------------------------
; Visual inspection of t3amp
; now do t3amp
inc=where(rawt3_all.calib eq 1)
int=where(rawt3_all.target eq 1)
intc=where(rawt3_all.calib eq 1 or rawt3_all.target eq 1)


smoothlength=3

for psn=0,num_norms-1 do begin  ;

   rawt3_all=rawt3_all_original ; will get reset each time we change normalization since each method may have different
   ;bad files..... unfortunately.

   ; The closure phase set might have marked some files as bad and so
   ; we need to reflect that here. Unfortunatley, this requires a case statement.
   if norm_labels[psn] eq 'CHOP' then begin
      rawt3_all.badflag_t3amp = t3_chop.badflag_t3amp
      rawt3_all.badflag_t3phi = t3_chop.badflag_t3phi
      calt3_temp = calt3_chop
   endif
   if norm_labels[psn] eq 'FIBER' then begin
      rawt3_all.badflag_t3amp = t3_fiber.badflag_t3amp
      rawt3_all.badflag_t3phi = t3_fiber.badflag_t3phi
      calt3_temp = calt3_fiber
   endif
   if norm_labels[psn] eq 'SHUTTER' then begin
      rawt3_all.badflag_t3amp = t3_shut.badflag_t3amp
      rawt3_all.badflag_t3phi = t3_shut.badflag_t3phi
      calt3_temp = calt3_shut
   endif
   if norm_labels[psn] eq 'DAQSYNC' then begin
      rawt3_all.badflag_t3amp = t3_daqsync.badflag_t3amp
      rawt3_all.badflag_t3phi = t3_daqsync.badflag_t3phi
      calt3_temp = calt3_daqsync
   endif
   if norm_labels[psn] eq 'XCHAN' then begin
      rawt3_all.badflag_t3amp = t3_xchan.badflag_t3amp
      rawt3_all.badflag_t3phi = t3_xchan.badflag_t3phi
      calt3_temp = calt3_xchan
   endif

   skip1:

   for tri=0,ntriangles-1 do begin

      zoomup=0
      skipsetup_t3amp:
      ;smoothlength = float(dialog_input(prompt="Smoothing (1/3) Length (hrs):",initial=smoothlength))
      skipsetup1_t3amp:

      inc0=where(rawt3_all.calib eq 1 and rawt3_all.badflag_t3amp(tri,spec) eq 0)
      int0=where(rawt3_all.target eq 1 and rawt3_all.badflag_t3amp(tri,spec) eq 0)
      intc0=where( (rawt3_all.calib eq 1 or rawt3_all.target eq 1) and $
                   rawt3_all.badflag_t3amp(tri,spec) eq 0,ct_tc)
      if ct_tc eq 0 then goto,skipthis3


      mintimes=min(rawt3_all(intc0).uthrs)
      maxtimes=max(rawt3_all(intc0).uthrs)
      num_interp = fix((maxtimes-mintimes+.4)/(smoothlength/5.) ) > 20
      uthrs=pvector([mintimes-.2,maxtimes+.2,num_interp])
      numsrc=n_elements(int0)

      !p.multi=0
      titlabel='('+norm_labels(psn)+' NORM)'

      if norm_labels[psn] eq 'CHOP' then begin

         ;  result_fun = caltables_t3amp(uthrs,rawt3_all(intc0).uthrs,$
         ;     rawt3_all(intc0).t3amp(tri,spec),rawt3_all(intc0).t3amp_err(tri,spec),$
         ;     rawt3_all(intc0).t3phi(tri,spec),$
         ;     calt3_chop(intc0).t3amp(tri,spec),smooth_fun,booterror=booterror_fun,nsamps=5)


         t3amp=rawt3_all.t3amp(tri,spec) / calt3_chop.t3amp(tri,spec)
         t3amp_err=rawt3_all.t3amp_err(tri,spec) / calt3_chop.t3amp(tri,spec)
      endif
      if norm_labels[psn] eq 'FIBER' then begin
         ;result_fun = caltables_t3amp(uthrs,rawt3_all(intc0).uthrs,$
         ;     rawt3_all(intc0).t3amp(tri,spec),rawt3_all(intc0).t3amp_err(tri,spec),$
         ;     rawt3_all(intc0).t3phi(tri,spec),$
         ;     calt3_fiber(intc0).t3amp(tri,spec),smooth_fun,booterror=booterror_fun,nsamps=5)

         t3amp=rawt3_all.t3amp(tri,spec) / calt3_fiber.t3amp(tri,spec)
         t3amp_err=rawt3_all.t3amp_err(tri,spec) / calt3_fiber.t3amp(tri,spec)
      endif
      if norm_labels[psn] eq 'SHUTTER' then begin
         ;result_fun = caltables_t3amp(uthrs,rawt3_all(intc0).uthrs,$
         ;     rawt3_all(intc0).t3amp(tri,spec),rawt3_all(intc0).t3amp_err(tri,spec),$
         ;     rawt3_all(intc0).t3phi(tri,spec),$
         ;     calt3_shut(intc0).t3amp(tri,spec),smooth_fun,booterror=booterror_fun,nsamps=5)

         t3amp=rawt3_all.t3amp(tri,spec) / calt3_shut.t3amp(tri,spec)
         t3amp_err=rawt3_all.t3amp_err(tri,spec) / calt3_shut.t3amp(tri,spec)
      endif
      if norm_labels[psn] eq 'DAQSYNC' then begin
         ;result_fun = caltables_t3amp(uthrs,rawt3_all(intc0).uthrs,$
         ;     rawt3_all(intc0).t3amp(tri,spec),rawt3_all(intc0).t3amp_err(tri,spec),$
         ;     rawt3_all(intc0).t3phi(tri,spec),$
         ;     calt3_daqsync(intc0).t3amp(tri,spec),smooth_fun,booterror=booterror_fun,nsamps=5)

         t3amp=rawt3_all.t3amp(tri,spec) / calt3_daqsync.t3amp(tri,spec)
         t3amp_err=rawt3_all.t3amp_err(tri,spec) / calt3_daqsync.t3amp(tri,spec)
      endif
      if norm_labels[psn] eq 'XCHAN' then begin
         ;result_fun = caltables_t3amp(uthrs,rawt3_all(intc0).uthrs,$
         ;     rawt3_all(intc0).t3amp(tri,spec),rawt3_all(intc0).t3amp_err(tri,spec),$
         ;     rawt3_all(intc0).t3phi(tri,spec),$
         ;     calt3_xchan(intc0).t3amp(tri,spec),smooth_fun,booterror=booterror_fun,nsamps=5)

         t3amp=rawt3_all.t3amp(tri,spec) / calt3_xchan.t3amp(tri,spec)
         t3amp_err=rawt3_all.t3amp_err(tri,spec) / calt3_xchan.t3amp(tri,spec)
      endif


      plotsym,0,/fill

      ut=rawt3_all.uthrs
      if (zoomup eq 0) then begin ; standard scaling.
         yr=[  min([0,t3amp(intc0)]),max(t3amp(intc0))]
         xr=[mintimes-.2,maxtimes+.2]
      endif else begin
         yr=zoomyr
         xr=zoomxr              ;defined elsewheres
      endelse

      ploterror,ut(intc0),t3amp(intc0),t3amp_err(intc0), psym=4,$
                xtit="Time (hrs)",ytit="Triple Amplitude",xr=xr,xst=1,yr=yr,$
                title="Triangle Num:"+string(fix(tri)) + '('+get_triname(tri,mircarray=mircarray)+')'+titlabel
      plotsym,0,/fill
      oplot,[-1e6,1e6],[0,0],line=1

     ; Draw all boundaries between blocks and label.
      for bp=0,nblocks-1 do begin
         inbp=where(rawvis2_all.block eq headers(bp).block,ctbp)
         bptimes=rawvis2_all(inbp).uthrs
         oplot,[min(bptimes),min(bptimes)],!y.crange,line=2
         xyouts,min(bptimes),!y.crange(1),'!C'+targets(bp),orient=-90
      endfor


      if (zoomup eq 1) then legend,['ZOOM UP'],box=0

      ;oplot,uthrs,result_fun,thick=2
      ;oplot,uthrs,result_fun+booterror_fun,thick=.5
      ;oplot,uthrs,result_fun-booterror_fun,thick=.5


      print,' '
      print,' Click left of axis to continue (or to goback to unzoomed view)'
      print,' Click right of axis to change smoothing length'
      print,' Click near datapoint to remove a time sample.'
      print,' Click below axis twice to specify a range to remove.'
      print,' Click above plot twice to specify a range to ZOOM up'
      ;stop
      if noclick eq 1 or fasttrack eq 1 then x0=!x.crange(0)-1 else $
         cursor,x0,y0
      wait,.2
      if (x0 ge !x.crange(1)) then begin
         smoothlength = float(dialog_input(prompt="Smoothing (1/3) Length (hrs):",initial=smoothlength))
         goto,skipsetup1_t3amp
      endif

      if (x0 le !x.crange(1) and x0 ge !x.crange(0)) then begin
         ;  identify bad thing:
         if (y0 ge !y.crange(0) and y0 lt !y.crange(1)) then begin

            temp=convert_coord(x0,y0 ,/data,/to_device)
            temp1=convert_coord(ut[intc0],t3amp[intc0],/data,/to_device)
            ri2at, (temp[0])-reform(temp1[0,*]) , (temp[1])-reform(temp1[1,*]),$
                   tdiff,angle

            ; ri2at,(x0-ut(intc0))/(xr[1]-xr[0]), (y0-t3amp(intc0))/(yr(1)-yr(0)), tdiff,angle
            in0=where(tdiff eq min(tdiff))
            rawt3_all(intc0(in0)).badflag_t3amp(tri,*)=1 ; mark all waves
         endif
         if (y0 lt !y.crange(0)) then begin
            print,'Click 2nd point for range'
            cursor,x1,y1
            wait,.2
            xmin=min([x0,x1])
            xmax=max([x0,x1])

            in0=where(  ut(intc0) ge xmin and ut(intc0) le xmax,ctbad )
            for cb=0,ctbad-1 do begin
               rawt3_all(intc0(in0(cb))).badflag_t3amp(tri,*)=1 ; mark all waves
            endfor

         endif                      ;kill-range
         if (y0 ge !y.crange(1)) then begin ; new zoomwindow
            print,'Click 2nd point for zoom range'
            cursor,x1,y1
            wait,.2
            xmin=min([x0,x1])
            xmax=max([x0,x1])
            zoomxr=[xmin,xmax]

            in0=where(  ut(intc0) ge xmin and ut(intc0) le xmax,ctbad )
            zoomyr=[  min([0,t3amp(intc0(in0))]),max(t3amp(intc0(in0)))]
            zoomup=1

         endif                  ; zoom kill

         goto,skipsetup1_t3amp
      endif
      if (zoomup eq 1) then begin
         zoomup=0
         print,' Back to normal view of all data'
         goto,skipsetup1_t3amp  ; jdm should rewrite this using routines... or gui.
      endif

      ; where do I save this calibrated info?
      ;t3_chop.. etc.
      ;----Save calibration factors and the badflag arrays but do not apply yet. This will be done
      ; in last step after correlation analysis.
      ; This is a good time to deal with the correlated errors from the calibrator....but how? -jdm
      skipthis3:
      print,' Calculating Calibration Table for All Spectral Channels...'

      for s=0,dimy-1 do begin

         inc0=where(rawt3_all.calib eq 1 and rawt3_all.badflag_t3amp(tri,s) eq 0)
         int0=where(rawt3_all.target eq 1 and rawt3_all.badflag_t3amp(tri,s) eq 0)
         intc0=where( (rawt3_all.calib eq 1 or rawt3_all.target eq 1) and $
                      rawt3_all.badflag_t3amp(tri,s) eq 0)

         if norm_labels[psn] eq 'CHOP' then begin
            t3_chop.badflag_t3amp(tri,s)=rawt3_all.badflag_t3amp(tri,s)
         endif
         if norm_labels[psn] eq 'FIBER' then begin
            t3_fiber.badflag_t3amp(tri,s)=rawt3_all.badflag_t3amp(tri,s)
         endif
         if norm_labels[psn] eq 'SHUT' then begin
            t3_shut.badflag_t3amp(tri,s)=rawt3_all.badflag_t3amp(tri,s)
         endif
         if norm_labels[psn] eq 'DAQSYNC' then begin
            t3_daqsync.badflag_t3amp(tri,s)=rawt3_all.badflag_t3amp(tri,s)
         endif
         if norm_labels[psn] eq 'XCHAN' then begin
            t3_xchan.badflag_t3amp(tri,s)=rawt3_all.badflag_t3amp(tri,s)
         endif


      endfor                    ; loop over spectral channels.

   endfor                       ; triangles
   skipclicks2:
endfor                          ; loop over psnorm methods.


;stop
if restoreflags eq 0 then begin
vis2_xchan_badflag = vis2_xchan.badflag
t3_xchan_badflag_t3phi= t3_xchan.badflag_t3phi
t3_xchan_badflag_t3amp= t3_xchan.badflag_t3amp
save,vis2_xchan_badflag,t3_xchan_badflag_t3phi,t3_xchan_badflag_t3amp,file='lastflags_all.idlvar'
endif else begin; RESTORE !!
if n_elements(vis2_xchan.badflag) ne n_elements(vis2_xchan_badflag) then begin
  print,'You appear to be using flags that are not compatible...'
  stop
endif
vis2_xchan.badflag=vis2_xchan_badflag
t3_xchan.badflag_t3phi=t3_xchan_badflag_t3phi
t3_xchan.badflag_t3amp=t3_xchan_badflag_t3amp
endelse


;stop
jdmskip:
; We are essentially done I think.....
;---------------------------------------------------------
; probably should do some kind of closure amplitude average
; and save the averaged info in a custom table or idlvar file/table.
;---------------------------------------------------------

;---------------------------------------------------------------
; To do:  1. decide how to average the blocks...
;         2. save data
; 	  3. separately track errors due to calibrator size uncertainties?
;         4. closure amp stuff would be good nice too.
;--------------------------------------------------------------------

;-------------------------------
; OIFITS PREP
;-------------------------------
; Ask for unique identifier for saving data.

identifier = 'XYZ'
datestring=get_fluordate()
identifier =dialog_input(title='Choose SaveName',$
 initial=identifier, prompt = "Enter your Initials (e.g, JDM)" );  Identifier for Saving:")
identifier0=identifier+'_'+datestring

skipahead_mark:
;------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------
; RAWOUTPUT MODULE
;---------------------------------------------------------
use_scatter=0
skipahead_mark2:

if use_scatter eq 0 then identifier = identifier0; +'_'+mode2

; probably don't need to set these up again, but just in case.
; setup other common things.
specname = instructure(0).mode ; all are supposed to be same mode.
lambda=bestlambda
lambdaerr=abs(deriv(lambda))

; The spectrometer info will be setup by the
; setup_combiner script in the future....

telindex=intarr(ntel)
for tt=0,ntel-1 do begin
 if tt eq 0 then tel=mircarray.b1
 if tt eq 1 then tel=mircarray.b2
 if tt eq 2 then tel=mircarray.b3
 if tt eq 3 then tel=mircarray.b4
 if tt eq 4 then tel=mircarray.b5
 if tt eq 5 then tel=mircarray.b6

; TEL ORDER: s1,s2,e1,w1,w2,e2
  if tel eq 'S1' then teli=0
  if tel eq 'S2' then teli=1
  if tel eq 'E1' then teli=2
  if tel eq 'W1' then teli=3
  if tel eq 'W2' then teli=4
  if tel eq 'E2' then teli=5
 telindex(tt)=teli
endfor





;----------------------------------------


; When saving date, will loop over the norm methods and loop over the
; targets.
;
; also have to restore the original values for calibrators since the
; structures have been 'calibrated' by diameters.
;

for psn=0,num_norms -1 do begin

   normlabel=norm_labels(psn)
   if normlabel eq 'CHOP' then begin ;chop
      vis2_norm = vis2_chop          ;contains info on calibration error
      calvis2_norm = calvis2_chop    ; contains info the measurement error
      t3_norm = t3_chop
      calt3_norm = calt3_chop
   endif
   if normlabel eq 'FIBER' then begin ;chop
      vis2_norm = vis2_fiber
      calvis2_norm = calvis2_fiber
      t3_norm = t3_fiber
      calt3_norm = calt3_fiber
   endif
   if normlabel eq 'SHUTTER' then begin ;chop
      vis2_norm = vis2_shut
      calvis2_norm = calvis2_shut
      t3_norm = t3_shut
      calt3_norm =calt3_shut
   endif
   if normlabel eq 'DAQSYNC' then begin ;chop
      vis2_norm = vis2_daqsync
      calvis2_norm = calvis2_daqsync
      t3_norm = t3_daqsync
      calt3_norm =calt3_daqsync
   endif
   if normlabel eq 'XCHAN' then begin ;chop
      vis2_norm = vis2_xchan
      calvis2_norm = calvis2_xchan
      t3_norm = t3_xchan
      calt3_norm =calt3_xchan

      ; print do xtalk check. this could in principle be done for all norm types,since
      ; ps_norms doesn't enter in that much. that means could be used for old (non-xchan]
      ; data by just dropping in this section of the code !! but right now only 6-T data but
      ; this part could be easily changed in mirc_xtalkfit
      ;stop

      hey=dialog_checklist(['Yes. apply xtalk correction (for bright targets with low vis)','No'], tit='Do you want to apply crosstalk correction?',$
      /exclusive, initial=0)

      xtalk_flag=reform(where(hey eq 0)) ;1 =yes, 0=no
      if xtalk_flag eq 1 then begin ; do xtalk
      cutoff_vis=.1
      cutoff_vis = float(dialog_input(prompt="Only apply for visibilites less than what cutoff?",initial=cutoff_vis))
      smooth_num=20000 ; only oprate on the "split" # of chunks, which is meant to be the
      ; minimum scale that tings might change.
      ;smooth_num = float(dialog_input(prompt="The true visibility might change over how many blocks?",initial=smooth_num))
      endif ; xtalk .

   endif ;xchan


   ; Loop through all the target AND CALIBRATOR blocks FOR NOCAL !!!
   ; each block will have its own OI-file for later merging.
   tbin=where(headers.target eq 1 or headers.calib eq 1,ntb)
   filelist = strarr(ntb)       ; to help with merging blocks of same files at end
   for i=0,ntb-1 do begin


      ;----- OIFITS -----
      ; temporary just for setting up target info.
      dins=where(rawvis2_all.block eq headers(tbin(i)).block,nd)
      uvins=rawvis2_all(dins).uvin

      obs_date0 = uvinfo( uvins(0) ).date_obs
      ra        = uvinfo( uvins(0) ).ra
      dec       = uvinfo( uvins(0) ).dec
      target    = targets(tbin(i))
      targetname= targets_file(tbin(i)) ; for creating name of oifits file

      ; to proceed we have to next lopop over baselines, because
      ; the number of bad-files is strongly dependent on this.

      ; Config the data arrays for oifits-writing
      ; Do we need different ones for the closure phases since sometimes
      ; the # of bad files are different?

      ;-----------------------
      ;  SETUP FOR VIS2
      ;-----------------------

      dummy = {      uu_meters : dblarr(nbaselines),$
                     vv_meters : dblarr(nbaselines),$
                     obs_time  : dblarr(nbaselines),$ ; might be diff per baseline
                     int_time  : dblarr(nbaselines),$ ; due to data selection
                     mjd       : dblarr(nbaselines),$ ;
                     vis2      : fltarr(nbaselines,dimy),$
                     vis2_err  : fltarr(nbaselines,dimy),$
                     badflag   : intarr(nbaselines,dimy),$
                     telindex1 : intarr(nbaselines), $
                     telindex2 : intarr(nbaselines), $
                     nbaselines: nbaselines,$
                     dimy      : dimy}

      dummy.badflag(*)=1        ; must overwrite with correct values to be valid
      vis2_oifits=replicate(dummy,nsplit_max)

      ;stop
      for base=0,nbaselines-1 do begin
         dins=where(rawvis2_all.block eq headers(tbin(i)).block and $
                    vis2_norm.badflag(base,0) eq 0,nd) ; only count 'good files'
	 ;base on spectral channel 0 (all channels have same badflag here)
         if nd eq 1 then begin
            dins=[dins,dins]
            nd=2
         endif
         if nd ne 0 then begin  ; do anaysis.. otherwise skip.
            uvins=rawvis2_all(dins).uvin
            nuv=n_elements(uvins)
            if nuv ge 2 then median_time = median(uvinfo[uvins].int_time) else median_time=uvinfo[uvins].int_time

            ; nsplit=nint((nd*median_time*1.0/(avg_time*60.))+0.499)
            ; mircblocksplitter, nd, nsplit, minfiles, i0,i1
            temp_times = uvinfo[uvins].time
            mirctimesplitter, temp_times, avg_time*60., minfiles, i0,i1,nsplit=nsplit0

            ; There is a problem if NOTHING survivse thes cut... not sure what happens.
            ; JDM

            if nsplit0 gt 2 then $
               print,'typical number of files per point: ',median(i1-i0),median(i1-i0)*median_time,' secs GOAL: ',avg_Time*60.,' secs'



            ;---------------

            for ns=0,nsplit0-1 do begin
               din=  dins(i0(ns):i1(ns))                   ; pick out the ns-th subset of this block
               uvin=uvins(i0(ns):i1(ns))                   ; pick out the uv coverage for these files
               ; setup information for saving data.
               vis2_oifits(ns).obs_time(base) = mean( uvinfo(uvin).time) ;
               vis2_oifits(ns).int_time(base) = max(uvinfo(uvin).time) - min(uvinfo(uvin).time)
               vis2_oifits(ns).mjd(base)  = mean(uvinfo(uvin).mjd) ;

               vis2_oifits(ns).uu_meters(base) = $
                  mean( uvinfo(uvin).ucoord(ihelp(base),jhelp(base) ) ) ;
               vis2_oifits(ns).vv_meters(base) = $
                  mean( uvinfo(uvin).vcoord(ihelp(base),jhelp(base) ) )
               vis2_oifits(ns).telindex1(base) = telindex(ihelp(base))
               vis2_oifits(ns).telindex2(base) = telindex(jhelp(base))


               ; Do the correlation analysi here using mirc_xtalkfit !!!
               ; JDM 2014Sep

               ;check if xchan..
               ;print,'do xtalk check'
               ;JDM2014
               if xtalk_flag eq 1 then begin ; do xtalk.
                  nd0=n_elements(din)
                  if nd0 gt 5 then begin ; skip this if less than 6 good files!

                     temp_cvis2_base = dblarr(nbaselines,dimy,nd0,2)
                     temp_norms=dblarr(nbaselines,dimy,nd0)
                     temp_cvis2_base[*,*,*,0] = original_rawvis2_all(din).vis2(*,*)
                     temp_cvis2_base[*,*,*,1] = original_rawvis2_all(din).vis2_err(*,*)
                     temp_norms[*,*,*] = calvis2_norm[din].vis2[*,*]
                     ;cvis2_base=fltarr(

                     temp_new_cvis2_base=mirc_xtalkfit(temp_cvis2_base, temp_norms, coefficients, flux_fiber[tbin[i]], comboinfo, cutoff_vis=cutoff_vis, smooth_num=smooth_num,chi2s=chi2s,/single,b_index=[base])
                     ;if base eq 0 then stop ; check if smooth_num works well.
                     if chi2s[base,4] gt 15 then begin
                        print,chi2s[base,4]
                        ploterror,temp_norms[base,4,*],temp_cvis2_base[base,4,*,0],temp_cvis2_base[base,4,*,1],psym=4
                        oploterror,temp_norms[base,4,*],temp_new_cvis2_base[base,4,*,0],temp_cvis2_base[base,4,*,1],psym=5
                        print,'Some issue with the cross talk correction. if you dont know what to do, enter .c to continue'
                        stop
                     endif
                     original_rawvis2_all(din).vis2(base,*)=temp_new_cvis2_base[base,*,*,0] ;technically I shouldn't
                     ;change this array because I am sequentially doing fits..
                     ;but should be ok since corrections only affect weak baselines.

                  endif ; at least 5 good points.
               endif ; do the xtalk correction!




               for j=0,dimy-1 do begin ; wavelength channels
                  v2_numerator     = original_rawvis2_all(din).vis2(base,j)
                  v2_numerator_err = original_rawvis2_all(din).vis2_err(base,j)
                  v2_denominator   =   calvis2_norm(din).vis2(base,j)

                  ; Suppressing weights since individual files show too much variation due to
                  ; the fluctuating coupling.
                  dum = boot_ratio(v2_numerator,v2_denominator,$
                                   result=result)
                  if result[1] eq 0 then begin
                     ;print,'note from john: Check to see if my fix works...'
                     print,' Problem with distributions... crude -fix- '
                     ;stop
                     result[1]=max(original_rawvis2_all[din].vis2_err(base,j)/original_rawvis2_all(din).vis2(base,j))*result[0]
                  endif

                  data_pe = abs(result(1)/result(0))
                  total_pe =  data_pe ; no calibration error here

                  vis2_oifits(ns).vis2(base,j) = result(0)
                  vis2_oifits(ns).vis2_err(base,j) = abs(result(0))*total_pe
                  vis2_oifits(ns).badflag(base,j) = 0

                  abc=where(result ne result,ctbadd)
                  if ctbadd ne 0 then stop


               endfor           ; j dimy

            endfor              ; nsplit
         endif                  ; skip the splits if no files!!

      endfor                    ; baselines

      ; I should average the closure amps here!! and put in infostring
      ; kind of tricky -- with all the badflagging etc.

      ;-----------------------
      ;  SETUP FOR T3
      ;-----------------------

      dummy = {      uu1_meters : dblarr(ntriangles),$
                     vv1_meters : dblarr(ntriangles),$
                     uu2_meters : dblarr(ntriangles),$
                     vv2_meters : dblarr(ntriangles),$
                     obs_time  : dblarr(ntriangles),$ ; might be diff per baseline
                     int_time  : dblarr(ntriangles),$ ; due to data selection
                     mjd       : dblarr(ntriangles),$
                     t3amp      : fltarr(ntriangles,dimy),$
                     t3amp_err  : fltarr(ntriangles,dimy),$
                     t3phi      : fltarr(ntriangles,dimy),$
                     t3phi_err  : fltarr(ntriangles,dimy),$
                     badflag_t3amp   : intarr(ntriangles,dimy),$
                     badflag_t3phi   : intarr(ntriangles,dimy),$
                     telindex1    : intarr(ntriangles), $
                     telindex2    : intarr(ntriangles), $
                     telindex3    : intarr(ntriangles), $
                     ntriangles: ntriangles,$
                     dimy      : dimy}
      dummy.badflag_t3amp(*)=1  ; must overwrite with correct values to be valid
      dummy.badflag_t3phi(*)=1  ; must overwrite with correct values to be valid

      t3_oifits=replicate(dummy,nsplit_max)


      for tri=0,ntriangles-1 do begin
         dins=where(rawt3_all.block eq headers(tbin(i)).block and $
                    t3_norm.badflag_t3phi(tri,0) eq 0,nd) ; only count 'good closure phasefiles'
         if nd ne 0 then begin                            ; skip if no good files in the block
            ;  the triple amps will be averaged based on their own
	    ;  badflags, but the NSPLITS will be determined by the
	    ;  closure phase files and the uv coverages will reflect
	    ;  the uv footprint for good closure phases.
	    ;  This means the mean uv for the triple amps will be slighlty
	    ;  incorrect.

            uvins=rawt3_all(dins).uvin
            nuv=n_elements(uvins)
            if nuv ge 2 then median_time = median(uvinfo[uvins].int_time) else median_time=uvinfo[uvins].int_time

	    ; nsplit=nint((nd*median_time*1.0/(avg_time*60.))+0.499)
	    ; mircblocksplitter, nd, nsplit, minfiles, i0,i1
            temp_times = uvinfo[uvins].time
            mirctimesplitter, temp_times, avg_time*60., minfiles, i0,i1,nsplit=nsplit0

            if nsplit0 gt 2 then $
               print,'typical number of files per point: ',median(i1-i0),median(i1-i0)*median_time,' secs GOAL: ',avg_Time*60.,' secs'


	    ;---------------

            for ns=0,nsplit0-1 do begin
               din=  dins(i0(ns):i1(ns)) ; pick out the ns-th subset of this block
               uvin=uvins(i0(ns):i1(ns)) ; pick out the uv coverage for these files
	       ; setup information for saving data.
               t3_oifits(ns).obs_time(tri) = mean( uvinfo(uvin).time)
               t3_oifits(ns).int_time(tri) = max(uvinfo(uvin).time) - min(uvinfo(uvin).time)
               t3_oifits(ns).mjd(tri)  = mean(uvinfo(uvin).mjd)

               t3_oifits(ns).uu1_meters(tri) = $
                  mean( uvinfo(uvin).ucoord(itrihelp(tri),jtrihelp(tri) ) )
               t3_oifits(ns).vv1_meters(tri) = $
                  mean( uvinfo(uvin).vcoord(itrihelp(tri),jtrihelp(tri) ) )
               t3_oifits(ns).uu2_meters(tri) = $
                  mean( uvinfo(uvin).ucoord(jtrihelp(tri),ktrihelp(tri) ) )
               t3_oifits(ns).vv2_meters(tri) = $
                  mean( uvinfo(uvin).vcoord(jtrihelp(tri),ktrihelp(tri) ) )
               t3_oifits(ns).telindex1(tri) = telindex(itrihelp(tri))
               t3_oifits(ns).telindex2(tri) = telindex(jtrihelp(tri))
               t3_oifits(ns).telindex3(tri) = telindex(ktrihelp(tri))


               for j=0,dimy-1 do begin ; wavelength channels

                  ; NEVER write a t3phi unless there is a good t3amp -- this is a debateable choice since
                  ; we can measure closure phsaes weven when we can't calibrate t3amps, but
                  ; this leads to confusion in the datafiles and also will require a level of data selection
                  ; in imaging codes.

                  goodamp=where(t3_norm(din).badflag_t3amp(tri,j) eq 0 $
                                and t3_norm(din).badflag_t3phi(tri,j) eq 0,ct)
                  if ct ne 0 then begin ; at least one amp + good cp.
                     din1=din(goodamp)

                     ;First Closure Phase.

                     t3phi     = rawt3_all(din1).t3phi(tri,j)
                     t3phi_err = rawt3_all(din1).t3phi_err(tri,j)

                    ; JDM: may need to try different techniques to do this averaging.
                    ;      the approach I'm using is not very outlier resistant so it is relying
                    ;      on good data selection

                     result = boot_angles(t3phi,weights=1./t3phi_err^2)

                     if result[1] eq 0 then begin
                        print,'note from john: Check to see if my fix works...t3phi'
                        print,' Problem with distributions... crude -fix- '

                        ;stop
                        result[1]=max(t3phi_err)
                     endif


                     data_err = result(1)

                     total_err = data_err

                     t3_oifits(ns).t3phi(tri,j) = result(0)
                     t3_oifits(ns).t3phi_err(tri,j) = total_err
                     t3_oifits(ns).badflag_t3phi(tri,j) = 0

                  endif else begin
                     t3_oifits(ns).t3phi(tri,j) =  !values.f_nan
                     t3_oifits(ns).t3phi_err(tri,j) =  !values.f_nan
                     t3_oifits(ns).badflag_t3phi(tri,j) = 1
                  endelse


                  ; NOW do the triple AMP -- we have to  check if the
                  ; triple amp for these files  were flagged as bad.
                  ; Also, we need to USE the Closure phase for the averaging...
                  ; Thus, it will be impossible to have a t3amp results w/o
                  ; a t3phi -- although not the other way around.....

                  goodamp=where(t3_norm(din).badflag_t3amp(tri,j) eq 0 $
                                and t3_norm(din).badflag_t3phi(tri,j) eq 0,ct)
                  if ct ne 0 then begin ; at least one amp + good cp.
                     din1=din(goodamp)
                     t3amp_numerator     = original_rawt3_all(din1).t3amp(tri,j)
                     t3amp_numerator_err = original_rawt3_all(din1).t3amp_err(tri,j)
                     t3amp_phase         = rawt3_all(din1).t3phi(tri,j)
                     t3amp_denominator   =   calt3_norm(din1).t3amp(tri,j)

 	             ;Biases due to weighting are problematic...
                     eq_wts=replicate(1.0,n_elements(t3amp_numerator_err))
                     dum = boot_ratio_t3amp_wt(t3amp_numerator,$
                                               t3amp_denominator,t3amp_phase,eq_wts,$
                                               result=result)
                     if result[1] eq 0 then begin
                        print,'note from john: Check to see if my fix works... t3amp'
                        print,' Problem with distributions... crude -fix- '

                        result[1] = max(t3amp_numerator_err/t3amp_numerator)*result[0]
                     endif



                     data_pe = abs(result(1)/result(0))
                     total_pe = data_pe


                     t3_oifits(ns).t3amp(tri,j) = result(0)
                     t3_oifits(ns).t3amp_err(tri,j) = abs(result(0))*total_pe
                     t3_oifits(ns).badflag_t3amp(tri,j) = 0
                  endif else begin ; no t3amp data.
                     t3_oifits(ns).t3amp(tri,j) = !values.f_nan
                     t3_oifits(ns).t3amp_err(tri,j) = !values.f_nan
                     t3_oifits(ns).badflag_t3amp(tri,j) = 1
                     ;stop
                     ;set NULL
                  endelse


               endfor           ; j dimy

            endfor              ; nsplit
         endif                  ; nd>0
      endfor                    ; triangles


      ; I should average the closure amps here!! and put in infostring
      ; kind of tricky -- with all the badflagging etc.


      ;SAVE FILE
      ;p0=rstrpos(infile,'/')
      ;p1=rstrpos(infile,'_info.idlvar')
      ;rootname = strmid(infile,p0+1,p1-p0-1)
      p0=rstrpos(reducefile,'/')
      p1=rstrpos(reducefile,'_reduced.idlvar')
      rootname = strmid(reducefile,p0+1,p1-p0-1)


      ; Replace '(' or ')' with period.
      ; I might need to do search for illegal characters.

      rootname2='MIRC_L1.'+rootname+'.'+targetname+normlabel+'.'+identifier
      oifits_name=rootname2+'.oifits'

      info_name  =rootname2+'.info.txt'

      oifits_file = oifits_path + oifits_name
      info_file = oifits_path + info_name
      filelist(i) = oifits_name
      tempq=(vis2_oifits.badflag)
      tempin=where(tempq[0:4,4,*] ne 1,ct)
      print,oifits_name,'number of baseline 0-*: ',ct

      ;;stop
      mircfull2oidata, obs_date0,ra,dec, target,$
                       vis2_oifits, t3_oifits,$
                       specname, lambda,lambdaerr,telindex,file=oifits_file


   endfor                       ; i - looping through target blocks.
   ; Do Merging.
   ; Find unique

   allbin = where(headers.target eq 1 or headers.calib eq 1)
   tnames=headers(allbin).source

   ;JDM2015 somewhere here we can flag blocks that have the etalons.
   ; give new names

   ;ETALON section


   tuniq=tnames(uniq(tnames,sort(tnames)))
   for uqt=0,n_elements(tuniq)-1 do begin

      allin=where(tnames eq tuniq(uqt),ct)
      print,tnames[allin]

      infiles=oifits_path+filelist(allin)
      outfile=oifits_path+'MIRC_L1.'+rootname+'.'+tuniq(uqt)+'.ALL.'+normlabel+'.'+identifier+'.oifits'
      outps = oifits_path+'MIRC_L1.'+rootname+'.'+tuniq(uqt)+'.ALL.'+normlabel+'.'+identifier+'.ps'


      merge_oidata,infile=infiles,outfile=outfile
      ; now remove the individual files to reduce clutter in the directory.
      for d=0,n_elements(infiles)-1 do begin
         textstring = '\rm '+infiles[d]
         spawn,textstring
      endfor

     ; after this process, sometimes there are files with no data. I will do a check now
      read_oidata,outfile,a,b,c,d,e,f
      TOTALDATA=n_elements(d) +n_elements(e)+n_elements(f)
      if totaldata eq 0 then begin
         print,'no data in ' + outfile
         textstring = '\rm '+outfile
         spawn,textstring
      endif else begin          ; make ps file and move into position
      print,' in here'

         mirc_reduce_ps,oifits=outfile,outps=outps
         newdir=oifits_path+rootname+'_'+identifier
         spawn,'mkdir '+newdir
         spawn,'\mv '+outfile+' '+newdir+'/.'
         spawn,'\mv '+outps+  ' '+newdir+'/.'

      endelse
     ;now move  files into directory structure.

   endfor                       ; loop through all the uniq targets in the file.

endfor                          ; psn norm methods

if use_scatter eq 1 then begin
   use_scatter=0
   goto,skipahead_mark2
endif


print,'ALL FINISHED for REAL!'
print,'Congratulations!'

end
