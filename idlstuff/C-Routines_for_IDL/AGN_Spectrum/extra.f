c---
      subroutine alf(cbuf, lbuf, kp, ctok, ltok)
      character cbuf*(*), ctok*(*)
      integer   lbuf, kp, ltok
c---
c strip preceeding spaces, copy token from cbuf to ctok, and
c then leave kp pointing one character before the start of the
c next token.  note: if ltok=0 then ctok is not modified.
c---
c cbuf      i    character array to be parsed
c lbuf      i    the last valid character in cbuf
c kp+1      i/o  points at next character to be parsed
c ctok        o  the parsed token
c ltok        o  number of valid characters in ctok
c---
c [aft]
c---
      integer   itab
      parameter (itab=9)
      character chr
      integer   iquote
c---
      ltok=0
      call alfsks(cbuf,lbuf,kp)
      if(kp.ge.lbuf) goto 900
      kp=kp+1
      iquote=0
      if(cbuf(kp:kp).eq.'"') then
         kp=kp+1
         iquote=1
      end if
c---
  150 chr=cbuf(kp:kp)
      if(iquote.eq.0) then
         if(chr.eq.' ' ) goto 200
         if(chr.eq.char(itab)) goto 200
         if(chr.eq.',' ) goto 900
      else
         if(chr.eq.'"') then
            if(cbuf(kp+1:kp+1).ne.'"') goto 200
            kp=kp+1
         end if
      end if
      if(ltok.lt.len(ctok)) ltok=ltok+1
      ctok(ltok:ltok)=chr
      if(kp.ge.lbuf) goto 200
      kp=kp+1
      goto 150
c---
c- strip trailing spaces and tabs and up to one comma.
  200 call alfsks(cbuf,lbuf,kp)
      if(cbuf(kp+1:kp+1).eq.',') kp=kp+1
c---
  900 if(ltok.lt.len(ctok)) ctok(ltok+1:)=' '
      return
      end
c*********
      subroutine alfsks(cbuf, lbuf, kp)
      character cbuf*(*)
      integer   lbuf, kp
c---
c skip spaces and tabs.  upon return either cbuf(kp+1:kp+1).ne.' '
c or kp.ge.lbuf
c---
c cbuf    i
c lbuf    i
c kp      i/o
c---
c [aft]
c---
      integer   itab
      parameter (itab=9)
      character chr
c---
  100 if(kp.ge.lbuf) return
      chr=cbuf(kp+1:kp+1)
      if(chr.eq.' ' .or. chr.eq.char(itab)) then
         kp=kp+1
         goto 100
      end if
      return
      end

c*********
c*********
      subroutine edicom(cbuf, lbuf)
      character cbuf*(*)
      integer   lbuf
c*********
c---
c this command scans the %ed% string and executes any sub-commands
c found.
c---
c      include 'edicmn.inc'
c edicmn.inc
c
c iup..ieof are the user-defined control keys.
c iftype is the number returned by fortyp.
c icedit is <>0 if line editing has been switched on.
c
      integer   iup, idown, ileft, iright, ibeg, iend
      integer   iwrite, ierase, ieof, iftype, icedit
c
      common/edicmn/ iup, idown, ileft, iright, ibeg, iend,
     :          iwrite, ierase, ieof, iftype, icedit
      integer   icomm(9)
      equivalence (iup,icomm(1))
c
      character ctok*64
      character coms(9)*1
      integer   i, itmp, j, kp, ltok
      integer   ifirst, isave
      save      ifirst, isave
      data coms/'u','d','l','r','b','e','w','x','z'/
      data ifirst/1/
c---
      if(ifirst.ne.0) then
         ifirst=0
         icedit=0
         iup   =  2
         idown = 14
         iright=  6
         ileft =  4
         ibeg  =  8
         iend  =  5
         iwrite= 18
         ierase= 21
         ieof  = 26
         call fortyp(iftype)
      end if
      kp=0
c
  100 call alf(cbuf, lbuf, kp, ctok, ltok)
      if(ltok.le.0) goto 900
      call upc(ctok)
      if(ctok(1:2).eq.'of') then
c %ed% off
c         if(icedit.ne.0) call ttrset()
         icedit=0
      else if(ctok(1:2).eq.'on') then
c %ed% on
         if(iftype.ne.0) then
c            call ttinit()
            icedit=1
         else
            write(6,*) 'edicom--line editing has not been implemented.'
         end if
      else if(ctok(1:3).eq.'res') then
c %ed% res
         if(isave.ne.0) then
c            call ttinit()
            icedit=1
         end if
      else if(ctok(1:3).eq.'sav') then
c %ed% sav
         isave=icedit
c         if(icedit.ne.0) call ttrset()
         icedit=0
      else if(ctok(1:3).eq.'ver') then
c %ed% ver
         write(6,*)  'edicom--version: 1990-mar-25'
      else if(ctok(1:1).eq.'?') then
c %ed% ?
         if(icedit.eq.0) then
            write(6,*) 'off'
         else
            write(6,151) (coms(j),char(icomm(j)+64),j=1,9)
  151       format(1x,'on ',20(1x,a,a))
         end if
      else if(ltok.eq.2) then
c assignment
         if(ctok(2:2).lt.'a' .or. ctok(2:2).gt.'z') goto 100
         itmp=ichar(ctok(2:2))-ichar('a')+1
         do 180 i=1,9
            if(ctok(1:1).eq.coms(i)) then
               do 170 j=1,9
                  if(j.eq.i) goto 170
                  if(icomm(j).ne.itmp) goto 170
                  write(6,*) 'edicom--error, ^',ctok(2:2),
     :               ' was used in ',coms(j),ctok(2:2),'.'
                  goto 100
  170          continue
               icomm(i)=itmp
               goto 100
            end if
  180    continue
      end if
      goto 100
c---
  900 continue
      return
      end


c*********
      subroutine conc(cbuf)
      character cbuf*(*)
c---
c convert case of cbuf to default case of system.
c for unix this is lower case.
c---
c cbuf    i/o  the file name to be converted.
c---
c      call locase(cbuf)
      return
      end
c*********

c*********
      subroutine upc(cstr)
      character cstr*(*)
c---
c converts character array in cstr to upper case.
c---
c cstr    i/o  the string to convert to upper case
c---
      integer   i, icstr
c---
      do 120 i= 1, len(cstr)
         icstr= ichar(cstr(i:i))
         if(icstr.gt.96 .and. icstr.lt.123) cstr(i:i)= char(icstr-32)
  120 continue
      return
      end
c*********
c*********
      subroutine openwr(lun,cfile,cstat,cacc,car,lrecl,iread,ier)
      character cfile*(*), cstat*(*), cacc*(*), car*(*)
      integer   lun, lrecl, iread, ier
c---
c wrapup routine for the fortran open statement.  this version will
c force the filename to be case insensitive, unless the the first
c character of cfile is a backslash.  if the first character is a
c backslash, then it should not be considered part of the filename
c on any system.
c   cfile can be of the form 'disk,dir/sub,file' in which case, openwr
c will create 3 strings and call ptend to construct the (system-
c dependent) file name.  this provides a system independent way
c to specify file names.
c---
c lun     i    the logical unit number to be opened
c cfile   i/o  the name of the file
c cstat   i    ='old', 'new' or 'unknown'
c cacc    i    ='d' for direct access, 'u' for unformatted, 'e' for append
c car     i    carriage control, 'l' for 'list' or 'f'
c lrecl   i    record length (required for direct access files)
c iread   i    <>0 for readonly, shared
c ier       o  <>0 if file not opened
c---
c 1989-aug-10 - do not call conc if file name starts with backslash [aft]
c 1989-jul-07 - allow file names of form 'disk,dir/sub,file' [aft]
c 1989-feb-13 - latest mod [aft]
c---
      integer   lenact
c
      character ctmp*72,tmpstr*72
      character cdisk*32, cdir*32
      character cdum
      integer   ifile, ios, is, kp, lfile, ltmp
c---
   11 format(a)
c---
c      write (6,*)'in openwr',lun,cfile,cstat,cacc,car,lrecl,iread,ier
c
      lfile=lenact(cfile)
c---
c special treatment for terminal device.
c---
      if(lfile.eq.0) then
c needs additional check to prevent trashing a terminal
         if(car(1:1).eq.'l'.or.car(1:1).eq.'l') then
            open(unit=lun, file='/dev/tty', status='old',
     :        iostat=ios)
            goto 140
         end if
      end if
c---
c sort out file type.
c  0 = sequential access, unformatted         " " only
c  1 = direct access => unformatted           "d" + anything
c  2 = sequential access, unformatted         "u" only
c  3 = sequential access, formatted, append   "e" only
c  4 = sequential access, unformatted, append "ue" only
c
      ctmp=cacc
      call upc(ctmp)
      if(index(ctmp,'D').ne.0) then
         ifile=1
      else if(index(ctmp,'U').ne.0 .and. index(ctmp,'e').ne.0) then
         ifile=4
      else if(index(ctmp,'U').ne.0) then
         ifile=2
      else if(index(ctmp,'E').ne.0) then
         ifile=3
      else
         ifile=0
      end if
c---
c check for a "+" as the first character of the first name in which case
c force an append
      if(cfile(1:1) .eq. '+') then
         if(ifile.eq.0) ifile=3
         if(ifile.eq.2) ifile=4
         cfile=cfile(2:)
         lfile=lfile-1
      end if
c
      kp=0
      call alf(cfile, lfile, kp, ctmp, ltmp)
      is=1
      if(ltmp.lt.lfile) then
c discovered a 'system independent' file specification.  decode.
         cdisk=ctmp
         call alf(cfile, lfile, kp, ctmp, ltmp)
         cdir=ctmp
         call alf(cfile, lfile, kp, ctmp, ltmp)
         cfile=ctmp
         call ptend(cdisk, cdir, cfile)
         lfile=lenact(cfile)
         call conc(cfile)
      else
c if first chacter is a backslash, then ignore it
         if(cfile(1:1).eq.char(92)) then
            is=2
         else
c otherwise convert case
            call conc(cfile)
         end if
      end if
c
c      write (6,*)'ine openwr, ifile=',ifile
c---
  130 continue
      if(ifile.eq.1) then
         open(unit=lun, file=cfile(is:lfile), status=cstat,
     :      access='direct', recl=lrecl, iostat=ios)
c         write (6,*)'direct access',is,lfile,
c     $                cfile(is:lfile),cstat,lrecl,ios
      else if(ifile.eq.2) then
         open(unit=lun, file=cfile(is:lfile), status=cstat,
     :      form='unformatted',iostat=ios)
      else if(ifile.eq.3) then
         open(unit=lun, file=cfile(is:lfile), status=cstat,
     :      iostat=ios)
c     :      fileopt='eof',iostat=ios)
      else if(ifile.eq.4) then
         open(unit=lun, file=cfile(is:lfile), status=cstat,
     :      form='unformatted',iostat=ios)
c     :      fileopt='eof',form='unformatted',iostat=ios)
      else
         open(unit=lun, file=cfile(is:lfile), status=cstat,
     :      iostat=ios)
      end if
c
  140 if(ios.eq.0) then
         inquire(unit=lun,name=cfile)
      else
c special treatment under unix, to allow user to overwrite an
c existing file
         if(ios.eq.117) then
            call edicom('sav', 3)
            tmpstr=cfile(is:lfile)//' exists, reuse it?'
            call prompt(tmpstr,0)
            read(*,11,err=900,end=900) cdum
            call edicom('res', 3)
            if(cdum.eq.'y' .or. cdum.eq.'y') then
               open(unit=lun, file=cfile(is:lfile), status='unknown')
               close(unit=lun, status='delete')
               goto 130
            end if
         end if
      end if   
  900 ier=ios
      return
      end

c********
      subroutine prompt(cbuf, lbuf)
      character cbuf*(*)
      integer   lbuf
c---
c write cbuf to the terminal leaving the cursor at the end of cbuf
c when finished.
c---
c cbuf    i    the prompt string
c lbuf    i    the number of valid characters in cbuf (can be zero)
c---
      integer   lenact
      integer   lb
c---
      lb= lbuf
      if(lb.eq.0) lb= lenact(cbuf)
      if(lb.gt.0) then
         write(*,111) cbuf(:lb)
  111    format(1x, a, ' ', $)
      else
c due to a sun funny, writing nothing to the line should space the
c cursor over one space, to allign with other text.
         write(*,121)
  121    format(' ',$)
      end if
      return
      end


c********
      subroutine ptend(cdisk, cdir, cfile)
      character cdisk*(*), cdir*(*), cfile*(*)
c---
c prefix extension.  add 'prefix' to file name.
c---
c cdisk     i    the disk name
c cdir      i    the directory name
c cfile     i/o  the file name
c---
c  2-aug-1988 - [aft]
c---
      integer   lenact
      integer   i, ldisk, ldir, lfile, lout
c---
      ldisk=lenact(cdisk)+2
      ldir =lenact(cdir)+1
      lfile=lenact(cfile)
      lout=ldisk+ldir+lfile
      if(lout.gt.len(cfile)) then
         write(*,*) 'ptend error, not enough room to create filename.'
         stop
      end if
c---
c- move characters in cfile
      do 190 i=0,lfile-1
        cfile(lout-i:lout-i)=cfile(lfile-i:lfile-i)
  190 continue
c---
c- finish the job
      cfile(1:ldisk)='/'//cdisk(:ldisk-2)//'/'
      cfile(ldisk+1:ldisk+ldir)=cdir(:ldir-1)//'/'
      call conc(cfile)
      return
      end
c*********
      subroutine fortyp(iftype)
      integer   iftype
c---
c iftype should return
c  0  if single character io (rdchr, putstr) has not been implemented.
c -1  if fortran i/o preceeds write statements with a cr//lf.
c +1  if cr//lf follow fortran write operations.
c---
      iftype=+1
      return
      end
c*********
      subroutine frelun(lun)
      integer   lun
c---  
c release a logical unit number.
c---
c lun       i    the logical unit number
c---
c 1989-feb-13 - latest mod [aft]
c---  
      return
      end
c*********
c*********
      subroutine getlun(lun)
      integer   lun
c---
c returns a free logical unit number.  note, units 1-9 should be
c avoided, since many programs use these units without warning.
c---
c lun       o  an unopened logical unit number.
c---
      integer   i
      logical   qopen
c---
      do 180 i= 99, 10, -1
         inquire(unit= i,  opened= qopen)
         if(.not.qopen) then
            lun= i
            return
         end if
  180 continue
      write(*,*) 'getlun, error, out of units.'
      return
      end
c*********

      subroutine xwrite(cstr, nonimp)
      character cstr*(*)
      integer   nonimp
c entry xwrtpr
      character cout*(*)
c---
c xparse subroutine to write a string on the terminal and/or the log
c file, depending on the values of the chattyness flags
c---
c cstr    i    string to be written.  n.b. the first character is
c              assumed to contain fortran carriage control infor-
c              mation
c nonimp  i    the unimportance of the string.  if =0, the string
c              is never written, else if 1 < nonimp < trmcht (terminal
c              chattyness) the string is written to the terminal
c              and if 1 < abs(nonimp) < logcht the string is
c              appended to the log file (if currently open).
c---
c 1986-mar-08 - rashafer
c---
c     include 'xparinc.inc'
c xparinc.inc
c include file for xparse to define the xparsc common block
c
c xprmsg   c256 message work space, used for input and output
c               -of messages
c xprmsw   c256 additional string work space
c ixpioe   i4   last io status value
c qxcase   l4   if ture, case is not significant in matches
c qxpart   l4   if true, the entire string need not be matched
c qxpfir   l4   if ture, the first partial match is accepted, without
c               -raising the error condition
c xpreof   c4   special eof string
c lxpeof   i4   lenact of xpreof
c blank    c1   blank
c commnt   c1   comment character
c opnstr   c1   open string character
c clsstr   c1   close string character
c comma    c1   ordinary break character
c skpall   c1   special break character (end of parse)
c spcbr1   c1   special break 1
c spcbr2   c1   special break 2
c indrct   c1   indirect file indicator
c contin   c1   continuation char
c tab      c1   tab char
c command  c1   command char
c opsys    c1   operating systems char
c inquiry  c1   the inquiry char (for help purposes)
c xxxspr   c2   spare
c lgunit   i4   no longer used, should be deleted next time file is modified.
c trmcht   i4   the terminal chattyness
c logcht   i4   the logfile chattyness
c request_inquiry  c4   the request prompt
c require_inquiry  l4   if true, then prompts are required
c return_inquiry   l4   give prompts on the current line
c
      character*512 xprmsg,xprmsw
      character*4 xpreof, request_inquiry
      character*2 xxxspr
      character*1 blank,commnt,opnstr,clsstr,comma,skpall,spcbr1,
     &      spcbr2,indrct,contin,tab,command,opsys,inquiry
      integer   ixpioe,lxpeof,lgunit,trmcht,logcht
      logical   qxcase,qxpart,qxpfir,require_inquiry,return_inquiry
      common/xparsc/lxpeof,lgunit,trmcht,logcht,
     &      xprmsg,xprmsw,ixpioe,qxcase,qxpart,qxpfir,xpreof,
     &      blank,commnt,opnstr,clsstr,comma,skpall,spcbr1,
     &      spcbr2,indrct,contin,tab,command,opsys,inquiry,xxxspr,
     &      request_inquiry,require_inquiry,return_inquiry
      save   /xparsc/
c
      integer   lenact
c
      real      rbuf
      integer   lout, lstr, nbuf
c---
   11 format(a)
c---
      if(nonimp.eq.0) return
      lstr=max(1,lenact(cstr))
      if((nonimp.gt.0).and.(nonimp.le.trmcht))
     :      write(*,11) cstr(:lstr)
      if(abs(nonimp).le.logcht)
     :    call logger(5, rbuf, nbuf, cstr, lstr)
      return
c*********
      entry xwrtpr(cout)
c---
c xwrtpr is used to write a string to the terminal/log file in
c exactly the same way as for a prompt, but without requesting
c any input.  useful for creating multiple line prompts.
c---
      lout=lenact(cout)
      write(*,11) cout(:max(1,lout))
      call logger(5, rbuf, nbuf, cout, lout)
      return
      end
      subroutine logger(ifunc, rbuf, nbuf, cbuf, lbuf)
      integer   ifunc, nbuf, lbuf
      real      rbuf(*)
      character cbuf*(*)
c---
c device driver for the log file.
c---
c 1990-sep-16 - [aft]
c---
      integer   lenact
c
      character cdnam*64, ctmp*64
      character cacc*10
      character cstat*8
      integer   ier, itmp, lunloc
      save      lunloc, cdnam
      data      lunloc/0/, cdnam/'log.log'/
c---
   11 format(a)
c---
      goto( 10, 20, 30, 40, 50) ifunc
  900 write(*,901) ifunc
  901 format('unimplemented function in logger device driver:',i5)
      nbuf = -1
      return
c
c--- ifunc = 1, set default file name-----------------------------------
c
   10 continue
      if(lenact(cbuf).ne.0) then
         cdnam=cbuf
      end if
      return
c
c--- ifunc = 2, get default name and unit number -----------------------
c
   20 continue
      rbuf(1)=lunloc
      nbuf=1
      if(lunloc.ne.0) then
         cbuf=cdnam
      else
         inquire(unit=lunloc,name=cbuf)
      end if
      lbuf=lenact(cdnam)
      return
c
c--- ifunc = 3, open the file ------------------------------------------
c
   30 continue
      if(lbuf.gt.0) then
         ctmp=cbuf
         call xtend(ctmp,'log')
      else
         ctmp=cdnam
      end if
c
c if already open, close current file, and open a new version.
      if(lunloc.ne.0) close(unit=lunloc)
c user can force the selection of a logical unit number.
      lunloc=rbuf(1)
      if(lunloc.ne.0) then
         close(unit=lunloc)
      else
         call getlun(lunloc)
      end if
c
      if(rbuf(2).ne.0.0) then
         cacc='e'
         cstat='unknown'
      else
         cacc=' '
         cstat='new'
      end if
c
      call openwr(lunloc, ctmp, cstat, cacc, 'l', 0, 0, ier)
      if(ier.eq.0) then
         nbuf=0
         rbuf(1)=lunloc
      else
         nbuf=-1
         rbuf(1)=0.0
      end if
      return
c
c--- ifunc = 4, close the file -----------------------------------------
c
   40 continue
      if(lunloc.ne.0) then
         close(unit=lunloc, status=cbuf)
         call frelun(lunloc)
         lunloc=0
      end if
      return
c
c--- ifunc = 5, write a string to the file (if active) -----------------
c
   50 continue
      if(lunloc.ne.0) then
         itmp=lenact(cbuf)
         write(lunloc,11) cbuf(:itmp)
      end if
      return
      end
      subroutine xtend(cnam, cext)
      character cnam*(*), cext*(*)
c---
c if cext contains a dot then this routine will force cext to be the
c extension.  if cext does not contain a dot then then cext will be
c appended only if cnam does not already have an extension.  for
c example, if cnam='fun' and cext='dat' then on output cnam='fun.dat'.
c---
c cnam    i/o  the file name.
c cext      o  the extension.
c---
c 1989-feb-13 - don't xtend zero length file names [aft]
c---
      integer   lenact
c
      character ctmp*10
      integer   idot, iend, ista, lext, lnam
c---
      lnam=lenact(cnam)
      if(lnam.eq.0) return
      lext=lenact(cext)
      ctmp=cext
      call conc(ctmp)
c---
c find first dot after the directory
      call dirpos(cnam, ista, iend)
      idot=iend+index(cnam(iend+1:),'.')
c exit if file name contains an extension, and requested extension is
c not forcing.
      if(idot.gt.iend .and. ctmp(1:1).ne.'.') return
      if(idot.eq.iend) then
         idot=lnam+1
         cnam(idot:idot)='.'
      end if
      if(ctmp(1:1).eq.'.') then
         cnam(idot+1:idot+lext-1)=ctmp(2:lext)
         cnam(idot+lext:)=' '
      else if(ctmp(1:1).ne.'.') then
         cnam(idot+1:idot+lext)=ctmp(:lext)
         cnam(idot+lext+1:)=' '
      end if
      return
      end
c*********
      subroutine dirpos(cfile, ista, iend)
      character cfile*(*)
      integer   ista, iend
c---
c cfile is a full filespec containing the disk and directory names.
c upon return cfile(ista:iend) contains only the directory spec.
c---
c cfile   i    the full file name
c ista      o  first valid character in directory spec, can be zero
c iend      o  last valid character in directory spec, can be zero
c---
c 1989-jul-08 - [aft]
c---
      integer   lenact
      integer   i, lfile
c---
      lfile=lenact(cfile)
      iend=0
      do 190 i=1,lfile
         if(cfile(i:i).eq.',') iend=i
         if(cfile(i:i).eq.'/') iend=i
  190 continue
      ista=min(1,iend)
      return
      end
