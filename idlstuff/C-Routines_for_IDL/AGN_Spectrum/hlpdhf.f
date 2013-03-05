      program hlpdhf
c---
c xhelp utility program to convert a sequential help format file
c to the direct access form used by the actual help routines.
c the program accomplishes this through two passes through the .hlp
c file, the first to generate the sizes and addresses of the various
c topics, the second to actually write out the dhf file.
c    this version can read the same files that can be installed into
c a vms help library.  all rules for creating a vms .hlp file should
c be followed, with the following restriction:
c the first line of the .hlp file must contain the "1 help" topic,
c followed by the top level message.  this is because the current
c version of gtxhlp assumes that the first topic in the .dhf file
c contains the top-level message.
c---
c 30-may-1985 - rashafer
c  3-aug-1988 - modified to read .hlp files [aft]
c---
      integer mxchr, mxtop, mxbuf
      parameter (mxchr=20000)
      parameter (mxtop=1000)
      parameter (mxbuf=512)
      integer lenact, newtop
c
      character cpics*(mxchr)
      character cline*256
      character ctok*64
      character cin*40, cout*40, cscr*12
      character cbuf*(mxbuf), ctmp*1
      character c4tmp*4
      integer i4tmp
      equivalence (c4tmp, i4tmp)
      integer ichstr(mxtop), ichstp(mxtop), inisub(mxtop), iprtop(mxtop)
     &        , inxtop(mxtop), itxsiz(mxtop), islsiz(mxtop)
      integer i, icssiz, icsize, ictsiz, idepth, ieqtop, ier
      integer irhelp, ios, ipage, iptop, irec, isubt, itlen, itop
      integer jdepth, jtop, kp, lbuf, lin, lout
      integer lline, lscr, ltok, ltops, ntop
      logical qequiv(mxtop)
c---
c** ichstr   first char of topic name
c** ichstp   pointer to last char of topic name
c** inisub   topic number of first sub-topic or if qequiv is true,
c      indicating that the topic has been equivalenced to
c      a previously defined topic, it is the topic number
c      of that topic.
c** iprtop   topic number of the topic that this topic belongs to.
c      the root topic is topic number 1, which has by
c      definition iprtop= 0.
c** inxtop   topic number of the next topic in the same level.  if
c      zero, then there is currently no next topic at this
c      level.
c** itxsiz   total no. of bytes in the topic text (including the
c      byte for the string lengths)
c** islsiz   total no. of bytes in the selection page, including
c      the byte for the string lengths and addresses.
c
c---
 11   format (a)
c---
c** first get the name of the input file
 100  call gtbuf('input .hlp file, [output .dhf file]:', ier)
      if (ier.lt.0) then
         goto 900
      endif
      call gtchar(cin, lin)
      call gtchar(cout, lout)
      if (lout.le.0) cout = cin
      call xtend(cin, 'hlp')
      open (unit=1, file=cin, status='old', iostat=ios)
      if (ios.ne.0) then
         write (6, *) 'unable to open sequential file ',
     &                cin(:lenact(cin))
         goto 100
      endif
c** create output file name from optional second argument or from
c** first argument
      call xtend(cout, 'dhf')
      call openwr(2, cout, 'new', 'd', ' ', mxbuf, 0, ios)
      if (ios.ne.0) then
         write (6, *) 'unable to open direct access output file ',
     &                cout(:lenact(cout))
         goto 100
      endif
      irhelp=0
c---
c- first scan the file, making the index
      itop = 0
      idepth = 0
      ntop = 0
      ltops = 0
      itxsiz(1) = 0
      irec = 0
 120  continue
      read (1, 11, iostat=ios) cline
      if (ios.ne.0) then
         goto 200
      endif
      irec = irec + 1
      if (newtop(cline).eq.0) then
c- repeat until next topic is found (counting bytes)
         if (ntop.eq.0) then
            write (6, *) 'error--file must begin with a topic.'
            stop
         endif
         itxsiz(ntop) = itxsiz(ntop) + max(1, lenact(cline)) + 1
         goto 120
      endif
c** the end of a block's text
      jtop = ntop + 1
      if (jtop.gt.mxtop) then
         write (6, *) '*error*: exceeded max no. of topics'
         goto 900
      endif
      itxsiz(jtop) = 0
c** find out the depth of the block.  for level 1, the topic 'help'
c** generates the xhelp level 0 text.
      jdepth = ichar(cline(1:1)) - 48
      lscr = lenact(cline)
      if (jdepth.eq.1 .and. lscr.lt.12) then
         cscr = cline
         call upc(cscr)
         if (lscr.gt.4) then
         if (cscr(lscr-4:lscr).eq.' help') then
            if(irec.eq.1) then
               irhelp=1
               jdepth = jdepth - 1
            else
               if(irhelp.eq.0) then
                  write (6, *) '"1 help" must be first line in file.'
                  stop
               end if
            endif
         endif
         endif
      endif
      if (ntop.le.0) then
c- first time, set initial depth.
         idepth = jdepth
      else
         if (jdepth.gt.idepth) then
            idepth = idepth + 1
            if (jdepth.gt.idepth) then
               write (6, *) 'topic ', jtop, ' has incorrect nesting'
            endif
c** the depth has increased by one, so check if the last topic was
c** not an equivalenced one
            if (qequiv(ntop)) then
               write (6, *) 'an equivalenced topic can not have a',
     &                      ' sub-topic'
               goto 900
            endif
            inisub(ntop) = jtop
            iprtop(jtop) = ntop
            inxtop(jtop) = 0
         elseif (jdepth.eq.idepth) then
c** the same depth as the last topic
            inxtop(ntop) = jtop
            iprtop(jtop) = iprtop(ntop)
            inxtop(jtop) = 0
         else
c** the topic depth must be reduced
            do 160 i = idepth - 1, jdepth, -1
               ntop = iprtop(ntop)
 160        continue
            idepth = jdepth
            inxtop(ntop) = jtop
            iprtop(jtop) = iprtop(ntop)
            inxtop(jtop) = 0
         endif
      endif
c** jtop is currently the most recent topic, while ntop is the topic
c** that it is subsidiary to
      inxtop(jtop) = 0
      inisub(jtop) = 0
      lline = lenact(cline)
      kp = 1
      call alf(cline, lline, kp, ctok, ltok)
c** check for an equivalenced topic
      qequiv(jtop) = .false.
      if (ltok.ge.2) then
         if (ctok(ltok-1:ltok).eq.'==') then
            call fndtop(ctok, ltok-2, ntop, ieqtop, cpics, 
     @                  ichstr, ichstp)
            if (ieqtop.ne.0) then
               qequiv(jtop) = .true.
               inisub(jtop) = ieqtop
               ltok = ltok - 2
            else
               write (6, *) 'unable to make an equivalence match'
            endif
         endif
      endif
c---
      ichstr(jtop) = ltops + 1
      ltops = ltops + ltok
      cpics(ichstr(jtop):ltops) = ctok(:ltok)
      ichstp(jtop) = ltops
      iptop = iprtop(jtop)
      if (iptop.gt.0) islsiz(iptop) = islsiz(iptop) + ltok + 5
      islsiz(jtop) = 0
      ntop = jtop
      goto 120
c---
 200  continue
c---
c- a few debug statements (to print the index)
c   write(*,*) '   i,               topic,star,stop,isub,prev,'//
c     :      'next,txsz,slsiz'
c   do 220 i=1,ntop
c     write(*,211) i,cpics(ichstr(i) :ichstp(i)), ichstr(i),
c     :      ichstp(i), inisub(i), iprtop(i), inxtop(i),
c     :       itxsiz(i), islsiz(i)
c211     format(i5,1x,a20,7i5,l2)
c220   continue
c---
c** having accumulated the information itxsiz will become the block
c** address of the topic, and islsiz will be the address of the
c** selection information.
      icsize = ichstp(1) + 1
      do 270 itop = 1, ntop
         ictsiz = itxsiz(itop) + 9
c** the term 9 is added to include 4 bytes for the default higher
c** level block add., the selection page add. (4 bytes), and a
c** single byte for the terminating text string (indicated by a
c** single length byte= 0)
         if ((.not.qequiv(itop)) .and. (inisub(itop).ne.0)) then
            icssiz = islsiz(itop) + 1
         else
            icssiz = 0
         endif
         itxsiz(itop) = icsize
         icsize = icsize + ictsiz
         if (icssiz.ne.0) then
            islsiz(itop) = icsize
            icsize = icsize + icssiz
         else
            islsiz(itop) = 0
         endif
 270  continue
c---
c- rewind the input file, go through again copying the text to output.
      rewind (1)
      ctmp = char(ichstp(1))
      lbuf = 0
      ipage = 0
      call chtobl(ctmp, 1, cbuf, lbuf, ipage, mxbuf, 2)
      call chtobl(cpics(ichstr(1):), ichstp(1), cbuf, lbuf, ipage,
     &            mxbuf, 2)
      itop = 0
      ios = 0
 290  if (ios.eq.0) then
         read (1, 11, iostat=ios) cline
         if (itop.eq.0 .or. newtop(cline).eq.1 .or. ios.ne.0) then
c** a new topic has begun
            if (itop.gt.0) then
c** close out the previous topic's text page(s)
               ctmp = char(0)
               call chtobl(ctmp, 1, cbuf, lbuf, ipage, mxbuf, 2)
               if (islsiz(itop).gt.0) then
c** write out previous topic's selection page
                  isubt = inisub(itop)
 310              if (isubt.gt.0) then
                     itlen = ichstp(isubt) - ichstr(isubt) + 1
                     ctmp = char(itlen)
                     call chtobl(ctmp, 1, cbuf, lbuf, ipage, mxbuf, 2)
                     call chtobl(cpics(ichstr(isubt):), itlen, cbuf,
     &                           lbuf, ipage, mxbuf, 2)
                     if (qequiv(isubt)) then
                        i4tmp = itxsiz(inisub(isubt))
                        call chtobl(c4tmp, 4, cbuf, lbuf, ipage, mxbuf,
     &                              2)
                     else
                        i4tmp = itxsiz(isubt)
                        call chtobl(c4tmp, 4, cbuf, lbuf, ipage, mxbuf,
     &                              2)
                     endif
                     isubt = inxtop(isubt)
                     goto 310
                  endif
                  ctmp = char(0)
                  call chtobl(ctmp, 1, cbuf, lbuf, ipage, mxbuf, 2)
               endif
            endif
c
            if (ios.eq.0) then
c** the old topic was terminated by reading in a new topic card
               itop = itop + 1
c** the beginning of a topic has the default previous topic and the
c** select page's address
               if (iprtop(itop).eq.0) then
                  i4tmp = 0
                  call chtobl(c4tmp, 4, cbuf, lbuf, ipage, mxbuf, 2)
               else
                  i4tmp = itxsiz(iprtop(itop))
                  call chtobl(c4tmp, 4, cbuf, lbuf, ipage, mxbuf, 2)
               endif
               i4tmp = islsiz(itop)
               call chtobl(c4tmp, 4, cbuf, lbuf, ipage, mxbuf, 2)
            endif
         else
c** an ordinary addition to the current topic text page
            itlen = max(lenact(cline), 1)
            ctmp = char(itlen)
            call chtobl(ctmp, 1, cbuf, lbuf, ipage, mxbuf, 2)
            call chtobl(cline, itlen, cbuf, lbuf, ipage, mxbuf, 2)
         endif
         goto 290
      endif
c---
c flush the buffer.
      call chtobl(cline, -1, cbuf, lbuf, ipage, mxbuf, 2)
 900  continue
      c4tmp = 'off'
      call edicom(c4tmp,3)
      end
c*********
      integer function newtop(cline)
      character cline*(*)
c---
c searches the the start of a new topic.  for vms style, the topic
c lines start with a single digit followed by a space.  in the future,
c the topic level may require two digits, in which case this code
c will to be enhanced.
c---
c 29-jul-1988 - [aft]
c---
      integer ix
c---
      newtop = 0
      if (cline(2:2).ne.' ') return
      ix = ichar(cline(1:1))
      if (49.le.ix .and. ix.le.57) newtop = 1
      return
      end
      subroutine gtxhlp(iunit, cfile, ctopic)
      integer iunit
      character cfile*(*), ctopic*(*)
c---
c main subroutine to manipulate a dhf format help file interactively.
c---
c iunit   i    unit for reading the help file
c cfile   i    name of help file
c ctopic  i/o  initial command string on input, on output the
c              -remainder of the prompt string when a @ special
c              -character is read.
c---
c  6-jul-1985 rashafer
c version 1.1: modified to allow return with update information in
c      the ctopic field when a @ special character is read.
c 29-jul-1988 - use gtbuf and some reformatting of output [aft]
c 11-jan-1989 fwjhaberl
c      some cosmetic changes to make it more similar to vms help.
c      at the lowest (sub-)topic level the prompt for the
c      next (sub-)topic asks you for (sub-)topics at this
c      lowest level still.
c---
c** the block size (size of a single page)
      integer mblsiz
      parameter (mblsiz=512)
c** the number of pages to keep in the internal buffer.
      integer mxpage, iclsiz, ilnsiz
      parameter (mxpage=5)
      parameter (iclsiz=13, ilnsiz=78)
c** mxtop = maximum depth of topics
      integer mxtop
      parameter (mxtop=32)
      integer lenact
c
      character cbuf*(mxpage*mblsiz)
      character ccom*256, chead*256, cout*256
      character ctok*64, ctok1*64
      character comch*1
      character ctmp*1
      integer itab(6+2*mxpage)
      integer icadd, icmode, ier, iend, ilen, imatch
      integer ipage, isel, ista, itmp
      integer jcadd, jchar, jcslnm, jmatch, jpage, kp
      integer lcom, lenstr, lentop, lentot, ltok
      integer mode, nchar, ntop
      logical qcom, qmatch, qndone
      logical qlow
c
      character c4tmp
      integer i4tmp
      equivalence (c4tmp, i4tmp)
c
      integer intpad(mxtop), inslad(mxtop), ichstr(mxtop), ichstp(mxtop)
     &        , icnxsl(mxtop), icslnm(mxtop)
c**
c** intpad   initial topic address
c** inslad   initial selection page address
c** ichstr   first char of topic name
c** ichstp   last car of topic name
c** icnxsl   next selection address for this topic
c** icslnm   current selection no. for this topic, i.e.,
c**          -icslnm(i) is the sub-topic number current for topic i.
c---
c- set gtbuf mode to 1 (ignore @, $ and ! characters)
      mode = 1
      call gtmode(mode)
c
      ier = 0
      call xsquez(cfile, ' ', 0, ilen)
      ccom = cfile
      call xtend(ccom, 'dhf')
      call opnche(iunit, ccom, itab, mxpage*mblsiz, .true., mblsiz, ier)
      if (ier.ne.0) then
         goto 900
      endif
      ntop = 1
      ccom = ctopic
      call rdche(itab, cbuf, ctmp, 0, 1, ier)
      lenstr = ichar(ctmp)
      call rdche(itab, cbuf, chead, 1, lenstr, ier)
      call dirpos(cfile, ista, iend)
      chead = cfile(iend+1:)
      itmp = index(chead, '.')
      if (itmp.gt.0) chead(itmp:) = ' '
      ichstr(1) = 1
      ichstp(1) = lenact(chead)
      icslnm(1) = 0
c---
c** icadd is the current address
c** icmode is the current mode
c** 0 - at the start of the text page
c** 1 - in the midst of the text page
c** 2 - beginning of the selection page
c** 3 - beginning of selection page (continuing from text)
c** 4 - continuing the selection page
c** 5 - at the end of the current topic (= eotext if no sel.)
c** 6 - move to the beginning of the 'next' topic.
      icmode = 0
      icadd = lenstr + 1
      icadd = icadd + 4
      call rdche(itab, cbuf, c4tmp, icadd, 4, ier)
      inslad(1) = i4tmp
      icadd = icadd + 4
      intpad(1) = icadd
c** come from for a new command line to process
 100  continue
      lcom = lenact(ccom)
      if (lcom.ne.0) then
c** some command line
         qcom = .true.
         kp = 0
 130     if ((qcom) .and. (kp.lt.lcom)) then
            kp = kp + 1
            comch = ccom(kp:kp)
            if (comch.eq.'?') then
c** check if the next char also ?, indicating that xhelp
c** documentation is needed
               if (ccom(kp+1:kp+1).eq.'?') then
                  kp = kp + 1
                  call shfdoc
               else
c** goto the selection page of current topic otherwise, popup
c** to prev. topic sel page
                  icadd = inslad(ntop)
                  if (icadd.eq.0) then
                     write (6, '(1x,2a)') 'no selections for ',
     &                            chead(ichstr(ntop):ichstp(ntop))
                     ntop = ntop - 1
                     icadd = inslad(ntop)
                  endif
                  icmode = 2
               endif
            elseif (comch.eq.'^') then
               if (ntop.gt.1) then
                  ntop = ntop - 1
                  icmode = 2
                  icadd = inslad(ntop)
               else
                  goto 900
               endif
            elseif (comch.eq.'<') then
               if (ntop.gt.1) then
c** move to the previous topic if possible
                  jcslnm = max(icslnm(ntop-1)-1, 1)
                  icslnm(ntop-1) = jcslnm
                  icadd = inslad(ntop-1)
c** skip over the intervening ones
                  do 160 isel = 1, jcslnm - 1
                     call rdche(itab, cbuf, ctmp, icadd, 1, ier)
                     lentop = ichar(ctmp)
                     icadd = icadd + 1
                     call rdche(itab, cbuf, cout, icadd, lentop, ier)
                     icadd = icadd + 4 + lentop
  160             continue
c** now install the new topic name
                  call rdche(itab, cbuf, ctmp, icadd, 1, ier)
                  lentop = ichar(ctmp)
                  icadd = icadd + 1
                  call rdche(itab, cbuf, chead(ichstr(ntop):), icadd,
     &                       lentop, ier)
                  ichstp(ntop) = ichstr(ntop) + lentop - 1
                  icadd = icadd + lentop
                  call rdche(itab, cbuf, c4tmp, icadd, 4, ier)
                  jcadd = i4tmp
                  icnxsl(ntop-1) = icadd + 4
                  icadd = jcadd + 4
                  call rdche(itab, cbuf, c4tmp, icadd, 4, ier)
                  inslad(ntop) = i4tmp
                  icadd = icadd + 4
                  intpad(ntop) = icadd
               else
                  icadd = intpad(ntop)
               endif
               icmode = 0
            elseif (comch.eq.'>') then
               if (ntop.gt.1) then
                  icadd = icnxsl(ntop-1)
                  call rdche(itab, cbuf, ctmp, icadd, 1, ier)
                  lentop = ichar(ctmp)
                  if (lentop.gt.0) then
                     icslnm(ntop-1) = icslnm(ntop-1) + 1
                     icadd = icadd + 1
                     call rdche(itab, cbuf, chead(ichstr(ntop):), icadd,
     &                          lentop, ier)
                     ichstp(ntop) = ichstr(ntop) + lentop - 1
                     icadd = icadd + lentop
                     call rdche(itab, cbuf, c4tmp, icadd, 4, ier)
                     jcadd = i4tmp
                     icnxsl(ntop-1) = icadd + 4
                     icadd = jcadd + 4
                     call rdche(itab, cbuf, c4tmp, icadd, 4, ier)
                     inslad(ntop) = i4tmp
                     icadd = icadd + 4
                     intpad(ntop) = icadd
                  else
                     icadd = intpad(ntop)
                  endif
               else
                  icadd = intpad(ntop)
               endif
               icmode = 0
            elseif (comch.eq.'/' .and. lcom.eq.1) then
c** skip to the next topic in the heirarchy
               if (inslad(ntop).ne.0) then
                  icadd = inslad(ntop)
                  call rdche(itab, cbuf, ctmp, icadd, 1, ier)
                  lentop = ichar(ctmp)
               else
                  lentop = 0
 220              if (lentop.eq.0) then
                     ntop = ntop - 1
                     if (ntop.le.0) then
                        lentop = -1
                     else
                        icadd = icnxsl(ntop)
                        call rdche(itab, cbuf, ctmp, icadd, 1, ier)
                        lentop = ichar(ctmp)
                     endif
                     goto 220
                  endif
               endif
c---
               if (ntop.le.0) then
c** have moved all the way to the end
                  ntop = 1
                  icadd = intpad(1)
               else
                  ntop = ntop + 1
                  icadd = icadd + 1
                  ichstr(ntop) = ichstp(ntop-1) + 1
                  call rdche(itab, cbuf, chead(ichstr(ntop):), icadd,
     &                       lentop, ier)
                  icadd = icadd + lentop
                  ichstp(ntop) = ichstp(ntop-1) + lentop
                  call rdche(itab, cbuf, c4tmp, icadd, 4, ier)
                  jcadd = i4tmp
                  icnxsl(ntop-1) = icadd + 4
                  icadd = jcadd + 4
                  call rdche(itab, cbuf, c4tmp, icadd, 4, ier)
                  inslad(ntop) = i4tmp
                  icadd = icadd + 4
                  intpad(ntop) = icadd
               endif
               icmode = 0
            elseif (comch.eq.'@') then
c** return, with the ctopic reset
               ctopic = ccom(kp+1:)
               goto 990
            else
               qcom = .false.
               kp = kp - 1
            endif
            goto 130
         endif
c---
 250     if (kp.lt.lcom) then
c** there is still stuff on the command string, which you treat as
c** sub-topic name
            if (inslad(ntop).eq.0) then
               kp = lcom + 1
               write (6, '(1x,2a)') 'no selections for ',
     &                      chead(ichstr(ntop):ichstp(ntop))
            else
               call alf(ccom, lcom, kp, ctok, ltok)
               ctok(ltok+1:) = ' '
               imatch = 0
               jcadd = inslad(ntop)
               jcslnm = 0
c** loop for finding matches with topic names
 120           continue
               call rdche(itab, cbuf, ctmp, jcadd, 1, ier)
               lentop = ichar(ctmp)
               if (lentop.eq.0) then
                  goto 140
               endif
               jcslnm = jcslnm + 1
               jcadd = jcadd + 1
               call rdche(itab, cbuf, cout, jcadd, lentop, ier)
               jcadd = jcadd + lentop
               call smatch(ctok, cout(:lentop), .true., .true., qmatch)
               if (qmatch) then
                  call rdche(itab, cbuf, c4tmp, jcadd, 4, ier)
                  jmatch = i4tmp
                  if ((lentop.ne.ltok) .and. imatch.ne.0) then
c** there is a cross match
                     if (imatch.gt.0) then
                        write (6, '(1x,a)') ' multiple choices:'
                        write (6, '(1x,a)') 
     &                     chead(ichstr(ntop):ichstp(ntop))
c** cancel the selector indicator
                        ntop = ntop - 1
                        imatch = -1
                     endif
                     write (6, '(1x,a)') cout(:lentop)
                  elseif (lentop.eq.ltok) then
c** an exact match
                     if (imatch.gt.0) then
                        ntop = ntop - 1
                     elseif (imatch.lt.0) then
                        write(*, '(1x,a)') '...oops, exact match'
                     endif
                     imatch = jmatch
                  else
                     imatch = jmatch
                  endif
c---
                  if (imatch.gt.0) then
c** install the name
                     nchar = ichstp(ntop)
                     icslnm(ntop) = jcslnm
                     icnxsl(ntop) = jcadd + 4
                     ntop = ntop + 1
                     ichstr(ntop) = nchar + 1
                     ichstp(ntop) = nchar + lentop
                     chead(ichstr(ntop):ichstp(ntop)) = cout(:lentop)
                     if (lentop.eq.ltok) then
                        goto 140
                     endif
                  endif
               endif
               jcadd = jcadd + 4
               goto 120
c** come from at the end of the loop
 140           continue
c** check if a match is found, otherwise insert the correction
c** and continue
               if (imatch.eq.0) then
                  write (6, '(1x,2a)') 'sorry, there is no sub-topic ',
     &                         ctok(:ltok)
                  ccom = '?'
                  goto 100
               elseif (imatch.lt.0) then
                  ctok = 'please re-enter:'
                  goto 280
               else
                  icadd = imatch
                  icmode = 0
                  icadd = icadd + 4
                  call rdche(itab, cbuf, c4tmp, icadd, 4, ier)
                  inslad(ntop) = i4tmp
                  icadd = icadd + 4
                  intpad(ntop) = icadd
               endif
            endif
            goto 250
         endif
      endif
c** now having processed the command line, output some help begin with
c** the header
      if (icmode.le.1) then
         if (icmode.eq.0) then
            call wrttop(ntop, chead, ichstr, ichstp, cout)
            ctok1 = cout
         endif
         ipage = 0
         qndone = .true.
 320     if (qndone) then
            call rdche(itab, cbuf, ctmp, icadd, 1, ier)
            lenstr = ichar(ctmp)
            if (lenstr.gt.0) then
               ipage = ipage + 1
               icadd = icadd + 1
               call rdche(itab, cbuf, cout, icadd, lenstr, ier)
               write (6, '(1x,a)') cout(:lenstr)
               qndone = ipage.lt.18
               icadd = icadd + lenstr
            else
               qndone = .false.
            endif
            goto 320
         endif
         if (lenstr.eq.0) then
c** the end of the topic text was reached
            if (inslad(ntop).gt.0) then
c** there are sub-tropic selections, so see if there is room for
c** sub-topic names
               jcadd = inslad(ntop)
               icadd = jcadd
               jchar = 0
               qndone = .true.
               jpage = 0
 350           if (qndone) then
                  call rdche(itab, cbuf, ctmp, jcadd, 1, ier)
                  lenstr = ichar(ctmp)
                  if (lenstr.gt.0) then
                     lentot = (((lenstr/iclsiz)+1)*iclsiz)
                     if (lentot+jchar.gt.ilnsiz) then
                        jpage = jpage + 1
                        qndone = jpage.lt.18
                        jchar = 0
                     endif
                     if (qndone) then
                        jcadd = jcadd + lenstr + 5
                        jchar = jchar + lentot
                     endif
                  else
                     qndone = .false.
                  endif
                  goto 350
               endif
c---
               write (6, *)
               if (ntop.gt.1) write (6, '(1x,a)')
     &                             '  additional information available:'
               write (6, *)
               if ((ipage+jpage+1.lt.18) .or. (jpage.ge.18)) then
                  icmode = 2
                  goto 400
               else
                  icmode = 3
               endif
            else
c** there were no-sub topics, so set the pop up get
               icmode = 5
            endif
         else
c** not the end
            icmode = 1
         endif
      endif
c---
      if ((icmode.eq.2) .or. (icmode.eq.3) .or. (icmode.eq.4)) then
         ipage = 0
         if (icmode.eq.2) then
            call wrttop(ntop, chead, ichstr, ichstp, cout)
            ctok1 = cout
c** come from when the selection page is to be written directly
c** after the sub topic text.
 400        continue
            icadd = inslad(ntop)
         endif
         qndone = .true.
         jchar = 0
         cout = ' '
 420     if (qndone) then
            call rdche(itab, cbuf, ctmp, icadd, 1, ier)
            lenstr = ichar(ctmp)
            if (lenstr.gt.0) then
               lentot = (((lenstr/iclsiz)+1)*iclsiz)
               if (lentot+jchar.gt.ilnsiz) then
c** write out the line
                  write (6, '(1x,a)') cout(:jchar)
                  ipage = ipage + 1
                  qndone = ipage.lt.18
                  jchar = 0
                  cout = ' '
               endif
               if (qndone) then
                  icadd = icadd + 1
                  call rdche(itab, cbuf, cout(jchar+1:), icadd, lenstr,
     &                       ier)
                  jchar = jchar + lentot
                  icadd = icadd + lenstr + 4
               endif
            else
               qndone = .false.
            endif
            goto 420
         endif
         if (lenstr.gt.0) then
c** there are more sub-topic selections
            icmode = 4
         else
            icmode = 5
         endif
         if (jchar.gt.0) then
c** write out the last line
            write (6, '(1x,a)') cout(:jchar)
         endif
      endif
c---
      if (icmode.eq.1 .or. icmode.eq.4) then
         ctok = 'press return to continue ... '
      elseif (icmode.eq.5) then
         write (6, *)
         ilen = lenact(ctok1)
         if (inslad(ntop).eq.0) then
c ** reached lowest sub-topic level
            ilen = ilen - (ichstp(ntop)-ichstr(ntop)+2)
            ntop = ntop - 1
            qlow = .true.
         else
            qlow = .false.
         endif
         if (ntop.eq.1) then
            ctok = ctok1(1:ilen)//' topic? '
         else
            ctok = ctok1(1:ilen)//' sub-topic? '
         endif
      else
         ctok = 'xhelp>'
      endif
c---
 280  call gtbuf(ctok, ier)
      if (ier.lt.0) then
         goto 900
      endif
      call gtrest(ccom, lcom)
      ccom(lcom+1:) = ' '
      if (icmode.eq.5) then
         if (lcom.le.0) then
            ccom(1:1) = '^'
            if (qlow) ntop = ntop + 1
         elseif (ccom.eq.'/' .or. ccom.eq.'>' .or. ccom.eq.'<' .or.
     &           ccom.eq.'??') then
            if (qlow) ntop = ntop + 1
         endif
      endif
      if (icmode.eq.1 .and. ccom.eq.'/') ccom = ' '
      goto 100
c** come from for the end of the line
 900  continue
      ctopic = ' '
 990  continue
      call gtmode(mode)
      return
      end
      subroutine opnche(iunit, cfile, itab, mxbuf, qro, lrecl, ierr)
      character cfile*(*)
      integer   iunit, mxbuf, lrecl, ierr
      integer   itab(mxbuf)
      logical   qro
c---
c dache subroutine to open the cached da file.
c---
c iunit   i    fortran lun to use
c cfile   i    name of file to open
c itab    i/r  the buffer used to hold the cache and associated info
c              -(see details below) equivalenced to the array itab
c mxbuf   i    the size of the buffer
c qro     i    if true, file is to be open readonly
c              -before they are overwritten in the cache (obsolete)
c lrecl   i    number of bytes per record in file to be opened
c ierr      r  error flag, if zero then no error
c              -if 1, then unable to open the da file
c              -2 - unable to determine the recordlength
c              -3 - mxbuf too small to hold a single cache page.
c---
c 28-jun-1985 - rashafer
c  5-oct-1986 - sun unix requires recl for direct access opens [kaa]
c 30-jul-1988 - add call to openwr for general vax/sun usage [aft]
c---
c****** the buffer is divided into a header of 20 bytes:
c ipsize        i4      the size of a page in bytes
c npage         i4      the no. of pages in the cache
c lastpg        i4      the last page accessed
c lastcp        i4      the cache ptr of the last page accessed
c junit         i4      the unit to use for read and writes
c ****** followed by the cache pointers
c icadd(npage)  i4      the page no. currently in the ith position of
c                       the cache
c icuse(npage)  i4      the use index.  the lower the absolute value
c                       of the no., the more recently it has been
c                       accessed.  if a page is empty, the value
c                       is 0.  if the page has been written into, the
c                       value is <0.
c****** followed by the cached pages
c---
 
      integer ios, lenact
      integer npage, ipage
 
      if (lrecl.le.0) then
         if (ierr.eq.0) write (6, *) 'opnche: illegal record length:',
     &                               lrecl
         ierr = 4
         return
      endif
c---
      close (iunit)
      if (qro) then
         call openwr(iunit, cfile, 'old', 'd', ' ', lrecl, 1, ios)
      else
         call openwr(iunit, cfile, 'old', 'd', ' ', lrecl, 0, ios)
      endif
c
      if (ios.ne.0) then
         if (ierr.eq.0) then
 160        write (6, 161) cfile(:lenact(cfile)), ios
 161        format (' opnche:  unable to open file ', a, ', iostat=',
     &              i5)
         endif
         ierr = 1
         return
      endif
c
      npage = mxbuf/lrecl
      if (npage.le.0) then
         if (ierr.eq.0) write (6, *) 'opnche: no room for one page of',
     &                               lrecl, ' bytes in buffer of size',
     &                               mxbuf
         ierr = 3
         return
      endif
      itab(1) = lrecl
      itab(2) = npage
      itab(3) = -1
      itab(4) = 0
      itab(5) = iunit
      do 190 ipage = 1, npage
         itab(5+ipage) = -1
         itab(5+npage+ipage) = 0
 190  continue
      ierr = 0
      return
      end
      subroutine chtobl(cray, lray, cbuf, lbuf, ipage, mxbuf, iunit)
      character cray*(*), cbuf*(*)
      integer lray, lbuf, ipage, mxbuf, iunit
c---
c shf subroutine the accumulates bytes into a buffer array and then
c write out sequential blocks into the indicated direct access file.
c---
c cray    i    string of bytes to be accumulated
c lray    i    no. of bytes to so accumulate, when < 0 then dump
c              -the buffer in any case
c cbuf    i/o  accumulation buffer
c lbuf    i/o  current location in buffer (must be set to zero before
c              -the first call to c1tobl
c ipage   i/o  current d.a. record in file
c mxbuf   i    size of buffer (bytyes)
c iunit   i    fortran unit no.  file must have been previously opened
c              -as a fortran unformatted direct access file.  be
c              -sure that the record length in the open statement
c              -was a number of full-words (4 byte quantitites).
c              -(thus mxbuf must be a multiple of 4)
c---
c 27 jun 1985 - rashafer
c 25 feb 1987 to do a record write when lray is a negative number.
c 10 mar 1987 to properly handle multiple buffer length records
c---
      integer i
c---
      if (lray.le.0) then
         if (lbuf.gt.0) then
c** non empty buffer so dump it
            ipage = ipage + 1
            write (iunit, rec=ipage) cbuf
            lbuf = 0
         endif
         return
      endif
c---
      do 190 i = 1, lray
         lbuf = lbuf + 1
         cbuf(lbuf:lbuf) = cray(i:i)
         if (lbuf.eq.mxbuf) then
            ipage = ipage + 1
            write (iunit, rec=ipage) cbuf
            lbuf = 0
         endif
  190 continue
      return
      end
      integer function ichbeg(itab, cbuf, iadd, nbyte, lstadd, ierr)
      integer   itab(*)
      character cbuf*(*)
      integer   iadd, nbyte, lstadd, ierr
c---
c dache function call to return the buffer address of a
c given da file address
c---
c ichbeg    r  the offset in the cache buffer that points
c              -to the first byte needed (byte 0 is the beginning
c              -of the buffer and the file)
c itab    i/r  see opnche for a description
c cbuf         -cache buffer
c iadd    i    byte file address
c nbyte   i    last byte needed
c lstadd    r  last buffer address actually in cache block in the
c              -needed range
c ierr      r  error flag, 4 - rdpag error
c---
c 28-jun-85 - rashafer
c---
      integer ipsize, npage
      integer maxoff, icpage, iccp, i, iunit, imax
      integer jmax, iuse, jpage, ipage, iusoff, icoff
      integer juse, jpgoff
c---
c** see opnche for layout of the buffer
      ipsize = itab(1)
      npage = itab(2)
      iusoff = 5 + npage
      icpage = (iadd/ipsize)
      icoff = iadd - icpage*ipsize
      maxoff = min(icoff+nbyte, ipsize) - 1
      if (icpage.eq.itab(3)) then
         iccp = itab(4)
      else
c** find if the current page is in the cache
         iccp = 0
         i = 0
 110     if ((iccp.eq.0) .and. (i.lt.npage)) then
            i = i + 1
            if (itab(5+i).eq.icpage) iccp = i
            goto 110
         endif
         if (iccp.eq.0) then
c** the requested page not in the cache, find free page
            iunit = itab(5)
            i = 0
            imax = 0
            jmax = 0
 150        if ((iccp.eq.0) .and. (i.lt.npage)) then
               i = i + 1
               iuse = abs(itab(i+iusoff))
               if (iuse.eq.0) then
c** a pristine page
                  iccp = i
               elseif (iuse.gt.imax) then
                  imax = iuse
                  jmax = i
               endif
               goto 150
            endif
            if (iccp.eq.0) then
               iccp = jmax
               if (iccp.le.0) iccp = 1
               if (itab(iusoff+iccp).lt.0) then
c** the cache page has been updated, so first write it back out
                  jpage = itab(5+iccp)
                  call wrtpag(iunit, cbuf, iccp, ipsize, jpage, ierr)
               endif
            endif
            do 180 ipage = 1, npage
               iuse = itab(ipage+iusoff)
               if (iuse.ne.0) itab(ipage+iusoff)
     &              = sign(abs(iuse)+1, iuse)
 180        continue
            itab(iccp+iusoff) = 1
            call rdpag(iunit, cbuf, iccp, ipsize, icpage, ierr)
            itab(5+iccp) = icpage
         else
c** mark some other page as more current
            iuse = abs(itab(iccp+iusoff))
            do 220 ipage = 1, npage
               juse = itab(ipage+iusoff)
               if ((abs(juse).lt.iuse) .and. (juse.ne.0)) then
                  itab(ipage+iusoff) = sign(abs(juse)+1, juse)
               endif
 220        continue
         endif
         itab(3) = icpage
         itab(4) = iccp
      endif
c** calculate the byte address of the 0th byte of the selected page
      jpgoff = ipsize*(iccp-1)
      ichbeg = jpgoff + icoff
      lstadd = jpgoff + maxoff
      return
      end
 
      subroutine fndtop(ctok, ltok, ntop, ieqtop, cpics, ichstr, ichstp)
      character ctok*(*), cpics*(*)
      integer ltok, ntop, ieqtop, ichstr(*), ichstp(*)
c---
c shf subroutine to find the topic corresponding to a given topic
c ctok, as input on unit one.  see shftod for argument details.
c---
c ctok    i    the topic to match
c ieqtop    r  the topic number of the indicated topic ctok
c              -if zero, then no match was found.
c---
c 21-dec-1988 - cleaned up [aft]
c  1-jun-1985 - rashafer
c---
      integer i
      logical qmatch
c---
      do 190 i = 1, ntop
         call smatch(ctok(1:ltok), cpics(ichstr(i):ichstp(i)),
     :      .false., .true., qmatch)
         if (qmatch) then
            ieqtop = i
            return
         endif
 190  continue
      ieqtop = 0
      return
      end
      subroutine rdche(itab, cbuf, cout, iadd, nbyte, ier)
      integer   itab(*)
      character cbuf*(*), cout*(*)
      integer   iadd, nbyte, ier
c---
c dache subroutine to read from a direct access file
c---
c itab    i/o  see opnche for description
c cbuf    i/o  cache buffer
c cout    o    returned values
c iadd    i    starting byte to read (starting with 0)
c nbyte   i    no. of cout needed
c ier     i/o  dache error flag
c---
c 28-jun-85 - rashafer
c  1-aug-88 - better use of character array [aft]
c---
      integer icur, jadd, jbyte
      integer jer, jbeg, jend, ichbeg, jtrans
c---
      ier = 0
      icur = 0
      jadd = iadd
      jbyte = nbyte
      jer = ier
c*** repeat
 110  continue
      jbeg = ichbeg(itab, cbuf, jadd, jbyte, jend, jer)
      jtrans = jend - jbeg + 1
      if (jer.ne.0) then
         ier = jer
         icur = nbyte
      else
         cout(icur+1:icur+jtrans) = cbuf(jbeg+1:jend+1)
         icur = icur + jtrans
      endif
      jadd = jadd + jtrans
      jbyte = jbyte - jtrans
c** until current number in cout equals total requested
      if (icur.lt.nbyte) then
         goto 110
      endif
      return
      end
      subroutine rdpag(iunit, buf, icpage, ipsize, ifpage, ier)
      integer   iunit, icpage, ipsize, ifpage, ier
      character*1 buf(ipsize, icpage)
c---
c dache subroutine to readin a page of a da file into the buffer.
c---
c 1985-jun-29 - rashafer
c---
      integer   ios, i
c---
      read (iunit, rec=ifpage+1, iostat=ios) 
     :      (buf(i,icpage), i=1, ipsize)
      if (ios.ne.0) then
         if (ier.eq.0) write (6, 101) ios, ifpage + 1, iunit
 101     format (' *error*:rdpag: read error ', i3, ' for record', i4,
     &           ' on unit', i3)
         ier = 4
      endif
      return
      end
      subroutine shfdoc
c---
c dache subroutine that writes some useful information to the
c default output unit
c---
c 24-jul-1985 - rashafer
c---
c       write(*,*)'    xhelp - a standard help facility:  brief notes:'
      write (6, *) 'at any prompt you can type the following:'
      write (6, *) '<ctrl-z>         - to exit from facility'
      write (6, *)
     &          'one or more <cr> - pop up one or more levels until you'
      write (6, *)
     &            '                   exit from facility (at the lowest'
      write (6, *) '                   level the first <cr> gives the '
      write (6, *)
     &           '                   (sub-)topic list at current level)'
      write (6, *) '?   - to see (sub-)topic list at current level'
      write (6, *)
     &          '/   - skip to next (sub-)topic (you can go through the'
      write (6, *) '      complete help file by repeating /)'
      write (6, *) '>   - skip to following topic at current level'
      write (6, *) '<   - go to prev. topic at current level'
      write (6, *) '^   - pop up one level'
      write (6, *) 'or any of the sub-topic names'
      return
      end
      subroutine wrtpag(iunit, c1buf, icpage, ipsize, ifpage, ierr)
      integer   iunit, icpage, ipsize, ifpage, ierr
      character*1 c1buf(ipsize, icpage)
c---
c dache subroutine to readin a page of a da file into the buffer
c---
c 1985-jun-29 - rashafer
c---
      integer i, ios
c---
      write (iunit, rec=ifpage+1, iostat=ios) 
     :         (c1buf(i,icpage), i=1, ipsize)
      if (ios.ne.0) then
         if (ierr.eq.0) write (6, 111) ios, ifpage + 1, iunit
 111     format (' *error*:wrtpag: write error', i4, ' for record', i4,
     &           ' on unit', i2)
         ierr = 5
      endif
      return
      end
      subroutine wrttop(ntop, chead, ichstr, ichstp, cout)
      character chead*(*), cout*(*)
      integer ntop, ichstr(*), ichstp(*)
c---
c dache facitility subroutine to write out the current topic chain
c---
c ntop      i    no of topic levels currently active
c chead     i    contains topic names
c ichstr    i    start character for the topic name
c ichstp    i    stop character for the topic name
c cout        o  work space for output
c---
c 21-jul-1985 - rashafer
c---
      integer lenact
c
      integer itop, nchar, lenc
c---
      cout = chead(ichstr(1):ichstp(1))
      nchar = ichstp(1) + 1
      write (6, *)
      do 130 itop = 2, ntop
         lenc = ichstp(itop) - ichstr(itop) + 1
         if (lenc+nchar.gt.78) then
            nchar = lenact(cout(:nchar-1))
            write (6, 111) cout(:nchar)
 111        format (1x, a)
            nchar = 0
            cout(:max0(nchar-1,1)) = ' '
         endif
         cout(nchar+1:nchar+lenc) = chead(ichstr(itop):ichstp(itop))
         nchar = nchar + lenc + 1
 130  continue
      if (nchar.gt.0) then
         nchar = lenact(cout(:nchar-1))
         write (6, 111) cout(:nchar)
      endif
      write (6, *)
      return
      end
      subroutine smatch(test, string, qpart, qcase, qmatch)
      character test*(*), string*(*)
      logical   qpart, qcase, qmatch
c---
c search for single string matches in xhelp files.
c---
c test    i    test string
c string  i    base string
c qpart   i    if true, partial matches (where test is shorter than
c                  string) is allowed, if false, exact matches are needed
c qcase   i    if true, then upper and lower case are not signific.
c qmatch    r  if true, a match was found, if false, then not
c---
c 30-may-1985 - rashafer
c---
      integer   lenact

      integer length, ilen, ictest, icstr
c---
      qmatch = .false.
      length = lenact(test)
      if (qpart) then
         if (length.gt.lenact(string)) return
      else
         if (length.ne.lenact(string)) return
      endif
      if (qcase) then
         do 190 ilen = 1, length
            ictest = ichar(test(ilen:ilen))
            if ( ictest.ge.ichar('a') .and. ictest.le.ichar('z'))
     :         ictest = ictest + ichar('a')-ichar('a')
            icstr = ichar(string(ilen:ilen))
            if ( icstr.ge.ichar('a') .and. icstr.le.ichar('z'))
     :         icstr = icstr + ichar('a')-ichar('a')
            if (icstr.ne.ictest) return
 190     continue
      else
         if (string(:length).ne.test(:length)) return
      endif
      qmatch = .true.
      return
      end
      subroutine xgtnum (string, iparse, nreq, desc, ndesc,
     &                 	 valmin, valmax, nrange, numtype, retval, 
     &                   nret, iflag,idelim, * , *, *)
c		rashafer 9 march 1986
c		modified from xgtstr 14 april 1986
c	xparse subroutine to peel off a requested number of arguements as
c	general numbers.
c	string	c*	i/r: parse string
c	iparse	c*	i/r: parse position
c	nreq	i4	i: no. of strings requested
c	desc	c*(ndesc)	i: description of the string requested
c	ndesc	i4	i: no. of descriptions passed
c	valmin	***	i: the minimum allowed value for the ith string
c	valmax	***	i: the maximum allowed value for the ith string
c	nrange	i4	i: the number of min-max ranges (if <=0 then no
c			range checking is performed)
c	numtype	i4	i: index indicating the type of value requested:
c			= 3 for integer.
c			= 4 for real
c	retval	***	r: the values actually picked up (no change
c				on skips)
c	nret	i4	r: no. of strings actually processed (where skips
c				and infinite skips are included as processed)
c				if nret not = nreq then there was a fall off
c				the end of the string, or a special delimeter
c				met (see input value of idelim)
c	iflag	i4	r: value of condition flag at last call to xgtarg (or
c				to xinfix for eofs)
c	idelim	i4	i/r: if <= -1  there is no checking for special
c			delimeters, else on return idelim contains the value
c			of the delimeter that triggered the return (if 0
c			then a comma, if 1 then special del 1, etc.)
c	alternate returns:
c 		*1	same as first alternate return of xgtarg (fell off
c				line) or eof
c		*2	same as second alternate return of xgtarg (infinite
c				skip)
c		*3	special delimeter (see idelim)
c
c **	** n.b. values in the above list indicated as type *** are actually
c **	** arrays of the kind indicated by numtype, although in this routine
c **	** they are declared as byte arrays.
      character*(*) string,desc(*)
      integer iparse,nreq,ndesc,nret,iflag,idelim,nrange,numtype
c      byte valmin(0:),valmax(0:),retval(0:)
      character valmin(*),valmax(*),retval(*)
c
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
c
      integer lenact, jdelim, ibeg, iend, destmp, ierr
      integer offset,curdesc,range,rangeoff
      logical qskip,isnum,inrange
      integer ntype
      parameter(ntype=4)
      integer multi(ntype)
      logical xtsti4,xpari4,xtstr4,xparr4
      data multi/-1,-1,4,4/

      if((numtype.le.0).or.(numtype.gt.ntype).or.
     &    (multi(ntype).le.0))then
         write(6,*)
     &	  '*error*:xparse:xgtnum: unsupported numerical type',numtype
         return
      end if
c      write (6,*)'in xgtnum'
      nret=0
      do while(nret.lt.nreq)
c         write (6,*)'nret=',nret
c         write (6,*)string,iparse,ibeg,iend,qskip,iflag,jdelim
         call xgtarg(string,iparse,ibeg,iend,qskip,iflag,jdelim,*910,
     &               *920)
c         write (6,*)string,iparse,ibeg,iend,qskip,iflag,jdelim
         offset=multi(numtype)*nret
         nret=nret+1
         curdesc = min(ndesc, nret)
         if(.not.qskip)then
            goto (100,100,130,140) numtype
100         continue
            return
c **        ** integer
130         continue
            isnum = xpari4(string(ibeg:iend), curdesc,
     &              desc(curdesc), retval(offset+1))
c            write (6,*)'integer--',isnum
c            write (6,*)xprsmg
            goto 200
c **        ** real
 140        continue
            isnum = xparr4(string(ibeg:iend), curdesc,
     &              desc(curdesc), retval(offset+1))
c            write (6,*)'real--',isnum
c            write (6,*)xprsmg
            goto 200
c **        ** end of the case (numtype)
200         continue
         else
            isnum=.true.
         end if
         if(isnum.and.(nrange.gt.0))then
            range = min(nret,nrange)
            rangeoff = (range-1)*multi(numtype)
            goto (200,200,230,240),numtype
c **        ** integer
230         inrange = xtsti4(retval(offset+1),valmin(rangeoff+1),
     &                valmax(rangeoff+1),curdesc, desc(curdesc))
c            write (6,*)'integer, inrange=',inrange
            isnum = inrange
            goto 300
c **        ** real
240         continue
            inrange = xtstr4(retval(offset+1),valmin(rangeoff+1),
     &                valmax(rangeoff+1), curdesc, desc(curdesc))
c            write (6,*)'real, inrange=',inrange
            isnum = inrange
            goto 300
c **        ** end of case (numtype)
300         continue
         else
            inrange = .true.
         end if
         if(.not.isnum)then
c            write (6,*)'not isnum'
            call xinfix(xprmsg(:lenact(xprmsg)+1),string,iparse,ierr)
            nret=nret-1
            if(ierr.ne.0)then
               if(.not.inrange) then
c                  write (6,*)'not in range'
                  goto (300,300,330,340),numtype
c **              ** i4
330               continue
                  call xseti4(valmin(rangeoff+1),retval(offset+1))
c                  write (6,*)'after xseti4--',
c     $                valmin(rangeoff+1),retval(offset+1)
                  goto 400
c **              ** r4
340               continue
                  call xsetr4(valmin(rangeoff+1),retval(offset+1))
c                  write (6,*)'after xsetr4--',
c     $                valmin(rangeoff+1),retval(offset+1)
                  goto 400
c **              ** end case numtype
400               continue
               end if
               iflag=-1
c               write (6,*)'return 1'
               return 1
            end if
         elseif((idelim.ge.0).and.(jdelim.gt.0))then
            idelim=jdelim
c            write (6,*)'return 3'
            return 3
         end if
      end do
      if((idelim.ge.0))idelim=jdelim
c      write (6,*)'plain return'
      return
c **  ** alternate returns from xgtarg
c **  ** fell off the end
910   continue
c      write (6,*)'return 1 910'
      return 1
c **  ** infinite skip condition
920   continue
      nret = nreq
c      write (6,*)'return 2 920'
      return 2
      end
      subroutine xgtmch(string, iparse, matchs, nmatch, type, imatch,
     &                  iflag, idelim, *, *)
c		rashafer 86 mar 5
c	xparse subroutine to match the next argument on the parse string
c	with a list of allowed matches.  if the argument does not match
c	the argument is discarded, and the user is prompted for a
c	replacement.
c	string	c*	i/r: parse string
c	iparse	i4	i/r: parse position
c	matchs	c*(nmatch)	i: array of allowed matches
c	nmatch	i4	i: no. of allowed matches
c	type	c*	i: description of the entries in the arrray
c	imatch	i4	i/r: the index in the array of the next argument
c			if the argument is skipped over, then imatch is
c			unchanged.
c			(if <= 0 no successful match, see iflag)
c	iflag	i4	r: unusual circumstance flag
c			1 - fell off the end of the argument string
c			-1 - eof while prompting for replacement string
c	idelim	i4	r: the delimeter flag for what follows the argument
c	alternate returns:
c		*1 - fell off the string
c		*2 - eof while prompting for replacement string
	character*(*) string,matchs(*),type
        character*80 strtmp1,strtmp2
	integer iparse,nmatch,imatch,iflag,idelim
c
	integer lenact
	integer ibeg,iend,iflag2,ierr, jmatch
	logical qskip
c **	** come from for a successful insertion after a bad match
100	continue
	iflag=0
	call xgtarg(string,iparse,ibeg,iend,qskip,iflag2,idelim)
	if((qskip).and.(imatch.gt.0))then
	    return
	    end if
	if(iflag2.ne.0)then
	    iflag=1
	    return 1
	    end if
	jmatch=imatch
	if(.not.qskip)call xmatch(string(ibeg:iend),matchs,nmatch,jmatch)
	if(jmatch.le.0)then
	    call xunids(string(ibeg:iend),matchs,nmatch,jmatch,type)
	    if(jmatch.lt.0) imatch = min(nmatch,abs(jmatch))
	    if(imatch.gt.0) then
                strtmp1='insert selection: ('//matchs(imatch)
     &		 (:lenact(matchs(imatch)))//') '
                 strtmp2=string
		call xinfix(strtmp1,strtmp2,iparse,ierr)
	    else
		call xinfix('insert selection (no default):',string,
     &		 iparse,ierr)
		end if
	    if(ierr.lt.0)then
		imatch=jmatch
		iflag=-1
		return 2
		end if
c **	    ** now go back and try again
	    goto 100
	    end if
	imatch=jmatch
	return
	end
	subroutine xnxtrd(prompt, string, iparse, ierr, *, *)
c		xnxtrd - rashafer 27 feb 1986
c	xparse subroutine to flush out the current parse string and
c	read in a new one
c	prompt	c*	i: prompt string
c	string	c*	i/r: parse string
c	iparse	i4	i/r: parse position (on return zero, indicating
c			a fresh position)
c	ierr	i4	r: error flag. 0 - success
c				<0	- eof
c				>0	- i/o error of some sort on read
c	alternate returns:
c		*1 - eof
c		*2 - i/o error
	character*(*) prompt,string
	integer iparse,ierr
c
c	include 'xparinc.inc'
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
        external xprsbd
c
c **	** come from when xcheck detected a special input line and thus
c **	** the ordinary command line had to be read in
c        write (6,*)'in xnxtrd'
100	continue
	return_inquiry = .false.  
c  turn of inquiry mode for remainder
	call xpflsh(string,iparse,ierr)
	if(ierr.ne.0) return 1
	call xcread(prompt,string,ierr)
	return_inquiry = require_inquiry 
c  re-instate inquiry mode, if selected
       if(ierr.eq.0)then
          iparse=0
          call xcheck(string,iparse,ierr)
          if(ierr.ne.0) goto 100 
c     try another read
c       write (6,*)'xnxtrd ierr=',i4
       elseif(ierr.gt.1)then
          return 2
       end if
       if(ierr.lt.0)then
          return 1
       end if
       return
       end

      subroutine xcread(cprom, cbuf, ier)
      character cprom*(*), cbuf*(*)
      integer   ier
c---
c xparse subroutine to read in a command line.
c---
c cprom   i    the prompt string, if blank then there is no prompt.
c cbuf      o  the returned string after input (blank only if there
c              is an eof or an error, or if an empty string input)
c ier       o  =0 then something was read otherwise the io error flag
c              (if <0 then an eof was raised)
c---
c 1989-aug-08 - new gtbuf version - [aft]
c---
      integer   lbuf
c---
c      write (6,*)'in xcread'
      call gtbuf(cprom, ier)
      if(ier.ne.0) return
      call gtrest(cbuf, lbuf)
      return
      end
c*********
      block data xprsbd
c---
c xparse block data routine to handle the initial values
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
c
      data qxcase,qxpart,qxpfir/3*.true./
      data xpreof,lxpeof/'/*  ',2/
      data blank,opnstr,clsstr,comma,commnt,skpall,spcbr1,spcbr2,indrct
     &          /' ','"','"',',','!','/','|','|','@'/
            data contin,tab/'-','	'/
      data command,opsys,inquiry/'#','$','?'/
      data request_inquiry,require_inquiry,return_inquiry
     &       /'??',2*.false./
      data lgunit,trmcht,logcht/0,2*10/
      end

	subroutine xcheck(string, iparse, ierr)
c
c		rashafer 29 april 1987
c	xparse subroutine to check a line for any of the special characters
c	and then act on them if they are:
c		command (default #) pass on to the xparse routine
c		opsys (default $) pass on to the xopsys routine
	implicit none
c
	character*(*) string	
c  i/r: the string to be parsed
	integer iparse	
c  i/r: the parse position (moved to point to the
				
c  position before the next non-blank char if not
				
c  one of the special charcters: 
	integer ierr		
c  r: status flag: 0 - no special char found
				
c 		1 = a special char found and
				
c  		successfully handled (should
				
c 		probably initiate a re-read
				
c 		-1 an eof was generated somewhere
				
c 		along the way
				
c 		other:  some error condition
c
c	include 'xparinc.inc'
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
        external xprsbd
c
	integer lenstr
	logical xchkbl
c
	lenstr = len(string)
	do while(xchkbl(string(iparse+1:iparse+1)).and.(iparse.lt.lenstr))
	    iparse=iparse+1
	    end do
	ierr=0
	if(iparse.lt.lenstr)then
	    if(string(iparse+1:iparse+1) .eq. command) then
		iparse=iparse+1
		call xcomnd(string,iparse,ierr)
		if(ierr.eq.1)then
		    ierr=2
		elseif(ierr.eq.0)then
		    ierr=1
		    end if
	    elseif(string(iparse+1:iparse+1) .eq. opsys) then
		iparse=iparse+1
		call xopsys(string,iparse)
		ierr=1
		end if
	    end if
	return
	end
	subroutine xopsys(instrg, lenn)
	character*(*) instrg
	integer lenn
c---
c spawn a process to execute operating system commands
c---
c instrg  i/r  parse string
c lenn	  i/r  current parse position
c---
c rashafer 4 dec 1984
c---
	integer lcopy
	integer lenact, iend, ibeg, iflag, iret, idelim
c
c	include 'xparinc.inc'
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
        external xprsbd
	character*100 copy
	logical qskip, request_inquire

	if(lenn.lt.len(instrg))then
	    copy = instrg(lenn+1:)
	    if(lenact(copy).gt.0)then
		request_inquire = .false.
c  turn off inquiry mode for blank lines
		call xgtarg(instrg,lenn,ibeg,iend,qskip,iflag,idelim,*200,*200
     &		 ,*200)
		if((.not.qskip).and.instrg(ibeg:iend).eq.inquiry)then
		    call xinfix(' input an empty line or a single command:'
     &		     ,instrg,lenn,iflag)
		    if(iflag.ne.0)then
			lenn=len(instrg)
			return
			end if
		    copy=instrg(lenn+1:)
		    if(lenact(copy).le.0)goto 200
		    end if
		lcopy = lenact(copy)
		write(*,*) 'spawning...'
		call spawn(copy, lcopy, iret)
		lenn=len(instrg)
		return
		end if
	    end if
c **	** come from for an empty line on the parse string
200	continue
	write(*,*)' type ''lo'' to return to program'
	call spawn(' ',0,iret)
	lenn=len(instrg)
	return
	end

	subroutine xinfix(prompt,string,iparse,ierr,*,*)
c		rashafer 1986 mar 5
c	xparse subroutine to prompt the user at the terminal for a fixed
c	parameter to be inserted before the current parse position
c	prompt	c*	i: prompt string
c	string	c*	i/r: parse string
c	iparse	i4	i/r: parse position
c	ierr	i4	r: error flag returned during the read (n.b. that
c			if the insertion fails no error is currently generated
c			although it is signaled by the secondary return)
c	alternate returns:
c		*1 - eof
c		*2 - the insertion failed for space (or someother i/o error)
	character*(*) prompt,string
	integer iparse,ierr,ierr2
c	include 'xparinc.inc'
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
        integer lenact
	integer oldlen,newlen
	call xnewcm('*',.true.,ierr)
	if(iparse.lt.len(string))then
c	    ** xprmsw is a workspace area
	    xprmsw=' '//string(iparse+1:)
	    oldlen=lenact(xprmsw)
	else
	    oldlen=0
	    end if
	call xcread(prompt,string,ierr)
	if(ierr.ne.0)then
	    if(oldlen.gt.1)then
	        string=xprmsw(2:oldlen)
	    else
		string=' '
		end if
	    iparse=0
	    if(ierr.lt.0)then
c		** n.b. that there was no need to remove the terminal from
c		** the command file when the eof was generated (ierr<0) as it
c		** was already done by xcread
		return 1
	    else
		call xrmvcm(ierr2)
		return 2
		end if
	    end if
	call xrmvcm(ierr2)
	newlen=lenact(string)
	if(newlen.le.0)then
	    newlen=1
	    string(1:1)=comma
	    end if
	iparse=0
	string(newlen+1:)=xprmsw(:oldlen)
	if(newlen+oldlen.gt.len(string))then
	    write(*,*)'*warning*: during an insertion some data was lost:',
     &		'"',xprmsw(len(string)-newlen+1:oldlen),'"'
	    return 2
	    end if
	return
c **	** handler for bad return on file open... equivalent to an eof
900	continue
	ierr=-1
	return 1
	end

      subroutine xgtarg(string, iparse, ibeg, iend, qskip, iflag,
     :   idelim, *, *, *)
      character string*(*)
      integer   iparse, ibeg, iend, iflag, idelim
      logical   qskip
c---
c version 2.0:  allows for the request_prompt string (default ??
c      to be checked for.  also, if the stirng is in prompt mode
c      then instead of the falling off the end condition being
c      returned, a ? is returned (inquiry character)
c---
c string  i    the string to be parsed
c iparse  i/r  the parse position, at the call the value
c              should be at the position just before where parsing
c              is to begin.  therefore, for a string just begun it
c              should be 0.  on return it will be set at the position
c              appropriate to obtain the next argument.  when
c              iparse is > len(string), then there was no
c              argument found before the end of the string reached
c              (see the iflag = 1 condition below) while if
c              iparse = len(string) the next call to xgtarg
c              will run off the end.  (not all run-off-the-ends
c              can be predicted this way when continuation chars
c              are used).
c              n.b. as
c              each argument is parsed, its form may be modified,
c              particularly arguments in character strings, so it
c              is difficult to backup over previously parsed strings.
c              thus, except for new strings (where iparse must
c              be initially set to zero) the value should not
c              be adjusted by the calling program.
c ibeg      r  the first character of the returned argument
c iend      r  the last character of the returned argument
c qskip     r  if true, the argument field has been 'skipped over'
c              (i.e., it contains an empty field).
c iflag     r  the condition flag:
c              -1 - an eof was generated while processing a
c                   continuation.  iparse is set > len(string).
c               0 - nothing unusual
c               1 - the parse fell off the end (also indicated
c                   by iparse > len(string)).
c               2 - the field reached the terminal skip character.
c                   this and subsequent calls will raise the
c                   qskip=.true. condition.
c idelim    r  the delimeter flag:
c               0 - the argument was normally delimeted by a coma,
c                   end of line, or the infinite skip char.
c               1 - the argument was delimeted by the
c                   first special delimeter
c               2 - the argument was delimeted by the
c                   second special delimeter
c alternate returns:
c      alternate returns are based on the value of abs(iflag)
c---
c 20 october 1985 - rashafer 
c---
      external xprsbd
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
      logical xchkdl,xchkbl,qdel,qdone,qcopy,qquote
      integer lenstr,icopy,ival,jval
c---
      lenstr=len(string)
      iflag=0
      idelim=0
      iparse=iparse+1
      ibeg=1
      iend=0
c** come from for case of continuation char processing or the
c** request_inquiry string has tobggled off return_inquiry
  100 continue
c** skip to first non-blank character
      do while((iparse.le.lenstr).and.xchkbl(string(iparse:iparse)))
          iparse=iparse+1
      end do
      if(iparse.gt.lenstr) then
         qskip=.false.
         if(return_inquiry) then
            iparse=lenstr-1
            ibeg=lenstr-1
            iend=lenstr-1
            string(ibeg:lenstr)=inquiry
            return
         else
            iparse=lenstr+1
            iflag=1
            return 1
         end if
      end if
      qquote=.false.
      if(string(iparse:iparse).eq.opnstr) then
         qquote=.true.
         qskip=.false.
         ibeg=iparse+1
         qdone=.false.
         qcopy=.false.
         do while((iparse.lt.lenstr).and.(.not.qdone))
            iparse=iparse+1
            if(string(iparse:iparse).eq.clsstr) then
               if(iparse.eq.lenstr) then
                  qdone=.true.
               else
                  qdone=string(iparse+1:iparse+1).ne.clsstr
                  if(.not.qdone) then
                     if(.not.qcopy) then
                        qcopy=.true.
                        icopy=iparse
                     end if
                     iparse=iparse+1
                  end if
               end if
            end if
            if(qcopy) then
               string(icopy:icopy)=string(iparse:iparse)
               icopy=icopy+1
            end if
         end do
         if(qdone) then
            if(qcopy) then
               iend=icopy-2
            else
               iend=max(iparse-1,ibeg)
            end if
            iparse=iparse+1
         else
            if(qcopy) then
               iend=icopy-1
            else
               iend=lenstr
            end if
         end if
c** as a quoted string ends on the end-of-string character
c** then a false qdel means that the end was not found by
c** a delimeter (such as comma or the specials), so the
c** existence of the specials must still be checked
         qdel=.false.
      elseif(xchkdl(string(iparse:),ival)) then
c** the first argument of the string was a delimeter, indicating
c** that either there was an empty field (i.e. just to a comment)
c** or a skipped field to a comma, infinite-skip, or some special
c** delimeter
         if(ival.eq.1) then
            if(return_inquiry) then
               iparse=lenstr-1
               ibeg=lenstr-1
               iend=lenstr-1
               string(ibeg:lenstr)=inquiry
               return
            else
               iparse=lenstr+1
               iflag=1
               return 1
            end if
         else
            qskip=.true.
            if(ival.gt.0) then
               if(ival.eq.2) then
                  iparse=iparse-1
                  iflag=2
                  return 2
               else
c** a special delimeter
                  idelim=ival-2
               end if
            end if
            return
         end if
      else
c** some ordinary character
         qskip=.false.
         ibeg=iparse
         qdel=.false.
         qdone=.false.
         do while((.not.qdone).and.(iparse.lt.lenstr))
            iparse=iparse+1
            if(xchkdl(string(iparse:),ival)) then
               qdone=.true.
               qdel=.true.
            else
               qdone=xchkbl(string(iparse:iparse))
            end if
         end do
         if(qdone) then
            iend=iparse-1
         else
c** iparse must be the end of the string
            iend=iparse
c** set iparse to point past the end of the string (no more
c** processing will then probably be necessary)
            iparse=iparse+1
         end if
      end if
c** if not already at a delimeter, then advance until a non-blank
c** character is found
      if(.not.qdel) then
         do while((iparse.le.lenstr).and.xchkbl(string(iparse:iparse)))
            iparse=iparse+1
         end do
         if(iparse.le.lenstr) then
            if(xchkdl(string(iparse:),jval)) then
               if(jval.eq.1) then
                  iparse=lenstr
               elseif(jval.eq.2) then
                  iparse=iparse-1
               elseif(jval.ge.3) then
c** as a special delimeter ends the argument, qdel is
c** set to allow for the appropriate special return
                  qdel=.true.
                  ival=jval
               end if
            else
               iparse=iparse-1
            end if
         else
c** reset iparse to point to the last char in the string, so
c** that iparse > lenstr is reserved for the condition that
c** this call to xgtarg caused a run-off-the-end, while
c** iparse = lenstr means that the next call to xgtarg will
c** generate that condition
            iparse = lenstr
         end if
      end if
      if(qdel) then
c** the argument ended on a delimeter
         if(ival.eq.1) then
            iparse=lenstr
         elseif(ival.eq.2) then
            iparse=iparse-1
c** here the infinite-skip condition char was found.  the
c** special condition will not be raised until the next
c** argument is parsed.
            return
         elseif(ival.ge.2) then
            idelim=ival-2
            return
         end if
      end if
c** now check for continuation condition
      if((iparse.eq.lenstr).and.(.not.qquote).and.(string(ibeg:iend).eq.
     &    contin)) then
         call xcread('->',string,ival)
         if(ival.eq.0) then
            iparse=1
            goto 100
         end if
         iflag=-1
         return 1
      end if
c** check for inquiry request
      if((.not.qquote).and.(string(ibeg:iend).eq.request_inquiry)) then
         return_inquiry=.not. return_inquiry
         if(.not. return_inquiry) then
            iparse=iparse+1
            goto 100
         end if
         ibeg=iend
         string(ibeg:iend)=inquiry
      end if
      return
      end

	subroutine xcomnd(string,iparse,ierr)
c		rashafer 30 april 1987
c		subroutine to look for special xparse commands
c
	character*(*) string	
c  i/r: parse string
	integer iparse	
c  i/r: parse position
	integer ierr		
c  r: condition flag 0 = success
				
c 	-1 = eof during process
				
c 	positive values:  error during command
c
c	include 'xparinc.inc'
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
        external xprsbd
	integer ncom
	parameter (ncom=12)
	character*10 xcom(ncom)
	integer icom,idelim,jerr
	integer lenact
	logical xgtonf
	data xcom/'chatter','debug',
     &	 'exit','fix*','help*','inquire','log','pause','quit',
     &	 'set*','suggest*','xparsehelp'/
	data icom/0/

	call xgtmch(string,iparse,xcom,ncom,
     &	 'xparse supported commands (# commands)',icom,ierr,idelim)
	if(ierr.lt.0)return
	ierr=0
	if(icom.le.0) return
c **	** do case for command index
	goto (650,500,200,100,100,300,600,400,200,100,100,100),icom
100	continue
        call xwrite(' xparse command "'//xcom(icom)(:lenact(xcom(icom)))
     &	            //'" is not currently implemented.',1)
	goto 9000
c **	** exit or quit handler
200	continue
	call xwrite(' xparse generated exit',11)
	call exit
c **	** inquire handler
300	continue
	require_inquiry=xgtonf(string,iparse,'set inquiry mode (??)',
     &	     require_inquiry,jerr)
	goto 9000
c **	** pause handler
400	continue
	call xwrite(' program paused, use continue to resume',-1)
	pause 'program paused, use continue to resume'
	goto 9000
c **	** debug handler
500	continue
	call xdebug
	goto 9000
c **	** log handler
600	continue
	call xlogit(string,iparse,ierr)
	goto 9000
c **	** chatter handler
650	continue
	call xchatr(string,iparse,ierr)
	goto 9000
c **	** end case for command index
9000	continue
	return
	end
	subroutine xchatr(string,lenn,iflag)
c	setcht		rashafer	16 march 1985
c		xspec subroutine to set the chattyness level
c
c	version 2.0 (xchatr)  xparse version
	implicit none
	character*(*)string 
c 	i/r: parse string
	integer lenn		
c  i/r: parse position
	integer iflag 
c  r: the status condition (-1 is an eof)
c
	character*14 descr(2)
	integer chatvl(2)
	integer nret
	data descr/'chattyness','log file chat.'/
	call xgtcht(chatvl(1),chatvl(2))
	call xgtint(string,lenn,2,descr,2,chatvl,nret,iflag,-1)
	call xchaty(chatvl(1),chatvl(2))
	return
	end
	subroutine xgtcht (inter, log)
c		rashafer 22 march 1986
c	xparse subroutine to return the current interactive and log file
c	chattyness levels
	implicit none
c	inter	i4	r: the current interactive chattyness
c	log	i4	r: the current log file chattyness
	integer inter,log
c	include 'xparinc.inc'
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
        inter=trmcht
	log=logcht
	return
	end

	subroutine xgtint(string, iparse, nreq, descr, ndesc, retval,
     &	                  nret, iflag, idelim, *, *, *)
c		rashafer 14 april 1986
c	xparse subroutine to return a list of i*4 (using xgtnum)
	character*(*) string,descr(*)
	integer retval(*),iparse,nreq,ndesc,nret,iflag,idelim
	call xgtnum(string,iparse,nreq,descr,ndesc,retval,retval,0,3,
     &		retval, nret,
     &		iflag, idelim, *910, *920, *930)
	return
910	continue
	return 1
920	continue
	return 2
930	continue
	return 3
	end

	subroutine xchaty(inter,log)
c		rashafer 22 mar 1986
c	xparse subroutine to set the chattyness levels used by xwrite (et al.)
	implicit none
c	inter	i4	i: the chattyness level for the interactive terminal
c	log	i4	i: the chattyness level for the log file writes
c		n.b. if inter or log are < 0 then the current values are not
c		modified
	integer inter,log
c	include 'xparinc.inc'
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
        if(inter.ge.0)trmcht=inter
	if(log.ge.0)logcht=log
	return
	end

	subroutine xdebug
c---
c xparse subroutine to handle the #debug condition.
c---
c may 14 1987 - rashafer
c---

	call debugh

c for vax case runs code :
c	external ss$_debug
c	call lib$signal(ss$_debug)

	return
	end

      subroutine xlogit(cstr, iparse, iflag)
      character cstr*(*)
      integer   iparse, iflag
c---
c xspec subroutine to open a log file... if the file name
c is none then the file is closed.  also sets the logging
c of command file info.
c { [ <filename> <base log> <step size> ] |
c        [   none   ] |
c        [   set     [ [ default <default name> ] |
c              [ append [ on | off ] ]    |
c              [ stamp  [ on | off ] ]    |
c              [ id <id string> ] ] ] }
c---
c cstr    i    parse string
c iparse  i/o  parse position string
c iflag     o  status flag 0=success, -1=eof during parse
c---
c 1990-sep-14 - modified to support log file driver [aft]
c version 2.0 (xlogit)  new xparse version for log file opening syntax
c setlog   rashafer 7 aug 1984
c---
c     include  'xparinc.inc'
      integer   nset_opt
      parameter (nset_opt=4)
      integer   lenact
c
      character*60 filename
      character*60 keyword
      character*60 id
      character*40 default_file
      character*35 descr(2)
      character    ctmp
      character*7 set_opt(nset_opt)
      real      rbuf
      integer ios, lstr, nbuf
      integer nret,inqunit,ierr,idelim
      integer chatvl(2)
      integer iset_opt
      logical exist, opnlog
      logical qpart,xqmtch,xgtonf,qjunk
      logical append, stamp
      data default_file/'log.log'/
      data filename/'log.log'/
      data set_opt/'append','default','id','stamp'/
      data descr/'chattyness level to log commands',
     &      'increment for indirect command files'/

      keyword=filename
      call xgtstr( cstr, iparse, 1,
     &    'log file name, `none'', or `set''', 1, keyword, nret,
     &    iflag, -1)
      call logger(2, rbuf, nbuf, ctmp, lstr)
      opnlog = rbuf.ne.0.
      if(nret.le.0) then
c** if the file is already open, or if an eof occured during the
c** handling of the '?', return, otherwise use the current filename
         if((opnlog).or.(iflag.lt.0)) return
      end if
      if(xqmtch('none',keyword,qpart)) then
         if(opnlog) call xclslg('keep')
      else if(xqmtch('set',keyword,qpart)) then
         call xgtmch( cstr, iparse, set_opt, nset_opt,'set option',
     &        iset_opt,
     &        iflag, idelim)
         if(iflag.ne.0) return
         goto (110,120,130,140) iset_opt
c** do case iset_opt
c**  append handler
  110    continue
         append = xgtonf(cstr,iparse,'append log file status',append,
     &       iflag)
         goto 199
c**  default handler
  120    continue
         call xgtstr(cstr,iparse,1,'default log file name',1,
     &       default_file,nret,iflag,idelim)
         goto 199
c**  id handler
  130    continue
         call xgtstr(cstr,iparse,1,'log file id string',1,id,nret,
     &       iflag,idelim)
         goto 199
c**  stamp handler
  140    continue
         qjunk = xgtonf(cstr, iparse, 'log time stamp status', stamp,
     &       iflag)
         goto 199
c** end case iset_opt
  199    continue
      else
         inquire(file=keyword,name=filename,exist=exist,number=inqunit
     &        ,iostat=ios)
         if(.not.(opnlog.and.exist.and.(inqunit.eq.nint(rbuf)))) then
c** the file is not opened currently to the input file
            if(opnlog) call xclslg('keep')
            inqunit=rbuf
            call xopnlg(filename, append, id, stamp,
     &       inqunit, iflag)
            if(iflag.ne.0) then
               call xwrite(' unable to open the log file `'
     &           //filename(:lenact(filename))//'''',5)
               return
            end if
         end if
c** get the command file chattyness levels
         call xgtint(cstr,iparse,2,descr,2,chatvl,nret,iflag,-1)
      end if
      return
      end

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

      subroutine xopnlg(cfile, qpend, cstr, qtime, iunit, ierr)
      character cfile*(*), cstr*(*)
      logical   qpend, qtime
      integer   iunit, ierr
c---
c xparse subroutine to open the log file for output.
c---
c cfile     i/o  the file to be opened
c qpend     i    if true, the file is oppend append, else new
c cstr      i    string used as an initial line(with form feed)
c qtime     i    if true, a date/time stamp is placed on the
c                -initial line
c iunit     i/o  if iunit>0 then it is the fortran unit number to be used.
c                -if iunit<=0 then on output it is the fortran unit
c                -number of the log file, as generated by getlun
c ierr        o  error flag from open
c---
c 1990-sep-12 - new version to use logger driver [aft]
c---
      integer   lenact
c
      character ctmp*80
      character ctime*18
      real      rbuf(2)
      integer   ltmp, nbuf
c---
      rbuf(1)=iunit
      if(qpend) then
         rbuf(2)=1.0
      else
         rbuf(2)=0.0
      end if
      nbuf=2
      ltmp=lenact(cfile)
      call logger(3, rbuf, nbuf, cfile, ltmp)
      if(nbuf.lt.0) return
      iunit=rbuf(1)
c---
c if requested write current date and time
      if(qtime) then
c         call getdat(ctime(1:9))
c         ctime(10:10)=' '
c         call gettim(ctime(11:18))
      else
          ctime=' '
      end if
c---
c write startup string as a comment.
      ltmp=lenact(cstr)
      if(ltmp.gt.0 .or. qtime) then
         write(ctmp,'(a,1x,a)') '!'//cstr(:max(1,ltmp)),
     :      ctime(:max(1,lenact(ctime)))
         call logger(5, rbuf, nbuf, ctmp, ltmp)
      end if
      return
      end

      subroutine xgtstr (string, iparse, nreq, desc, ndesc, retstr, 
     &                   nret, iflag, idelim, * , *, *)
c		rashafer 9 march 1986
c	xparse subroutine to peel off a requested number of arguements as
c	strings
c	string	c*	i/r: parse string
c	iparse	c*	i/r: parse position
c	nreq	i4	i: no. of strings requested
c	desc	c*(ndesc)	i: description of the string requested
c	ndesc	i4	i: no. of descriptions passed (if = 0 then '?' will
c			not trigger a description)
c	retstr	c*(nreq)	r: the strings actually picked up (no change
c				on skips)
c	nret	i4	r: no. of strings actually processed (where skips
c				and infinite skips are included as processed)
c				if nret not = nreq then there was a fall off
c				the end of the string, or a special delimeter
c				met (see input value of idelim)
c	iflag	i4	r: value of condition flag at last call to xgtarg
c	idelim	i4	i/r: if <= -1  there is no checking for special
c			delimeters, else on return idelim contains the value
c			of the delimeter that triggered the return (if 0
c			then a comma, if 1 then special del 1, etc.)
c	alternate returns:
c 		*1	same as first alternate return of xgtarg (fell off
c				line)
c		*2	same as second alternate return of xgtarg (infinite
c				skip)
c		*3	special delimeter (see idelim)
      character*(*) string,desc(*),retstr(*)
      integer iparse,nreq,ndesc,nret,iflag,idelim
c
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
c
      integer lenact
      integer jdelim,ibeg,iend,destmp,ierr
      logical qskip
      nret=0
      do while(nret.lt.nreq)
         call xgtarg(string,iparse,ibeg,iend,qskip,iflag,jdelim,*910,
     &               *920)
         nret=nret+1
         if(.not.qskip)then
            if((ndesc.gt.0).and.(string(ibeg:iend).eq.inquiry))then
               destmp=min(ndesc,nret)
               xprmsg=desc(destmp)(:lenact(desc(destmp)))//' : ('//
     &                retstr(nret)(:max(1,lenact(retstr(nret))))//')'
               call xinfix(xprmsg(:lenact(xprmsg)+1),string,iparse,ierr)
               nret=nret-1
               if(ierr.eq.1)then
                  iflag=-1
                  return 1
               end if
            else
               retstr(nret)=string(ibeg:iend)
               if((idelim.ge.0).and.(jdelim.ne.0))then
                  idelim=jdelim
                  return 3
               end if
            end if
         end if
      end do
      if((idelim.ge.0))idelim=jdelim
      return
c **  ** alternate returns from xgtarg
c **  ** fell off the end
910   continue
      return 1
c **  ** infinite skip condition
920   continue
      nret = nreq
      return 2
      end

	logical function xgtonf(string, iparse, type, default, ierr)
c		rashafer 30 april 1987
c	xparse subroutine to look for the string on or off on the parse string
	implicit none
c
c	xgtonf	returns true if the next string is "on"
c		returns false if the next string is "off"
	character*(*) string
	integer iparse
	character*(*) type	
c  i: a descriptor string
	logical default	
c  i: the value returned if skipped over or
				
c 	run off the end
	integer ierr		
c  r: status -- 0 = success
				
c 		-1 = eof
c
	character*3 opt(2)
	integer iopt, idelim
	data opt/'on','off'/
	if(default)then
	    iopt=1
	else
	    iopt=2
	    end if
	call xgtmch(string,iparse,opt,2,type,iopt,ierr,idelim)
	xgtonf = iopt.eq.1
	return
	end

      logical function xqmtch(test, base, qpart)
      character*(*) test,base
      logical   qpart
c---
c xparse function to see if a string is the same as a base string
c---
c xqmtch       function result: if true, then a match exists
c test    i    the string to be tested (this must have already been
c              -lenact-ed so that its length has removed padding
c base    i    the string to be compared against
c qpart     o  if true, then only a partial match was made (only
c              -allowed if the flag qxpart is true), otherwise
c              -a complete match.
c---
c 1985-oct-18 - rashafer
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
      integer   lenact, ixplwr
      integer   ltest, lbase, itest, ibase, i
c---
      ltest=len(test)
      lbase=lenact(base)
      qpart=.false.
      if(lbase.lt.ltest) then
c** the base string was too short
         xqmtch=.false.
         return
      elseif(lbase.gt.ltest) then
c** the base string was longer than the test string
         if(.not.qxpart) then
c** partial matches are not allowed
            xqmtch=.false.
            return
         else
            qpart=.true.
         end if
      end if
c---
      if(qxcase)then
c** case is not significant so do the conversion
         do i = 1,ltest
            itest=ichar(test(i:i))
            ibase=ichar(base(i:i))
            if(itest.ne.ibase) then
               itest=ixplwr(itest)
               ibase=ixplwr(ibase)
               if(itest.ne.ibase) then
                  xqmtch=.false.
                  return
               end if
            end if
         end do
         xqmtch=.true.
      else
         xqmtch=test.eq.base(:ltest)
      end if
      return
      end

      logical function xchkbl(string)
      character*(*) string
c---
c function to check if a character is one of the xparse blank
c characters (blank or horizontal tab).
c---
c xchkbl    o  if true, then char is a delimeter
c string  i    the first character is the one to be checked
c---
c 1985-oct-22 - rashafer
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
c---
      xchkbl=(string(1:1).eq.blank).or.(string(1:1).eq.tab)
      return
      end

      logical function xchkdl(string, ival)
      character*(*) string
      integer   ival
c---
c function to check if a character is one of the xparse delimeters.
c ival
c 0    = comma
c 1    = commnt - the comment begin char
c 2    = skpall - the terminal skip character
c 3    = spcbr1 - the first 'special break' char
c 4    = spcbr2 - the second 'special break' char
c---
c xchkdl    o  if true, then char is a delimeter
c string  i    the first character is the one to be checked
c ival      o  if xchkdl is true, this flag indicates which
c              delimeter was matched.
c---
c 1985-oct-22 - rashafer
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
      character*1 char
      integer   ipoint
c---
      ival=-1
      char=string(1:1)
      if(char.eq.comma) then
         ival=0
      elseif(char.eq.commnt) then
         ival=1
      elseif(char.eq.skpall) then
c** check if there are no other non-blank characters before the end
c** of line (or the comment char)
         ipoint=2
         do while((ipoint.le.len(string)).and.
     :            (string(ipoint:ipoint).eq.blank))
            ipoint=ipoint+1
         end do
         if((ipoint.gt.len(string)).or.
     :      (string(ipoint:ipoint).eq.commnt)) then
             ival=2
         end if
      elseif(char.eq.spcbr1) then
         ival=3
      elseif(char.eq.spcbr2) then
         ival=4
      end if
      xchkdl=ival.ge.0
      return
      end

	subroutine xpflsh (string, iparse, ierr, *)
c		xpflsh    rashafer feb 27 86
c	xparse subroutine to `flush' out the parse string (skipping over
c	any continuation cards developed)
c	string	c*	i/r: parse string
c	iparse	i4	i/r: parse posistion (on return equal to the
c			length of the string)
c	ierr	i4	r: error flag.  0 - success
c			-1 = eof during the reading of a continuation.
c	alternate returns:
c		*1	eof (ierr = -1)
	character*(*) string
	integer iparse,ierr
c
	integer jerr,ibeg,iend,iflag,idelim,lenstr
	logical qskip
c
	lenstr = len(string)
	do while(iparse.lt.lenstr)
		call xgtarg(string,iparse,ibeg,iend,qskip,iflag,idelim)
c		** the end has come if an infinite skip or an eof under
c		** a continuation processing is reached
		if(iflag.eq.2)iparse=lenstr
		end do
	if(iflag.eq.-1)then
		ierr=-1
	else
		ierr= 0
		return
		end if
	return
	end

	subroutine xunids(chose,choice,nch,jch,type)
c		rashafer 19 oct 1985
c	    xparse subroutine to print a message when a string was unidentified
c	chose	c*	i: the string tested (special case if '?')
c	choice	c*(nch)	i: the possible values
c	nch	i4	i: the no. of possible values
c	jch	i4	i: the error indicator: if 0 then ther had been
c			no (allowed) match found.  if < 0 then
c			multiple matches were indicated starting with
c			choice (-hcg(
c	type	c*	i: the type of value being looked for
	implicit none
	integer nch,jch
	character*(*) chose,choice(nch),type
	external xprsbd
c	include 'xparinc.inc'
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
        integer linew
	parameter (linew=72)
	logical qpart,xqmtch,qpart2
	integer lenact
	integer lc,ncol,icol,ic,kch,ich,lentyp
	qpart=jch.lt.0
	lentyp=lenact(type)
	if(qpart)then
	    write(*,*)' a unique example of ',type(:lentyp),
     &	     ' must be chosen:'
	    ich= -jch
	else
	    ich = 1
	    if(chose.ne.inquiry)then
c **		** does not match the enquire mode
		write(*,*)'''',chose(:lenact(chose)),''' does not match.'
		end if
	    write(*,*)' choose from the following ',type(:lentyp),':'
	    end if
	lc= len(choice(1))
	ncol=linew/(lc+1)
	icol=0
	ic=0
	xprmsw(:linew)= ' '
	do kch = ich,nch
	    if((.not.qpart).or.(xqmtch(chose,choice(kch),qpart2)))then
		icol=icol+1
		xprmsw(ic+1:ic+lc) = choice(kch)
		if(icol.eq.ncol)then
		    write(*,'(1x,a)')xprmsw(:linew)
		    xprmsw(:linew)= ' '
		    icol=0
		    ic=0
		else
		    ic=ic+lc+1
		    end if
		end if
	    end do
	if(icol.gt.0)then
	    write(*,'(1x,a)')xprmsw(:linew)
	    end if
	return
	end

	subroutine xmatch(string,array,narray,iret,*)
c		rashafer 17 oct 1985
c	xparse subroutine to test if a string matches one of the strings
c	contained in an array of strings
c	string	c*	i: test string
c	array	c*(narray) i: array of strins to be matched against
c	narray	i4	i: no of eelements in the array
c	iret	i4	r: if > 0 then a match was fournd with element iret
c			   if = 0 no match was found
c			   if < 0 then a non-unique match was found
c			   where -iret is the first such none unique match and
c			   the qxpfir flag was not set
c	alternate returns:
c		1 -	iret <= 0
	implicit none
	character*(*) string,array(*)
	integer iret,narray
c	include 'xparinc.inc'
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
        external xprsbd
	integer lenact
	logical qpart,xqmtch
	integer jret,lenst,iarray
	jret=0
	lenst=lenact(string)
	if((.not.qxpart).and.lenst.gt.len(array(1)))then
	    iret=0
	    return 1
	    end if
	do iarray = 1,narray
	    if(xqmtch(string(:lenst),array(iarray),qpart))then
		if(qxpfir.or.(.not.qpart))then
		    iret=iarray
		    return
		    end if
		if(jret.eq.0)then
		    jret=iarray
		    end if
		end if
	    end do
	if(jret.eq.0)then
c	    ** ther was never any match for the string
	    iret=0
	else
	    iret=-jret
	    end if
	return 1
	end

	subroutine xseti4(inval,outval)
c		rashafer 28 may
c	xparse subroutine to equate two i*4s
	implicit none
	integer inval,outval
	outval=inval
	return
	end
	subroutine xsetr4(inval,outval)
c		rashafer 28 may
c	xparse subroutine to equate two r*4s
	implicit none
	real inval,outval
	outval=inval
	return
	end

	logical function xtstr4(value,min,max,idescr,descr)
c
c		rashafer 28 may 1986
c	xparse subroutine to test if a value is in range
c	real version
	implicit none
	real value,min,max
	integer idescr
	character*(*) descr
	character*15 valstr(7)
c	include 'xparinc.inc'
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
        xtstr4 = (value.ge.min).and.(value.le.max)
	if(xtstr4)return
c	** now create the error message
	write(valstr,'(3e12.4)')value,min,max
	if(idescr.gt.0)then
	    call xpmsg3(descr,valstr)
	else
	    call xpmsg3('an integer',valstr)
	    end if
	if(value.lt.min)then
	    value=min
	else
	    value=max
	    end if
	return
	end
	logical function xtsti4(value,min,max,idescr,descr)
c
c		rashafer 28 may 1986
c	xparse subroutine to test if a value is in range
c	integer version
	implicit none
	integer value,min,max
	integer idescr
	character*(*) descr
	character*12 valstr(7)
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
        xtsti4 = (value.ge.min).and.(value.le.max)
	if(xtsti4)return
c	** now create the error message
	write(valstr,'(3i10)')value,min,max
	if(idescr.gt.0)then
	    call xpmsg3(descr,valstr)
	else
	    call xpmsg3('an integer',valstr)
	    end if
	if(value.lt.min)then
	    value=min
	else
	    value=max
	    end if
	return
	end

	logical function xparr4( string, idescr, descr, retval)
c		rashafer 14 april 1986
c	xparse routine to parse out an integer
	implicit none
c	xparr4	l4	function: if true, an integer was parsed and
c			returned in retval
c	string	c*	i: string to hold value
c	idescr	i4	i: if > 0 a descriptor is passed
c	descr	c*	i: a descriptor of the variable
c	retval	r4	i/r: the current value
	character*(*) string, descr
	integer idescr
	real  retval
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
        external xprsbd
	logical xisint,xisnum,isnum
	integer sqzlen
	character*15 tmpstr
	real*8 dval
	isnum=xisnum(string,dval,xisint,sqzlen)
	if(isnum)then
	    xparr4=.true.
	    retval=dval
	else
	    xparr4=.false.
	    write(tmpstr,'(1pe12.4)')retval
	    call xsquez(tmpstr,' ',0,sqzlen)
	    if(string.eq.inquiry)then
		    if(idescr.gt.0)then
			call xpmsg1(descr,tmpstr(:sqzlen))
		    else
			call xpmsg1('number',tmpstr(:sqzlen))
			end if
	    else
		    if(idescr.gt.0)then
			call xpmsg2('a number',string,descr,tmpstr(:sqzlen))
		    else
			call xpmsg2('a number',string,' ',tmpstr(:sqzlen))
			end if
		    end if
	    end if
	return
	end
	logical function xpari4( string, idescr, descr, retval)
c		rashafer 14 april 1986
c	xparse routine to parse out an integer
	implicit none
c	xpari4	l4	function: if true, an integer was parsed and
c			returned in retval
c	string	c*	i: string to hold value
c	idescr	i4	i: if > 0 a descriptor is passed
c	descr	c*	i: a descriptor of the variable
c	retval	i4	i/r: the current value
	character*(*) string, descr
	integer idescr, retval
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
        external xprsbd
	logical xisint,qint
	integer sqzlen
	character*12 tmpstr
	real*8 dval
	call xisnum(string,dval,xisint,retval)
	if(xisint)then
	    xpari4=.true.
	else
	    xpari4=.false.
	    write(tmpstr,'(i10)')retval
	    call xsquez(tmpstr,' ',0,sqzlen)
	    if(string.eq.inquiry)then
		    if(idescr.gt.0)then
			call xpmsg1(descr,tmpstr(:sqzlen))
		    else
			call xpmsg1('integer',tmpstr(:sqzlen))
			end if
	    else
		    if(idescr.gt.0)then
			call xpmsg2('an integer',string,descr,tmpstr(:sqzlen))
		    else
			call xpmsg2('an integer',string,' ',tmpstr(:sqzlen))
			end if
		    end if
	    end if
	return
	end
	subroutine xpmsg2 (type, string, descr, value)
c		rashafer 14 april 1986
c	xparse message subroutine, to place an information message in
c	the xprmsg string (in the common block)
	character*(*) descr, value, string, type
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
        integer lenact
	xprmsg='unable to parse '//type//' from "'//
     &         string(:lenact(string))//'"'
	if(lenact(descr).gt.0)then
		xprmsg(lenact(xprmsg)+2:)='for '//descr
		end if
	xprmsg(lenact(xprmsg)+1:)= ': ('//value(:lenact(value))//')'
	return
	end

	subroutine xpmsg1 (descr, value)
c		rashafer 14 april 1986
c	xparse message subroutine, to place an information message in
c	the xprmsg string (in the common block)
	character*(*) descr, value
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
	integer lenact
	xprmsg = descr(:lenact(descr)) // ': (' //
     &	 value(:lenact(value)) // ')'
	return
	end

      subroutine gtbuf(cprom, ier)
      character cprom*(*)
      integer   ier
c entry gtchar
      character ctoken*(*)
      integer   ltoken
c entry gtdble
      real*8 dnum
c entry gtint
      integer   inum
c entry gtlev
      integer   lev
c entry gtmode
      integer   inmode
c entry gtreal
      real      rnum
c entry ldbuf
      character cmd(*)*(*)
      integer   ncmd, ierld
c entry xgetin
      character cext1*(*)
c entry xsetin
      character cext2*(*)
c---
c when first called this routine tries to read the string on the
c command line.  if no string was given or on future calls, this
c routine will prompt the user for input.
c---
c special first characters (must be first character on line)
c !     ignore the entire line (comment line)
c $     spawn to vms
c @file open next level of indirect file and read buffer from there.
c @     force buffer to be read from terminal (even when reading a file).
c special character can appear anywhere in line.
c !     (not enclosed in ") means rest of line is a comment.
c---
c cprom     i    used to prompt the user.
c ier         o  =0 no error, =-1 if ^z was pressed.
c---
c requires:         for:
c alf.for           alf, fpnum, isnum, alfsks
c lenact.for        lenact
c logger.for        logger
c script.for        script
c xtend.for         xtend
c xwrcmd            xwrcmd
c .xanlib]sys.for   rdforn, spawn, prompt
c---
c 1990-feb-02 - added parameters to indirect files [aft]
c 1989-dec-13 - added command recall [aft]
c 1989-jul-06 - system independent indirect command files [aft]
c---
c implementation of indirect file parameters:  currently mxargs
c determines the maximum total number of parameters used at all
c indirect levels.  the total number of characters in all parameter
c strings is the len(cargs) currently set to 1024.  since some levels
c may use no parameters and others many, a two level pointer system
c is used.  istack(1,level) determines the offset into ipos when
c reading commands at depth of level.  istack(2,level) is the number
c of valid parameters at that level.  let ioff=istack(1,level) then
c ipos(1,ioff+ipar) gives the start position in cargs for parameter
c number ipar.  ipos(2,ioff+ipar) is the end position.
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
c
      integer  mxlev
      parameter (mxlev=60)
      integer   mxprev
      parameter (mxprev=20)
      integer   mxargs
      parameter (mxargs=32)
      real      fpnum
      integer   lenact, isnum
c
      character cprev(mxprev)*132, cargs*1024
      save      cprev,             cargs
      character cbuf*256
      save      cbuf
      character cnam*64, ctmp*64, ctok*32
      character cdef*10
      save      cdef
      real      rbuf(2)
      real      tmp
      integer   lunsav(mxlev), istack(2,0:mxlev), ipos(2,mxargs)
      save      lunsav,        istack,            ipos
      integer   i, ibas, ic, ie, idspawn, ifnb, iload
      integer   ios, ipar, iper, iq, is, itmp, itread
      integer   lnam, lprom, ltok, ltmp, luntry
      integer   icnt, iecho, ifirst, itop, kp, lbuf
      save      icnt, iecho, ifirst, itop, kp, lbuf
      integer   laspos, lbase, level, mode, nargs, nload
      save      laspos, lbase, level, mode, nargs, nload
c
      data cdef/'xco'/
      data idspawn/0/, iecho/1/, ifirst/1/, iload/0/, itop/0/
      data laspos/0/, lbase/0/, level/0/
      data mode/0/, nload/0/, nargs/0/
c---
c- first time through, search for foreign arguments.
      lprom=lenact(cprom)
c      write (6,*)'in gtbuf',ifirst
      if(ifirst.ne.0) then
c fill recall buffer with blanks, ignoring any commands the user
c has previously loaded.
         do 80 i=nload+1,mxprev
            cprev(i)=' '
   80    continue
         do 90 i=1, mxargs
            istack(1,i)=0
            istack(2,i)=0
   90    continue
         icnt=0
c         call trlog('spawn_disable', 13, cbuf, lbuf)
         if(lbuf.gt.0) idspawn=1
c         call trlog('gtbuf_edit', 10, cbuf, lbuf)
c         call edicom(cbuf, lbuf)
c
c         call rdforn(cbuf, lbuf)
         ier=0
         if(lbuf.gt.0) then
c- read something, go parse it.
            itop=itop+1
            icnt=icnt+1
            if(nload.gt.0) then
c cause a previously loaded command to be run first.
               nload=nload+1
               cprev(nload)=cbuf(:lbuf)
               cbuf=cprev(1)
               lbuf=lenact(cbuf)
            else
c otherwise, just load the command into the command stack.
               cprev(1)=cbuf(:lbuf)
            end if
            ios=0
            goto 150
         else
c- if no foreign argument and no prompt, then exit quietly.
            if(lprom.eq.0) goto 900
         end if
      end if
      itread=0
c---
  100 ios=0
      ier=0
      kp=0
      if(level.gt.lbase) then
         read(lunsav(level),'(a)',iostat=ios) cbuf
         lbuf=lenact(cbuf)
      else
cxxx this is now obsolete code
         if(iload.ne.0) then
            cbuf=cprom
            lbuf=lprom
            ier=1
         else
cxxx end obsolete section
            if(nload.eq.0 .or. level.lt.0) then
               if(icedit.ne.0) then
c use super-duper single character io
c                  call qline (cprom, lprom, cprev, icnt, itop,
c     :                  mxprev, cbuf, lbuf, ios)
               else
c use standard fortran io
                  call prompt(cprom,0)
                  read(5,'(a)',iostat=ios) cbuf
                  lbuf=lenact(cbuf)
               end if
               if(lbuf.gt.0) then
c only push non-blank commands onto the stack
                  icnt=icnt+1
                  itop = mod(itop, mxprev) + 1
                  cprev(itop) = cbuf
               end if
            else
c read from a 'loaded' command
               icnt=icnt+1
               itop=itop+1
               if(itop.gt.mxprev) itop=1
               cbuf=cprev(itop)
               lbuf=lenact(cbuf)
               if(itop.eq.nload) nload=0
            end if
         end if
      end if
c---
c- use /* to denote eof
  150 if(cbuf(1:2).eq.'/*') ios=-1
      if(ios.ne.0) then
         lbuf=0
         if(level.gt.0) then
            ier= 1
            close(unit=lunsav(level))
            call frelun(lunsav(level))
            nargs=nargs-istack(2,level)
            istack(2,level)=0
            level=level-1
            laspos=ipos(2,level)
            if(iload.eq.0) goto 100
         else
            ier=-1
            if(level.lt.0) then
c user has entered an eof and hence is unable to correct an error that
c occurred while reading an indirect file.  log this response and
c unwind the indirect file stack.
               call xwrcmd(level, cprom, lprom, '/*', 2)
               level=abs(level)
               itmp=level
               do 160 i=1,itmp
                  close(unit=lunsav(level))
                  call frelun(lunsav(level))
                  nargs=nargs-istack(2,level)
                  istack(2,level)=0
                  level=level-1
  160          continue
               laspos=ipos(2,level)
            end if
         end if
      else
c---
c if the first character is '!' then ignore the complete line.
         if(mode.eq.0 .and. cbuf(1:1).eq.'!') goto 100
c otherwise, strip off comments, this is done by scaning the input
c line looking for a ! not enclosed in ".
  170    iq=0
         iper=0
         ic=lbuf
         do 180 i=1,ic
            if(iq.eq.0) then
c- even number of " marks.
               if(cbuf(i:i).eq.'"') then
                  iq=1
               else if(cbuf(i:i).eq.'!') then
                  lbuf=i-1
                  lbuf=lenact(cbuf(:lbuf))
                  goto 190
               else if(cbuf(i:i).eq.'%') then
                  if(iper.eq.0) iper=i
               end if
            else
c- odd number of " marks.
               if(cbuf(i:i).eq.'"') iq=0
            end if
  180    continue
  190    continue
c---
c got an gtbuf command, deal with it
         if(iper.ne.0) then
            ltok=index(cbuf(iper+1:),'%')+1
            if(ltok.le.0) goto 240
            ctok=cbuf(iper:iper+ltok-1)
            kp=iper+ltok-1
            call upc(ctok)
            if(ctok(1:3).eq.'%ed') then
c %edit% for command line editing
               call edicom(cbuf(kp+1:), lbuf-kp)
               goto 100
            else if(ctok(1:2).eq.'%e') then
c %echo% to control the printing of command read from indirect files.
               call alf(cbuf, lbuf, kp, ctok, ltok)
               if(ltok.le.0) goto 100
               call upc(ctok)
               if(ctok(1:2).eq.'of') then
                  iecho=0
               else
                  iecho=1
               end if
               goto 100
            else if(ctok(1:2).eq.'%l') then
c %log%
               call alf(cbuf, lbuf, kp, ctok, ltok)
               ctmp=ctok(1:3)
               call upc(ctmp(1:3))
               if(ltok.ge.3 .and. ctmp(1:3).eq.'%of') then
                  call logger(4, rbuf, itmp, ctok, ltok)
               else
                  rbuf(1)=0.0
                  rbuf(2)=0.0
                  itmp=2
                  call logger(3, rbuf, itmp, ctok, ltok)
               end if
               goto 100
            else if(ctok(1:2).eq.'%r') then
c %recall%all
               call alf(cbuf, lbuf, kp, ctok, ltok)
               if(ltok.le.0) goto 100
               call upc(ctok)
               if(ctok(1:1).eq.'a') then 
                  if(icnt.le.mxprev) then
                     itmp=0
                     ibas=0
                  else
                     itmp=itop
                     ibas=icnt-mxprev
                  end if
                  do 220 i=1,mxprev
                     itmp=itmp+1
                     if(itmp.gt.mxprev) itmp=1
                     ibas=ibas+1
                     write(*,211) ibas,cprev(itmp)(:lenact(cprev(itmp)))
  211                format(1x,i5,': ',a)
                     if(itmp.eq.itop) goto 100
  220             continue
                  goto 100
               end if
            else if(ctok(1:2).eq.'%s') then
c %script%
               call alf(cbuf, lbuf, kp, ctok, ltok)
               ctmp=ctok(1:3)
               call upc(ctmp(1:3))
               if(ltok.ge.3 .and. ctmp(1:3).eq.'%of') then
                  call script(4, tmp, itmp, ctok, ltok)
               else
                  itmp=0
                  call script(3, tmp, itmp, ctok, ltok)
               end if
               goto 100
            else if(isnum(ctok(2:),ltok-2).ne.0) then
               ipar=fpnum(ctok(2:),ltok-2,itmp)
               if(ipar.gt.0 .and. ipar.le.istack(2,level)) then
                  itmp=istack(1,level)+ipar
                  is=ipos(1,itmp)
                  ie=ipos(2,itmp)
                  ctmp=cargs(is:ie)//cbuf(iper+ltok:)
               else
                  ctmp=cbuf(iper+ltok:)
               end if
               cbuf(iper:)=ctmp
               iper=index(cbuf(iper:),'%')
               lbuf=lenact(cbuf)
               if(iper.gt.0) goto 170
            end if
         end if
c---
c write (unloaded) info to log file, with comments stripped.  this is
c because the prompt is written as a comment.
         if(iload.eq.0) then
            call xwrcmd(level, cprom, lprom, cbuf, lbuf)
c if reading an indirect file, then echo response.
            if(level.gt.lbase .and. iecho.ne.0) then
               write(*,201) cprom(:lprom),cbuf(:lbuf)
  201          format(1x,a,' ',a)
            end if
         end if
c---
c- check for special characters in first column.
  240    if(mode.eq.0) then
            if(cbuf(1:1).eq.'$') then
c- $ to spawn
               if(idspawn.eq.0) then
                  lbuf=lbuf-1
                  write(*,*) 'spawning...'
                  call edicom('sav',3)
                  call spawn(cbuf(2:),lbuf,ier)
                  call edicom('res',3)
               else
                  write(*,*) 'spawning has been disabled.'
               end if
               goto 100
            else if(cbuf(1:1).eq.'@') then
               if(lbuf.eq.1) then
c- @ with no file name means read next line from terminal
                  level=-level
                  itread=1
                  goto 100
               end if
c- indirect command file.
               if(level.ge.mxlev) then
                  write(*,261)  mxlev
  261             format(1x,'gtbuf--too many indirect files.  ',
     :               'max depth is',i6,'.')
                  ier=1
                  goto 900
               end if
               kp=1
               call alf(cbuf, lbuf, kp, cnam, lnam)
               call getlun(luntry)
               if(lenact(cdef).gt.0) call xtend(cnam,cdef)
c search current directory
               call openwr(luntry,cnam,'old',' ',' ',0,1,ios)
  270          if(ios.eq.0) then
                  level=abs(level)+1
                  lunsav(level)=luntry
                  istack(1,level)=istack(1,level-1)+istack(2,level-1)
                  istack(2,level)=0
                  if(kp.lt.lbuf) then
  280                call alf(cbuf, lbuf, kp, ctmp, ltmp)
                     if(nargs.ge.mxargs) then
                        write(*,281) mxargs,ctmp(:ltmp)
  281 format(' gtbuf--maximum number of arguments is',i5,'.'/
     :       ' gtbuf--no room for "',a,'".')
                        goto 300
                     end if
                     if(laspos+ltmp.gt.len(cargs)) then
                        write(*,291) ctmp(:ltmp)
  291 format(' gtbuf--arguments too long.  no room for "',a,'".')
                        goto 300
                     end if
                     istack(2,level)=istack(2,level)+1
                     nargs=nargs+1
                     ipos(1,nargs)=laspos+1
                     ipos(2,nargs)=laspos+ltmp
                     cargs(laspos+1:laspos+ltmp)=ctmp
                     laspos=laspos+ltmp
                     if(kp.lt.lbuf) goto 280
                  end if
c don't start reading indirect file on first call.
  300             if(ifirst.eq.0 .or. lprom.gt.0) goto 100
                  lbuf=0
                  goto 900
               end if
c search user defined directory (if it exists)
c               call trlog('my_xcoms',8,ctmp,ltmp)
               if(ltmp.gt.0) then
                  ctmp(ltmp+1:)=cnam
                  call openwr(luntry,ctmp,'old',' ',' ',0,1,ios)
                  if(ios.eq.0) goto 270
               end if
c search system directory
               ctmp=cnam
               call ptend('xanadu','lib/xcoms',ctmp)
               call openwr(luntry,ctmp,'old',' ',' ',0,1,ios)
               if(ios.eq.0) goto 270
c nowhere to be found
               call frelun(luntry)
               write(*,311) cnam(:lenact(cnam)),ios
  311          format(' gtbuf--unable to open ',a,', ios=',i5)
               if(ifirst.ne.0 .and. lprom.eq.0) goto 900
               if(iload.eq.0 .or. level.gt.lbase) goto 100
               ier=1
            end if
         end if
c---
c- skip leading blanks.
         call alfsks(cbuf,lbuf,kp)
         ifnb=kp
      end if
c---
c- gtrest should return preceeding blanks (to make select work)
c- hence, we must reset kp=0.
  900 iload=0
      lbase=0
      ifirst=0
      kp=0
      if(itread.ne.0) level=abs(level)
      return
c---
      entry gtchar(ctoken, ltoken)
c---
c return the next token as a character string.
c---
      call alf(cbuf,lbuf,kp,ctoken,ltoken)
      return
c---
      entry gtdble(dnum, ier)
c---
c return the next token as a real*8 number.  note, the
c parser is not used to decode the number, hence things like no or
c embedded arithmetic operators are not allowed.
c---
      call alf(cbuf,lbuf,kp,ctok,ltok)
      if(isnum(ctok,ltok).eq.0) goto 920
      read(ctok(:ltok),911,err=920) dnum
  911 format(d20.0)
      ier=0
      return
  920 ier=1
      return
c---
      entry gtint(inum, ier)
c---
c return the next token as an integer number.
c---
      call alf(cbuf,lbuf,kp,ctok,ltok)
      ier=0
      if(ltok.gt.0) then
         tmp=fpnum(ctok,ltok,ier)
         if(ier.eq.0) then
            inum=nint(tmp)
         else
            ier=1
         end if
      end if
      return
c---
      entry gtlev(lev)
c---
c return the current level of indirectness.  this assumes that the
c gtbuf parser has already been called.  in order to allow plt to
c be called by routines that do not use the gtbuf parser, set the
c 'already been called' flag now.
c---
      lev=level
      ifirst=0
      return
c---
      entry gtmode(inmode)
c---
c toggle the mode to inmode.  note, original mode is returned in
c inmode.
c---
      itmp=mode
      mode=inmode
      inmode=itmp
      return
c---
      entry gtprom(lev)
c---
c for the next call to gtbuf, stuff the prompt into cbuf.  this
c provides a tricky way to use gtbuf to parse lines generated by
c an external routine.
c---
      iload=1
      lbase=lev
      return
c---
      entry gtreal(rnum, ier)
c---
c return the next token as a real number.
c---
      call alf(cbuf,lbuf,kp,ctok,ltok)
      ier=0
      if(ltok.gt.0) then
         tmp=fpnum(ctok,ltok,ier)
         if(ier.eq.0) then
            rnum=tmp
         else
            ier=1
         end if
      end if
      return
c---
      entry gtpeek(ctoken, ltoken)
c---
c peek at the next token, i.e., return the token but do not advance
c the character pointer.
c---
      call alfsks(cbuf,lbuf,kp)
      if(kp.ge.lbuf) then
         ctoken=' '
         ltoken=0
      else
         ctoken(1:1)=cbuf(kp+1:kp+1)
         ltoken=1
         if(cbuf(kp+2:kp+2).ne.' ' .and. cbuf(kp+2:kp+2).ne.',') then
            ctoken(2:2)=cbuf(kp+2:kp+2)
            ltoken=2
         end if
      end if
      return
c---
      entry gtrest(ctoken, ltoken)
c---
c get rest of string, excluding comments.
c---
      if(kp.lt.lbuf) then
         ctoken=cbuf(kp+1:lbuf)
         ltoken=max(lbuf-kp,0)
      else
         ctoken=' '
         ltoken=0
      end if
      kp=lbuf
      return
c---
      entry gtseco
c---
c someone called xinird.  this will keep gtbuf from also reading the
c command line.
c---
      ifirst=0
      return
c*********
      entry ldbuf(cmd, ncmd, ierld)
c---
c load commands into the queue without prompting the user.
c---
c cmd(*)  i    the commands to be loaded
c ncmd    i    the number of commands to load
c ierld     o  =0 all commands loaded, otherwise the number of
c              -unloaded commands
c---
      if(nload.eq.0) nload=itop
      ierld=0
      do 820 i=1,ncmd
         itmp=nload+1
         if(itmp.gt.mxprev) itmp=1
         if(itmp.eq.itop) then
            ierld=ierld+1
         else
            nload=itmp
            cprev(nload)=cmd(i)
         end if
  820 continue
      return
c*********
      entry xgetin(cext1)
c---
c return the default extension for indirect files.
c---
      cext1=cdef
      return
c*********
      entry xnewcm
c---
c this routine should be used when the software wishes to force
c the next few lines to be read from the terminal.
c---
      level=-abs(level)
      return
c*********
      entry xrmvcm
c---
c return control to the current indirect file (if it exists).
c---
      level= abs(level)
      return
c*********
      entry xsetin(cext2)
c---
c set the default extension for indirect files.
c---
      cdef=cext2
      return
      end

c- contains entry points for
c alf
c alfsks
c fpnum
c isnum
c irange
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
      real function fpnum(ctok, ltok, ier)
      character ctok*(*)
      integer   ltok, ier
c---
c decodes number found starting at ctok(1:1).
c anything which cannot be read as part of a number terminates
c the decoding.
c---
c ctok    i
c ltok    i
c ier       o
c---
c [aft]
c---
      real       no
      parameter (no=-1.2e-34)
      character chr, cdum*20
      real      onum
      integer   icomb, kp, ibeg, ios
c---
c- check for no-data flag
      if(ltok.eq.2) then
         if(ctok(1:1).ne.'n' .and. ctok(1:1).ne.'n') goto 100
         if(ctok(2:2).ne.'o' .and. ctok(2:2).ne.'o') goto 100
         fpnum=no
         ier=0
         return
      end if
c---
  100 onum=0.
      icomb=0
      kp=1
      chr=ctok(kp:kp)
c---
  120 fpnum=0.
      ibeg=kp
      if(chr.eq.'-' .or. chr.eq.'+') then
         if(kp.ge.ltok) goto 900
         kp=kp+1
         chr=ctok(kp:kp)
      end if
c---
  130 if(chr.ge.'0' .and. chr.le.'9') then
         kp=kp+1
         chr=ctok(kp:kp)
         goto 130
      end if
c---
c- check for decimal point followed by more digits.
      if(chr.eq.'.') then
  150    kp=kp+1
         chr=ctok(kp:kp)
         if(chr.lt.'0' .or. chr.gt.'9') goto 160
         goto 150
      end if
c---
c- check for exponent.
  160 if(kp.eq.ibeg) goto 190
      call upc(chr)
      if(chr.ne.'e' .and. chr.ne.'d') goto 190
      kp=kp+1
      chr=ctok(kp:kp)
      if(chr.eq.'+' .or. chr.eq.'-') then
         kp=kp+1
         chr=ctok(kp:kp)
      end if
  180 if(chr.ge.'0' .and. chr.le.'9') then
         kp=kp+1
         chr=ctok(kp:kp)
         goto 180
      end if
c---
  190 if(kp.eq.1) goto 900
      cdum=ctok(ibeg:kp-1)
      if(kp-ibeg.lt.len(cdum)) cdum(kp-ibeg+1:)=','
      read(cdum,191,iostat=ios) fpnum
  191 format(d20.0)
      if(ios.ne.0) goto 900
c---
      if(icomb.ne.0) then
         if(icomb.eq.1) fpnum=onum+fpnum
         if(icomb.eq.2) fpnum=onum-fpnum
         if(icomb.eq.3) fpnum=onum*fpnum
         if(icomb.eq.4) fpnum=onum/fpnum
      end if
c- ignore / as last character.
      icomb=0
      if(chr.eq.'/' .and. kp.ne.ltok) goto 540
      if(chr.eq.'*') goto 530
      if(chr.eq.'-') goto 520
      if(chr.eq.'+') goto 510
      ier=0
      return
c---
  540 icomb=icomb+1
  530 icomb=icomb+1
  520 icomb=icomb+1
  510 icomb=icomb+1
      if(kp.ge.ltok) goto 900
      kp=kp+1
      chr=ctok(kp:kp)
      onum=fpnum
      fpnum=0.
      goto 120
c---
c- error within a number.
  900 ier=1
      return
      end
c*********
      integer function isnum(ctok, ltok)
      character ctok*(*)
      integer   ltok
c---
c return -1 if ctok is a number, 0 otherwise.
c---
c ctok    i
c ltok    i
c---
c [aft]
c---
      character chr
c---
      isnum=0
      if(ltok.le.0) return
      chr=ctok(1:1)
      if(chr.eq.'+' .or. chr.eq.'-' .or. chr.eq.'.' .or.
     :     (chr.ge.'0' .and. chr.le.'9')) then
c- first digit ok, check last digit
         chr=ctok(ltok:ltok)
         if(chr.ge.'0' .and. chr.le.'9') goto 900
         if(ltok.gt.1 .and. chr.eq.'.') goto 900
      end if
c- now check for no keyword.
      if(ltok.ne.2) return
      if(chr.ne.'n' .and. chr.ne.'n') return
      chr=ctok(2:2)
      if(chr.ne.'o' .and. chr.ne.'o') return
c---
  900 isnum=-1
      return
      end
c*********
      subroutine irange(ctok, ltok, ilodef, ihidef, ilo, ihi, ier)
      character ctok*(*)
      integer   ltok, ilodef, ihidef, ilo, ihi, ier
c---
c parse a range from the token.
c---
c ctok    i
c ltok    i
c ilodef  i
c ihidef  i
c ilo       o
c ihi       o
c ier       o
c---
c [aft]
c---
      real      fpnum
      integer   itmp
c---
      ilo=ilodef
      ihi=ihidef
      itmp=index(ctok,'..')
      if(itmp.gt.0) then
         if(itmp.gt.1) then
            ilo=fpnum(ctok,itmp-1,ier)
         end if
         if(itmp.lt.ltok-1) then
            ihi=fpnum(ctok(itmp+2:),ltok-itmp-1,ier)
         end if
      else if(ltok.gt.0) then
         ilo=fpnum(ctok,ltok,ier)
         ihi=ilo
      end if
      return
      end
      integer function lenact(cbuf)
      character cbuf*(*)
c---
c function to return the active length of a character string, not
c counting any trailing blanks.  n.b. an all blank string will
c return zero as the length.
c---
c cbuf    i    string whose length is to be measured.
c---
c 1988-jun-13 - standard fortran version [aft]
c---
      integer   i
c---
      do 190 i=len(cbuf),1,-1
         if(cbuf(i:i).ne.' ') then
            lenact=i
            return
         end if
  190 continue
      lenact=0
      return
      end
      subroutine script(ifunc, rbuf, nbuf, cbuf, lbuf)
      integer   ifunc, nbuf, lbuf
      real      rbuf(*)
      character cbuf*(*)
c---
c device driver for the script file.
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
      data      lunloc/0/, cdnam/'script.xco'/
c---
   11 format(a)
c---
      goto( 10, 20, 30, 40, 50) ifunc
  900 write(*,901) ifunc
  901 format('unimplemented function in script device driver:',i5)
      nbuf = -1
      return
      
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
         call xtend(ctmp,'xco')
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

      subroutine xwrcmd(level, cprom, lprom, cbuf, lbuf)
      integer   level, lprom, lbuf
      character cprom*(*), cbuf*(*)
c---
c write the command prompt and response (if logging).
c---
c 1990-sep-25 - use log and script file device drivers [aft]
c---
      character ctmp*132, chriek*16
      real      rbuf
      integer   ltmp, nbuf, newl
      data chriek/'!!!!!!!!!!!!!!!!'/
c---
c prepare and write the log file entry.
      write(ctmp,111) cprom(:lprom),cbuf(:lbuf)
  111 format(a,1x,a)
      ltmp=0
      call logger(5, rbuf, nbuf, ctmp, ltmp)
c---
c create script string.
      ctmp=' '
      ltmp=0
      if(level.gt.0) then
         ltmp=min(max(1,level),len(chriek))
         ctmp(:ltmp)=chriek(:ltmp)
      end if
      if(lbuf.gt.0) ctmp(ltmp+1:ltmp+lbuf)=cbuf(:lbuf)
      ltmp=ltmp+lbuf
      if(ltmp.lt.30) then
         ltmp=30
      else
         ltmp=ltmp+1
      end if
      ctmp(ltmp+1:ltmp+2)='! '
      ltmp=ltmp+2
      if(lprom.gt.0) then
         newl=min(len(ctmp),ltmp+lprom)
         ctmp(ltmp+1:newl)=cprom(:lprom)
         ltmp=newl
      end if
      call script(5, rbuf, nbuf, ctmp, ltmp)
      return
      end

	subroutine xparse(delstr)
c		xparse	rashafer	2 mar 86
c	xparse initialization subroutine - used to modify the delimeter values
	implicit none
c	delstr	c*	i: the string that modifies the delimeters used
c			by the various xparse subroutines, a non-blank char
c			indicates that that particular char should be
c			modified
	external xprsbd
	character*(*) delstr
	integer lenact,lenstr
	return
	end

      subroutine edicom(cbuf, lbuf)
      character cbuf*(*)
      integer   lbuf
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

c--- system dependent routines (unix version)
c--- contains entry points for:
c conc    (sub)  converts filename to system preferred case
c debugh
c dirpos  (sub)  return the number of characters in directory spec
c frelun  (sub)  free up a logical unit number
c getdat  (sub)  get current date
c getlun  (sub)  get a free logical unit number
c getran  (fun)  get a random number
c gettim  (sub)  get the current time
c getts
c getuse  (sub)  get the current user id
c gtrcln  (fun)  get the record length of a direct access file
c i2conv  (fun)  convert i*2 from ibm byte order to processor byte order
c i4conv  (fun)  convert i*4 from ibm byte order to processor byte order
c initts
c errset
c iand
c ior
c ishft
c locase  (sub)  convert to lower case
c openwr  (sub)  open wrapup
c pltter         toggle graphics/alpha mode on terminal
c prompt  (sub)  prompt
c ptend   (sub)  add prefix to file names
c rdforn  (sub)  read a foreign command, i.e., the command line
c spawn   (sub)  spawn to the operating system
c trlog   (sub)  translate logical name
c upc     (sub)  convert to upper case
c xcopy   (sub)  copy file
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
      subroutine debugh
c---
c no unix equivalent to vms command set up
c---
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
      real function getran(iseed)
      integer iseed
c---
c returns a real uniformly distributed between 0 and 1
c---
c iseed   i/o  random seed (must be integer).
c---
      integer   itmp
c---
c the following line works on the sun because the sun ignores
c integer overflows.
      iseed=1+iseed*69069
c prepare to convert to a 24-bit floating point number.  do the division
c with integers which truncate, rather than floating point, which on
c the sun would round the number.
      itmp=iseed/256
c remove the twos-complement sign.
      if(itmp.lt.0) itmp=itmp+2**24
c getran should now never be rounded to 1.0.
      getran=float(itmp)/(2.**24)
      return
      end
c*********
      subroutine gettim (ctime)
      character*8 ctime
c---
c return time in format hh:mm:ss
c---
c ctime     o  the current time
c---
      character udate*24
c---
c      call fdate(udate)
      ctime=udate(12:19)
      return
      end
c*********
c*********
      subroutine getuse(cuser, ierr)
      character cuser*(*)
      integer   ierr
c---
c return the name of the current user
c---
c cuser     o  the current user id
c ierr      o  =1 if valid, <>0 could not generate user id
c---
      integer   lenact
      character*4 user
      data user/'user'/
c---
      cuser = ' '
c      call getenv(user, cuser)
      if(lenact(cuser) .eq. 0) then
         ierr=1
      else
         ierr=0
      end if
      return
      end
c*********
      integer function gtrcln(cfile)
      character cfile*(*)
c---
c get the record length (in bytes) of a direct access unformatted file
c assuming that the first sixteen bytes are :
c		c*4	machine id ('sun ' or 'vax ')
c		i*4	check integer
c		r*4	check real
c		i*4	number of bytes per record
c---
c cfile   i    the filename
c gtrcln  o    if successful = # bytes per record
c              if failed to open file = 0
c              if failed integer check = -1
c              if failed real check = -2
c---
      integer fidint
      real    fidrl
      parameter (fidint=123456, fidrl=1.23456e+07)
      real    chkrl
      integer chkint, nbytes, ilun, ierr
      character cmach*4

      call getlun(ilun)
      call openwr(ilun, cfile, 'old', 'd', ' ', 16, 1, ierr)
      if (ierr .ne. 0) then
         gtrcln=0
         return
      end if
      read(ilun,rec=1,iostat=ierr) cmach, chkint, chkrl, nbytes
      close(ilun)
      call frelun(ilun)
      if (chkint .ne. fidint) then
         gtrcln=-1
         return
      end if
      if (chkrl .ne. fidrl) then
         gtrcln=-2
         return
      end if
      gtrcln=nbytes
      return
      end
c
c*********
      integer function i4conv(ishf)
      integer ishf
c---
c convert an integer from ibm byte order to processor byte order.
c---
c ishf    i    the number to be converted
c---
      i4conv=ishf
      return
      end
c*********
      subroutine initts
c---
c initialize the runtime statistics
c---
      integer   time
      real      tarray(2), etime, begcpu, rstart
      integer   istart, itbeg
      save istart, rstart
c---
c      istart = time()
c      rstart = etime(tarray)
      return
c
      entry sttime(itbeg, begcpu)
c return the system time when the program was started.
      itbeg =istart
      begcpu=rstart
      return
      end
c*********
      subroutine errset(ier, q1, q2, q3, q4, imxlim)
      integer   ier, imxlim
      logical   q1, q2, q3, q4
c---
c wrapup routine to load a dummy vms errset routine.
c---
      return
      end
c********
c*********
      subroutine locase(cbuf)
      character cbuf*(*)
c---
c converts cbuf to lower case.
c---
c cbuf    i/o  the character array to convert to lower case
c---
      integer   i, itmp
c---
      do 120 i=1,len(cbuf)
         itmp=ichar(cbuf(i:i))
         if(itmp.gt.64 .and. itmp.lt.91) cbuf(i:i)=char(itmp+32)
  120 continue
      return
      end
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
      subroutine pltter(ctype)
      character ctype*(*)
c---
c ctype='a' switches terminal into alpha mode.
c ctype='g' switches terminal into graphics mode.
c note: caps.
c---
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
      subroutine spawn(cbuf, lbuf, ier)
      character cbuf*(*)
      integer   lbuf, ier
c---
c spawn to operating system.  if lbuf=0 then this routine should
c spawn a shell and leave the user in the shell until the user logs
c out or exits the shell.  if lbuf<>0 then the system should only
c execute that one command and immediately return to the calling
c routine.
c---
c cbuf    i    the system command to execute
c lbuf    i    the number of valid characters in cbuf (can be zero)
c ier       o  =0 spawn was successful, <>0 otherwise
c---

      if(lbuf.gt.0) then
c         call system(cbuf(:lbuf))
      else
         write(*,*) 'type exit to return.'
c         call system('/bin/csh')
      end if
      return
      end
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
      subroutine xcopy(ifile, ofile)
c     fwj haberl   25-oct-1990 
c     copy file 'ifile' to 'ofile'
      character*(*) ifile, ofile
      
      character*512 string
      integer lstring, lenact, ierr

      string = 'cp '//ifile(1:lenact(ifile))//' '//
     &                      ofile(1:lenact(ifile))
      lstring = lenact(string)
      call spawn(string, lstring, ierr)
      return
      end

      subroutine xsquez(string, rmchar, nsave, length)
      character string*(*)
      character rmchar*1
      integer   nsave, length
c---
c subroutine to take a string and xsquez out any blanks
c** this is pretty grungy code.  surely we can do better than this
c---
c string  i/o  string to have chars removed removed
c rmchar  i    the character to be removed
c nsave   i    all strings of char of length nsave are to be
c              retained.  to remove all occurances nsave <=0
c length    o  actual length of xsquezed string
c---
c 1984-jul-29 - rashafer
c 1985-mar-08 - to allow arbitrary characters to be removed, and
c               to allow an arbitrary no. of them to be retained.
c---
      integer   irm, initl, lenact, initb, ival, icount
      logical   qon
c---
      irm=ichar(rmchar(1:1))
      initl=lenact(string)
      initb=index(string(:initl),rmchar(1:1))
      if(initb.eq.0) then
         length=initl
      else
         qon=.false.
         length=initb-1
         initb=initb-1
  100    if(initb.lt.initl) then
            initb=initb+1
            ival=ichar(string(initb:initb))
            if(ival.ne.irm) then
               qon=.false.
               length=length+1
               string(length:length)=char(ival)
            else
               if(.not.qon) then
                  qon=.true.
                  icount=1
               else
                  icount=icount+1
               end if
               if(icount.le.nsave) then
                  length=length+1
                  string(length:length)=rmchar(1:1)
               end if
            end if
            goto 100
         end if
         if(length.lt.initl) then
             string(length+1:initl)= ' '
         end if
      end if
      return
      end
	subroutine xpmsg3(descr,valstr)
c		subroutine to write a message that a number is out of range
c		in the xprmsg string in the common block
	implicit none
	character*(*) descr, valstr(3)
c	include 'xparinc.inc'
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
	integer i,retlen(3),lenact
	do i=1,3
	    call xsquez(valstr(i),' ',0,retlen(i))
	    end do
	xprmsg = 'the value, '//valstr(1)(:retlen(1))//
     &	 ', of '//descr(:lenact(descr))//
     &	 ' is outside the allowed range ('
     &	//valstr(2)(:retlen(2))//','//valstr(3)(:retlen(3))//'):'
	return
	end
      subroutine xclslg(cdisp)
      character cdisp*(*)
c---
c xparse subroutine to close (and dispose) of the logfile.
c---
c cdisp  i    the disposition string 'save'/'delete'/' '
c---
c 1990-sep-12 - new version to use logger driver - [aft]
c---
      real     rbuf
      integer  ltmp, nbuf
c---
      call logger(4, rbuf, nbuf, cdisp, ltmp)
      return
      end
	logical function xisnum(string, dval, qint, ival)
c	xisnum	rashafer 6 april 1984
c	modified from retnum (rick utility) 8 march 1986 for xparse inclusion
c		subroutine to return the numerical value of a string using
c		the fortran translation routines
	implicit none
c	** returned value
c	xisnum	l4	-: if true, a real number was parsed
c
c	string	c*	i: string to be translated
c	dval	r8	r: the returned value (modified only if xisnum is true)
c	qint	l4	r: if true, the string was an integer
c	ival	i4	r: the returned integer value (modified if qint is true)
	character*(*) string
	logical qint
	real*8 dval
	integer ival
c
	integer lenact,len,ios
	integer ivaltm
	real*8 dvaltm
	len=lenact(string)
	if(len.le.0)then
		xisnum=.false.
		qint=.false.
		return
		end if
	read(string,'(f15.0)',iostat=ios)dvaltm
	xisnum=ios.eq.0
	if(xisnum)then
		dval=dvaltm
		read(string,'(i8)',iostat=ios)ivaltm
		qint=ios.eq.0
		if(qint)then
		    ival=ivaltm
		    end if
	else
		qint=.false.
		end if
	return
	end
	integer function ixplwr (ich)
c		rashafer 19 oct 1985
c	xparse subroutine to convert upper to lower case by the value of a
c	singel character
c ***********************************************************************
c *               special vax version (ascii)              		*
c ***********************************************************************
	implicit none
	integer ich
	integer iua,iuz,iconv
	iua=ichar('A')
        iuz=ichar('Z')
	iconv=ichar('a')-ichar('A')
	if((ich.ge.iua).and.(ich.le.iuz))then
	    ixplwr = ich + iconv
	else
	    ixplwr = ich
	    end if
	return
	end

