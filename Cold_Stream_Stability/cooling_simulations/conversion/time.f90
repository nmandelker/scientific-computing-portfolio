program main
implicit none

	integer :: noutput1, noutput2, i, Noutputs
	character(len=1024) :: filename, dir_name, command, input_arg, temp_label
	real(8) :: t

	if( iargc() .ge. 1 ) then
		call getarg(1,input_arg)
		read(input_arg,*) noutput1
		if( iargc() .ge. 2 ) then
			call getarg(2,input_arg)
			read(input_arg,*) noutput2
		else
			noutput2 = 1
		end if
	else
		noutput1 = 1
	end if
	Noutputs = noutput2 - noutput1 + 1
	print *, 'noutput1, noutput2, Noutputs'
	print *, noutput1, noutput2, Noutputs

	open(unit=16,file='./output/time.txt',form='formatted')
	write(16,'(i3.3)') Noutputs

	do i=noutput1,noutput2
		write(filename,'(a,i5.5,a,i5.5,a)') './output/output_',i,'/info_',i,'.txt'
		open(unit=17,file=filename,form='formatted')
		read(17,*)
		read(17,*)
		read(17,*)
		read(17,*)
		read(17,*)
		read(17,*)
		read(17,*)
		read(17,*)
		read(17,'(A13,E23.15)')temp_label,t
		close(unit=17)
		write(16,'(F6.4)') t
	end do
	close(unit=16)
end program main

