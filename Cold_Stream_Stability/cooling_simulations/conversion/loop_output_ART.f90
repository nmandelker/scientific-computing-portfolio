program main
implicit none

	integer :: noutput1, noutput2, i
	character(len=1024) :: dir_name1, command, input_arg

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
	print *, 'noutput'
	print *, noutput1, noutput2

	do i=noutput1,noutput2
		write(dir_name1,'(a,i5.5)') './output/output_',i

		write(command,'(a,a,a)') './make_ART_format.exe -inp ',trim(dir_name1),' -out ART_output -xmi 0.0 -xma 1.0 '
		call system(trim(command))
	end do
end program main
	
