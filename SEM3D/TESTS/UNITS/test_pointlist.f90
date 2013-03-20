program test

  use point_list

  type(PointTable) :: table
  integer :: id, i
  real, allocatable, dimension(:,:) :: res

  call init_point_table(table, 100)

  call add_point(table, 0., 0., 0., id)
  write(*,*) "0,0,0 : ", id

  call add_point(table, 0., 1., 0., id)
  write(*,*) "0,1,0 : ", id

  call add_point(table, 0., 0., 1., id)
  write(*,*) "0,0,1 : ", id

  call add_point(table, 0., 0., 0., id)
  write(*,*) "0,0,0 : ", id

  call add_point(table, 0., 1., 0., id)
  write(*,*) "0,1,0 : ", id

  do i=1,2000
      call add_point(table, i/23., i/17.+12, i/12. + 17, id)
  end do

  call fill_table_array(table, res)

  write(*,*) "TABLE OF POINTS:", size(res,2)
  do i=0,min(size(res,2)-1,100)
      write(*,*) res(:,i)
  end do
end program
