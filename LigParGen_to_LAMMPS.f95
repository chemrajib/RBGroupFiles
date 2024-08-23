! This code use the output from LigParGen and convert the lammps data 
! file to have unique atoms, bonds, etc.
! This code in fact can also be used to merge multiple lammps data files.
! Written by Rajib Biswas (2024).
module param
  implicit none
  integer, parameter :: inpf = 100

  integer :: i, j, k, l  
  
  integer :: natoms, nbonds, ndihs, nimps, nangles
  integer :: oatmtps, obtps, oangtps, odihtps, oimptps
  integer :: ntps, ntc 
  integer :: atomtps, bondtps, angletps, dihedtps, improptps, mol_indx
  integer, allocatable, dimension(:) :: angle_indx, bond_indx, dih_indx, imp_indx 
  integer, allocatable, dimension(:) :: batom1, batom2, btps
  integer, allocatable, dimension(:) :: aatom1, aatom2, aatom3, antps
  integer, allocatable, dimension(:) :: datom1, datom2, datom3, datom4, dhtps
  integer, allocatable, dimension(:) :: iatom1, iatom2, iatom3, iatom4, iptps
  integer, allocatable, dimension(:) :: a_index, m_index, oa_tps, a_tps
  integer, allocatable, dimension(:) :: oangtp_index, obondtp_indx
  integer, allocatable, dimension(:) :: angtp_indx, bondtp_indx, dihtp_indx, imptp_indx

  double precision, allocatable, dimension(:) :: charge, rx, ry, rz
  double precision, allocatable, dimension(:,:) :: omass, mass 
  double precision, allocatable, dimension(:,:) :: opaircff, oangcff, obondcff
  double precision, allocatable, dimension(:,:) :: paircff, angcff, bondcff
  double precision, allocatable, dimension(:,:) :: oimpcff, impcff
  double precision, allocatable, dimension(:,:) :: odihcff, dihcff
  !double precision, allocatable, dimension(:,:) :: dyna1, dyna2 
  
  double precision :: xlo, xhi, ylo, yhi, zlo, zhi, dr1, dr2, tol = 1.0e-6 

  logical :: check 
end module param

program lammps_out
  use param
  implicit none
  character (len=8) :: x1, fmt
  character (len=1024) :: file1, file2
  
  file2='H_096CAD'
  file1=trim(file2)//'_edited.data'
  open(unit = inpf, file=trim(file2)//'.lmp', action='read')

  ! Read the lammps data file 
  call read_lammps  

  print*, 'LAMMPS data file ', trim(file2)//'.lmp', ' Reading complete'


  ! Find the unique atomtps and pair cfficient
  call find_atomtps
  ! Find unique masses
  allocate(mass(2,atomtps))
  ntc = 1; i = 0
  mass = 0.d0
  check = .true.
  do while (check)
     i = i + 1
     if ((a_tps(i).eq.ntc).and.(mass(2,ntc).eq.0.d0)) then
        mass(1, ntc) = ntc 
        mass(2, ntc) = omass(2,i)
        !print*, i, ntc, omass(2,i)
        ntc = ntc + 1
     end if
     if (ntc.gt.atomtps) check = .false.
  end do

  ! Find unique angcff
  !call find_angtps
  !call find_bondtps
  !call find_dihtps
  call find_imptps
end program lammps_out
!=======================================================================================
! Find the unique atomtps and pair cfficient
subroutine find_atomtps
  use param
  implicit none 
  double precision, allocatable, dimension(:,:) :: dyna1, dyna2
  
  a_tps = 0; i = 1
  do  
     if (i.eq.1) then 
        ntps = 1 
        a_tps(i) = ntps
        allocate(dyna1(2,ntps))
        dyna1(1,1) = opaircff(2,i)
        dyna1(2,1) = opaircff(3,i)
     end if 
     do j = i+1, oatmtps
        do k = 1, ntps
           if ((dyna1(1,k).eq.opaircff(2,j)).and.(dyna1(2,k).eq.opaircff(3,j))) then 
              a_tps(j) = k
           end if
        end do   
     end do
        check = .true. 
        j = i
        if (j.lt.oatmtps) then 
           do while (check)
           j = j+1
           ntc = 0 
           do k = 1, ntps 
              if ((dyna1(1,k).ne.opaircff(2,j)).or.(dyna1(2,k).ne.opaircff(3,j))) then 
                 ntc = ntc + 1
              end if
           end do
           !if (j.eq.oatmtps) print*, 1, ntc, ntps 
           if (ntc.eq.ntps) then
              check = .false.
              ntc = ntps + 1
           end if 
           if (j.eq.oatmtps) exit 
        end do
        end if 
           if (.not.check) then 
           a_tps(j) = ntc 
           allocate(dyna2(2,ntc))
           do l = 1, ntc - 1
             dyna2(1,l) = dyna1(1,l)
             dyna2(2,l) = dyna1(2,l)
          end do
          dyna2(1,ntc) = opaircff(2,j)
          dyna2(2,ntc) = opaircff(3,j)
          deallocate(dyna1)
          allocate(dyna1(2,ntps))
          dyna1 = dyna2
          deallocate(dyna2)
          ntps = ntc
          end if 
          i = j
     if (i .eq. oatmtps) exit
  end do
  atomtps = ntps
  allocate(paircff(2,atomtps))
  paircff = dyna1
  deallocate(dyna1)
 end subroutine find_atomtps
 !========================================================================
 subroutine find_bondtps
  use param
  implicit none
  double precision, allocatable, dimension(:,:) :: dyna1, dyna2
  allocate(bondtp_indx(obtps))
  bondtp_indx = 0; i = 1
  !do j = 1, obtps
  !   print*, obondcff(:,j)
  !end do 
  !stop 
  do  
     if (i.eq.1) then 
        ntps = 1 
        bondtp_indx(i) = ntps
        allocate(dyna1(2,ntps))
        dyna1(1,1) = obondcff(2,i)
        dyna1(2,1) = obondcff(3,i)
     end if 
     do j = i+1, obtps
        do k = 1, ntps
           if ((dyna1(1,k).eq.obondcff(2,j)).and.(dyna1(2,k).eq.obondcff(3,j))) then 
              bondtp_indx(j) = k 
           end if
        end do   
     end do
        check = .true. 
        j = i
        if (j.lt.obtps) then 
           do while (check)
           j = j+1
           ntc = 0 
           do k = 1, ntps 
              if ((dyna1(1,k).ne.obondcff(2,j)).or.(dyna1(2,k).ne.obondcff(3,j))) then 
                 ntc = ntc + 1
              end if
           end do
           if (ntc.eq.ntps) then
              check = .false.
              ntc = ntps + 1
           end if 
           if (j.eq.obtps) exit 
        end do
        end if 
           if (.not.check) then 
           bondtp_indx(j) = ntc 
           allocate(dyna2(2,ntc))
           do l = 1, ntc - 1
             dyna2(1,l) = dyna1(1,l)
             dyna2(2,l) = dyna1(2,l)
          end do
          dyna2(1,ntc) = obondcff(2,j)
          dyna2(2,ntc) = obondcff(3,j)
          deallocate(dyna1)
          allocate(dyna1(2,ntps))
          dyna1 = dyna2
          deallocate(dyna2)
          ntps = ntc
          end if 
          i = j
     if (i .eq. obtps) exit
  end do
  bondtps = ntps
  allocate(bondcff(2, bondtps))
  bondcff = dyna1
  deallocate(dyna1)
  print*, bondtps 
  do j = 1, obtps
     print*, j, bondtp_indx(j)
  end do 
  do j = 1, bondtps
     print*, j, bondcff(1,j), bondcff(2,j)
  end do 
 end subroutine find_bondtps
!=======================================================================================
! Find unique angles
 !========================================================================
 subroutine find_angtps
  use param
  implicit none
  double precision, allocatable, dimension(:,:) :: dyna1, dyna2
  allocate(angtp_indx(oangtps))
  angtp_indx = 0; i = 1; ntps = 1
  do  
     if (i.eq.1) then 
        ntps = 1 
        angtp_indx(i) = ntps
        allocate(dyna1(2,ntps))
        dyna1(1,1) = oangcff(2,i)
        dyna1(2,1) = oangcff(3,i)
     end if 
     do j = i+1, oangtps
        do k = 1, ntps
           if ((dyna1(1,k).eq.oangcff(2,j)).and.(dyna1(2,k).eq.oangcff(3,j))) then 
              angtp_indx(j) = k 
           end if
        end do   
     end do
        check = .true. 
        j = i
        if (j.lt.oangtps) then 
           do while (check)
           j = j+1
           ntc = 0 
           do k = 1, ntps 
              if ((dyna1(1,k).ne.oangcff(2,j)).or.(dyna1(2,k).ne.oangcff(3,j))) then 
                 ntc = ntc + 1
              end if
           end do
           if (ntc.eq.ntps) then
              check = .false.
              ntc = ntps + 1
           end if 
           if (j.eq.oangtps) exit 
        end do
        end if 
           if (.not.check) then 
           angtp_indx(j) = ntc 
           allocate(dyna2(2,ntc))
           do l = 1, ntc - 1
             dyna2(1,l) = dyna1(1,l)
             dyna2(2,l) = dyna1(2,l)
          end do
          dyna2(1,ntc) = oangcff(2,j)
          dyna2(2,ntc) = oangcff(3,j)
          deallocate(dyna1)
          allocate(dyna1(2,ntps))
          dyna1 = dyna2
          deallocate(dyna2)
          ntps = ntc
          end if 
          i = j
     if (i .eq. oangtps) exit
  end do
  angletps = ntps
  allocate(angcff(2,angletps))
  angcff = dyna1
  deallocate(dyna1)
  print*, angletps 
  do j = 1, oangtps
     print*, j, angtp_indx(j)
  end do 
  do j = 1, angletps
     print*, j, angcff(1,j), angcff(2,j)
  end do 
 end subroutine find_angtps
!=============================================================================
! Find unique dihedrals
 !========================================================================
 subroutine find_dihtps
  use param
  implicit none
  double precision, allocatable, dimension(:,:) :: dyna1, dyna2
  allocate(dihtp_indx(odihtps))
  dihtp_indx = 0; i = 1; ntps = 1
  do  
     if (i.eq.1) then 
        ntps = 1 
        dihtp_indx(i) = ntps
        allocate(dyna1(4,ntps))
        dyna1(1,1) = odihcff(2,i)
        dyna1(2,1) = odihcff(3,i)
        dyna1(3,1) = odihcff(4,i)
        dyna1(4,1) = odihcff(5,i)
     end if 
     do j = i+1, odihtps
        do k = 1, ntps
           if ((dyna1(1,k).eq.odihcff(2,j)).and.(dyna1(2,k).eq.odihcff(3,j)).and. &
              (dyna1(3,k).eq.odihcff(4,j)).and.(dyna1(4,k).eq.odihcff(5,j))) then 
              dihtp_indx(j) = k 
           end if
        end do   
     end do
        check = .true. 
        j = i
        if (j.lt.odihtps) then 
           do while (check)
           j = j+1
           ntc = 0 
           do k = 1, ntps 
              if ((dyna1(1,k).ne.odihcff(2,j)).or.(dyna1(2,k).ne.odihcff(3,j)).or. &
              (dyna1(3,k).ne.odihcff(4,j)).or.(dyna1(4,k).ne.odihcff(5,j))) then 
                 ntc = ntc + 1
              end if
           end do
           if (ntc.eq.ntps) then
              check = .false.
              ntc = ntps + 1
           end if 
           if (j.eq.odihtps) exit 
        end do
        end if 
           if (.not.check) then 
           dihtp_indx(j) = ntc 
           allocate(dyna2(4,ntc))
           do l = 1, ntc - 1
              do k = 1, 4 
                 dyna2(k,l) = dyna1(k,l)
             end do
          end do
          dyna2(1,ntc) = odihcff(2,i)
          dyna2(2,ntc) = odihcff(3,i)
          dyna2(3,ntc) = odihcff(4,i)
          dyna2(4,ntc) = odihcff(5,i)
          deallocate(dyna1)
          allocate(dyna1(4,ntps))
          dyna1 = dyna2
          deallocate(dyna2)
          ntps = ntc
          end if 
          i = j
     if (i .eq. odihtps) exit
  end do
  dihedtps = ntps
  allocate(dihcff(4,dihedtps))
  dihcff = dyna1
  deallocate(dyna1)
  print*, dihedtps 
  do j = 1, odihtps
     print*, j, dihtp_indx(j)
  end do 
  do j = 1, dihedtps
     print*, j, dihcff(:,j)
  end do 
 end subroutine find_dihtps
 !=============================================================================
! Find unique imrpopers 
 !========================================================================
 subroutine find_imptps
  use param
  implicit none
  double precision, allocatable, dimension(:,:) :: dyna1, dyna2
  integer, parameter :: nprm = 3 
  allocate(imptp_indx(oimptps))
  imptp_indx = 0; i = 1; ntps = 1
  do  
     if (i.eq.1) then 
        ntps = 1 
        imptp_indx(i) = ntps
        allocate(dyna1(nprm,ntps))
        do l = 1, nprm
           dyna1(l,1) = oimpcff(l+1,i)
        end do 
     end if 
     do j = i+1, oimptps
        do k = 1, ntps
           if ((dyna1(1,k).eq.oimpcff(2,j)).and.(dyna1(2,k).eq.oimpcff(3,j)).and. &
              (dyna1(3,k).eq.oimpcff(4,j))) then 
              imptp_indx(j) = k 
           end if
        end do   
     end do
        check = .true. 
        j = i
        if (j.lt.oimptps) then 
           do while (check)
           j = j+1
           ntc = 0 
           do k = 1, ntps 
              if ((dyna1(1,k).ne.oimpcff(2,j)).or.(dyna1(2,k).ne.oimpcff(3,j)).or. &
              (dyna1(3,k).ne.oimpcff(4,j))) then 
                 ntc = ntc + 1
              end if
           end do
           if (ntc.eq.ntps) then
              check = .false.
              ntc = ntps + 1
           end if 
           if (j.eq.oimptps) exit 
        end do
        end if 
           if (.not.check) then 
           imptp_indx(j) = ntc 
           allocate(dyna2(nprm,ntc))
           do l = 1, ntc - 1
              dyna2(:,l) = dyna1(:,l)
          end do
          dyna2(:,ntc) = oimpcff(2:4,i)
          deallocate(dyna1)
          allocate(dyna1(nprm,ntps))
          dyna1 = dyna2
          deallocate(dyna2)
          ntps = ntc
          end if 
          i = j
     if (i .eq. oimptps) exit
  end do
  improptps = ntps
  allocate(impcff(nprm,improptps))
  impcff = dyna1
  deallocate(dyna1)
  print*, improptps
  do j = 1, oimptps
     print*, j, imptp_indx(j)
  end do 
  do j = 1, improptps
     print*, j, impcff(:,j)
  end do 
 end subroutine find_imptps
!=============================================================================
subroutine read_lammps
  use param
  implicit none
  read(inpf,*)
  read(inpf,*)
  
  read(inpf,*) natoms
  read(inpf,*) nbonds
  read(inpf,*) nangles
  read(inpf,*) ndihs
  read(inpf,*) nimps
  
  read(inpf,*)

  read(inpf,*) oatmtps 
  read(inpf,*) obtps 
  read(inpf,*) oangtps 
  read(inpf,*) odihtps
  read(inpf,*) oimptps 

  read(inpf,*)
  
  read(inpf,*) xlo, xhi 
  read(inpf,*) ylo, yhi
  read(inpf,*) zlo, zhi 
  
  allocate(charge(natoms), rx(natoms), ry(natoms), rz(natoms))
  allocate(a_index(natoms), m_index(natoms), a_tps(natoms),oa_tps(natoms))
  allocate(omass(2,oatmtps), opaircff(3,oatmtps))
  allocate(obondcff(3,obtps),oangcff(3,oangtps))
  allocate(odihcff(5,odihtps))
  allocate(oimpcff(4,oimptps))
  allocate(bond_indx(nbonds),angle_indx(nangles),dih_indx(ndihs),imp_indx(nimps))
  allocate(btps(nbonds),batom1(nbonds),batom2(nbonds))
  allocate(antps(nangles),aatom1(nangles),aatom2(nangles),aatom3(nangles))
  allocate(dhtps(ndihs),datom1(ndihs),datom2(ndihs),datom3(ndihs),datom4(ndihs))
  allocate(iptps(nimps),iatom1(nimps),iatom2(nimps),iatom3(nimps),iatom4(nimps))
  
  read(inpf,*)
  read(inpf,*)
  read(inpf,*)
  
  do i = 1, oatmtps
     read(inpf,*) omass(:,i)
  end do
  
  read(inpf,*)
  read(inpf,*)
  read(inpf,*)
  
  do i = 1, oatmtps
     read(inpf,*) opaircff(:,i)
  end do

  read(inpf,*)
  read(inpf,*)
  read(inpf,*)
  
  do i = 1, obtps
     read(inpf,*) obondcff(:,i)
  end do
  
  read(inpf,*)
  read(inpf,*)
  read(inpf,*)
  
  do i = 1, oangtps
     read(inpf,*) oangcff(:,i)
  end do
  
  read(inpf,*)
  read(inpf,*)
  read(inpf,*)
  
  do i = 1, odihtps
     read(inpf,*) odihcff(:,i)
  end do 
  
  read(inpf,*)
  read(inpf,*)
  read(inpf,*)
  
  do i = 1, oimptps
     read(inpf,*) oimpcff(:,i)
  end do

  read(inpf,*)
  read(inpf,*)
  read(inpf,*)
   
  do i = 1, natoms
     read(inpf,*) a_index(i), m_index(i), oa_tps(i), charge(i), rx(i), ry(i), rz(i) 
  end do
  
  read(inpf,*)
  read(inpf,*)
  read(inpf,*)
  
  do i = 1, nbonds
     read(inpf,*) bond_indx(i), btps(i), batom1(i), batom2(i)
  end do
  
  read(inpf,*)
  read(inpf,*)
  read(inpf,*)
  
  do i = 1, nangles
     read(inpf,*) angle_indx(i), antps(i), aatom1(i), aatom2(i), aatom3(i)
  end do
  
  read(inpf,*)
  read(inpf,*)
  read(inpf,*)
  
  do i = 1, ndihs
     read(inpf,*) dih_indx(i), dhtps(i), datom1(i), datom2(i), datom3(i), datom4(i)
  end do
  
  read(inpf,*)
  read(inpf,*)
  read(inpf,*)
  
  do i = 1, nimps
     read(inpf,*) imp_indx(i), iptps(i), iatom1(i), iatom2(i), iatom3(i), iatom4(i)
  end do
end subroutine 
!=====================================================================================
subroutine write_lammps(file1)
  use param
  implicit none
  character (len=1024) :: x1, file1
  open(unit=11,file=file1, form='formatted', action='write')
1 format(A55)
2 format(I6,2x,A16)
3 format(4x,F12.6,4x,F12.6,2x,2(A3,x))
4 format(A6)
5 format(I4,2x,F8.4)
  write(11,1) 'Generated with tools_lammpsout written by R. Biswas'
  write(11,*)
  write(11,2)          natoms, 'atoms'
  write(11,2)          nbonds, 'bonds'
  write(11,2)          nangles, 'angles'
  write(11,2)          ndihs, 'dihedrals'
  write(11,2)          nimps, 'impropers'
  write(11,*)
  write(11,2)          atomtps, 'atom tps'
  write(11,2)          bondtps, 'bond tps'
  write(11,2)          angletps, 'angle tps'
  write(11,2)          dihedtps, 'dihedral tps'
  write(11,2)          improptps, 'improper tps'
  write(11,*)
  write(11,3)    xlo, xhi,  'xlo', 'xhi'
  write(11,3)    ylo, yhi,  'ylo', 'yhi'
  write(11,3)    zlo, zhi,  'zlo', 'zhi'
  write(11,*)
  write(11,'(A6)')     'Masses'
  write(11,*)
  !do j = 1, atomtps
  !   write(11,'(I4,2x,F8.4,2x,A1,2x,A2)')    j, Mass(j)!, '#', a_tps(j) 
  !end do
  write(11,*)
  write(11,'(A5)') 'Atoms'
  write(11,*)
  do j = 1, natoms
     write(11,'(3(I6,2x),4(F10.5,2x))') a_index(j), m_index(j), a_tps(j), charge(a_tps(j)), rx(j), ry(j), rz(j)
  end do
  write(11,*)
  write(11,'(A5)') 'Bonds'
  write(11,*)
  do j = 1, nbonds
     write(11,'(4(I6,2x))') j, btps(j), batom1(j), batom2(j)
  end do
  write(11,*)
  write(11,'(A6)') 'Angles'
  write(11,*)
  do j = 1, nangles
     write(11,'(5(I6,2x))') j, antps(j), aatom1(j), aatom2(j), aatom3(j)
  end do
  close(11)
end subroutine write_lammps
