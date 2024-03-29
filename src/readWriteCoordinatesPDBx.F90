subroutine writeCoordinatesPDBx(uout)
  use moduleSystem 
  use moduleElements 
  implicit none
  integer, intent(inout) :: uout

  integer :: nb, ib
  integer :: iatm, jatm, imol
  character(cp) :: l

  write(uout,'("data_cell")')
  write(uout,'("_cell.length_a",1x,f10.4)')frame % cell(1)
  write(uout,'("_cell.length_b",1x,f10.4)')frame % cell(2)
  write(uout,'("_cell.length_c",1x,f10.4)')frame % cell(3)
  write(uout,'("_cell.angle_alpha",1x,f10.4)')frame % cell(4)
  write(uout,'("_cell.angle_beta",1x,f10.4)')frame % cell(5)
  write(uout,'("_cell.angle_gamma",1x,f10.4)')frame % cell(6)
  write(uout,'("#")')

  write(uout,'("loop_")')
  write(uout,'("_struct_conn.id")')
  write(uout,'("_struct_conn.conn_type_id")')
  write(uout,'("_struct_conn.ptnr1_label_asym_id")')
  write(uout,'("_struct_conn.ptnr1_label_comp_id")')
  write(uout,'("_struct_conn.ptnr1_label_seq_id")')
  write(uout,'("_struct_conn.ptnr1_label_atom_id")')
  write(uout,'("_struct_conn.ptnr2_label_asym_id")')
  write(uout,'("_struct_conn.ptnr2_label_comp_id")')
  write(uout,'("_struct_conn.ptnr2_label_seq_id")')
  write(uout,'("_struct_conn.ptnr2_label_atom_id")')

  do ib=1,numberOfUniqueBonds
    iatm = listOfUniqueBonds(1,ib)
    jatm = listOfUniqueBonds(2,ib)
    imol = atomToMoleculeIndex (iatm)
    write(uout,'("bond",i0," covale",2(1x,"A",1x,a,1x,i5,1x,a))') ib, &
      listOfMolecules(imol) % resname, imol, frame % lab(iatm), &
      listOfMolecules(imol) % resname, imol, frame % lab(jatm) 
  end do
  write(uout,'("#")')

  write(uout,'("loop_")')
  write(uout,'("_atom_site.group_PDB")')
  write(uout,'("_atom_site.id")')
  write(uout,'("_atom_site.type_symbol")')
  write(uout,'("_atom_site.label_atom_id")')
  write(uout,'("_atom_site.label_comp_id")')
  write(uout,'("_atom_site.label_asym_id")')
  write(uout,'("_atom_site.label_seq_id")')
  write(uout,'("_atom_site.Cartn_x")')
  write(uout,'("_atom_site.Cartn_y")')
  write(uout,'("_atom_site.Cartn_z")')
  
  do iatm=1,frame%natoms
    imol = atomToMoleculeIndex (iatm)
    l = frame%lab(iatm)
    write(uout,'("ATOM",1x,i6,3(1x,a6),1x,"A",1x,i6,3(1x,f10.4))') &
      iatm, getElement(l), l, listOfMolecules(imol) % resname, imol, frame%pos(:,iatm)
  end do

end subroutine writeCoordinatesPDBx
