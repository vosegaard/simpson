#   Miscellaneous Tcl routines
#   Copyright (C) 1999 Mads Bak
#
#   This file is part of the SIMPSON General NMR Simulation Package
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version. 
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# Given the three principal elements of the chemical shift tensor
# in ppm defined on the shift/spectrometer/deshielding scale,
# returns the isotropic component, the anisotropic component and the asymmetry. 

proc csapar {s1 s2 s3} {
  set iso [expr ($s1+$s2+$s3)/3.0]

  set v [list [expr abs($s1-$iso)] [expr abs($s2-$iso)] [expr abs($s3-$iso)]]
  set a [list $s1 $s2 $s3]
    
  if [expr [lindex $v 0] > [lindex $v 1] ] {
    set v [list [lindex $v 1] [lindex $v 0] [lindex $v 2]]
    set a [list [lindex $a 1] [lindex $a 0] [lindex $a 2]]
  }

  if [expr [lindex $v 1] > [lindex $v 2] ] {
    set v [list [lindex $v 0] [lindex $v 2] [lindex $v 1]]
    set a [list [lindex $a 0] [lindex $a 2] [lindex $a 1]]
  }
  if [expr [lindex $v 0] > [lindex $v 1] ] {
    set v [list [lindex $v 1] [lindex $v 0] [lindex $v 2]]
    set a [list [lindex $a 1] [lindex $a 0] [lindex $a 2]]
  }

  if [expr [lindex $v 1] > [lindex $v 2] ] {
    set v [list [lindex $v 0] [lindex $v 2] [lindex $v 1]]
    set a [list [lindex $a 0] [lindex $a 2] [lindex $a 1]]
  }
  set z [lindex $a 2]
  set y [lindex $a 0]
  set x [lindex $a 1]
  
  set aniso [expr $z-$iso ]
  set eta [expr double($y-$x)/$aniso ];
  return [list $iso $aniso $eta]
}

proc csaprinc {iso aniso eta} {
  set zz [expr $aniso + $iso]
  set xx [expr $iso-$aniso*(1.0+$eta)/2.0]
  set yy [expr $xx + $eta*$aniso]
  return [list $xx $yy $zz]
}


proc putmatrix {m {fm "%9.3g"}} {
   foreach i $m {
     foreach j $i {
        if {[llength $j] == 2} {
          puts -nonewline [format "($fm,$fm) " [lindex $j 0] [lindex $j 1]]
        } else {
          puts -nonewline [format $fm $j]
        }
     }
     puts ""
   }
}



proc contourplot {file xlabel ylabel} {
  set f [open $file.gnu w]
  puts $f "
  set param 
  set view 0,0,1
  set cntrparam bspline
  set cntrparam levels 10
  set nosurface
  set xlabel '$xlabel'
  set ylabel '$ylabel'
  set contour
  set term post
  set output '$file.ps'
  splot '$file' w l
  "
  close $f
  exec gnuplot $file.gnu
  puts "Generated: $file.ps"
}

proc 2dplot {file xlabel ylabel {title {}}} {

  set f [open $file.gnu w]
  puts $f "
  set term post
  set param
  set view 75,20,1
  set contour
  set title  '$title'
  set xlabel '$xlabel'
  set ylabel '$ylabel'
  set output '${file}.ps'
  plot '$file' w l
  "
  close $f
  exec gnuplot $file.gnu
  puts "Generated: ${file}.ps"
}


proc 3dplot {file xlabel ylabel {zrange {}}} {

  set zrng {}
  if {[llength $zrange] == 2} {
     set zrng "set zrange \[[join $zrange :]\]"
  }
  set f [open $file.gnu w]
  puts $f "
  set term post
  set param
  $zrng
  set view 75,20,1
  set contour
  set xlabel '$xlabel'
  set ylabel '$ylabel'
  set output '${file}-3d.ps'
  splot '$file' w l
  "
  close $f
  exec gnuplot $file.gnu
  puts "Generated: ${file}-3d.ps"
}
