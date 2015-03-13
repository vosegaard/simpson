#   Data manipulation Tcl-routines
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
#
#   Tcl routines in addition to the ones implemented in ftools.c
#

proc fexpr {desc reexpr imexpr} {  
   set reex [concat expr $reexpr]
   set imex [concat expr $imexpr]
   set np [fget $desc -np]
   set ni [fget $desc -ni]
   if {$ni > 1} {set np [expr $np*$ni]}
   for {set i 1} {$i <= $np} {incr i} {
      set c [findex $desc $i]
      set re [lindex $c 0]
      set im [lindex $c 1]
      set kre [eval $reex]
      set kim [eval $imex]
      fsetindex $desc $i $kre $kim
   }
}

