#   Relaxation related Tcl routines
#   Copyright (C) 2006-2009 Zdenek Tosner
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
#   Tcl routines for relaxation module
#      - evaluation of the 'relax' section in the input file
#      

proc relax { data } {
   global relax par
   
   set gm_ok { none spherical_top symmetric_top asymmetric_top }
   set lm_ok { rigid model_free model_free_ext diffusion_on_a_cone diffusion_in_a_cone 3_sites_jump }
   set cc_ok { shift dipole quadrupole }
   set auto 0
   set cross 0
   
   set data [split $data "\n"]
   foreach lst $data {
     if [string match #* [string trimleft $lst]] continue
     set nam [lindex $lst 0]
     if ![string length $nam] continue

     switch -exact $nam {
         global_motion {
	                set dum [lindex $lst 1]
                        if {[lsearch -exact $gm_ok $dum]==-1} {
			   puts stderr "error: unknown global motion '$dum' in relax section, must be one" 
			   set f [join $gm_ok {, }]
			   puts stderr "       of $f"
			   exit
			}
			if [info exists relax($nam)] {
                          puts stderr "error: '$nam' already exists in array relax"
                          exit
                        }
	                set relax($nam) [lrange $lst 1 9999]
		       }
	 shift {
	        set dum "$nam\_[lindex $lst 1]"
		set dum2 [lindex $lst 2]
		if {[lsearch -exact $lm_ok $dum2]==-1} {
		   set dum1 [regsub -all {_} $dum { }]
	           puts stderr "error: unknown local motion '$dum2' for '$dum1' in relax section, " 
	           puts stderr "       must be one of [join [lrange $lm_ok 0 3] {, }],"
		   puts stderr "       [join [lrange $lm_ok 4 10] {, }]"
		   exit
		}
		if [info exists relax($dum)] {
                  puts stderr "error: '$dum' already exists in array relax"
                  exit
                }
		set relax($dum) [lrange $lst 2 9999]
		incr auto
	       }
	 dipole {
	         set dum "$nam\_[lindex $lst 1]_[lindex $lst 2]"
		 set dum2 [lindex $lst 3]
		 if {[lsearch -exact $lm_ok $dum2]==-1} {
		    set dum1 [regsub -all {_} $dum { }]
	            puts stderr "error: unknown local motion '$dum2' for '$dum1' in relax section, " 
		    puts stderr "       must be one of [join [lrange $lm_ok 0 3] {, }],"
		    puts stderr "       [join [lrange $lm_ok 4 10] {, }]"
		    exit
		 }
		 if [info exists relax($dum)] {
                    puts stderr "error: '$dum' already exists in array relax"
                    exit
                 }
		 set relax($dum) [lrange $lst 3 9999]
		 incr auto
	        }
         quadrupole {
	             set dum "$nam\_[lindex $lst 1]"
		     set dum2 [lindex $lst 2]
		     if {[lsearch -exact $lm_ok $dum2]==-1} {
		        set dum1 [regsub -all {_} $dum { }]
	                puts stderr "error: unknown local motion '$dum2' for '$dum1' in relax section, " 
	                puts stderr "       must be one of [join [lrange $lm_ok 0 3] {, }],"
		        puts stderr "       [join [lrange $lm_ok 4 10] {, }]"
		        exit
		     }
		     if [info exists relax($dum)] {
                       puts stderr "error: '$dum' already exists in array relax"
                       exit
                     }
		     set relax($dum) [lrange $lst 2 9999]
		     incr auto
	            }
         random_field {
	               set dum "$nam\_[lindex $lst 1]"
		       if [info exists relax($dum)] {
                          puts stderr "error: '$dum' already exists in array relax"
                          exit
                       }
		       set relax($dum) [lrange $lst 2 9999]
		       incr auto
		      }
         cross_correlation {
                            puts stderr "error: cross correlation not implemented yet!!!"
			    exit
	                    if [info exists relax($nam)] {
                               puts stderr "error: '$nam' already exists in array relax"
                               exit
                            }
	                    set relax($nam) [lrange $lst 1 9999]
			    incr cross
			   }
         default {
	          puts stderr "error: unknown name '$nam' in relax section, must be one"
		  puts "       of global_motion, shift, dipole, quadrupole, random_field,"
		  puts "       cross_correlation"
		  exit
		 }
     
     }
   }
   
   # default setting for global motion (none)
   if {[lsearch -exact [array names relax] global_motion] == -1} {
     set relax(global_motion) none
   }
   set relax(auto)  $auto
   set relax(cross) $cross

}
