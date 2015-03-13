
#   SIMPSON Tcl routines
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
#   Tcl routines central to the SIMPSON simulation
#      - evaluation of the 'spinsys' and 'par' sections in the input file
#      - default values for the 'main' and 'pulseq' procedures
#      - fsimpson and simpson procedures that performs the simulation
#      

proc markvars {} {
      uplevel #0 {
        set _simpson_omitvars(_simpson_omitvars) 1
        set _simpson_omitvars(_slave_state) 1
        set _simpson_omitvars(_globs) 1
        set _simpson_omitvars(__name) 1
        foreach __name [info vars] {
           set _simpson_omitvars($__name) 1
        }
      }
} 
    
proc savestate {} {
      uplevel #0 {
      if [info exists _slave_state] {unset _slave_state}
      set _slave_state {}
      set _globs "global"
        foreach __name [info globals]  {
          if [info exists _simpson_omitvars($__name)] continue 
          lappend _globs $__name
          if [array exists $__name] {
             lappend _slave_state  [list array set $__name [array get $__name]]
          } else {
               lappend _slave_state  [list set $__name [set $__name]]
          }
        }
        return [join [list [join $_globs " "] [join $_slave_state \n] ] \n]
      }
}

proc signalhandler {code message} {
  global stop sigwascalled par

  if [info exists sigwascalled] {
     puts stderr "Program is aborted."
     #if [file exists $par(name).lock] {
     #   file delete $par(name).lock
     #}
     exit
  }
  set stop 1
  puts stderr "Received '$message', setting flag to terminate."
  puts stderr "Repeat action to abort program immediately."
  set sigwascalled 1
}

# Sets array 'varname' to the contents of data
# List elements in data are
#      elemname value
# or   elemname {value value value}

proc setvar {varname data {overwrite 0}} {
  global $varname
  set var $varname

  set data [split $data "\n"]

  foreach lst $data {
    if [string match #* [string trimleft $lst]] continue
    set nam [lindex $lst 0]
    if ![string length $nam] continue
    if !$overwrite {
      if [info exists ${var}($nam)] {
        puts stderr "error: '$nam' already exists in array ${var}"
        exit
      }
    }
    set ${var}($nam) [lrange $lst 1 9999]
  }
}

# This is the template for specifying the spin interactions
# It defines the names of the parameters, N is a number

setvar ssnam {
  shift N iso aniso eta alpha beta gamma
  jcoupling N N iso aniso eta alpha beta gamma
  quadrupole N order aniso eta alpha beta gamma
  dipole N N aniso alpha beta gamma
  dipole_ave N N aniso eta alpha beta gamma
}

# Substitutes variable names with the corresponding values
# found in 'ssval'
# resfreq is the resonance-frequency (in Hz) for the nucleus in question
# or zero if convertion (e.g. 100p ==> 100*resfreq)

proc ssSubstExpr {ex resfreq} {
  global ssval par

  if [regexp {^[0-9.Ee+-]+$} $ex] {
    return $ex
  }
  set origex $ex
  if [regsub {([0-9.Ee+-]+)p} $ex "(\\1*($resfreq))" dummy] {
     if {$resfreq == 0} {
        puts stderr "error in expression '$origex'. Can only use ppm to hz convertion for chemical shift."
        exit
     }
     if ![info exists par(proton_frequency)] {
        puts stderr "error in expression '$origex'. 'proton_frequency' must be set when converting ppm to Hz values."
        exit
     }
     set freq [expr abs(($par(proton_frequency)/1.0e6)*($resfreq/1.0e6))]
     if ![regsub {([0-9.Ee+-]+)p} $ex "(\\1*($freq))" ex] {
        puts stderr "error: illegal value of 'proton_frequency' = $par(proton_frequency) "
        exit
     }
  }
  set i 0
  while {[regsub -all {([a-z]+(_[0-9]+)+_[a-z]+)} $ex {$ssval(\1)} ex] != 0} {
    set ex [subst $ex]
    if {[incr i] > 100} {
     puts stderr "error: substitution over 100 times af the expression '$origex'"
     puts stderr "       have you made any circular references ?"
     exit
    }
  }    
  if [catch {set ex [expr $ex]} res] {
    puts stderr "error: $res"
    exit  
  }
  return $ex
}

# Sets variables in array 'ssval' corresponding to values read from 'ss'
# with the names taken from the namelist 'ssnam'

proc ssSetValues {} {
  global spinsys ssval ssnam

  foreach iact [array names spinsys] {
    set val $spinsys($iact)
    if {[string compare $iact "nuclei"] == 0} {
      set ssval(nuclei) $val
      continue
    }
    if {[string compare $iact "channels"] == 0} {
      set ssval(channels) $val
      continue
    }
    regsub {[X]+$} $iact {} iact

    if ![info exists ssnam($iact)] {
      puts stderr "error: in spinsys section at '$iact', correct syntax is:\n  spinsys {"
      foreach i [array names ssnam] {
        puts stderr "    $i $ssnam($i)"
      }
      puts stderr "  }"
      exit
    }
    set nam $ssnam($iact)

    if {[llength $val] != [llength $nam]} {
       puts stderr "error: in spinsys: input list    '$iact $val'"
       puts stderr "       does not match definition '$iact $nam'"
       exit
    }
    set nuc ""
    set k 1
    foreach j $nam {
      if {[string compare $j "N"] == 0} {
        set nuc "${nuc}_[lindex $val [expr $k - 1]]"
      } else {
        set ssval($iact${nuc}_$j) [lindex $val [expr $k - 1]]
      }
      incr k
    }
  }
}

# The inverse operation of ssSetNames
# Sets the variables in array 'spinsysres' corresponding to values read from 'ssval'
# with the names taken from the namelist 'ssnam'

proc ssSetSpinsys {} {
  global spinsys spinsysres ssval ssnam

  foreach iact [array names spinsys] {
    set orignam $iact
    set val $spinsys($orignam)

    if {[string compare $iact "nuclei"] == 0} {
      set spinsysres(nuclei) $ssval(nuclei)
      continue
    }
    if {[string compare $iact "channels"] == 0} {
      set spinsysres(channels) $ssval(channels)
      continue
    }
    regsub {[X]+$} $iact {} iact
    set nam $ssnam($iact)
    if {[llength $val] != [llength $nam]} {
       puts stderr "error: input list '$val' for interaction '$iact' doesn't match definition '$nam'"
       exit
    }
    set lst {}
    set nuc ""
    set k 1
    foreach j $nam {
      set ival [lindex $val [expr $k - 1]]
      if {[string compare $j "N"] == 0} {
        set nuc "${nuc}_$ival"
        lappend lst $ival
      } else {
        lappend lst $ssval($iact${nuc}_$j)
      }
      incr k
    }
    set spinsysres($orignam) $lst
  }
}


# Sets the array 'spinsys' with the values specified in the 'data' list
# References to each other and expressions is resolved
# Warning: circular references are allowed because

proc spinsys_resolve { { fitval {} } } {
  global spinsys ssval spinsysres par

  ssSetValues
  foreach i $fitval {
    set nam [lindex $i 0]
    set val [lindex $i 1]
    if ![info exists ssval($nam)] {
      puts stderr "error: variable '$nam' does not correspond to any variable in section 'spinsys'"
      exit
    }
    set ssval($nam) $val      
  }

  ssSetSpinsys  
  if ![info exists spinsysres(nuclei)] {
    puts stderr "error: spinsys has to include 'nuclei' definition"
    exit
  }
  set k 1
  foreach j $spinsysres(nuclei) {
    set nuc($k) $j
    incr k
  }
  foreach i [array names spinsysres] {  
    if ![string compare $i nuclei] {
      continue
    }
    if ![string compare $i channels] {
      continue
    }
    set lst {}
    set resfreq 0
    if ![string first shift $i ] {
       set resfreq [resfreq $nuc([lindex $spinsysres($i) 0]) 1e6]
    }
    set k 1
    foreach j $spinsysres($i) {
      if {$k == 2 || $k == 3} {
        lappend lst [ssSubstExpr $j $resfreq]
      } else {
        lappend lst [ssSubstExpr $j 0]
      }
      incr k
    }
    set spinsysres($i) $lst
  }
}

proc spinsys { data } {
  global spinsys

  set okpar { channels nuclei dipole quadrupole shift jcoupling mixingterm dipole_ave }

  set data [split $data "\n"]
  foreach lst $data {
    if [string match #* [string trimleft $lst]] continue
    set nam [lindex $lst 0]
    if ![string length $nam] continue

    if {[lsearch -exact $okpar $nam] == -1} {
      puts stderr "error: unknown name '$nam' in spinsys section, must be one"
      set f [join $okpar {, }]
      puts stderr "       of $f"
      exit
    }
    
    while {[info exists spinsys($nam)]} {
      set nam "${nam}X"
    }
    set spinsys($nam) [lrange $lst 1 9999]
  }
}



proc cmp_length {a b} {
  return [expr [string length $a] < [string length $b]]
}

proc par {data} {
  global par

  set okpar {
    proton_frequency spin_rate sw sw1 np ni method rotor_angle
    gamma_angles fixed_rep real_spec block_diag detect_operator
    crystal_file start_operator name verbose various variable pulse_sequence
    conjugate_fid dipole_check gamma_zero use_cluster new_cluster cluster_port
    inner_rotor_angle outer_rotor_angle inner_spin_rate outer_spin_rate dor
    string oc_tol_cg oc_tol_ls oc_mnbrak_step oc_max_iter oc_cutoff 
    oc_cutoff_iter oc_var_save_iter oc_var_save_proc oc_cg_min_step oc_max_brack_eval
    oc_max_brent_eval oc_verbose oc_grad_level oc_method
    rfprof_file use_3_angle_set acq_adjoint
    zprofile zvals relax prop_method use_sparse num_cores averaging_file
    points_per_cycle ED_symmetry
    oc_lbfgs_eps oc_lbfgs_tol_ls oc_lbfgs_max_ls_eval oc_lbfgs_m
    sparsity sparse_tol maxfulldim maxdimdiagonalize
  }
  # ZT: added 'string', 'rfprof_file' and OC_ related to allowed parameters (above)
  # ZT: relax decides wheather relaxation is invoked
  # ZT: prop_method is for which method to use for propagator calculation

  set allowsubst {
    proton_frequency spin_rate sw sw1 np ni rotor_angle
    gamma_angles fixed_rep real_spec block_diag gamma_zero
    inner_rotor_angle outer_rotor_angle inner_spin_rate outer_spin_rate
    variable
    oc_tol_cg oc_tol_ls oc_mnbrak_step oc_max_iter oc_cutoff 
    oc_cutoff_iter oc_var_save_iter oc_cg_min_step oc_max_brack_eval oc_max_brent_eval
    points_per_cycle
  }
  set data [split $data "\n"]
  foreach lst $data {
    if [string match #* [string trimleft $lst]] continue
    set nam [lindex $lst 0]
    set orignam $nam
    if ![string length $nam] continue
    if {[lsearch -exact $okpar $nam] == -1} {
      puts stderr "error: unknown name '$nam' in parameter section"
      exit
    }
    if ![string compare $nam variable] {
      set nam [lindex $lst 1]
      set lst [lrange $lst 1 9999]    
    }
    # ZT: added handling of 'string'
    if ![string compare $nam string] {
      set nam [lindex $lst 1]
      if {![string length $nam]} {
      	puts stderr "error: incomplete definition of string variable"
	exit
      }
      set lst [join [lrange $lst 1 9999]]    
    }
    # ZT: 'string' end of modification 
    if [info exists par($nam)] {
      puts stderr "error: '$nam' already exists in array par"
      exit
    }
    set ex [lrange $lst 1 9999]
    if {[lsearch -exact $allowsubst $orignam] != -1} {
      set origex $ex
      regsub -all {([A-Za-z_][A-Za-z0-9_]*\(?)} $ex { \1 } ex
      regsub -all {([0-9\.]+) ([Ee])} $ex {\1\2} ex
      regsub -all {([0-9\.]+[Ee]) ([+-][0-9\.]+)} $ex {\1\2} ex
      set ex "$ex "
      foreach n [lsort -command cmp_length [array names par]] {
        set i 0
        if {[string first " $n " $ex] != -1} {        
          if ![regsub -all " \($n\) " $ex { $par(\1) } ex] {
             puts stderr "error: unable to match $n in $ex"
             exit
          }
          set ex [subst "$ex"]
        }
      }
      if [catch {set ex [expr $ex]} res] {
        puts stderr "error: $res"
        exit  
      }
    }
    set par($nam) $ex
  }
}

proc fsimpson {{fitpar {}}} {
  global par _fd

  if [info exists par(verbose)] {
    if {[string index $par(verbose) 0] == 1}  {
      puts "Parameters"
      foreach k [lsort [array names par]] {
         puts [format "  %-20s %s" $k $par($k)]
      }
    }
  }
  spinsys_resolve $fitpar
  set par(tcalc) [lindex [time {uplevel #0 {set _fd [internalsimpson]}}] 0]
  set fd $_fd
  unset _fd
  return $fd
}

proc main {} {
   puts "There was no main section in the input file!"
}

proc pulseq {} {
   puts "There was no pulseq section in the input file!"
}







