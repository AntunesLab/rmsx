# grid_color_and_scale_NewTube_env_centered_palette_Xaxis.tcl
# Usage:
#   vmd -dispdev win -e grid_color_and_scale_NewTube_env_centered_palette_Xaxis.tcl \
#       -args file1.pdb file2.pdb ... [palette]
#
# Example:
#   vmd -dispdev win -e grid_color_and_scale_NewTube_env_centered_palette_Xaxis.tcl -args 1ubq.pdb 1ubq_copy.pdb turbo

# ----------------------------------------------------------------------
# identical normalize_list definition
proc normalize_list {input_list} {
    if {[llength $input_list] == 0} {return {}}
    set sorted_list [lsort -real $input_list]
    set list_length [llength $sorted_list]
    if {$list_length == 1} {return {5}}
   # set p5_idx  [expr {int($list_length * 0.05)}]
   # set p95_idx [expr {int($list_length * 0.95)}]
    #foreach v {p5_idx p95_idx} {
    #    if {[set $v] >= $list_length} {set $v [expr {$list_length - 1}]}
    #    if {[set $v] < 0} {set $v 0}
    #}
   # set scale_min [lindex $sorted_list $p5_idx]
    #set scale_max [lindex $sorted_list $p95_idx]
    # ... inside normalize_list ...
    set scale_min [lindex $sorted_list 0]
    set scale_max [lindex $sorted_list end]
    if {$scale_min == $scale_max} {return [lrepeat $list_length 5]}
    set range [expr {$scale_max - $scale_min}]
    set normalized_list {}
    #foreach val $input_list {
    #    if {$val < $scale_min} {set val $scale_min}
    #    if {$val > $scale_max} {set val $scale_max}
    #    set norm [expr {round(($val - $scale_min) / $range * 10.0)}]
    #    if {$norm < 0} {set norm 0}
    #    if {$norm > 10} {set norm 10}
    #    lappend normalized_list $norm
    #}
    # --- This is the NEW (continuous) code ---
    foreach val $input_list {
        if {$val < $scale_min} {set val $scale_min}
        if {$val > $scale_max} {set val $scale_max}
        
        # Calculate the continuous value from 0.0 to 10.0
        set norm [expr {(($val - $scale_min) / $range) * 10.0}]
        
        lappend normalized_list $norm
    }

    return $normalized_list
}

# ===== Representation + decoupled mapping controls ====================
# Geometry rep
set ::REP       "NewTube"   ;# or "NewCartoon"
set ::THICK      0.30       ;# baseline thickness (Å) [geom only]
set ::RES        120
set ::ASPECT     1.00
set ::SPLINE     0

# Thickness modulation pipeline (independent of color):
# We'll store per-residue thickness driver into 'user2'
set ::USERSCALE  1.0        ;# multiply 0..10 values
set ::USEROFFSET 2.0        ;# add offset before clamp
set ::MASK_OPACITY 0.30
if {[info exists env(RMSX_VMD_MASK_OPACITY)] && $env(RMSX_VMD_MASK_OPACITY) ne ""} {
    set ::MASK_OPACITY $env(RMSX_VMD_MASK_OPACITY)
}
set ::MASK_MATERIAL "Transparent"

# Color pipeline (independent of thickness):
# We'll color by 'user2' (can switch to Beta via hotkey if desired)
set ::COLORMETHOD "User2"    ;# "User2" | "Beta" | "User" | etc.
set ::COLMIN       0.0
set ::COLMAX      10.0

# Helper: apply only geometry
proc getRepIndices {molid} {
    set rep_ids {}
    set numreps [molinfo $molid get numreps]
    for {set repid 0} {$repid < $numreps} {incr repid} {
        lappend rep_ids $repid
    }
    return $rep_ids
}

proc applyGeom {} {
    foreach molid $::molList {
        foreach repid [getRepIndices $molid] {
            mol modstyle $repid $molid $::REP $::THICK $::RES $::ASPECT $::SPLINE
        }
    }
    puts [format {Geometry -> %s  THICK=%.2f  RES=%d  ASPECT=%.2f  SPLINE=%d} \
          $::REP $::THICK $::RES $::ASPECT $::SPLINE]
}
# Helper: apply only color/range
proc applyColor {} {
    if {$::COLMAX <= $::COLMIN} { set ::COLMAX [expr {$::COLMIN + 0.1}] }
    foreach molid $::molList {
        foreach repid [getRepIndices $molid] {
            mol modcolor $repid $molid $::COLORMETHOD
            mol scaleminmax $molid $repid $::COLMIN $::COLMAX
        }
    }
    puts [format {Color -> %s  range=[%0.2f..%0.2f]} $::COLORMETHOD $::COLMIN $::COLMAX]
}

proc maskValueIsTrue {value} {
    set lowered [string tolower [string trim $value]]
    return [expr {$lowered eq "1" || $lowered eq "true" || $lowered eq "yes" || $lowered eq "y" || $lowered eq "on"}]
}

proc loadMaskClauses {mask_csv_path} {
    if {![file exists $mask_csv_path]} {
        return {}
    }

    array unset chain_residues
    set handle [open $mask_csv_path r]
    gets $handle header_line

    while {[gets $handle line] >= 0} {
        if {[string trim $line] eq ""} {
            continue
        }
        set fields [split $line ","]
        if {[llength $fields] < 3} {
            continue
        }

        set residue_id [string trim [lindex $fields 0]]
        set chain_id [string trim [lindex $fields 1]]
        set masked_flag [string trim [lindex $fields 2]]
        if {$residue_id eq "" || ![maskValueIsTrue $masked_flag]} {
            continue
        }

        if {[info exists chain_residues($chain_id)]} {
            lappend chain_residues($chain_id) $residue_id
        } else {
            set chain_residues($chain_id) [list $residue_id]
        }
    }
    close $handle

    set clauses {}
    foreach chain_id [array names chain_residues] {
        set residue_ids [lsort -integer -unique $chain_residues($chain_id)]
        if {[llength $residue_ids] == 0} {
            continue
        }
        set residue_clause [join $residue_ids " "]
        if {$chain_id eq ""} {
            lappend clauses "resid $residue_clause"
        } else {
            lappend clauses "((segid $chain_id and resid $residue_clause) or (chain $chain_id and resid $residue_clause))"
        }
    }

    return $clauses
}

proc configureMaskMaterial {} {
    material change opacity $::MASK_MATERIAL $::MASK_OPACITY
    material change ambient $::MASK_MATERIAL 0.35
    material change diffuse $::MASK_MATERIAL 0.55
    material change specular $::MASK_MATERIAL 0.05
    material change shininess $::MASK_MATERIAL 0.03
    material change mirror $::MASK_MATERIAL 0.00
}

proc applyMaskedTransparency {mask_csv_path} {
    set mask_clauses [loadMaskClauses $mask_csv_path]
    if {[llength $mask_clauses] == 0} {
        return
    }

    set mask_selection [join $mask_clauses " or "]
    configureMaskMaterial
    catch {display depthsort on}

    foreach molid $::molList {
        if {[molinfo $molid get numreps] < 1} {
            continue
        }

        set current_style [molinfo $molid get "{representation 0}"]
        set current_selection [molinfo $molid get "{selection 0}"]
        set current_color [molinfo $molid get "{color 0}"]

        mol modselect 0 $molid "($current_selection) and not ($mask_selection)"
        mol addrep $molid
        set last_rep [expr {[molinfo $molid get numreps] - 1}]

        mol modselect $last_rep $molid "($current_selection) and ($mask_selection)"
        eval "mol modstyle $last_rep $molid $current_style"
        mol modcolor $last_rep $molid $current_color
        mol scaleminmax $molid $last_rep $::COLMIN $::COLMAX
        mol modmaterial $last_rep $molid $::MASK_MATERIAL
    }

    puts [format {Masked transparency applied: %s} $mask_csv_path]
}

# --- NEW PROCEDURE ---
# This procedure lays out the molecules and is called by the new hotkeys
proc setGridSpacing {new_spacing} {
    global spacing molList g_current_offsets
    if {[info exists g_current_offsets]} {
        puts "Undoing previous layout..."
        set idx 0
        foreach molid $molList {
            set last_dx [lindex $g_current_offsets $idx]
            set sel [atomselect $molid "all"]
            $sel moveby [list [expr {-$last_dx}] 0.0 0.0]
            $sel delete
            incr idx
        }
    }
    set spacing $new_spacing
    set g_current_offsets {}
    puts [format {Positioning molecules (spacing = %0.2f Å)} $spacing]

    set positions {}
    set count 0
    foreach molid $molList {
        set dx [expr {$count * $spacing}]
        lappend positions [list $dx 0.0 0.0]
        incr count
    }

    set totalX 0.0
    set npos [llength $positions]
    foreach pos $positions { set totalX [expr {$totalX + [lindex $pos 0]}] }
    set meanX [expr {$totalX / double($npos)}]
    set offsetX [expr {0.0 - $meanX}]

    set idx 0
    foreach molid $molList pos $positions {
        set dx [expr {[lindex $pos 0] + $offsetX}]
        lappend g_current_offsets $dx
        set sel [atomselect $molid "all"]
        $sel moveby [list $dx 0.0 0.0]
        $sel delete
        incr idx
    }
}
# --- END NEW PROCEDURE ---

# ----------------------------------------------------------------------
# Argument parsing and optional palette selection
if {$argc < 1} {
    puts "Usage: vmd -dispdev win -e grid_color_and_scale_NewTube_env_centered_palette_Xaxis.tcl -args pdb1 pdb2 ... [palette]"
    error "Need at least one PDB"
}
set lastArg [lindex $argv end]
if {[file exists $lastArg]} {
    set pdbFiles $argv
    set palette "viridis"
} else {
    set pdbFiles [lrange $argv 0 end-1]
    set palette $lastArg
    if {$palette eq ""} { set palette "viridis" }
}
puts [format {Palette: %s} $palette]

set N [llength $pdbFiles]
puts [format {Loading %d PDB files} $N]

# ----------------------------------------------------------------------
# Load and store molecule IDs
global molList
set molList {}
foreach f $pdbFiles {
    mol new $f type pdb waitfor all
    lappend molList [molinfo top get id]
}

# ----------------------------------------------------------------------
# Estimate spacing based on molecular size
global spacing
set sel0 [atomselect [lindex $molList 0] "all"]
set c0 [measure center $sel0]
set maxdist 0.0
foreach pos [$sel0 get {x y z}] {
    set d [vecdist $c0 $pos]
    if {$d > $maxdist} {set maxdist $d}
}
$sel0 delete
set spacing [expr {$maxdist * 2.2}]
puts [format {Auto spacing = %0.2f Å} $spacing]

# Arrange initially along X and center:
setGridSpacing $spacing

# ----------------------------------------------------------------------
# Per-residue average B-factors
puts "Computing per-residue average B-factors"
set allAvgs {}
set resIndexMap {}
foreach molid $molList {
    set selAll [atomselect $molid "all"]
    foreach resid [lsort -unique [$selAll get resid]] {
        set selR [atomselect $molid "resid $resid"]
        set b [$selR get beta]
        if {[llength $b]==0} {$selR delete; continue}
        set sum 0.0
        foreach v $b {set sum [expr {$sum + $v}]}
        set meanB [expr {$sum / double([llength $b])}]
        lappend allAvgs $meanB
        lappend resIndexMap [list $molid $resid]
        $selR delete
    }
    $selAll delete
}

# ----------------------------------------------------------------------
# Normalize (0..10)
set normValues [normalize_list $allAvgs]
set minNorm [lindex [lsort -real $normValues] 0]
set maxNorm [lindex [lsort -real $normValues] end]
puts [format {Color/scale normalized range (raw): %0.2f – %0.2f} $minNorm $maxNorm]

# ----------------------------------------------------------------------
# Choose per-atom fields:
#  - 'user'  -> thickness (independent)
#  - 'user2' -> color (independent)
# Also set env knobs so reps read thickness from user where supported.
set env(VMDMODULATERIBBON) user
set env(VMDMODULATENEWTUBE) user
set env(VMDMODULATENEWCARTOON) user

# Assign fields
for {set i 0} {$i<[llength $resIndexMap]} {incr i} {
    lassign [lindex $resIndexMap $i] molid resid
    set vraw [lindex $normValues $i]          ;# 0..10

    # thickness driver -> user (scaled independently)
    set t [expr {$::USERSCALE*$vraw + $::USEROFFSET}]
    if {$t < 0.0}  { set t 0.0 }
    if {$t > 10.0} { set t 10.0 }

    # color driver -> user2 (independent)
    set c $vraw

    set sel [atomselect $molid "resid $resid"]
    set n [$sel num]
    $sel set user  [lrepeat $n $t]
    $sel set user2 [lrepeat $n $c]
    $sel delete
}

# ----------------------------------------------------------------------
# Apply representation + color mapping (decoupled)
foreach molid $molList {
    mol delrep 0 $molid
    mol representation $::REP $::THICK $::RES $::ASPECT $::SPLINE
    mol addrep $molid
    mol modmaterial 0 $molid AOChalky

}
applyGeom
applyColor

# ----------------------------------------------------------------------
# Display style and palette handling
if {[catch {color scale method $palette} err]} {
    puts "Could not apply palette '$palette' ($err). Falling back to viridis."
    if {[catch {color scale method viridis}]} {
        color scale method BWR
        puts "Fallback palette: BWR"
    } else {
        puts "Fallback palette: viridis"
    }
} else {
    puts [format {Palette applied: %s} $palette]
}

display projection Orthographic
display rendermode GLSL
display shadows on
display depthcue off
display ambientocclusion on
axes location Off
stage location Off
color Display Background white
rotate x by 90

puts "Layout ready: X-axis centered, color=user2, size=user"

# ===============================================================
# rotateAllMols_3Dkeys.tcl
# ---------------------------------------------------------------
# Rotates all loaded molecules about their centers of mass
# along the X, Y, or Z axes using keyboard shortcuts.
#
# Keys:
#   u / i → rotate about X (+5° / -5°)
#   n / m → rotate about Y (+5° / -5°)
#   j / k → rotate about Z (+5° / -5°)
# ===============================================================

proc rotateAllMols {axis angle} {
    foreach mid [molinfo list] {
        if {[lsearch -exact [molinfo list] $mid] == -1} { continue }
        if {[molinfo $mid get numatoms] == 0 || [molinfo $mid get fixed]} { continue }
        set sel [atomselect $mid "all"]
        set com [measure center $sel weight mass]
        set R [transaxis $axis $angle]
        $sel moveby [vecscale -1 $com]
        $sel move $R
        $sel moveby $com
        $sel delete
    }
}
proc rotateAllMolsX {angle} { rotateAllMols x $angle }
proc rotateAllMolsY {angle} { rotateAllMols y $angle }
proc rotateAllMolsZ {angle} { rotateAllMols z $angle }

# ---------------------------------------------------------------
# Key bindings
# ---------------------------------------------------------------
# X-axis rotation
user add key {u} { rotateAllMolsX 5 }
user add key {i} { rotateAllMolsX -5 }
# Y-axis rotation
user add key {n} { rotateAllMolsY 5 }
user add key {m} { rotateAllMolsY -5 }
# Z-axis rotation
user add key {j} { rotateAllMolsZ 5 }
user add key {k} { rotateAllMolsZ -5 }
puts {Hotkeys: u/i ±X, n/m ±Y, j/k ±Z, +/- spacing, [ and ] thickness, 9/0/_/; size, ,/. color range, c color toggle}

# --- Grid spacing hotkeys ---
user add key {=} { global spacing ; setGridSpacing [expr {$spacing + 2.0}] }
user add key {+} { global spacing ; setGridSpacing [expr {$spacing + 2.0}] }
user add key {-} { global spacing ; setGridSpacing [expr {max($spacing - 2.0, 0.5)}] }

# --- Style hotkeys (decoupled) ---------------------------------------
# Base thickness (geometry only): [  ]
user add key {[} { set ::THICK  [expr {$::THICK + 0.05}] ; applyGeom }
user add key {]} { set ::THICK  [expr {max(0.05, $::THICK - 0.05)}] ; applyGeom }

# Optional: thickness modulation amplitude/offset (user) — geometry only
user add key {9} { set ::USERSCALE [expr {$::USERSCALE * 0.9}] ; puts [format {USERSCALE=%0.3f} $::USERSCALE] }
user add key {0} { set ::USERSCALE [expr {$::USERSCALE * 1.1}] ; puts [format {USERSCALE=%0.3f} $::USERSCALE] }
user add key {_} { set ::USEROFFSET [expr {$::USEROFFSET - 0.2}] ; puts [format {USEROFFSET=%0.3f} $::USEROFFSET] }
user add key {;} { set ::USEROFFSET [expr {$::USEROFFSET + 0.2}] ; puts [format {USEROFFSET=%0.3f} $::USEROFFSET] }

# Color range (color only): ,  .
user add key {,} { set ::COLMIN [expr {$::COLMIN + 0.5}] ; applyColor }
user add key {.} { set ::COLMAX [expr {$::COLMAX + 0.5}] ; applyColor }

# Toggle color method between User2 and Beta (color only)
user add key {c} {
    if {$::COLORMETHOD eq "User2"} { set ::COLORMETHOD "Beta" } else { set ::COLORMETHOD "User2" }
    applyColor
    puts [format {Coloring by %s} $::COLORMETHOD]
}

# ----------------------------------------------------------------------
